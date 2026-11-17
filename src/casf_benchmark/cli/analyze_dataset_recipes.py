"""Run the configured pharmacophore/diversity dataset-construction study."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import importlib.metadata
import json
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from casf_benchmark.analysis.dataset_recipes import (
    allocation_variants, bootstrap_interval, greedy_distance_clusters,
    paired_gains_losses, recipe_grid, recipe_statistics, stable_seed,
)
from casf_benchmark.analysis.pharmacophore import (
    atomic_json_dump, export_candidate_library, load_pickle, molecule_coordinate_hash,
    sha256, stable_smiles_id,
)
from casf_benchmark.cli.analyze_pharmacophores import (
    _match_group_exact_features, backend, load_config, radius_tag, write_frame_atomic,
)

DEFAULT_CONFIG = Path(__file__).resolve().parents[1] / "config/dataset_recipe_experiment.yaml"


def save_npz(path: Path, **values):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temporary.open("wb") as handle:
        np.savez_compressed(handle, **values)
    os.replace(temporary, path)


def context(config, panel):
    source = load_config(Path(config["panels"][panel]["source_config"]))
    old = Path(source["output_dir"])
    manifest = pd.read_csv(old / "query_manifest.csv", keep_default_na=False)
    protocol = Path(config['output_dir'])/'protocol.json'
    if protocol.is_file():
        frozen = json.loads(protocol.read_text())['panels'][panel]
        if (sha256(Path(source['_path'])) != frozen['source_config_sha256'] or
                sha256(old/'query_manifest.csv') != frozen['query_manifest_sha256'] or
                sha256(old/'query_freeze_manifest.json') != frozen['query_freeze_sha256']):
            raise ValueError(f'{panel}: frozen source configuration or query provenance changed')
    root = Path(config["output_dir"]) / panel
    root.mkdir(parents=True, exist_ok=True)
    return source, old, manifest, root


def pool_path(config, panel, name, source):
    replacement = config["panels"][panel].get("completed_pools", {}).get(name)
    if replacement:
        path = Path(config["output_dir"]) / replacement
        if path.is_file():
            return path
    method = config["panels"][panel]["methods"][name]
    if method not in source["methods"]:
        raise FileNotFoundError(f"Complete the {panel}/{name} generation pool first")
    return Path(source["methods"][method]["path"])


def freeze(config, config_path):
    root = Path(config["output_dir"])
    root.mkdir(parents=True, exist_ok=True)
    protocol_path = root / "protocol.json"
    previous = json.loads(protocol_path.read_text()) if protocol_path.is_file() else None
    if previous and previous['config_sha256'] != sha256(config_path):
        raise ValueError("A frozen protocol already exists with a different configuration")
    rows, provenance = [], {}
    for panel in config["panels"]:
        source, old, manifest, _ = context(config, panel)
        provenance[panel] = {"source_config_sha256": sha256(Path(source["_path"])),
                             "query_manifest_sha256": sha256(old / "query_manifest.csv"),
                             "query_freeze_sha256": sha256(old / "query_freeze_manifest.json"),
                             "molecule_count": int(manifest.expected_smiles.nunique()),
                             "eligible_queries": int(manifest.eligible.sum())}
        for name, method in config["panels"][panel]["methods"].items():
            path = source["methods"].get(method, {}).get("path")
            if path:
                rows.append({"panel": panel, "source": name, "method": method,
                             "path": path, "sha256": sha256(Path(path))})
    inventory = pd.DataFrame(rows)
    inventory_path = root / "immutable_source_inventory.csv"
    if inventory_path.is_file() and not pd.read_csv(inventory_path).equals(inventory):
        raise ValueError('Original source inventory changed after freezing')
    if previous and previous['panels'] != provenance:
        raise ValueError('Frozen query/reference provenance changed')
    write_frame_atomic(inventory, inventory_path)
    recipes = [{"recipe": r.name, **{s: dict(r.weights).get(s, 0) for s in
                                   [*config["sources"], "rdkit_raw", "rdkit_minimized"]}} for r in recipe_grid(config['sources'])]
    write_frame_atomic(pd.DataFrame(recipes), root / "recipes.csv")
    payload = {"status": "frozen_exploratory_protocol", "configuration": config,
               "config_sha256": sha256(config_path), "panels": provenance,
               "recipe_count": len(recipes), "primary_gate": "posebusters",
               "rejected_candidates_consume_budget": True,
               "reference_weighting_sensitivity": "equal reference-geometry/contact-pattern groups",
               "reference_group_definition": "greedy native heavy-atom RMSD < radius, identical feature types, and an injective feature assignment within radius after the native atom fit",
               "uncertainty_unit": "molecule; repeated library subsets are not independent molecules",
               "generation_shortfalls": "quota capped at source availability; equal-availability sensitivity includes every molecule",
               "selection": "uniform subsets; no query-informed candidate selection"}
    atomic_json_dump(protocol_path, payload)
    print(json.dumps({"recipe_count": len(recipes), "panels": provenance}, indent=2), flush=True)


def prepare_source(config, panel, name, smiles, query_group, source, old, root):
    from rdkit import Chem
    from casf_benchmark.generation.conformer_sets import posebusters_geometry_passes

    path = pool_path(config, panel, name, source)
    mols = load_pickle(path)[smiles]
    sid = stable_smiles_id(smiles)
    method = config["panels"][panel]["methods"][name]
    hashes = [molecule_coordinate_hash(mol) for mol in mols]
    destination = root / "prepared" / sid / f"{name}.npz"
    query_ids = query_group.query_id.astype(str).tolist()
    radii = config["radii"]
    previous_count, previous = 0, None
    if destination.is_file():
        with np.load(destination, allow_pickle=True) as stored:
            previous = {k: stored[k] for k in stored.files}
        old_hashes = list(map(str, previous["coordinate_hashes"]))
        if old_hashes != hashes[:len(old_hashes)] or list(previous["query_ids"]) != query_ids:
            raise ValueError(f"Prepared source prefix/query changed: {panel}/{name}/{sid}")
        previous_count = len(old_hashes)
        if previous_count == len(mols):
            metadata_path = destination.with_suffix('.json')
            metadata = json.loads(metadata_path.read_text()) if metadata_path.is_file() else {
                'candidate_count': len(mols), 'reused_prefix_count': len(mols),
                'new_candidates_assessed': 0, 'pb_passes': int(previous['passes'].sum()),
                'matrix_sha256': sha256(destination), 'new_pb_checks': [],
                'query_manifest_sha256': sha256(old/'query_manifest.csv')}
            metadata.update(source_path=str(path.resolve()), source_sha256=sha256(path))
            atomic_json_dump(metadata_path, metadata)
            return mols
    elif method in source["methods"]:
        records_path = old / "candidate_libraries" / method / f"{sid}.records.json"
        pb_path = old / "posebusters" / method / f"{sid}.npz"
        if pb_path.is_file() and (records_path.is_file() or not len(query_group)):
            records = json.loads(records_path.read_text()) if records_path.is_file() else []
            old_hashes = ([r["coordinate_sha256"] for r in records] if records else
                          [molecule_coordinate_hash(m) for m in load_pickle(Path(source['methods'][method]['path']))[smiles]])
            if old_hashes != hashes[:len(old_hashes)]:
                raise ValueError(f"Archived source coordinates disagree: {panel}/{name}/{sid}")
            previous_count = len(old_hashes)
            previous = {"passes": np.load(pb_path)["passes"].astype(bool)}
            for radius in radii:
                if not len(query_group):
                    previous[f"matches_{radius_tag(radius)}"] = np.zeros((0,previous_count),bool)
                    continue
                p = old / "matches" / method / sid / f"{radius_tag(radius)}.npz"
                with np.load(p, allow_pickle=True) as stored:
                    if list(stored["query_ids"]) != query_ids or stored["matches"].shape != (len(query_ids), previous_count):
                        raise ValueError("Archived match matrix correspondence changed")
                    if list(stored['candidate_ids']) != [r['candidate_id'] for r in records]:
                        raise ValueError("Archived candidate ordering changed")
                    previous[f"matches_{radius_tag(radius)}"] = stored["matches"].astype(bool)
    added_count = len(mols) - previous_count
    if added_count < 0:
        raise ValueError("A completed parent pool lost source candidates")
    passes = np.zeros(len(mols), dtype=bool)
    if previous_count:
        if len(previous["passes"]) != previous_count:
            raise ValueError("Archived PoseBusters mask does not align with candidates")
        passes[:previous_count] = previous["passes"]
    checks = []
    reference = Chem.MolFromSmiles(smiles)
    check_mols, check_indices = [], []
    for index in range(previous_count, len(mols)):
        candidate = Chem.RemoveHs(Chem.Mol(mols[index]))
        match = candidate.GetSubstructMatch(reference)
        if len(match) != reference.GetNumAtoms() or candidate.GetNumAtoms() != reference.GetNumAtoms():
            checks.append({"candidate_index": index, "identity_mapping": False})
            continue
        check_mols.append(Chem.RenumberAtoms(candidate, list(match)))
        check_indices.append(index)
    if check_mols:
        result = posebusters_geometry_passes(check_mols, reference,
                                            max_workers=int(config["workers"]), energy_num_threads=1)
        passes[check_indices] = result.passes
        checks.extend({"candidate_index": i, **row} for i, row in zip(check_indices, result.check_rows, strict=True))
    values = {"passes": passes, "coordinate_hashes": np.asarray(hashes, dtype=object),
              "query_ids": np.asarray(query_ids, dtype=object)}
    new_matrices = {radius: np.zeros((len(query_ids), added_count), dtype=bool) for radius in radii}
    if added_count and len(query_ids):
        working = root / "feature_work"
        working.mkdir(parents=True, exist_ok=True)
        query_link = working / "queries"
        if not query_link.exists():
            query_link.symlink_to((old / "queries").resolve(), target_is_directory=True)
        settings = dict(source)
        settings["output_dir"] = str(working)
        settings["backend"] = {**source["backend"], "workers": config["workers"]}
        records = export_candidate_library(method=method, smiles=smiles, mols=mols,
                                           sdf_path=working / "candidate_libraries" / method / f"{sid}.sdf")
        _match_group_exact_features(config=settings, engine=backend(settings), method=method,
                                    smiles=smiles, mols=mols, candidate_records=records[previous_count:],
                                    query_group=query_group, radii=radii)
        for radius in radii:
            p = working / "matches" / method / sid / f"{radius_tag(radius)}.npz"
            with np.load(p, allow_pickle=True) as stored:
                if (list(stored['query_ids']) != query_ids or
                        list(stored['candidate_ids']) != [r['candidate_id'] for r in records[previous_count:]] or
                        stored['matches'].shape != (len(query_ids),added_count)):
                    raise ValueError('Appended match matrix correspondence changed')
                new_matrices[radius] = stored["matches"].astype(bool)
    for radius in radii:
        key = f"matches_{radius_tag(radius)}"
        prefix = previous[key] if previous_count and key in previous else np.zeros((len(query_ids), previous_count), bool)
        values[key] = np.concatenate([prefix, new_matrices[radius]], axis=1)
    for a, b in zip(radii, radii[1:]):
        if np.any(values[f"matches_{radius_tag(a)}"] & ~values[f"matches_{radius_tag(b)}"]):
            raise ValueError("Increasing tolerance lost matches")
    if hashes != [molecule_coordinate_hash(mol) for mol in mols]:
        raise ValueError("Analysis mutated a source conformer")
    save_npz(destination, **values)
    atomic_json_dump(destination.with_suffix(".json"), {
        "source_path": str(path.resolve()), "source_sha256": sha256(path),
        "candidate_count": len(mols), "reused_prefix_count": previous_count,
        "new_candidates_assessed": added_count, "pb_passes": int(passes.sum()),
        "new_pb_checks": checks, "matrix_sha256": sha256(destination),
        "query_manifest_sha256": sha256(old / "query_manifest.csv")})
    print(f"prepared {panel}/{name}/{sid}: {len(mols)} candidates, {passes.sum()} PB, {previous_count} reused", flush=True)
    return mols


def reference_groups(config, group, old):
    """Reference-only redundancy grouping, retaining different contact patterns."""
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign
    from scipy.optimize import linear_sum_assignment

    radius = config["reference_cluster_radius"]
    records, centers, assignments = [], [], []
    for row in group.itertuples(index=False):
        mol = next(m for m in Chem.SDMolSupplier(row.raw_sdf_path, removeHs=True) if m is not None)
        query = json.loads((old / "queries" / radius_tag(config["radii"][0]) / f"{row.query_id}.json").read_text())
        types = np.array([point["name"] for point in query["points"]])
        xyz = np.array([[point[k] for k in 'xyz'] for point in query["points"]])
        record = (mol, types, xyz)
        assignment = None
        for index, center in enumerate(centers):
            cmol, ctypes, cxyz = center
            if sorted(types) != sorted(ctypes):
                continue
            rmsd, transform, _ = rdMolAlign.GetBestAlignmentTransform(mol, cmol)
            if rmsd >= radius:
                continue
            moved = xyz @ transform[:3,:3].T + transform[:3,3]
            bad = (types[:,None] != ctypes[None,:]) | (np.linalg.norm(moved[:,None,:]-cxyz[None,:,:],axis=2) > radius)
            a,b = linear_sum_assignment(bad.astype(int))
            if not np.any(bad[a,b]):
                assignment = index
                break
        if assignment is None:
            assignment = len(centers)
            centers.append(record)
        assignments.append(assignment)
        mapping_path = old / "reference_maps" / f"{row.query_id}.json"
        mapping = json.loads(mapping_path.read_text()) if mapping_path.is_file() else {}
        native = mapping.get("native_assembly", {})
        occupancies = [a.get("occupancy", 1) for a in native.get("raw_to_assembly_atom", [])]
        records.append({"query_id": row.query_id, "reference_group": assignment,
                        "pdb_id": row.pdb_id, "name": row.name, "feature_count": len(types),
                        "minimum_ligand_occupancy": min(occupancies) if occupancies else np.nan,
                        "reference_geometry_status": row.reference_geometry_failure_reasons or 'passed_archived_checks',
                        "identity_matches": row.identity_matches_audit, "stereo_matches": row.stereo_matches_audit})
    counts = np.bincount(assignments) if assignments else np.array([])
    for row in records:
        row["query_weight"] = 1 / len(records)
        row["group_weight"] = 1 / (len(counts) * counts[row["reference_group"]])
    return pd.DataFrame(records)


def prepare(config, panel, index, sources=None):
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.*')
    source, old, manifest, root = context(config, panel)
    smiles_list = sorted(manifest.expected_smiles.unique())
    smiles = smiles_list[index]
    group = manifest[(manifest.expected_smiles == smiles) & manifest.eligible].sort_values('reference_index')
    sid = stable_smiles_id(smiles)
    for name in (sources or config['panels'][panel]['methods']):
        if name not in config['panels'][panel]['methods']:
            raise ValueError(f'Unknown source {name}')
        prepare_source(config, panel, name, smiles, group, source, old, root)
    path = root / 'prepared' / sid / 'references.csv'
    if len(group) and not path.is_file():
        write_frame_atomic(reference_groups(config, group, old), path)
    elif not len(group):
        write_frame_atomic(pd.DataFrame(columns=['query_id','reference_group','query_weight','group_weight']), path)
    missing = [name for name in config['panels'][panel]['methods']
               if not (root/'prepared'/sid/f'{name}.npz').is_file()]
    atomic_json_dump(root / 'prepared' / sid / 'complete.json',
                     {'status':'partial' if missing else 'complete','missing_sources':missing,
                      'smiles':smiles,'index':index,'eligible_queries':len(group),'config_sha256':config['_sha256']})


def load_prepared(config, panel, smiles):
    root = Path(config['output_dir']) / panel / 'prepared' / stable_smiles_id(smiles)
    data = {}
    for name in config['panels'][panel]['methods']:
        with np.load(root / f'{name}.npz', allow_pickle=True) as archive:
            data[name] = {k:archive[k] for k in archive.files}
    return data, pd.read_csv(root / 'references.csv')


def distance_matrix(mols, path, signature, workers):
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign
    if path.is_file():
        with np.load(path) as archive:
            if str(archive['signature']) == signature:
                return archive['distances']
    started = time.perf_counter()
    if not mols:
        distances = np.empty((0,0),np.float32)
    else:
        template = Chem.RemoveHs(Chem.Mol(mols[0]))
        combined = Chem.Mol(template)
        combined.RemoveAllConformers()
        for mol in mols:
            mol = Chem.RemoveHs(Chem.Mol(mol))
            atom_map = mol.GetSubstructMatch(template)
            if len(atom_map) != template.GetNumAtoms() or mol.GetNumAtoms() != template.GetNumAtoms():
                raise ValueError('Diversity parent pools have incompatible molecular graphs')
            ordered = Chem.RenumberAtoms(mol, atom_map)
            combined.AddConformer(Chem.Conformer(ordered.GetConformer()),assignId=True)
        lower = np.asarray(rdMolAlign.GetAllConformerBestRMS(combined, numThreads=int(workers)), dtype=np.float32)
        distances = np.zeros((len(mols),len(mols)),np.float32)
        rows,cols = np.tril_indices(len(mols),-1)
        distances[rows,cols] = lower
        distances[cols,rows] = lower
    if not np.isfinite(distances).all():
        raise ValueError('Nonfinite diversity distance')
    save_npz(path, distances=distances, signature=np.array(signature))
    print(f'distance cache {path.name}: {len(mols)} conformers in {time.perf_counter()-started:.1f}s',flush=True)
    return distances


def diversity(config, panel, index):
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.*')
    source, old, manifest, root = context(config,panel)
    smiles = sorted(manifest.expected_smiles.unique())[index]
    sid = stable_smiles_id(smiles)
    data, references = load_prepared(config,panel,smiles)
    populations = {name:len(value['passes']) for name,value in data.items()}
    pools = {name:load_pickle(pool_path(config,panel,name,source))[smiles] for name in data}
    for name in data:
        if [molecule_coordinate_hash(m) for m in pools[name]] != list(data[name]['coordinate_hashes']):
            raise ValueError('Diversity pool changed after preparation')
    rows = []
    for family in [config['sources'],['rdkit_raw'],['rdkit_minimized']]:
        offsets, sequence, all_hashes = {}, [], []
        for name in family:
            offsets[name] = len(sequence)
            sequence.extend(pools[name])
            all_hashes.extend(data[name]['coordinate_hashes'])
        signature = hashlib.sha256('|'.join(all_hashes).encode()).hexdigest()
        matrix = distance_matrix(sequence,root/'distances'/sid/('_'.join(family)+'.npz'),signature,config['workers'])
        for recipe in recipe_grid(config['sources']):
            if not set(dict(recipe.weights)).issubset(family):
                continue
            for scenario in ['requested','equal_availability']:
                for budget in config['budgets']:
                    actual_budget = budget if scenario=='requested' else min(budget,min(populations.values()))
                    variants = allocation_variants(recipe,actual_budget)
                    for draw in range(config['sampling_repeats'] * len(variants)):
                        repeat = draw // len(variants)
                        # Shared per-source permutations provide nested draws for all recipes.
                        selected = []
                        variant_index = draw % len(variants)
                        variant = variants[variant_index]
                        for name,n in variant.items():
                            rng = np.random.default_rng(stable_seed(config['seed'],panel,smiles,name,repeat,'selection'))
                            indices = rng.permutation(populations[name])[:min(n,populations[name])]
                            selected.extend((offsets[name]+int(i),bool(data[name]['passes'][i])) for i in indices)
                        rng = np.random.default_rng(stable_seed(config['seed'],panel,smiles,repeat,'order'))
                        order = rng.permutation(len(selected))
                        for gate in ['raw','posebusters']:
                            indices = np.array([selected[i][0] for i in order if gate=='raw' or selected[i][1]],int)
                            for radius in config['cluster_radii']:
                                metrics = greedy_distance_clusters(matrix,indices,radius)
                                rows.append({'panel':panel,'smiles':smiles,'recipe':recipe.name,'scenario':scenario,
                                             'budget':budget,'sampling_budget':actual_budget,'repeat':repeat,
                                             'allocation_variant':variant_index,'gate':gate,
                                             'cluster_radius':radius,'realized_draws':len(selected),'valid_count':len(indices),
                                             'eligible_queries':len(references),**metrics})
    write_frame_atomic(pd.DataFrame(rows),root/'diversity'/f'{sid}.csv')
    atomic_json_dump(root/'diversity'/f'{sid}.metadata.json', {
        'config_sha256': config['_sha256'], 'smiles': smiles,
        'prepared_sha256': {name: sha256(root/'prepared'/sid/f'{name}.npz') for name in data},
        'diversity_csv_sha256': sha256(root/'diversity'/f'{sid}.csv')})
    print(f'diversity {panel}/{index} complete: {len(rows)} rows',flush=True)


def probabilities(config,panel):
    source,old,manifest,root = context(config,panel)
    rows,query_rows,gain_rows = [],[],[]
    for smiles, group in manifest[manifest.eligible].groupby('expected_smiles',sort=True):
        data, refs = load_prepared(config,panel,smiles)
        if list(refs.query_id) != group.sort_values('reference_index').query_id.tolist():
            raise ValueError('Reference weighting/order mismatch')
        populations = {name:len(d['passes']) for name,d in data.items()}
        valid_counts = {name:int(d['passes'].sum()) for name,d in data.items()}
        for gate in ['raw','posebusters']:
            counts = populations if gate=='raw' else valid_counts
            for radius in config['radii']:
                hits = {name:d[f'matches_{radius_tag(radius)}'][:,np.ones(populations[name],bool) if gate=='raw' else d['passes']].sum(axis=1) for name,d in data.items()}
                for scenario in ['requested','equal_availability']:
                    for budget in config['budgets']:
                        n = budget if scenario=='requested' else min(budget,min(populations.values()))
                        for recipe in recipe_grid(config['sources']):
                            result = recipe_statistics(recipe,n,populations,hits,counts)
                            for weighting,column in [('queries','query_weight'),('reference_groups','group_weight')]:
                                weights = refs[column].to_numpy(float)
                                rows.append({'panel':panel,'smiles':smiles,'recipe':recipe.name,'gate':gate,'radius':radius,
                                             'scenario':scenario,'budget':budget,'sampling_budget':n,'weighting':weighting,
                                             'query_count':len(refs),'reference_groups':refs.reference_group.nunique(),
                                             'coverage':float(weights@result['coverage']),'support':float(weights@result['support']),
                                             'realized_draws':result['realized_draws'],'expected_valid':result['expected_valid'],
                                             'complete_requested_budget':result['realized_draws']==budget})
                            if gate=='posebusters' and radius==config['radii'][0] and budget==max(config['budgets']) and scenario=='requested':
                                for i,row in enumerate(refs.itertuples(index=False)):
                                    query_rows.append({'panel':panel,'smiles':smiles,'recipe':recipe.name,'query_id':row.query_id,
                                                       'coverage':result['coverage'][i],'support':result['support'][i],
                                                       'feature_count':row.feature_count,'query_weight':row.query_weight,
                                                       'group_weight':row.group_weight})
                                for comparator in recipe_grid(config['sources']):
                                    if comparator.name not in [*config['sources'],'rdkit_minimized'] or comparator.name==recipe.name:
                                        continue
                                    gained,lost = paired_gains_losses(recipe,comparator,n,populations,hits)
                                    gain_rows.append({'panel':panel,'smiles':smiles,'recipe':recipe.name,'comparator':comparator.name,
                                                      'gained_queries':float(gained.sum()),'lost_queries':float(lost.sum()),
                                                      'weighted_gain':float(refs.query_weight.to_numpy()@gained),
                                                      'weighted_loss':float(refs.query_weight.to_numpy()@lost)})
    write_frame_atomic(pd.DataFrame(rows),root/'per_molecule_coverage.csv')
    write_frame_atomic(pd.DataFrame(query_rows),root/'per_query_support.csv')
    write_frame_atomic(pd.DataFrame(gain_rows),root/'gains_losses_per_molecule.csv')
    return pd.DataFrame(rows)


def diagnostic_summaries(config, panel, coverage):
    """Exploratory strata; no molecule/query filtering or recipe tuning."""
    from rdkit import Chem
    from rdkit.Chem import Lipinski

    source, old, manifest, root = context(config, panel)
    overlap_path = old / 'training_overlap/exact_identity_overlap.csv'
    overlap = pd.read_csv(overlap_path).set_index('expected_smiles') if overlap_path.is_file() else None
    if overlap is not None and not overlap.index.is_unique:
        raise ValueError('Training-overlap audit has duplicate identities')
    descriptors, references = [], []
    for smiles, group in manifest.groupby('expected_smiles', sort=True):
        mol = Chem.MolFromSmiles(smiles)
        rotatable = Lipinski.NumRotatableBonds(mol)
        status = 'unknown'
        if overlap is not None and smiles in overlap.index:
            status = 'exact_overlap' if bool(overlap.loc[smiles, 'exact_identity_overlap']) else 'no_exact_overlap'
        descriptors.append({'smiles': smiles, 'name': group.iloc[0]['name'],
                            'heavy_atoms': mol.GetNumHeavyAtoms(), 'rotatable_bonds': rotatable,
                            'flexibility': '0-3' if rotatable <= 3 else '4-7' if rotatable <= 7 else '8+',
                            'eligible_queries': int(group.eligible.sum()), 'qwen_training_overlap': status})
        path = root / 'prepared' / stable_smiles_id(smiles) / 'references.csv'
        refs = pd.read_csv(path)
        if len(refs):
            references.append(refs.assign(smiles=smiles, panel=panel))
    descriptors = pd.DataFrame(descriptors)
    write_frame_atomic(descriptors, root / 'molecule_descriptors.csv')
    refs = pd.concat(references, ignore_index=True)
    write_frame_atomic(refs, root / 'reference_group_audit.csv')
    atomic_json_dump(root / 'reference_quality_summary.json', {
        'eligible_queries': len(refs), 'molecules': int(descriptors.eligible_queries.gt(0).sum()),
        'molecules_without_queries': int(descriptors.eligible_queries.eq(0).sum()),
        'redundancy_groups': int(refs.groupby('smiles').reference_group.nunique().sum()),
        'queries_with_minimum_occupancy_below_one': int(refs.minimum_ligand_occupancy.lt(1).sum()),
        'queries_with_unknown_occupancy': int(refs.minimum_ligand_occupancy.isna().sum()),
        'reference_geometry_status_counts': refs.reference_geometry_status.value_counts().to_dict(),
        'training_overlap_audit': str(overlap_path) if overlap is not None else None,
        'training_overlap_audit_sha256': sha256(overlap_path) if overlap is not None else None,
        'training_overlap_caution': 'Exact Qwen identity audit only; no_exact_overlap does not mean scaffold-disjoint or unseen to every generator.'})
    primary = coverage[(coverage.gate == 'posebusters') & (coverage.budget == 1000) &
                       (coverage.scenario == 'requested') & (coverage.weighting == 'queries')]
    strata = []
    joined = primary.merge(descriptors, on='smiles', validate='many_to_one')
    for column in ['flexibility', 'qwen_training_overlap']:
        for (recipe, radius, label), group in joined.groupby(['recipe', 'radius', column], sort=False):
            lo, hi = bootstrap_interval(group.coverage.to_numpy(), replicates=config['bootstrap_replicates'], seed=config['seed'])
            strata.append({'stratifier': column, 'stratum': label, 'recipe': recipe, 'radius': radius,
                           'coverage': group.coverage.mean(), 'ci_low': lo, 'ci_high': hi, 'molecules': len(group)})
    queries = pd.read_csv(root / 'per_query_support.csv')
    queries = queries.merge(refs[['query_id', 'reference_geometry_status', 'minimum_ligand_occupancy']],
                            on='query_id', validate='many_to_one')
    queries['feature_stratum'] = pd.cut(queries.feature_count, [0, 4, 6, np.inf], labels=['3-4', '5-6', '7+'])
    queries['native_geometry'] = np.select(
        [queries.reference_geometry_status.eq('passed_archived_checks'),
         queries.reference_geometry_status.eq('not_reference_posebusters_screened')],
        ['passed_archived_checks', 'not_screened'], default='flagged_archived_checks')
    queries['native_occupancy'] = np.where(queries.minimum_ligand_occupancy.isna(), 'unknown',
        np.where(queries.minimum_ligand_occupancy.lt(1), 'below_one', 'at_least_one'))
    # Query strata are still macro-averaged over molecules, not pooled over queries.
    for column in ['feature_stratum', 'native_geometry', 'native_occupancy']:
        grouped = queries.groupby(['smiles', 'recipe', column], observed=True).agg(
            coverage=('coverage', 'mean'), support=('support', 'mean'), queries=('query_id', 'size')).reset_index()
        write_frame_atomic(grouped, root / f'{column}_per_molecule.csv')
        for (recipe, label), group in grouped.groupby(['recipe', column], observed=True, sort=False):
            lo, hi = bootstrap_interval(group.coverage.to_numpy(), replicates=config['bootstrap_replicates'], seed=config['seed'])
            strata.append({'stratifier': column, 'stratum': str(label), 'recipe': recipe, 'radius': .5,
                           'coverage': group.coverage.mean(), 'ci_low': lo, 'ci_high': hi, 'molecules': len(group)})
    write_frame_atomic(pd.DataFrame(strata), root / 'stratified_coverage.csv')


def summarize(config,panel,allow_partial=False):
    source,old,manifest,root = context(config,panel)
    expected = [stable_smiles_id(s) for s in sorted(manifest.expected_smiles.unique())]
    missing = [sid for sid in expected if not (root/'prepared'/sid/'complete.json').is_file()
               or json.loads((root/'prepared'/sid/'complete.json').read_text()).get('status')!='complete']
    if missing:
        raise RuntimeError(f'{panel}: {len(missing)} molecules lack prepared results')
    pool_hashes = {name: sha256(pool_path(config,panel,name,source)) for name in config['panels'][panel]['methods']}
    for sid in expected:
        for name in pool_hashes:
            metadata = json.loads((root/'prepared'/sid/f'{name}.json').read_text())
            if metadata['source_sha256'] != pool_hashes[name]:
                raise RuntimeError(f'{panel}/{sid}/{name}: parent pool changed; rerun prepare')
    coverage = probabilities(config,panel)
    diagnostic_summaries(config, panel, coverage)
    summaries = []
    keys = ['panel','recipe','gate','radius','scenario','budget','weighting']
    for key,group in coverage.groupby(keys,sort=False):
        lo,hi = bootstrap_interval(group.coverage.to_numpy(),replicates=config['bootstrap_replicates'],seed=config['seed'])
        summaries.append(dict(zip(keys,key))|{'coverage':group.coverage.mean(),'ci_low':lo,'ci_high':hi,
                                            'molecule_count':len(group),'query_count':group.query_count.sum(),
                                            'support':group.support.mean(),'realized_draws':group.realized_draws.mean(),
                                            'expected_valid':group.expected_valid.mean(),
                                            'complete_budget_molecules':int(group.complete_requested_budget.sum())})
    aggregate = pd.DataFrame(summaries)
    write_frame_atomic(aggregate,root/'coverage_summary.csv')
    primary = coverage[(coverage.gate=='posebusters')&(coverage.radius==.5)&(coverage.budget==1000)&
                       (coverage.scenario=='requested')&(coverage.weighting=='queries')]
    pivot = primary.pivot(index='smiles',columns='recipe',values='coverage')
    paired = []
    for left,right in itertools.combinations(pivot.columns,2):
        difference = (pivot[left]-pivot[right]).to_numpy()
        lo,hi = bootstrap_interval(difference,replicates=config['bootstrap_replicates'],seed=config['seed'])
        paired.append({'left':left,'right':right,'difference':float(difference.mean()),'ci_low':lo,'ci_high':hi,
                       'molecules':len(difference),'intervals':'exploratory, unadjusted for recipe search'})
    write_frame_atomic(pd.DataFrame(paired),root/'paired_recipe_differences.csv')
    gains = pd.read_csv(root/'gains_losses_per_molecule.csv').groupby(['recipe','comparator']).agg(
        gained_queries=('gained_queries','sum'),lost_queries=('lost_queries','sum'),
        weighted_gain=('weighted_gain','mean'),weighted_loss=('weighted_loss','mean')).reset_index()
    write_frame_atomic(gains,root/'gains_losses_summary.csv')
    diversity_missing = []
    for sid in expected:
        meta_path = root/'diversity'/f'{sid}.metadata.json'
        metadata = json.loads(meta_path.read_text()) if meta_path.is_file() else {}
        current = {name: sha256(root/'prepared'/sid/f'{name}.npz') for name in pool_hashes}
        if (not (root/'diversity'/f'{sid}.csv').is_file() or metadata.get('config_sha256') != config['_sha256']
                or metadata.get('prepared_sha256') != current):
            diversity_missing.append(sid)
    if diversity_missing and not allow_partial:
        raise RuntimeError(f'{panel}: coverage finished, but {len(diversity_missing)} diversity results remain')
    paths = [root/'diversity'/f'{sid}.csv' for sid in expected if sid not in diversity_missing]
    diversity_summary = pd.DataFrame()
    if paths:
        draws = pd.concat([pd.read_csv(path) for path in paths],ignore_index=True)
        molkeys = ['panel','smiles','recipe','scenario','budget','gate','cluster_radius','eligible_queries']
        per_molecule = draws.groupby(molkeys,dropna=False).agg(clusters=('clusters','mean'),
            largest_fraction=('largest_fraction','mean'),effective_clusters=('effective_clusters','mean'),
            valid_count=('valid_count','mean'),realized_draws=('realized_draws','mean')).reset_index()
        write_frame_atomic(per_molecule,root/'per_molecule_diversity.csv')
        # Include no-query molecules here only, not in the joint coverage endpoint.
        support_keys = ['recipe','scenario','budget','gate','cluster_radius']
        support = per_molecule.groupby(support_keys).agg(
            clusters=('clusters','mean'), largest_fraction=('largest_fraction','mean'),
            effective_clusters=('effective_clusters','mean'), valid_count=('valid_count','mean'),
            molecule_count=('smiles','nunique')).reset_index()
        write_frame_atomic(support, root/'diversity_all_molecules_supporting.csv')
        divkeys = ['panel','recipe','scenario','budget','gate','cluster_radius']
        joint = per_molecule[per_molecule.eligible_queries>0]
        divrows = []
        for key,group in joint.groupby(divkeys,sort=False):
            lo,hi = bootstrap_interval(group.clusters.to_numpy(),replicates=config['bootstrap_replicates'],seed=config['seed'])
            divrows.append(dict(zip(divkeys,key))|{'clusters':group.clusters.mean(),'clusters_ci_low':lo,'clusters_ci_high':hi,
                'largest_fraction':group.largest_fraction.mean(),'effective_clusters':group.effective_clusters.mean(),
                'valid_count':group.valid_count.mean(),'molecule_count':len(group)})
        diversity_summary = pd.DataFrame(divrows)
        write_frame_atomic(diversity_summary,root/'diversity_summary.csv')
        comparisons = []
        primary_diversity = joint[(joint.gate=='posebusters') & (joint.budget==1000) &
                                 (joint.scenario=='requested') & (joint.cluster_radius==.5)]
        for metric in ['clusters', 'largest_fraction', 'effective_clusters', 'valid_count']:
            values = primary_diversity.pivot(index='smiles', columns='recipe', values=metric)
            for left,right in itertools.combinations(values.columns,2):
                difference = (values[left]-values[right]).to_numpy()
                lo,hi = bootstrap_interval(difference,replicates=config['bootstrap_replicates'],seed=config['seed'])
                comparisons.append({'left':left, 'right':right, 'metric':metric,
                                    'difference':float(np.nanmean(difference)), 'ci_low':lo, 'ci_high':hi,
                                    'molecules':int(np.isfinite(difference).sum()),
                                    'intervals':'exploratory, unadjusted for recipe search'})
        write_frame_atomic(pd.DataFrame(comparisons), root/'paired_diversity_differences.csv')
    endpoint = aggregate[(aggregate.gate=='posebusters')&(aggregate.budget==1000)&
                         (aggregate.scenario=='requested')&(aggregate.weighting=='queries')]
    table = endpoint.pivot(index='recipe',columns='radius',values='coverage').rename(columns={.5:'coverage_0p5',.75:'coverage_0p75',1.:'coverage_1p0'})
    meta = endpoint[endpoint.radius==.5].set_index('recipe')[['ci_low','ci_high','expected_valid','realized_draws','complete_budget_molecules','support']]
    table = table.join(meta)
    # Never combine partial-cohort diversity with full-cohort coverage in a figure/table.
    if len(diversity_summary) and not diversity_missing:
        d = diversity_summary[(diversity_summary.gate=='posebusters')&(diversity_summary.budget==1000)&(diversity_summary.scenario=='requested')]
        for radius in config['cluster_radii']:
            subset = d[d.cluster_radius==radius].set_index('recipe')
            table[f'clusters_{radius_tag(radius)}'] = subset.clusters
            table[f'largest_fraction_{radius_tag(radius)}'] = subset.largest_fraction
    write_frame_atomic(table.reset_index(),root/'recipe_endpoint_table.csv')
    if 'clusters_r0p5' in table:
        # Descriptive two-endpoint frontier; not a significance test or a tuned score.
        scores = table[['coverage_0p5', 'clusters_r0p5']].to_numpy()
        frontier = []
        for score in scores:
            dominated = np.any(np.all(scores >= score, axis=1) & np.any(scores > score, axis=1))
            frontier.append(not dominated)
        write_frame_atomic(table.reset_index().assign(descriptive_pareto_frontier=frontier), root/'joint_frontier.csv')
    render(config,panel,aggregate,table,diversity_summary)
    completions_pending = [name for name,path in config['panels'][panel].get('completed_pools',{}).items()
                           if not (Path(config['output_dir'])/path).is_file()]
    status = ('analysis_complete_on_current_pools' if completions_pending else 'complete') if not diversity_missing else 'coverage_complete_diversity_partial'
    recipe_count = len(recipe_grid(config['sources']))
    atomic_json_dump(root/'run_manifest.json',{'status':status,'recipes':recipe_count,'eligible_molecules':len(pivot),
        'source_completions_pending':completions_pending,
        'diversity_molecules_missing':diversity_missing,'config_sha256':config['_sha256'],
        'generation_pools': {name: {'path': str(pool_path(config,panel,name,source)),
                                  'sha256': sha256(pool_path(config,panel,name,source))}
                             for name in config['panels'][panel]['methods']},
        'cli_source_sha256':sha256(Path(__file__)),
        'statistics_source_sha256':sha256(Path(__file__).resolve().parents[1]/'analysis/dataset_recipes.py'),
        'analysis_package_versions': {name:importlib.metadata.version(name) for name in ['rdkit','posebusters','numpy','scipy']},
        'summary_sha256':sha256(root/'recipe_endpoint_table.csv')})
    report = [f'# {panel}: dataset recipe exploration','',f'Status: {status}. All comparisons are exploratory.',
              f'Parent-pool completions still pending: {", ".join(completions_pending) or "none"}.',
              '', 'Primary: PoseBusters-passing hits; failures consume budget; all frozen eligible molecules retained.',
              'The realized candidate counts and complete-budget molecule counts must accompany coverage rankings.',
              '', '```text', table.sort_values('coverage_0p5',ascending=False).to_string(), '```',
              '', f'Use paired_recipe_differences.csv for paired uncertainty. Intervals are unadjusted for searching {recipe_count} recipes.',
              'Use equal_availability rows in coverage_summary.csv to distinguish mixture effects from source shortfalls.',
              'Reference-group weighting and 0.75/1.0 Å tolerances are sensitivities, not replacement endpoints.']
    report_path = root/'FINDINGS.md'
    temporary = report_path.with_name(f'.{report_path.name}.{os.getpid()}.tmp')
    temporary.write_text('\n'.join(report)+'\n')
    os.replace(temporary, report_path)
    print(table.sort_values('coverage_0p5',ascending=False).to_string(),flush=True)


def render(config,panel,aggregate,table,diversity):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    root = Path(config['output_dir'])/panel/'figures'
    root.mkdir(exist_ok=True)
    fig,ax = plt.subplots(figsize=(10,8))
    ordered = table.sort_values('coverage_0p5')
    colors = ['#777777' if name.startswith('rdkit') else '#246a93' for name in ordered.index]
    ax.barh(ordered.index,ordered.coverage_0p5*100,color=colors)
    ax.errorbar(ordered.coverage_0p5*100,np.arange(len(ordered)),
                xerr=np.maximum(0,np.vstack([ordered.coverage_0p5-ordered.ci_low,ordered.ci_high-ordered.coverage_0p5])*100),
                fmt='none',color='black',capsize=2)
    ax.set_xlabel('Experimental pharmacophore coverage (%) · 0.5 Å · PB-valid hits')
    ax.set_title(f'{panel}: requested 1,000 candidates; source shortfalls reported in table')
    fig.tight_layout()
    fig.savefig(root/'recipe_coverage.png',dpi=180)
    fig.savefig(root/'recipe_coverage.pdf')
    plt.close(fig)
    if 'clusters_r0p5' in table:
        fig,axes = plt.subplots(1,2,figsize=(12,6))
        focus = {'qwen': 'Qwen', 'loqi': 'LoQI', 'flowr': 'FlowR',
                 'qwen75_loqi25': 'Qwen/LoQI 75/25', 'qwen50_loqi50': 'Qwen/LoQI 50/50',
                 'qwen25_loqi75': 'Qwen/LoQI 25/75',
                 'qwen50_loqi25_flowr25': 'Qwen/LoQI/FlowR 50/25/25',
                 'rdkit_raw': 'RDKit raw', 'rdkit_minimized': 'RDKit minimized'}
        focus = {name: label for name, label in focus.items() if name in table.index}
        palette = plt.get_cmap('tab10')
        for radius,ax in zip([.5,.75],axes):
            coverage_column = f'coverage_{"0p5" if radius==.5 else "0p75"}'
            other = table.loc[~table.index.isin(focus)]
            ax.scatter(other.clusters_r0p5,other[coverage_column]*100,s=25,color='#b8b8b8',alpha=.7,
                       label='Other mixture ratios',zorder=1)
            for number,(name,label) in enumerate(focus.items()):
                row = table.loc[name]
                marker = 'X' if name.startswith('rdkit') else '^' if 'flowr25' in name else 's' if '_' in name else 'o'
                ax.scatter(row.clusters_r0p5,row[coverage_column]*100,s=65,color=palette(number),
                           marker=marker,label=label,zorder=3)
            ax.set_xlabel('Mean PB-valid geometric clusters · 0.5 Å')
            ax.set_ylabel(f'Pharmacophore coverage (%) · {radius:g} Å')
            ax.margins(.1)
            ax.grid(alpha=.2)
        fig.suptitle(f'{panel}: recipe diversity and experimental feature coverage (exploratory)')
        handles,labels = axes[0].get_legend_handles_labels()
        fig.legend(handles,labels,loc='lower center',ncol=3,fontsize=9,frameon=False)
        fig.tight_layout(rect=[0,.18,1,.95])
        fig.savefig(root/'recipe_diversity_coverage.png',dpi=200)
        fig.savefig(root/'recipe_diversity_coverage.pdf')
        plt.close(fig)
    fig,axes = plt.subplots(1,2,figsize=(13,5))
    focus = aggregate[(aggregate.gate=='posebusters')&(aggregate.radius==.5)&(aggregate.weighting=='queries')]
    for scenario,ax in zip(['requested','equal_availability'],axes):
        for name,group in focus[focus.scenario==scenario].groupby('recipe'):
            group=group.sort_values('budget')
            ax.plot(group.budget,group.coverage*100,marker='.',label=name,linewidth=1)
        ax.set_title(scenario.replace('_',' '))
        ax.set_xlabel('Requested candidate allowance')
        ax.set_ylabel('Coverage (%) · 0.5 Å')
        ax.grid(alpha=.2)
    axes[1].legend(fontsize=6,bbox_to_anchor=(1.02,1),loc='upper left')
    fig.tight_layout()
    fig.savefig(root/'recipe_budget_curves.png',dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['freeze','prepare','diversity','summarize'])
    parser.add_argument('--config',type=Path,default=DEFAULT_CONFIG)
    parser.add_argument('--panel',choices=['plinder22','atc48'])
    parser.add_argument('--index',type=int)
    parser.add_argument('--workers',type=int)
    parser.add_argument('--sources',nargs='+',help='Prepare a subset while another generator is pending')
    parser.add_argument('--allow-partial',action='store_true')
    args=parser.parse_args()
    config=yaml.safe_load(args.config.read_text())
    config['_sha256']=sha256(args.config)
    if args.workers:
        config['workers']=args.workers
    if args.stage=='freeze':
        freeze(config,args.config)
        return
    protocol=Path(config['output_dir'])/'protocol.json'
    if not protocol.is_file() or json.loads(protocol.read_text())['config_sha256']!=config['_sha256']:
        raise ValueError('Freeze the unchanged protocol before analysis')
    if args.panel is None:
        parser.error('--panel required')
    if args.stage=='summarize':
        summarize(config,args.panel,args.allow_partial)
    else:
        if args.index is None:
            parser.error('--index required')
        if args.stage=='prepare':
            prepare(config,args.panel,args.index,args.sources)
        else:
            diversity(config,args.panel,args.index)


if __name__=='__main__':
    main()
