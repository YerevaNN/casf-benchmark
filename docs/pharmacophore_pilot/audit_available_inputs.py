"""Read-only inventory and feature-extraction smoke check; no query matching or optimization."""
from pathlib import Path
from collections import Counter
import hashlib
import importlib.util
import json
import pickle
import shutil

import numpy as np
import pandas as pd
from rdkit import Chem, RDConfig, rdBase
from rdkit.Chem import ChemicalFeatures

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
REF = Path('/mnt/weka/mbedrosian/data/plinder/cleaned/split_by_molecule_dedup_0.05_anchor_triplet')
PAPER = ROOT.parent / 'casf-benchmark-paper/jcim-overleaf'
POOLS = {
    'qwen_1p7b_v2_coordtoken_step41054': Path('/mnt/weka/mbedrosian/pharma_generation_analysis/v2_1p7b_coordtoken_step41054/gen_results/qwen_1p7b_v2_coordtoken_step41054_plinder/generation_results.pickle'),
    'loqi_druglike': Path('/mnt/weka/mbedrosian/codex_dir/druglike_off_the_shelf/gen_results/loqi_druglike/generation_results.pickle'),
    'rdkit_random_raw_druglike': Path('/mnt/weka/vtarasov/outputs/casf_benchmark_runs/druglike_baselines/gen_results/rdkit_random_raw_druglike/generation_results.pickle'),
    'rdkit_random_minimized_22_druglike': Path('/mnt/weka/vtarasov/outputs/casf_benchmark_runs/druglike_baselines/gen_results/rdkit_random_minimized_22_druglike/generation_results.pickle'),
}

def digest(path):
    with path.open('rb') as handle:
        h = hashlib.sha256()
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()

fdef = Path(RDConfig.RDDataDir) / 'BaseFeatures.fdef'
factory = ChemicalFeatures.BuildFeatureFactory(str(fdef))
families = {'Donor', 'Acceptor', 'Aromatic', 'Hydrophobe', 'LumpedHydrophobe', 'PosIonizable', 'NegIonizable'}
references, identity, hashes = [], {}, {}
for path in sorted(REF.glob('*.sdf')):
    hashes[str(path)] = digest(path)
    for record_index, mol in enumerate(Chem.SDMolSupplier(str(path), removeHs=False)):
        assert mol is not None
        before = mol.GetConformer().GetPositions().copy()
        smiles = Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)
        identity[smiles] = mol.GetProp('molecule_name')
        feats = [f for f in factory.GetFeaturesForMol(mol) if f.GetFamily() in families]
        centers = np.array([tuple(f.GetPos()) for f in feats])
        # Only an upper-bound eligibility screen: co-located feature labels are one center.
        centers = np.unique(np.round(centers, 6), axis=0)
        counts = Counter(f.GetFamily() for f in feats)
        np.testing.assert_array_equal(before, mol.GetConformer().GetPositions())
        references.append(dict(ccd=path.stem, name=identity[smiles], smiles=smiles,
            record_index=record_index, conf_index=mol.GetProp('conf_index'),
            original_reference_index=mol.GetProp('reference_index'),
            pdb_id=mol.GetProp('pdb_id'), system_id=mol.GetProp('system_id'),
            n_features=len(feats), n_distinct_centers=len(centers),
            center_rank=int(np.linalg.matrix_rank(centers-centers.mean(axis=0))) if len(centers) else 0,
            feature_counts=json.dumps(counts, sort_keys=True)))
refs = pd.DataFrame(references)
assert len(refs) == 1067 and len(identity) == 22
refs.to_csv(OUT/'reference-feature-inventory.csv', index=False)
paper = pd.read_csv(PAPER/'plinder22-data/per-molecule.csv')
pool_rows, input_files = [], {}
for method, path in POOLS.items():
    input_files[method] = dict(path=str(path), sha256=digest(path), size_bytes=path.stat().st_size)
    with path.open('rb') as handle:
        pool = pickle.load(handle)
    assert set(identity) <= set(pool)
    current = paper[paper.label.eq(method)].set_index('smiles')
    for smiles, name in identity.items():
        candidates = pool[smiles]
        readable = [m for m in candidates if m is not None and m.GetNumConformers()]
        assert len(readable) == current.loc[smiles, 'n_generated']
        m = Chem.Mol(readable[0])
        before = m.GetConformer().GetPositions().copy()
        feats = factory.GetFeaturesForMol(m)
        np.testing.assert_array_equal(before, m.GetConformer().GetPositions())
        pool_rows.append(dict(method=method, smiles=smiles, name=name,
            stored_entries=len(candidates), readable_conformers=len(readable),
            paper_pb_pass=int(current.loc[smiles, 'n_pb_pass']),
            first_conformer_features=len(feats),
            all_have_one_conformer=all(m.GetNumConformers()==1 for m in readable)))
    input_files[method]['extra_identities_excluded'] = len(set(pool)-set(identity))
    del pool
pd.DataFrame(pool_rows).to_csv(OUT/'pool-inventory.csv', index=False)
annotations = pd.read_csv(ROOT/'docs/plinder22_audit/reference_audit.csv', keep_default_na=False)
assert len(annotations) == 1067
cache = ROOT/'docs/plinder23_audit/cache'
required_pdbs = sorted(refs.pdb_id.unique())
missing_pdbs = [p for p in required_pdbs if not (cache/f'{p}.cif.gz').exists()]
modules={}
for name in ['scipy','gemmi','plip','prolif','openbabel','pmapper']:
    spec=importlib.util.find_spec(name)
    modules[name]=spec.origin if spec else None
result=dict(rdkit_version=rdBase.rdkitVersion,feature_definition=str(fdef),feature_definition_sha256=digest(fdef),
    reference_files=hashes, pools=input_files, modules=modules,
    executables={name:shutil.which(name) for name in ['pharmitserver','pharmit','pharmer','plip','obabel']},
    n_reference_molecules=22,n_references=1067,n_required_pdbs=len(required_pdbs),missing_cached_pdbs=missing_pdbs,
    n_unique_raw_sdf_paths=int(annotations.raw_sdf_paths.nunique()),
    all_raw_sdf_paths_exist=all(Path(p).is_file() for p in annotations.raw_sdf_paths.unique()),
    n_single_raw_matches=int(annotations.raw_match_count.eq(1).sum()),
    n_covalent_flags=int(annotations.covalent_flag.eq(True).sum()),
    n_reference_geometry_failures=int(annotations.geometry_pass.eq(False).sum()),
    all_feature_centers_at_least_four=int((refs.n_distinct_centers>=4).sum()),
    feature_smoke_check='All 1,067 references and the first conformer of each of 88 method-molecule pools; coordinates unchanged.',
    limits=['All-feature counts are not interaction-selected query counts.',
            'Cached PDB availability does not verify assembly-copy or ligand atom mapping.',
            'PB counts are aggregate records, not recovered candidate masks.',
            'No conformer matching, PB rerun, minimization or new generation performed.'])
(OUT/'input-audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['rdkit_version','n_references','n_required_pdbs','missing_cached_pdbs','all_feature_centers_at_least_four']}))
print(pd.DataFrame(pool_rows).groupby('method').readable_conformers.agg(['sum','min','max']).to_string())
print(refs.groupby('name').n_distinct_centers.min().sort_values().to_string())
