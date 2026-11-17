"""PLINDER pharmacophore-query extraction, matching, and budget statistics.

The module treats cached coordinates as immutable. Pharmit supplies feature
definitions and interaction-query recognition. Matching exhaustively enumerates
typed feature correspondences before a proper rigid Kabsch check, avoiding
false negatives from Pharmit's indexed triplet search.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import pickle
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stable_smiles_id(smiles: str) -> str:
    return hashlib.sha256(smiles.encode()).hexdigest()[:16]


def stable_query_id(row: Any) -> str:
    raw = f"{row.ccd}|{int(row.conf_index)}|{int(row.reference_index)}|{row.system_id}"
    return f"q_{hashlib.sha256(raw.encode()).hexdigest()[:20]}"


def atomic_json_dump(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(temporary, path)


def load_pickle(path: Path) -> Any:
    with path.open("rb") as handle:
        return pickle.load(handle)


def canonical_isomeric_smiles(mol: Any) -> str:
    from rdkit import Chem

    return Chem.MolToSmiles(Chem.RemoveHs(Chem.Mol(mol)), isomericSmiles=True)


def molecule_coordinate_hash(mol: Any) -> str:
    """Hash topology and exact double-precision coordinates without properties."""
    from rdkit import Chem

    flags = Chem.PropertyPickleOptions.NoProps | Chem.PropertyPickleOptions.CoordsAsDouble
    return hashlib.sha256(mol.ToBinary(flags)).hexdigest()


def proper_kabsch(source: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Return the proper rotation/translation mapping source onto target."""
    source = np.asarray(source, dtype=float)
    target = np.asarray(target, dtype=float)
    if source.shape != target.shape or source.ndim != 2 or source.shape[1] != 3:
        raise ValueError("Kabsch inputs must have identical (n, 3) shapes")
    source_center = source.mean(axis=0)
    target_center = target.mean(axis=0)
    covariance = (source - source_center).T @ (target - target_center)
    left, _, right_t = np.linalg.svd(covariance)
    rotation = left @ right_t
    if np.linalg.det(rotation) < 0:
        left[:, -1] *= -1
        rotation = left @ right_t
    translation = target_center - source_center @ rotation
    fitted = source @ rotation + translation
    rmsd = float(np.sqrt(np.mean(np.sum((fitted - target) ** 2, axis=1))))
    if np.linalg.det(rotation) < 1.0 - 1e-8:
        raise RuntimeError("Kabsch produced an improper rotation")
    return rotation, translation, rmsd


def explicit_atom_mapping(cleaned: Any, native: Any, max_matches: int = 10000) -> dict[str, Any]:
    """Find the stereo-aware atom map and proper rigid transform with lowest RMSD."""
    from rdkit import Chem

    clean = Chem.RemoveHs(Chem.Mol(cleaned))
    raw = Chem.RemoveHs(Chem.Mol(native))
    if clean.GetNumAtoms() != raw.GetNumAtoms():
        raise ValueError("Cleaned/native heavy-atom counts differ")
    matches = raw.GetSubstructMatches(
        clean,
        useChirality=True,
        uniquify=False,
        maxMatches=max_matches,
    )
    if not matches:
        raise ValueError("No stereo-aware cleaned/native atom correspondence")
    clean_xyz = np.asarray(clean.GetConformer().GetPositions(), dtype=float)
    raw_xyz = np.asarray(raw.GetConformer().GetPositions(), dtype=float)
    best: tuple[float, tuple[int, ...], np.ndarray, np.ndarray] | None = None
    for match in matches:
        rotation, translation, rmsd = proper_kabsch(clean_xyz, raw_xyz[list(match)])
        if best is None or rmsd < best[0]:
            best = (rmsd, match, rotation, translation)
    assert best is not None
    rmsd, match, rotation, translation = best
    return {
        "cleaned_to_native_atom_map": list(enumerate(map(int, match))),
        "rotation": rotation.tolist(),
        "translation": translation.tolist(),
        "rmsd": rmsd,
        "n_candidate_maps": len(matches),
    }


def _mmcif_rows(block: Any, category: str) -> list[dict[str, str]]:
    data = block.get_mmcif_category(category)
    return [dict(zip(data, values)) for values in zip(*data.values())] if data else []


def _original_subchain(subchain: str) -> str:
    match = re.match(r"^orig_(.*)_\d+$", subchain)
    if match:
        return match.group(1)
    if subchain.startswith("orig_") and subchain.endswith("_"):
        return subchain[5:-1]
    return subchain


def _atom_xyz(atom: Any) -> np.ndarray:
    return np.asarray([atom.pos.x, atom.pos.y, atom.pos.z], dtype=float)


def _best_assembly_ligand_copy(model: Any, ligand_asym: str, native_mol: Any) -> dict[str, Any]:
    """Resolve the assembly ligand copy and alternate locations by native coordinates."""
    from rdkit import Chem
    from scipy.optimize import linear_sum_assignment

    raw = Chem.RemoveHs(Chem.Mol(native_mol))
    raw_xyz = np.asarray(raw.GetConformer().GetPositions(), dtype=float)
    raw_elements = [atom.GetSymbol().upper() for atom in raw.GetAtoms()]
    candidates = []
    for chain_index, chain in enumerate(model):
        for residue_index, residue in enumerate(chain):
            if not residue.subchain.startswith(f"orig_{ligand_asym}_"):
                continue
            atoms = [atom for atom in residue if not atom.element.is_hydrogen]
            if len(atoms) < len(raw_elements):
                continue
            coords = np.asarray([_atom_xyz(atom) for atom in atoms])
            elements = [atom.element.name.upper() for atom in atoms]
            cost = np.linalg.norm(raw_xyz[:, None, :] - coords[None, :, :], axis=2)
            for i, element in enumerate(raw_elements):
                cost[i, np.asarray(elements) != element] = 1e6
            raw_indices, assembly_indices = linear_sum_assignment(cost)
            if len(raw_indices) != len(raw_elements):
                continue
            selected_costs = cost[raw_indices, assembly_indices]
            if np.any(selected_costs >= 1e5):
                continue
            rmsd = float(np.sqrt(np.mean(selected_costs**2)))
            candidates.append(
                {
                    "chain_index": chain_index,
                    "residue_index": residue_index,
                    "chain_name": chain.name,
                    "residue_name": residue.name,
                    "residue_seqid": str(residue.seqid),
                    "subchain": residue.subchain,
                    "rmsd_without_superposition": rmsd,
                    "max_atom_distance": float(selected_costs.max()),
                    "raw_to_assembly_atom": [
                        {
                            "raw_atom_index": int(raw_index),
                            "assembly_atom_index": int(assembly_index),
                            "assembly_atom_name": atoms[assembly_index].name,
                            "altloc": str(atoms[assembly_index].altloc),
                            "occupancy": float(atoms[assembly_index].occ),
                        }
                        for raw_index, assembly_index in zip(
                            raw_indices, assembly_indices, strict=True
                        )
                    ],
                }
            )
    if not candidates:
        raise ValueError(f"No assembly ligand candidates for asym {ligand_asym}")
    candidates.sort(key=lambda item: item["rmsd_without_superposition"])
    if len(candidates) > 1 and math.isclose(
        candidates[0]["rmsd_without_superposition"],
        candidates[1]["rmsd_without_superposition"],
        abs_tol=1e-6,
    ):
        raise ValueError("Assembly ligand copy is coordinate-ambiguous")
    result = candidates[0]
    result["n_assembly_candidates"] = len(candidates)
    return result


def _choose_altloc_atoms(residue: Any) -> list[Any]:
    grouped: dict[str, list[Any]] = {}
    for atom in residue:
        if atom.element.is_hydrogen:
            continue
        grouped.setdefault(atom.name, []).append(atom)
    chosen = []
    for name in grouped:
        options = grouped[name]
        options.sort(
            key=lambda atom: (
                -float(atom.occ),
                str(atom.altloc) not in {"\x00", " ", "."},
                str(atom.altloc),
            )
        )
        chosen.append(options[0])
    return chosen


def reconstruct_native_receptor(
    *,
    mmcif_path: Path,
    assembly_name: str,
    ligand_asym: str,
    native_mol: Any,
    receptor_path: Path,
    contact_cutoff: float,
) -> dict[str, Any]:
    """Select the exact assembly copy and export every contacting protein chain."""
    import gemmi
    from rdkit import Chem
    from scipy.spatial.distance import cdist

    block = gemmi.cif.read(str(mmcif_path)).sole_block()
    structure = gemmi.make_structure_from_block(block)
    entity_poly = {
        row["entity_id"]: row.get("type", "")
        for row in _mmcif_rows(block, "_entity_poly.")
    }
    asym_entity = {
        row["id"]: row["entity_id"]
        for row in _mmcif_rows(block, "_struct_asym.")
    }
    for chain in structure[0]:
        chain.name = f"orig_{chain.name}_"
        for residue in chain:
            residue.subchain = f"orig_{residue.subchain}_"
    for assembly in structure.assemblies:
        for generator in assembly.generators:
            generator.chains = [f"orig_{name}_" for name in generator.chains]
            generator.subchains = [f"orig_{name}_" for name in generator.subchains]
    assembly = next(
        (item for item in structure.assemblies if item.name == str(assembly_name)),
        None,
    )
    if assembly is None:
        raise ValueError(f"Biological assembly {assembly_name} not found")
    model = gemmi.make_assembly(
        assembly,
        structure[0],
        gemmi.HowToNameCopiedChain.AddNumber,
    )
    ligand_copy = _best_assembly_ligand_copy(model, ligand_asym, native_mol)
    raw = Chem.RemoveHs(Chem.Mol(native_mol))
    ligand_xyz = np.asarray(raw.GetConformer().GetPositions(), dtype=float)

    protein_chain_indices: list[int] = []
    chain_records = []
    context = {"nucleic_acids": [], "metals": [], "waters": [], "cofactors": []}
    metal_symbols = {
        "LI", "NA", "K", "RB", "CS", "MG", "CA", "SR", "BA", "MN", "FE",
        "CO", "NI", "CU", "ZN", "CD", "HG",
    }
    for chain_index, chain in enumerate(model):
        chain_atoms = []
        polymer_types = set()
        for residue in chain:
            original_asym = _original_subchain(residue.subchain)
            polymer_type = entity_poly.get(asym_entity.get(original_asym, ""), "")
            if polymer_type:
                polymer_types.add(polymer_type)
            chain_atoms.extend(_atom_xyz(atom) for atom in _choose_altloc_atoms(residue))
        if not chain_atoms:
            continue
        minimum = float(cdist(ligand_xyz, np.asarray(chain_atoms)).min())
        is_protein = any(value.startswith("polypeptide") for value in polymer_types)
        if is_protein and minimum <= contact_cutoff:
            protein_chain_indices.append(chain_index)
            chain_records.append(
                {
                    "assembly_chain": chain.name,
                    "original_subchains": sorted(
                        {_original_subchain(residue.subchain) for residue in chain}
                    ),
                    "polymer_types": sorted(polymer_types),
                    "minimum_ligand_distance": minimum,
                }
            )
        elif minimum <= contact_cutoff and any(
            "ribonucleotide" in value or "deoxyribonucleotide" in value
            for value in polymer_types
        ):
            context["nucleic_acids"].append(chain.name)

    # Record nearby non-protein context, but do not give it to the extractor.
    for chain in model:
        for residue in chain:
            original_asym = _original_subchain(residue.subchain)
            polymer_type = entity_poly.get(asym_entity.get(original_asym, ""), "")
            if polymer_type:
                continue
            atoms = _choose_altloc_atoms(residue)
            if not atoms:
                continue
            minimum = float(cdist(ligand_xyz, np.asarray([_atom_xyz(a) for a in atoms])).min())
            if minimum > contact_cutoff:
                continue
            if residue.name in {"HOH", "DOD"}:
                context["waters"].append(residue.name)
            elif len(atoms) == 1 and atoms[0].element.name.upper() in metal_symbols:
                context["metals"].append(atoms[0].element.name.upper())
            elif not residue.subchain.startswith(f"orig_{ligand_asym}_"):
                context["cofactors"].append(residue.name)
    if not protein_chain_indices:
        raise ValueError("No directly contacting protein chain")

    receptor_path.parent.mkdir(parents=True, exist_ok=True)
    chain_alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    if len(protein_chain_indices) > len(chain_alphabet):
        raise ValueError("Too many contacting protein chains for PDB export")
    lines = []
    serial = 1
    for output_chain_index, model_chain_index in enumerate(protein_chain_indices):
        chain = model[model_chain_index]
        output_chain = chain_alphabet[output_chain_index]
        for residue in chain:
            original_asym = _original_subchain(residue.subchain)
            polymer_type = entity_poly.get(asym_entity.get(original_asym, ""), "")
            if not polymer_type.startswith("polypeptide"):
                continue
            residue_number = int(residue.seqid.num)
            insertion_code = str(residue.seqid.icode).strip()[:1]
            for atom in _choose_altloc_atoms(residue):
                x, y, z = _atom_xyz(atom)
                atom_name = atom.name[:4]
                element = atom.element.name.upper()[:2]
                lines.append(
                    f"ATOM  {serial:5d} {atom_name:>4s} {residue.name[:3]:>3s} "
                    f"{output_chain}{residue_number:4d}{insertion_code:1s}   "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}{float(atom.occ):6.2f}"
                    f"{float(atom.b_iso):6.2f}          {element:>2s}\n"
                )
                serial += 1
        lines.append("TER\n")
    lines.append("END\n")
    receptor_path.write_text("".join(lines), encoding="ascii")
    ligand_copy["protein_chains"] = chain_records
    ligand_copy["context"] = {key: sorted(values) for key, values in context.items()}
    ligand_copy["receptor_sha256"] = sha256(receptor_path)
    ligand_copy["receptor_path"] = str(receptor_path)
    return ligand_copy


def distinct_centers(points: Sequence[dict[str, Any]], separation: float) -> list[np.ndarray]:
    centers: list[np.ndarray] = []
    for point in points:
        center = np.asarray([point[axis] for axis in "xyz"], dtype=float)
        if all(np.linalg.norm(center - existing) >= separation for existing in centers):
            centers.append(center)
    return centers


def maximum_triangle_area(centers: Sequence[np.ndarray]) -> float:
    maximum = 0.0
    for first, second, third in combinations(centers, 3):
        area = 0.5 * float(np.linalg.norm(np.cross(second - first, third - first)))
        maximum = max(maximum, area)
    return maximum


@dataclass(frozen=True)
class PharmitBackend:
    executable: Path
    timeout_seconds: int = 3600
    skip_triplet_fingerprints: bool = True
    skip_geometric_pruning: bool = True

    def _run(self, arguments: Sequence[str], *, cwd: Path | None = None) -> subprocess.CompletedProcess:
        command = [str(self.executable), *map(str, arguments)]
        result = subprocess.run(
            command,
            cwd=cwd,
            check=False,
            text=True,
            capture_output=True,
            timeout=self.timeout_seconds,
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"Pharmit command failed with exit {result.returncode}: "
                f"stdout={result.stdout[-4000:]!r}; stderr={result.stderr[-4000:]!r}"
            )
        return result

    def extract_interaction_query(
        self,
        ligand_path: Path,
        receptor_path: Path,
        output_path: Path,
        radius: float,
    ) -> dict[str, Any]:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        self._run(
            [
                "--cmd", "pharma", "--in", ligand_path, "--receptor", receptor_path,
                "--out", output_path,
            ]
        )
        query = json.loads(output_path.read_text())
        # The CLI writes JSON null when no interaction features are found.
        # Treat that as an empty feature set so reference eligibility can record
        # a scientific exclusion instead of a backend-processing failure.
        points = [point for point in (query.get("points") or []) if point.get("enabled", True)]
        for point in points:
            point["radius"] = float(radius)
            point["vector_on"] = False
            point["enabled"] = True
            point["requirement"] = "required"
        query["points"] = points
        atomic_json_dump(output_path, query)
        return query

    def create_database(self, sdf_path: Path, database_path: Path) -> str:
        if database_path.exists():
            shutil.rmtree(database_path)
        database_path.parent.mkdir(parents=True, exist_ok=True)
        result = self._run(
            [
                "--cmd", "dbcreate", "--dbdir", database_path, "--in", sdf_path,
                "--reduceconfs", "0",
            ]
        )
        return result.stdout + result.stderr

    def recognize_features(self, sdf_path: Path, output_path: Path) -> list[list[dict[str, Any]]]:
        """Recognize every conformer's features using Pharmit's definitions."""
        output_path.parent.mkdir(parents=True, exist_ok=True)
        self._run(["--cmd", "pharma", "--q", "--in", sdf_path, "--out", output_path])
        serialized = output_path.read_text(encoding="utf-8")
        decoder = json.JSONDecoder()
        position = 0
        feature_sets: list[list[dict[str, Any]]] = []
        while True:
            while position < len(serialized) and serialized[position].isspace():
                position += 1
            if position == len(serialized):
                break
            payload, position = decoder.raw_decode(serialized, position)
            feature_sets.append(
                [point for point in (payload.get("points") or []) if point.get("enabled", True)]
            )
        return feature_sets

    def search_many(
        self,
        *,
        database_path: Path,
        query_paths: Sequence[Path],
        output_paths: Sequence[Path],
        max_hits: int,
        workers: int,
    ) -> str:
        if len(query_paths) != len(output_paths):
            raise ValueError("Every query needs an output path")
        arguments: list[Any] = [
            "--cmd", "dbsearch", "--dbdir", database_path,
            "--max-hits", max_hits, "--max-orient", 1, "--reduceconfs", 0,
            "--nthreads", workers, "--print=false",
        ]
        if self.skip_triplet_fingerprints:
            # Pharmit's compressed triplet fingerprint is an acceleration
            # heuristic and can reject exact geometric matches, particularly
            # when query radii change.  Disable it for complete membership.
            arguments.append("--skip-fingers")
        if self.skip_geometric_pruning:
            # Pharmit's higher-order triplet bounds are not conservative for
            # all repeated-feature/radius combinations. Independent triplet
            # ranges and the exact final rigid correspondence check remain on.
            arguments.append("--skip-geometric-pruning")
        for path in query_paths:
            arguments.extend(["--in", path])
        for path in output_paths:
            arguments.extend(["--out", path])
        result = self._run(arguments)
        return result.stdout + result.stderr


def transform_query(query: dict[str, Any], rotation: np.ndarray, translation: np.ndarray) -> dict[str, Any]:
    transformed = json.loads(json.dumps(query))
    rotation = np.asarray(rotation, dtype=float)
    translation = np.asarray(translation, dtype=float)
    if not math.isclose(float(np.linalg.det(rotation)), 1.0, abs_tol=1e-7):
        raise ValueError("Query transform must be a proper rotation")
    for point in transformed["points"]:
        center = np.asarray([point[axis] for axis in "xyz"], dtype=float)
        moved = center @ rotation + translation
        point.update(dict(zip("xyz", map(float, moved), strict=True)))
    return transformed


def export_candidate_library(
    *,
    method: str,
    smiles: str,
    mols: Sequence[Any],
    sdf_path: Path,
) -> list[dict[str, Any]]:
    from rdkit import Chem

    sdf_path.parent.mkdir(parents=True, exist_ok=True)
    original_hashes = [molecule_coordinate_hash(mol) for mol in mols]
    records = []
    with Chem.SDWriter(str(sdf_path)) as writer:
        for index, mol in enumerate(mols):
            candidate_id = f"{method}__{stable_smiles_id(smiles)}__{index:04d}"
            record = {
                "method": method,
                "smiles": smiles,
                "smiles_id": stable_smiles_id(smiles),
                "candidate_index": index,
                "candidate_id": candidate_id,
                "coordinate_sha256": original_hashes[index],
                "exported": False,
                "export_error": "",
            }
            try:
                candidate = Chem.Mol(mol)
                if candidate.GetNumConformers() == 0:
                    raise ValueError("no_conformer")
                candidate.SetProp("_Name", candidate_id)
                candidate.SetProp("candidate_id", candidate_id)
                writer.write(candidate)
                record["exported"] = True
            except Exception as exc:  # keep the supplied denominator and a false match
                record["export_error"] = f"{type(exc).__name__}: {exc}"
            records.append(record)
    after_hashes = [molecule_coordinate_hash(mol) for mol in mols]
    if after_hashes != original_hashes:
        raise RuntimeError("Candidate export mutated cached coordinates")
    return records


def read_match_titles(path: Path) -> list[str]:
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    if not path.is_file() or path.stat().st_size == 0:
        return []
    titles = []
    supplier = Chem.SDMolSupplier(str(path), removeHs=False)
    for mol in supplier:
        if mol is None:
            raise ValueError(f"Unreadable Pharmit result in {path}")
        titles.append(mol.GetProp("_Name"))
    if len(titles) != len(set(titles)):
        raise ValueError(f"Duplicate candidate membership in {path}")
    return titles


def export_normalized_feature_library(
    *,
    smiles: str,
    mols: Sequence[Any],
    candidate_ids: Sequence[str],
    sdf_path: Path,
) -> dict[str, Any]:
    """Export molecules in a common heavy-atom order for feature recognition."""
    from rdkit import Chem

    if len(mols) != len(candidate_ids):
        raise ValueError("Every feature-library molecule needs a candidate ID")
    reference = Chem.MolFromSmiles(smiles)
    if reference is None:
        raise ValueError(f"Cannot parse feature-library SMILES: {smiles}")
    sdf_path.parent.mkdir(parents=True, exist_ok=True)
    with Chem.SDWriter(str(sdf_path)) as writer:
        for mol, candidate_id in zip(mols, candidate_ids, strict=True):
            candidate = Chem.RemoveHs(Chem.Mol(mol))
            match = candidate.GetSubstructMatch(reference, useChirality=True)
            if len(match) != candidate.GetNumAtoms():
                # Top-up acceptance already froze molecular identity. This map
                # only establishes a stable atom permutation for feature output;
                # it must not become a second stereochemical acceptance gate.
                match = candidate.GetSubstructMatch(reference, useChirality=False)
            if len(match) != candidate.GetNumAtoms() or len(match) != reference.GetNumAtoms():
                raise ValueError(
                    f"Cannot put {candidate_id} into the coordinate-free reference atom order"
                )
            normalized = Chem.RenumberAtoms(candidate, list(match))
            normalized.SetProp("_Name", candidate_id)
            writer.write(normalized)
    return {
        "record_count": len(candidate_ids),
        "reference_atom_count": reference.GetNumAtoms(),
        "atom_order": "coordinate-free isomeric SMILES heavy-atom order",
        "sdf_sha256": sha256(sdf_path),
    }


def exact_feature_correspondence_transform(
    query_coordinates: np.ndarray,
    query_types: Sequence[str],
    candidate_coordinates: np.ndarray,
    candidate_types: Sequence[str],
    radius: float,
    *,
    numerical_tolerance: float = 2e-4,
) -> dict[str, Any] | None:
    """Return one exact typed correspondence and its candidate-to-query transform."""
    query = np.asarray(query_coordinates, dtype=float)
    candidate = np.asarray(candidate_coordinates, dtype=float)
    if query.shape != (len(query_types), 3):
        raise ValueError("Query feature coordinates/types disagree")
    if candidate.shape != (len(candidate_types), 3):
        raise ValueError("Candidate feature coordinates/types disagree")
    domains = [
        [index for index, feature_type in enumerate(candidate_types) if feature_type == qtype]
        for qtype in query_types
    ]
    if any(not domain for domain in domains):
        return None
    query_distances = np.linalg.norm(query[:, None, :] - query[None, :, :], axis=2)
    candidate_distances = np.linalg.norm(
        candidate[:, None, :] - candidate[None, :, :], axis=2
    )
    assigned: dict[int, int] = {}
    used: set[int] = set()
    pair_limit = 2.0 * float(radius) + 2.0 * numerical_tolerance

    def search() -> dict[str, Any] | None:
        if len(assigned) == len(query):
            query_indices = sorted(assigned)
            candidate_indices = [assigned[index] for index in query_indices]
            selected = candidate[candidate_indices]
            rotation, translation, _ = proper_kabsch(selected, query[query_indices])
            residuals = np.linalg.norm(
                selected @ rotation + translation - query[query_indices], axis=1
            )
            if np.all(residuals <= float(radius) + numerical_tolerance):
                return {
                    "query_indices": query_indices,
                    "candidate_indices": candidate_indices,
                    "rotation": rotation,
                    "translation": translation,
                    "residuals": residuals,
                }
            return None

        next_query = -1
        next_values: list[int] | None = None
        for query_index, domain in enumerate(domains):
            if query_index in assigned:
                continue
            values = [
                candidate_index
                for candidate_index in domain
                if candidate_index not in used
                and all(
                    abs(
                        query_distances[query_index, assigned_query]
                        - candidate_distances[candidate_index, assigned_candidate]
                    )
                    <= pair_limit
                    for assigned_query, assigned_candidate in assigned.items()
                )
            ]
            if not values:
                return None
            if next_values is None or len(values) < len(next_values):
                next_query = query_index
                next_values = values
        assert next_values is not None
        for candidate_index in next_values:
            assigned[next_query] = candidate_index
            used.add(candidate_index)
            result = search()
            if result is not None:
                return result
            used.remove(candidate_index)
            del assigned[next_query]
        return None

    return search()


def exact_feature_correspondence_match(
    query_coordinates: np.ndarray,
    query_types: Sequence[str],
    candidate_coordinates: np.ndarray,
    candidate_types: Sequence[str],
    radius: float,
    *,
    numerical_tolerance: float = 2e-4,
) -> bool:
    """Exhaustively test injective typed correspondences under proper rotation."""
    return (
        exact_feature_correspondence_transform(
            query_coordinates,
            query_types,
            candidate_coordinates,
            candidate_types,
            radius,
            numerical_tolerance=numerical_tolerance,
        )
        is not None
    )


def exact_feature_membership(
    query_points: Sequence[dict[str, Any]],
    candidate_coordinates: Sequence[np.ndarray] | np.ndarray,
    candidate_types: Sequence[str] | Sequence[Sequence[str]],
    radius: float,
) -> np.ndarray:
    """Exhaustively match one query against all candidate conformers."""
    query_coordinates = np.asarray(
        [[point[axis] for axis in "xyz"] for point in query_points], dtype=float
    )
    query_types = [str(point["name"]) for point in query_points]
    coordinate_sets = list(candidate_coordinates)
    if candidate_types and isinstance(candidate_types[0], str):
        type_sets: Sequence[Sequence[str]] = [candidate_types] * len(coordinate_sets)  # type: ignore[list-item]
    else:
        type_sets = candidate_types  # type: ignore[assignment]
    if len(coordinate_sets) != len(type_sets):
        raise ValueError("Candidate feature coordinates/types have different record counts")
    return np.asarray(
        [
            exact_feature_correspondence_match(
                query_coordinates,
                query_types,
                np.asarray(conformer_coordinates, dtype=float),
                conformer_types,
                radius,
            )
            for conformer_coordinates, conformer_types in zip(
                coordinate_sets, type_sets, strict=True
            )
        ],
        dtype=bool,
    )


def hypergeometric_no_hit_probability(population: int, hits: int, draws: int) -> float:
    """Exact uniform-without-replacement no-hit probability."""
    if population < 0 or not 0 <= hits <= population:
        raise ValueError("Invalid population/hit count")
    if not 0 <= draws <= population:
        raise ValueError("Draw count exceeds population")
    misses = population - hits
    if draws > misses:
        return 0.0
    if draws == 0 or hits == 0:
        return 1.0
    log_probability = (
        math.lgamma(misses + 1)
        - math.lgamma(draws + 1)
        - math.lgamma(misses - draws + 1)
        - math.lgamma(population + 1)
        + math.lgamma(draws + 1)
        + math.lgamma(population - draws + 1)
    )
    return float(math.exp(log_probability))


def hit_probability(population: int, hits: int, budget: int) -> float:
    draws = min(population, budget)
    return 1.0 - hypergeometric_no_hit_probability(population, hits, draws)


def mixture_hit_probability(
    qwen_population: int,
    qwen_hits: int,
    loqi_population: int,
    loqi_hits: int,
    total_budget: int,
) -> float:
    if total_budget % 2:
        raise ValueError("Equal-allocation mixture requires an even budget")
    each = total_budget // 2
    if each > qwen_population or each > loqi_population:
        raise ValueError("Mixture budget exceeds a source pool")
    no_qwen = hypergeometric_no_hit_probability(qwen_population, qwen_hits, each)
    no_loqi = hypergeometric_no_hit_probability(loqi_population, loqi_hits, each)
    return 1.0 - no_qwen * no_loqi


def paired_bootstrap_difference(
    left: np.ndarray,
    right: np.ndarray,
    *,
    replicates: int,
    seed: int,
) -> dict[str, float]:
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    if left.shape != right.shape or left.ndim != 1 or len(left) == 0:
        raise ValueError("Paired bootstrap inputs must be nonempty equal vectors")
    difference = left - right
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, len(difference), size=(replicates, len(difference)))
    samples = difference[indices].mean(axis=1)
    return {
        "difference": float(difference.mean()),
        "ci_low": float(np.quantile(samples, 0.025)),
        "ci_high": float(np.quantile(samples, 0.975)),
        "bootstrap_replicates": int(replicates),
        "bootstrap_seed": int(seed),
    }


def backend_smoke_test(backend: PharmitBackend, output_dir: Path) -> dict[str, Any]:
    """Run late-sentinel, result-limit, and rigid-invariance checks."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    output_dir.mkdir(parents=True, exist_ok=True)
    molecule = Chem.AddHs(Chem.MolFromSmiles("CN1CCC[C@H]1c1cccnc1"))
    conformer_ids = list(AllChem.EmbedMultipleConfs(molecule, numConfs=11, randomSeed=7123))
    if len(conformer_ids) != 11:
        raise RuntimeError("Could not generate smoke-test conformers")
    query_mol = Chem.Mol(molecule)
    query_mol.RemoveAllConformers()
    query_mol.AddConformer(Chem.Conformer(molecule.GetConformer(conformer_ids[0])))
    query_mol.SetProp("_Name", "SMOKE_QUERY")
    query_ligand = output_dir / "query_ligand.sdf"
    with Chem.SDWriter(str(query_ligand)) as writer:
        writer.write(query_mol)
    library = output_dir / "library.sdf"
    with Chem.SDWriter(str(library)) as writer:
        for index, conformer_id in enumerate(conformer_ids[1:]):
            candidate = Chem.Mol(molecule)
            candidate.RemoveAllConformers()
            candidate.AddConformer(Chem.Conformer(molecule.GetConformer(conformer_id)))
            candidate.SetProp("_Name", f"SMOKE_{index:03d}")
            writer.write(candidate)
        sentinel = Chem.Mol(query_mol)
        sentinel.SetProp("_Name", "SMOKE_LATE_SENTINEL")
        writer.write(sentinel)
    query_path = output_dir / "query.json"
    backend._run(["--cmd", "pharma", "--in", query_ligand, "--out", query_path])
    query = json.loads(query_path.read_text())
    for point in query["points"]:
        point["radius"] = 0.1
        point["vector_on"] = False
    atomic_json_dump(query_path, query)
    database = output_dir / "database"
    create_log = backend.create_database(library, database)
    matches = output_dir / "matches.sdf"
    search_log = backend.search_many(
        database_path=database,
        query_paths=[query_path],
        output_paths=[matches],
        max_hits=1000,
        workers=2,
    )
    titles = read_match_titles(matches)
    if "SMOKE_LATE_SENTINEL" not in titles:
        raise RuntimeError("Pharmit did not recover the late exact-match sentinel")

    angle = 0.731
    rotation = np.asarray(
        [
            [math.cos(angle), -math.sin(angle), 0.0],
            [math.sin(angle), math.cos(angle), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    transformed = transform_query(query, rotation, np.asarray([7.0, -3.0, 2.0]))
    transformed_path = output_dir / "query_transformed.json"
    atomic_json_dump(transformed_path, transformed)
    transformed_matches = output_dir / "matches_transformed.sdf"
    backend.search_many(
        database_path=database,
        query_paths=[transformed_path],
        output_paths=[transformed_matches],
        max_hits=1000,
        workers=2,
    )
    transformed_titles = read_match_titles(transformed_matches)
    if set(transformed_titles) != set(titles):
        raise RuntimeError("Rigid query transform changed Pharmit membership")
    report = {
        "status": "passed",
        "late_sentinel_recovered": True,
        "candidate_titles_preserved": True,
        "result_limit": 1000,
        "library_size": 11,
        "matches": titles,
        "rigid_transform_membership_invariant": True,
        "proper_rotation_determinant": float(np.linalg.det(rotation)),
        "create_log": create_log,
        "search_log": search_log,
        "backend_executable": str(backend.executable),
        "backend_sha256": sha256(backend.executable),
        "triplet_fingerprint_pruning_disabled": backend.skip_triplet_fingerprints,
        "higher_order_triplet_pruning_disabled": backend.skip_geometric_pruning,
    }
    atomic_json_dump(output_dir / "smoke_report.json", report)
    return report

