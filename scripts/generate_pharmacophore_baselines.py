#!/usr/bin/env python3
"""Generate coordinate-free RDKit raw and force-field baseline pools."""

from __future__ import annotations

import argparse
import json
import os
import pickle
import tempfile
from collections import OrderedDict
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem


def load_smiles(path: Path) -> list[str]:
    with path.open("rb") as handle:
        payload = pickle.load(handle)
    if not isinstance(payload, dict) or not payload:
        raise TypeError("Input must be a nonempty SMILES-keyed mapping")
    return list(payload)


def atomic_pickle(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="wb", dir=path.parent, prefix=f".{path.name}.", delete=False
    ) as handle:
        temporary = Path(handle.name)
        pickle.dump(value, handle, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(temporary, path)


def single_conformer(molecule: Chem.Mol, conformer_id: int) -> Chem.Mol:
    result = Chem.Mol(molecule)
    conformer = Chem.Conformer(molecule.GetConformer(conformer_id))
    result.RemoveAllConformers()
    result.AddConformer(conformer, assignId=True)
    return Chem.RemoveHs(result)


def generate(args: argparse.Namespace) -> None:
    smiles = load_smiles(args.input)
    if not 0 <= args.index < len(smiles):
        raise ValueError(f"Index {args.index} outside 0..{len(smiles) - 1}")
    molecule_smiles = smiles[args.index]
    molecule = Chem.AddHs(Chem.MolFromSmiles(molecule_smiles))
    if molecule is None:
        raise ValueError(f"Unreadable SMILES: {molecule_smiles}")
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = int(args.seed + args.index)
    parameters.useRandomCoords = True
    parameters.pruneRmsThresh = -1.0
    parameters.numThreads = int(args.threads)
    conformer_ids = list(
        AllChem.EmbedMultipleConfs(
            molecule,
            numConfs=args.num_conformers,
            params=parameters,
        )
    )
    if len(conformer_ids) != args.num_conformers:
        raise RuntimeError(
            f"Embedded {len(conformer_ids)}/{args.num_conformers} conformers"
        )
    raw = [single_conformer(molecule, conformer_id) for conformer_id in conformer_ids]
    if not AllChem.MMFFHasAllMoleculeParams(molecule):
        raise RuntimeError(f"MMFF94s parameters unavailable for {molecule_smiles}")
    statuses = AllChem.MMFFOptimizeMoleculeConfs(
        molecule,
        numThreads=args.threads,
        maxIters=args.max_iterations,
        mmffVariant="MMFF94s",
    )
    if len(statuses) != args.num_conformers:
        raise RuntimeError("MMFF result count differs from conformer count")
    minimized = [
        single_conformer(molecule, conformer_id) for conformer_id in conformer_ids
    ]
    payload = {
        "smiles": molecule_smiles,
        "raw": raw,
        "minimized": minimized,
        "mmff_status_counts": {
            str(status): sum(observed == status for observed, _ in statuses)
            for status in sorted({observed for observed, _ in statuses})
        },
        "mmff_energy_min": min(float(energy) for _, energy in statuses),
        "mmff_energy_max": max(float(energy) for _, energy in statuses),
    }
    atomic_pickle(args.output / "parts" / f"{args.index:03d}.pickle", payload)
    print(
        f"[{args.index:03d}] conformers={len(raw)} statuses={payload['mmff_status_counts']}",
        flush=True,
    )


def merge(args: argparse.Namespace) -> None:
    smiles = load_smiles(args.input)
    raw: OrderedDict[str, list[Chem.Mol]] = OrderedDict()
    minimized: OrderedDict[str, list[Chem.Mol]] = OrderedDict()
    metadata = []
    for index, expected_smiles in enumerate(smiles):
        part_path = args.output / "parts" / f"{index:03d}.pickle"
        with part_path.open("rb") as handle:
            part = pickle.load(handle)
        if part["smiles"] != expected_smiles:
            raise ValueError(f"Part {index} SMILES/order mismatch")
        raw[expected_smiles] = part["raw"]
        minimized[expected_smiles] = part["minimized"]
        metadata.append(
            {
                "index": index,
                "smiles": expected_smiles,
                "raw_count": len(part["raw"]),
                "minimized_count": len(part["minimized"]),
                "mmff_status_counts": part["mmff_status_counts"],
                "mmff_energy_min": part["mmff_energy_min"],
                "mmff_energy_max": part["mmff_energy_max"],
            }
        )
    atomic_pickle(args.output / "rdkit_etkdg_raw.pickle", raw)
    atomic_pickle(args.output / "rdkit_etkdg_mmff94s.pickle", minimized)
    (args.output / "generation_manifest.json").write_text(
        json.dumps(
            {
                "input": str(args.input.resolve()),
                "molecule_count": len(smiles),
                "num_conformers_per_molecule": len(next(iter(raw.values()))),
                "embedding": "RDKit ETKDGv3, coordinate-free molecular identity",
                "minimization": "MMFF94s, 200 iterations; nonconvergence retained",
                "failed_embeddings_retained": False,
                "molecules": metadata,
            },
            indent=2,
        )
        + "\n"
    )


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("command", choices=("generate", "merge"))
    result.add_argument("--input", type=Path, required=True)
    result.add_argument("--output", type=Path, required=True)
    result.add_argument("--index", type=int, default=0)
    result.add_argument("--num-conformers", type=int, default=1000)
    result.add_argument("--threads", type=int, default=8)
    result.add_argument("--max-iterations", type=int, default=200)
    result.add_argument("--seed", type=int, default=20261116)
    return result


if __name__ == "__main__":
    arguments = parser().parse_args()
    (generate if arguments.command == "generate" else merge)(arguments)
