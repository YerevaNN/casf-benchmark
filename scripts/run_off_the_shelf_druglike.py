#!/usr/bin/env python3
"""Run an off-the-shelf conformer model on the 23-molecule druglike set.

The historical model adapters live under ``/mnt/weka/mbedrosian/codex_dir`` and
were written against the old ``molgen3D.pharmacophore`` package. This driver
installs narrow import aliases to the equivalent ``casf_benchmark`` modules,
loads only the model-specific sampler from each adapter, and writes the
``generation_results.pickle`` contract used by the druglike evaluators.

Array jobs normally process one molecule each into ``parts/``. Run with
``--merge-parts`` after the array completes to validate and combine all parts.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import pickle
import sys
import tempfile
import types
from collections import OrderedDict
from pathlib import Path
from typing import Any

DEFAULT_DRUGLIKE_PICKLE = Path("/mnt/weka/vtarasov/druglike_smi.pickle")
CODEX_ROOT = Path("/mnt/weka/mbedrosian/codex_dir")

MODEL_ADAPTERS = {
    "loqi": CODEX_ROOT / "loqi/codex/loqi_casf_generate.py",
    "nextmol_dmt_l": CODEX_ROOT / "nextmol_dmt_l/codex/nextmol_dmt_l_casf_generate.py",
    "torsional_diffusion": (
        CODEX_ROOT / "torsional_diffusion/codex/torsional_diffusion_casf_generate.py"
    ),
    "mcf_drugs_l": CODEX_ROOT / "mcf_drugs_l/codex/mcf_casf_generate.py",
}


def load_druglike_smiles(path: Path) -> list[str]:
    with path.open("rb") as handle:
        data = pickle.load(handle)
    if not isinstance(data, dict):
        raise TypeError(f"Expected a SMILES-keyed mapping in {path}, got {type(data).__name__}")
    smiles = list(data)
    if not smiles or not all(isinstance(value, str) and value for value in smiles):
        raise ValueError(f"{path} does not contain non-empty SMILES keys")
    return smiles


def part_path(output_dir: Path, molecule_index: int) -> Path:
    return output_dir / "parts" / f"{molecule_index:03d}.pickle"


def _atomic_pickle_dump(value: object, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="wb", dir=destination.parent, prefix=f".{destination.name}.", delete=False
    ) as handle:
        temporary = Path(handle.name)
        pickle.dump(value, handle, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(temporary, destination)


def normalize_conformers(mols: list) -> list:
    """Remove explicit hydrogens so PB and COV/MAT use the shared heavy topology."""
    from rdkit import Chem

    normalized = []
    for mol in mols:
        if mol is None or mol.GetNumConformers() == 0:
            continue
        normalized.append(Chem.RemoveHs(mol))
    return normalized


def merge_parts(druglike_pickle: Path, output_dir: Path) -> Path:
    smiles = load_druglike_smiles(druglike_pickle)
    merged: OrderedDict[str, list] = OrderedDict()
    metadata = []
    missing = []
    for index, expected_smiles in enumerate(smiles):
        path = part_path(output_dir, index)
        if not path.is_file():
            missing.append(path)
            continue
        with path.open("rb") as handle:
            part = pickle.load(handle)
        if not isinstance(part, dict) or list(part) != [expected_smiles]:
            raise ValueError(
                f"{path} must contain exactly molecule {index} ({expected_smiles!r}); "
                f"found keys {list(part) if isinstance(part, dict) else type(part).__name__}"
            )
        mols = part[expected_smiles]
        if not isinstance(mols, list):
            raise TypeError(f"{path}: conformers must be a list")
        mols = normalize_conformers(mols)
        merged[expected_smiles] = mols
        metadata.append(
            {"molecule_index": index, "smiles": expected_smiles, "n_conformers": len(mols)}
        )
    if missing:
        preview = ", ".join(str(path) for path in missing[:5])
        raise FileNotFoundError(f"Missing {len(missing)} inference part(s): {preview}")

    destination = output_dir / "generation_results.pickle"
    _atomic_pickle_dump(merged, destination)
    (output_dir / "generation_summary.json").write_text(
        json.dumps(
            {
                "n_molecules": len(merged),
                "total_conformers": sum(len(mols) for mols in merged.values()),
                "molecules": metadata,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return destination


def install_legacy_casf_aliases() -> None:
    """Expose moved CASF helpers under the import names expected by old adapters."""
    import molgen3D
    from casf_benchmark import paths
    from casf_benchmark.chembl3d import loader
    from casf_benchmark.generation import conformer_sets

    package_name = "molgen3D.pharmacophore"
    package = types.ModuleType(package_name)
    package.__path__ = []  # type: ignore[attr-defined]
    package.chembl3d_conformer_loader = loader
    package.conformer_sets_layout = paths
    package.generate_casf_smiles_conformer_sets = conformer_sets

    setattr(molgen3D, "pharmacophore", package)
    sys.modules[package_name] = package
    sys.modules[f"{package_name}.chembl3d_conformer_loader"] = loader
    sys.modules[f"{package_name}.conformer_sets_layout"] = paths
    sys.modules[f"{package_name}.generate_casf_smiles_conformer_sets"] = conformer_sets


def load_adapter(model: str):
    install_legacy_casf_aliases()
    path = MODEL_ADAPTERS[model]
    if not path.is_file():
        raise FileNotFoundError(f"Missing {model} adapter: {path}")
    module_name = f"_casf_external_{model}"
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load adapter module from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def build_sampler(model: str, args: argparse.Namespace) -> Any:
    adapter = load_adapter(model)
    if model == "loqi":
        previous_cwd = Path.cwd()
        os.chdir(CODEX_ROOT / "loqi")
        try:
            return adapter.LoqiSampler(
                config_path=CODEX_ROOT / "loqi/scripts/conf/loqi/loqi.yaml",
                ckpt_path=CODEX_ROOT / "loqi/data/loqi.ckpt",
                batch_size=args.batch_size,
                n_steps=args.inference_steps,
                atom_aware_batching=True,
                target_molecule_size=50,
                add_hs=True,
                use_stereo_bonds=True,
            )
        finally:
            os.chdir(previous_cwd)
    if model == "nextmol_dmt_l":
        return adapter.NextMolDmtLGenerator(
            repo_root=CODEX_ROOT / "nextmol_dmt_l/repo",
            checkpoint=CODEX_ROOT / "nextmol_dmt_l/checkpoints/drugs_dmt_l_e2999.ckpt",
            llm_model=str(CODEX_ROOT / "nextmol_dmt_l/checkpoints/mollama"),
            infer_batch_size=args.batch_size,
            sampling_steps=args.inference_steps,
            seed=args.seed,
            precision="bf16-mixed",
            device=args.device,
        )
    if model == "torsional_diffusion":
        return adapter.TorsionalDiffusionGenerator(
            model_dir=CODEX_ROOT / "torsional_diffusion/checkpoints/workdir/drugs_default",
            ckpt="best_model.pt",
            batch_size=args.batch_size,
            inference_steps=args.inference_steps,
            seed=args.seed,
            device=args.device,
            ode=False,
            no_random=False,
            pre_mmff=False,
            post_mmff=False,
        )
    if model == "mcf_drugs_l":
        return adapter.McfSampler(
            ckpt_path=CODEX_ROOT / "mcf_drugs_l/repo/ckpts/mcf_drugs_l.ckpt",
            batch_size=args.batch_size,
            n_eigenfuncs=32,
            seed=args.seed,
            sampling_fn="standard",
        )
    raise ValueError(f"Unknown model: {model}")


def generate_one(
    model: str, sampler: Any, smiles: str, num_conformers: int, molecule_index: int
) -> tuple[list, int, list[str]]:
    if model == "loqi":
        mols, errors, batches = sampler.generate(smiles, num_conformers)
        return mols, batches, list(errors)
    if model == "nextmol_dmt_l":
        mols, batches = sampler.generate(smiles, num_conformers)
        return mols, batches, []
    if model == "torsional_diffusion":
        mols, batches = sampler.generate(smiles, num_conformers, f"druglike_{molecule_index:03d}")
        return mols, batches, []
    if model == "mcf_drugs_l":
        mols, errors, batches = sampler.generate(smiles, num_conformers)
        return mols, batches, list(errors)
    raise ValueError(f"Unknown model: {model}")


def default_batch_size(model: str) -> int:
    return {
        "loqi": 512,
        "nextmol_dmt_l": 250,
        "torsional_diffusion": 128,
        "mcf_drugs_l": 250,
    }[model]


def default_inference_steps(model: str) -> int:
    return {"loqi": 25, "nextmol_dmt_l": 100, "torsional_diffusion": 20, "mcf_drugs_l": 1}[
        model
    ]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", choices=tuple(MODEL_ADAPTERS), required=True)
    parser.add_argument("--druglike-pickle", type=Path, default=DEFAULT_DRUGLIKE_PICKLE)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--num-conformers", type=int, default=1000)
    parser.add_argument("--molecule-offset", type=int, default=0)
    parser.add_argument("--limit-molecules", type=int, default=1)
    parser.add_argument("--batch-size", type=int)
    parser.add_argument("--inference-steps", type=int)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--device", choices=("cuda", "cpu"), default="cuda")
    parser.add_argument("--merge-parts", action="store_true")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.merge_parts:
        destination = merge_parts(args.druglike_pickle, args.output_dir)
        print(f"Merged inference parts -> {destination}", flush=True)
        return
    if args.num_conformers <= 0 or args.limit_molecules <= 0 or args.molecule_offset < 0:
        raise SystemExit("Conformer count/limit must be positive and offset must be non-negative")

    smiles = load_druglike_smiles(args.druglike_pickle)
    selected = list(
        enumerate(
            smiles[args.molecule_offset : args.molecule_offset + args.limit_molecules],
            start=args.molecule_offset,
        )
    )
    if not selected:
        raise SystemExit(f"No molecules selected at offset {args.molecule_offset}")
    args.batch_size = args.batch_size or default_batch_size(args.model)
    args.inference_steps = args.inference_steps or default_inference_steps(args.model)

    sampler = build_sampler(args.model, args)
    for index, molecule_smiles in selected:
        mols, batches, errors = generate_one(
            args.model, sampler, molecule_smiles, args.num_conformers, index
        )
        usable = normalize_conformers(mols)
        if errors:
            print(f"[{index:03d}] warnings={errors}", flush=True)
        if not usable:
            raise RuntimeError(f"{args.model} generated no usable conformers for molecule {index}")
        destination = part_path(args.output_dir, index)
        _atomic_pickle_dump(OrderedDict([(molecule_smiles, usable)]), destination)
        print(
            f"[{index:03d}] model={args.model} conformers={len(usable)}/"
            f"{args.num_conformers} batches={batches} -> {destination}",
            flush=True,
        )


if __name__ == "__main__":
    main()
