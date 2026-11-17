#!/usr/bin/env python3
"""Top up the PLINDER-22 Qwen cache without modifying its original artifacts.

This deliberately mirrors the acceptance path in
``molgen3D.evaluation.inference.process_batch``.  Every batch is an atomic resume
unit and receives a distinct recorded seed, so an interrupted job never replays
an already committed random stream.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import pickle
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


DEFAULT_SOURCE_DIR = Path(
    "/mnt/weka/mbedrosian/pharma_generation_analysis/"
    "v2_1p7b_coordtoken_step41054/gen_results/"
    "qwen_1p7b_v2_coordtoken_step41054_plinder"
)
DEFAULT_OUTPUT_DIR = Path(
    "/mnt/weka/mbedrosian/code/casf-benchmark/data/results/"
    "pharmacophore_plinder22/qwen_topup"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_pickle_dump(path: Path, payload: object) -> None:
    import cloudpickle

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temporary.open("wb") as handle:
        cloudpickle.dump(payload, handle, protocol=4)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def atomic_json_dump(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def load_pickle(path: Path) -> Any:
    import cloudpickle

    with path.open("rb") as handle:
        return cloudpickle.load(handle)


def molecule_fingerprint(mol: Any) -> str:
    """Hash graph and stored coordinates without normalizing or changing them."""
    from rdkit import Chem

    flags = Chem.PropertyPickleOptions.AllProps | Chem.PropertyPickleOptions.CoordsAsDouble
    return hashlib.sha256(mol.ToBinary(flags)).hexdigest()


def smiles_id(smiles: str) -> str:
    return hashlib.sha256(smiles.encode()).hexdigest()[:16]


def load_committed_batches(batch_dir: Path) -> tuple[dict[str, list], list[dict]]:
    replacements: dict[str, list] = defaultdict(list)
    records: list[dict] = []
    for path in sorted(batch_dir.glob("batch_*.pickle")):
        payload = load_pickle(path)
        if payload["batch_file"] != path.name:
            raise ValueError(f"Batch identity mismatch in {path}")
        records.extend(payload["records"])
        for smiles, mols in payload["accepted"].items():
            replacements[smiles].extend(mols)
    return dict(replacements), records


def generate_batch(
    *,
    model: Any,
    tokenizer: Any,
    fsq_decoder: Any,
    prompts: list[tuple[str, str, str]],
    generation_config: Any,
    seed: int,
    batch_index: int,
    original_counts: dict[str, int],
    accepted_before: dict[str, int],
) -> tuple[dict[str, list], list[dict]]:
    """Generate and apply the original molecular-identity/parsing checks."""
    import torch
    from molgen3D.evaluation.inference import set_seed
    from molgen3D.evaluation.utils import extract_between, same_molecular_graph
    from molgen3D.data_processing.smiles_encoder_decoder import strip_smiles

    set_seed(seed)
    texts = [prompt for _, prompt, _ in prompts]
    tokenized = tokenizer(
        texts,
        return_tensors="pt",
        padding=True,
        pad_to_multiple_of=8,
    )
    tokenized = {
        key: value.to(model.device, non_blocking=True)
        for key, value in tokenized.items()
    }
    tokenized["attention_mask"] = tokenized["attention_mask"].contiguous()
    eos_token_id = tokenizer.convert_tokens_to_ids("[/CONFORMER]")
    if eos_token_id is None:
        eos_token_id = tokenizer.eos_token_id
    with torch.inference_mode():
        outputs = model.generate(
            input_ids=tokenized["input_ids"],
            attention_mask=tokenized["attention_mask"],
            max_new_tokens=1000,
            eos_token_id=eos_token_id,
            generation_config=generation_config,
            use_cache=True,
            return_dict_in_generate=True,
            output_attentions=False,
            output_hidden_states=False,
        )
        sequences = outputs.sequences.detach().cpu()
    decoded_outputs = tokenizer.batch_decode(sequences, skip_special_tokens=False)

    accepted: dict[str, list] = defaultdict(list)
    records: list[dict] = []
    for position, (decoded, (target_smiles, _, attempt_id)) in enumerate(
        zip(decoded_outputs, prompts, strict=True)
    ):
        canonical_smiles = extract_between(decoded, "[SMILES]", "[/SMILES]")
        generated_conformer = extract_between(decoded, "[CONFORMER]", "[/CONFORMER]")
        generated_smiles = ""
        status = "rejected"
        reason = "no_eos"
        mol_obj = None
        if generated_conformer:
            generated_smiles = strip_smiles(generated_conformer)
            if not same_molecular_graph(canonical_smiles, generated_smiles):
                reason = "smiles_mismatch"
            else:
                try:
                    mol_obj = fsq_decoder.decode_mol(generated_conformer)
                except Exception as exc:  # identical accept/reject boundary; detail is provenance
                    reason = "mol_parse_fail"
                    parse_error = f"{type(exc).__name__}: {exc}"
                else:
                    status = "accepted"
                    reason = "accepted"
        current_offset = accepted_before.get(target_smiles, 0) + len(accepted[target_smiles])
        candidate_index = original_counts[target_smiles] + current_offset
        candidate_id = (
            f"qwen_topup:{smiles_id(target_smiles)}:{candidate_index:04d}"
            if status == "accepted"
            else ""
        )
        if mol_obj is not None:
            accepted[target_smiles].append(mol_obj)
        record = {
            "attempt_id": attempt_id,
            "batch_index": batch_index,
            "batch_position": position,
            "seed": seed,
            "molecule_smiles": target_smiles,
            "smiles_id": smiles_id(target_smiles),
            "candidate_id": candidate_id,
            "candidate_index": candidate_index if candidate_id else None,
            "status": status,
            "rejection_reason": "" if status == "accepted" else reason,
            "canonical_smiles": canonical_smiles,
            "generated_smiles": generated_smiles,
            "decoded_output_sha256": hashlib.sha256(decoded.encode()).hexdigest(),
            "generated_conformer": generated_conformer,
        }
        if reason == "mol_parse_fail":
            record["parse_error"] = parse_error
        records.append(record)
    return dict(accepted), records


def write_attempt_tables(output_dir: Path, records: list[dict]) -> None:
    jsonl_path = output_dir / "attempt_log.jsonl"
    temporary = jsonl_path.with_name(f".{jsonl_path.name}.{os.getpid()}.tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record, sort_keys=True) + "\n")
    os.replace(temporary, jsonl_path)

    fields = [
        "attempt_id",
        "batch_index",
        "batch_position",
        "seed",
        "molecule_smiles",
        "smiles_id",
        "candidate_id",
        "candidate_index",
        "status",
        "rejection_reason",
        "canonical_smiles",
        "generated_smiles",
        "decoded_output_sha256",
    ]
    csv_path = output_dir / "attempt_log.csv"
    temporary = csv_path.with_name(f".{csv_path.name}.{os.getpid()}.tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fields)
        writer.writeheader()
        writer.writerows({key: record.get(key, "") for key in fields} for record in records)
    os.replace(temporary, csv_path)


def finalize(
    *,
    source_dir: Path,
    output_dir: Path,
    original: dict[str, list],
    metadata: dict,
    replacements: dict[str, list],
    records: list[dict],
    source_hashes_before: dict[str, str],
    base_seed: int,
    max_attempts_per_molecule: int,
    batch_size: int,
) -> None:
    replacement_only = {smiles: list(replacements.get(smiles, [])) for smiles in original}
    combined = {
        smiles: list(mols) + list(replacement_only[smiles])
        for smiles, mols in original.items()
    }
    final_counts = {smiles: len(mols) for smiles, mols in combined.items()}
    if any(count != 1000 for count in final_counts.values()):
        raise RuntimeError(f"Top-up incomplete: {final_counts}")

    original_fingerprints = {
        smiles: [molecule_fingerprint(mol) for mol in mols]
        for smiles, mols in original.items()
    }
    atomic_pickle_dump(output_dir / "qwen_replacements_only.pickle", replacement_only)
    atomic_pickle_dump(output_dir / "qwen_completed_1000.pickle", combined)
    reloaded = load_pickle(output_dir / "qwen_completed_1000.pickle")
    for smiles, fingerprints in original_fingerprints.items():
        observed = [
            molecule_fingerprint(mol)
            for mol in reloaded[smiles][: len(fingerprints)]
        ]
        if observed != fingerprints:
            raise RuntimeError(f"Original coordinates changed in combined cache: {smiles}")

    write_attempt_tables(output_dir, records)
    attempts = Counter(record["molecule_smiles"] for record in records)
    rejections: dict[str, Counter] = defaultdict(Counter)
    for record in records:
        if record["status"] != "accepted":
            rejections[record["molecule_smiles"]][record["rejection_reason"]] += 1
    rows = []
    for smiles, original_mols in original.items():
        rows.append(
            {
                "smiles_id": smiles_id(smiles),
                "molecule_smiles": smiles,
                "original_count": len(original_mols),
                "additional_attempts": attempts[smiles],
                "accepted_additions": len(replacement_only[smiles]),
                "smiles_mismatch": rejections[smiles]["smiles_mismatch"],
                "mol_parse_fail": rejections[smiles]["mol_parse_fail"],
                "no_eos": rejections[smiles]["no_eos"],
                "surplus_generated": 0,
                "final_count": final_counts[smiles],
            }
        )
    with (output_dir / "topup_counts.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    source_hashes_after = {
        "generation_results.pickle": sha256(source_dir / "generation_results.pickle"),
        "generation_metadata.json": sha256(source_dir / "generation_metadata.json"),
    }
    if source_hashes_after != source_hashes_before:
        raise RuntimeError("An immutable source artifact changed during top-up")
    manifest = {
        "status": "complete",
        "source_dir": str(source_dir.resolve()),
        "source_hashes_before": source_hashes_before,
        "source_hashes_after": source_hashes_after,
        "source_metadata": metadata,
        "base_seed": base_seed,
        "batch_seed_rule": "base_seed + batch_index",
        "batch_size": batch_size,
        "max_attempts_per_molecule": max_attempts_per_molecule,
        "original_counts": {smiles: len(mols) for smiles, mols in original.items()},
        "final_counts": final_counts,
        "attempts": len(records),
        "accepted_additions": sum(len(mols) for mols in replacement_only.values()),
        "rejections_by_reason": dict(
            Counter(
                record["rejection_reason"]
                for record in records
                if record["status"] != "accepted"
            )
        ),
        "surplus_samples_generated_not_retained": 0,
        "batches": [
            {
                "batch_index": batch_index,
                "seed": next(
                    record["seed"]
                    for record in records
                    if record["batch_index"] == batch_index
                ),
                "attempts": sum(
                    record["batch_index"] == batch_index for record in records
                ),
                "accepted": sum(
                    record["batch_index"] == batch_index
                    and record["status"] == "accepted"
                    for record in records
                ),
                "surplus_generated_not_retained": 0,
            }
            for batch_index in sorted({record["batch_index"] for record in records})
        ],
        "combined_cache_sha256": sha256(output_dir / "qwen_completed_1000.pickle"),
        "replacement_cache_sha256": sha256(output_dir / "qwen_replacements_only.pickle"),
        "original_prefix_coordinate_hashes_verified": True,
        "acceptance_logic": "original molecular identity and CoordToken parsing only",
        "posebusters_used_for_acceptance": False,
        "sampling": {
            "temperature": metadata["temperature"],
            "top_p": metadata["top_p"],
            "max_new_tokens": 1000,
            "generation_config": "molgen3D sampling_configs['top_p_sampling1']",
        },
    }
    atomic_json_dump(output_dir / "topup_manifest.json", manifest)


def run(args: argparse.Namespace) -> None:
    from molgen3D.config.sampling_config import sampling_configs
    from molgen3D.data_processing.coordtoken import CoordTokenCodec, FORMAT
    from molgen3D.evaluation.inference import load_model_tokenizer

    source_dir = args.source_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    batch_dir = output_dir / "batches"
    batch_dir.mkdir(exist_ok=True)
    source_hashes_before = {
        "generation_results.pickle": sha256(source_dir / "generation_results.pickle"),
        "generation_metadata.json": sha256(source_dir / "generation_metadata.json"),
    }
    metadata = json.loads((source_dir / "generation_metadata.json").read_text())
    original = load_pickle(source_dir / "generation_results.pickle")
    if metadata["temperature"] != 1.0 or metadata["top_p"] != 1.0:
        raise ValueError("Unexpected source sampling settings")
    if metadata["coordtoken_format"] != FORMAT:
        raise ValueError("Unexpected CoordToken format")
    if metadata["batch_size"] != args.batch_size:
        raise ValueError(
            f"Batch size must preserve source setting {metadata['batch_size']}"
        )
    if any(len(mols) > 1000 for mols in original.values()):
        raise ValueError("Source pool exceeds the requested target")

    replacements, records = load_committed_batches(batch_dir)
    original_counts = {smiles: len(mols) for smiles, mols in original.items()}
    accepted_counts = {smiles: len(replacements.get(smiles, [])) for smiles in original}
    remaining = {
        smiles: 1000 - original_counts[smiles] - accepted_counts[smiles]
        for smiles in original
    }
    if all(value == 0 for value in remaining.values()):
        finalize(
            source_dir=source_dir,
            output_dir=output_dir,
            original=original,
            metadata=metadata,
            replacements=replacements,
            records=records,
            source_hashes_before=source_hashes_before,
            base_seed=args.base_seed,
            max_attempts_per_molecule=args.max_attempts_per_molecule,
            batch_size=args.batch_size,
        )
        return

    model, tokenizer = load_model_tokenizer(
        model_path=metadata["derived_model_path"],
        tokenizer_path=metadata["derived_model_path"],
        torch_dtype="bfloat16",
        attention_imp="sdpa",
        device=args.device,
    )
    fsq_decoder = CoordTokenCodec(
        metadata["coordtoken_checkpoint"],
        expected_sha256=metadata["coordtoken_checkpoint_sha256"],
        device=model.device,
    )

    batch_index = len(list(batch_dir.glob("batch_*.pickle")))
    attempt_counts = Counter(record["molecule_smiles"] for record in records)
    while True:
        replacements, records = load_committed_batches(batch_dir)
        accepted_counts = {
            smiles: len(replacements.get(smiles, [])) for smiles in original
        }
        remaining = {
            smiles: 1000 - original_counts[smiles] - accepted_counts[smiles]
            for smiles in original
        }
        if all(value == 0 for value in remaining.values()):
            break
        exhausted = [
            smiles
            for smiles, missing in remaining.items()
            if missing > 0 and attempt_counts[smiles] >= args.max_attempts_per_molecule
        ]
        if exhausted:
            atomic_json_dump(
                output_dir / "topup_incomplete.json",
                {"exhausted": exhausted, "remaining": remaining},
            )
            raise RuntimeError(f"Maximum attempts reached for {len(exhausted)} molecules")

        planned: list[tuple[str, str, str]] = []
        for smiles in sorted(original, key=len):
            allowed = args.max_attempts_per_molecule - attempt_counts[smiles]
            count = min(remaining[smiles], allowed)
            for offset in range(count):
                ordinal = attempt_counts[smiles] + offset
                attempt_id = f"{smiles_id(smiles)}:{ordinal:04d}"
                planned.append((smiles, f"[SMILES]{smiles}[/SMILES]", attempt_id))
        planned = planned[: args.batch_size]
        if not planned:
            raise RuntimeError("No attempts planned while top-up remains incomplete")
        seed = args.base_seed + batch_index
        accepted, batch_records = generate_batch(
            model=model,
            tokenizer=tokenizer,
            fsq_decoder=fsq_decoder,
            prompts=planned,
            generation_config=sampling_configs["top_p_sampling1"],
            seed=seed,
            batch_index=batch_index,
            original_counts=original_counts,
            accepted_before=accepted_counts,
        )
        batch_path = batch_dir / f"batch_{batch_index:04d}.pickle"
        atomic_pickle_dump(
            batch_path,
            {
                "batch_file": batch_path.name,
                "batch_index": batch_index,
                "seed": seed,
                "prompts": planned,
                "records": batch_records,
                "accepted": accepted,
            },
        )
        attempt_counts.update(record["molecule_smiles"] for record in batch_records)
        batch_index += 1

    replacements, records = load_committed_batches(batch_dir)
    finalize(
        source_dir=source_dir,
        output_dir=output_dir,
        original=original,
        metadata=metadata,
        replacements=replacements,
        records=records,
        source_hashes_before=source_hashes_before,
        base_seed=args.base_seed,
        max_attempts_per_molecule=args.max_attempts_per_molecule,
        batch_size=args.batch_size,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--base-seed", type=int, default=20261116)
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--max-attempts-per-molecule", type=int, default=500)
    parser.add_argument("--device", choices=["cuda", "cpu"], default="cuda")
    return parser


if __name__ == "__main__":
    cli = build_parser().parse_args()
    started = time.perf_counter()
    status = 'failed'
    try:
        run(cli)
        status = 'complete'
    finally:
        runtime_path = cli.output_dir / 'runtime_sessions.json'
        sessions = json.loads(runtime_path.read_text()) if runtime_path.is_file() else []
        sessions.append({'status': status, 'device': cli.device, 'host': os.uname().nodename,
                         'slurm_job_id': os.environ.get('SLURM_JOB_ID'),
                         'elapsed_seconds': time.perf_counter() - started,
                         'batch_size': cli.batch_size, 'base_seed': cli.base_seed})
        atomic_json_dump(runtime_path, sessions)
