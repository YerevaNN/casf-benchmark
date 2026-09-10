"""The druglike ingest must accept any generator's label, not just Qwen checkpoints.

Menua runs the same COV/MAT eval for RDKit, loqi and others and pushes the rows into
the same two dashboard tables, so adding a method has to be a catalog entry plus its
eval CSVs. This drives `scripts/build_druglike_covmat.py` over a synthetic run root
holding one Qwen checkpoint and one non-Qwen generator, and checks both land.
"""

from __future__ import annotations

import importlib.util
import sqlite3
import sys
from pathlib import Path
from types import ModuleType

import pandas as pd
import pytest

from casf_benchmark.catalog import GenerationRun
from casf_benchmark.paths import REPO_ROOT

COVMAT_REPORT = """\
================================================================================
COVMAT EVALUATION RESULTS
================================================================================

EVALUATION SUMMARY
----------------------------------------
Total molecules in ground truth: 23
Total conformers in ground truth: 2450
Missing molecules (no conformers): {missing}

COVERAGE AND RECALL METRICS
----------------------------------------
Threshold: 0.75
Molecule success rate: {success}
Coverage-Recall (COV-R):
  Mean:   {cov_r}
  Median: 0.7787
Coverage-Precision (COV-P):
  Mean:   0.5210
  Median: 0.5310
Matching-Recall (MAT-R):
  Mean:   1.2700
  Median: 1.1000
Matching-Precision (MAT-P):
  Mean:   1.4000
  Median: 1.3000
"""

RMSD_COLUMNS = (
    "geom_smiles,min_rmsd,max_rmsd,avg_rmsd,cov_r_075,cov_p_075,"
    "mat_r,mat_p,num_true_confs,num_gen_confs,num_valid_rmsd_pairs,sub_smiles"
)


def load_script() -> ModuleType:
    path = REPO_ROOT / "scripts" / "build_druglike_covmat.py"
    spec = importlib.util.spec_from_file_location("build_druglike_covmat", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_eval_outputs(run_dir: Path, smiles: str, *, cov_r: str) -> None:
    eval_dir = run_dir / "eval_20260901_000000"
    eval_dir.mkdir(parents=True)
    (eval_dir / "covmat_results.txt").write_text(
        COVMAT_REPORT.format(cov_r=cov_r, missing="1", success="0.9130")
    )
    (eval_dir / "rmsd_matrix.csv").write_text(
        f"{RMSD_COLUMNS}\n{smiles},0.3051,0.3395,0.3223,1,1,0.3051,0.3223,1,2,2,{smiles}\n"
    )


@pytest.fixture
def ingest_env(tmp_path: Path) -> dict[str, object]:
    smiles = "Cc1ccc(N)c2c1OCCC2=O"
    gen_root = tmp_path / "gen_results"
    write_eval_outputs(gen_root / "qwen_run_dir", smiles, cov_r="0.6859")
    write_eval_outputs(gen_root / "rdkit_run_dir", smiles, cov_r="0.4210")

    # Stands in for the druglike_eval.sqlite that build_druglike_eval_db.py writes.
    druglike_db = tmp_path / "druglike_eval.sqlite"
    labels = ["qwen_1p7b_fsq_bigdata_step47023", "rdkit_etkdg"]
    con = sqlite3.connect(druglike_db)
    pd.DataFrame(
        {"label": labels, "n_molecules": [23, 23], "overall_pb_pass_rate": [0.77, 0.61]}
    ).to_sql("summary", con, index=False)
    pd.DataFrame(
        {"label": labels, "smiles": [smiles, smiles], "pb_pass_rate": [0.8, 0.6]}
    ).to_sql("per_molecule", con, index=False)
    con.close()

    return {
        "gen_root": gen_root,
        "druglike_db": druglike_db,
        "extended_db": tmp_path / "extended.sqlite",
        "smiles": smiles,
    }


def run_ingest(env: dict[str, object], entries: tuple[GenerationRun, ...], monkeypatch) -> None:
    module = load_script()
    monkeypatch.setattr(module, "load_generation_run_entries", lambda: entries)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "build_druglike_covmat.py",
            "--generation-results-root",
            str(env["gen_root"]),
            "--druglike-db",
            str(env["druglike_db"]),
            "--extended-db",
            str(env["extended_db"]),
        ],
    )
    module.main()


def test_publishes_qwen_and_non_qwen_rows_side_by_side(ingest_env, monkeypatch) -> None:
    entries = (
        GenerationRun(
            label="qwen_1p7b_fsq_bigdata_step47023",
            cohorts={"druglike": "qwen_run_dir"},
            descriptors={},
        ),
        GenerationRun(
            label="rdkit_etkdg",
            cohorts={"druglike": "rdkit_run_dir"},
            descriptors={"generator": "RDKit", "display_label": "RDKit ETKDG"},
        ),
    )
    run_ingest(ingest_env, entries, monkeypatch)

    con = sqlite3.connect(ingest_env["extended_db"])
    published = pd.read_sql("select * from extended_druglike_summary", con)
    per_molecule = pd.read_sql("select * from extended_druglike_per_molecule", con)
    con.close()

    # Descriptors lead the published table, so the tab reads generator-first.
    assert list(published.columns[:6]) == [
        "generator",
        "display_label",
        "model_size",
        "tokenizer",
        "recipe",
        "step",
    ]

    summary = published.set_index("label")
    assert set(summary.index) == {"qwen_1p7b_fsq_bigdata_step47023", "rdkit_etkdg"}

    # The Qwen label carries its own metadata and is parsed from the label alone.
    qwen = summary.loc["qwen_1p7b_fsq_bigdata_step47023"]
    assert qwen["generator"] == "Qwen"
    assert qwen["model_size"] == "1.7B"
    assert qwen["tokenizer"] == "FSQ"
    assert qwen["step"] == 47023
    assert qwen["cov_r_mean"] == pytest.approx(0.6859)

    # The non-Qwen label gets its descriptors from the catalog and stays otherwise blank
    # rather than being mangled into Qwen-shaped fields.
    rdkit = summary.loc["rdkit_etkdg"]
    assert rdkit["generator"] == "RDKit"
    assert rdkit["display_label"] == "RDKit ETKDG"
    assert pd.isna(rdkit["model_size"])
    assert pd.isna(rdkit["step"])
    assert rdkit["cov_r_mean"] == pytest.approx(0.4210)
    # Metrics common to both generators are populated for both.
    assert rdkit["mat_r_mean"] == pytest.approx(1.27)
    assert rdkit["molecule_success_rate"] == pytest.approx(0.9130)

    # Per-molecule COV/MAT merged for both labels.
    assert set(per_molecule["label"]) == {"qwen_1p7b_fsq_bigdata_step47023", "rdkit_etkdg"}
    assert bool(per_molecule["mat_r"].notna().all())


def test_rerunning_does_not_duplicate_descriptor_columns(ingest_env, monkeypatch) -> None:
    entries = (
        GenerationRun(label="rdkit_etkdg", cohorts={"druglike": "rdkit_run_dir"}, descriptors={}),
    )
    run_ingest(ingest_env, entries, monkeypatch)
    # The script writes its enrichment back into druglike_eval.sqlite, so a second pass
    # reads its own output; columns must not accumulate _x/_y suffixes.
    run_ingest(ingest_env, entries, monkeypatch)

    con = sqlite3.connect(ingest_env["extended_db"])
    summary = pd.read_sql("select * from extended_druglike_summary", con)
    con.close()

    assert not [c for c in summary.columns if c.endswith(("_x", "_y"))]
    assert list(summary.columns).count("generator") == 1
    assert summary.loc[summary["label"] == "rdkit_etkdg", "generator"].iloc[0] == "rdkit"
