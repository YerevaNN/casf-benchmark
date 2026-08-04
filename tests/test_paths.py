from __future__ import annotations

from pathlib import Path

from casf_benchmark.paths import (
    DEFAULT_RUNS_ROOT,
    GENERATION_SCRIPT,
    resolve_generation_dir,
    resolve_geometric_paths,
    resolve_paths,
)


def test_generation_script_exists_in_repo() -> None:
    assert GENERATION_SCRIPT.is_file()
    assert GENERATION_SCRIPT.name == "conformer_sets.py"


def test_resolve_paths_for_bundled_run_layout() -> None:
    root = DEFAULT_RUNS_ROOT / "qwen_core"
    paths = resolve_paths(root)
    assert paths.generation_dir == root / "generation"
    assert paths.analysis_dir == root / "analysis"
    assert paths.tables_dir == paths.analysis_dir / "tables"
    assert paths.cache_dir == paths.analysis_dir / "cache"
    assert paths.manifest_path.name == "manifest.tsv"
    assert (paths.tables_dir / "geometric_per_ligand_long.csv").exists()


def test_resolve_paths_for_weka_generation_layout(tmp_path: Path) -> None:
    root = tmp_path / "core_pb_full_dynamic_chembl_count"
    generation = root / "generation"
    generation.mkdir(parents=True)
    (generation / "manifest.tsv").write_text("mol_id\tgeneration_method\tstatus\n", encoding="utf-8")

    paths = resolve_paths(root)
    assert paths.generation_dir == generation
    assert paths.analysis_dir == root / "analysis"
    assert paths.manifest_path.exists()


def test_resolve_generation_dir_accepts_root_or_generation(tmp_path: Path) -> None:
    root = tmp_path / "conformer_sets_all"
    generation = root / "generation"
    generation.mkdir(parents=True)
    assert resolve_generation_dir(root) == generation
    assert resolve_generation_dir(generation) == generation


def test_resolve_generation_dir_legacy_flat_layout(tmp_path: Path) -> None:
    root = tmp_path / "legacy"
    (root / "rdkit_random_raw_fixed").mkdir(parents=True)
    assert resolve_generation_dir(root) == root


def test_resolve_geometric_paths_v2_artifacts(tmp_path: Path) -> None:
    root = tmp_path / "dataset"
    generation = root / "generation"
    generation.mkdir(parents=True)
    (generation / "rdkit_random_raw_fixed").mkdir()
    paths = resolve_geometric_paths(root)
    assert paths.generation_dir == generation
    assert paths.analysis_dir == root / "analysis"
    assert paths.report_md == root / "analysis" / "geometric_report.md"
    assert paths.ligand_metrics_parts_dir == root / "analysis" / "cache" / "geometric_generation_parts"
    assert paths.chembl_mol_cache_dir == root / "analysis" / "cache" / "chembl3d_mols"
