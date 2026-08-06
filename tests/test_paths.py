from __future__ import annotations

from pathlib import Path

from casf_benchmark.paths import (
    BUNDLED_CHEMBL_MAP_CSV,
    BUNDLED_CHEMBL_REF_MAP_CSV,
    DEFAULT_CORE_PHARMA_ROOT,
    DEFAULT_REF_PHARMA_ROOT,
    DEFAULT_RUNS_ROOT,
    DEFAULT_QWEN_GENERATION_ROOT,
    DEFAULT_QWEN_REF_GENERATION_ROOT,
    DEFAULT_WEKA_ANALYSIS_SOURCES_CONFIG,
    GENERATION_SCRIPT,
    WEKA_CHEMBL_MAP_CSV,
    WEKA_CHEMBL_REF_MAP_CSV,
    WEKA_DATA_ROOT,
    WEKA_PHARMA_ROOT,
    _first_existing_path,
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


def test_weka_mapping_csv_paths_match_cluster_layout() -> None:
    assert WEKA_CHEMBL_MAP_CSV == WEKA_DATA_ROOT / "casf16" / "casf16_core_chembl3d_exact_intersection.csv"
    assert WEKA_CHEMBL_REF_MAP_CSV == WEKA_DATA_ROOT / "casf16" / "casf16_ref_chembl3d_exact_intersection.csv"


def test_weka_pharma_roots_match_submit_casf_defaults() -> None:
    assert DEFAULT_CORE_PHARMA_ROOT == WEKA_PHARMA_ROOT / "core_pb_full_dynamic_chembl_count"
    assert DEFAULT_REF_PHARMA_ROOT == WEKA_PHARMA_ROOT / "ref_pb_full_dynamic_chembl_count"


def test_weka_qwen_roots_match_analysis_config() -> None:
    assert DEFAULT_QWEN_GENERATION_ROOT.name == "qwen_gens"
    assert DEFAULT_QWEN_REF_GENERATION_ROOT.name == "casf16_ref_qwen_1k"
    assert DEFAULT_WEKA_ANALYSIS_SOURCES_CONFIG.name == "casf_analysis_sources.weka.yaml"


def test_first_existing_path_prefers_weka_mapping_when_present(tmp_path: Path, monkeypatch) -> None:
    bundled = tmp_path / "bundled.csv"
    weka = tmp_path / "weka.csv"
    bundled.write_text("bundled\n", encoding="utf-8")
    assert _first_existing_path(weka, bundled) == bundled
    weka.write_text("weka\n", encoding="utf-8")
    assert _first_existing_path(weka, bundled) == weka


def test_bundled_mapping_csvs_exist_in_repo() -> None:
    assert BUNDLED_CHEMBL_MAP_CSV.is_file()
    assert BUNDLED_CHEMBL_REF_MAP_CSV.is_file()
