"""Directory layout for CASF16 conformer_sets datasets."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

METHODS = (
    "rdkit_random_raw_fixed",
    "rdkit_random_minimized_fixed",
    "torsion_raw_fixed",
    "torsion_minimized_fixed",
    "rdkit_random_raw_dynamic",
    "rdkit_random_minimized_dynamic",
    "torsion_raw_dynamic",
    "torsion_minimized_dynamic",
    "rdkit_random_raw_chembl_count",
    "rdkit_random_minimized_chembl_count",
    "torsion_raw_chembl_count",
    "torsion_minimized_chembl_count",
)

GENERATION_SUBDIR = "generation"
ANALYSIS_SUBDIR = "analysis"
TABLES_SUBDIR = "tables"
CACHE_SUBDIR = "cache"
GEOMETRIC_LIGAND_PARTS_DIRNAME = "geometric_generation_parts"
CHEMBL3D_MOL_CACHE_DIRNAME = "chembl3d_mols"

REPO_ROOT = Path(__file__).resolve().parents[2]
PACKAGE_ROOT = Path(__file__).resolve().parent
DATA_ROOT = REPO_ROOT / "data"
RESULTS_ROOT = DATA_ROOT / "results"
MAPPING_ROOT = DATA_ROOT / "mapping"

GENERATION_SCRIPT = PACKAGE_ROOT / "generation" / "conformer_sets.py"

WEKA_ROOT = Path(os.environ.get("CASF_BENCHMARK_DATA_ROOT", "/mnt/weka/mbedrosian"))
WEKA_DATA_ROOT = WEKA_ROOT / "data"
WEKA_PHARMA_ROOT = WEKA_ROOT / "pharma_generation_analysis"
WEKA_CODEX_ROOT = WEKA_ROOT / "codex_dir"

DEFAULT_DASHBOARD_DB = RESULTS_ROOT / "casf_analysis_dashboard.sqlite"
DEFAULT_EXTENDED_DB = RESULTS_ROOT / "extended_casf_analysis.sqlite"
DEFAULT_MASTER_CSV = RESULTS_ROOT / "casf_analysis_master.csv"
DEFAULT_PER_LIGAND_LONG_CSV = RESULTS_ROOT / "casf_per_ligand_long.csv"
DEFAULT_RUNS_ROOT = RESULTS_ROOT / "runs"

DEFAULT_CASF16_DATA = WEKA_DATA_ROOT / "casf16" / "CASF16"
DEFAULT_CASF16_REF_DATA = WEKA_DATA_ROOT / "casf16" / "CASF16_REF"
DEFAULT_CASF_LIGAND_DIR = DEFAULT_CASF16_DATA / "ligands"
DEFAULT_CASF_OPT_LIGAND_DIR = DEFAULT_CASF16_DATA / "ligands_opt"
DEFAULT_CASF_REF_LIGAND_DIR = DEFAULT_CASF16_REF_DATA / "ligands"
DEFAULT_CORE_LIGAND_DIR = DEFAULT_CASF16_DATA / "core_chembl3d_exact_intersection_ligands"
DEFAULT_REF_INTERSECTION_LIGAND_DIR = DEFAULT_CASF16_REF_DATA / "ref_chembl3d_exact_intersection_ligands"

DEFAULT_CHEMBL_MAP_CSV = MAPPING_ROOT / "casf16_core_chembl3d_exact_intersection.csv"
DEFAULT_CHEMBL_REF_MAP_CSV = MAPPING_ROOT / "casf16_ref_chembl3d_exact_intersection.csv"
DEFAULT_CHEMBL_DATASET_ROOT = WEKA_DATA_ROOT / "chembl3d"
DEFAULT_CHEMBL3D_INDEX_DIR = WEKA_DATA_ROOT / "chembl3d_index"
DEFAULT_CHEMBL3D_INDEX_CSV = DEFAULT_CHEMBL3D_INDEX_DIR / "chembl3d_topology_smiles_index.csv"

DEFAULT_CORE_PHARMA_ROOT = WEKA_PHARMA_ROOT / "core_pb_full_dynamic_chembl_count"
DEFAULT_REF_PHARMA_ROOT = WEKA_PHARMA_ROOT / "ref_pb_full_dynamic_chembl_count"
DEFAULT_REFERENCE_DATASETS_CORE_ROOT = WEKA_PHARMA_ROOT / "reference_datasets" / "core"
DEFAULT_REFERENCE_DATASETS_REF_ROOT = WEKA_PHARMA_ROOT / "reference_datasets" / "ref"

DEFAULT_QWEN_GENERATION_ROOT = WEKA_CODEX_ROOT / "qwen_gens"
DEFAULT_QWEN_REF_GENERATION_ROOT = WEKA_CODEX_ROOT / "qwen" / "generations" / "casf16_ref_qwen_1k"
DEFAULT_WEKA_DASHBOARD_DB = WEKA_PHARMA_ROOT / "casf_analysis_dashboard.sqlite"
DEFAULT_WEKA_MASTER_CSV = WEKA_PHARMA_ROOT / "casf_analysis_master.csv"
DEFAULT_WEKA_PER_LIGAND_LONG_CSV = WEKA_PHARMA_ROOT / "casf_per_ligand_long.csv"
DEFAULT_WEKA_ANALYSIS_SOURCES_CONFIG = REPO_ROOT / "config" / "casf_analysis_sources.weka.yaml"
WEKA_LOQI_CORE_ROOT = WEKA_CODEX_ROOT / "loqi" / "generations" / "casf16_core_loqi_1k"
WEKA_LOQI_REF_ROOT = WEKA_CODEX_ROOT / "loqi" / "generations" / "casf16_ref_loqi_1k"
WEKA_NEXTMOL_DMT_L_CORE_ROOT = (
    WEKA_CODEX_ROOT / "nextmol_dmt_l" / "generations" / "casf16_core_nextmol_dmt_l_1k"
)
WEKA_NEXTMOL_DMT_L_REF_ROOT = (
    WEKA_CODEX_ROOT / "nextmol_dmt_l" / "generations" / "casf16_ref_nextmol_dmt_l_1k"
)
WEKA_TORSIONAL_DIFFUSION_CORE_ROOT = (
    WEKA_CODEX_ROOT / "torsional_diffusion" / "generations" / "casf16_core_torsional_diffusion_1k"
)
WEKA_TORSIONAL_DIFFUSION_REF_ROOT = (
    WEKA_CODEX_ROOT / "torsional_diffusion" / "generations" / "casf16_ref_torsional_diffusion_1k"
)
WEKA_MCF_DRUGS_L_CORE_ROOT = (
    WEKA_CODEX_ROOT / "mcf_drugs_l" / "generations" / "casf16_core_mcf_drugs_l_1k"
)
WEKA_MCF_DRUGS_L_REF_ROOT = (
    WEKA_CODEX_ROOT / "mcf_drugs_l" / "generations" / "casf16_ref_mcf_drugs_l_1k"
)

DEFAULT_CONFORMER_SETS_ROOT = DEFAULT_CORE_PHARMA_ROOT
DEFAULT_GENERATION_DIR = DEFAULT_CORE_PHARMA_ROOT / GENERATION_SUBDIR
DEFAULT_CORE_INTERSECTION_ROOT = DEFAULT_CORE_PHARMA_ROOT


def _has_method_dirs(directory: Path) -> bool:
    return any((directory / method).is_dir() for method in METHODS)


def resolve_generation_dir(output_dir: Path) -> Path:
    """Return the directory where SDFs and manifest files are written."""
    output_dir = output_dir.resolve()
    if output_dir.name == GENERATION_SUBDIR:
        return output_dir
    if _has_method_dirs(output_dir):
        return output_dir
    return output_dir / GENERATION_SUBDIR


@dataclass(frozen=True)
class ConformerSetsPaths:
    root_dir: Path
    generation_dir: Path
    analysis_dir: Path
    casf_ligand_dir: Path
    tables_dir: Path
    cache_dir: Path

    @property
    def manifest_path(self) -> Path:
        return self.generation_dir / "manifest.tsv"

    @property
    def manifest_parts_dir(self) -> Path:
        return self.generation_dir / "manifest_parts"


def resolve_paths(root_dir: Path, casf_ligand_dir: Path | None = None) -> ConformerSetsPaths:
    """Resolve generation/ and analysis/ paths under a conformer_sets root."""
    root_dir = root_dir.resolve()
    if (root_dir / GENERATION_SUBDIR).is_dir():
        generation_dir = root_dir / GENERATION_SUBDIR
        analysis_dir = root_dir / ANALYSIS_SUBDIR
    elif _has_method_dirs(root_dir):
        generation_dir = root_dir
        analysis_dir = root_dir
    else:
        generation_dir = root_dir / GENERATION_SUBDIR
        analysis_dir = root_dir / ANALYSIS_SUBDIR

    if casf_ligand_dir is None:
        casf_ligand_dir = root_dir.parent / "ligands"
    casf_ligand_dir = casf_ligand_dir.resolve()

    if analysis_dir == root_dir:
        tables_dir = analysis_dir
        cache_dir = analysis_dir
    else:
        tables_dir = analysis_dir / TABLES_SUBDIR
        cache_dir = analysis_dir / CACHE_SUBDIR

    return ConformerSetsPaths(
        root_dir=root_dir,
        generation_dir=generation_dir,
        analysis_dir=analysis_dir,
        casf_ligand_dir=casf_ligand_dir,
        tables_dir=tables_dir,
        cache_dir=cache_dir,
    )


@dataclass(frozen=True)
class GeometricAnalysisPaths:
    base: ConformerSetsPaths
    report_md: Path
    ligand_metrics_parts_dir: Path
    chembl_mol_cache_dir: Path

    @property
    def root_dir(self) -> Path:
        return self.base.root_dir

    @property
    def generation_dir(self) -> Path:
        return self.base.generation_dir

    @property
    def analysis_dir(self) -> Path:
        return self.base.analysis_dir

    @property
    def casf_ligand_dir(self) -> Path:
        return self.base.casf_ligand_dir

    @property
    def manifest_path(self) -> Path:
        return self.base.manifest_path

    @property
    def manifest_parts_dir(self) -> Path:
        return self.base.manifest_parts_dir

    @property
    def tables_dir(self) -> Path:
        return self.base.tables_dir

    @property
    def cache_dir(self) -> Path:
        return self.base.cache_dir


def resolve_geometric_paths(
    root_dir: Path,
    casf_ligand_dir: Path | None = None,
) -> GeometricAnalysisPaths:
    base = resolve_paths(root_dir, casf_ligand_dir)
    cache_dir = base.cache_dir
    return GeometricAnalysisPaths(
        base=base,
        report_md=base.analysis_dir / "geometric_report.md",
        ligand_metrics_parts_dir=cache_dir / GEOMETRIC_LIGAND_PARTS_DIRNAME,
        chembl_mol_cache_dir=cache_dir / CHEMBL3D_MOL_CACHE_DIRNAME,
    )
