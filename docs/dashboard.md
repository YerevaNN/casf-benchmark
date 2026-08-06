# Dashboard

Read-only Streamlit UI over precomputed analysis SQLite.

**App:** [`apps/dashboard/streamlit_app.py`](../apps/dashboard/streamlit_app.py)  
**Build scripts:** [`scripts/build_casf_analysis_master_csv.py`](../scripts/build_casf_analysis_master_csv.py) · [`scripts/build_casf_analysis_dashboard_db.py`](../scripts/build_casf_analysis_dashboard_db.py)  
**Config:** [`src/casf_benchmark/config/casf_analysis_sources.yaml`](../src/casf_benchmark/config/casf_analysis_sources.yaml) (repo-relative) · [`src/casf_benchmark/config/casf_analysis_sources.weka.yaml`](../src/casf_benchmark/config/casf_analysis_sources.weka.yaml) (cluster)

## Data flow

```
data/results/runs/{run_id}/analysis/tables/geometric_per_ligand_long.csv
        →  build_casf_analysis_master_csv.py
        →  build_casf_analysis_dashboard_db.py
        →  data/results/casf_analysis_dashboard.sqlite
        →  streamlit run apps/dashboard/streamlit_app.py
```

Bundled DB for clone-and-run: [`data/results/casf_analysis_dashboard.sqlite`](../data/results/casf_analysis_dashboard.sqlite).

Override: `CASF_DASHBOARD_DB`, `CASF_EXTENDED_DB` (extended tab — [`data/results/extended_casf_analysis.sqlite`](../data/results/extended_casf_analysis.sqlite)).

## SQLite tables

| Table | Contents |
| --- | --- |
| `comparison_rows` | Cohort aggregates by family/tier/method |
| `comparison_strata` | Same, split by rotatable-bond / heavy-atom bins |
| `per_ligand_long` | Full per-ligand metrics |
| `analysis_sources` | Run roots from YAML |
| `catalog_families` | From [`casf_generation_families.yaml`](../src/casf_benchmark/config/casf_generation_families.yaml) |

## UI

Sidebar: ligand set, tier, family, stratum breakdown. Tabs: Overview, Clustering, Energy, CASF hits, Funnel. Extended section loads optional `extended_*` tables.

Per-table **Columns & rows** expander: hide columns; filter by method, stratum, or `mol_id` substring. `method` column stays pinned on horizontal scroll.

## Commands

```bash
# View (minimal env — see extras.md)
streamlit run apps/dashboard/streamlit_app.py

# Rebuild DB from bundled per-run CSVs (after pip install -e ".[dev]")
casf-build-master-csv
casf-build-dashboard-db

# Full cluster rebuild (all sources + reference runs)
sbatch scripts/run_pb_once_casf_analysis.sbatch
```

Adding a generator: entry in `casf_generation_families.yaml` + analyzed run root in `casf_analysis_sources.yaml` → rebuild scripts above.

Public hosting: [extras.md](extras.md#streamlit-cloud).

See also: [analyzer.md](analyzer.md) · [extras.md](extras.md)
