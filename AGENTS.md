# AGENTS.md

## Cursor Cloud specific instructions

This is a single-product Python repo (`casf-benchmark`): a CASF–ChEMBL3D conformer
benchmark pipeline plus an interactive **Streamlit dashboard**, which is the only
long-running service. Persistence is embedded **SQLite** (no DB server); precomputed
results are bundled under `data/results/` (real files in this repo, so `git lfs pull`
is not needed to run the dashboard).

### Environment
- Dependencies are installed into a virtualenv at `.venv` (gitignored). The startup
  update script (re)creates it and runs `pip install -e ".[dev]"`. The upstream README
  uses conda, but conda is not required here — pip covers the dashboard + tests.
- `.[dev]` pulls in `streamlit`, `pytest`, and (via `posebusters`) `rdkit` from PyPI,
  so the rdkit-dependent tests run without conda. `zarr` (conda-only, analysis
  pipeline) is NOT installed; the cluster generation/analysis pipeline is not runnable
  here anyway (needs Slurm + the `/mnt/weka/mbedrosian` dataset).
- Activate the venv before running anything: `source .venv/bin/activate`.

### Test / run
- Test: `pytest` (config in `pyproject.toml`; 68 pass, 1 skipped).
  - Gotcha: `tests/test_standalone_wiring.py` invokes console scripts (e.g.
    `casf-analyze-conformer-sets`) via `subprocess` and relies on them being on `PATH`.
    Run tests with the venv **activated** (not `.venv/bin/python -m pytest`), otherwise
    `.venv/bin` is off `PATH` and those two tests fail with `FileNotFoundError`.
- Run the dashboard (the product): `streamlit run apps/dashboard/streamlit_app.py`.
  It reads `data/results/casf_analysis_dashboard.sqlite` by default; override with
  `CASF_DASHBOARD_DB` (and `CASF_EXTENDED_DB` for the extended tab, which uses
  `data/results/extended_casf_analysis.sqlite`).
- No linter/formatter is configured in this repo.
