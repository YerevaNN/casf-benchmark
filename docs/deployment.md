# Deploy the dashboard (Streamlit Community Cloud)

YerevaNN org repos cannot always authorize the Streamlit GitHub App. Use a **personal mirror** under [menuab/casf-benchmark](https://github.com/menuab/casf-benchmark) for public web hosting.

## Is the mirror always in sync?

**Not automatically.** Git remotes are independent:

| Remote | Repo | Role |
| --- | --- | --- |
| `origin` | `YerevaNN/casf-benchmark` | Source of truth (development) |
| `personal` | `menuab/casf-benchmark` | Streamlit Cloud deploy target |

After you push only to YerevaNN, **menuab does not update** until you also push to `personal`.

### One-command sync (recommended)

```bash
./scripts/push_with_mirror.sh
```

This runs `git push origin main` then `git push personal main`. Streamlit Cloud watches **menuab** and redeploys within a few minutes.

### Manual two-step push

```bash
git push origin main
git push personal main
```

## One-time setup

### 1. Personal mirror remote (on your machine)

```bash
cd casf-benchmark
git remote add personal git@github.com:menuab/casf-benchmark.git
```

Create an empty public repo `menuab/casf-benchmark` on GitHub first if it does not exist.

Initial mirror:

```bash
git push -u personal main
```

### 2. Streamlit Community Cloud

1. Sign in at [share.streamlit.io](https://share.streamlit.io) with **your** GitHub (`menuab`).
2. **New app** → Repository: `menuab/casf-benchmark`, branch: `main`.
3. **Main file path:** `apps/dashboard/streamlit_app.py`
4. **Python version:** 3.10
5. Deploy.

Streamlit reads `requirements.txt` at the repo root (dashboard deps + `-e .` for `casf_benchmark`).

Your public URL will look like:

```
https://casf-benchmark-menuab.streamlit.app
```

(you choose the app subdomain during setup)

## Optional: fully automatic mirror (GitHub Actions)

If you can add secrets to `YerevaNN/casf-benchmark`, a workflow can push to `menuab` on every `main` push. Requires a Personal Access Token with `repo` scope on `menuab`, stored as `PERSONAL_MIRROR_PAT`. See `.github/workflows/mirror-personal.yml` (disabled until the secret exists).

Without that secret, use `./scripts/push_with_mirror.sh` after each YerevaNN push.

## Updating the live dashboard

1. Commit changes locally.
2. Run `./scripts/push_with_mirror.sh`.
3. Wait ~1–3 minutes for Streamlit Cloud to rebuild.

Dashboard-only data updates (new SQLite under `data/results/`) follow the same flow — the bundled DB in git is what the hosted app reads.
