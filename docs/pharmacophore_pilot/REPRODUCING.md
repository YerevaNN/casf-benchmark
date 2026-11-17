# Reproducing the PLINDER-22 pharmacophore pilot

The pilot uses Pharmit commit `49d2c0a573f6e58bca62f674b33874bdd42c2274`
for protein-interaction query extraction and candidate feature recognition.
The backend is isolated from the benchmark `.venv`; build it with:

```bash
cd /mnt/weka/mbedrosian/code/casf-benchmark
bash scripts/build_pharmit_pilot.sh
```

The build script pins the source commits for Pharmit, smina, Open Babel,
FastCGI and BitMagic, installs its dependency environment, applies documented
modern-toolchain compatibility patches, and applies the recorded opt-in patch
that can disable non-conservative indexed-search prefilters. The patch does not
change pharmacophore recognition.
The executable used for this run has SHA-256
`e4da894249f78c26e938374f5acf37006a9215792149af9c430f798bf77bd0d3`.
`data/results/pharmacophore_plinder22/backend_build_manifest.json` records the
resolved source commits, installed packages, linked libraries and compatibility
patches from the working installation.

## Qwen completion

The top-up runner uses the original model bundle, decoder contract, prompt,
`top_p_sampling1`, temperature 1.0, top-p 1.0, 1,000-token limit, molecular
identity check and CoordToken parser. It uses a new seed for each atomic batch.

```bash
mkdir -p data/results/pharmacophore_plinder22/qwen_topup
/opt/slurm/bin/sbatch --output=data/results/pharmacophore_plinder22/qwen_topup/slurm-%j.out \
  scripts/topup_pharmacophore_qwen.sbatch
```

The completed run was Slurm job `390125`. It retained 53 replacements from 59
attempts (six molecular-identity rejections), producing exactly 1,000 accepted
records for every molecule. The original cache and metadata SHA-256 values are
unchanged before and after; see `qwen_topup/topup_manifest.json`.

## Pilot stages

Activate the repository environment so its console scripts and test dependencies
are on `PATH`:

```bash
cd /mnt/weka/mbedrosian/code/casf-benchmark
source .venv/bin/activate
python -m casf_benchmark.cli.analyze_pharmacophores inventory
python -m casf_benchmark.cli.analyze_pharmacophores smoke
python -m casf_benchmark.cli.analyze_pharmacophores queries
python -m casf_benchmark.cli.analyze_pharmacophores match
python -m casf_benchmark.cli.analyze_pharmacophores posebusters
python -m casf_benchmark.cli.analyze_pharmacophores rmsd --workers 12
python -m casf_benchmark.cli.analyze_pharmacophores summarize
pytest
```

Equivalently, after the Qwen top-up:

```bash
python -m casf_benchmark.cli.analyze_pharmacophores all --rmsd-workers 12
```

The PoseBusters replay removes hydrogens and renumbers each candidate into the
coordinate-free reference-SMILES atom order before applying the common checks,
matching the archived evaluator. Failed candidates retain their original index,
consume the sampling denominator and are never replenished. The RMSD stage
recalculates deterministic row/column anchors before accepting any existing
cache; it computes only the Qwen replacement columns absent upstream.

### Feature matching after least-squares alignment

Pharmit's default indexed triplet search was rejected for final membership. In
a discarded diagnostic, its compressed triplet fingerprint lost 4,555
memberships as radii loosened. Disabling that fingerprint still left 1,462
cross-radius losses, showing that other triplet decomposition/index bounds were
not conservative for this panel. The pilot therefore uses Pharmit only to
recognize the unchanged feature definitions for every supplied conformer.

For final membership, the Python stage exhaustively enumerates all injective,
type-compatible query/candidate feature correspondences. Pair-distance checks
are used only as necessary-condition pruning; every surviving correspondence is
validated by a proper rigid Kabsch fit and a per-feature positional-radius check.
This positional implementation does not use feature vectors, matching the
frozen `vector_on=false` queries. It supports conformation-dependent Pharmit
feature counts (for example, merged nearby hydrophobes), records the recognized
feature caches and signatures, and requires candidate-level set inclusion at
0.5, 0.75 and 1.0 Å. The final run had zero violations. Pharmit's patched
indexed mode remains part of the backend smoke and reference self/invariance
checks, but is explicitly not used for final candidate membership.

The correspondence search is exhaustive, but the continuous alignment test is
not a complete rigid-feasibility solver. Kabsch minimizes the sum of squared
point errors, not the largest point error, so a feasible match can be rejected.
The audit at `expert_audit_2026-11-17/pharmacophore_methods_en_ru.html` explains
the limitation, and `expert_audit_2026-11-17/matching_limitation_check.py`
reproduces a synthetic false negative. Its impact on molecular results and
rankings is unknown; no fix or result recalculation is implied by this report.

The configuration, paths, radii, budget grid, seeds and backend pin are in
`src/casf_benchmark/config/pharmacophore_plinder22.yaml`. All generated files
are under `data/results/pharmacophore_plinder22/`; no manuscript or dashboard
artifact is read as an output destination or overwritten.

## PLINDER-ATC48 extension

The broader panel is frozen independently of generator results. It retains all
remaining PLINDER 2024-06/v2 ATC/DrugBank identities after excluding identities
in the original PLINDER shortlist and one molecule containing an element not
supported by the common method panel. Up to eight reference-valid, noncovalent
complexes from distinct PDB entries are selected by a seeded hash ordering. The
result contains 48 molecules and 379 references; 311 queries from 44 molecules
pass the same reference-only eligibility rules. See
`data/pharmacophore_atc48/freeze_manifest.json`.

Build and validate the frozen reference panel with:

```bash
source .venv/bin/activate
python scripts/build_pharmacophore_atc48.py
python -m casf_benchmark.cli.analyze_pharmacophores \
  --config src/casf_benchmark/config/pharmacophore_plinder_atc48.yaml inventory
python -m casf_benchmark.cli.analyze_pharmacophores \
  --config src/casf_benchmark/config/pharmacophore_plinder_atc48.yaml smoke
python -m casf_benchmark.cli.analyze_pharmacophores \
  --config src/casf_benchmark/config/pharmacophore_plinder_atc48.yaml queries
```

`scripts/generate_pharmacophore_baselines.py` produces the two reference-blind
ETKDG pools. Qwen and LoQI use their original generation and acceptance paths.
The Qwen source is the merged accepted pool (47,943 conformers; 992--1,000 per
molecule), so requested single-library budgets are capped at availability. The
50:50 endpoint needs only 500 Qwen conformers per molecule and is unaffected by
this shortfall. LoQI uses deterministic CPU shards because the GPU queue was
unavailable. Thirty-eight molecules completed as four 250-conformer shards. Ten
molecules (zero-based indices 3, 4, 7, 8, 11, 13, 22, 23, 24, and 38) used eight
125-conformer shards to stay within wall-time limits. All shards use seed
`20261117 + 10000 * molecule_index + shard_index` and are merged in fixed shard
order. The mixed shard counts are explicit in the filenames and merged metadata.
The following clean submission reproduces the final layout without launching the
superseded four-shard tasks for those ten molecules:

```bash
LOQI_SHARD_COUNT=4 sbatch \
  --time=03:00:00 \
  --array=0-11,20-27,36-43,48-51,56-87,100-151,156-191%128 \
  --output=data/results/pharmacophore_atc48/generation/loqi/cpu-shard-%A_%a.out \
  scripts/run_pharmacophore_loqi_cpu_shard.sbatch
LOQI_SHARD_COUNT=8 sbatch \
  --time=02:00:00 \
  --array=24-39,56-71,88-95,104-111,176-199,304-311%80 \
  --output=data/results/pharmacophore_atc48/generation/loqi/cpu8-shard-%A_%a.out \
  scripts/run_pharmacophore_loqi_cpu_shard.sbatch
source /mnt/weka/mbedrosian/codex_dir/loqi/codex/env.sh
python scripts/run_off_the_shelf_druglike.py \
  --model loqi \
  --druglike-pickle data/pharmacophore_atc48/generation_input.pickle \
  --output-dir data/results/pharmacophore_atc48/generation/loqi \
  --merge-parts
python scripts/audit_pharmacophore_loqi_generation.py \
  --generation-input data/pharmacophore_atc48/generation_input.pickle \
  --generation-dir data/results/pharmacophore_atc48/generation/loqi \
  --output data/results/pharmacophore_atc48/generation/loqi/run_manifest.json \
  --base-seed 20261117
```

The array task is `molecule_index * shard_count + shard_index`. If a shard is
interrupted, rerun that exact task with the same shard count; complete outputs are
written atomically. In the executed run, interrupted four-shard tasks for molecule
21 (shard 1), molecule 33 (all shards, after undersized-memory and wall-time
attempts), and molecule 46 (shard 3, after reaching the original wall-time limit)
were rerun with the same task IDs and seeds. Do not mix four- and eight-shard
files for one molecule. Five unsharded files from the abandoned monolithic attempt
remain as ignored artifacts; complete shard sets take precedence, and the run
manifest records their hashes.

Generation artifacts remain separate from query and matching outputs. Once all
four method pools exist, run:

```bash
python -m casf_benchmark.cli.analyze_pharmacophores \
  --config src/casf_benchmark/config/pharmacophore_plinder_atc48.yaml match
python -m casf_benchmark.cli.analyze_pharmacophores \
  --config src/casf_benchmark/config/pharmacophore_plinder_atc48.yaml posebusters
python -m casf_benchmark.cli.analyze_pharmacophores \
  --config src/casf_benchmark/config/pharmacophore_plinder_atc48.yaml rmsd --workers 80
python -m casf_benchmark.cli.analyze_pharmacophores \
  --config src/casf_benchmark/config/pharmacophore_plinder_atc48.yaml summarize
```

The exact-identity training audit is reproducible with
`scripts/audit_pharmacophore_training_overlap.py`. It checks isomeric InChIKey
identity only and must not be interpreted as excluding scaffold or near-neighbor
overlap. The extension broadens structural coverage but is not an unseen-molecule
generalization benchmark.
