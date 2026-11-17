PHARMACOPHORE QUERY EXPERIMENT

This directory documents the PLINDER22 pilot and its ATC48 extension. The code
tests fixed conformer libraries against pharmacophore queries extracted from
experimental protein-ligand complexes. It does not train a model or run docking.

CURRENT EXPERIMENT
src/casf_benchmark/config/dataset_recipe_experiment.yaml is the active protocol:
Qwen, LoQI, their 75/25, 50/50 and 25/75 mixtures, and RDKit raw/minimized controls.
FlowR is excluded. dataset_recipe_experiment_flowr_v1.yaml is retained only as
the historical configuration for the superseded study, not the active protocol.
The original PLINDER22 pilot config also lists the broader initial method panel.

CODE
src/casf_benchmark/analysis/pharmacophore.py:
    Protein context, feature extraction, rigid matching and coverage statistics.
src/casf_benchmark/cli/analyze_pharmacophores.py:
    Inventory, backend checks, queries, matching, PoseBusters and summaries.
src/casf_benchmark/analysis/dataset_recipes.py:
    Fixed-budget allocations, sampling probabilities and diversity statistics.
src/casf_benchmark/cli/analyze_dataset_recipes.py:
    Freeze the protocol, prepare source pools, calculate diversity and summarize.
src/casf_benchmark/cli/report_dataset_recipes.py:
    Validate completed results and generate a report.
scripts/build_pharmacophore_atc48.py:
    Select ATC48 and download its reference ligand/complex structures.
scripts/build_pharmit_pilot.sh and scripts/pharmit_exact_search.patch:
    Build the pinned Pharmit feature backend and its indexed-search patch.
scripts/generate_pharmacophore_baselines.py:
    Generate paired RDKit ETKDG raw and MMFF94s-minimized pools.
scripts/topup_pharmacophore_qwen.py:
    Complete Qwen pools using the original identity/parsing acceptance rules.
scripts/run_off_the_shelf_druglike.py:
    Existing model adapter, extended with deterministic conformer sharding.
    LoQI launch and audit scripts record the shard layout and random seeds.

ENVIRONMENT AND INPUTS
Activate .venv before running tests or entry points. Besides the existing dev
dependencies, the complete pipeline needs SciPy, Gemmi, requests and a pandas
Parquet engine such as pyarrow. The recorded PoseBusters version is 0.6.5.
Pharmit is a separate pinned installation; Qwen and LoQI require their original
model environments, adapters and checkpoints. Slurm is needed only for the
provided cluster launchers, not for unit tests.

The configs retain the original cluster paths and input hashes. Adapt paths in
COPIES of those configs for another installation; do not alter the saved frozen
protocol in place. You must supply the audited PLINDER references and source
mmCIF/SDF files, the historical ATC candidate/index tables, and the conformer
libraries or generation environments. The upstream ATC annotation-join script
and explicit ATC/DrugBank ID mapping are still unavailable.

The branch contains code, configs, tests and selected reports. It intentionally
does not contain checkpoints, source datasets, generated conformer pools,
prepared match matrices, large result caches, or the complete review bundle.
References in reports to the review bundle refer to that separate local artifact.
The English/Russian HTML report itself is self-contained.

STARTING POINTS
source .venv/bin/activate
python -m casf_benchmark.cli.analyze_pharmacophores --help
python -m casf_benchmark.cli.analyze_dataset_recipes --help
python -m casf_benchmark.cli.report_dataset_recipes --help
python scripts/build_pharmacophore_atc48.py --help
pytest

REPRODUCING.md describes the original pilot and extension stages. For the current
two-source study, prepare the reference queries and source conformer libraries
first, then freeze the active config in a new output directory, run prepare and
diversity for each molecule in each panel, and finally summarize both panels:

python -m casf_benchmark.cli.analyze_dataset_recipes freeze --config CONFIG.yaml
python -m casf_benchmark.cli.analyze_dataset_recipes prepare --config CONFIG.yaml --panel plinder22 --index 0
python -m casf_benchmark.cli.analyze_dataset_recipes diversity --config CONFIG.yaml --panel plinder22 --index 0
python -m casf_benchmark.cli.analyze_dataset_recipes summarize --config CONFIG.yaml --panel plinder22
python -m casf_benchmark.cli.report_dataset_recipes --config CONFIG.yaml

The prepare/diversity commands above show one molecule only. Repeat them for all
22 PLINDER22 and 48 ATC48 molecules, and summarize ATC48 too before reporting.
No jobs are submitted unless the submission helper is called with --submit.

KNOWN MATCHING LIMITATION
The search enumerates one-to-one typed feature assignments, then checks the
least-squares Kabsch alignment for each. It is NOT a complete test of whether
any rigid alignment can place all features within tolerance. A feasible match
can be missed. Its effect on the real molecule results and rankings is unknown.
No algorithm fix or recalculation is included in this commit.

Run the small synthetic diagnostic without the datasets or Pharmit:
python docs/pharmacophore_pilot/expert_audit_2026-11-17/matching_limitation_check.py

READABLE METHODS REPORT
docs/pharmacophore_pilot/expert_audit_2026-11-17/pharmacophore_methods_en_ru.html
Open it in a browser and switch English / Русский at the top. Exact parameters
and unresolved expert-review questions are included in both languages.
