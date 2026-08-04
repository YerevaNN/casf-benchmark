#!/usr/bin/env python3
"""Build a wide CASF analysis master CSV from analyzer outputs."""

from __future__ import annotations

import argparse
from pathlib import Path

from casf_benchmark.analysis.dataset import (
    DEFAULT_GLOBAL_LONG_OUTPUT,
    build_comparison_master,
    build_global_per_ligand_long,
    expand_sources,
)

from casf_benchmark.paths import DEFAULT_MASTER_CSV, DEFAULT_PER_LIGAND_LONG_CSV

DEFAULT_CONFIG = Path(__file__).resolve().parents[1] / "config" / "casf_analysis_sources.yaml"
DEFAULT_OUTPUT = DEFAULT_MASTER_CSV


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output-csv", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--per-ligand-long-csv", type=Path, default=DEFAULT_GLOBAL_LONG_OUTPUT)
    args = parser.parse_args()
    config = args.config.resolve()
    global_long = build_global_per_ligand_long(config, args.per_ligand_long_csv.resolve())
    master = build_comparison_master(global_long)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    master.to_csv(args.output_csv.resolve(), index=False)
    print(f"Per-ligand long: {len(global_long)} row(s): {args.per_ligand_long_csv}")
    print(f"Comparison master: {len(master)} row(s) from {len(expand_sources(config))} configured run(s): {args.output_csv}")


if __name__ == "__main__":
    main()
