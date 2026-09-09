from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from casf_benchmark.catalog import load_generation_runs
from casf_benchmark.paths import DEFAULT_QWEN_GENERATION_RUNS_CONFIG


def write_config(tmp_path: Path, body: str) -> Path:
    path = tmp_path / "runs.yaml"
    path.write_text(textwrap.dedent(body))
    return path


def test_bundled_config_covers_every_cohort() -> None:
    assert DEFAULT_QWEN_GENERATION_RUNS_CONFIG.is_file()
    for cohort in ("casf16_core", "casf16_ref", "druglike"):
        runs = load_generation_runs(cohort)
        assert runs, f"no runs configured for {cohort}"
        for label, dirname in runs:
            assert label.startswith("qwen_")
            assert dirname


def test_bundled_config_labels_are_unique_per_cohort() -> None:
    for cohort in ("casf16_core", "casf16_ref", "druglike"):
        labels = [label for label, _ in load_generation_runs(cohort)]
        assert len(labels) == len(set(labels))


def test_casf16_cohorts_are_a_subset_of_druglike() -> None:
    # Every checkpoint benchmarked on CASF16 was also run on the druglike set, so the
    # dashboard can join the two; the reverse does not hold.
    druglike = {label for label, _ in load_generation_runs("druglike")}
    for cohort in ("casf16_core", "casf16_ref"):
        assert {label for label, _ in load_generation_runs(cohort)} <= druglike


def test_preserves_config_order(tmp_path: Path) -> None:
    config = write_config(
        tmp_path,
        """
        runs:
          - label: b_second
            cohorts: {druglike: dir_b}
          - label: a_first
            cohorts: {druglike: dir_a}
        """,
    )
    assert load_generation_runs("druglike", config) == (
        ("b_second", "dir_b"),
        ("a_first", "dir_a"),
    )


def test_omits_checkpoints_without_the_requested_cohort(tmp_path: Path) -> None:
    config = write_config(
        tmp_path,
        """
        runs:
          - label: both
            cohorts: {casf16_core: core_dir, druglike: druglike_dir}
          - label: druglike_only
            cohorts: {druglike: only_dir}
        """,
    )
    assert load_generation_runs("casf16_core", config) == (("both", "core_dir"),)
    assert load_generation_runs("unknown_cohort", config) == ()


def test_rejects_duplicate_labels(tmp_path: Path) -> None:
    config = write_config(
        tmp_path,
        """
        runs:
          - label: dupe
            cohorts: {druglike: first}
          - label: dupe
            cohorts: {druglike: second}
        """,
    )
    with pytest.raises(ValueError, match="duplicate label"):
        load_generation_runs("druglike", config)


def test_rejects_non_mapping_config(tmp_path: Path) -> None:
    config = write_config(tmp_path, "- not\n- a\n- mapping\n")
    with pytest.raises(ValueError, match="must be a mapping"):
        load_generation_runs("druglike", config)
