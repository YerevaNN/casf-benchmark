"""Canonical CASF generation family catalog."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

import yaml

DEFAULT_CATALOG_PATH = Path(__file__).resolve().parents[2] / "config" / "casf_generation_families.yaml"
TIERS = ("fixed", "dynamic", "chembl_count")
REFERENCE_TIER = "reference"


@dataclass(frozen=True)
class FamilySpec:
    id: str
    display_label: str
    core_root: Path | None
    ref_root: Path | None
    tiers: tuple[str, ...] = TIERS
    method_prefix: str | None = None
    generator: str = ""
    variant: str = ""


@dataclass(frozen=True)
class ParsedMethod:
    family: str
    tier: str
    variant: str
    method: str


def _path_or_none(value: object) -> Path | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text or text.lower() in {"none", "null"}:
        return None
    return Path(text)


def _load_catalog_data(catalog_path: Path = DEFAULT_CATALOG_PATH) -> dict[str, Any]:
    with catalog_path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ValueError(f"CASF catalog must be a mapping: {catalog_path}")
    return data


@lru_cache(maxsize=8)
def load_families(catalog_path: Path = DEFAULT_CATALOG_PATH) -> list[FamilySpec]:
    data = _load_catalog_data(Path(catalog_path))
    default_tiers = tuple(str(value) for value in data.get("tiers", TIERS))
    families = []
    for raw in data.get("families", []):
        if not isinstance(raw, dict):
            raise ValueError(f"CASF family entry must be a mapping: {raw!r}")
        family_id = str(raw["id"])
        families.append(
            FamilySpec(
                id=family_id,
                display_label=str(raw.get("display_label", family_id)),
                generator=str(raw.get("generator", "")),
                method_prefix=str(raw.get("method_prefix", family_id)),
                core_root=_path_or_none(raw.get("core_root")),
                ref_root=_path_or_none(raw.get("ref_root")),
                tiers=tuple(str(value) for value in raw.get("tiers", default_tiers)),
                variant=str(raw.get("variant") or (family_id.removeprefix("qwen_") if family_id.startswith("qwen_") else "")),
            )
        )
    ids = [family.id for family in families]
    if len(ids) != len(set(ids)):
        raise ValueError("CASF catalog contains duplicate family id(s)")
    return families


@lru_cache(maxsize=8)
def load_reference_config(catalog_path: Path = DEFAULT_CATALOG_PATH) -> dict[str, Any]:
    data = _load_catalog_data(Path(catalog_path))
    reference = data.get("reference", {})
    if not isinstance(reference, dict):
        raise ValueError("CASF catalog reference section must be a mapping")
    return reference


def family_by_id(family_id: str, catalog_path: Path = DEFAULT_CATALOG_PATH) -> FamilySpec:
    for family in load_families(catalog_path):
        if family.id == family_id:
            return family
    raise KeyError(f"Unknown CASF generation family: {family_id}")


def family_method_name(family_id: str, tier: str) -> str:
    if tier not in TIERS:
        raise ValueError(f"Unknown CASF generation tier {tier!r}; expected one of {TIERS}")
    return f"{family_id}_{tier}"


def parse_method_name(method: str) -> ParsedMethod:
    text = str(method)
    for tier in TIERS:
        suffix = f"_{tier}"
        if not text.endswith(suffix):
            continue
        family_id = text[: -len(suffix)]
        family = family_by_id(family_id)
        return ParsedMethod(
            family=family.id,
            tier=tier,
            variant=family.variant if family.id.startswith("qwen_") else "",
            method=text,
        )
    raise ValueError(f"Unknown or untiered CASF generation method: {method}")


def maybe_parse_method_name(method: object) -> ParsedMethod | None:
    try:
        return parse_method_name(str(method))
    except (KeyError, ValueError):
        return None


def families_for_ligand_set(ligand_set: str) -> list[FamilySpec]:
    if ligand_set not in {"core", "ref"}:
        raise ValueError(f"Unknown ligand_set {ligand_set!r}; expected 'core' or 'ref'")
    families = []
    for family in load_families():
        if ligand_set == "ref" and family.ref_root is None:
            continue
        families.append(family)
    return families


def resolve_generation_root(family_id: str, ligand_set: str) -> Path | None:
    family = family_by_id(family_id)
    if ligand_set == "core":
        return family.core_root
    if ligand_set == "ref":
        return family.ref_root
    raise ValueError(f"Unknown ligand_set {ligand_set!r}; expected 'core' or 'ref'")


def generation_method_order(ligand_set: str = "core") -> list[str]:
    return [
        family_method_name(family.id, tier)
        for family in families_for_ligand_set(ligand_set)
        for tier in family.tiers
    ]


def display_label_for_method(method: str) -> str:
    parsed = parse_method_name(method)
    family = family_by_id(parsed.family)
    return family.display_label
