"""Fetch the precomputed dashboard SQLite files from a GitHub Release.

Both dashboard DBs are larger than GitHub's 100MB hard blob limit (and Git LFS is
disabled for this repository), so they are published as release assets instead of
being committed -- see the `.gitignore` entries for them. A Weka checkout, or any
machine that rebuilt them locally, already has them on disk and never downloads.
A fresh clone (Streamlit Community Cloud, most notably) has neither, and this
module is what lets the dashboard serve the final analysis there with no Weka or
SSH access.

Which results the dashboard serves is pinned by environment, so republishing is a
new release plus an env change rather than a code change:

  CASF_DASHBOARD_RELEASE       release tag to pull from (default below)
  CASF_DASHBOARD_RELEASE_REPO  `owner/name` holding the release (default below);
                               the app itself may be served from a mirror, while
                               the assets stay on the canonical repository

Only the two files the app actually opens are fetchable. `casf_per_ligand_long.csv`
is a rebuild input for the extended analysis, never read by the dashboard, and is
deliberately not listed here.
"""

from __future__ import annotations

import os
import urllib.request
from pathlib import Path
from typing import Callable

from casf_benchmark.paths import DEFAULT_DASHBOARD_DB, DEFAULT_EXTENDED_DB

DEFAULT_RELEASE_TAG = "dashboard-data-qwen-druglike"
DEFAULT_RELEASE_REPO = "YerevaNN/casf-benchmark"

#: The only paths this module will ever write. Keyed by location rather than by
#: bare filename so that a DB the operator pointed us at elsewhere is never
#: silently overwritten by release contents.
RELEASE_ASSET_PATHS: tuple[Path, ...] = (DEFAULT_DASHBOARD_DB, DEFAULT_EXTENDED_DB)

CHUNK_BYTES = 1 << 20
_USER_AGENT = "casf-benchmark-dashboard"

#: Called with (bytes_downloaded, total_bytes); total is 0 when unknown.
ProgressCallback = Callable[[int, int], None]


def release_tag() -> str:
    """Release tag to fetch dashboard data from."""
    return (os.environ.get("CASF_DASHBOARD_RELEASE") or "").strip() or DEFAULT_RELEASE_TAG


def release_repo() -> str:
    """`owner/name` of the repository whose release holds the dashboard data."""
    return (os.environ.get("CASF_DASHBOARD_RELEASE_REPO") or "").strip() or DEFAULT_RELEASE_REPO


def asset_url(asset_name: str, tag: str | None = None, repo: str | None = None) -> str:
    """Public download URL for one asset of the pinned release."""
    return (
        f"https://github.com/{repo or release_repo()}"
        f"/releases/download/{tag or release_tag()}/{asset_name}"
    )


def _normalize(path: Path) -> Path:
    return Path(os.path.expanduser(str(path))).resolve()


def is_release_asset(path: Path) -> bool:
    """True when `path` is one of the default DB locations published as an asset.

    Strict on purpose: a path typed into the dashboard sidebar, or set via
    `CASF_DASHBOARD_DB` / `CASF_EXTENDED_DB`, belongs to whoever set it and is
    left alone even if it happens to share a filename with a release asset.
    """
    candidate = _normalize(path)
    return any(candidate == _normalize(default) for default in RELEASE_ASSET_PATHS)


def fetch_release_asset(
    path: Path,
    *,
    tag: str | None = None,
    repo: str | None = None,
    progress: ProgressCallback | None = None,
    timeout: float = 60.0,
) -> bool:
    """Download `path` from the pinned release; return True if a download happened.

    A no-op (returning False) when the file is already on disk or when `path` is not
    one of `RELEASE_ASSET_PATHS`, which is what keeps Weka and local rebuilds working
    untouched.

    Writes to a sibling `.part` file and renames only on success, so an interrupted
    or failed fetch never leaves a truncated SQLite file behind that a later run
    would mistake for a complete download.
    """
    if path.exists():
        return False
    if not is_release_asset(path):
        return False

    url = asset_url(path.name, tag=tag, repo=repo)
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(path.name + ".part")
    request = urllib.request.Request(url, headers={"User-Agent": _USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            total = int(response.headers.get("Content-Length") or 0)
            downloaded = 0
            if progress is not None:
                progress(downloaded, total)
            with part.open("wb") as handle:
                while True:
                    chunk = response.read(CHUNK_BYTES)
                    if not chunk:
                        break
                    handle.write(chunk)
                    downloaded += len(chunk)
                    if progress is not None:
                        progress(downloaded, total)
        part.replace(path)
    except BaseException:
        part.unlink(missing_ok=True)
        raise
    return True
