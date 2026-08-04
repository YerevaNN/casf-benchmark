#!/usr/bin/env bash
# Push to YerevaNN (origin) and the personal mirror (menuab) in one step.
#
# Usage:
#   ./scripts/push_with_mirror.sh              # push current branch
#   ./scripts/push_with_mirror.sh main         # push explicit branch
#
# The personal mirror triggers Streamlit Community Cloud redeploy.

set -euo pipefail

BRANCH="${1:-$(git branch --show-current)}"
PERSONAL_REMOTE="${PERSONAL_REMOTE:-personal}"

if ! git remote get-url "$PERSONAL_REMOTE" &>/dev/null; then
  echo "Missing remote '$PERSONAL_REMOTE'. Add it once:" >&2
  echo "  git remote add personal git@github.com:menuab/casf-benchmark.git" >&2
  exit 1
fi

echo "→ Pushing $BRANCH to origin (YerevaNN)..."
git push origin "$BRANCH"

echo "→ Pushing $BRANCH to $PERSONAL_REMOTE (menuab mirror)..."
git push "$PERSONAL_REMOTE" "$BRANCH"

echo "Done. Streamlit Cloud will redeploy from menuab/casf-benchmark when configured."
