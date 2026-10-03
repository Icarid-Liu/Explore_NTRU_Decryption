#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DIST="${ROOT}/artifact-dist"

if ! git -C "${ROOT}" diff --quiet || ! git -C "${ROOT}" diff --cached --quiet; then
    echo "error: tracked working tree has uncommitted changes; archive a reviewed commit only" >&2
    exit 1
fi

commit="$(git -C "${ROOT}" rev-parse HEAD)"
short="$(git -C "${ROOT}" rev-parse --short=12 HEAD)"
name="Explore_NTRU_Decryption-asiacrypt2026-${short}.tar.gz"

mkdir -p "${DIST}"
git -C "${ROOT}" archive \
    --format=tar.gz \
    --prefix="Explore_NTRU_Decryption/" \
    --output="${DIST}/${name}" \
    "${commit}"

echo "Commit: ${commit}"
echo "Archive: ${DIST}/${name}"
sha256sum "${DIST}/${name}"
