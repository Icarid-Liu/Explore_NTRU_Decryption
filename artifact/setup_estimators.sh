#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPTS="${ROOT}/scripts"

clone_at_commit() {
    local url="$1"
    local dest="$2"
    local commit="$3"

    if [[ -e "${dest}" && ! -d "${dest}/.git" ]]; then
        echo "error: ${dest} exists but is not a Git checkout" >&2
        exit 1
    fi

    if [[ ! -d "${dest}/.git" ]]; then
        git clone --no-tags "${url}" "${dest}"
    fi

    git -C "${dest}" fetch --no-tags origin "${commit}"
    git -C "${dest}" checkout --detach "${commit}"

    local actual
    actual="$(git -C "${dest}" rev-parse HEAD)"
    if [[ "${actual}" != "${commit}" ]]; then
        echo "error: expected ${commit}, got ${actual} in ${dest}" >&2
        exit 1
    fi

    echo "${dest}: ${actual}"
}

clone_at_commit \
    "https://github.com/malb/lattice-estimator.git" \
    "${SCRIPTS}/lattice-estimator" \
    "3e48ef421ec256afddb3e7d2249a77eab6e9ba12"

clone_at_commit \
    "https://github.com/identitymapping/enhanced_lattice-estimator.git" \
    "${SCRIPTS}/enhanced_lattice-estimator" \
    "876b66173f4354a96ddafc0ce3a79767ec43c6d4"

clone_at_commit \
    "https://github.com/yonghaason/PrimalMeetLWE.git" \
    "${SCRIPTS}/PrimalMeetLWE" \
    "61115115830c909e42758f2774606074bf98afb1"
