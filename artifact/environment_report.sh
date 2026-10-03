#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

section() {
    printf '\n===== %s =====\n' "$1"
}

section "Artifact"
echo "UTC time: $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
echo "Repository: ${ROOT}"
if command -v git >/dev/null 2>&1 && git -C "${ROOT}" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    echo "Git commit: $(git -C "${ROOT}" rev-parse HEAD)"
    echo "Git branch: $(git -C "${ROOT}" rev-parse --abbrev-ref HEAD)"
    echo "Git status:"
    git -C "${ROOT}" status --short
fi

section "Operating system"
uname -a
if [[ -r /etc/os-release ]]; then
    cat /etc/os-release
fi

section "CPU"
if command -v lscpu >/dev/null 2>&1; then
    lscpu
else
    grep -m1 -E 'model name|Hardware' /proc/cpuinfo || true
fi

section "Toolchain"
for tool in gcc clang make valgrind openssl git; do
    if command -v "${tool}" >/dev/null 2>&1; then
        echo "--- ${tool} ---"
        "${tool}" --version 2>&1 | head -n 3 || true
    else
        echo "${tool}: not found"
    fi
done

section "Sage/Python"
if command -v sage >/dev/null 2>&1; then
    sage --version
    sage -python - <<'PY'
import platform
print("Python:", platform.python_version())
for name in ("mpmath", "gmpy2", "tqdm", "jupyter", "nbconvert"):
    try:
        module = __import__(name)
        print(f"{name}: {getattr(module, '__version__', 'version unavailable')}")
    except Exception as exc:
        print(f"{name}: unavailable ({exc})")
PY
else
    echo "sage: not found"
fi

section "Pinned estimator checkouts"
for d in lattice-estimator enhanced_lattice-estimator PrimalMeetLWE; do
    p="${ROOT}/scripts/${d}"
    if [[ -d "${p}/.git" ]]; then
        echo "${d}: $(git -C "${p}" rev-parse HEAD)"
    else
        echo "${d}: not fetched"
    fi
done
