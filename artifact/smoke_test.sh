#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${ROOT}/artifact-output/smoke"
mkdir -p "${OUT}"

echo "[1/5] Recording environment"
bash "${ROOT}/artifact/environment_report.sh" | tee "${OUT}/environment.txt"

echo "[2/5] Running a numerical decryption-failure case"
(
    cd "${ROOT}"
    sage -python scripts/estimate_dec_failures.py --case end-t-512
) | tee "${OUT}/estimate_dec_failures_end_t_512.txt"

build_end() {
    local dir="$1"
    local label="$2"
    echo "Building ${label}"
    (
        cd "${ROOT}/${dir}"
        make main -j"$(nproc)"
    ) >"${OUT}/build_${label}.txt" 2>&1
}

echo "[3/5] Building portable END implementations"
build_end "END KEM/END-512/END-512-257-C" "end512_c"
build_end "END KEM/END-1024/END-1024-257-C" "end1024_c"

echo "[4/5] Building FCL-ML-KEM correctness test"
(
    cd "${ROOT}/FCL-ML-KEM"
    make test/test_kyber512 -j"$(nproc)"
    ./test/test_kyber512
) | tee "${OUT}/fcl_ml_kem_512_test.txt"

echo "[5/5] AVX2 build check (when supported)"
if grep -qw avx2 /proc/cpuinfo 2>/dev/null; then
    build_end "END KEM/END-512/END-512-257-AVX2" "end512_avx2"
    build_end "END KEM/END-1024/END-1024-257-AVX2" "end1024_avx2"
else
    echo "CPU does not advertise AVX2; skipping AVX2 build." | tee "${OUT}/avx2_skipped.txt"
fi

echo "Smoke test completed. Logs: ${OUT}"
