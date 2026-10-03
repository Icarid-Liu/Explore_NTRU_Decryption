#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${ROOT}/artifact-output/full"
mkdir -p "${OUT}"

bash "${ROOT}/artifact/environment_report.sh" | tee "${OUT}/environment.txt"
bash "${ROOT}/artifact/setup_estimators.sh" | tee "${OUT}/setup_estimators.txt"

(
    cd "${ROOT}"
    sage -python scripts/estimate_dec_failures.py --case all
) | tee "${OUT}/estimate_dec_failures.txt"

(
    cd "${ROOT}"
    sage -python artifact/reproduce_dfr.py --case all
) | tee "${OUT}/dfr_frontier.txt"

bash "${ROOT}/artifact/run_notebooks.sh" all | tee "${OUT}/notebooks.txt"

run_end() {
    local dir="$1"
    local label="$2"
    (
        cd "${ROOT}/${dir}"
        make main -j"$(nproc)"
        ./main
    ) | tee "${OUT}/${label}.txt"
}

run_end "END KEM/END-512/END-512-257-C" "end512_c"
run_end "END KEM/END-1024/END-1024-257-C" "end1024_c"

if grep -qw avx2 /proc/cpuinfo 2>/dev/null; then
    run_end "END KEM/END-512/END-512-257-AVX2" "end512_avx2"
    run_end "END KEM/END-1024/END-1024-257-AVX2" "end1024_avx2"
else
    echo "CPU does not advertise AVX2; AVX2 runs skipped." | tee "${OUT}/avx2_skipped.txt"
fi

(
    cd "${ROOT}/FCL-ML-KEM"
    make test/test_kyber512 nistkat/PQCgenKAT_kem512 -j"$(nproc)"
    ./test/test_kyber512
    ./nistkat/PQCgenKAT_kem512
) | tee "${OUT}/fcl_ml_kem_512_correctness_kat.txt"

echo
echo "Full reproduction workflow completed."
echo "Outputs: ${ROOT}/artifact-output/"
echo "For FCL-ML-KEM cycle counts, run the speed target documented in README.md on the benchmark machine."
