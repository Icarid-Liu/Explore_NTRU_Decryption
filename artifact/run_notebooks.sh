#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPTS="${ROOT}/scripts"
OUT="${ROOT}/artifact-output/notebooks"
MODE="${1:-all}"

mkdir -p "${OUT}"

if ! command -v sage >/dev/null 2>&1; then
    echo "error: SageMath is required" >&2
    exit 1
fi

run_nb() {
    local input="$1"
    local output="$2"
    (
        cd "${SCRIPTS}"
        export PYTHONPATH="${SCRIPTS}/lattice-estimator:${PYTHONPATH:-}"
        sage -python -m jupyter nbconvert \
            --to notebook \
            --execute "${input}" \
            --ExecutePreprocessor.timeout=-1 \
            --output "${output}" \
            --output-dir "${OUT}"
    )
}

case "${MODE}" in
    dependencies)
        run_nb "dependencies_impact_on_wrap_errors.ipynb" "dependencies_impact_on_wrap_errors.executed.ipynb"
        ;;
    security)
        for d in lattice-estimator enhanced_lattice-estimator PrimalMeetLWE; do
            if [[ ! -d "${SCRIPTS}/${d}/.git" ]]; then
                echo "error: missing scripts/${d}; run: bash artifact/setup_estimators.sh" >&2
                exit 1
            fi
        done
        run_nb "concrete_security_estimates.ipynb" "concrete_security_estimates.executed.ipynb"
        ;;
    all)
        run_nb "dependencies_impact_on_wrap_errors.ipynb" "dependencies_impact_on_wrap_errors.executed.ipynb"
        for d in lattice-estimator enhanced_lattice-estimator PrimalMeetLWE; do
            if [[ ! -d "${SCRIPTS}/${d}/.git" ]]; then
                echo "error: missing scripts/${d}; run: bash artifact/setup_estimators.sh" >&2
                exit 1
            fi
        done
        run_nb "concrete_security_estimates.ipynb" "concrete_security_estimates.executed.ipynb"
        ;;
    *)
        echo "usage: $0 [dependencies|security|all]" >&2
        exit 2
        ;;
esac

echo "Executed notebook(s) written to ${OUT}"
