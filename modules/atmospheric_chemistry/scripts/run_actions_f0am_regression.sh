#!/usr/bin/env bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../" && pwd)"
ACTIONS_DIR="${REPO_ROOT}/modules/atmospheric_chemistry/test/tests/actions"
APP="${REPO_ROOT}/modules/atmospheric_chemistry/atmospheric_chemistry-opt"
CSVDIFF="${REPO_ROOT}/python/mooseutils/csvdiff.py"

echo "=== SCRIPT_DIR: ${SCRIPT_DIR} ==="
echo "=== REPO_ROOT: ${REPO_ROOT} ==="
echo "=== ACTIONS_DIR: ${ACTIONS_DIR} ==="
echo "=== APP: ${APP} ==="

if [[ ! -x "${APP}" ]]; then
  echo "ERROR: app not found: ${APP}" >&2
  false
fi

if [[ ! -f "${CSVDIFF}" ]]; then
  echo "ERROR: csvdiff not found: ${CSVDIFF}" >&2
  false
fi

cd "${ACTIONS_DIR}"

run_case() {
  local input="$1"
  local log="$2"
  echo "=== RUN ${input} ==="
  "${APP}" -i "${input}" > "${log}" 2>&1
}

check_csv() {
  local actual="$1"
  local gold="$2"
  echo "=== DIFF ${actual} vs ${gold} ==="
  python3 "${CSVDIFF}" "${actual}" "${gold}" --relative-tolerance 5e-2 --abs-zero 1e-6 --ignore-fields ONE RO2 CH3ONO
}

run_case "vs_F0AM_chamber_S1_box.i" "vs_F0AM_chamber_S1_box.repro.log"
run_case "vs_F0AM_chamber_S2_box.i" "vs_F0AM_chamber_S2_box.repro.log"
run_case "vs_F0AM_chamber_S3_box.i" "vs_F0AM_chamber_S3_box.repro.log"
run_case "vs_F0AM_dielcycle_box.i" "vs_F0AM_dielcycle_box.repro.log"
run_case "vs_atchem2_transport_building.i" "vs_atchem2_transport_building.repro.log"

check_csv "vs_F0AM_chamber_S1_box.csv" "gold/vs_F0AM_chamber_S1_box.csv" > vs_F0AM_chamber_S1_box.csv-check.log
check_csv "vs_F0AM_chamber_S2_box.csv" "gold/vs_F0AM_chamber_S2_box.csv" > vs_F0AM_chamber_S2_box.csv-check.log
check_csv "vs_F0AM_chamber_S3_box.csv" "gold/vs_F0AM_chamber_S3_box.csv" > vs_F0AM_chamber_S3_box.csv-check.log
cd -
echo "All requested regressions completed."
