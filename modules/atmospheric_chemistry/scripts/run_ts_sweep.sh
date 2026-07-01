#!/usr/bin/env bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
ACTIONS_DIR="${REPO_ROOT}/modules/atmospheric_chemistry/test/tests/actions"

APP=""
for c in \
  "${REPO_ROOT}/atmospheric_chemistry-opt" \
  "${REPO_ROOT}/modules/atmospheric_chemistry/atmospheric_chemistry-opt"; do
  if [[ -x "${c}" ]]; then
    APP="${c}"
    break
  fi
done

if [[ -z "${APP}" ]]; then
  echo "ERROR: atmospheric_chemistry-opt not found" >&2
  false
fi

cd "${ACTIONS_DIR}"

run_case() {
  local input="$1"
  local log="$2"
  echo "Running ${input}"
  "${APP}" -i "${input}" | tee "${log}"
}

run_case "vs_F0AM_chamber_S1_box_ts_bdf_rtol1e-2.i" "ts_bdf_rtol1e-2.log"
run_case "vs_F0AM_chamber_S1_box_ts_bdf_rtol1e-3.i" "ts_bdf_rtol1e-3.log"
run_case "vs_F0AM_chamber_S1_box_ts_arkimex_rtol1e-2.i" "ts_arkimex_rtol1e-2.log"
cd -
echo "Sweep complete"