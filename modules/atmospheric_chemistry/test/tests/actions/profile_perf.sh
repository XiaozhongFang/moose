#!/bin/bash
set -euo pipefail

# Run PETSc TS performance sweep for chamber S1 (10800 s, dt=100).
cd modules/atmospheric_chemistry/test/tests/actions

APP=../../../../atmospheric_chemistry-opt

run_case() {
	local input_file="$1"
	local log_file="$2"
	echo "===== Running ${input_file} ====="
	${APP} -i "${input_file}" -PerformanceLog | tee "${log_file}"
	echo "----- Summary (${log_file}) -----"
	grep -E "TOTAL RUN TIME IS|TS|Step|Solve" "${log_file}" || true
}

run_case vs_F0AM_chamber_S1_box_ts_bdf_rtol1e-2.i ts_bdf_rtol1e-2.log
run_case vs_F0AM_chamber_S1_box_ts_bdf_rtol1e-3.i ts_bdf_rtol1e-3.log
run_case vs_F0AM_chamber_S1_box_ts_arkimex_rtol1e-2.i ts_arkimex_rtol1e-2.log

echo "===== TS performance sweep finished ====="
