#!/bin/bash
# Run a short performance profile of the box model
cd modules/atmospheric_chemistry/test/tests/actions
../../../../atmospheric_chemistry-opt -i vs_F0AM_chamber_S1_box.i --end_time 500 --time_step_interval 1000000 -PerformanceLog
