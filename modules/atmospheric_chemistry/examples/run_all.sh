#!/bin/bash
# Run all F0AM box model examples
# Usage: ./run_all.sh [atmospheric_chemistry-opt path]

EXE="${1:-./atmospheric_chemistry-opt}"
DIR="$(cd "$(dirname "$0")" && pwd)"

echo "=== F0AM Box Model Examples ==="
echo "Executable: $EXE"
echo ""

# Chamber S1 (low NOx)
echo "--- S1 Chamber (NO2=0.1ppb, 610 sp, 3h) ---"
$EXE -i "$DIR/box_chamber/S1_chamber.i" --no-gui 2>&1 | tail -1
echo ""

# Chamber S2 (mid NOx)
echo "--- S2 Chamber (NO2=1ppb, 610 sp, 3h) ---"
$EXE -i "$DIR/box_chamber/S2_chamber.i" --no-gui 2>&1 | tail -1
echo ""

# Chamber S3 (high NOx)
echo "--- S3 Chamber (NO2=10ppb, 610 sp, 3h) ---"
$EXE -i "$DIR/box_chamber/S3_chamber.i" --no-gui 2>&1 | tail -1
echo ""

# Diel cycle (2908 sp)
echo "--- Diel Cycle (2908 sp, 24h) ---"
$EXE -i "$DIR/box_dielcycle/box_dielcycle.i" --no-gui 2>&1 | tail -1
echo ""

# EKMA tutorial (5 sp)
echo "--- EKMA Tutorial (5 sp) ---"
$EXE -i "$DIR/box_ekma/box_ekma.i" --no-gui 2>&1 | tail -1
echo ""

echo "=== All examples completed ==="
