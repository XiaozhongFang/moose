# F0AM Chamber Box Model Examples

These examples reproduce the F0AM `ExampleSetup_Chamber` scenarios:
- **S1**: Isoprene (10ppb) + NO₂ (0.1ppb) + H₂O₂ (200ppb) — low NOx
- **S2**: Isoprene (10ppb) + NO₂ (1ppb) + H₂O₂ (200ppb) — mid NOx
- **S3**: Isoprene (10ppb) + NO₂ (10ppb) + H₂O₂ (200ppb) — high NOx
- **S2b**: Restart from S2 end state with JFAC×10 (light intensity ×10)

## Mechanism

All scenarios use **MCMv331_Inorg_Isoprene** (610 species, 1974 reactions) with
**BottomUp photolysis** (lab lamp spectrum integration).

## Running

```bash
# Run all scenarios
cd modules/atmospheric_chemistry
./atmospheric_chemistry-opt -i examples/box_chamber/S1_chamber.i
./atmospheric_chemistry-opt -i examples/box_chamber/S2_chamber.i
./atmospheric_chemistry-opt -i examples/box_chamber/S3_chamber.i
./atmospheric_chemistry-opt -i examples/box_chamber/S2b_chamber.i
```

## Required Files

- `doc/content/modules/atmospheric_chemistry/database/MCMv331_Inorg_Isoprene.fac`
- `doc/content/modules/atmospheric_chemistry/database/photolysis/bottomup/` (CrossSections + QuantumYields + lamp flux)
