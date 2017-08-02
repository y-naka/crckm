# crckm
Calculate Reaction Coverage in KEGG MODULE

## How to use
1. Download reference KEGG MODULE DEFINITION.
```bash
python3 src/download.py
```
2. Format KEGG MODULE DEFINITION and calculate Reaction Coverage.
```bash
python3 src/format_and_calculation.py $KO_MATRIX
```