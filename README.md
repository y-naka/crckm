# crckm
Calculate Reaction Coverage in KEGG MODULE

## How to use
1. Install all dependencies in CRCKM.
```bash
pip install -r requirements.txt
```
2. Download reference KEGG MODULE DEFINITION.
```bash
python3 src/download.py
```
3. Format KEGG MODULE DEFINITION and calculate Reaction Coverage.
```bash
python3 src/format_and_calculation.py $KO_MATRIX
```