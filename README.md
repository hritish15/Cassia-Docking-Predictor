# Cassia Docking Therapeutic Predictor

This project implements molecular docking to predict the therapeutic potential of Cassia phytochemicals. Binding affinities (kcal/mol) are estimated against curated protein targets for each disease area, then mapped to a 0-1 therapeutic score and predictive rate.

**Important:** The output is a docking-predicted research screening score. It is not a clinically proven cure probability. Estimated affinities should be replaced with real AutoDock Vina results for publication.

## Compound classes

- Alkaloids
- Flavonoids
- Saponins
- Anthraquinones
- Tannins

## Disease areas & target proteins

| Disease | Protein Targets | PDB IDs |
|---|---|---|
| Inflammation | COX-2, TNF-alpha, IL-6 | 5KIR, 2AZ5, 1ALU |
| Diabetes | alpha-glucosidase, DPP-4 | 5KZW, 1X70 |
| Antimicrobial | DHFR, CYP51 | 3SRW, 5EQB |
| Cancer Research | EGFR, PI3K, VEGFR2 | 1M17, 4L23, 2P2H |
| Antiviral | 3CLpro | 6LU7 |
| Neuroprotective | AChE | 4EY7 |

## Files

- `index.html` - browser app (loads docking results via fetch)
- `styles.css` - app styling
- `app.js` - app logic with docking-aware rendering
- `data/compounds.json` - Cassia compound list with SMILES
- `data/docking_results.json` - docking binding affinities per compound per target
- `data/proteins/` - target PDB structures (download with --prepare-targets)
- `scripts/dock_all_compounds.py` - main docking pipeline

ML model files (legacy) moved to `models/ml/`.

## Run the browser app

Open `index.html` in a browser. No server required.

## Run the Streamlit GUI

```bash
pip install streamlit rdkit-pypi
streamlit run streamlit_app.py
```

## Generate docking results (estimated mode)

```bash
python scripts/dock_all_compounds.py
```

This uses RDKit-based QSAR estimation (no Vina required). Replace with real docking for publications.

## Generate docking results (AutoDock Vina mode)

1. Install [AutoDock Vina](https://vina.scripps.edu/) and OpenBabel.
2. Download and prepare target proteins:

```bash
python scripts/dock_all_compounds.py --prepare-targets
```

3. Run docking:

```bash
python scripts/dock_all_compounds.py --vina
```

## Binding affinity → therapeutic score mapping

| Affinity (kcal/mol) | Therapeutic Score | Label |
|---|---|---|
| ≤ -9.0 | ≥ 0.80 | Very High |
| -7.5 to -9.0 | 0.65 - 0.80 | High |
| -6.5 to -7.5 | 0.50 - 0.65 | Moderate |
| -5.5 to -6.5 | 0.30 - 0.50 | Low |
| > -5.5 | < 0.30 | Very Low |

The **predictive rate** is the therapeutic score expressed as a percentage.

