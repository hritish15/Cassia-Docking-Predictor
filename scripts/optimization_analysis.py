import csv
import json
from pathlib import Path
from collections import defaultdict

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "data"

# ---------------------------------------------------------------------------
# 1. Load docking results & group by compound class per disease target
# ---------------------------------------------------------------------------
with open(DATA / "docking_results.json") as f:
    docking = json.load(f)

# compound_class mapping → csv columns (with split weights for mixed classes)
CLASS_TO_CSV = {
    "Flavonoid":       {"FLAVANOIDS": 1.0},
    "Anthraquinone":   {"ANTHRA": 1.0},
    "Tannin":          {"TANNINS": 1.0},
    "Tannin / Phenolic": {"TANNINS": 0.5, "PHENOLS": 0.5},
    "Tannin / Flavonoid": {"TANNINS": 0.5, "FLAVANOIDS": 0.5},
}

disease_labels = {
    "inflammation": "Inflammation",
    "diabetes": "Diabetes",
    "antimicrobial": "Antimicrobial",
    "cancer_research": "Cancer Research",
    "antiviral": "Antiviral",
    "neuroprotective": "Neuroprotective",
}

# Collect avg affinity per (disease, compound_class)
class_affinities = defaultdict(list)  # key: (disease_key, class_name) -> list of best_affinities

for cid, cdata in docking.items():
    if cdata["docking_status"] != "complete":
        continue
    cls = cdata["compound_class"]
    for dk, dd in cdata["disease_predictions"].items():
        if dd["best_affinity_kcal_mol"] is not None:
            class_affinities[(dk, cls)].append(dd["best_affinity_kcal_mol"])

# Map mean affinity to 0-1 importance using absolute biological thresholds.
# -5.0 kcal/mol → 0.0 (very weak / non-binder)
# -11.0 kcal/mol → 1.0 (very strong binder)
def affinity_to_importance(avg_affinity: float) -> float:
    raw = (avg_affinity - (-5.0)) / (-11.0 - (-5.0))
    return round(max(0.0, min(1.0, raw)), 3)

disease_class_importance = {}  # {disease_key: {class_name: importance}}
for dk in disease_labels:
    scores = {}
    for (ddk, cls), vals in class_affinities.items():
        if ddk != dk:
            continue
        avg = sum(vals) / len(vals)
        scores[cls] = affinity_to_importance(avg)
    disease_class_importance[dk] = scores

# ---------------------------------------------------------------------------
# 2. Load optimizations.csv
# ---------------------------------------------------------------------------
opt_rows = []
with open(DATA / "optimizations.csv", newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        plant_part = row.get("plant part", "").strip()
        # The last column is the solvent (no header)
        solvent = list(row.values())[-1].strip()
        opt_rows.append({
            "plant_part": plant_part,
            "solvent": solvent,
            "values": {
                k: float(v) if v.strip() else 0.0
                for k, v in row.items()
                if k != "plant part" and k.strip()
            }
        })

# Rename the value columns that might include the solvent column
# Actually the DictReader has keys from the header row. The header is:
# "plant part\tALKALOIDS\tFLAVANOIDS\tTANNINS\tTERPENOIDS\tSTEROIDS\tPHENOLS\tSAPONINS\tANTHRA\t(solvent - no header)"
# So last key would be something like "" or the solvent name
# Let me re-read more carefully
opt_rows = []
with open(DATA / "optimizations.csv", newline="") as f:
    content = f.read().strip().splitlines()

# Columns: 0=plant part, 1-8=phytochemicals, 9=solvent (headerless)
PHYTO_COLS = ["ALKALOIDS", "FLAVANOIDS", "TANNINS", "TERPENOIDS",
              "STEROIDS", "PHENOLS", "SAPONINS", "ANTHRA"]

for line in content[1:]:
    parts = [p.strip() for p in line.split("\t")]
    if len(parts) < 2:
        continue
    plant_part = parts[0]
    solvent = parts[-1]
    values = {}
    for i, col in enumerate(PHYTO_COLS):
        raw = parts[i + 1] if i + 1 < len(parts) else ""
        values[col] = float(raw) if raw.strip() else 0.0
    opt_rows.append({
        "plant_part": plant_part,
        "solvent": solvent,
        "values": values,
    })

# ---------------------------------------------------------------------------
# 3. Score each (plant part, solvent) per disease target
# ---------------------------------------------------------------------------
print("=" * 90)
print("EXTRACTION OPTIMIZATION RECOMMENDATIONS PER DISEASE TARGET")
print("=" * 90)

for dk, dlabel in disease_labels.items():
    weights = disease_class_importance.get(dk, {})
    if not weights:
        print(f"\n{dlabel}: No docking data available")
        continue

    print(f"\n{'─' * 90}")
    print(f"  {dlabel}")
    print(f"{'─' * 90}")
    print(f"  Compound class importance (0-1):")
    for cls, w in sorted(weights.items(), key=lambda x: -x[1]):
        print(f"    {cls:25s}  {w:.3f}")
    print()

    # Score each extraction method
    scored = []
    for row in opt_rows:
        score = 0.0
        classes_covered = []
        for cls, w in weights.items():
            csv_map = CLASS_TO_CSV.get(cls, {})
            for csv_col, split_w in csv_map.items():
                presence = row["values"].get(csv_col, 0.0)
                contribution = w * split_w * presence
                if contribution > 0 and csv_col not in classes_covered:
                    classes_covered.append(csv_col)
                score += contribution
        scored.append((score, row["plant_part"], row["solvent"], classes_covered))

    scored.sort(key=lambda x: -x[0])
    top = scored[:5]

    print(f"  Top extraction methods:")
    print(f"  {'Score':>6s}  {'Plant Part':12s}  {'Solvent':20s}  {'Classes Covered'}")
    print(f"  {'─'*6}  {'─'*12}  {'─'*20}  {'─'*30}")
    for score, part, solvent, covered in top:
        covered_s = ", ".join(sorted(set(covered))) if covered else "none"
        print(f"  {score:6.3f}  {part:12s}  {solvent:20s}  {covered_s}")
    print()

# ---------------------------------------------------------------------------
# 4. Also show best per phytochemical class directly
# ---------------------------------------------------------------------------
print("=" * 90)
print("BEST EXTRACTION PER PHYTOCHEMICAL CLASS")
print("=" * 90)

for col in PHYTO_COLS:
    scored = []
    for row in opt_rows:
        val = row["values"].get(col, 0.0)
        scored.append((val, row["plant_part"], row["solvent"]))
    scored.sort(key=lambda x: -x[0])
    top = scored[:3]
    print(f"\n  {col}:")
    for val, part, solvent in top:
        presence = "Present" if val >= 1.0 else ("Mild" if val >= 0.75 else "No")
        print(f"    {part:12s} + {solvent:20s}  → {presence} ({val})")
