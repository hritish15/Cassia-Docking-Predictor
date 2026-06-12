import csv
import json
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import py3Dmol

import streamlit as st


ROOT = Path(__file__).resolve().parent
DOCKING_PATH = ROOT / "data" / "docking_results.json"
COMPOUNDS_PATH = ROOT / "data" / "compounds.json"
OPTIMIZATIONS_PATH = ROOT / "data" / "optimizations.csv"

DISEASE_LABELS = {
    "inflammation": "Inflammation",
    "diabetes": "Diabetes",
    "antimicrobial": "Antimicrobial",
    "cancer_research": "Cancer Research",
    "antiviral": "Antiviral",
    "neuroprotective": "Neuroprotective",
}


@st.cache_data
def load_compounds() -> list[dict]:
    return json.loads(COMPOUNDS_PATH.read_text(encoding="utf-8"))


@st.cache_data
def load_docking() -> dict:
    docking = json.loads(DOCKING_PATH.read_text(encoding="utf-8"))
    return docking


def score_color(score) -> str:
    if score is None:
        return "#a8c3b3"
    if score >= 0.8:
        return "#94d84e"
    if score >= 0.65:
        return "#4fd0b0"
    if score >= 0.5:
        return "#e3b44b"
    return "#e87171"


def affinity_color(affinity) -> str:
    if affinity is None:
        return "#a8c3b3"
    if affinity <= -9.0:
        return "#94d84e"
    if affinity <= -7.0:
        return "#4fd0b0"
    if affinity <= -6.0:
        return "#e3b44b"
    return "#e87171"


def render_docking_card(disease_key: str, prediction: dict) -> None:
    score = prediction.get("therapeutic_score")
    affinity = prediction.get("best_affinity_kcal_mol")
    rate = prediction.get("predictive_rate", "N/A")
    color = score_color(score)
    aff_color = affinity_color(affinity)
    width = 0 if score is None else int(round(float(score) * 100))

    targets_hit = prediction.get("targets_hit", [])
    hit_count = prediction.get("targets_hit_count", 0)
    target_str = ", ".join(targets_hit[:3])
    if hit_count > 3:
        target_str += f" +{hit_count - 3}"

    aff_text = f"{affinity} kcal/mol" if affinity is not None else "N/A"

    st.markdown(
        f"""
        <div class="prediction-card">
          <div class="disease-title">{DISEASE_LABELS.get(disease_key, disease_key)}</div>
          <div class="affinity" style="color: {aff_color};">{aff_text}</div>
          <div class="score" style="color: {color};">{rate}</div>
          <div class="bar"><span style="width: {width}%; background: {color};"></span></div>
          <div class="label">{prediction.get('therapeutic_label', '')}</div>
          <div class="target">Best target: {target_str}</div>
        </div>
        """,
        unsafe_allow_html=True,
    )


# ---------------------------------------------------------------------------
# Optimization analysis
# ---------------------------------------------------------------------------
PHYTO_COLS = ["ALKALOIDS", "FLAVANOIDS", "TANNINS", "TERPENOIDS",
              "STEROIDS", "PHENOLS", "SAPONINS", "ANTHRA"]

CLASS_TO_CSV = {
    "Flavonoid":           {"FLAVANOIDS": 1.0},
    "Anthraquinone":       {"ANTHRA": 1.0},
    "Tannin":              {"TANNINS": 1.0},
    "Tannin / Phenolic":   {"TANNINS": 0.5, "PHENOLS": 0.5},
    "Tannin / Flavonoid":  {"TANNINS": 0.5, "FLAVANOIDS": 0.5},
}


def affinity_to_importance(avg_affinity: float) -> float:
    raw = (avg_affinity - (-5.0)) / (-11.0 - (-5.0))
    return max(0.0, min(1.0, raw))


@st.cache_data
def load_optimization_csv() -> list[dict]:
    lines = OPTIMIZATIONS_PATH.read_text().strip().splitlines()
    rows = []
    for line in lines[1:]:
        parts = [p.strip() for p in line.split("\t")]
        if len(parts) < 2:
            continue
        values = {}
        for i, col in enumerate(PHYTO_COLS):
            raw = parts[i + 1] if i + 1 < len(parts) else ""
            values[col] = float(raw) if raw.strip() else 0.0
        rows.append({
            "plant_part": parts[0],
            "solvent": parts[-1],
            "values": values,
        })
    return rows


@st.cache_data
def compute_class_importance(docking: dict) -> dict:
    class_affinities = defaultdict(list)
    for cdata in docking.values():
        if cdata.get("docking_status") != "complete":
            continue
        cls = cdata["compound_class"]
        for dk, dd in cdata.get("disease_predictions", {}).items():
            aff = dd.get("best_affinity_kcal_mol")
            if aff is not None:
                class_affinities[(dk, cls)].append(aff)

    importance = {}
    for dk in DISEASE_LABELS:
        scores = {}
        for (ddk, cls), vals in class_affinities.items():
            if ddk != dk:
                continue
            avg = sum(vals) / len(vals)
            scores[cls] = affinity_to_importance(avg)
        importance[dk] = scores
    return importance


def render_optimization_page(docking, optimisation_rows) -> None:
    importance = compute_class_importance(docking)

    st.markdown('<div class="main-title">Extraction Optimization</div>', unsafe_allow_html=True)
    st.markdown(
        '<div class="subtitle">Best plant part & solvent per disease target based on docking results</div>',
        unsafe_allow_html=True,
    )

    for dk, dlabel in DISEASE_LABELS.items():
        weights = importance.get(dk, {})
        if not weights:
            continue

        st.markdown(f"### {dlabel}")
        cols = st.columns([1, 2])
        with cols[0]:
            st.write("**Compound class importance**")
            for cls, w in sorted(weights.items(), key=lambda x: -x[1]):
                pct = int(w * 100)
                bar_width = max(4, pct)
                st.markdown(
                    f'<div style="margin-bottom:4px;font-size:0.85rem">'
                    f'{cls}: {pct}%'
                    f'<div style="background:#10291f;border-radius:999px;height:6px;width:100%">'
                    f'<div style="background:#4fd0b0;width:{bar_width}%;height:6px;border-radius:999px"></div>'
                    f'</div></div>',
                    unsafe_allow_html=True,
                )

        with cols[1]:
            scored = []
            for row in optimisation_rows:
                score = 0.0
                covered = set()
                for cls, w in weights.items():
                    for csv_col, sw in CLASS_TO_CSV.get(cls, {}).items():
                        presence = row["values"].get(csv_col, 0.0)
                        if presence >= 0.75:
                            covered.add(csv_col)
                        score += w * sw * presence
                scored.append((score, row["plant_part"], row["solvent"], covered))

            scored.sort(key=lambda x: -x[0])
            top = scored[:5]

            st.write("**Top extraction methods**")
            st.markdown(
                f'<div style="display:grid;grid-template-columns:60px 100px 130px 1fr;gap:4px;font-size:0.8rem;'
                f'color:#a8c3b3;margin-bottom:4px">'
                f'<div>Score</div><div>Part</div><div>Solvent</div><div>Classes</div></div>',
                unsafe_allow_html=True,
            )
            for score, part, solvent, covered in top:
                covered_s = ", ".join(sorted(covered)) if covered else "—"
                pct = int(score / max(s[0] for s in scored) * 100) if scored else 0
                st.markdown(
                    f'<div style="display:grid;grid-template-columns:60px 100px 130px 1fr;gap:4px;font-size:0.85rem;'
                    f'align-items:center;margin-bottom:2px">'
                    f'<div style="color:#4fd0b0;font-weight:700">{score:.2f}</div>'
                    f'<div>{part}</div>'
                    f'<div>{solvent}</div>'
                    f'<div style="font-size:0.8rem;color:#a8c3b3">{covered_s}</div>'
                    f'</div>',
                    unsafe_allow_html=True,
                )
        st.divider()

    # Per-phytochemical-class extraction table
    st.markdown("### Best Extraction per Phytochemical Class")
    for col in PHYTO_COLS:
        best = []
        for row in optimisation_rows:
            val = row["values"].get(col, 0.0)
            best.append((val, row["plant_part"], row["solvent"]))
        best.sort(key=lambda x: -x[0])
        top = best[:3]
        items = []
        for val, part, solvent in top:
            label = "Present" if val >= 1.0 else ("Mild" if val >= 0.75 else "No")
            items.append(f"{part} + {solvent} → {label}")
        st.write(f"**{col}:** {' | '.join(items)}")


# ---------------------------------------------------------------------------
# 3D Docking Visualisation
# ---------------------------------------------------------------------------
DISEASE_TARGETS = {
    "inflammation": {"targets": [
        {"name": "COX-2",      "pdb_id": "5KIR", "center": [23.5, 20.8, 19.1], "size": [22, 22, 22]},
        {"name": "TNF-alpha",  "pdb_id": "2AZ5", "center": [-13.9, 71.4, 26.9], "size": [26, 26, 26]},
        {"name": "IL-6",       "pdb_id": "1ALU", "center": [-5.0, -12.0, 1.0], "size": [26, 26, 26]},
    ]},
    "diabetes": {"targets": [
        {"name": "alpha-glucosidase", "pdb_id": "5KZW", "center": [-4.4, -16.1, -13.5], "size": [30, 30, 30]},
        {"name": "DPP-4",            "pdb_id": "1X70", "center": [17.5, 57.7, 35.1], "size": [30, 30, 30]},
    ]},
    "antimicrobial": {"targets": [
        {"name": "DHFR",  "pdb_id": "3SRW", "center": [-5.1, -32.1, 6.1], "size": [22, 22, 22]},
        {"name": "CYP51", "pdb_id": "5EQB", "center": [16.9, 11.7, 17.8], "size": [22, 22, 22]},
    ]},
    "cancer_research": {"targets": [
        {"name": "EGFR",   "pdb_id": "1M17", "center": [22.0, 0.3, 52.8], "size": [22, 22, 22]},
        {"name": "PI3K",   "pdb_id": "4L23", "center": [31.8, 45.8, 42.2], "size": [26, 26, 26]},
        {"name": "VEGFR2", "pdb_id": "2P2H", "center": [7.6, 19.1, 45.1], "size": [22, 22, 22]},
    ]},
    "antiviral": {"targets": [
        {"name": "3CLpro", "pdb_id": "6LU7", "center": [-11.6, 14.6, 65.2], "size": [22, 22, 22]},
    ]},
    "neuroprotective": {"targets": [
        {"name": "AChE", "pdb_id": "4EY7", "center": [-1.6, -50.2, 2.1], "size": [26, 26, 26]},
    ]},
}

VINA_BIN = ROOT / "scripts" / "vina"
PROTEINS_DIR = ROOT / "data" / "proteins"

VINA_TIMEOUT = 300


def smiles_to_pdb_block(smiles: str) -> str | None:
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    result = AllChem.EmbedMolecule(mol, params)
    if result != 0:
        params.useRandomCoords = True
        result = AllChem.EmbedMolecule(mol, params)
    if result != 0:
        return None
    AllChem.MMFFOptimizeMolecule(mol)
    return Chem.MolToPDBBlock(mol)


def run_vina_pair(smiles: str, disease_key: str, target_info: dict) -> str | None:
    target_dir = PROTEINS_DIR / disease_key
    receptor_pdbqt = target_dir / f"{target_info['name']}.pdbqt"
    receptor_pdb = target_dir / f"{target_info['pdb_id']}.pdb"
    if not receptor_pdbqt.exists() or not receptor_pdb.exists():
        return None

    pdb_block = smiles_to_pdb_block(smiles)
    if pdb_block is None:
        return None

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        ligand_pdb = tmp_path / "ligand.pdb"
        ligand_pdbqt = tmp_path / "ligand.pdbqt"
        out_pdbqt = tmp_path / "out.pdbqt"

        ligand_pdb.write_text(pdb_block)

        r = subprocess.run(
            ["obabel", str(ligand_pdb), "-O", str(ligand_pdbqt), "--gen3D"],
            capture_output=True, text=True, timeout=60,
        )
        if r.returncode != 0:
            return None

        center = target_info["center"]
        size = target_info["size"]
        r = subprocess.run(
            [
                str(VINA_BIN), "--receptor", str(receptor_pdbqt),
                "--ligand", str(ligand_pdbqt),
                "--out", str(out_pdbqt),
                "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
                "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
                "--exhaustiveness", "2",
            ],
            capture_output=True, text=True, timeout=VINA_TIMEOUT,
        )
        if r.returncode != 0:
            return None

        # Read the best pose
        return out_pdbqt.read_text()

    return None


def show_3d_docking(compound_name: str, smiles: str, disease_key: str, target_info: dict) -> None:
    # Run Vina
    with st.spinner("Running docking for 3D view (~30s)..."):
        try:
            vina_output = run_vina_pair(smiles, disease_key, target_info)
        except Exception:
            vina_output = None

    if vina_output is None:
        st.error("Docking failed. Ensure Vina is available and the ligand can be processed.")
        return

    # Convert Vina output PDBQT → PDB
    docked_pdb = ""
    try:
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            pdbqt_f = tmp / "out.pdbqt"
            pdb_f = tmp / "out.pdb"
            pdbqt_f.write_text(vina_output)
            subprocess.run(["obabel", str(pdbqt_f), "-O", str(pdb_f)],
                           capture_output=True, text=True, timeout=30)
            if pdb_f.exists():
                docked_pdb = pdb_f.read_text()
    except Exception:
        pass

    if not docked_pdb:
        st.error("Could not convert docked pose for 3D view.")
        return

    # Generate HTML for the three views
    pdb_id = target_info["pdb_id"]
    complex_html = None
    receptor_html = None
    ligand_html = None

    st.markdown(f"#### {target_info['name']} — Docking Pose (exhaustiveness=2)")

    # Complex view: docked ligand (PDB) + receptor fetched from RCSB
    try:
        v = py3Dmol.view(width=500, height=400, viewergrid=(1, 1))
        v.addModel(docked_pdb, "pdb")
        v.setStyle({"stick": {"color": "#e3b44b"}})
        v.zoomTo()
        complex_html = v._make_html()
    except Exception:
        pass

    # Receptor view: fetch via PDB ID from RCSB
    try:
        v = py3Dmol.view(query=f"pdb:{pdb_id}", width=400, height=350)
        v.setStyle({"cartoon": {"color": "spectrum"}})
        v.zoomTo()
        receptor_html = v._make_html()
    except Exception:
        pass

    # Ligand-only view: from RDKit
    lig_pdb = smiles_to_pdb_block(smiles)
    if lig_pdb:
        try:
            v = py3Dmol.view(width=400, height=350, viewergrid=(1, 1))
            v.addModel(lig_pdb, "pdb")
            v.setStyle({"stick": {"colorscheme": "Jmol"}})
            v.zoomTo()
            ligand_html = v._make_html()
        except Exception:
            pass

    # Display the three viewers
    st.markdown("**Docked Ligand Pose**")
    if complex_html:
        st.components.v1.html(complex_html, height=420)
    else:
        st.caption("Not available")

    col1, col2 = st.columns(2)
    with col1:
        st.markdown(f"**Receptor (RCSB: {pdb_id})**")
        if receptor_html:
            st.components.v1.html(receptor_html, height=370)
        else:
            st.caption("Not available")
    with col2:
        st.markdown(f"**Ligand: {compound_name}**")
        if ligand_html:
            st.components.v1.html(ligand_html, height=370)
        else:
            st.caption("Not available")


def main() -> None:
    st.set_page_config(
        page_title="Cassia Docking Predictor",
        page_icon="AI",
        layout="wide",
    )

    st.markdown(
        """
        <style>
          .stApp {
            background: radial-gradient(circle at top left, #153d28, #07130d 50%);
            color: #edf8f0;
          }
          section[data-testid="stSidebar"] {
            background: #0c1f17;
            border-right: 1px solid #24543f;
          }
          .main-title {
            font-size: 2.4rem;
            font-weight: 800;
            margin-bottom: 0.2rem;
          }
          .subtitle {
            color: #a8c3b3;
            margin-bottom: 1.4rem;
          }
          .notice {
            border-left: 4px solid #e3b44b;
            background: rgba(227, 180, 75, 0.11);
            padding: 0.85rem 1rem;
            border-radius: 6px;
            color: #f4d99a;
            margin: 1rem 0 1.4rem;
          }
          .smiles-box {
            border: 1px solid #24543f;
            border-radius: 8px;
            background: #06140e;
            padding: 1rem;
            overflow-wrap: anywhere;
            font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
            color: #d7ffe8;
          }
          .prediction-card {
            border: 1px solid #24543f;
            border-radius: 8px;
            background: rgba(12, 31, 23, 0.9);
            padding: 1rem;
            min-height: 200px;
          }
          .disease-title {
            font-weight: 700;
            color: #edf8f0;
            margin-bottom: 0.3rem;
          }
          .affinity {
            font-size: 1rem;
            line-height: 1.2;
            font-weight: 600;
            margin-bottom: 0.2rem;
          }
          .score {
            font-size: 2.1rem;
            line-height: 1;
            font-weight: 800;
            margin-bottom: 0.65rem;
          }
          .bar {
            height: 8px;
            background: #10291f;
            border-radius: 999px;
            overflow: hidden;
            margin-bottom: 0.75rem;
          }
          .bar span {
            display: block;
            height: 100%;
            border-radius: 999px;
          }
          .label {
            color: #edf8f0;
            font-weight: 700;
            margin-bottom: 0.35rem;
          }
          .target {
            color: #a8c3b3;
            line-height: 1.35;
            font-size: 0.85rem;
          }
          div[data-testid="column"] {
            padding: 0 6px;
          }
        </style>
        """,
        unsafe_allow_html=True,
    )

    if not DOCKING_PATH.exists():
        st.error("Docking results not found. Run: python scripts/dock_all_compounds.py")
        return

    compounds = load_compounds()
    docking = load_docking()
    optimisation_rows = load_optimization_csv()

    st.sidebar.title("Cassia Docking Predictor")
    view = st.sidebar.radio("View", ["Compound Docking", "Extraction Optimization"])

    if view == "Extraction Optimization":
        render_optimization_page(docking, optimisation_rows)
        return

    def compound_id(name: str) -> str:
        return name.lower().replace(" ", "-").replace("/", "-")

    classes = sorted({c["compound_class"] for c in compounds})
    selected_class = st.sidebar.selectbox("Compound class", classes)
    class_cmpds = [c for c in compounds if c["compound_class"] == selected_class]
    selected_name = st.sidebar.selectbox("Compound", [c["compound_name"] for c in class_cmpds])
    compound = next(c for c in class_cmpds if c["compound_name"] == selected_name)
    cid = compound_id(compound["compound_name"])
    output = docking.get(cid, {})

    st.markdown('<div class="main-title">Cassia Docking Therapeutic Predictor</div>', unsafe_allow_html=True)
    st.markdown(
        '<div class="subtitle">Molecular Docking (AutoDock Vina / Estimated)</div>',
        unsafe_allow_html=True,
    )

    docking_status = output.get("docking_status", "")

    if docking_status == "skipped_pending_structure":
        st.markdown(
            '<div class="notice">This compound has a pending structure. '
            'Docking will be available once the SMILES is confirmed experimentally.</div>',
            unsafe_allow_html=True,
        )
    else:
        st.markdown(
            '<div class="notice">Research-screening prediction based on molecular docking (binding affinity kcal/mol). '
            'Not a clinically proven cure probability. Replace estimated affinities with real Vina results: '
            '<code>python scripts/dock_all_compounds.py --vina</code></div>',
            unsafe_allow_html=True,
        )

    left, right = st.columns([1.2, 1])
    with left:
        st.subheader(compound["compound_name"])
        st.write(f"**Class:** {compound.get('compound_class', '')}")
        st.write(f"**Source:** {compound.get('plant_source', '')}")
        st.write(f"**Known activity:** {compound.get('known_activity', '')}")
        st.write(f"**Target proteins:** {compound.get('target_protein', '')}")
    with right:
        st.write("**Predictor info**")
        st.write("Method: `Molecular Docking`")
        st.write(f"Compound status: `{docking_status}`")
        disease_count = len(output.get("disease_predictions", {}))
        st.write(f"Disease areas assessed: `{disease_count}`")

    st.write("**SMILES input**")
    st.markdown(f'<div class="smiles-box">{compound.get("SMILES", "")}</div>', unsafe_allow_html=True)

    if docking_status == "complete":
        st.subheader("Docking-Based Therapeutic Potential")

        predictions_data = output.get("disease_predictions", {})
        items = list(predictions_data.items())
        cols_per_row = 3
        for start in range(0, len(items), cols_per_row):
            columns = st.columns(cols_per_row)
            for column, (disease_key, prediction) in zip(columns, items[start: start + cols_per_row]):
                with column:
                    render_docking_card(disease_key, prediction)

        with st.expander("Per-target binding details"):
            for disease_key, prediction in predictions_data.items():
                st.write(f"**{DISEASE_LABELS.get(disease_key, disease_key)}**")
                for t in prediction.get("all_targets", []):
                    aff = t.get("binding_affinity_kcal_mol", "N/A")
                    tname = t.get("name", "")
                    st.write(f"  {tname} (PDB: {t.get('pdb_id', '')}): {aff} kcal/mol — {t.get('binding_label', '')}")

                    if docking_status == "complete" and VINA_BIN.exists():
                        # Find matching target info for 3D
                        matched = None
                        for dt in DISEASE_TARGETS.get(disease_key, {}).get("targets", []):
                            if dt["name"] == tname:
                                matched = dt
                                break
                        if matched:
                            btn_key = f"btn_3d_{cid}_{tname}"
                            if st.button(f"View 3D: {tname}", key=btn_key):
                                show_3d_docking(
                                    compound["compound_name"],
                                    compound.get("SMILES", ""),
                                    disease_key,
                                    matched,
                                )


if __name__ == "__main__":
    main()
