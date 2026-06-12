#!/usr/bin/env python3
"""CLI for viewing docking-based therapeutic predictions."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
COMPOUNDS_JSON = ROOT / "data" / "compounds.json"
DOCKING_JSON = ROOT / "data" / "docking_results.json"
MODEL_OUTPUTS_JSON = ROOT / "data" / "model_outputs.json"


def slugify(value: str) -> str:
    return value.lower().replace(" / ", "-").replace(" ", "-")


def main() -> None:
    parser = argparse.ArgumentParser(description="View docking predictions for a Cassia compound")
    parser.add_argument("compound", help="Compound name (e.g. Quercetin)")
    parser.add_argument("--vina", action="store_true", help="Run real Vina docking (requires vina installed)")
    args = parser.parse_args()

    compounds = json.loads(COMPOUNDS_JSON.read_text(encoding="utf-8"))
    query = args.compound.lower()
    compound = next(
        (
            item
            for item in compounds
            if slugify(item["compound_name"]) == query or item["compound_name"].lower() == query
        ),
        None,
    )
    if compound is None:
        available = ", ".join(item["compound_name"] for item in compounds)
        raise SystemExit(f"Compound not found. Available: {available}")

    cid = slugify(compound["compound_name"])

    if DOCKING_JSON.exists():
        outputs = json.loads(DOCKING_JSON.read_text(encoding="utf-8"))
        result = outputs.get(cid)
        if result and result["docking_status"] == "complete":
            print(f"{compound['compound_name']} ({compound['compound_class']})")
            print(f"SMILES: {compound['SMILES']}")
            print(f"Method: Molecular Docking")
            print()
            for disease_key, prediction in result["disease_predictions"].items():
                aff = prediction["best_affinity_kcal_mol"]
                score = prediction["therapeutic_score"]
                rate = prediction["predictive_rate"]
                label = prediction["therapeutic_label"]
                targets = ", ".join(prediction["targets_hit"])
                print(f"{disease_key}:")
                print(f"  Binding affinity: {aff} kcal/mol")
                print(f"  Therapeutic score: {score} ({label})")
                print(f"  Predictive rate: {rate}")
                print(f"  Targets hit: {targets}")
                print()
            return

    if MODEL_OUTPUTS_JSON.exists():
        outputs = json.loads(MODEL_OUTPUTS_JSON.read_text(encoding="utf-8"))
        result = outputs.get(cid)
        if result:
            print(f"{compound['compound_name']} ({compound['compound_class']})")
            print(f"SMILES: {compound['SMILES']}")
            print(f"Method: Legacy ML seed output")
            print()
            for disease, prediction in result["predictions"].items():
                print(f"- {disease}: {prediction['score']:.2f} ({prediction['label']}) target={prediction['target']}")
            return

    raise SystemExit(
        f"No predictions found for {compound['compound_name']}. "
        f"Run: python scripts/dock_all_compounds.py"
    )


if __name__ == "__main__":
    main()
