#!/usr/bin/env python3
"""
ASCO/CAP HER2 IHC & FISH Interpretation Algorithm

Implements the 2018 ASCO/CAP guidelines for HER2 status determination in breast cancer
using IHC scoring (0-3+) and FISH (fluorescence in situ hybridization) interpretation.

References:
  - Wolff AC et al. J Clin Oncol. 2018;36(20):2135-2151 (ASCO/CAP 2018 update)
  - Wolff AC et al. J Clin Oncol. 2007;25(1):118-145 (ASCO/CAP 2007 original)

Zero-dependency Python implementation (stdlib only).
License: MIT
"""

import argparse
import csv
import json
import sys
from typing import Dict, Any, List, Optional


# ---------------------------------------------------------------------------
# IHC Interpretation
# ---------------------------------------------------------------------------

def interpret_ihc(ihc_score: int, percent_staining: float = 0.0,
                  membrane_completeness: str = "complete",
                  staining_intensity: str = "strong") -> Dict[str, Any]:
    """
    Interpret HER2 IHC staining result per ASCO/CAP 2018.

    Parameters
    ----------
    ihc_score : int
        IHC score: 0, 1, 2, or 3.
    percent_staining : float
        Percentage of invasive tumor cells showing staining (0-100).
    membrane_completeness : str
        'complete' or 'incomplete' membrane staining pattern.
    staining_intensity : str
        'faint', 'weak_moderate', or 'strong' staining intensity.

    Returns
    -------
    dict with keys: ihc_score, percent_staining, her2_status, category,
                    reflex_fish_required, description, treatment_implications
    """
    if ihc_score not in (0, 1, 2, 3):
        raise ValueError(f"IHC score must be 0, 1, 2, or 3; got {ihc_score}")
    if not 0.0 <= percent_staining <= 100.0:
        raise ValueError(f"percent_staining must be 0-100; got {percent_staining}")

    result = {
        "ihc_score": ihc_score,
        "percent_staining": percent_staining,
        "membrane_completeness": membrane_completeness,
        "staining_intensity": staining_intensity,
    }

    if ihc_score == 0:
        result["her2_status"] = "Negative"
        result["category"] = "Negative"
        result["reflex_fish_required"] = False
        result["description"] = (
            "No staining or incomplete, faint/barely perceptible membrane "
            "staining in <=10% of invasive tumor cells."
        )
        result["treatment_implications"] = (
            "HER2-targeted therapy (trastuzumab, pertuzumab, T-DXd) is NOT indicated."
        )

    elif ihc_score == 1:
        result["her2_status"] = "Negative"
        result["category"] = "Negative"
        result["reflex_fish_required"] = False
        result["description"] = (
            "Faint/barely perceptible incomplete membrane staining in >10% "
            "of invasive tumor cells."
        )
        result["treatment_implications"] = (
            "HER2-targeted therapy is NOT indicated. Consider reflex FISH "
            "if strong clinical suspicion."
        )

    elif ihc_score == 2:
        result["her2_status"] = "Equivocal"
        result["category"] = "Equivocal"
        result["reflex_fish_required"] = True
        result["description"] = (
            "Weak to moderate complete membrane staining in >10% of invasive "
            "tumor cells, OR intense complete membrane staining in <=10%."
        )
        result["treatment_implications"] = (
            "Reflex FISH testing is REQUIRED to determine HER2 status. "
            "Do not initiate HER2-targeted therapy based on IHC 2+ alone."
        )

    elif ihc_score == 3:
        result["her2_status"] = "Positive"
        result["category"] = "Positive"
        result["reflex_fish_required"] = False
        result["description"] = (
            "Uniform intense complete membrane staining in >10% of invasive "
            "tumor cells."
        )
        result["treatment_implications"] = (
            "Patient is eligible for HER2-targeted therapy: trastuzumab, "
            "pertuzumab, and/or T-DXd (trastuzumab deruxtecan)."
        )

    return result


# ---------------------------------------------------------------------------
# FISH Interpretation
# ---------------------------------------------------------------------------

def interpret_fish(her2_copy_number: float, cep17_copy_number: float,
                   her2_ceph_ratio: Optional[float] = None,
                   num_cells_counted: int = 20,
                   guideline_year: int = 2018) -> Dict[str, Any]:
    """
    Interpret HER2 FISH results per ASCO/CAP 2018 guidelines.

    Parameters
    ----------
    her2_copy_number : float
        Average HER2 signals per cell.
    cep17_copy_number : float
        Average CEP17 (chromosome 17 centromere) signals per cell.
    her2_ceph_ratio : float, optional
        HER2/CEP17 ratio. If None, calculated from copy numbers.
    num_cells_counted : int
        Number of cells counted (default 20).
    guideline_year : int
        2007 or 2018 (default 2018). The 2018 update reclassifies groups 4 and 5.

    Returns
    -------
    dict with keys: her2_copy_number, cep17_copy_number, ratio, fish_status,
                    fish_group, description, treatment_implications
    """
    if her2_copy_number < 0:
        raise ValueError(f"her2_copy_number must be >= 0; got {her2_copy_number}")
    if cep17_copy_number <= 0:
        raise ValueError(f"cep17_copy_number must be > 0; got {cep17_copy_number}")
    if guideline_year not in (2007, 2018):
        raise ValueError(f"guideline_year must be 2007 or 2018; got {guideline_year}")

    if her2_ceph_ratio is None:
        her2_ceph_ratio = her2_copy_number / cep17_copy_number

    ratio = round(her2_ceph_ratio, 2)
    her2_cn = round(her2_copy_number, 1)

    result = {
        "her2_copy_number": her2_cn,
        "cep17_copy_number": round(cep17_copy_number, 1),
        "ratio": ratio,
        "num_cells_counted": num_cells_counted,
        "guideline_year": guideline_year,
    }

    # ASCO/CAP 2018 FISH groups
    if ratio >= 2.0 and her2_cn >= 4.0:
        result["fish_status"] = "Positive"
        result["fish_group"] = "Group 1"
        result["description"] = (
            f"HER2/CEP17 ratio >=2.0 ({ratio}) AND HER2 copy number >=4.0 ({her2_cn}). "
            "Positive for HER2 amplification."
        )

    elif ratio >= 2.0 and her2_cn < 4.0:
        if guideline_year == 2018:
            result["fish_status"] = "Positive"
            result["fish_group"] = "Group 4 (2018: Positive)"
            result["description"] = (
                f"HER2/CEP17 ratio >=2.0 ({ratio}) but HER2 copy number <4.0 ({her2_cn}). "
                "Per 2018 ASCO/CAP update, reclassified as POSITIVE (was equivocal in 2007)."
            )
        else:
            result["fish_status"] = "Equivocal"
            result["fish_group"] = "Group 4 (2007: Equivocal)"
            result["description"] = (
                f"HER2/CEP17 ratio >=2.0 ({ratio}) but HER2 copy number <4.0 ({her2_cn}). "
                "Per 2007 ASCO/CAP, this was equivocal."
            )

    elif ratio < 2.0 and her2_cn >= 6.0:
        result["fish_status"] = "Positive"
        result["fish_group"] = "Group 2"
        result["description"] = (
            f"HER2/CEP17 ratio <2.0 ({ratio}) but HER2 copy number >=6.0 ({her2_cn}). "
            "Positive for HER2 amplification (ratio-independent criterion)."
        )

    elif ratio < 2.0 and 4.0 <= her2_cn < 6.0:
        if guideline_year == 2018:
            result["fish_status"] = "Negative"
            result["fish_group"] = "Group 5 (2018: Negative)"
            result["description"] = (
                f"HER2/CEP17 ratio <2.0 ({ratio}) AND HER2 copy number 4.0-5.9 ({her2_cn}). "
                "Per 2018 ASCO/CAP update, reclassified as NEGATIVE (was equivocal in 2007)."
            )
        else:
            result["fish_status"] = "Equivocal"
            result["fish_group"] = "Group 5 (2007: Equivocal)"
            result["description"] = (
                f"HER2/CEP17 ratio <2.0 ({ratio}) AND HER2 copy number 4.0-5.9 ({her2_cn}). "
                "Per 2007 ASCO/CAP, this was equivocal."
            )

    else:  # ratio < 2.0 and her2_cn < 4.0
        result["fish_status"] = "Negative"
        result["fish_group"] = "Group 3"
        result["description"] = (
            f"HER2/CEP17 ratio <2.0 ({ratio}) AND HER2 copy number <4.0 ({her2_cn}). "
            "Negative for HER2 amplification."
        )

    # Treatment implications
    if result["fish_status"] == "Positive":
        result["treatment_implications"] = (
            "Patient is eligible for HER2-targeted therapy: trastuzumab, "
            "pertuzumab, and/or T-DXd (trastuzumab deruxtecan)."
        )
    elif result["fish_status"] == "Negative":
        result["treatment_implications"] = (
            "HER2-targeted therapy is NOT indicated."
        )
    else:
        result["treatment_implications"] = (
            "Equivocal result. Additional testing or repeat FISH recommended."
        )

    return result


# ---------------------------------------------------------------------------
# Combined IHC + FISH Assessment
# ---------------------------------------------------------------------------

def assess_her2_status(ihc_score: int, percent_staining: float = 0.0,
                       her2_copy_number: Optional[float] = None,
                       cep17_copy_number: Optional[float] = None,
                       her2_ceph_ratio: Optional[float] = None,
                       guideline_year: int = 2018) -> Dict[str, Any]:
    """
    Combined HER2 assessment using IHC and optional FISH data.

    Parameters
    ----------
    ihc_score : int
        IHC score (0, 1, 2, or 3).
    percent_staining : float
        Percentage of invasive tumor cells with membrane staining.
    her2_copy_number : float, optional
        Average HER2 FISH copy number per cell.
    cep17_copy_number : float, optional
        Average CEP17 copy number per cell.
    her2_ceph_ratio : float, optional
        Pre-computed HER2/CEP17 ratio.
    guideline_year : int
        2007 or 2018 (default 2018).

    Returns
    -------
    dict with combined assessment including final_status, ihc_result,
    fish_result (if FISH data provided), and treatment recommendations.
    """
    ihc_result = interpret_ihc(ihc_score, percent_staining)

    result = {
        "ihc_result": ihc_result,
        "final_status": ihc_result["her2_status"],
        "method": "IHC only",
    }

    # If IHC is 0 or 3+, FISH is not needed for determination
    if ihc_score in (0, 3):
        result["fish_result"] = None
        result["notes"] = (
            f"IHC {ihc_score}+ is definitive. FISH testing is not required."
        )

    # If IHC is 1+, generally negative but may reflex
    elif ihc_score == 1:
        if her2_copy_number is not None and cep17_copy_number is not None:
            fish_result = interpret_fish(her2_copy_number, cep17_copy_number,
                                         her2_ceph_ratio, guideline_year=guideline_year)
            result["fish_result"] = fish_result
            result["method"] = "IHC + FISH"
            result["final_status"] = fish_result["fish_status"]
            result["notes"] = (
                "IHC 1+ is negative. FISH was performed per clinical concern."
            )
        else:
            result["fish_result"] = None
            result["notes"] = (
                "IHC 1+ is negative. FISH may be considered if clinical suspicion."
            )

    # If IHC is 2+ (equivocal), FISH is required
    elif ihc_score == 2:
        if her2_copy_number is not None and cep17_copy_number is not None:
            fish_result = interpret_fish(her2_copy_number, cep17_copy_number,
                                         her2_ceph_ratio, guideline_year=guideline_year)
            result["fish_result"] = fish_result
            result["method"] = "IHC + FISH"
            result["final_status"] = fish_result["fish_status"]
            result["notes"] = (
                "IHC 2+ is equivocal. FISH result determines final HER2 status."
            )
        else:
            result["fish_result"] = None
            result["final_status"] = "Equivocal"
            result["notes"] = (
                "IHC 2+ is equivocal. FISH testing is REQUIRED to determine final status."
            )

    # Treatment recommendations based on final status
    if result["final_status"] == "Positive":
        result["treatment_recommendations"] = {
            "trastuzumab": "Indicated (Herceptin)",
            "pertuzumab": "Indicated in combination with trastuzumab (Perjeta)",
            "t_dx": "Indicated if prior therapy (Enhertu, trastuzumab deruxtecan)",
            "t_dm1": "Indicated for metastatic disease (Kadcyla, T-DM1)",
        }
    elif result["final_status"] == "Negative":
        result["treatment_recommendations"] = {
            "trastuzumab": "NOT indicated",
            "pertuzumab": "NOT indicated",
            "t_dx": "NOT indicated",
            "t_dm1": "NOT indicated",
        }
    else:
        result["treatment_recommendations"] = {
            "note": "Equivocal status - complete FISH testing before treatment decisions."
        }

    return result


# ---------------------------------------------------------------------------
# Batch Processing
# ---------------------------------------------------------------------------

def process_batch(input_csv: str, output_csv: str) -> int:
    """
    Process a CSV of HER2 cases and write results.

    Expected CSV columns: ihc_score, percent_staining,
    Optional: her2_copy_number, cep17_copy_number, her2_ceph_ratio, guideline_year
    """
    with open(input_csv, mode="r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)

    out_fields = fieldnames + ["final_status", "method", "notes"]
    out_rows = []

    for row in rows:
        ihc = int(row.get("ihc_score", 0))
        pct = float(row.get("percent_staining", 0))

        her2_cn = row.get("her2_copy_number")
        cep17 = row.get("cep17_copy_number")
        ratio = row.get("her2_ceph_ratio")
        year = int(row.get("guideline_year", 2018))

        fish_cn = float(her2_cn) if her2_cn else None
        fish_cep = float(cep17) if cep17 else None
        fish_ratio = float(ratio) if ratio else None

        assessment = assess_her2_status(
            ihc_score=ihc,
            percent_staining=pct,
            her2_copy_number=fish_cn,
            cep17_copy_number=fish_cep,
            her2_ceph_ratio=fish_ratio,
            guideline_year=year,
        )

        row_dict = dict(row)
        row_dict["final_status"] = assessment["final_status"]
        row_dict["method"] = assessment["method"]
        row_dict["notes"] = assessment.get("notes", "")
        out_rows.append(row_dict)

    with open(output_csv, mode="w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=out_fields)
        writer.writeheader()
        writer.writerows(out_rows)

    return len(out_rows)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(
        description="ASCO/CAP HER2 IHC & FISH Interpretation Algorithm"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # IHC subcommand
    p_ihc = subparsers.add_parser("ihc", help="Interpret IHC score")
    p_ihc.add_argument("--score", type=int, required=True, choices=[0, 1, 2, 3],
                        help="IHC score (0, 1, 2, or 3)")
    p_ihc.add_argument("--percent", type=float, default=0.0,
                        help="Percent of cells with staining (0-100)")

    # FISH subcommand
    p_fish = subparsers.add_parser("fish", help="Interpret FISH results")
    p_fish.add_argument("--her2-cn", type=float, required=True,
                        help="Average HER2 copy number per cell")
    p_fish.add_argument("--cep17-cn", type=float, required=True,
                        help="Average CEP17 copy number per cell")
    p_fish.add_argument("--ratio", type=float, default=None,
                        help="HER2/CEP17 ratio (auto-calculated if omitted)")
    p_fish.add_argument("--guideline-year", type=int, default=2018,
                        choices=[2007, 2018], help="ASCO/CAP guideline year")

    # Combined assessment
    p_combined = subparsers.add_parser("assess", help="Combined IHC + FISH assessment")
    p_combined.add_argument("--ihc", type=int, required=True, choices=[0, 1, 2, 3],
                            help="IHC score")
    p_combined.add_argument("--percent", type=float, default=0.0,
                            help="Percent staining")
    p_combined.add_argument("--her2-cn", type=float, default=None,
                            help="HER2 copy number (optional)")
    p_combined.add_argument("--cep17-cn", type=float, default=None,
                            help="CEP17 copy number (optional)")
    p_combined.add_argument("--ratio", type=float, default=None,
                            help="HER2/CEP17 ratio (optional)")
    p_combined.add_argument("--guideline-year", type=int, default=2018,
                            choices=[2007, 2018])

    # Batch
    p_batch = subparsers.add_parser("batch", help="Batch process CSV")
    p_batch.add_argument("-i", "--input", required=True, help="Input CSV path")
    p_batch.add_argument("-o", "--output", default="results.csv", help="Output CSV path")

    args = parser.parse_args(argv)

    if args.command == "ihc":
        result = interpret_ihc(args.score, args.percent)
        print(json.dumps(result, indent=2))

    elif args.command == "fish":
        result = interpret_fish(
            her2_copy_number=args.her2_cn,
            cep17_copy_number=args.cep17_cn,
            her2_ceph_ratio=args.ratio,
            guideline_year=args.guideline_year,
        )
        print(json.dumps(result, indent=2))

    elif args.command == "assess":
        result = assess_her2_status(
            ihc_score=args.ihc,
            percent_staining=args.percent,
            her2_copy_number=args.her2_cn,
            cep17_copy_number=args.cep17_cn,
            her2_ceph_ratio=args.ratio,
            guideline_year=args.guideline_year,
        )
        # Simplify output (remove nested for CLI readability)
        output = {
            "final_status": result["final_status"],
            "method": result["method"],
            "ihc_status": result["ihc_result"]["her2_status"],
            "fish_status": result["fish_result"]["fish_status"] if result["fish_result"] else "N/A",
            "notes": result["notes"],
            "treatment_recommendations": result["treatment_recommendations"],
        }
        print(json.dumps(output, indent=2))

    elif args.command == "batch":
        count = process_batch(args.input, args.output)
        print(f"Processed {count} records -> {args.output}")


if __name__ == "__main__":
    main()
