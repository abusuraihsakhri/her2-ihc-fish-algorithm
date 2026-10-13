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
import math
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
    if type(ihc_score) is not int or ihc_score not in (0, 1, 2, 3):
        raise ValueError(f"IHC score must be 0, 1, 2, or 3; got {ihc_score}")
    if isinstance(percent_staining, bool) or not isinstance(percent_staining, (int, float)) or not math.isfinite(percent_staining) or not 0.0 <= percent_staining <= 100.0:
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
            "This result does not independently determine treatment eligibility."
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
            "HER2 is negative by IHC. Further evaluation may be considered "
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
            "HER2 protein overexpression detected; clinical eligibility for "
            "specific treatments requires clinical assessment."
        )

    return result


# ---------------------------------------------------------------------------
# FISH Interpretation
# ---------------------------------------------------------------------------

def _fish_group(ratio: float, copies: float) -> int:
    """ASCO/CAP 2018 dual-probe ISH group (reaffirmed in 2023)."""
    if ratio >= 2:
        return 1 if copies >= 4 else 2
    if copies >= 6:
        return 3
    return 4 if copies >= 4 else 5


def _finite_number(value: float, name: str, *, strictly_positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{name} must be a finite number")
    if (value <= 0) if strictly_positive else (value < 0):
        raise ValueError(f"{name} must be {'> 0' if strictly_positive else '>= 0'}; got {value}")
    return float(value)


def interpret_fish(her2_copy_number: float, cep17_copy_number: float,
                   her2_ceph_ratio: Optional[float] = None,
                   num_cells_counted: int = 20,
                   guideline_year: int = 2018,
                   ihc_score: Optional[int] = None,
                   recount_her2_copy_number: Optional[float] = None,
                   recount_cep17_copy_number: Optional[float] = None,
                   recount_ratio: Optional[float] = None,
                   recount_cells_counted: int = 20) -> Dict[str, Any]:
    """Classify dual-probe ISH, with required IHC/recount adjudication for groups 2-4.

    The 2023 update reaffirms the 2018 algorithm. Historical 2007 rules
    are deliberately not offered because this function cannot validate
    the full legacy pathway. Do not use for single-probe ISH.
    """
    if guideline_year not in (2018, 2023):
        raise ValueError("Only the 2018/2023 dual-probe ISH guideline is supported")
    her2 = _finite_number(her2_copy_number, "her2_copy_number")
    cep17 = _finite_number(cep17_copy_number, "cep17_copy_number", strictly_positive=True)
    ratio = her2 / cep17 if her2_ceph_ratio is None else _finite_number(her2_ceph_ratio, "her2_ceph_ratio")
    if her2_ceph_ratio is not None and abs(ratio - her2 / cep17) > 0.1:
        raise ValueError("Provided HER2/CEP17 ratio disagrees with the supplied copy numbers")
    if type(num_cells_counted) is not int or num_cells_counted < 20:
        raise ValueError("num_cells_counted must be an integer >= 20")
    if ihc_score is not None and (type(ihc_score) is not int or ihc_score not in (0, 1, 2, 3)):
        raise ValueError("ihc_score must be 0, 1, 2, or 3")

    group = _fish_group(ratio, her2)  # NEVER classify using rounded display values.
    recount_inputs = (recount_her2_copy_number, recount_cep17_copy_number)
    if (recount_inputs[0] is None) != (recount_inputs[1] is None):
        raise ValueError("Both recount HER2 and CEP17 copy numbers are required")
    if recount_ratio is not None and recount_inputs[0] is None:
        raise ValueError("Recount copy numbers are required with recount_ratio")

    recount_group = None
    if recount_inputs[0] is not None:
        r_her2 = _finite_number(recount_inputs[0], "recount_her2_copy_number")
        r_cep17 = _finite_number(recount_inputs[1], "recount_cep17_copy_number", strictly_positive=True)
        r_ratio = r_her2 / r_cep17 if recount_ratio is None else _finite_number(recount_ratio, "recount_ratio")
        if recount_ratio is not None and abs(r_ratio - r_her2 / r_cep17) > 0.1:
            raise ValueError("Recount ratio disagrees with recount copy numbers")
        if type(recount_cells_counted) is not int or recount_cells_counted < 20:
            raise ValueError("recount_cells_counted must be an integer >= 20")
        recount_group = _fish_group(r_ratio, r_her2)

    if group == 1:
        status, note = "Positive", "Dual-probe ISH group 1: amplified."
    elif group == 5:
        status, note = "Negative", "Dual-probe ISH group 5: not amplified."
    elif ihc_score is None:
        status, note = "Pending IHC", "Concurrent IHC is required for ISH groups 2-4."
    elif ihc_score == 3:
        status, note = "Positive", "Concurrent IHC 3+ supports HER2 positivity."
    elif ihc_score in (0, 1):
        status, note = "Negative", "Concurrent IHC 0/1+ supports HER2 negativity."
    elif recount_group is None:
        status, note = "Pending recount", (
            "IHC 2+ requires a blinded independent recount of at least 20 cells "
            "for ISH groups 2-4."
        )
    elif recount_group != group:
        status, note = "Indeterminate", (
            "ISH group changed on blinded recount; adjudicate according to laboratory procedures."
        )
    else:
        status = "Positive" if group == 3 else "Negative"
        note = (
            "Blinded independent recount remains in ISH group "
            f"{group}; combined IHC/ISH interpretation is {status.lower()}."
        )

    return {
        "her2_copy_number": round(her2, 2),
        "cep17_copy_number": round(cep17, 2),
        "ratio": round(ratio, 3),
        "num_cells_counted": num_cells_counted,
        "guideline_year": guideline_year,
        "fish_group": f"Group {group}",
        "fish_status": status,
        "description": note,
        "recount_group": f"Group {recount_group}" if recount_group is not None else None,
        "treatment_implications": "This test does not independently determine treatment eligibility.",
    }


# ---------------------------------------------------------------------------
# Combined IHC + FISH Assessment
# ---------------------------------------------------------------------------

def assess_her2_status(ihc_score: int, percent_staining: float = 0.0,
                       her2_copy_number: Optional[float] = None,
                       cep17_copy_number: Optional[float] = None,
                       her2_ceph_ratio: Optional[float] = None,
                       guideline_year: int = 2018,
                       recount_her2_copy_number: Optional[float] = None,
                       recount_cep17_copy_number: Optional[float] = None) -> Dict[str, Any]:
    """Combine pathologist-assigned IHC and dual-probe ISH results.

    Input IHC score is assumed to have been validated by the pathologist;
    staining percentage alone cannot establish the score.
    """
    ihc = interpret_ihc(ihc_score, percent_staining)
    if (her2_copy_number is None) != (cep17_copy_number is None):
        raise ValueError("Both HER2 and CEP17 copy numbers must be supplied together")
    if guideline_year not in (2018, 2023):
        raise ValueError("Only the 2018/2023 dual-probe ISH guideline is supported")
    result = {
        "ihc_result": ihc,
        "final_status": ihc["her2_status"],
        "method": "IHC only",
        "fish_result": None,
    }
    if her2_copy_number is not None:
        fish = interpret_fish(
            her2_copy_number, cep17_copy_number, her2_ceph_ratio,
            guideline_year=guideline_year, ihc_score=ihc_score,
            recount_her2_copy_number=recount_her2_copy_number,
            recount_cep17_copy_number=recount_cep17_copy_number,
        )
        result["fish_result"] = fish
        result["method"] = "IHC + dual-probe ISH"
        if ihc_score == 3:
            result["final_status"] = "Positive"
            result["notes"] = ("IHC 3+ is positive. Review discordant ISH results in context."
                               if fish["fish_group"] == "Group 5" else
                               "IHC 3+ confirms HER2 protein overexpression.")
        else:
            result["final_status"] = fish["fish_status"]
            result["notes"] = fish["description"]
    elif ihc_score == 2:
        result["notes"] = "IHC 2+ is equivocal: reflex ISH is required."
    else:
        result["notes"] = "IHC result reported without dual-probe ISH."
    if her2_copy_number is None and (recount_her2_copy_number is not None or recount_cep17_copy_number is not None):
        raise ValueError("A recount cannot be supplied without initial ISH measurements")
    result["treatment_recommendations"] = {
        "note": (
            "HER2 reporting does not establish treatment eligibility. "
            "Therapy depends on disease setting, complete pathology, prior treatment, "
            "drug-specific indications and multidisciplinary evaluation. "
            "IHC 1+ or IHC 2+/ISH-negative disease may qualify for certain "
            "antibody-drug conjugates in appropriate clinical settings."
        )
    }
    return result


# ---------------------------------------------------------------------------
# Batch Processing
# ---------------------------------------------------------------------------

def process_batch(input_csv: str, output_csv: str) -> int:
    """Process CSV rows without silently inventing missing IHC results.

    Each invalid row raises a descriptive error before the output file is opened.
    """
    with open(input_csv, mode="r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        fields = list(reader.fieldnames or [])
        if "ihc_score" not in fields:
            raise ValueError("Missing required CSV column: ihc_score")
        rows = list(reader)
    out_rows = []
    for line, row in enumerate(rows, start=2):
        try:
            if not row.get("ihc_score", "").strip():
                raise ValueError("ihc_score is required")
            kwargs = {
                "ihc_score": int(row["ihc_score"]),
                "percent_staining": float(row.get("percent_staining") or 0),
                "guideline_year": int(row.get("guideline_year") or 2018),
            }
            for key in ("her2_copy_number", "cep17_copy_number", "her2_ceph_ratio",
                        "recount_her2_copy_number", "recount_cep17_copy_number"):
                if row.get(key) not in (None, ""):
                    kwargs[key] = float(row[key])
            assessment = assess_her2_status(**kwargs)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"CSV line {line}: {exc}") from exc
        safe_row = {}
        for key in fields:
            val = str(row.get(key) or "")
            # Prevent spreadsheet formula execution when CSV is opened.
            safe_row[key] = "'" + val if val.lstrip().startswith(("=", "+", "-", "@")) else val
        safe_row.update(
            final_status=assessment["final_status"],
            method=assessment["method"],
            notes=assessment["notes"],
        )
        out_rows.append(safe_row)
    out_fields = list(dict.fromkeys(fields + ["final_status", "method", "notes"]))
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
                        choices=[2018, 2023], help="ASCO/CAP guideline year")

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

    # Serve (FastAPI)
    p_serve = subparsers.add_parser("serve", help="Start FastAPI REST server")
    p_serve.add_argument("--host", default="0.0.0.0", help="Host to bind (default: 0.0.0.0)")
    p_serve.add_argument("--port", type=int, default=8000, help="Port to bind (default: 8000)")

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

    elif args.command == "serve":
        import uvicorn
        from agents.api import app
        uvicorn.run(app, host=args.host, port=args.port)


if __name__ == "__main__":
    main()
