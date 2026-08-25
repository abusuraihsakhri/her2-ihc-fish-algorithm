# HER2 IHC/FISH Interpretation Algorithm

> **ASCO/CAP HER2 Status Determination for Breast Cancer**
> Reference: Wolff AC et al. J Clin Oncol. 2018;36(20):2135-2151

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10+-3776AB.svg)

## Overview

Real implementation of the ASCO/CAP 2018 HER2 interpretation algorithm for breast cancer, including:

- **IHC scoring** (0-3+) with proper classification and reflex FISH requirements
- **FISH interpretation** with all 5 ASCO/CAP groups (2007 and 2018 guidelines)
- **Combined IHC + FISH assessment** with treatment eligibility determination
- **Treatment implications** for trastuzumab, pertuzumab, T-DXd, T-DM1

## Quick Start

```bash
# IHC interpretation
python her2_algorithm.py ihc --score 3 --percent 90

# FISH interpretation
python her2_algorithm.py fish --her2-cn 8.0 --cep17-cn 2.0

# Combined assessment
python her2_algorithm.py assess --ihc 2 --percent 30 --her2-cn 10.0 --cep17-cn 2.0

# Batch processing
python her2_algorithm.py batch -i cases.csv -o results.csv
```

## IHC Scoring (ASCO/CAP 2018)

| Score | Staining Pattern | Status | Action |
|-------|-----------------|--------|--------|
| 0 | No staining or incomplete, faint in ≤10% | **Negative** | No HER2 therapy |
| 1+ | Faint/barely perceptible incomplete in >10% | **Negative** | Reflex FISH if clinical concern |
| 2+ | Weak to moderate complete in >10% | **Equivocal** | Reflex FISH **required** |
| 3+ | Strong complete membrane in >10% | **Positive** | HER2 therapy eligible |

## FISH Interpretation (ASCO/CAP 2018)

| Group | Ratio | HER2 CN | 2018 Status |
|-------|-------|---------|-------------|
| 1 | ≥2.0 | ≥4.0 | **Positive** |
| 2 | <2.0 | ≥6.0 | **Positive** |
| 3 | <2.0 | <4.0 | **Negative** |
| 4 | ≥2.0 | <4.0 | **Positive** (was equivocal in 2007) |
| 5 | <2.0 | 4.0-5.9 | **Negative** (was equivocal in 2007) |

## Treatment Implications

- **HER2 Positive**: Eligible for trastuzumab (Herceptin), pertuzumab (Perjeta), T-DXd (Enhertu), T-DM1 (Kadcyla)
- **HER2 Negative**: HER2-targeted therapy NOT indicated

## Python API

```python
from her2_algorithm import interpret_ihc, interpret_fish, assess_her2_status

# IHC only
result = interpret_ihc(3, percent_staining=90.0)
print(result["her2_status"])  # "Positive"

# FISH only
result = interpret_fish(her2_copy_number=8.0, cep17_copy_number=2.0)
print(result["fish_status"])  # "Positive"

# Combined
result = assess_her2_status(ihc_score=2, percent_staining=30.0,
                             her2_copy_number=10.0, cep17_copy_number=2.0)
print(result["final_status"])  # "Positive"
```

## Running Tests

```bash
python -m pytest test_her2_algorithm.py -v
```

## License

MIT License. See [LICENSE](LICENSE) for details.
