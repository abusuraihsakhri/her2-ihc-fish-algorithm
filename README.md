# Her2 IHC FISH Algorithm

> **Domain:** Digital Pathology & Quantitative Histopathology  
> **Reference Guidelines & Standards:** `College of American Pathologists (CAP) Synoptic Protocols & DICOM WSI`

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

</div>

---

## 📖 What It Does

ASCO/CAP HER2 IHC & FISH Interpretation Algorithm

Implements the 2018 ASCO/CAP guidelines for HER2 status determination in breast cancer
using IHC scoring (0-3+) and FISH (fluorescence in situ hybridization) interpretation.

References:
  - Wolff AC et al. J Clin Oncol. 2018;36(20):2135-2151 (ASCO/CAP 2018 update)
  - Wolff AC et al. J Clin Oncol. 2007;25(1):118-145 (ASCO/CAP 2007 original)

Zero-dependency Python implementation (stdlib only).
License: MIT

---

## ⚙️ Key Capabilities & Algorithmic Modules

### 🔬 Analytical Functions

- **`interpret_ihc()`**: Interpret HER2 IHC staining result per ASCO/CAP 2018.

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
- **`interpret_fish()`**: Interpret HER2 FISH results per ASCO/CAP 2018 guidelines.

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
- **`assess_her2_status()`**: Combined HER2 assessment using IHC and optional FISH data.

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
- **`process_batch()`**: Process a CSV of HER2 cases and write results.

Expected CSV columns: ihc_score, percent_staining,
Optional: her2_copy_number, cep17_copy_number, her2_ceph_ratio, guideline_year
- **`main()`** — calculates and validates main parameters.

---

## 📐 Mathematical Formulation & Logic

```text
  if ihc_score == 0:
  elif ihc_score == 1:
  elif ihc_score == 2:
  elif ihc_score == 3:
  HER2/CEP17 ratio. If None, calculated from copy numbers.
```

---

## 💻 CLI Quickstart & Usage

### 1. Guided Interactive Mode
```bash
python cli.py
```

### 2. Direct Parameterized Evaluation
```bash
python cli.py --input data.csv
```

### Parameter Reference
- `--interactive`: Launch guided terminal interactive wizard.
- `--input <path>`: Evaluate input from JSON or CSV specification.
- `--json`: Output deterministic structured results in JSON format.

### Input Data Schema

| Field | Description | Requirement |
|:------|:------------|:------------|
| `Patient_ID` | Parameter / observation metric | Required |
| `v1` | Parameter / observation metric | Required |
| `v2` | Parameter / observation metric | Required |
| `v3` | Parameter / observation metric | Required |

---

## 🛡️ Security & Enterprise Architecture

* **Zero-PHI Outbound Interceptor:** Active AST and regex inspection blocking SSNs, MRNs, phone numbers, and patient identifiers.
* **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation and state transition.
* **Air-Gapped LLM Reasoning Adapter:** Agnostic integration for local Ollama instances (`llama3`, `mistral`), Claude 3.5 Sonnet, GPT-4o, and deterministic test mocks.
* **Active Learning Bayesian Calibration:** Dynamic tracker updating worker reliability weights and monitoring Brier calibration drift.
* **FastAPI & Prometheus Telemetry:** Exposes OpenAPI 3.1 REST endpoints and operational Prometheus metrics (`/metrics`).

---

## 🧪 Testing & Verification

Run the automated test suite:

```bash
pytest -v
```

Execute high-throughput batch simulation benchmarks:

```bash
python simulator.py --tasks 1000 --concurrency 8
```

---

## 🐳 Container Deployment

```bash
docker build -t her2-ihc-fish-algorithm .
docker run -p 8000:8000 her2-ihc-fish-algorithm
```
