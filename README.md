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

### 1. IHC Interpretation
```bash
python cli.py ihc --score 2 --percent 30.0
```

### 2. FISH Interpretation
```bash
python cli.py fish --her2-cn 8.0 --cep17-cn 2.0
python cli.py fish --her2-cn 3.0 --cep17-cn 1.0 --guideline-year 2018
```

### 3. Combined Assessment
```bash
python cli.py assess --ihc 2 --percent 30.0 --her2-cn 10.0 --cep17-cn 2.0
```

### 4. Batch CSV Processing
```bash
python cli.py batch -i sample.csv -o results.csv
```

### 5. Start REST API Server
```bash
python cli.py serve --host 0.0.0.0 --port 8000
```

### Input Data Schema (for batch CSV)

| Field | Description | Requirement |
|:------|:------------|:------------|
| `ihc_score` | IHC score (0, 1, 2, or 3) | Required |
| `percent_staining` | Percentage of cells with staining (0-100) | Required |
| `her2_copy_number` | Average HER2 FISH copy number | Optional |
| `cep17_copy_number` | Average CEP17 copy number | Optional |
| `her2_ceph_ratio` | Pre-computed HER2/CEP17 ratio | Optional |
| `guideline_year` | 2007 or 2018 (default: 2018) | Optional |

---

## 🛡️ Security & Enterprise Architecture

* **Zero-PHI Outbound Interceptor:** Active regex inspection blocking SSNs, MRNs, phone numbers, emails, and patient identifiers from outbound audit logs.
* **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation. **Requires `AUDIT_SECRET_KEY` environment variable** - generate with: `python -c "import secrets; print(secrets.token_hex(32))"`
* **Air-Gapped LLM Reasoning Adapter:** Deterministic mock integration with Zero-PHI protection on all prompts.
* **FastAPI REST API:** Exposes health, audit, and chat endpoints with PHI protection.

### Environment Variables

| Variable | Required | Description |
|:---------|:---------|:------------|
| `AUDIT_SECRET_KEY` | Yes (for agents/API) | HMAC-SHA256 secret key for audit trail. Generate a secure random key. |
| `MODEL_PROVIDER` | No | LLM provider: `mock` (default), `ollama`, `claude`, `openai` |

---

## 🧪 Testing & Verification

### Run the automated test suite:

```bash
pytest -v
```

### Run high-throughput simulation benchmarks:

```bash
python simulator.py 1000
```

### Test CLI commands directly:

```bash
python cli.py ihc --score 3 --percent 90.0
python cli.py fish --her2-cn 8.0 --cep17-cn 2.0
python cli.py assess --ihc 2 --percent 30.0 --her2-cn 10.0 --cep17-cn 2.0
python cli.py batch -i sample.csv -o results.csv
```

---

## 🐳 Container Deployment

### Docker

```bash
# Generate a secure audit key
export AUDIT_SECRET_KEY=$(python -c "import secrets; print(secrets.token_hex(32))")

# Build and run
docker build -t her2-ihc-fish-algorithm .
docker run -e AUDIT_SECRET_KEY=$AUDIT_SECRET_KEY -p 8000:8000 her2-ihc-fish-algorithm
```

### Docker Compose

```bash
# Set your audit key in .env or environment
export AUDIT_SECRET_KEY=$(python -c "import secrets; print(secrets.token_hex(32))")
docker-compose up -d
```

> **SECURITY NOTE:** Always provide a secure `AUDIT_SECRET_KEY` at runtime. Never hardcode secrets in Docker images or compose files.
