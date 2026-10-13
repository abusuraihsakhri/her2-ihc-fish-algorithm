# HER2 IHC / ISH Interpreter

A Python command-line, CSV and REST interface for interpreting breast cancer HER2 immunohistochemistry and **dual-probe in situ hybridization (ISH/FISH)** measurements under the **2018 ASCO/CAP guideline reaffirmed in 2023**. A static browser interface at \`index.html\` runs the same Python implementation locally with Pyodide.

> **Educational / research reference only.** Not validated as a clinical diagnostic device or a substitute for pathology review, assay quality control, multidisciplinary interpretation or current institutional procedures. **Do not enter patient identifiers.**

## Features

- Interpret an externally assigned IHC 0, 1+, 2+ or 3+ score. The percentage is descriptive; IHC scoring must be performed by a qualified observer.
- Classify dual-probe ISH into **groups 1–5** using unrounded HER2/CEP17 ratio and mean HER2 signals per cell.
- Require concurrent IHC for less common ISH **groups 2–4**; IHC 2+ cases in these groups remain pending until an independent blinded recount of at least 20 cells is supplied. Changed groups require adjudication.
- Support CSV batch interpretation, JSON CLI results, local REST endpoint, and downloadable browser JSON results.
- Reject nonfinite values, inconsistent ratios, incomplete pairs and invalid CSV rows; neutralize spreadsheet-formula prefixes in CSV output.
- Retain an optional, separate sample audit/supervisor API. This audit API is **not** the clinical interpretation algorithm.

## Browser application

The static entry point is \`index.html\`, with logic in \`app.js\` and \`her2_algorithm.py\`. The page loads Pyodide v0.26.4 from jsDelivr, then executes calculations entirely in the browser. JavaScript/WebAssembly and network access to the CDN are required on first load. Large runtime downloads can take longer on slow connections.

**Privacy:** Inputs are not sent to an application API, stored in cookies/localStorage, or persisted by this application. Third-party CDN requests fetch the runtime; no entered measurement values are included in those requests. Do not interpret the browser's local-only workflow as a guarantee of HIPAA compliance.

If deployed using GitHub Pages from the repository root, \`index.html\` is the app entry point. GitHub Pages does not execute the FastAPI service.

## Python CLI

Only Python 3.10+ standard library is needed for the clinical interpretation CLI.

\`\`\`bash
python cli.py ihc --score 2 --percent 30
python cli.py fish --her2-cn 8 --cep17-cn 2 --guideline-year 2023
python cli.py fish --her2-cn 3 --cep17-cn 1 --ihc 2
python cli.py fish --her2-cn 3 --cep17-cn 1 --ihc 2 --recount-her2-cn 3.1 --recount-cep17-cn 1
python cli.py assess --ihc 2 --percent 30 --her2-cn 5 --cep17-cn 3
python cli.py batch -i sample.csv -o results.csv
\`\`\`

For **group 2, 3 or 4** with IHC 2+, a missing independent recount returns \`Pending recount\`, not a definitive HER2-positive or HER2-negative result. A changed ISH group on recount returns \`Indeterminate\`, requiring procedural adjudication. An isolated dual-probe ISH group 2–4 result without IHC returns \`Pending IHC\`.

The historical **2007** interpretation pathway is not supported; requests with \`guideline_year=2007\` are rejected rather than applying an incomplete legacy algorithm.

## CSV schema

\`ihc_score\` is mandatory; optional columns include \`percent_staining\`, \`her2_copy_number\`, \`cep17_copy_number\`, \`her2_ceph_ratio\`, \`guideline_year\`, \`recount_her2_copy_number\`, \`recount_cep17_copy_number\`. Initial and recount signal means must each be supplied as HER2/CEP17 pairs. A row error stops the batch without writing a partial output file.

## REST API (self-hosted)

\`\`\`bash
python -m pip install -r requirements.txt
export AUDIT_SECRET_KEY="$(python -c 'import secrets; print(secrets.token_hex(32))')"
python cli.py serve --host 127.0.0.1 --port 8000
\`\`\`

- \`POST /api/her2\`: interpret structured measurements using the same Python clinical algorithm.
- \`GET /health\`: application health.
- \`POST /api/audit\`, \`GET /api/audit/logs\`, \`POST /api/chat\`, \`GET /metrics\`: retained **nonclinical demonstration** audit/supervisor endpoints. The chat endpoint uses a deterministic mock; it does **not** perform validated clinical inference or call external language models.

The demonstration endpoints have no authentication and the PHI regular-expression filter is **not a comprehensive de-identification mechanism**. Do **not** expose the FastAPI server on a public network or send patient information without implementing appropriate authentication, access controls, threat modeling and compliance review. \`AUDIT_SECRET_KEY\` is required for the audit subsystem; never commit a production key. Signed audit entries are process-local and not a persistent ledger.

## Development and testing

\`\`\`bash
python -m pip install -r requirements-dev.txt
pytest -v
python -m compileall -q her2_algorithm.py agents tests
node --check app.js
\`\`\`

CI covers Python 3.10–3.12, regression tests and audit-integrity checks. \`Dockerfile\` and \`docker-compose.yml\` provide a local server option; supply \`AUDIT_SECRET_KEY\` securely at runtime.

## Clinical scope and sources

The [CAP 2023 HER2 guideline update](https://www.cap.org/cap-guidelines/her2-testing-in-breast-cancer-2023-guideline-update/) reaffirms the 2018 dual-probe ISH interpretation framework. The 2023 guideline does **not** introduce \`HER2-low\` as a distinct HER2 interpretive category, although IHC 1+ and IHC 2+/ISH-negative results may be relevant to therapy selection in specific settings. **This tool does not recommend a drug or evaluate clinical eligibility.**

Key sources:

- Wolff AC et al. *J Clin Oncol*. 2018;36:2105–2122. doi:[10.1200/JCO.2018.77.8738](https://doi.org/10.1200/JCO.2018.77.8738).
- Wolff AC et al. *Arch Pathol Lab Med*. 2023. doi:[10.5858/arpa.2023-0950-SA](https://doi.org/10.5858/arpa.2023-0950-SA).

## Browser compatibility

Modern browsers supporting WebAssembly and ES2020 are expected to work; real-device and cross-browser regression testing have not been completed. Large Pyodide downloads may be blocked by restrictive networks. Python execution works without a browser.

## License

[MIT](LICENSE).
