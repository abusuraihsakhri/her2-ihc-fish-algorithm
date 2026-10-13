"use strict";

// Python is the sole source of interpretation logic. No inputs are sent to an API.
const form = document.getElementById("case-form");
const runButton = document.getElementById("calculate");
const runtimeStatus = document.getElementById("runtime-status");
const errorEl = document.getElementById("error");
const resultBox = document.getElementById("result");
const emptyEl = document.getElementById("empty");
const exportButton = document.getElementById("export");
let runtime = null;
let lastResult = null;

function renderText(id, value) {
  document.getElementById(id).textContent = value == null ? "—" : String(value);
}
function optionalNumber(id) {
  const value = document.getElementById(id).value.trim();
  if (!value) return null;
  const num = Number(value);
  if (!Number.isFinite(num)) throw Error("Measurements must be finite numbers.");
  return num;
}
function collectInput() {
  const values = {
    ihc_score: Number(document.getElementById("ihc").value),
    percent_staining: optionalNumber("percent"),
    guideline_year: 2023,
  };
  for (const [inputId, pythonKey] of [
    ["her2", "her2_copy_number"],
    ["cep17", "cep17_copy_number"],
    ["recount-her2", "recount_her2_copy_number"],
    ["recount-cep17", "recount_cep17_copy_number"],
  ]) {
    const value = optionalNumber(inputId);
    if (value !== null) values[pythonKey] = value;
  }
  if ((values.her2_copy_number === undefined) !== (values.cep17_copy_number === undefined)) {
    throw Error("Enter both initial HER2 and CEP17 means, or leave both blank.");
  }
  if ((values.recount_her2_copy_number === undefined) !==
      (values.recount_cep17_copy_number === undefined)) {
    throw Error("Enter both independent recount means, or leave both blank.");
  }
  if (!Number.isFinite(values.percent_staining)) throw Error("Enter a valid staining percentage.");
  return values;
}
function showResult(data) {
  lastResult = data;
  errorEl.textContent = "";
  const final = data.final_status;
  const statusClass = final === "Positive" ? "positive" : final === "Negative" ? "negative" : "pending";
  const state = document.getElementById("result-status");
  state.className = "result-state " + statusClass;
  state.textContent = final;
  renderText("result-note", data.notes);
  renderText("result-ihc", data.ihc_result.ihc_score + (data.ihc_result.ihc_score ? "+" : "") + " (" + data.ihc_result.her2_status + ")");
  renderText("result-group", data.fish_result && data.fish_result.fish_group);
  renderText("result-fish", data.fish_result && data.fish_result.fish_status);
  renderText("result-ratio", data.fish_result && data.fish_result.ratio);
  renderText("result-method", data.method);
  emptyEl.hidden = true;
  resultBox.hidden = false;
}
form.addEventListener("submit", (event) => {
  event.preventDefault();
  errorEl.textContent = "";
  if (!runtime) return;
  try {
    const data = collectInput();
    runtime.globals.set("input_json", JSON.stringify(data));
    const json = runtime.runPython(
      "json.dumps(her2_algorithm.assess_her2_status(**json.loads(input_json)), allow_nan=False)"
    );
    showResult(JSON.parse(json));
  } catch (error) {
    lastResult = null;
    resultBox.hidden = true;
    emptyEl.hidden = false;
    errorEl.textContent = String(error?.message || error).replace(/^.*ValueError:\s*/, "");
  }
});
form.addEventListener("reset", () => {
  lastResult = null;
  resultBox.hidden = true;
  emptyEl.hidden = false;
  errorEl.textContent = "";
});
exportButton.addEventListener("click", () => {
  if (!lastResult) return;
  const blob = new Blob([JSON.stringify(lastResult, null, 2) + "\n"], {type: "application/json"});
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = "her2-interpretation.json";
  link.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
});
async function initializePython() {
  try {
    if (typeof window.loadPyodide !== "function") throw Error("Python runtime script unavailable.");
    runtime = await window.loadPyodide({ indexURL: "https://cdn.jsdelivr.net/pyodide/v0.26.4/full/" });
    const response = await fetch("./her2_algorithm.py");
    if (!response.ok) throw Error("Python algorithm could not be loaded from this site.");
    const source = await response.text();
    runtime.FS.writeFile("/home/pyodide/her2_algorithm.py", source);
    runtime.runPython("import json, sys; sys.path.insert(0, '/home/pyodide'); import her2_algorithm");
    runButton.disabled = false;
    runtimeStatus.textContent = "Python runtime ready · local execution";
  } catch (error) {
    runtime = null;
    runButton.disabled = true;
    runtimeStatus.textContent = "Interpreter unavailable";
    errorEl.textContent = String(error?.message || error) +
      " Check connectivity or open the Python CLI in the repository.";
  }
}
initializePython();
