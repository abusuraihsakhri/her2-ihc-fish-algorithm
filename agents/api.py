"""
FastAPI REST API Server for Her2 Ihc Fish Algorithm.
"""
from typing import Dict, Any, List, Optional
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from .base import AuditLogger, PHIGuard
from .models import SystemTaskPayload, ConsensusDossier
from .supervisor import SystemSupervisor
from her2_algorithm import assess_her2_status

supervisor = SystemSupervisor(model_provider="mock")

app = FastAPI(
    title="Her2 Ihc Fish Algorithm API",
    description="Enterprise Distributed Component Platform (Digital Pathology & Histology Systems)",
    version="3.0.0-ENTERPRISE",
)


class ChatRequest(BaseModel):
    query: str


@app.get("/health")
def health():
    return {"status": "HEALTHY", "service": "her2-ihc-fish-algorithm", "domain": "Digital Pathology & Histology Systems", "standard": "CAP Cancer Protocols / DICOM WSI PS3.16", "version": "3.0.0-ENTERPRISE"}


@app.get("/metrics")
def metrics():
    return {
        "dossiers_processed_total": len(supervisor.dossier_registry),
        "audit_blocks_total": len(AuditLogger.get_trail()),
        "system_status": "NOMINAL_OPTIMAL"
    }


@app.post("/api/audit")
def api_audit(payload: SystemTaskPayload):
    dossier = supervisor.process_task(payload)
    return dossier.to_dict()


@app.post("/api/chat")
def api_chat(req: ChatRequest):
    try:
        ans = supervisor.query_supervisory_chat(req.query)
        return {"response": ans}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.get("/api/audit/logs")
def api_audit_logs():
    return {"audit_trail": AuditLogger.get_trail(), "verified": AuditLogger.verify_integrity()}


class HER2AssessmentRequest(BaseModel):
    """De-identified breast cancer HER2 measurements; no patient identifiers."""
    ihc_score: int
    percent_staining: float = 0.0
    her2_copy_number: Optional[float] = None
    cep17_copy_number: Optional[float] = None
    her2_ceph_ratio: Optional[float] = None
    guideline_year: int = 2023
    recount_her2_copy_number: Optional[float] = None
    recount_cep17_copy_number: Optional[float] = None


@app.post("/api/her2")
def api_her2_assessment(request: HER2AssessmentRequest):
    """Execute the same 2018/2023 algorithm used by the CLI and browser."""
    try:
        return assess_her2_status(**request.model_dump(exclude_none=True))
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
