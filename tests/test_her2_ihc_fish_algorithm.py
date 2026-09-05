"""
Automated Pytest Test Suite for Her2 Ihc Fish Algorithm.
Domain: Digital Pathology & Histology Systems
Standard: CAP Cancer Protocols / DICOM WSI PS3.16
"""
import os
import sys
from pathlib import Path

# Set required environment variable before importing agents
os.environ.setdefault("AUDIT_SECRET_KEY", "test-audit-key-for-unit-tests-only-do-not-use-in-production")

sys.path.insert(0, str(Path(__file__).parent.parent))

import pytest
from agents.base import PHIGuard, AuditLogger, SecurityException
from agents.models import SystemTaskPayload, UrgencyLevel, SystemIntegrityStatus
from agents.workers import InvariantQCWorker, SafetyEscalationWorker, ProtocolConformanceWorker
from agents.supervisor import SystemSupervisor
from cli import main


def test_phi_guard_enforcement():
    with pytest.raises(SecurityException):
        PHIGuard.assert_no_phi("Patient MRN-994827 blood culture positive for Staphylococcus")

    # Clean text passes
    PHIGuard.assert_no_phi("Analytical assay specimen KEY-001 optimal")


def test_phi_guard_redaction():
    """Test that PHI is properly redacted."""
    text = "Patient MRN-123456 has SSN 123-45-6789"
    redacted = PHIGuard.redact_phi(text)
    assert "MRN-123456" not in redacted
    assert "123-45-6789" not in redacted
    assert "[REDACTED_IDENTIFIER]" in redacted


def test_audit_trail_requires_secret_key():
    """Test that AuditTrail requires a secret key."""
    import os
    # Save and clear the env var
    original_key = os.environ.pop("AUDIT_SECRET_KEY", None)
    try:
        from agents.base import AuditTrail
        with pytest.raises(RuntimeError, match="AUDIT_SECRET_KEY"):
            AuditTrail()
    finally:
        # Restore the env var
        if original_key:
            os.environ["AUDIT_SECRET_KEY"] = original_key


def test_specialized_workers():
    # Worker 1: QC Invariant
    p1 = SystemTaskPayload(task_id="T1", target_identifier="KEY-01", primary_metric=35.0)
    alerts1 = InvariantQCWorker.evaluate(p1)
    assert len(alerts1) == 1
    assert alerts1[0].urgency == UrgencyLevel.ELEVATED

    # Worker 2: Safety
    p2 = SystemTaskPayload(task_id="T2", target_identifier="KEY-02", primary_metric=10.0, is_critical_flag=True)
    alerts2 = SafetyEscalationWorker.evaluate(p2)
    assert len(alerts2) == 1
    assert alerts2[0].urgency == UrgencyLevel.CRITICAL_STAT

    # Worker 3: Protocol Conformance
    p3 = SystemTaskPayload(task_id="T3", target_identifier="KEY-03", primary_metric=10.0, status_descriptor="DISCORDANT_ANOMALY")
    alerts3 = ProtocolConformanceWorker.evaluate(p3)
    assert len(alerts3) == 1


def test_supervisor_consensus_and_audit():
    supervisor = SystemSupervisor(model_provider="mock")
    payload = SystemTaskPayload(
        task_id="TASK-PROD-01",
        target_identifier="KEY-PROD-01",
        primary_metric=12.0,
        secondary_metric=4.0,
        status_descriptor="NOMINAL"
    )
    dossier = supervisor.process_task(payload)
    assert dossier.overall_urgency == UrgencyLevel.ROUTINE
    assert dossier.integrity_status == SystemIntegrityStatus.VALIDATED
    assert dossier.audit_hash != ""

    # Verify cryptographic audit trail
    assert AuditLogger.verify_integrity() is True

    # CLI tests - verify existing subcommands work without error
    main(["ihc", "--score", "3"])
    main(["fish", "--her2-cn", "8.0", "--cep17-cn", "2.0"])
    main(["assess", "--ihc", "2", "--percent", "30.0", "--her2-cn", "10.0", "--cep17-cn", "2.0"])
