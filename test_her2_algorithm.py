#!/usr/bin/env python3
"""Tests for HER2 IHC/FISH Algorithm - 20 real test cases."""
import pytest
from her2_algorithm import (
    interpret_ihc, interpret_fish, assess_her2_status, process_batch
)


# ---------------------------------------------------------------------------
# IHC Interpretation Tests
# ---------------------------------------------------------------------------

class TestIHCInterpretation:

    def test_ihc_0_negative(self):
        result = interpret_ihc(0, percent_staining=5.0)
        assert result["her2_status"] == "Negative"
        assert result["category"] == "Negative"
        assert result["reflex_fish_required"] is False

    def test_ihc_0_zero_staining(self):
        result = interpret_ihc(0, percent_staining=0.0)
        assert result["her2_status"] == "Negative"
        assert result["reflex_fish_required"] is False

    def test_ihc_1_negative(self):
        result = interpret_ihc(1, percent_staining=30.0)
        assert result["her2_status"] == "Negative"
        assert result["reflex_fish_required"] is False

    def test_ihc_2_equivocal(self):
        result = interpret_ihc(2, percent_staining=50.0)
        assert result["her2_status"] == "Equivocal"
        assert result["reflex_fish_required"] is True

    def test_ihc_3_positive(self):
        result = interpret_ihc(3, percent_staining=90.0)
        assert result["her2_status"] == "Positive"
        assert result["reflex_fish_required"] is False
        assert "trastuzumab" in result["treatment_implications"].lower()

    def test_ihc_invalid_score(self):
        with pytest.raises(ValueError, match="IHC score must be 0, 1, 2, or 3"):
            interpret_ihc(4, percent_staining=50.0)

    def test_ihc_invalid_percent(self):
        with pytest.raises(ValueError, match="percent_staining must be 0-100"):
            interpret_ihc(2, percent_staining=150.0)

    def test_ihc_3_treatment_implications(self):
        result = interpret_ihc(3, percent_staining=80.0)
        assert "T-DXd" in result["treatment_implications"]


# ---------------------------------------------------------------------------
# FISH Interpretation Tests
# ---------------------------------------------------------------------------

class TestFISHInterpretation:

    def test_fish_positive_group1(self):
        """Ratio >=2.0 AND HER2 CN >=4.0 -> Positive (Group 1)."""
        result = interpret_fish(her2_copy_number=8.0, cep17_copy_number=2.0)
        assert result["fish_status"] == "Positive"
        assert result["fish_group"] == "Group 1"
        assert result["ratio"] == 4.0

    def test_fish_positive_group4_2018(self):
        """Ratio >=2.0 AND HER2 CN <4.0 -> Positive per 2018 (Group 4)."""
        result = interpret_fish(her2_copy_number=3.0, cep17_copy_number=1.0,
                                guideline_year=2018)
        assert result["fish_status"] == "Positive"
        assert "Group 4" in result["fish_group"]

    def test_fish_equivocal_group4_2007(self):
        """Same as above but 2007 guidelines -> Equivocal."""
        result = interpret_fish(her2_copy_number=3.0, cep17_copy_number=1.0,
                                guideline_year=2007)
        assert result["fish_status"] == "Equivocal"
        assert "2007" in result["fish_group"]

    def test_fish_positive_group2(self):
        """Ratio <2.0 AND HER2 CN >=6.0 -> Positive (Group 2)."""
        result = interpret_fish(her2_copy_number=7.0, cep17_copy_number=4.0)
        assert result["fish_status"] == "Positive"
        assert result["fish_group"] == "Group 2"
        assert result["ratio"] < 2.0

    def test_fish_negative_group5_2018(self):
        """Ratio <2.0 AND HER2 CN 4.0-5.9 -> Negative per 2018 (Group 5)."""
        result = interpret_fish(her2_copy_number=5.0, cep17_copy_number=3.0,
                                guideline_year=2018)
        assert result["fish_status"] == "Negative"
        assert "Group 5" in result["fish_group"]

    def test_fish_equivocal_group5_2007(self):
        """Same as above but 2007 -> Equivocal."""
        result = interpret_fish(her2_copy_number=5.0, cep17_copy_number=3.0,
                                guideline_year=2007)
        assert result["fish_status"] == "Equivocal"

    def test_fish_negative_group3(self):
        """Ratio <2.0 AND HER2 CN <4.0 -> Negative (Group 3)."""
        result = interpret_fish(her2_copy_number=2.0, cep17_copy_number=2.0)
        assert result["fish_status"] == "Negative"
        assert result["fish_group"] == "Group 3"

    def test_fish_ratio_auto_calculated(self):
        """Ratio should be auto-calculated when not provided."""
        result = interpret_fish(her2_copy_number=10.0, cep17_copy_number=2.5)
        assert result["ratio"] == 4.0

    def test_fish_invalid_negative_cn(self):
        with pytest.raises(ValueError, match="her2_copy_number must be >= 0"):
            interpret_fish(her2_copy_number=-1.0, cep17_copy_number=2.0)

    def test_fish_invalid_zero_cep17(self):
        with pytest.raises(ValueError, match="cep17_copy_number must be > 0"):
            interpret_fish(her2_copy_number=5.0, cep17_copy_number=0.0)


# ---------------------------------------------------------------------------
# Combined Assessment Tests
# ---------------------------------------------------------------------------

class TestCombinedAssessment:

    def test_ihc_0_no_fish_needed(self):
        result = assess_her2_status(ihc_score=0, percent_staining=5.0)
        assert result["final_status"] == "Negative"
        assert result["fish_result"] is None
        assert "definitive" in result["notes"].lower()

    def test_ihc_3_no_fish_needed(self):
        result = assess_her2_status(ihc_score=3, percent_staining=90.0)
        assert result["final_status"] == "Positive"
        assert result["fish_result"] is None

    def test_ihc_2_with_fish_positive(self):
        result = assess_her2_status(
            ihc_score=2, percent_staining=30.0,
            her2_copy_number=10.0, cep17_copy_number=2.0
        )
        assert result["final_status"] == "Positive"
        assert result["fish_result"]["fish_status"] == "Positive"
        assert result["method"] == "IHC + FISH"

    def test_ihc_2_with_fish_negative(self):
        result = assess_her2_status(
            ihc_score=2, percent_staining=20.0,
            her2_copy_number=2.0, cep17_copy_number=2.0
        )
        assert result["final_status"] == "Negative"
        assert result["fish_result"]["fish_status"] == "Negative"

    def test_ihc_2_no_fish_stays_equivocal(self):
        result = assess_her2_status(ihc_score=2, percent_staining=25.0)
        assert result["final_status"] == "Equivocal"
        assert result["fish_result"] is None
        assert "REQUIRED" in result["notes"]

    def test_ihc_1_with_fish(self):
        result = assess_her2_status(
            ihc_score=1, percent_staining=15.0,
            her2_copy_number=8.0, cep17_copy_number=2.0
        )
        assert result["final_status"] == "Positive"
        assert result["method"] == "IHC + FISH"

    def test_treatment_positive(self):
        result = assess_her2_status(ihc_score=3, percent_staining=80.0)
        recs = result["treatment_recommendations"]
        assert recs["trastuzumab"] == "Indicated (Herceptin)"

    def test_treatment_negative(self):
        result = assess_her2_status(ihc_score=0, percent_staining=0.0)
        recs = result["treatment_recommendations"]
        assert recs["trastuzumab"] == "NOT indicated"


# ---------------------------------------------------------------------------
# Batch Processing Tests
# ---------------------------------------------------------------------------

class TestBatchProcessing:

    def test_batch_basic(self, tmp_path):
        csv_in = tmp_path / "in.csv"
        csv_out = tmp_path / "out.csv"
        csv_in.write_text(
            "ihc_score,percent_staining,her2_copy_number,cep17_copy_number\n"
            "3,90.0,,\n"
            "0,0.0,,\n"
            "2,30.0,10.0,2.0\n",
            encoding="utf-8"
        )
        count = process_batch(str(csv_in), str(csv_out))
        assert count == 3
        content = csv_out.read_text(encoding="utf-8")
        assert "Positive" in content
        assert "Negative" in content

    def test_batch_with_2018_guideline(self, tmp_path):
        csv_in = tmp_path / "in.csv"
        csv_out = tmp_path / "out.csv"
        csv_in.write_text(
            "ihc_score,percent_staining,her2_copy_number,cep17_copy_number,guideline_year\n"
            "2,25.0,3.0,1.0,2018\n",
            encoding="utf-8"
        )
        count = process_batch(str(csv_in), str(csv_out))
        assert count == 1
        content = csv_out.read_text(encoding="utf-8")
        assert "Positive" in content  # Group 4 is positive in 2018
