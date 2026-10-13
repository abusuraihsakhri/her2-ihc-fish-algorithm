"""ASCO/CAP 2018 dual-probe ISH regression cases (reaffirmed in 2023)."""
import csv
import math
import pytest
from her2_algorithm import interpret_fish, assess_her2_status, interpret_ihc, process_batch


@pytest.mark.parametrize("her2,cep17,group", [
    (8.0, 2.0, "Group 1"),
    (3.0, 1.0, "Group 2"),
    (7.0, 4.0, "Group 3"),
    (5.0, 3.0, "Group 4"),
    (2.0, 2.0, "Group 5"),
])
def test_correct_ish_group_numbering(her2, cep17, group):
    assert interpret_fish(her2, cep17)["fish_group"] == group


@pytest.mark.parametrize("group_values,ihc,expected", [
    ((3.0, 1.0), 0, "Negative"),
    ((3.0, 1.0), 1, "Negative"),
    ((3.0, 1.0), 3, "Positive"),
    ((7.0, 4.0), 0, "Negative"),
    ((7.0, 4.0), 1, "Negative"),
    ((7.0, 4.0), 3, "Positive"),
    ((5.0, 3.0), 0, "Negative"),
    ((5.0, 3.0), 1, "Negative"),
    ((5.0, 3.0), 3, "Positive"),
])
def test_ihc_adjudication(group_values, ihc, expected):
    assert interpret_fish(*group_values, ihc_score=ihc)["fish_status"] == expected


@pytest.mark.parametrize("her2,cep17,recount_h,recount_c,expected", [
    (3.0, 1.0, 3.1, 1.0, "Negative"),   # Group 2 remains group 2.
    (7.0, 4.0, 7.2, 4.0, "Positive"),   # Group 3 remains group 3.
    (5.0, 3.0, 5.1, 3.0, "Negative"),   # Group 4 remains group 4.
    (5.0, 3.0, 6.5, 4.0, "Indeterminate"),  # Recount changes group.
])
def test_ihc2_recount(her2, cep17, recount_h, recount_c, expected):
    pending = interpret_fish(her2, cep17, ihc_score=2)
    assert pending["fish_status"] == "Pending recount"
    result = interpret_fish(
        her2, cep17, ihc_score=2,
        recount_her2_copy_number=recount_h,
        recount_cep17_copy_number=recount_c,
    )
    assert result["fish_status"] == expected


def test_group2_not_automatically_positive():
    assessment = assess_her2_status(2, her2_copy_number=3, cep17_copy_number=1)
    assert assessment["final_status"] == "Pending recount"


def test_group3_not_automatically_positive_without_ihc_or_recount():
    assert interpret_fish(7, 4)["fish_status"] == "Pending IHC"
    assert assess_her2_status(2, her2_copy_number=7, cep17_copy_number=4)["final_status"] == "Pending recount"


def test_group4_2023_follows_2018():
    result = interpret_fish(5, 3, guideline_year=2023, ihc_score=2,
                            recount_her2_copy_number=5.1, recount_cep17_copy_number=3.0)
    assert result["fish_status"] == "Negative"


def test_thresholds_use_unrounded_values():
    # Displayed 2.000 MUST NOT accidentally cross 2.0 threshold.
    result = interpret_fish(3.999999, 2.0)
    assert result["fish_group"] == "Group 4"
    assert interpret_fish(4.0, 2.0)["fish_group"] == "Group 1"
    assert interpret_fish(6.0, 4.0)["fish_group"] == "Group 3"


@pytest.mark.parametrize("a,b", [(float("nan"), 2), (float("inf"), 2),
                                  (2, float("nan")), (2, float("inf")),
                                  (-1, 2), (3, 0)])
def test_nonfinite_or_invalid_copy_numbers(a, b):
    with pytest.raises(ValueError):
        interpret_fish(a, b)


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), -1, 101])
def test_invalid_ihc_percent(bad):
    with pytest.raises(ValueError):
        interpret_ihc(2, bad)


def test_inconsistent_ratio_rejected():
    with pytest.raises(ValueError, match="disagrees"):
        interpret_fish(10.0, 2.0, her2_ceph_ratio=2.0)


def test_incomplete_fish_data_rejected():
    with pytest.raises(ValueError, match="Both HER2 and CEP17"):
        assess_her2_status(2, her2_copy_number=5.0)


def test_csv_requires_ihc_and_is_atomic_on_invalid_row(tmp_path):
    source, output = tmp_path / "cases.csv", tmp_path / "results.csv"
    source.write_text("ihc_score,percent_staining\n2,40\n,10\n", encoding="utf-8")
    with pytest.raises(ValueError, match="CSV line 3"):
        process_batch(str(source), str(output))
    assert not output.exists()


def test_csv_formula_injection_escaped(tmp_path):
    source, output = tmp_path / "cases.csv", tmp_path / "results.csv"
    source.write_text("ihc_score,case_note\n0,=HYPERLINK(42)\n", encoding="utf-8")
    process_batch(str(source), str(output))
    with output.open(encoding="utf-8", newline="") as handle:
        row = next(csv.DictReader(handle))
    assert row["case_note"].startswith("'=")


def test_ihc3_and_ish5_discordance_flag():
    res = assess_her2_status(3, her2_copy_number=2, cep17_copy_number=2)
    assert res["final_status"] == "Positive"
    assert "discordant" in res["notes"].lower()
