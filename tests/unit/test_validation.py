"""Tests for validation module."""
import pytest
from main_workflow.reporting.core.validation import (
    validate_threshold_values,
    validate_gene_count,
    validate_custom_gene_list,
    validate_contrasts_selected
)


def test_validate_threshold_values_valid():
    """Test validation with valid threshold values."""
    lfc_val, pval_val, has_error = validate_threshold_values("1.5", "0.01")

    assert lfc_val == 1.5
    assert pval_val == 0.01
    assert has_error is False


def test_validate_threshold_values_invalid():
    """Test validation with invalid threshold values."""
    lfc_val, pval_val, has_error = validate_threshold_values("invalid", "0.05")

    assert lfc_val == 1.0  # Default
    assert pval_val == 0.05  # Default
    assert has_error is True


def test_validate_threshold_values_partial_invalid():
    """Test validation when one threshold is invalid."""
    lfc_val, pval_val, has_error = validate_threshold_values("2.0", "not_a_number")

    assert lfc_val == 1.0  # Default
    assert pval_val == 0.05  # Default
    assert has_error is True


def test_validate_threshold_values_negative():
    """Test validation allows negative values."""
    lfc_val, pval_val, has_error = validate_threshold_values("-1.5", "0.05")

    assert lfc_val == -1.5
    assert pval_val == 0.05
    assert has_error is False


def test_validate_gene_count_valid():
    """Test gene count validation with valid input."""
    gene_count, has_error = validate_gene_count("50")

    assert gene_count == 50
    assert has_error is False


def test_validate_gene_count_invalid():
    """Test gene count validation with invalid input."""
    gene_count, has_error = validate_gene_count("not_a_number")

    assert gene_count == 50  # Default
    assert has_error is True


def test_validate_gene_count_zero():
    """Test gene count validation with zero."""
    gene_count, has_error = validate_gene_count("0")

    assert gene_count == 50  # Default
    assert has_error is True


def test_validate_gene_count_negative():
    """Test gene count validation with negative number."""
    gene_count, has_error = validate_gene_count("-10")

    assert gene_count == 50  # Default
    assert has_error is True


def test_validate_gene_count_large():
    """Test gene count validation with large number."""
    gene_count, has_error = validate_gene_count("10000")

    assert gene_count == 10000
    assert has_error is False


def test_validate_custom_gene_list_valid():
    """Test custom gene list validation with non-empty list."""
    result = validate_custom_gene_list(['GENE1', 'GENE2', 'GENE3'])

    assert result is True


def test_validate_custom_gene_list_empty():
    """Test custom gene list validation with empty list."""
    result = validate_custom_gene_list([])

    assert result is False


def test_validate_custom_gene_list_single_gene():
    """Test custom gene list validation with single gene."""
    result = validate_custom_gene_list(['GENE1'])

    assert result is True


def test_validate_contrasts_selected_valid():
    """Test contrasts validation with selected contrasts."""
    result = validate_contrasts_selected([('GSE123', 'contrast1'), ('GSE456', 'contrast2')])

    assert result is True


def test_validate_contrasts_selected_empty():
    """Test contrasts validation with no contrasts."""
    result = validate_contrasts_selected([])

    assert result is False


def test_validate_contrasts_selected_single():
    """Test contrasts validation with single contrast."""
    result = validate_contrasts_selected([('GSE123', 'contrast1')])

    assert result is True