"""Tests for organism_utils module."""
import pytest
from core.organism_utils import (
    group_datasets_by_organism,
    group_contrasts_by_organism,
    filter_genes_by_organism,
    filter_genes_by_organism_datasets,
    get_organism_display_name,
    get_organisms_from_datasets
)


def test_group_datasets_by_organism(sample_analysis_info):
    """Test grouping datasets by organism."""
    selected = ['GSE123456', 'GSE789012']
    result = group_datasets_by_organism(sample_analysis_info, selected)

    assert len(result) == 2
    assert 'Homo sapiens' in result
    assert 'Mus musculus' in result
    assert result['Homo sapiens'] == ['GSE123456']
    assert result['Mus musculus'] == ['GSE789012']


def test_group_datasets_by_organism_with_unknown(sample_analysis_info):
    """Test grouping with dataset not in analysis_info."""
    selected = ['GSE123456', 'GSE_UNKNOWN']
    result = group_datasets_by_organism(sample_analysis_info, selected)

    assert len(result) == 1
    assert 'Homo sapiens' in result
    assert 'GSE_UNKNOWN' not in str(result)


def test_group_contrasts_by_organism(sample_analysis_info):
    """Test grouping contrasts by organism."""
    contrasts = [
        ('GSE123456', 'Treatment_vs_Control'),
        ('GSE789012', 'KO_vs_WT')
    ]
    result = group_contrasts_by_organism(sample_analysis_info, contrasts)

    assert len(result) == 2
    assert result['Homo sapiens'] == [('GSE123456', 'Treatment_vs_Control')]
    assert result['Mus musculus'] == [('GSE789012', 'KO_vs_WT')]


def test_filter_genes_by_organism(sample_deg_data_dict, sample_analysis_info):
    """Test filtering genes by organism."""
    genes = ['GENE1', 'GENE2', 'GENE_MISSING']
    contrasts = [('GSE123456', 'Treatment_vs_Control')]

    result = filter_genes_by_organism(
        deg_data=sample_deg_data_dict,
        cpm_data={},
        analysis_info=sample_analysis_info,
        genes=genes,
        organism='Homo sapiens',
        selected_contrasts=contrasts
    )

    assert 'GENE1' in result
    assert 'GENE2' in result
    assert 'GENE_MISSING' not in result


def test_filter_genes_by_organism_with_cpm(sample_deg_data_dict, sample_cpm_data, sample_analysis_info):
    """Test filtering genes considering CPM data."""
    genes = ['GENE1', 'GENE2', 'GENE3', 'GENE_MISSING']
    contrasts = [('GSE123456', 'Treatment_vs_Control')]

    result = filter_genes_by_organism(
        deg_data=sample_deg_data_dict,
        cpm_data={'GSE123456': sample_cpm_data},
        analysis_info=sample_analysis_info,
        genes=genes,
        organism='Homo sapiens',
        selected_contrasts=contrasts
    )

    # Should find genes from both DEG and CPM data
    assert 'GENE1' in result
    assert 'GENE2' in result
    assert 'GENE3' in result
    assert 'GENE_MISSING' not in result


def test_filter_genes_by_organism_datasets(sample_cpm_data, sample_analysis_info):
    """Test filtering genes by organism for datasets."""
    genes = ['GENE1', 'GENE2', 'GENE_MISSING']
    datasets = ['GSE123456']

    result = filter_genes_by_organism_datasets(
        cpm_data={'GSE123456': sample_cpm_data},
        analysis_info=sample_analysis_info,
        genes=genes,
        organism='Homo sapiens',
        selected_datasets=datasets
    )

    assert 'GENE1' in result
    assert 'GENE2' in result
    assert 'GENE_MISSING' not in result


def test_get_organism_display_name():
    """Test organism display name formatting."""
    assert get_organism_display_name('homo sapiens') == 'Homo Sapiens'
    assert get_organism_display_name('mus musculus') == 'Mus Musculus'
    assert get_organism_display_name('Unknown') == 'Unknown Species'
    assert get_organism_display_name('') == 'Unknown Species'
    assert get_organism_display_name('Some Long Species Name Here') == 'Some Long Species Name Here'


def test_get_organisms_from_datasets(sample_analysis_info):
    """Test getting unique organisms from datasets."""
    selected = ['GSE123456', 'GSE789012']
    result = get_organisms_from_datasets(sample_analysis_info, selected)

    assert len(result) == 2
    assert 'Homo sapiens' in result
    assert 'Mus musculus' in result
    # Should be sorted
    assert result == sorted(result)


def test_get_organisms_from_datasets_empty():
    """Test with no datasets."""
    result = get_organisms_from_datasets({}, [])
    assert result == []


def test_filter_genes_empty_contrasts():
    """Test filtering with no contrasts."""
    result = filter_genes_by_organism(
        deg_data={},
        cpm_data={},
        analysis_info={},
        genes=['GENE1'],
        organism='Homo sapiens',
        selected_contrasts=[]
    )
    assert result == []