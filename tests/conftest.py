"""Pytest fixtures for UORCA reporting tests."""
import pytest
import pandas as pd
import time
from contextlib import contextmanager
from typing import Dict, Any


@pytest.fixture
def sample_deg_data():
    """Sample DEG data for testing."""
    return pd.DataFrame({
        'Gene': ['GENE1', 'GENE2', 'GENE3', 'GENE4', 'GENE5'],
        'logFC': [2.5, -1.8, 0.5, 3.2, -2.1],
        'adj.P.Val': [0.001, 0.003, 0.15, 0.0001, 0.002],
        'AveExpr': [10.5, 8.3, 12.1, 9.8, 11.2]
    })


@pytest.fixture
def sample_cpm_data():
    """Sample CPM data for testing."""
    return pd.DataFrame({
        'Gene': ['GENE1', 'GENE2', 'GENE3'],
        'Sample1': [100.5, 200.3, 50.1],
        'Sample2': [150.2, 180.5, 60.3],
        'Sample3': [120.8, 210.1, 55.7]
    })


@pytest.fixture
def sample_analysis_info():
    """Sample analysis_info structure."""
    return {
        'GSE123456': {
            'accession': 'GSE123456',
            'organism': 'Homo sapiens',
            'analysis_success': True,
            'unique_groups': ['Control', 'Treatment'],
            'contrasts': [
                {
                    'name': 'Treatment_vs_Control',
                    'expression': 'Treatment - Control',
                    'description': 'Treatment compared to Control'
                }
            ]
        },
        'GSE789012': {
            'accession': 'GSE789012',
            'organism': 'Mus musculus',
            'analysis_success': True,
            'unique_groups': ['WT', 'KO'],
            'contrasts': [
                {
                    'name': 'KO_vs_WT',
                    'expression': 'KO - WT',
                    'description': 'Knockout vs Wild Type'
                }
            ]
        }
    }


@pytest.fixture
def sample_deg_data_dict(sample_deg_data):
    """Sample deg_data structure (nested dict)."""
    return {
        'GSE123456': {
            'Treatment_vs_Control': sample_deg_data.copy()
        },
        'GSE789012': {
            'KO_vs_WT': sample_deg_data.copy()
        }
    }


@pytest.fixture
def benchmark_timer():
    """Context manager for timing code execution."""
    @contextmanager
    def timer():
        times = {}
        start = time.perf_counter()
        yield times
        times['elapsed'] = time.perf_counter() - start
    return timer


@pytest.fixture
def performance_threshold():
    """Performance degradation threshold (5%)."""
    return 1.05  # New code must be within 5% of original