# UORCA Streamlit App Refactoring and Testing Implementation Plan

## Overview

This plan details a comprehensive refactoring of the UORCA Explorer Streamlit application to extract business logic from UI components, establish testing infrastructure, and improve code maintainability. The refactoring will create testable, standalone modules while maintaining all existing functionality.

## Current State Analysis

### Architecture
- **Main Entry**: `main_workflow/reporting/uorca_explorer.py` - 329 lines
- **Tab Modules**: 7 tab files in `main_workflow/reporting/streamlit_tabs/`
- **Core Logic**: `main_workflow/reporting/ResultsIntegration.py` - ~1600 lines
- **Helpers**: `main_workflow/reporting/streamlit_tabs/helpers/__init__.py` - ~1537 lines

### Key Issues Identified
1. **No Testing Infrastructure**: Zero test files, no pytest configuration, no test coverage
2. **Mixed UI and Business Logic**: Data transformations embedded in Streamlit rendering functions
3. **Code Duplication**:
   - Duplicate imports: pandas imported 7 times in ai_assistant_tab.py
   - Duplicate logic: GEO accession sorting repeated 5+ times
   - Similar table creation: 3 different implementations for contrast tables
4. **Unused Code**: ~150-200 lines of empty functions and duplicate imports

### Current Patterns (Working Well)
- **3-Layer Tab Architecture**: render → form → fragment pattern is solid
- **Caching Strategy**: `@st.cache_resource` and `@st.cache_data` effectively used
- **Session State**: Consistent patterns for dataset/contrast selection
- **Organism Handling**: Multi-species support well-implemented

### Key Discoveries

**Already Pure Functions** (just need moving):
- `group_datasets_by_organism()` - helpers/__init__.py:1180
- `filter_genes_by_organism()` - helpers/__init__.py:1224
- `_create_contrast_table_data_filtered()` - heatmap_tab.py:711
- `_create_dataset_info_dataframe()` - datasets_info_tab.py:85
- `_create_contrast_info_dataframe()` - contrasts_info_tab.py:90

**Mixed Functions** (need splitting):
- LFC matrix extraction - ResultsIntegration.py:565-700
- Heatmap clustering logic - ResultsIntegration.py:758-800
- Form validation - heatmap_tab.py:512-686
- Gene availability checking - expression_plots_tab.py:436-462

## Desired End State

### New Directory Structure
```
main_workflow/reporting/
├── business_logic/              # NEW: Pure business logic
│   ├── __init__.py
│   ├── data_formatters.py       # Table creation, sorting
│   ├── gene_selection.py        # Gene frequency analysis
│   ├── organism_utils.py        # Organism grouping/filtering
│   ├── validation.py            # Form validation logic
│   └── heatmap_data.py          # LFC matrix, clustering
├── tests/                       # NEW: Test suite
│   ├── __init__.py
│   ├── conftest.py              # Pytest fixtures
│   ├── unit/
│   │   ├── test_data_formatters.py
│   │   ├── test_gene_selection.py
│   │   ├── test_organism_utils.py
│   │   ├── test_validation.py
│   │   └── test_heatmap_data.py
│   ├── integration/
│   │   ├── test_results_integration.py
│   │   └── test_tab_workflows.py
│   └── fixtures/                # Sample data for tests
│       └── sample_deg_data.csv
├── streamlit_tabs/              # EXISTING: Refactored tabs
│   ├── helpers/
│   │   └── __init__.py          # Slimmed down, UI-only helpers
│   ├── heatmap_tab.py           # Refactored to use business_logic
│   ├── expression_plots_tab.py  # Refactored to use business_logic
│   ├── datasets_info_tab.py     # Refactored to use business_logic
│   ├── contrasts_info_tab.py    # Refactored to use business_logic
│   └── sidebar_controls.py      # Refactored to use business_logic
├── ResultsIntegration.py        # Slimmed down, uses business_logic
└── uorca_explorer.py            # Unchanged (just imports)
```

### Success Criteria

#### Automated Verification:
- [ ] All tests pass: `pytest tests/ -v`
- [ ] Test coverage ≥80%: `pytest --cov=main_workflow/reporting/business_logic --cov-report=term-missing`
- [ ] Type checking passes: `pyright main_workflow/reporting/business_logic/`
- [ ] No duplicate code: Manual review confirms consolidation
- [ ] Performance benchmarks pass: New code ≤ 5% slower than original

#### Manual Verification:
- [ ] Streamlit app loads without errors
- [ ] All 7 tabs render correctly
- [ ] Heatmap generation produces identical output to original
- [ ] Expression plots produce identical output to original
- [ ] Dataset/contrast tables show same data
- [ ] Download functionality works (PDF, ZIP)
- [ ] Multi-organism workflows function correctly
- [ ] AI Assistant tab unaffected (excluded from refactor)

## What We're NOT Doing

1. **NOT refactoring AI Assistant tab** - Complex, working, excluded per requirements
2. **NOT changing external APIs** - ResultsIntegrator public methods stay compatible for now
3. **NOT modifying data file formats** - Still reading same CSV/JSON structures
4. **NOT changing Streamlit UI patterns** - Keeping 3-layer architecture
5. **NOT adding new features** - Pure refactoring for testability
6. **NOT worrying about backwards compatibility** - Free to break internal APIs

## Implementation Approach

### Strategy
1. **Test-First for New Modules**: Write tests before extracting logic
2. **Validate Equivalence**: Compare outputs of original vs refactored code
3. **Incremental Extraction**: Move one function at a time, test immediately
4. **Performance Tracking**: Benchmark before/after each phase

### Risk Mitigation
- Keep original functions temporarily during transition
- Use Git branches for each phase
- Validate against real UORCA results at each step
- Performance benchmarks catch regressions

---

## Phase 1: Setup & Quick Wins

### Overview
Establish testing infrastructure, remove duplicate/unused code, and move already-pure functions to new business_logic module. This phase creates the foundation without breaking existing functionality.

### Changes Required

#### 1. Setup Testing Infrastructure

**File**: `pyproject.toml`
**Changes**: Add testing dependencies

```toml
[project.optional-dependencies]
test = [
    "pytest>=8.0.0",
    "pytest-cov>=5.0.0",
    "pytest-mock>=3.14.0",
]
dev = [
    "pyright>=1.1.350",
]
```

**File**: `pytest.ini` (NEW)
**Changes**: Create pytest configuration

```ini
[pytest]
testpaths = main_workflow/reporting/tests
python_files = test_*.py
python_classes = Test*
python_functions = test_*
addopts =
    -v
    --strict-markers
    --tb=short
    --cov=main_workflow/reporting/business_logic
    --cov-report=term-missing
    --cov-report=html
markers =
    unit: Unit tests
    integration: Integration tests
    slow: Slow-running tests
```

#### 2. Create Directory Structure

**Command**: Create new directories
```bash
mkdir -p main_workflow/reporting/business_logic
mkdir -p main_workflow/reporting/tests/{unit,integration,fixtures}
touch main_workflow/reporting/business_logic/__init__.py
touch main_workflow/reporting/tests/{__init__.py,conftest.py}
```

#### 3. Remove Duplicate Imports

**File**: `main_workflow/reporting/streamlit_tabs/ai_assistant_tab.py`
**Changes**: Remove redundant pandas/numpy imports

```python
# REMOVE lines 1404, 1460, 1516, 1554, 1616, 1700:
# import pandas as pd  # Delete these 6 occurrences

# REMOVE lines 1617, 1701:
# import numpy as np  # Delete these 2 occurrences

# Top-level imports (lines 14-15) are sufficient
```

**File**: `main_workflow/reporting/streamlit_tabs/helpers/__init__.py`
**Changes**: Remove duplicate Path import and unused shutil

```python
# REMOVE line 20:
# from pathlib import Path  # Duplicate of line 13

# REMOVE line 17:
# import shutil  # Unused
```

#### 4. Remove Unused Functions

**File**: `main_workflow/reporting/streamlit_tabs/heatmap_tab.py`
**Changes**: Delete empty/unused functions

```python
# DELETE lines 986-989:
# def _display_comprehensive_filtering_info():
#     """This function is now integrated into the form validation and no longer displays separately"""
#     pass

# DELETE lines 1024-1027:
# def _build_repro_script():
#     """Wrapper for backward compatibility but not used"""
#     return _build_repro_script_static()
```

**File**: `main_workflow/reporting/streamlit_tabs/contrasts_info_tab.py`
**Changes**: Delete empty function

```python
# DELETE lines 214-216:
# def _render_selection_controls():
#     """Removed info message about sidebar selection"""
#     pass
```

**File**: `main_workflow/reporting/streamlit_tabs/helpers/__init__.py`
**Changes**: Remove unnecessary wrapper

```python
# DELETE lines 1102-1104:
# def safe_rerun():
#     """Call streamlit rerun."""
#     st.rerun()

# REPLACE all calls to safe_rerun() with st.rerun() directly
```

#### 5. Create Initial Test Fixtures

**File**: `main_workflow/reporting/tests/conftest.py` (NEW)
**Changes**: Create pytest fixtures for common test data

```python
"""Pytest fixtures for UORCA reporting tests."""
import pytest
import pandas as pd
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
```

#### 6. Create Performance Benchmarking Utilities

**File**: `main_workflow/reporting/tests/conftest.py` (APPEND)
**Changes**: Add benchmarking fixtures

```python
import time
from contextlib import contextmanager

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
```

### Success Criteria

#### Automated Verification:
- [ ] pytest runs successfully: `pytest tests/ -v`
- [ ] Test infrastructure validates: `pytest tests/conftest.py -v`
- [ ] No import errors: `python -c "import main_workflow.reporting.business_logic"`
- [ ] Type checking passes: `pyright main_workflow/reporting/streamlit_tabs/`
- [ ] Duplicate imports removed: `grep -r "import pandas as pd" main_workflow/reporting/streamlit_tabs/ai_assistant_tab.py | wc -l` returns 1
- [ ] Unused functions removed: Verify deleted lines no longer exist

#### Manual Verification:
- [ ] Streamlit app still loads: `streamlit run main_workflow/reporting/uorca_explorer.py`
- [ ] All tabs render without errors
- [ ] No visible changes to UI behavior
- [ ] AI Assistant tab still functions (imports unchanged)

---

## Phase 2: Extract Core Data Logic

### Overview
Extract pure data transformation functions from helpers and tabs into new business_logic modules. These functions have no Streamlit dependencies and can be tested independently.

### Changes Required

#### 1. Create Organism Utilities Module

**File**: `main_workflow/reporting/business_logic/organism_utils.py` (NEW)
**Changes**: Extract organism-related logic

```python
"""Organism grouping and filtering utilities."""
from typing import Dict, List, Set
import pandas as pd

def group_datasets_by_organism(
    analysis_info: Dict[str, Dict],
    selected_datasets: List[str]
) -> Dict[str, List[str]]:
    """
    Group datasets by their organism.

    Args:
        analysis_info: Analysis metadata dict (from ResultsIntegrator)
        selected_datasets: List of dataset IDs

    Returns:
        Dict mapping organism names to lists of dataset IDs
    """
    organism_groups = {}
    for dataset_id in selected_datasets:
        if dataset_id in analysis_info:
            organism = analysis_info[dataset_id].get('organism', 'Unknown')
            if organism not in organism_groups:
                organism_groups[organism] = []
            organism_groups[organism].append(dataset_id)
    return organism_groups


def group_contrasts_by_organism(
    analysis_info: Dict[str, Dict],
    selected_contrasts: List[tuple]
) -> Dict[str, List[tuple]]:
    """
    Group contrasts by their dataset's organism.

    Args:
        analysis_info: Analysis metadata dict
        selected_contrasts: List of (analysis_id, contrast_id) tuples

    Returns:
        Dict mapping organism names to lists of contrast tuples
    """
    organism_groups = {}
    for analysis_id, contrast_id in selected_contrasts:
        if analysis_id in analysis_info:
            organism = analysis_info[analysis_id].get('organism', 'Unknown')
            if organism not in organism_groups:
                organism_groups[organism] = []
            organism_groups[organism].append((analysis_id, contrast_id))
    return organism_groups


def filter_genes_by_organism(
    deg_data: Dict[str, Dict[str, pd.DataFrame]],
    cpm_data: Dict[str, pd.DataFrame],
    analysis_info: Dict[str, Dict],
    genes: List[str],
    organism: str,
    selected_contrasts: List[tuple]
) -> List[str]:
    """
    Filter genes to only include those found in datasets of a specific organism.

    Args:
        deg_data: DEG data dict (analysis_id -> contrast_id -> DataFrame)
        cpm_data: CPM data dict (analysis_id -> DataFrame)
        analysis_info: Analysis metadata dict
        genes: List of gene names to filter
        organism: Target organism name
        selected_contrasts: List of (analysis_id, contrast_id) tuples

    Returns:
        List of genes present in the organism's datasets
    """
    organism_genes: Set[str] = set()

    # Collect all genes from datasets of this organism
    for analysis_id, contrast_id in selected_contrasts:
        if analysis_id in analysis_info and analysis_info[analysis_id].get('organism') == organism:
            # Check DEG data
            if analysis_id in deg_data and contrast_id in deg_data[analysis_id]:
                deg_df = deg_data[analysis_id][contrast_id]
                if 'Gene' in deg_df.columns:
                    organism_genes.update(deg_df['Gene'].tolist())

            # Check CPM data
            if analysis_id in cpm_data:
                cpm_df = cpm_data[analysis_id]
                if 'Gene' in cpm_df.columns:
                    organism_genes.update(cpm_df['Gene'].tolist())

    # Return only genes in input list AND found in organism's data
    return [gene for gene in genes if gene in organism_genes]


def get_organism_display_name(organism: str) -> str:
    """
    Get a user-friendly display name for an organism.

    Args:
        organism: Scientific or common name of organism

    Returns:
        Cleaned display name for UI
    """
    if not organism or organism == 'Unknown':
        return 'Unknown Species'

    organism = organism.strip()

    # Capitalize first letter of each word for scientific names
    if len(organism.split()) <= 2:
        return organism.title()
    else:
        return organism
```

**File**: `main_workflow/reporting/tests/unit/test_organism_utils.py` (NEW)
**Changes**: Create tests

```python
"""Tests for organism_utils module."""
import pytest
from business_logic.organism_utils import (
    group_datasets_by_organism,
    group_contrasts_by_organism,
    filter_genes_by_organism,
    get_organism_display_name
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

def test_get_organism_display_name():
    """Test organism display name formatting."""
    assert get_organism_display_name('homo sapiens') == 'Homo Sapiens'
    assert get_organism_display_name('Unknown') == 'Unknown Species'
    assert get_organism_display_name('') == 'Unknown Species'
```

#### 2. Create Data Formatters Module

**File**: `main_workflow/reporting/business_logic/data_formatters.py` (NEW)
**Changes**: Extract table creation and sorting logic

```python
"""Data formatting utilities for UI display."""
import re
import pandas as pd
from typing import List, Dict, Any, Tuple

def sort_by_geo_accession(df: pd.DataFrame, accession_col: str = "Accession") -> pd.DataFrame:
    """
    Sort DataFrame by GEO accession number.

    Args:
        df: DataFrame to sort
        accession_col: Name of accession column

    Returns:
        Sorted DataFrame
    """
    df = df.copy()
    try:
        df["_AccessionNum"] = (
            df[accession_col]
            .astype(str)
            .str.extract(r"(\d+)", expand=False)
            .astype(int)
        )
    except Exception:
        df["_AccessionNum"] = float("inf")

    return df.sort_values(
        ["_AccessionNum", accession_col],
        ascending=[True, True]
    ).drop(columns=["_AccessionNum"], errors="ignore")


def create_contrast_table_data(
    analysis_info: Dict[str, Dict],
    deg_data: Dict[str, Dict[str, pd.DataFrame]],
    contrast_info: Dict[str, Dict],
    selected_datasets: List[str]
) -> List[Dict[str, Any]]:
    """
    Create contrast table data using standardized validation logic.

    Args:
        analysis_info: Analysis metadata
        deg_data: DEG data
        contrast_info: Contrast metadata
        selected_datasets: List of dataset IDs to include

    Returns:
        List of dicts for contrast table display
    """
    valid_contrasts = []

    for analysis_id in selected_datasets:
        # Check if analysis was successful
        if analysis_id not in analysis_info:
            continue
        info = analysis_info[analysis_id]
        if not info.get("analysis_success", False):
            continue

        # Get contrasts from analysis_info
        contrasts = info.get("contrasts", [])
        accession = info.get("accession", analysis_id)

        for contrast in contrasts:
            original_name = contrast["name"]

            # Validate contrast has DEG data
            if (analysis_id in deg_data and
                original_name in deg_data[analysis_id] and
                not deg_data[analysis_id][original_name].empty):

                # Find consistent name
                consistent_name = original_name
                for contrast_key, contrast_data in contrast_info.items():
                    if (contrast_data.get('original_name') == original_name and
                        contrast_data.get('analysis_id') == analysis_id):
                        consistent_name = contrast_data.get('name', original_name)
                        break

                valid_contrasts.append({
                    "Select": False,
                    "Accession": accession,
                    "Contrast": consistent_name,
                    "Description": contrast.get("description", "")
                })

    # Sort by accession
    df = pd.DataFrame(valid_contrasts)
    if not df.empty:
        df = sort_by_geo_accession(df, "Accession")
        return df.to_dict('records')
    return []


def create_dataset_info_table(
    analysis_info: Dict[str, Dict],
    dataset_info: Dict[str, Dict],
    deg_data: Dict[str, Dict]
) -> pd.DataFrame:
    """
    Create dataset information table.

    Args:
        analysis_info: Analysis metadata
        dataset_info: Dataset metadata (titles, descriptions)
        deg_data: DEG data for counting contrasts

    Returns:
        DataFrame with dataset information
    """
    rows = []

    for analysis_id, info in analysis_info.items():
        # Only include successful analyses
        if not info.get("analysis_success", False):
            continue

        accession = info.get("accession", analysis_id)
        organism = info.get("organism", "Unknown")

        # Count samples
        unique_groups = info.get("unique_groups", [])
        sample_count = len(unique_groups) if unique_groups else 0

        # Count contrasts
        contrast_count = len(deg_data.get(analysis_id, {}))

        # Get title/summary
        title = ""
        summary = ""
        if analysis_id in dataset_info:
            title = dataset_info[analysis_id].get("title", "")
            if title.startswith("Title:"):
                title = title[6:].strip()
            summary = dataset_info[analysis_id].get("summary", "")
            if summary.startswith("Summary:"):
                summary = summary[8:].strip()

        rows.append({
            "Accession": accession,
            "Organism": organism,
            "Samples": sample_count,
            "Contrasts": contrast_count,
            "Title": title,
            "Description": summary
        })

    df = pd.DataFrame(rows)
    if not df.empty:
        df = sort_by_geo_accession(df, "Accession")
    return df


def create_contrast_info_table(
    analysis_info: Dict[str, Dict],
    deg_data: Dict[str, Dict[str, pd.DataFrame]],
    pvalue_thresh: float,
    lfc_thresh: float
) -> pd.DataFrame:
    """
    Create contrast information table with DEG counts.

    Args:
        analysis_info: Analysis metadata
        deg_data: DEG data
        pvalue_thresh: P-value threshold for significance
        lfc_thresh: Log fold-change threshold

    Returns:
        DataFrame with contrast information
    """
    rows = []

    for analysis_id, contrasts_dict in deg_data.items():
        if analysis_id not in analysis_info:
            continue

        info = analysis_info[analysis_id]
        accession = info.get("accession", analysis_id)

        for contrast_id, deg_df in contrasts_dict.items():
            # Count significant DEGs
            if 'adj.P.Val' in deg_df.columns and 'logFC' in deg_df.columns:
                sig_df = deg_df[
                    (deg_df['adj.P.Val'] < pvalue_thresh) &
                    (abs(deg_df['logFC']) > lfc_thresh)
                ]
                n_sig = len(sig_df)
                n_total = len(deg_df)
            else:
                n_sig = 0
                n_total = len(deg_df)

            # Get description from analysis_info
            description = ""
            for contrast in info.get("contrasts", []):
                if contrast.get("name") == contrast_id:
                    description = contrast.get("description", "")
                    break

            rows.append({
                "Accession": accession,
                "Contrast": contrast_id,
                "Description": description,
                "Significant DEGs": n_sig,
                "Total Genes": n_total
            })

    df = pd.DataFrame(rows)
    if not df.empty:
        df = sort_by_geo_accession(df, "Accession")
    return df
```

**File**: `main_workflow/reporting/tests/unit/test_data_formatters.py` (NEW)
**Changes**: Create tests

```python
"""Tests for data_formatters module."""
import pytest
import pandas as pd
from business_logic.data_formatters import (
    sort_by_geo_accession,
    create_contrast_table_data,
    create_dataset_info_table,
    create_contrast_info_table
)

def test_sort_by_geo_accession():
    """Test GEO accession sorting."""
    df = pd.DataFrame({
        "Accession": ["GSE100", "GSE50", "GSE200", "GSE10"],
        "Data": ["a", "b", "c", "d"]
    })

    result = sort_by_geo_accession(df)

    assert result["Accession"].tolist() == ["GSE10", "GSE50", "GSE100", "GSE200"]

def test_create_contrast_table_data(sample_analysis_info, sample_deg_data_dict):
    """Test contrast table creation."""
    result = create_contrast_table_data(
        analysis_info=sample_analysis_info,
        deg_data=sample_deg_data_dict,
        contrast_info={},
        selected_datasets=['GSE123456', 'GSE789012']
    )

    assert len(result) == 2
    assert all('Accession' in row for row in result)
    assert all('Contrast' in row for row in result)
    assert all('Select' in row for row in result)

def test_create_dataset_info_table(sample_analysis_info, sample_deg_data_dict):
    """Test dataset info table creation."""
    result = create_dataset_info_table(
        analysis_info=sample_analysis_info,
        dataset_info={},
        deg_data=sample_deg_data_dict
    )

    assert len(result) == 2
    assert 'Accession' in result.columns
    assert 'Organism' in result.columns
    assert 'Contrasts' in result.columns

def test_create_contrast_info_table(sample_analysis_info, sample_deg_data_dict):
    """Test contrast info table with DEG counts."""
    result = create_contrast_info_table(
        analysis_info=sample_analysis_info,
        deg_data=sample_deg_data_dict,
        pvalue_thresh=0.05,
        lfc_thresh=1.0
    )

    assert not result.empty
    assert 'Significant DEGs' in result.columns
    assert 'Total Genes' in result.columns
```

#### 3. Create Gene Selection Module

**File**: `main_workflow/reporting/business_logic/gene_selection.py` (NEW)
**Changes**: Extract gene identification logic

```python
"""Gene selection and identification utilities."""
import pandas as pd
from typing import Dict, List, Tuple, Set
import logging

logger = logging.getLogger(__name__)

def identify_frequent_degs(
    deg_data: Dict[str, Dict[str, pd.DataFrame]],
    contrasts: List[Tuple[str, str]],
    top_n: int,
    p_thresh: float,
    lfc_thresh: float
) -> List[str]:
    """
    Identify genes frequently differentially expressed across contrasts.

    Args:
        deg_data: DEG data dict (analysis_id -> contrast_id -> DataFrame)
        contrasts: List of (analysis_id, contrast_id) tuples to analyze
        top_n: Number of top genes to return
        p_thresh: P-value threshold for significance
        lfc_thresh: Log fold-change threshold

    Returns:
        List of gene names sorted by frequency
    """
    gene_counts = {}
    contrasts_processed = 0

    for analysis_id, contrast_id in contrasts:
        if analysis_id in deg_data and contrast_id in deg_data[analysis_id]:
            deg_df = deg_data[analysis_id][contrast_id]

            # Filter for significant genes
            if 'adj.P.Val' in deg_df.columns and 'logFC' in deg_df.columns and 'Gene' in deg_df.columns:
                significant_genes = deg_df[
                    (deg_df['adj.P.Val'] < p_thresh) &
                    (abs(deg_df['logFC']) > lfc_thresh)
                ]['Gene'].tolist()

                contrasts_processed += 1

                # Count occurrences
                for gene in significant_genes:
                    gene_counts[gene] = gene_counts.get(gene, 0) + 1

    # Sort by frequency (descending) then by gene name (ascending)
    sorted_genes = sorted(gene_counts.items(), key=lambda x: (-x[1], x[0]))

    logger.info(f"Processed {contrasts_processed} contrasts, found {len(gene_counts)} unique genes")

    # Return top N genes
    return [gene for gene, count in sorted_genes[:top_n]]


def get_all_genes_from_cpm(cpm_data: Dict[str, pd.DataFrame]) -> List[str]:
    """
    Extract all unique genes from CPM data.

    Args:
        cpm_data: CPM data dict (analysis_id -> DataFrame)

    Returns:
        Sorted list of unique gene names
    """
    all_genes = set()
    for cpm_df in cpm_data.values():
        if 'Gene' in cpm_df.columns:
            all_genes.update(cpm_df['Gene'].tolist())
    return sorted(all_genes)


def get_available_genes_for_contrasts(
    deg_data: Dict[str, Dict[str, pd.DataFrame]],
    contrasts: List[Tuple[str, str]]
) -> Set[str]:
    """
    Get all genes available in selected contrasts.

    Args:
        deg_data: DEG data dict
        contrasts: List of (analysis_id, contrast_id) tuples

    Returns:
        Set of gene names present in at least one contrast
    """
    available_genes = set()

    for analysis_id, contrast_id in contrasts:
        if analysis_id in deg_data and contrast_id in deg_data[analysis_id]:
            deg_df = deg_data[analysis_id][contrast_id]
            if 'Gene' in deg_df.columns:
                available_genes.update(deg_df['Gene'].tolist())

    return available_genes


def validate_custom_genes(
    custom_genes: List[str],
    available_genes: Set[str]
) -> Tuple[List[str], List[str]]:
    """
    Validate custom gene list against available genes.

    Args:
        custom_genes: User-provided gene list
        available_genes: Genes available in the data

    Returns:
        Tuple of (found_genes, missing_genes)
    """
    found = [g for g in custom_genes if g in available_genes]
    missing = [g for g in custom_genes if g not in available_genes]
    return found, missing
```

**File**: `main_workflow/reporting/tests/unit/test_gene_selection.py` (NEW)
**Changes**: Create tests

```python
"""Tests for gene_selection module."""
import pytest
import pandas as pd
from business_logic.gene_selection import (
    identify_frequent_degs,
    get_all_genes_from_cpm,
    get_available_genes_for_contrasts,
    validate_custom_genes
)

def test_identify_frequent_degs(sample_deg_data_dict):
    """Test identifying frequently DE genes."""
    contrasts = [
        ('GSE123456', 'Treatment_vs_Control'),
        ('GSE789012', 'KO_vs_WT')
    ]

    result = identify_frequent_degs(
        deg_data=sample_deg_data_dict,
        contrasts=contrasts,
        top_n=3,
        p_thresh=0.05,
        lfc_thresh=1.0
    )

    assert len(result) <= 3
    assert all(isinstance(gene, str) for gene in result)

def test_get_all_genes_from_cpm(sample_cpm_data):
    """Test extracting all genes from CPM data."""
    cpm_dict = {'dataset1': sample_cpm_data}
    result = get_all_genes_from_cpm(cpm_dict)

    assert len(result) == 3
    assert result == ['GENE1', 'GENE2', 'GENE3']  # Should be sorted

def test_get_available_genes_for_contrasts(sample_deg_data_dict):
    """Test getting available genes for contrasts."""
    contrasts = [('GSE123456', 'Treatment_vs_Control')]
    result = get_available_genes_for_contrasts(sample_deg_data_dict, contrasts)

    assert len(result) == 5
    assert 'GENE1' in result
    assert 'GENE5' in result

def test_validate_custom_genes():
    """Test custom gene validation."""
    custom = ['GENE1', 'GENE2', 'MISSING1', 'MISSING2']
    available = {'GENE1', 'GENE2', 'GENE3'}

    found, missing = validate_custom_genes(custom, available)

    assert found == ['GENE1', 'GENE2']
    assert missing == ['MISSING1', 'MISSING2']
```

#### 4. Update Tabs to Use Business Logic

**File**: `main_workflow/reporting/streamlit_tabs/helpers/__init__.py`
**Changes**: Update imports and delegate to business_logic

```python
# Add import at top
from business_logic import organism_utils, data_formatters, gene_selection

# Replace existing functions with wrappers
def group_datasets_by_organism(ri: ResultsIntegrator, selected_datasets: List[str]) -> Dict[str, List[str]]:
    """Group datasets by their organism."""
    return organism_utils.group_datasets_by_organism(ri.analysis_info, selected_datasets)

def filter_genes_by_organism(ri: ResultsIntegrator, genes: List[str], organism: str, selected_contrasts: List[Tuple[str, str]]) -> List[str]:
    """Filter genes to only include those found in datasets of a specific organism."""
    return organism_utils.filter_genes_by_organism(
        ri.deg_data, ri.cpm_data, ri.analysis_info, genes, organism, selected_contrasts
    )

# Similarly update other helper functions to delegate to business_logic modules
```

**File**: `main_workflow/reporting/streamlit_tabs/datasets_info_tab.py`
**Changes**: Use new data_formatters module

```python
# Update imports
from business_logic.data_formatters import create_dataset_info_table

# Replace _create_dataset_info_dataframe function (line 85)
def _create_dataset_info_dataframe(ri: ResultsIntegrator) -> pd.DataFrame:
    """Create dataset info dataframe using business logic."""
    return create_dataset_info_table(
        analysis_info=ri.analysis_info,
        dataset_info=ri.dataset_info,
        deg_data=ri.deg_data
    )
```

**File**: `main_workflow/reporting/streamlit_tabs/contrasts_info_tab.py`
**Changes**: Use new data_formatters module

```python
# Update imports
from business_logic.data_formatters import create_contrast_info_table

# Replace _create_contrast_info_dataframe function (line 90)
def _create_contrast_info_dataframe(ri: ResultsIntegrator, pvalue_thresh: float, lfc_thresh: float) -> pd.DataFrame:
    """Create contrast info dataframe using business logic."""
    return create_contrast_info_table(
        analysis_info=ri.analysis_info,
        deg_data=ri.deg_data,
        pvalue_thresh=pvalue_thresh,
        lfc_thresh=lfc_thresh
    )
```

### Success Criteria

#### Automated Verification:
- [ ] All unit tests pass: `pytest tests/unit/ -v`
- [ ] Test coverage ≥80% for new modules: `pytest --cov=main_workflow/reporting/business_logic --cov-report=term-missing`
- [ ] Type checking passes: `pyright main_workflow/reporting/business_logic/`
- [ ] No import errors in tabs: `python -m py_compile main_workflow/reporting/streamlit_tabs/*.py`
- [ ] Performance benchmark: New functions within 5% of original timing

#### Manual Verification:
- [ ] Dataset info tab shows same table as before refactor
- [ ] Contrast info tab shows same DEG counts
- [ ] Multi-organism workflows still work correctly
- [ ] No visual changes to any UI elements

---

## Phase 3: Extract Heatmap and Validation Logic

### Overview
Extract complex data preparation logic from ResultsIntegration.py and validation logic from form functions. Split mixed UI/logic functions into separate concerns.

### Changes Required

#### 1. Create Heatmap Data Module

**File**: `main_workflow/reporting/business_logic/heatmap_data.py` (NEW)
**Changes**: Extract LFC matrix and clustering logic

```python
"""Heatmap data preparation and clustering utilities."""
import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist

def extract_lfc_matrix(
    deg_data: Dict[str, Dict[str, pd.DataFrame]],
    genes: List[str],
    contrasts: List[Tuple[str, str]],
    contrast_labels: List[str],
    p_thresh: float,
    lfc_thresh: float
) -> pd.DataFrame:
    """
    Extract log fold-change matrix with significance filtering.

    Args:
        deg_data: DEG data dict (analysis_id -> contrast_id -> DataFrame)
        genes: List of gene names to extract
        contrasts: List of (analysis_id, contrast_id) tuples
        contrast_labels: Display labels for contrasts
        p_thresh: P-value threshold for significance
        lfc_thresh: Absolute LFC threshold for significance

    Returns:
        DataFrame with genes as rows, contrasts as columns, LFC values
        (values below threshold set to 0)
    """
    heatmap_data = []

    for gene in genes:
        row = {'Gene': gene}

        for (analysis_id, contrast_id), contrast_label in zip(contrasts, contrast_labels):
            if analysis_id not in deg_data or contrast_id not in deg_data[analysis_id]:
                row[contrast_label] = 0
                continue

            df = deg_data[analysis_id][contrast_id]

            # Find LFC column
            lfc_col = None
            for col in ['logFC', 'log2FoldChange', 'log2FC', 'LogFC']:
                if col in df.columns:
                    lfc_col = col
                    break

            if lfc_col is None:
                row[contrast_label] = 0
                continue

            # Find p-value column
            p_value_col = None
            if 'adj.P.Val' in df.columns:
                p_value_col = 'adj.P.Val'
            elif 'P.Value' in df.columns:
                p_value_col = 'P.Value'

            if p_value_col is None:
                row[contrast_label] = 0
                continue

            # Get value for this gene
            gene_row = df[df['Gene'] == gene]
            if not gene_row.empty:
                p_value = gene_row.iloc[0][p_value_col]
                abs_lfc = abs(gene_row.iloc[0][lfc_col])

                # Apply significance thresholds
                if p_value >= p_thresh or abs_lfc <= lfc_thresh:
                    row[contrast_label] = 0
                else:
                    row[contrast_label] = gene_row.iloc[0][lfc_col]
            else:
                row[contrast_label] = 0

        heatmap_data.append(row)

    return pd.DataFrame(heatmap_data)


def cluster_heatmap_data(
    df: pd.DataFrame,
    cluster_genes: bool = True,
    cluster_contrasts: bool = True
) -> Tuple[List[str], List[str]]:
    """
    Perform hierarchical clustering on heatmap data.

    Args:
        df: Heatmap DataFrame (genes x contrasts)
        cluster_genes: Whether to cluster genes (rows)
        cluster_contrasts: Whether to cluster contrasts (columns)

    Returns:
        Tuple of (ordered_genes, ordered_contrasts)
    """
    gene_col = 'Gene'
    genes = df[gene_col].tolist()
    contrasts = [col for col in df.columns if col != gene_col]

    # Default: no clustering
    ordered_genes = genes.copy()
    ordered_contrasts = contrasts.copy()

    # Cluster genes (rows)
    if cluster_genes and len(genes) > 1:
        try:
            data_matrix = df[contrasts].values
            # Use correlation distance
            gene_linkage = linkage(data_matrix, method='average', metric='correlation')
            gene_dendrogram = dendrogram(gene_linkage, no_plot=True)
            gene_order = gene_dendrogram['leaves']
            ordered_genes = [genes[i] for i in gene_order]
        except Exception:
            # Fallback to original order
            pass

    # Cluster contrasts (columns)
    if cluster_contrasts and len(contrasts) > 1:
        try:
            data_matrix = df[contrasts].values.T
            contrast_linkage = linkage(data_matrix, method='average', metric='correlation')
            contrast_dendrogram = dendrogram(contrast_linkage, no_plot=True)
            contrast_order = contrast_dendrogram['leaves']
            ordered_contrasts = [contrasts[i] for i in contrast_order]
        except Exception:
            # Fallback to original order
            pass

    return ordered_genes, ordered_contrasts


def simplify_contrast_labels(
    contrasts: List[Tuple[str, str]],
    analysis_info: Dict[str, Dict]
) -> List[str]:
    """
    Create simplified labels for contrasts.

    Args:
        contrasts: List of (analysis_id, contrast_id) tuples
        analysis_info: Analysis metadata

    Returns:
        List of simplified labels
    """
    accession_counts = {}
    for analysis_id, _ in contrasts:
        if analysis_id in analysis_info:
            accession = analysis_info[analysis_id].get('accession', analysis_id)
            accession_counts[accession] = accession_counts.get(accession, 0) + 1

    labels = []
    for analysis_id, contrast_id in contrasts:
        if analysis_id in analysis_info:
            accession = analysis_info[analysis_id].get('accession', analysis_id)
            if accession_counts[accession] > 1:
                labels.append(f"{accession}:{contrast_id}")
            else:
                labels.append(contrast_id)
        else:
            labels.append(f"{analysis_id}:{contrast_id}")

    return labels
```

**File**: `main_workflow/reporting/tests/unit/test_heatmap_data.py` (NEW)
**Changes**: Create tests

```python
"""Tests for heatmap_data module."""
import pytest
import pandas as pd
from business_logic.heatmap_data import (
    extract_lfc_matrix,
    cluster_heatmap_data,
    simplify_contrast_labels
)

def test_extract_lfc_matrix(sample_deg_data_dict):
    """Test LFC matrix extraction."""
    genes = ['GENE1', 'GENE2']
    contrasts = [('GSE123456', 'Treatment_vs_Control')]
    labels = ['Test Contrast']

    result = extract_lfc_matrix(
        deg_data=sample_deg_data_dict,
        genes=genes,
        contrasts=contrasts,
        contrast_labels=labels,
        p_thresh=0.05,
        lfc_thresh=1.0
    )

    assert 'Gene' in result.columns
    assert 'Test Contrast' in result.columns
    assert len(result) == 2

def test_cluster_heatmap_data():
    """Test clustering logic."""
    df = pd.DataFrame({
        'Gene': ['GENE1', 'GENE2', 'GENE3'],
        'Contrast1': [2.5, 1.2, -1.5],
        'Contrast2': [2.3, 1.1, -1.4]
    })

    ordered_genes, ordered_contrasts = cluster_heatmap_data(
        df, cluster_genes=True, cluster_contrasts=True
    )

    assert len(ordered_genes) == 3
    assert len(ordered_contrasts) == 2
    assert set(ordered_genes) == {'GENE1', 'GENE2', 'GENE3'}

def test_simplify_contrast_labels(sample_analysis_info):
    """Test contrast label simplification."""
    contrasts = [
        ('GSE123456', 'Treatment_vs_Control'),
        ('GSE789012', 'KO_vs_WT')
    ]

    result = simplify_contrast_labels(contrasts, sample_analysis_info)

    assert len(result) == 2
    # Single contrast per dataset should use short labels
    assert 'GSE' not in result[0] or ':' not in result[0]
```

#### 2. Create Validation Module

**File**: `main_workflow/reporting/business_logic/validation.py` (NEW)
**Changes**: Extract form validation logic

```python
"""Form validation utilities."""
from typing import Dict, List, Tuple, Optional, Any

class ValidationError(Exception):
    """Custom exception for validation failures."""
    pass

def validate_numeric_thresholds(
    lfc_thresh: str,
    pvalue_thresh: str
) -> Tuple[float, float]:
    """
    Validate and parse numeric thresholds.

    Args:
        lfc_thresh: LFC threshold as string
        pvalue_thresh: P-value threshold as string

    Returns:
        Tuple of (lfc_float, pval_float)

    Raises:
        ValidationError: If values cannot be parsed
    """
    try:
        lfc_val = float(lfc_thresh)
        pval_val = float(pvalue_thresh)

        if lfc_val < 0:
            raise ValidationError("LFC threshold must be positive")
        if not (0 < pval_val <= 1):
            raise ValidationError("P-value threshold must be between 0 and 1")

        return lfc_val, pval_val
    except ValueError:
        raise ValidationError("Please enter valid numeric values for thresholds")


def validate_gene_count(gene_count_input: str) -> int:
    """
    Validate and parse gene count input.

    Args:
        gene_count_input: Gene count as string

    Returns:
        Gene count as integer

    Raises:
        ValidationError: If value invalid
    """
    try:
        count = int(gene_count_input)
        if count <= 0:
            raise ValidationError("Gene count must be positive")
        return count
    except ValueError:
        raise ValidationError("Please enter a valid positive integer for gene count")


def validate_heatmap_params(
    gene_selection_method: str,
    gene_count_input: str,
    custom_genes_list: List[str],
    lfc_thresh: str,
    pvalue_thresh: str,
    selected_contrasts: List[Tuple[str, str]],
    available_genes: set
) -> Dict[str, Any]:
    """
    Validate all heatmap form parameters.

    Args:
        gene_selection_method: "Frequent DEGs" or "Custom"
        gene_count_input: Number of genes (if Frequent)
        custom_genes_list: Custom gene list (if Custom)
        lfc_thresh: LFC threshold
        pvalue_thresh: P-value threshold
        selected_contrasts: List of selected contrasts
        available_genes: Set of genes available in data

    Returns:
        Dict of validated parameters

    Raises:
        ValidationError: If validation fails
    """
    # Validate contrasts
    if not selected_contrasts:
        raise ValidationError("No contrasts selected")

    # Validate thresholds
    lfc_val, pval_val = validate_numeric_thresholds(lfc_thresh, pvalue_thresh)

    # Validate genes based on method
    if gene_selection_method == "Frequent DEGs":
        gene_count = validate_gene_count(gene_count_input)
        return {
            'method': 'frequent',
            'gene_count': gene_count,
            'lfc_thresh': lfc_val,
            'pval_thresh': pval_val,
            'contrasts': selected_contrasts
        }
    else:  # Custom
        if not custom_genes_list:
            raise ValidationError("No custom genes provided")

        # Check gene availability
        found_genes = [g for g in custom_genes_list if g in available_genes]
        missing_genes = [g for g in custom_genes_list if g not in available_genes]

        if not found_genes:
            raise ValidationError(
                f"None of the provided genes were found in the selected contrasts. "
                f"Missing genes: {', '.join(missing_genes[:10])}"
            )

        return {
            'method': 'custom',
            'genes': found_genes,
            'missing_genes': missing_genes,
            'lfc_thresh': lfc_val,
            'pval_thresh': pval_val,
            'contrasts': selected_contrasts
        }
```

**File**: `main_workflow/reporting/tests/unit/test_validation.py` (NEW)
**Changes**: Create tests

```python
"""Tests for validation module."""
import pytest
from business_logic.validation import (
    validate_numeric_thresholds,
    validate_gene_count,
    validate_heatmap_params,
    ValidationError
)

def test_validate_numeric_thresholds_valid():
    """Test valid numeric thresholds."""
    lfc, pval = validate_numeric_thresholds("1.5", "0.05")
    assert lfc == 1.5
    assert pval == 0.05

def test_validate_numeric_thresholds_invalid():
    """Test invalid numeric thresholds."""
    with pytest.raises(ValidationError):
        validate_numeric_thresholds("not_a_number", "0.05")

    with pytest.raises(ValidationError):
        validate_numeric_thresholds("-1.0", "0.05")

    with pytest.raises(ValidationError):
        validate_numeric_thresholds("1.0", "1.5")  # p-value > 1

def test_validate_gene_count_valid():
    """Test valid gene count."""
    count = validate_gene_count("50")
    assert count == 50

def test_validate_gene_count_invalid():
    """Test invalid gene count."""
    with pytest.raises(ValidationError):
        validate_gene_count("0")

    with pytest.raises(ValidationError):
        validate_gene_count("-10")

    with pytest.raises(ValidationError):
        validate_gene_count("not_a_number")

def test_validate_heatmap_params_frequent():
    """Test heatmap validation with frequent DEGs."""
    contrasts = [('GSE123', 'contrast1')]
    params = validate_heatmap_params(
        gene_selection_method="Frequent DEGs",
        gene_count_input="50",
        custom_genes_list=[],
        lfc_thresh="1.0",
        pvalue_thresh="0.05",
        selected_contrasts=contrasts,
        available_genes=set()
    )

    assert params['method'] == 'frequent'
    assert params['gene_count'] == 50

def test_validate_heatmap_params_custom():
    """Test heatmap validation with custom genes."""
    contrasts = [('GSE123', 'contrast1')]
    available = {'GENE1', 'GENE2', 'GENE3'}

    params = validate_heatmap_params(
        gene_selection_method="Custom",
        gene_count_input="",
        custom_genes_list=['GENE1', 'GENE2', 'MISSING'],
        lfc_thresh="1.0",
        pvalue_thresh="0.05",
        selected_contrasts=contrasts,
        available_genes=available
    )

    assert params['method'] == 'custom'
    assert params['genes'] == ['GENE1', 'GENE2']
    assert params['missing_genes'] == ['MISSING']
```

#### 3. Refactor Heatmap Tab to Use Business Logic

**File**: `main_workflow/reporting/streamlit_tabs/heatmap_tab.py`
**Changes**: Split _render_combined_heatmap_form to separate UI and validation

```python
# Add imports
from business_logic.validation import validate_heatmap_params, ValidationError
from business_logic.gene_selection import get_available_genes_for_contrasts
from business_logic.heatmap_data import extract_lfc_matrix, cluster_heatmap_data, simplify_contrast_labels

# Refactor _render_combined_heatmap_form (lines 377-707)
# Split into two functions:

def _render_heatmap_forms(
    ri: ResultsIntegrator,
    selected_datasets: List[str],
    contrast_data: List[Dict[str, Any]]
) -> Optional[Tuple[str, str, List[str], str, str, List[Tuple[str, str]]]]:
    """
    Render heatmap configuration forms (UI only).

    Returns:
        Tuple of form inputs or None if not submitted
    """
    # UI code for contrast selection form
    with st.form("heatmap_contrasts_form"):
        st.subheader("1. Select Contrasts")
        # ... existing UI code ...
        contrasts_submitted = st.form_submit_button("Confirm Contrast Selection")

    # UI code for gene selection form
    gene_selection_method = st.radio(...)
    with st.form("heatmap_gene_form"):
        # ... existing UI code ...
        submitted = st.form_submit_button("Generate Heatmap")

    if submitted:
        return (
            gene_selection_method,
            gene_count_input,
            custom_genes_list,
            lfc_thresh,
            pvalue_thresh,
            selected_contrasts
        )
    return None


def _validate_and_process_heatmap_form(
    ri: ResultsIntegrator,
    form_inputs: Tuple
) -> Optional[Dict[str, Any]]:
    """
    Validate heatmap form inputs and prepare parameters (pure logic).

    Returns:
        Validated parameters dict or None if validation fails
    """
    (gene_selection_method, gene_count_input, custom_genes_list,
     lfc_thresh, pvalue_thresh, selected_contrasts) = form_inputs

    # Get available genes (business logic)
    available_genes = get_available_genes_for_contrasts(
        ri.deg_data, selected_contrasts
    )

    try:
        # Validate using business logic
        params = validate_heatmap_params(
            gene_selection_method=gene_selection_method,
            gene_count_input=gene_count_input,
            custom_genes_list=custom_genes_list,
            lfc_thresh=lfc_thresh,
            pvalue_thresh=pvalue_thresh,
            selected_contrasts=selected_contrasts,
            available_genes=available_genes
        )

        # Display validation details
        with st.expander("Gene and contrast validation", expanded=False):
            st.info(f"✓ {len(selected_contrasts)} contrasts selected")
            if params['method'] == 'custom':
                st.success(f"✓ {len(params['genes'])} genes found")
                if params['missing_genes']:
                    st.warning(f"⚠ {len(params['missing_genes'])} genes not found")

        return params

    except ValidationError as e:
        st.error(str(e))
        return None
```

### Success Criteria

#### Automated Verification:
- [ ] All unit tests pass: `pytest tests/unit/ -v`
- [ ] Integration tests pass: `pytest tests/integration/ -v` (to be created)
- [ ] Test coverage ≥80%: `pytest --cov=main_workflow/reporting/business_logic --cov-report=html`
- [ ] Type checking passes: `pyright main_workflow/reporting/business_logic/`
- [ ] No regressions in existing tests
- [ ] Performance within 5% of original: Compare heatmap generation time

#### Manual Verification:
- [ ] Heatmap tab form validation shows same error messages
- [ ] Heatmap with frequent DEGs produces identical output
- [ ] Heatmap with custom genes produces identical output
- [ ] Clustering produces same gene/contrast ordering
- [ ] Multi-contrast heatmaps render correctly
- [ ] Missing gene warnings display correctly

---

## Phase 4: Consolidate and Optimize

### Overview
Consolidate remaining duplicate code, optimize performance, and establish comprehensive testing coverage including integration tests.

### Changes Required

#### 1. Create Integration Tests

**File**: `main_workflow/reporting/tests/integration/test_tab_workflows.py` (NEW)
**Changes**: Test complete workflows

```python
"""Integration tests for tab workflows."""
import pytest
import pandas as pd
from unittest.mock import Mock
from business_logic import (
    data_formatters,
    gene_selection,
    heatmap_data,
    validation,
    organism_utils
)

@pytest.fixture
def mock_results_integrator(sample_analysis_info, sample_deg_data_dict, sample_cpm_data):
    """Mock ResultsIntegrator for integration tests."""
    ri = Mock()
    ri.analysis_info = sample_analysis_info
    ri.deg_data = sample_deg_data_dict
    ri.cpm_data = {'GSE123456': sample_cpm_data}
    ri.dataset_info = {}
    ri.contrast_info = {}
    return ri

def test_heatmap_workflow_end_to_end(mock_results_integrator):
    """Test complete heatmap generation workflow."""
    # Step 1: Create contrast table
    contrast_data = data_formatters.create_contrast_table_data(
        analysis_info=mock_results_integrator.analysis_info,
        deg_data=mock_results_integrator.deg_data,
        contrast_info=mock_results_integrator.contrast_info,
        selected_datasets=['GSE123456']
    )
    assert len(contrast_data) > 0

    # Step 2: Identify frequent DEGs
    contrasts = [('GSE123456', 'Treatment_vs_Control')]
    genes = gene_selection.identify_frequent_degs(
        deg_data=mock_results_integrator.deg_data,
        contrasts=contrasts,
        top_n=3,
        p_thresh=0.05,
        lfc_thresh=1.0
    )
    assert len(genes) > 0

    # Step 3: Extract LFC matrix
    labels = heatmap_data.simplify_contrast_labels(
        contrasts, mock_results_integrator.analysis_info
    )
    matrix = heatmap_data.extract_lfc_matrix(
        deg_data=mock_results_integrator.deg_data,
        genes=genes,
        contrasts=contrasts,
        contrast_labels=labels,
        p_thresh=0.05,
        lfc_thresh=1.0
    )
    assert not matrix.empty
    assert 'Gene' in matrix.columns

    # Step 4: Cluster
    ordered_genes, ordered_contrasts = heatmap_data.cluster_heatmap_data(
        matrix, cluster_genes=True, cluster_contrasts=True
    )
    assert len(ordered_genes) == len(genes)

def test_expression_workflow_with_organisms(mock_results_integrator):
    """Test expression plot workflow with multi-organism handling."""
    # Step 1: Group by organism
    datasets = ['GSE123456', 'GSE789012']
    organism_groups = organism_utils.group_datasets_by_organism(
        mock_results_integrator.analysis_info, datasets
    )
    assert len(organism_groups) == 2

    # Step 2: Get genes and filter by organism
    all_genes = gene_selection.get_all_genes_from_cpm(mock_results_integrator.cpm_data)

    for organism, org_datasets in organism_groups.items():
        contrasts = [('GSE123456', 'Treatment_vs_Control')]
        filtered_genes = organism_utils.filter_genes_by_organism(
            deg_data=mock_results_integrator.deg_data,
            cpm_data=mock_results_integrator.cpm_data,
            analysis_info=mock_results_integrator.analysis_info,
            genes=all_genes,
            organism=organism,
            selected_contrasts=contrasts
        )
        # Should have filtered to organism-specific genes
        assert len(filtered_genes) <= len(all_genes)
```

#### 2. Create Performance Benchmarks

**File**: `main_workflow/reporting/tests/test_performance.py` (NEW)
**Changes**: Performance regression tests

```python
"""Performance benchmarking tests."""
import pytest
import time
from business_logic import gene_selection, heatmap_data

@pytest.mark.slow
def test_gene_identification_performance(sample_deg_data_dict, benchmark_timer, performance_threshold):
    """Benchmark gene identification performance."""
    contrasts = [('GSE123456', 'Treatment_vs_Control')] * 10  # Simulate many contrasts

    with benchmark_timer() as timing:
        genes = gene_selection.identify_frequent_degs(
            deg_data=sample_deg_data_dict,
            contrasts=contrasts,
            top_n=100,
            p_thresh=0.05,
            lfc_thresh=1.0
        )

    # Should complete in reasonable time (adjust threshold as needed)
    assert timing['elapsed'] < 1.0, f"Gene identification took {timing['elapsed']:.3f}s"

@pytest.mark.slow
def test_lfc_matrix_extraction_performance(sample_deg_data_dict, benchmark_timer):
    """Benchmark LFC matrix extraction."""
    genes = ['GENE1', 'GENE2', 'GENE3', 'GENE4', 'GENE5']
    contrasts = [('GSE123456', 'Treatment_vs_Control')] * 5
    labels = [f"C{i}" for i in range(5)]

    with benchmark_timer() as timing:
        matrix = heatmap_data.extract_lfc_matrix(
            deg_data=sample_deg_data_dict,
            genes=genes,
            contrasts=contrasts,
            contrast_labels=labels,
            p_thresh=0.05,
            lfc_thresh=1.0
        )

    assert timing['elapsed'] < 0.5, f"Matrix extraction took {timing['elapsed']:.3f}s"

@pytest.mark.slow
def test_clustering_performance(benchmark_timer):
    """Benchmark clustering performance."""
    import pandas as pd
    import numpy as np

    # Create large test matrix
    n_genes = 100
    n_contrasts = 20
    df = pd.DataFrame({
        'Gene': [f'GENE{i}' for i in range(n_genes)],
        **{f'C{j}': np.random.randn(n_genes) for j in range(n_contrasts)}
    })

    with benchmark_timer() as timing:
        genes, contrasts = heatmap_data.cluster_heatmap_data(
            df, cluster_genes=True, cluster_contrasts=True
        )

    assert timing['elapsed'] < 2.0, f"Clustering took {timing['elapsed']:.3f}s"
```

#### 3. Update ResultsIntegration to Use Business Logic

**File**: `main_workflow/reporting/ResultsIntegration.py`
**Changes**: Refactor create_lfc_heatmap to use extracted functions

```python
# Add imports at top
from business_logic.heatmap_data import (
    extract_lfc_matrix,
    cluster_heatmap_data,
    simplify_contrast_labels
)

# Refactor create_lfc_heatmap method (lines 508-963)
def create_lfc_heatmap(
    self,
    genes: List[str],
    contrasts: List[Tuple[str, str]],
    p_thresh: float = 0.05,
    lfc_thresh: float = 1.0,
    cluster_genes: bool = True,
    cluster_contrasts: bool = True,
    title: Optional[str] = None,
    width: int = 1200,
    height: int = 800
) -> go.Figure:
    """
    Create a heatmap of log fold-changes across contrasts.

    Now uses extracted business logic for data preparation.
    """
    # Simplify contrast labels using business logic
    contrast_labels = simplify_contrast_labels(contrasts, self.analysis_info)

    # Extract LFC matrix using business logic
    heatmap_df = extract_lfc_matrix(
        deg_data=self.deg_data,
        genes=genes,
        contrasts=contrasts,
        contrast_labels=contrast_labels,
        p_thresh=p_thresh,
        lfc_thresh=lfc_thresh
    )

    if heatmap_df.empty:
        # Return empty figure
        fig = go.Figure()
        fig.add_annotation(text="No data available", showarrow=False)
        return fig

    # Cluster using business logic
    if cluster_genes or cluster_contrasts:
        ordered_genes, ordered_contrasts = cluster_heatmap_data(
            heatmap_df,
            cluster_genes=cluster_genes,
            cluster_contrasts=cluster_contrasts
        )
        # Reorder dataframe
        heatmap_df = heatmap_df.set_index('Gene').loc[ordered_genes, ordered_contrasts].reset_index()

    # Create Plotly figure (keep existing visualization code)
    # Lines 819-963 remain largely unchanged...
    # [Existing Plotly figure creation code]
```

#### 4. Create Comparison Test for Equivalence

**File**: `main_workflow/reporting/tests/integration/test_equivalence.py` (NEW)
**Changes**: Verify refactored code produces identical results

```python
"""Equivalence tests comparing old vs new implementations."""
import pytest
import pandas as pd
import numpy as np

def test_lfc_matrix_equivalence():
    """
    Verify that extracted LFC matrix function produces identical results
    to original implementation.

    This test would need access to original implementation for comparison.
    For now, it validates output structure and values.
    """
    # This is a template - actual implementation would compare against
    # saved output from original code or run both versions side-by-side
    pass

def test_contrast_table_equivalence(sample_analysis_info, sample_deg_data_dict):
    """Verify contrast table creation matches original."""
    from business_logic.data_formatters import create_contrast_table_data

    result = create_contrast_table_data(
        analysis_info=sample_analysis_info,
        deg_data=sample_deg_data_dict,
        contrast_info={},
        selected_datasets=['GSE123456', 'GSE789012']
    )

    # Verify structure
    assert all('Accession' in row for row in result)
    assert all('Contrast' in row for row in result)
    assert all('Select' in row for row in result)

    # Verify sorting (GSE with lower number should come first)
    accessions = [row['Accession'] for row in result]
    assert accessions[0] == 'GSE123456'
    assert accessions[1] == 'GSE789012'
```

#### 5. Add Comprehensive Documentation

**File**: `main_workflow/reporting/business_logic/README.md` (NEW)
**Changes**: Document the business logic modules

```markdown
# UORCA Business Logic Modules

This directory contains pure business logic extracted from the UORCA Streamlit app.
All functions here are UI-agnostic and fully testable.

## Module Overview

### `organism_utils.py`
Organism grouping and filtering utilities for multi-species analyses.

**Key Functions:**
- `group_datasets_by_organism()` - Group datasets by species
- `filter_genes_by_organism()` - Filter genes to organism-specific sets
- `get_organism_display_name()` - Format organism names for display

### `data_formatters.py`
Data formatting and table creation for UI display.

**Key Functions:**
- `sort_by_geo_accession()` - Sort DataFrames by GEO accession number
- `create_contrast_table_data()` - Generate contrast selection tables
- `create_dataset_info_table()` - Generate dataset information tables
- `create_contrast_info_table()` - Generate contrast statistics tables

### `gene_selection.py`
Gene identification and selection logic.

**Key Functions:**
- `identify_frequent_degs()` - Find genes frequently DE across contrasts
- `get_all_genes_from_cpm()` - Extract all genes from expression data
- `validate_custom_genes()` - Check custom gene lists against available data

### `heatmap_data.py`
Heatmap data preparation and clustering.

**Key Functions:**
- `extract_lfc_matrix()` - Build log fold-change matrix with filtering
- `cluster_heatmap_data()` - Perform hierarchical clustering
- `simplify_contrast_labels()` - Create concise contrast labels

### `validation.py`
Form input validation logic.

**Key Functions:**
- `validate_numeric_thresholds()` - Parse and validate LFC/p-value thresholds
- `validate_gene_count()` - Validate gene count input
- `validate_heatmap_params()` - Comprehensive heatmap form validation

## Testing

All modules have corresponding test files in `tests/unit/`.

Run tests:
```bash
pytest tests/unit/test_organism_utils.py -v
pytest tests/unit/test_data_formatters.py -v
pytest tests/unit/test_gene_selection.py -v
pytest tests/unit/test_heatmap_data.py -v
pytest tests/unit/test_validation.py -v
```

## Performance

Performance benchmarks are in `tests/test_performance.py`.

Run benchmarks:
```bash
pytest tests/test_performance.py -v -m slow
```
```

### Success Criteria

#### Automated Verification:
- [ ] All tests pass: `pytest tests/ -v`
- [ ] Test coverage ≥80%: `pytest --cov=main_workflow/reporting/business_logic --cov-report=html`
- [ ] Integration tests pass: `pytest tests/integration/ -v`
- [ ] Performance tests pass: `pytest tests/test_performance.py -v`
- [ ] Type checking passes: `pyright main_workflow/reporting/`
- [ ] No code duplication: `pylint --disable=all --enable=duplicate-code main_workflow/reporting/`

#### Manual Verification:
- [ ] All tabs render correctly
- [ ] Heatmaps identical to pre-refactor output
- [ ] Expression plots identical to pre-refactor output
- [ ] Dataset/contrast tables show same data
- [ ] Error messages display correctly
- [ ] Performance is comparable (within 5%)
- [ ] Download functionality works
- [ ] Multi-organism workflows function correctly

---

## Testing Strategy

### Unit Tests
**Coverage**: ≥80% of business_logic/ modules

**Focus Areas:**
- Pure functions in all business_logic modules
- Edge cases (empty data, missing columns, invalid inputs)
- Error handling (ValidationError exceptions)
- Data transformations (sorting, filtering, grouping)

**Example Test Structure:**
```python
def test_function_name_scenario():
    """Test description."""
    # Arrange
    input_data = {...}

    # Act
    result = function(input_data)

    # Assert
    assert result == expected
```

### Integration Tests
**Coverage**: Key workflows across multiple modules

**Focus Areas:**
- Complete heatmap generation workflow
- Expression plot data preparation
- Multi-organism handling
- Form validation → data processing → output generation

### Performance Tests
**Coverage**: Critical performance-sensitive operations

**Focus Areas:**
- Gene identification with 100+ contrasts
- LFC matrix extraction with 1000+ genes
- Clustering with 100x20 matrices
- Table creation with 50+ datasets

**Acceptance Criteria**: New code ≤ 5% slower than original

### Manual Testing Steps

#### 1. Baseline Comparison
Before refactoring, capture baseline:
```bash
# Generate heatmap with known parameters
# Save output as baseline_heatmap.png
# Record generation time

# Generate expression plots
# Save output as baseline_expression.png
# Record generation time
```

#### 2. Post-Refactor Comparison
After each phase:
```bash
# Generate same heatmap with identical parameters
# Compare output to baseline (visual inspection)
# Compare generation time (should be within 5%)

# Generate same expression plots
# Compare output to baseline
# Compare generation time
```

#### 3. Edge Case Testing
- Empty gene list
- Single contrast selected
- All genes missing from data
- Invalid threshold values (negative, zero)
- Very large gene lists (1000+)
- Very large contrast lists (50+)
- Mixed organism selections

#### 4. UI Interaction Testing
- Select/deselect datasets
- Change form parameters
- Submit forms with errors
- Download PDF/ZIP files
- Navigate between tabs
- Refresh/reload page

---

## Performance Considerations

### Caching Strategy (Preserved)
- `@st.cache_resource` for ResultsIntegrator (expensive I/O)
- `@st.cache_data(ttl=3600)` for computations (gene identification, figures)
- Session state for UI selections

### Optimization Opportunities

**Gene Identification** (gene_selection.py):
- Current: O(n*m) where n=contrasts, m=genes per contrast
- Already optimized with dict-based counting

**LFC Matrix Extraction** (heatmap_data.py):
- Current: O(g*c) where g=genes, c=contrasts
- Potential optimization: Vectorized pandas operations instead of loops

**Clustering** (heatmap_data.py):
- Current: scipy.cluster.hierarchy (efficient)
- Already optimized for typical data sizes (100 genes x 20 contrasts)

### Performance Targets
- Gene identification: <1s for 100 contrasts
- LFC matrix extraction: <0.5s for 100 genes x 20 contrasts
- Clustering: <2s for 100x20 matrix
- Full heatmap generation: <5s total
- Expression plot generation: <3s for 10 genes x 5 datasets

---

## Migration Notes

### Deployment Steps

1. **Phase 1**: Can deploy immediately (only removals and infrastructure)
2. **Phase 2**: Test thoroughly with sample data before deployment
3. **Phase 3**: Requires validation against production data
4. **Phase 4**: Final optimization and monitoring

### Rollback Strategy

Each phase uses git branching:
```bash
git checkout -b phase-1-setup
# Implement Phase 1
git commit -m "Phase 1: Setup and quick wins"

git checkout -b phase-2-extract-logic
# Implement Phase 2
git commit -m "Phase 2: Extract core data logic"

# If issues arise, rollback:
git checkout main
git reset --hard <previous-commit>
```

### Data Compatibility

- **No changes** to input file formats (CSV, JSON)
- **No changes** to output file formats (PDF, ZIP)
- **No changes** to ResultsIntegrator public API (initially)
- Internal function signatures changed, but encapsulated

### Environment Setup

Before starting:
```bash
# Install testing dependencies
uv add --dev pytest pytest-cov pytest-mock pyright

# Verify environment
pytest --version
pyright --version

# Run tests (should fail initially, no tests exist yet)
pytest tests/
```

---

## References

- **Original Analysis**: Detailed findings document from research phase
- **Streamlit Docs**: https://docs.streamlit.io/
- **Pytest Docs**: https://docs.pytest.org/
- **Pyright Docs**: https://microsoft.github.io/pyright/
- **UORCA README**: `README.md` in project root