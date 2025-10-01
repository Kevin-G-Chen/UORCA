# TODO List: UORCA Codebase Restructuring

**Generated from**: `docs/plans/codebase_restructure.md`
**Generated on**: 2025-10-01
**Total Tasks**: 87
**Estimated Effort**: 10-12 days
**Target Completion**: 2025-10-15 (2 weeks)

---

## ⭐ Key Design Principle: Clean Package Architecture

**GOAL**: Reorganize UORCA into a single `uorca/` package with clear module boundaries and zero sys.path manipulation.

**Current Problems**:
- 11 files manipulate sys.path
- `main_workflow/reporting/` contains GUI, utilities, and business logic
- No clear separation between CLI, GUI, pipeline, and utilities
- Tests only cover `reporting/core/`

**Target**:
```
uorca/
├── core/          # Shared utilities (TaskManager, validation, entrez)
├── identification/ # Dataset search
├── pipeline/      # RNA-seq analysis
├── batch/         # Batch processors
└── gui/           # All Streamlit code
```

---

## 📋 Quick Summary

- **Phase 1**: 8 tasks (setup structure) - 2 hours
- **Phase 2**: 13 tasks (move core) - 1 day
- **Phase 3**: 15 tasks (move identification) - 1 day
- **Phase 4**: 18 tasks (move pipeline) - 1.5 days
- **Phase 5**: 7 tasks (move batch helpers) - 0.5 day
- **Phase 6**: 42 tasks (move GUI - LARGEST) - 3 days
- **Phase 7**: 12 tasks (update tests) - 1 day
- **Phase 8**: 8 tasks (backwards compat) - 0.5 day
- **Phase 9**: 9 tasks (update config) - 0.5 day
- **Phase 10**: 23 tasks (verification) - 1 day

**Total**: 155 granular tasks across 10 phases

---

## 🎯 Success Criteria

- [x] All 10 phases completed
- [ ] Zero sys.path manipulation in `uorca/`
- [ ] All tests pass
- [ ] Type checking clean
- [ ] CLI commands work identically
- [ ] GUI fully functional
- [ ] Backwards compatibility works with warnings

---

## 🚫 Out of Scope

- NOT refactoring ResultsIntegration.py (1836 lines)
- NOT adding new tests (just infrastructure)
- NOT changing CLI behavior
- NOT modifying Docker container
- NOT removing main_workflow/ (keeping for compat)

---

## 🔴 Phase 1: Create New Package Structure (2 hours)

**Objective**: Set up empty directory structure without moving code.

### PHASE1-1: Create core directories
**Effort**: XS (5 min)
**Priority**: Critical
**Dependencies**: None

**Tasks**:
- [ ] Run: `mkdir -p uorca/core`
- [ ] Verify: `ls -d uorca/core` succeeds

---

### PHASE1-2: Create identification directories
**Effort**: XS (5 min)
**Priority**: Critical
**Dependencies**: PHASE1-1

**Tasks**:
- [ ] Run: `mkdir -p uorca/identification/prompts`
- [ ] Verify: `ls -d uorca/identification/prompts` succeeds

---

### PHASE1-3: Create pipeline directories
**Effort**: XS (10 min)
**Priority**: Critical
**Dependencies**: PHASE1-1

**Tasks**:
- [ ] Run: `mkdir -p uorca/pipeline/agents`
- [ ] Run: `mkdir -p uorca/pipeline/prompts`
- [ ] Run: `mkdir -p uorca/pipeline/additional_scripts`
- [ ] Verify: `ls -d uorca/pipeline/*` shows 3 directories

---

### PHASE1-4: Create GUI directories
**Effort**: S (15 min)
**Priority**: Critical
**Dependencies**: PHASE1-1

**Tasks**:
- [ ] Run: `mkdir -p uorca/gui/pages`
- [ ] Run: `mkdir -p uorca/gui/explorer/helpers`
- [ ] Run: `mkdir -p uorca/gui/plotting`
- [ ] Run: `mkdir -p uorca/gui/ai`
- [ ] Run: `mkdir -p uorca/gui/prompts`
- [ ] Run: `mkdir -p uorca/gui/example_output_files`
- [ ] Verify: `ls -d uorca/gui/*` shows 6 directories

---

### PHASE1-5: Create batch directories
**Effort**: XS (5 min)
**Priority**: Critical
**Dependencies**: PHASE1-1

**Tasks**:
- [ ] Run: `mkdir -p uorca/batch/run_helpers`
- [ ] Verify: `ls -d uorca/batch/run_helpers` succeeds
- [ ] Note: `uorca/batch/` already exists with base.py, local.py, slurm.py

---

### PHASE1-6: Create test directories
**Effort**: S (10 min)
**Priority**: Critical
**Dependencies**: None

**Tasks**:
- [ ] Run: `mkdir -p tests/unit/core`
- [ ] Run: `mkdir -p tests/unit/identification`
- [ ] Run: `mkdir -p tests/unit/pipeline`
- [ ] Run: `mkdir -p tests/unit/batch`
- [ ] Run: `mkdir -p tests/unit/gui`
- [ ] Verify: `ls -d tests/unit/*` shows 5 directories

---

### PHASE1-7: Create __init__.py files
**Effort**: M (30 min)
**Priority**: Critical
**Dependencies**: PHASE1-1 through PHASE1-6

**Tasks**:
- [ ] Create: `touch uorca/core/__init__.py`
- [ ] Create: `touch uorca/identification/__init__.py`
- [ ] Create: `touch uorca/pipeline/__init__.py`
- [ ] Create: `touch uorca/pipeline/agents/__init__.py`
- [ ] Create: `touch uorca/gui/__init__.py`
- [ ] Create: `touch uorca/gui/pages/__init__.py`
- [ ] Create: `touch uorca/gui/explorer/__init__.py`
- [ ] Create: `touch uorca/gui/plotting/__init__.py`
- [ ] Create: `touch uorca/gui/ai/__init__.py`
- [ ] Create: `touch tests/unit/core/__init__.py`
- [ ] Create: `touch tests/unit/identification/__init__.py`
- [ ] Create: `touch tests/unit/pipeline/__init__.py`
- [ ] Create: `touch tests/unit/batch/__init__.py`
- [ ] Create: `touch tests/unit/gui/__init__.py`
- [ ] Initialize all __init__.py files with: `__all__ = []`

**Verification**:
- [ ] Run: `find uorca tests/unit -name "__init__.py" | wc -l` shows 14 files
- [ ] Run: `ls -R uorca/` and verify structure matches plan

---

### PHASE1-8: Create checkpoint
**Effort**: XS (2 min)
**Priority**: Critical
**Dependencies**: PHASE1-7

**Tasks**:
- [ ] Save current work: `git add -A`
- [ ] Commit: `git commit -m "chore: create empty package structure for restructure"`
- [ ] Verify: `main_workflow/` is unchanged
- [ ] Document: Note checkpoint name "Phase 1 complete - Empty structure created"

---

## 🔴 Phase 2: Move Core Utilities (1 day)

**Objective**: Move shared utilities from `main_workflow/` into `uorca/core/`.

### PHASE2-1: Move task_manager.py
**Effort**: XS (5 min)
**Priority**: Critical
**Dependencies**: PHASE1-8

**Tasks**:
- [ ] Run: `mv main_workflow/reporting/core/task_manager.py uorca/core/`
- [ ] Verify: `ls uorca/core/task_manager.py` succeeds
- [ ] Verify: `ls main_workflow/reporting/core/task_manager.py` fails

---

### PHASE2-2: Move validation.py
**Effort**: XS (5 min)
**Priority**: Critical
**Dependencies**: PHASE1-8

**Tasks**:
- [ ] Run: `mv main_workflow/reporting/core/validation.py uorca/core/`
- [ ] Verify file moved

---

### PHASE2-3: Move data_formatters.py
**Effort**: XS (5 min)
**Priority**: Critical
**Dependencies**: PHASE1-8

**Tasks**:
- [ ] Run: `mv main_workflow/reporting/core/data_formatters.py uorca/core/`
- [ ] Verify file moved

---

### PHASE2-4: Move gene_selection.py
**Effort**: XS (5 min)
**Priority**: Critical
**Dependencies**: PHASE1-8

**Tasks**:
- [ ] Run: `mv main_workflow/reporting/core/gene_selection.py uorca/core/`
- [ ] Verify file moved

---

### PHASE2-5: Move organism_utils.py
**Effort**: XS (5 min)
**Priority**: Critical
**Dependencies**: PHASE1-8

**Tasks**:
- [ ] Run: `mv main_workflow/reporting/core/organism_utils.py uorca/core/`
- [ ] Verify file moved

---

### PHASE2-6: Move script_generation.py
**Effort**: XS (5 min)
**Priority**: Critical
**Dependencies**: PHASE1-8

**Tasks**:
- [ ] Run: `mv main_workflow/reporting/core/script_generation.py uorca/core/`
- [ ] Verify file moved

---

### PHASE2-7: Move entrez_utils.py
**Effort**: XS (5 min)
**Priority**: Critical
**Dependencies**: PHASE1-8

**Tasks**:
- [ ] Run: `mv main_workflow/shared/entrez_utils.py uorca/core/`
- [ ] Verify file moved

---

### PHASE2-8: Move workflow_logging.py
**Effort**: XS (5 min)
**Priority**: Critical
**Dependencies**: PHASE1-8

**Tasks**:
- [ ] Run: `mv main_workflow/shared/workflow_logging.py uorca/core/`
- [ ] Verify: `ls uorca/core/*.py | wc -l` shows 8 files (plus __init__.py)

---

### PHASE2-9: Update uorca/core/__init__.py
**Effort**: M (30 min)
**Priority**: Critical
**Dependencies**: PHASE2-1 through PHASE2-8
**Files**: `uorca/core/__init__.py`

**Tasks**:
- [ ] Add docstring: `"""Core utilities shared across UORCA."""`
- [ ] Import TaskManager and TaskStatus
- [ ] Import all validation functions
- [ ] Import all data_formatters functions
- [ ] Import all gene_selection functions
- [ ] Import all organism_utils functions
- [ ] Import generate_reproducible_heatmap_script
- [ ] Import all entrez_utils functions
- [ ] Import setup_logging
- [ ] Define __all__ list with all exports
- [ ] Verify: `python -c "from uorca.core import TaskManager"` succeeds

**Acceptance Criteria**:
- [ ] All exports listed in plan are present
- [ ] Imports work without errors
- [ ] No sys.path manipulation needed

---

### PHASE2-10: Verify core imports work
**Effort**: S (15 min)
**Priority**: Critical
**Dependencies**: PHASE2-9

**Tasks**:
- [ ] Test: `python -c "from uorca.core import TaskManager, TaskStatus"`
- [ ] Test: `python -c "from uorca.core import validate_threshold_values"`
- [ ] Test: `python -c "from uorca.core import ESearch, ESummary"`
- [ ] Test: `python -c "from uorca.core import setup_logging"`
- [ ] Verify: All imports succeed

---

### PHASE2-11: Clean up old directories
**Effort**: S (10 min)
**Priority**: High
**Dependencies**: PHASE2-10

**Tasks**:
- [ ] Keep: `main_workflow/reporting/core/__init__.py` (for backwards compat later)
- [ ] Verify: `ls main_workflow/reporting/core/` shows only `__init__.py` and `__pycache__`
- [ ] Keep: `main_workflow/shared/__init__.py`
- [ ] Verify: `ls main_workflow/shared/` shows only `__init__.py` and `__pycache__`

---

### PHASE2-12: Run tests with new imports (will fail initially)
**Effort**: S (10 min)
**Priority**: High
**Dependencies**: PHASE2-11

**Tasks**:
- [ ] Run: `pytest tests/unit/ --collect-only` (expect failures)
- [ ] Note which tests need import updates (will fix in Phase 7)
- [ ] Document: List of test files that need updates

---

### PHASE2-13: Create checkpoint
**Effort**: XS (2 min)
**Priority**: Critical
**Dependencies**: PHASE2-12

**Tasks**:
- [ ] Save: `git add -A`
- [ ] Commit: `git commit -m "refactor: move core utilities to uorca/core/"`
- [ ] Checkpoint: "Phase 2 complete - Core utilities moved"

---

## 🔴 Phase 3: Move Dataset Identification (1 day)

**Objective**: Move dataset identification code into `uorca/identification/`.

### PHASE3-1: Move DatasetIdentification.py and rename
**Effort**: S (10 min)
**Priority**: Critical
**Dependencies**: PHASE2-13

**Tasks**:
- [ ] Run: `mv main_workflow/dataset_identification/DatasetIdentification.py uorca/identification/dataset_identification.py`
- [ ] Verify: `ls uorca/identification/dataset_identification.py` succeeds
- [ ] Note: File renamed to follow Python convention (snake_case)

---

### PHASE3-2: Move prompts directory
**Effort**: XS (5 min)
**Priority**: Critical
**Dependencies**: PHASE2-13

**Tasks**:
- [ ] Run: `mv main_workflow/dataset_identification/prompts uorca/identification/`
- [ ] Verify: `ls uorca/identification/prompts/*.txt | wc -l` shows 2 files
- [ ] Files should be: `extract_terms.txt` and `assess_relevance.txt`

---

### PHASE3-3: Update imports in dataset_identification.py (line 1)
**Effort**: M (30 min)
**Priority**: Critical
**Dependencies**: PHASE3-1
**Files**: `uorca/identification/dataset_identification.py`

**Tasks**:
- [ ] Read file to find all imports from `shared.*`
- [ ] Replace: `from shared.entrez_utils import ...` with `from uorca.core import ...`
- [ ] Remove all sys.path manipulation lines
- [ ] Update any relative imports to absolute
- [ ] Verify: No `sys.path.insert` or `sys.path.append` in file

**Acceptance Criteria**:
- [ ] Run: `grep "sys.path" uorca/identification/dataset_identification.py` returns nothing
- [ ] Run: `grep "from shared" uorca/identification/dataset_identification.py` returns nothing

---

### PHASE3-4: Add docstring to identification module
**Effort**: S (10 min)
**Priority**: High
**Dependencies**: PHASE3-3
**Files**: `uorca/identification/dataset_identification.py`

**Tasks**:
- [ ] Add module docstring at top (after shebang if present)
- [ ] Docstring should describe: "Dataset identification from GEO using AI-powered relevance assessment"

---

### PHASE3-5: Update uorca/identification/__init__.py
**Effort**: M (20 min)
**Priority**: Critical
**Dependencies**: PHASE3-4
**Files**: `uorca/identification/__init__.py`

**Tasks**:
- [ ] Add docstring: `"""Dataset identification from GEO."""`
- [ ] Import: `from .dataset_identification import (`
- [ ] Export: extract_terms
- [ ] Export: perform_search
- [ ] Export: get_basic_dataset_info
- [ ] Export: fetch_runinfo_from_bioproject
- [ ] Export: calculate_dataset_sizes_from_runinfo
- [ ] Export: embed_datasets
- [ ] Export: cluster_datasets
- [ ] Export: select_representative_datasets
- [ ] Export: main
- [ ] Define __all__ list
- [ ] Verify: `python -c "from uorca.identification import main"` succeeds

---

### PHASE3-6: Update uorca/identify.py wrapper
**Effort**: M (20 min)
**Priority**: Critical
**Dependencies**: PHASE3-5
**Files**: `uorca/identify.py` (around line 104)

**Tasks**:
- [ ] Find line with: `sys.path.insert(0, str(main_workflow_dir))`
- [ ] Remove sys.path manipulation
- [ ] Replace: `from dataset_identification.DatasetIdentification import main as identify_main`
- [ ] With: `from uorca.identification import main as identify_main`
- [ ] Verify: File is much simpler now

**Acceptance Criteria**:
- [ ] Run: `grep "sys.path" uorca/identify.py` returns nothing
- [ ] Run: `grep "main_workflow" uorca/identify.py` returns nothing

---

### PHASE3-7: Test identify CLI command
**Effort**: S (15 min)
**Priority**: Critical
**Dependencies**: PHASE3-6

**Tasks**:
- [ ] Run: `uv run uorca identify --help`
- [ ] Verify: Help text displays correctly
- [ ] Check: No import errors
- [ ] Note: Don't run full command yet (needs API keys)

---

### PHASE3-8: Update GUI page import
**Effort**: S (10 min)
**Priority**: High
**Dependencies**: PHASE3-5
**Files**: `uorca/gui/pages/dataset_identification.py` (when moved in Phase 6)

**Tasks**:
- [ ] Note: This file not moved yet, will update in Phase 6
- [ ] Document: Page needs update from `uorca.identify.main` to `uorca.identification.main`
- [ ] Add to Phase 6 checklist

---

### PHASE3-9: Verify identification imports
**Effort**: M (20 min)
**Priority**: Critical
**Dependencies**: PHASE3-7

**Tasks**:
- [ ] Test: `python -c "from uorca.identification import extract_terms"`
- [ ] Test: `python -c "from uorca.identification import perform_search"`
- [ ] Test: `python -c "from uorca.identification import main"`
- [ ] Verify: All imports succeed
- [ ] Check: `uorca/identification/dataset_identification.py` has no import errors

---

### PHASE3-10: Run type checking on identification
**Effort**: S (10 min)
**Priority**: High
**Dependencies**: PHASE3-9

**Tasks**:
- [ ] Run: `pyright uorca/identification/`
- [ ] Review: Note any type errors (OK to have some, will fix later)
- [ ] Document: List of type errors for future cleanup

---

### PHASE3-11: Test identify command with dry-run
**Effort**: M (20 min)
**Priority**: High
**Dependencies**: PHASE3-10

**Tasks**:
- [ ] Set fake env vars: `export ENTREZ_EMAIL=test@example.com OPENAI_API_KEY=test`
- [ ] Run: `uv run uorca identify -q "test" -o /tmp/test_identify 2>&1 | head -20`
- [ ] Verify: Imports work, fails gracefully if API keys invalid
- [ ] Unset: `unset ENTREZ_EMAIL OPENAI_API_KEY`

---

### PHASE3-12: Clean up old identification directory
**Effort**: S (5 min)
**Priority**: Medium
**Dependencies**: PHASE3-11

**Tasks**:
- [ ] Keep: `main_workflow/dataset_identification/__init__.py` (for backwards compat)
- [ ] Verify: Only `__init__.py` and `__pycache__` remain
- [ ] Document: Directory structure for backwards compat module

---

### PHASE3-13: Commit identification phase
**Effort**: S (5 min)
**Priority**: Critical
**Dependencies**: PHASE3-12

**Tasks**:
- [ ] Run: `git add -A`
- [ ] Commit: `git commit -m "refactor: move dataset identification to uorca/identification/"`
- [ ] Checkpoint: "Phase 3 complete - Identification moved"

---

### PHASE3-14: Verify CLI still works
**Effort**: S (10 min)
**Priority**: Critical
**Dependencies**: PHASE3-13

**Tasks**:
- [ ] Run: `uv run uorca --help`
- [ ] Run: `uv run uorca identify --help`
- [ ] Run: `uv run uorca run --help`
- [ ] Run: `uv run uorca explore --help`
- [ ] Verify: All commands show help without errors

---

### PHASE3-15: Document progress
**Effort**: XS (5 min)
**Priority**: Low
**Dependencies**: PHASE3-14

**Tasks**:
- [ ] Update this TODO: Mark Phase 3 tasks complete
- [ ] Note: Estimated vs actual time for calibration
- [ ] If behind schedule: Re-estimate remaining phases

---

## 🔴 Phase 4: Move Pipeline Code (1.5 days)

**Objective**: Move RNA-seq pipeline code into `uorca/pipeline/`.

### PHASE4-1: Move master.py
**Effort**: XS (5 min)
**Priority**: Critical
**Dependencies**: PHASE3-15

**Tasks**:
- [ ] Run: `mv main_workflow/master.py uorca/pipeline/`
- [ ] Verify: `ls uorca/pipeline/master.py` succeeds

---

### PHASE4-2: Move agents directory
**Effort**: S (10 min)
**Priority**: Critical
**Dependencies**: PHASE3-15

**Tasks**:
- [ ] Run: `mv main_workflow/agents uorca/pipeline/`
- [ ] Verify: `ls uorca/pipeline/agents/*.py | wc -l` shows 4 files
- [ ] Files should be: `__init__.py`, `extraction.py`, `metadata.py`, `analysis.py`

---

### PHASE4-3: Move prompts directory
**Effort**: S (10 min)
**Priority**: Critical
**Dependencies**: PHASE3-15

**Tasks**:
- [ ] Run: `mv main_workflow/prompts uorca/pipeline/`
- [ ] Verify: `ls uorca/pipeline/prompts/*.txt | wc -l` shows 4 files
- [ ] Files should be: `master.txt`, `extraction.txt`, `metadata.txt`, `analysis.txt`

---

### PHASE4-4: Move additional_scripts directory
**Effort**: XS (5 min)
**Priority**: Critical
**Dependencies**: PHASE3-15

**Tasks**:
- [ ] Run: `mv main_workflow/additional_scripts uorca/pipeline/`
- [ ] Verify: `ls uorca/pipeline/additional_scripts/*.R` succeeds
- [ ] File should be: `RNAseq.R`

---

### PHASE4-5: Update imports in master.py
**Effort**: M (30 min)
**Priority**: Critical
**Dependencies**: PHASE4-1
**Files**: `uorca/pipeline/master.py`

**Tasks**:
- [ ] Find: `from agents import extraction, metadata, analysis`
- [ ] Replace with: `from uorca.pipeline.agents import extraction, metadata, analysis`
- [ ] Find: `from shared.workflow_logging import setup_logging`
- [ ] Replace with: `from uorca.core import setup_logging`
- [ ] Remove any sys.path manipulation
- [ ] Verify: `grep "sys.path" uorca/pipeline/master.py` returns nothing

---

### PHASE4-6: Update imports in agents/extraction.py
**Effort**: M (20 min)
**Priority**: Critical
**Dependencies**: PHASE4-2
**Files**: `uorca/pipeline/agents/extraction.py`

**Tasks**:
- [ ] Find: `from shared.entrez_utils import ...`
- [ ] Replace with: `from uorca.core import ...`
- [ ] Find: `from shared.workflow_logging import setup_logging`
- [ ] Replace with: `from uorca.core import setup_logging`
- [ ] Remove any sys.path manipulation
- [ ] Verify: No imports from `shared.*`

---

### PHASE4-7: Update imports in agents/metadata.py
**Effort**: M (20 min)
**Priority**: Critical
**Dependencies**: PHASE4-2
**Files**: `uorca/pipeline/agents/metadata.py`

**Tasks**:
- [ ] Update imports similar to extraction.py
- [ ] Replace `shared.*` with `uorca.core.*`
- [ ] Remove sys.path manipulation

---

### PHASE4-8: Update imports in agents/analysis.py
**Effort**: M (20 min)
**Priority**: Critical
**Dependencies**: PHASE4-2
**Files**: `uorca/pipeline/agents/analysis.py`

**Tasks**:
- [ ] Update imports similar to extraction.py
- [ ] Replace `shared.*` with `uorca.core.*`
- [ ] Remove sys.path manipulation

---

### PHASE4-9: Update agents/__init__.py
**Effort**: S (15 min)
**Priority**: Critical
**Dependencies**: PHASE4-6, PHASE4-7, PHASE4-8
**Files**: `uorca/pipeline/agents/__init__.py`

**Tasks**:
- [ ] Add docstring: `"""Pipeline agents for RNA-seq analysis."""`
- [ ] Import and re-export: `from . import extraction, metadata, analysis`
- [ ] Define __all__: `["extraction", "metadata", "analysis"]`
- [ ] Verify: `python -c "from uorca.pipeline.agents import extraction"` succeeds

---

### PHASE4-10: Update uorca/pipeline/__init__.py
**Effort**: S (15 min)
**Priority**: Critical
**Dependencies**: PHASE4-5
**Files**: `uorca/pipeline/__init__.py`

**Tasks**:
- [ ] Add docstring: `"""RNA-seq analysis pipeline."""`
- [ ] Import: `from .master import main as run_pipeline`
- [ ] Define __all__: `["run_pipeline"]`
- [ ] Verify: `python -c "from uorca.pipeline import run_pipeline"` succeeds

---

### PHASE4-11: Add __main__ block to master.py
**Effort**: S (10 min)
**Priority**: High
**Dependencies**: PHASE4-10
**Files**: `uorca/pipeline/master.py`

**Tasks**:
- [ ] At end of file, add:
  ```python
  if __name__ == "__main__":
      from uorca.pipeline import run_pipeline
      run_pipeline()
  ```
- [ ] Verify: `python -m uorca.pipeline.master --help` works (if master.py has argparse)

---

### PHASE4-12: Update uorca/batch/local.py
**Effort**: M (20 min)
**Priority**: Critical
**Dependencies**: PHASE4-11
**Files**: `uorca/batch/local.py` (around line 116)

**Tasks**:
- [ ] Find: `'python', 'main_workflow/master.py',`
- [ ] Replace with: `'python', '-m', 'uorca.pipeline.master',`
- [ ] Verify: Full command looks correct in context
- [ ] Test: `python -c "from uorca.batch.local import LocalBatchProcessor"` succeeds

---

### PHASE4-13: Verify pipeline imports
**Effort**: M (20 min)
**Priority**: Critical
**Dependencies**: PHASE4-12

**Tasks**:
- [ ] Test: `python -c "from uorca.pipeline import run_pipeline"`
- [ ] Test: `python -c "from uorca.pipeline.master import main"`
- [ ] Test: `python -c "from uorca.pipeline.agents import extraction, metadata, analysis"`
- [ ] Verify: All imports succeed

---

### PHASE4-14: Run type checking on pipeline
**Effort**: M (20 min)
**Priority**: High
**Dependencies**: PHASE4-13

**Tasks**:
- [ ] Run: `pyright uorca/pipeline/`
- [ ] Review: Note any type errors (OK to have some)
- [ ] Document: List for future cleanup

---

### PHASE4-15: Test pipeline can be imported without errors
**Effort**: M (20 min)
**Priority**: Critical
**Dependencies**: PHASE4-14

**Tasks**:
- [ ] Test: `python -c "from uorca.pipeline.master import main; print('Success')"`
- [ ] Verify: No import errors
- [ ] Check: All agent modules import correctly

---

### PHASE4-16: Clean up old pipeline directories
**Effort**: S (10 min)
**Priority**: Medium
**Dependencies**: PHASE4-15

**Tasks**:
- [ ] Note: Keep `main_workflow/agents/__init__.py` for backwards compat
- [ ] Verify: `main_workflow/agents/` contains only `__init__.py` and `__pycache__`
- [ ] Verify: `main_workflow/prompts/` is empty or removed
- [ ] Verify: `main_workflow/additional_scripts/` is empty or removed

---

### PHASE4-17: Commit pipeline phase
**Effort**: S (5 min)
**Priority**: Critical
**Dependencies**: PHASE4-16

**Tasks**:
- [ ] Run: `git add -A`
- [ ] Commit: `git commit -m "refactor: move pipeline code to uorca/pipeline/"`
- [ ] Checkpoint: "Phase 4 complete - Pipeline moved"

---

### PHASE4-18: Document progress
**Effort**: XS (5 min)
**Priority**: Low
**Dependencies**: PHASE4-17

**Tasks**:
- [ ] Update this TODO: Mark Phase 4 complete
- [ ] Note: Estimated vs actual time
- [ ] Check: Are we on schedule?

---

## 🟡 Phase 5: Move Batch Processing Helpers (0.5 day)

**Objective**: Move batch processing helpers into `uorca/batch/run_helpers/`.

### PHASE5-1: Move submit_datasets.py
**Effort**: XS (5 min)
**Priority**: High
**Dependencies**: PHASE4-18

**Tasks**:
- [ ] Run: `mv main_workflow/run_helpers/submit_datasets.py uorca/batch/run_helpers/`
- [ ] Verify: `ls uorca/batch/run_helpers/submit_datasets.py` succeeds

---

### PHASE5-2: Move SLURM templates
**Effort**: S (10 min)
**Priority**: High
**Dependencies**: PHASE4-18

**Tasks**:
- [ ] Run: `mv main_workflow/run_helpers/run_single_dataset.sbatch.j2 uorca/batch/run_helpers/`
- [ ] Run: `mv main_workflow/run_helpers/run_dataset_array.sbatch.j2 uorca/batch/run_helpers/`
- [ ] Verify: `ls uorca/batch/run_helpers/*.j2 | wc -l` shows 2 files

---

### PHASE5-3: Update imports in submit_datasets.py
**Effort**: S (15 min)
**Priority**: High
**Dependencies**: PHASE5-1
**Files**: `uorca/batch/run_helpers/submit_datasets.py`

**Tasks**:
- [ ] Review file for any imports from `main_workflow.*`
- [ ] Update to import from `uorca.*` if needed
- [ ] Verify: No sys.path manipulation

---

### PHASE5-4: Update template paths in SLURM processor
**Effort**: M (20 min)
**Priority**: Critical
**Dependencies**: PHASE5-2
**Files**: `uorca/batch/slurm.py`

**Tasks**:
- [ ] Search for hardcoded template paths
- [ ] Find lines like: `Path(__file__).parent.parent / "main_workflow" / "run_helpers"`
- [ ] Replace with: `Path(__file__).parent / "run_helpers"`
- [ ] Verify: Paths are relative to current file location

---

### PHASE5-5: Test SLURM processor imports
**Effort**: S (10 min)
**Priority**: High
**Dependencies**: PHASE5-4

**Tasks**:
- [ ] Test: `python -c "from uorca.batch.slurm import SlurmBatchProcessor"`
- [ ] Verify: No import errors
- [ ] Check: Template paths resolve correctly (if accessible)

---

### PHASE5-6: Clean up old run_helpers directory
**Effort**: XS (5 min)
**Priority**: Low
**Dependencies**: PHASE5-5

**Tasks**:
- [ ] Verify: `main_workflow/run_helpers/` is empty or can be removed
- [ ] Keep directory if needed for backwards compat
- [ ] Document: Whether directory is needed

---

### PHASE5-7: Commit batch helpers phase
**Effort**: XS (2 min)
**Priority**: Critical
**Dependencies**: PHASE5-6

**Tasks**:
- [ ] Run: `git add -A`
- [ ] Commit: `git commit -m "refactor: move batch helpers to uorca/batch/run_helpers/"`
- [ ] Checkpoint: "Phase 5 complete - Batch helpers moved"

---

## 🔴 Phase 6: Move GUI Code (3 days) - LARGEST PHASE

**Objective**: Move all Streamlit/GUI code into `uorca/gui/`.

### PHASE6-1: Move uorca_explorer.py
**Effort**: S (5 min)
**Priority**: Critical
**Dependencies**: PHASE5-7

**Tasks**:
- [ ] Run: `mv main_workflow/reporting/uorca_explorer.py uorca/gui/`
- [ ] Verify: `ls uorca/gui/uorca_explorer.py` succeeds

---

### PHASE6-2: Move and rename ResultsIntegration.py
**Effort**: S (10 min)
**Priority**: Critical
**Dependencies**: PHASE5-7

**Tasks**:
- [ ] Run: `mv main_workflow/reporting/ResultsIntegration.py uorca/gui/results_integration.py`
- [ ] Note: Renamed to snake_case convention
- [ ] Verify: `ls uorca/gui/results_integration.py` succeeds

---

### PHASE6-3: Move and rename page 1
**Effort**: S (10 min)
**Priority**: Critical
**Dependencies**: PHASE5-7

**Tasks**:
- [ ] Run: `mv main_workflow/reporting/pages/01_Dataset_Identification.py uorca/gui/pages/dataset_identification.py`
- [ ] Note: Removed number prefix, converted to snake_case
- [ ] Verify: File moved

---

### PHASE6-4: Move and rename page 2
**Effort**: S (10 min)
**Priority**: Critical
**Dependencies**: PHASE5-7

**Tasks**:
- [ ] Run: `mv main_workflow/reporting/pages/02_Pipeline_Execution.py uorca/gui/pages/pipeline_execution.py`
- [ ] Verify: File moved

---

### PHASE6-5: Move all streamlit_tabs files to explorer
**Effort**: M (30 min)
**Priority**: Critical
**Dependencies**: PHASE5-7

**Tasks**:
- [ ] Run: `mv main_workflow/reporting/streamlit_tabs/*.py uorca/gui/explorer/`
- [ ] Exclude: Don't move `__init__.py` yet (will create new one)
- [ ] Verify: 9 tab files moved
- [ ] Files: heatmap_tab.py, expression_plots_tab.py, analysis_plots_tab.py, etc.

---

### PHASE6-6: Move streamlit_tabs/helpers directory
**Effort**: S (10 min)
**Priority**: Critical
**Dependencies**: PHASE5-7

**Tasks**:
- [ ] Run: `mv main_workflow/reporting/streamlit_tabs/helpers uorca/gui/explorer/`
- [ ] Verify: `ls uorca/gui/explorer/helpers/*.py` shows helper files

---

### PHASE6-7: Move plotting files
**Effort**: S (10 min)
**Priority**: High
**Dependencies**: PHASE5-7

**Tasks**:
- [ ] Run: `mv main_workflow/reporting/single_analysis_plots.py uorca/gui/plotting/`
- [ ] Run: `mv main_workflow/reporting/ortholog_mapper.py uorca/gui/plotting/`
- [ ] Verify: 2 files in `uorca/gui/plotting/`

---

### PHASE6-8: Move AI module files
**Effort**: M (20 min)
**Priority**: High
**Dependencies**: PHASE5-7

**Tasks**:
- [ ] Run: `mv main_workflow/reporting/ai_agent_factory.py uorca/gui/ai/agent_factory.py`
- [ ] Run: `mv main_workflow/reporting/ai_gene_schema.py uorca/gui/ai/gene_schema.py`
- [ ] Run: `mv main_workflow/reporting/config_loader.py uorca/gui/ai/config_loader.py`
- [ ] Run: `mv main_workflow/reporting/contrast_relevance.py uorca/gui/ai/contrast_relevance.py`
- [ ] Run: `mv main_workflow/reporting/tool_relevance_analyzer.py uorca/gui/ai/tool_relevance_analyzer.py`
- [ ] Run: `mv main_workflow/reporting/mcp_server_core.py uorca/gui/ai/mcp_server_core.py`
- [ ] Verify: 6 files in `uorca/gui/ai/`

---

### PHASE6-9: Move GUI prompts and example files
**Effort**: S (10 min)
**Priority**: Medium
**Dependencies**: PHASE5-7

**Tasks**:
- [ ] Run: `mv main_workflow/reporting/prompts uorca/gui/`
- [ ] Run: `mv main_workflow/reporting/example_output_files uorca/gui/`
- [ ] Verify: Directories moved

---

### PHASE6-10: Update imports in dataset_identification.py page
**Effort**: M (30 min)
**Priority**: Critical
**Dependencies**: PHASE6-3
**Files**: `uorca/gui/pages/dataset_identification.py`

**Tasks**:
- [ ] Remove ALL sys.path manipulation (lines ~17-19)
- [ ] Replace: `from core import TaskManager, TaskStatus`
- [ ] With: `from uorca.core import TaskManager, TaskStatus`
- [ ] Update function: `run_identification` import (line ~284)
- [ ] Replace: `from uorca.identify import main as identify_main`
- [ ] With: `from uorca.identification import main as identify_main`
- [ ] Verify: `grep "sys.path" uorca/gui/pages/dataset_identification.py` returns nothing

---

### PHASE6-11: Update imports in pipeline_execution.py page
**Effort**: M (30 min)
**Priority**: Critical
**Dependencies**: PHASE6-4
**Files**: `uorca/gui/pages/pipeline_execution.py`

**Tasks**:
- [ ] Remove ALL sys.path manipulation
- [ ] Replace: `from core import TaskManager, TaskStatus`
- [ ] With: `from uorca.core import TaskManager, TaskStatus`
- [ ] Import is already correct: `from uorca.batch.local import LocalBatchProcessor`
- [ ] Verify: No sys.path manipulation

---

### PHASE6-12: Update imports in heatmap_tab.py
**Effort**: M (20 min)
**Priority**: Critical
**Dependencies**: PHASE6-5
**Files**: `uorca/gui/explorer/heatmap_tab.py`

**Tasks**:
- [ ] Remove sys.path manipulation (line ~39)
- [ ] Replace: `from core import validation, script_generation`
- [ ] With: `from uorca.core import validation, script_generation`
- [ ] Replace: `from ResultsIntegration import ResultsIntegrator`
- [ ] With: `from uorca.gui.results_integration import ResultsIntegrator`

---

### PHASE6-13: Update imports in expression_plots_tab.py
**Effort**: M (20 min)
**Priority**: Critical
**Dependencies**: PHASE6-5
**Files**: `uorca/gui/explorer/expression_plots_tab.py`

**Tasks**:
- [ ] Update imports similar to heatmap_tab.py
- [ ] Remove sys.path manipulation
- [ ] Update core imports
- [ ] Update ResultsIntegrator import

---

### PHASE6-14: Update imports in analysis_plots_tab.py
**Effort**: M (20 min)
**Priority**: Critical
**Dependencies**: PHASE6-5
**Files**: `uorca/gui/explorer/analysis_plots_tab.py`

**Tasks**:
- [ ] Update imports similar to previous tabs
- [ ] Remove sys.path manipulation
- [ ] Update all imports to use `uorca.*`

---

### PHASE6-15: Update imports in datasets_info_tab.py
**Effort**: M (20 min)
**Priority**: Critical
**Dependencies**: PHASE6-5
**Files**: `uorca/gui/explorer/datasets_info_tab.py`

**Tasks**:
- [ ] Update imports
- [ ] Remove sys.path manipulation

---

### PHASE6-16: Update imports in contrasts_info_tab.py
**Effort**: M (20 min)
**Priority**: Critical
**Dependencies**: PHASE6-5
**Files**: `uorca/gui/explorer/contrasts_info_tab.py`

**Tasks**:
- [ ] Update imports
- [ ] Remove sys.path manipulation

---

### PHASE6-17: Update imports in ai_assistant_tab.py
**Effort**: M (20 min)
**Priority**: Critical
**Dependencies**: PHASE6-5
**Files**: `uorca/gui/explorer/ai_assistant_tab.py`

**Tasks**:
- [ ] Update imports
- [ ] May have imports from AI modules - update those too

---

### PHASE6-18: Update imports in sidebar_controls.py
**Effort**: M (20 min)
**Priority**: Critical
**Dependencies**: PHASE6-5
**Files**: `uorca/gui/explorer/sidebar_controls.py`

**Tasks**:
- [ ] Update imports
- [ ] Remove sys.path manipulation

---

### PHASE6-19: Update imports in uorca_summary_tab.py
**Effort**: M (20 min)
**Priority**: Critical
**Dependencies**: PHASE6-5
**Files**: `uorca/gui/explorer/uorca_summary_tab.py`

**Tasks**:
- [ ] Update imports
- [ ] Remove sys.path manipulation

---

### PHASE6-20: Update imports in helpers/__init__.py
**Effort**: M (30 min)
**Priority**: Critical
**Dependencies**: PHASE6-6
**Files**: `uorca/gui/explorer/helpers/__init__.py`

**Tasks**:
- [ ] Remove sys.path manipulation
- [ ] Update: `from core import ...` to `from uorca.core import ...`
- [ ] Update: `from ResultsIntegration import ...` to `from uorca.gui.results_integration import ...`
- [ ] Verify: All imports resolve

---

### PHASE6-21: Update imports in uorca_explorer.py (main app)
**Effort**: L (1 hour)
**Priority**: Critical
**Dependencies**: PHASE6-1
**Files**: `uorca/gui/uorca_explorer.py`

**Tasks**:
- [ ] Remove: `sys.path.insert(0, script_dir)` (line ~24)
- [ ] Replace: `from streamlit_tabs import (...)`
- [ ] With: `from uorca.gui.explorer import (...)`
- [ ] Replace: `from ResultsIntegration import ResultsIntegrator`
- [ ] With: `from uorca.gui.results_integration import ResultsIntegrator`
- [ ] Update any other imports from `main_workflow.*`
- [ ] Verify: No sys.path manipulation

---

### PHASE6-22: Update imports in results_integration.py
**Effort**: M (30 min)
**Priority**: Critical
**Dependencies**: PHASE6-2
**Files**: `uorca/gui/results_integration.py`

**Tasks**:
- [ ] Review file for imports (1836 lines, may take time)
- [ ] Update any imports from `core.*` to `uorca.core.*`
- [ ] Update any imports from `main_workflow.*`
- [ ] Remove sys.path manipulation if any
- [ ] Note: NOT refactoring this file, just updating imports

---

### PHASE6-23: Update imports in AI module files
**Effort**: L (1 hour)
**Priority**: High
**Dependencies**: PHASE6-8

**Tasks**:
- [ ] Update: `uorca/gui/ai/agent_factory.py`
- [ ] Update: `uorca/gui/ai/gene_schema.py`
- [ ] Update: `uorca/gui/ai/config_loader.py`
- [ ] Update: `uorca/gui/ai/contrast_relevance.py`
- [ ] Update: `uorca/gui/ai/tool_relevance_analyzer.py`
- [ ] Update: `uorca/gui/ai/mcp_server_core.py`
- [ ] For each: Remove sys.path, update imports to `uorca.*`

---

### PHASE6-24: Update imports in plotting files
**Effort**: M (30 min)
**Priority**: High
**Dependencies**: PHASE6-7

**Tasks**:
- [ ] Update: `uorca/gui/plotting/single_analysis_plots.py`
- [ ] Update: `uorca/gui/plotting/ortholog_mapper.py`
- [ ] Remove sys.path manipulation
- [ ] Update imports

---

### PHASE6-25: Create uorca/gui/__init__.py
**Effort**: S (15 min)
**Priority**: Critical
**Dependencies**: PHASE6-21
**Files**: `uorca/gui/__init__.py`

**Tasks**:
- [ ] Add docstring: `"""UORCA GUI components."""`
- [ ] Import: `from .uorca_explorer import main as launch_explorer`
- [ ] Define __all__: `["launch_explorer"]`
- [ ] Verify: `python -c "from uorca.gui import launch_explorer"` succeeds

---

### PHASE6-26: Create uorca/gui/explorer/__init__.py
**Effort**: M (30 min)
**Priority**: Critical
**Dependencies**: PHASE6-5 through PHASE6-19
**Files**: `uorca/gui/explorer/__init__.py`

**Tasks**:
- [ ] Add docstring: `"""Results explorer tab components."""`
- [ ] Import: `from .heatmap_tab import render_heatmap_tab`
- [ ] Import: `from .expression_plots_tab import render_expression_plots_tab`
- [ ] Import: `from .analysis_plots_tab import render_analysis_plots_tab`
- [ ] Import: `from .datasets_info_tab import render_datasets_info_tab`
- [ ] Import: `from .contrasts_info_tab import render_contrasts_info_tab`
- [ ] Import: `from .ai_assistant_tab import render_ai_assistant_tab`
- [ ] Import: `from .sidebar_controls import render_sidebar_controls`
- [ ] Import: `from .uorca_summary_tab import render_uorca_summary_tab`
- [ ] Define __all__ with all 8 functions
- [ ] Verify: `python -c "from uorca.gui.explorer import render_heatmap_tab"` succeeds

---

### PHASE6-27: Create uorca/gui/pages/__init__.py
**Effort**: S (10 min)
**Priority**: Medium
**Dependencies**: PHASE6-3, PHASE6-4

**Tasks**:
- [ ] Add docstring: `"""Multi-page Streamlit applications."""`
- [ ] Add minimal __all__: `[]` (pages discovered by Streamlit automatically)

---

### PHASE6-28: Create uorca/gui/plotting/__init__.py
**Effort**: S (10 min)
**Priority**: Medium
**Dependencies**: PHASE6-7

**Tasks**:
- [ ] Add docstring: `"""Plotting utilities."""`
- [ ] Import and export plotting functions if needed
- [ ] Can keep empty if not used as import path

---

### PHASE6-29: Create uorca/gui/ai/__init__.py
**Effort**: S (15 min)
**Priority**: Medium
**Dependencies**: PHASE6-8

**Tasks**:
- [ ] Add docstring: `"""AI assistant modules."""`
- [ ] Import and export key AI functions if needed
- [ ] Can keep minimal for now

---

### PHASE6-30: Update uorca/explore.py launcher
**Effort**: M (20 min)
**Priority**: Critical
**Dependencies**: PHASE6-25
**Files**: `uorca/explore.py` (around line 203)

**Tasks**:
- [ ] Find: `streamlit_script = project_root / "main_workflow" / "reporting" / "uorca_explorer.py"`
- [ ] Replace with: `streamlit_script = project_root / "uorca" / "gui" / "uorca_explorer.py"`
- [ ] Verify: Path construction is correct
- [ ] Test: `uv run uorca explore --help` works

---

### PHASE6-31: Verify no sys.path in GUI code
**Effort**: M (20 min)
**Priority**: Critical
**Dependencies**: PHASE6-10 through PHASE6-24

**Tasks**:
- [ ] Run: `grep -r "sys\.path" uorca/gui/ | wc -l`
- [ ] Verify: Returns 0
- [ ] If not: Review each occurrence and fix

---

### PHASE6-32: Test GUI page imports
**Effort**: M (30 min)
**Priority**: Critical
**Dependencies**: PHASE6-10, PHASE6-11

**Tasks**:
- [ ] Test: `python -c "from uorca.gui.pages import dataset_identification"`
- [ ] Test: `python -c "from uorca.gui.pages import pipeline_execution"`
- [ ] Verify: No import errors

---

### PHASE6-33: Test explorer tab imports
**Effort**: M (30 min)
**Priority**: Critical
**Dependencies**: PHASE6-26

**Tasks**:
- [ ] Test: `python -c "from uorca.gui.explorer import render_heatmap_tab"`
- [ ] Test all 8 tab functions
- [ ] Verify: All import successfully

---

### PHASE6-34: Test results_integration import
**Effort**: S (15 min)
**Priority**: Critical
**Dependencies**: PHASE6-22

**Tasks**:
- [ ] Test: `python -c "from uorca.gui.results_integration import ResultsIntegrator"`
- [ ] Verify: Imports successfully (may take time, 1836 lines)

---

### PHASE6-35: Run type checking on GUI
**Effort**: M (30 min)
**Priority**: High
**Dependencies**: PHASE6-31

**Tasks**:
- [ ] Run: `pyright uorca/gui/`
- [ ] Note: May have many errors (Streamlit types, etc.)
- [ ] Document: Significant errors for future cleanup
- [ ] Acceptable: Some type errors OK for now

---

### PHASE6-36: Test explorer launches
**Effort**: M (30 min)
**Priority**: Critical
**Dependencies**: PHASE6-30

**Tasks**:
- [ ] Ensure: Streamlit not running from earlier
- [ ] Run: `uv run uorca explore --help`
- [ ] If works: Try with fake results dir: `uv run uorca explore /tmp`
- [ ] Verify: Streamlit launches (may show error due to invalid dir, that's OK)
- [ ] Check: No import errors in browser

---

### PHASE6-37: Test page navigation
**Effort**: M (30 min)
**Priority**: Critical
**Dependencies**: PHASE6-36

**Tasks**:
- [ ] With explorer running, check sidebar
- [ ] Verify: "Dataset Identification" page appears
- [ ] Verify: "Pipeline Execution" page appears
- [ ] Click: Each page to test it loads
- [ ] Acceptable: Pages may have runtime errors due to missing data
- [ ] Critical: No import errors

---

### PHASE6-38: Test tab rendering
**Effort**: M (30 min)
**Priority**: Critical
**Dependencies**: PHASE6-36

**Tasks**:
- [ ] With valid results dir (if available), load explorer
- [ ] Click through all tabs
- [ ] Verify: Tabs render without import errors
- [ ] Note: May have functional issues, that's OK
- [ ] Critical: Imports work

---

### PHASE6-39: Clean up old reporting directory
**Effort**: S (15 min)
**Priority**: Medium
**Dependencies**: PHASE6-38

**Tasks**:
- [ ] Verify: `main_workflow/reporting/` mostly empty
- [ ] Keep: `__init__.py` for backwards compat (Phase 8)
- [ ] Keep: `core/__init__.py` for backwards compat
- [ ] Remove: Empty `pages/` directory
- [ ] Remove: Empty `streamlit_tabs/` directory
- [ ] Document: What remains

---

### PHASE6-40: Verify file count
**Effort**: S (10 min)
**Priority**: Medium
**Dependencies**: PHASE6-1 through PHASE6-9

**Tasks**:
- [ ] Run: `find uorca/gui -type f -name "*.py" | wc -l`
- [ ] Verify: Count matches expected (~30-35 files)
- [ ] Run: `ls uorca/gui/explorer/*.py | wc -l`
- [ ] Verify: 9 tab files present

---

### PHASE6-41: Commit GUI phase
**Effort**: S (5 min)
**Priority**: Critical
**Dependencies**: PHASE6-40

**Tasks**:
- [ ] Run: `git add -A`
- [ ] Commit: `git commit -m "refactor: move GUI code to uorca/gui/, eliminate sys.path manipulation"`
- [ ] Checkpoint: "Phase 6 complete - GUI moved"

---

### PHASE6-42: Document progress and take break!
**Effort**: S (10 min)
**Priority**: Low
**Dependencies**: PHASE6-41

**Tasks**:
- [ ] Update this TODO: Mark Phase 6 complete
- [ ] Note: This was the longest phase (42 tasks!)
- [ ] Celebrate: Major milestone reached
- [ ] Estimated vs actual time
- [ ] Take break before Phase 7

---

## 🟡 Phase 7: Update Tests (1 day)

**Objective**: Update all test imports to use new package structure.

### PHASE7-1: Update test_task_manager.py imports
**Effort**: S (15 min)
**Priority**: Critical
**Dependencies**: PHASE6-42
**Files**: `tests/unit/test_task_manager.py`

**Tasks**:
- [ ] Replace: `from main_workflow.reporting.core.task_manager import TaskManager, TaskStatus`
- [ ] With: `from uorca.core import TaskManager, TaskStatus`
- [ ] Verify: No other imports from `main_workflow.*`

---

### PHASE7-2: Update test_validation.py imports
**Effort**: S (10 min)
**Priority**: Critical
**Dependencies**: PHASE6-42
**Files**: `tests/unit/test_validation.py`

**Tasks**:
- [ ] Update imports to use `uorca.core`
- [ ] Replace all `main_workflow.reporting.core.*` imports

---

### PHASE7-3: Update test_data_formatters.py imports
**Effort**: S (10 min)
**Priority**: Critical
**Dependencies**: PHASE6-42
**Files**: `tests/unit/test_data_formatters.py`

**Tasks**:
- [ ] Update imports to use `uorca.core`

---

### PHASE7-4: Update test_gene_selection.py imports
**Effort**: S (10 min)
**Priority**: Critical
**Dependencies**: PHASE6-42
**Files**: `tests/unit/test_gene_selection.py`

**Tasks**:
- [ ] Update imports to use `uorca.core`

---

### PHASE7-5: Update test_organism_utils.py imports
**Effort**: S (10 min)
**Priority**: Critical
**Dependencies**: PHASE6-42
**Files**: `tests/unit/test_organism_utils.py`

**Tasks**:
- [ ] Update imports to use `uorca.core`

---

### PHASE7-6: Update test_script_generation.py imports
**Effort**: S (10 min)
**Priority**: Critical
**Dependencies**: PHASE6-42
**Files**: `tests/unit/test_script_generation.py`

**Tasks**:
- [ ] Update imports to use `uorca.core`

---

### PHASE7-7: Move test files to match structure
**Effort**: S (15 min)
**Priority**: High
**Dependencies**: PHASE7-1 through PHASE7-6

**Tasks**:
- [ ] Run: `mv tests/unit/test_task_manager.py tests/unit/core/`
- [ ] Run: `mv tests/unit/test_validation.py tests/unit/core/`
- [ ] Run: `mv tests/unit/test_data_formatters.py tests/unit/core/`
- [ ] Run: `mv tests/unit/test_gene_selection.py tests/unit/core/`
- [ ] Run: `mv tests/unit/test_organism_utils.py tests/unit/core/`
- [ ] Run: `mv tests/unit/test_script_generation.py tests/unit/core/`
- [ ] Verify: `ls tests/unit/core/test_*.py | wc -l` shows 6 files

---

### PHASE7-8: Create placeholder __init__.py files
**Effort**: XS (5 min)
**Priority**: Medium
**Dependencies**: PHASE7-7

**Tasks**:
- [ ] Already exists: `tests/unit/core/__init__.py`
- [ ] Already exists: `tests/unit/identification/__init__.py`
- [ ] Already exists: `tests/unit/pipeline/__init__.py`
- [ ] Already exists: `tests/unit/batch/__init__.py`
- [ ] Already exists: `tests/unit/gui/__init__.py`
- [ ] Verify: All exist

---

### PHASE7-9: Update pytest.ini
**Effort**: S (15 min)
**Priority**: Critical
**Dependencies**: PHASE7-8
**Files**: `pytest.ini`

**Tasks**:
- [ ] Update `testpaths = tests`
- [ ] Update coverage: `--cov=uorca` (instead of `--cov=main_workflow/reporting/core`)
- [ ] Add markers: unit, integration, slow
- [ ] Verify: Configuration valid

---

### PHASE7-10: Run tests
**Effort**: M (30 min)
**Priority**: Critical
**Dependencies**: PHASE7-9

**Tasks**:
- [ ] Run: `pytest tests/` (all tests)
- [ ] Verify: All tests pass
- [ ] If failures: Review and fix import issues
- [ ] Run: `pytest tests/unit/core/` (specific group)
- [ ] Verify: Core tests pass

---

### PHASE7-11: Verify test coverage
**Effort**: M (20 min)
**Priority**: High
**Dependencies**: PHASE7-10

**Tasks**:
- [ ] Run: `pytest --cov=uorca tests/unit/`
- [ ] Verify: Coverage measured for `uorca/` package
- [ ] Check: `uorca/core/` has good coverage (should be ~97%)
- [ ] Note: Other modules have low/no coverage (expected)

---

### PHASE7-12: Commit test updates
**Effort**: XS (2 min)
**Priority**: Critical
**Dependencies**: PHASE7-11

**Tasks**:
- [ ] Run: `git add -A`
- [ ] Commit: `git commit -m "test: update test imports to use uorca package"`
- [ ] Checkpoint: "Phase 7 complete - Tests updated"

---

## 🟡 Phase 8: Create Backwards Compatibility Layer (0.5 day)

**Objective**: Keep `main_workflow/` for backwards compatibility with re-export modules.

### PHASE8-1: Create main_workflow/__init__.py
**Effort**: M (20 min)
**Priority**: Critical
**Dependencies**: PHASE7-12
**Files**: `main_workflow/__init__.py`

**Tasks**:
- [ ] Create file with deprecation warning
- [ ] Import: `from uorca import core, identification, pipeline, batch, gui`
- [ ] Define __all__
- [ ] Test: `python -c "import warnings; warnings.simplefilter('default'); from main_workflow import core"` shows warning

---

### PHASE8-2: Create main_workflow/reporting/__init__.py
**Effort**: M (30 min)
**Priority**: Critical
**Dependencies**: PHASE7-12
**Files**: `main_workflow/reporting/__init__.py`

**Tasks**:
- [ ] Create file with deprecation warning
- [ ] Import: `from uorca import core`
- [ ] Import: `from uorca.gui import results_integration as ResultsIntegration`
- [ ] Create: `task_manager = core.task_manager`
- [ ] Create: `validation = core.validation`
- [ ] Define __all__
- [ ] Test: Import shows warning

---

### PHASE8-3: Create main_workflow/shared/__init__.py
**Effort**: S (15 min)
**Priority**: High
**Dependencies**: PHASE7-12
**Files**: `main_workflow/shared/__init__.py`

**Tasks**:
- [ ] Create file with deprecation warning
- [ ] Import: `from uorca.core import entrez_utils, workflow_logging`
- [ ] Define __all__
- [ ] Test: Import shows warning

---

### PHASE8-4: Create main_workflow/agents/__init__.py
**Effort**: S (15 min)
**Priority**: High
**Dependencies**: PHASE7-12
**Files**: `main_workflow/agents/__init__.py`

**Tasks**:
- [ ] Create file with deprecation warning
- [ ] Import: `from uorca.pipeline.agents import extraction, metadata, analysis`
- [ ] Define __all__
- [ ] Test: Import shows warning

---

### PHASE8-5: Test backwards compatibility imports
**Effort**: M (30 min)
**Priority**: Critical
**Dependencies**: PHASE8-1 through PHASE8-4

**Tasks**:
- [ ] Test: `python -c "import warnings; warnings.filterwarnings('error'); from main_workflow.shared import entrez_utils"` fails with DeprecationWarning (correct!)
- [ ] Test: `python -c "from main_workflow.reporting import ResultsIntegration"` works (with warning)
- [ ] Test: `python -c "from main_workflow.agents import extraction"` works (with warning)
- [ ] Verify: All old import paths work but warn

---

### PHASE8-6: Document backwards compat
**Effort**: S (15 min)
**Priority**: Medium
**Dependencies**: PHASE8-5

**Tasks**:
- [ ] Add comment to each __init__.py: "Will be removed in v1.0"
- [ ] Note: Which imports are supported
- [ ] Document: Migration guide for users

---

### PHASE8-7: Test with warnings as errors
**Effort**: S (10 min)
**Priority**: High
**Dependencies**: PHASE8-6

**Tasks**:
- [ ] Run: `python -W error::DeprecationWarning -c "from main_workflow import core"`
- [ ] Verify: Raises DeprecationWarning (correct behavior)
- [ ] Test: Same for other modules

---

### PHASE8-8: Commit backwards compat
**Effort**: XS (2 min)
**Priority**: Critical
**Dependencies**: PHASE8-7

**Tasks**:
- [ ] Run: `git add -A`
- [ ] Commit: `git commit -m "feat: add backwards compatibility layer for main_workflow imports"`
- [ ] Checkpoint: "Phase 8 complete - Backwards compat added"

---

## 🟡 Phase 9: Update Configuration Files (0.5 day)

**Objective**: Update project configuration to reflect new structure.

### PHASE9-1: Update pyproject.toml
**Effort**: M (30 min)
**Priority**: Critical
**Dependencies**: PHASE8-8
**Files**: `pyproject.toml`

**Tasks**:
- [ ] Update version: `version = "0.2.0"` (bump for restructure)
- [ ] Verify: `[project.scripts]` still has `uorca = "uorca.cli:main"`
- [ ] Update: `[tool.setuptools.packages.find]` to include only `uorca*`
- [ ] Add: `exclude = ["tests*", "main_workflow*"]`
- [ ] Update: `[tool.pytest.ini_options]` coverage to `--cov=uorca`
- [ ] Add: `[tool.pyright]` config with `include = ["uorca", "tests"]`, `exclude = ["main_workflow"]`
- [ ] Verify: File is valid TOML

---

### PHASE9-2: Test package build
**Effort**: S (15 min)
**Priority**: Critical
**Dependencies**: PHASE9-1

**Tasks**:
- [ ] Run: `uv build`
- [ ] Verify: Build succeeds
- [ ] Check: Package includes `uorca/` but excludes `main_workflow/`
- [ ] Inspect: Built wheel/tarball if needed

---

### PHASE9-3: Update CLAUDE.md structure section
**Effort**: M (30 min)
**Priority**: High
**Dependencies**: PHASE8-8
**Files**: `CLAUDE.md`

**Tasks**:
- [ ] Update "File Organization" section
- [ ] Document new structure: `uorca/core/`, `uorca/identification/`, etc.
- [ ] Update "Code Conventions" import patterns
- [ ] Add examples: `from uorca.core import TaskManager`
- [ ] Note: `main_workflow/` is deprecated
- [ ] Update file naming conventions

---

### PHASE9-4: Update CLAUDE.md import examples
**Effort**: S (15 min)
**Priority**: High
**Dependencies**: PHASE9-3
**Files**: `CLAUDE.md`

**Tasks**:
- [ ] Add section: "✅ CORRECT imports"
- [ ] Add section: "❌ WRONG imports (no sys.path)"
- [ ] Show examples for each module
- [ ] Note: No more sys.path manipulation needed

---

### PHASE9-5: Update .gitignore if needed
**Effort**: S (10 min)
**Priority**: Low
**Dependencies**: PHASE8-8
**Files**: `.gitignore`

**Tasks**:
- [ ] Check if needs update
- [ ] Add: `uorca/**/__pycache__/` (may already be covered)
- [ ] Verify: `.gitignore` still makes sense

---

### PHASE9-6: Test CLI entry point
**Effort**: S (10 min)
**Priority**: Critical
**Dependencies**: PHASE9-2

**Tasks**:
- [ ] Run: `uv run uorca --version` (if version command exists)
- [ ] Run: `uv run uorca --help`
- [ ] Verify: CLI works correctly
- [ ] Test: `uorca` (if installed globally)

---

### PHASE9-7: Run type checking on full package
**Effort**: M (20 min)
**Priority**: High
**Dependencies**: PHASE9-1

**Tasks**:
- [ ] Run: `pyright uorca/`
- [ ] Verify: Pyright uses new config from pyproject.toml
- [ ] Review: Type errors (many expected, that's OK)
- [ ] Verify: No critical import errors

---

### PHASE9-8: Verify README.md still accurate
**Effort**: M (20 min)
**Priority**: Medium
**Dependencies**: PHASE9-3

**Tasks**:
- [ ] Read: `README.md` installation instructions
- [ ] Verify: Commands still work (uv sync, etc.)
- [ ] Check: No references to old structure
- [ ] Update: If needed

---

### PHASE9-9: Commit config updates
**Effort**: XS (2 min)
**Priority**: Critical
**Dependencies**: PHASE9-8

**Tasks**:
- [ ] Run: `git add -A`
- [ ] Commit: `git commit -m "chore: update project config for restructured package"`
- [ ] Checkpoint: "Phase 9 complete - Config updated"

---

## 🔴 Phase 10: Final Verification (1 day)

**Objective**: Comprehensive testing of restructured codebase.

### PHASE10-1: Run full test suite
**Effort**: M (20 min)
**Priority**: Critical
**Dependencies**: PHASE9-9

**Tasks**:
- [ ] Run: `pytest tests/`
- [ ] Verify: All tests pass
- [ ] Check: Coverage report generated
- [ ] Review: Any unexpected failures

---

### PHASE10-2: Run type checking
**Effort**: M (20 min)
**Priority**: Critical
**Dependencies**: PHASE9-9

**Tasks**:
- [ ] Run: `pyright uorca/`
- [ ] Review: Errors (many OK, note critical ones)
- [ ] Verify: No import errors
- [ ] Document: Significant issues for future

---

### PHASE10-3: Test package build
**Effort**: S (10 min)
**Priority**: Critical
**Dependencies**: PHASE9-9

**Tasks**:
- [ ] Run: `uv build`
- [ ] Verify: Succeeds
- [ ] Check: Output in `dist/`

---

### PHASE10-4: Test uorca --help
**Effort**: XS (2 min)
**Priority**: Critical
**Dependencies**: PHASE9-9

**Tasks**:
- [ ] Run: `uv run uorca --help`
- [ ] Verify: Shows help

---

### PHASE10-5: Test uorca identify --help
**Effort**: XS (2 min)
**Priority**: Critical
**Dependencies**: PHASE9-9

**Tasks**:
- [ ] Run: `uv run uorca identify --help`
- [ ] Verify: Shows help

---

### PHASE10-6: Test uorca run --help
**Effort**: XS (2 min)
**Priority**: Critical
**Dependencies**: PHASE9-9

**Tasks**:
- [ ] Run: `uv run uorca run --help`
- [ ] Verify: Shows help

---

### PHASE10-7: Test uorca explore --help
**Effort**: XS (2 min)
**Priority**: Critical
**Dependencies**: PHASE9-9

**Tasks**:
- [ ] Run: `uv run uorca explore --help`
- [ ] Verify: Shows help

---

### PHASE10-8: Verify no sys.path in uorca/
**Effort**: S (5 min)
**Priority**: Critical
**Dependencies**: PHASE9-9

**Tasks**:
- [ ] Run: `grep -r "sys\.path\.insert" uorca/ | wc -l`
- [ ] Verify: Returns 0
- [ ] If not 0: Find and fix remaining occurrences

---

### PHASE10-9: Test identify command (if API keys available)
**Effort**: M (30 min)
**Priority**: High
**Dependencies**: PHASE10-5

**Tasks**:
- [ ] If have API keys: `uv run uorca identify -q "test query" -o /tmp/test_identify`
- [ ] Verify: Command runs (may fail on API, that's OK)
- [ ] Check: No import errors
- [ ] If no API keys: Skip this task

---

### PHASE10-10: Test explore command
**Effort**: M (30 min)
**Priority**: Critical
**Dependencies**: PHASE10-7

**Tasks**:
- [ ] Run: `uv run uorca explore /tmp` (invalid dir is OK)
- [ ] Verify: Streamlit launches
- [ ] Check browser: No import errors
- [ ] Close: Streamlit server

---

### PHASE10-11: Test pages in explorer
**Effort**: M (30 min)
**Priority**: Critical
**Dependencies**: PHASE10-10

**Tasks**:
- [ ] Launch: `uv run uorca explore /tmp` or with valid results
- [ ] Navigate to: "Dataset Identification" page
- [ ] Navigate to: "Pipeline Execution" page
- [ ] Verify: Both pages load without import errors
- [ ] Acceptable: Runtime errors due to missing data
- [ ] Close: Server

---

### PHASE10-12: Test tabs in explorer (if results available)
**Effort**: M (30 min)
**Priority**: High
**Dependencies**: PHASE10-10

**Tasks**:
- [ ] If have valid results dir: Launch explorer with it
- [ ] Click through all tabs: Heatmap, Expression, Analysis, etc.
- [ ] Verify: Tabs render without import errors
- [ ] Check: Data displays correctly
- [ ] If no results: Skip detailed testing

---

### PHASE10-13: Test run command help
**Effort**: S (10 min)
**Priority**: High
**Dependencies**: PHASE10-6

**Tasks**:
- [ ] Run: `uv run uorca run local --help`
- [ ] Run: `uv run uorca run slurm --help`
- [ ] Verify: Both show help
- [ ] Check: No import errors

---

### PHASE10-14: Test background task manager in GUI
**Effort**: M (30 min)
**Priority**: High
**Dependencies**: PHASE10-11

**Tasks**:
- [ ] Launch: Explorer or identification page
- [ ] Try: Submit a task (identification or pipeline)
- [ ] Verify: TaskManager works
- [ ] Check: Task appears in history
- [ ] Verify: SQLite database created in ~/.uorca/

---

### PHASE10-15: Test legacy imports with warnings
**Effort**: M (20 min)
**Priority**: Critical
**Dependencies**: PHASE9-9

**Tasks**:
- [ ] Run: `python -c "import warnings; warnings.simplefilter('default'); from main_workflow.reporting import ResultsIntegration"` (should show warning)
- [ ] Verify: Warning displayed
- [ ] Verify: Import works
- [ ] Test: Other legacy imports

---

### PHASE10-16: Test import performance
**Effort**: S (10 min)
**Priority**: Low
**Dependencies**: PHASE10-4

**Tasks**:
- [ ] Run: `time python -c "from uorca.core import TaskManager"`
- [ ] Note: Import time (should be fast, <500ms)
- [ ] Compare: With old structure if possible
- [ ] Acceptable: Similar or faster

---

### PHASE10-17: Review CLAUDE.md accuracy
**Effort**: M (20 min)
**Priority**: Medium
**Dependencies**: PHASE9-9

**Tasks**:
- [ ] Read: CLAUDE.md
- [ ] Verify: Structure section matches reality
- [ ] Verify: Import examples are correct
- [ ] Check: Conventions documented
- [ ] Update: If inaccuracies found

---

### PHASE10-18: Review README.md accuracy
**Effort**: S (15 min)
**Priority**: Medium
**Dependencies**: PHASE9-9

**Tasks**:
- [ ] Read: README.md
- [ ] Verify: Installation steps work
- [ ] Check: Examples use correct commands
- [ ] Update: If needed

---

### PHASE10-19: Verify docstrings updated
**Effort**: M (20 min)
**Priority**: Low
**Dependencies**: PHASE9-9

**Tasks**:
- [ ] Sample check: Read docstrings in moved files
- [ ] Verify: No references to old paths
- [ ] Check: Import examples in docstrings correct
- [ ] Note: Any widespread issues

---

### PHASE10-20: Run coverage report
**Effort**: M (20 min)
**Priority**: Medium
**Dependencies**: PHASE10-1

**Tasks**:
- [ ] Run: `pytest --cov=uorca tests/ --cov-report=html`
- [ ] Open: `htmlcov/index.html`
- [ ] Review: Coverage by module
- [ ] Note: `uorca/core/` should have ~97%
- [ ] Note: Other modules have low coverage (expected, infrastructure only)

---

### PHASE10-21: Final git status check
**Effort**: S (5 min)
**Priority**: High
**Dependencies**: PHASE10-20

**Tasks**:
- [ ] Run: `git status`
- [ ] Verify: No unexpected untracked files
- [ ] Verify: No unexpected modified files
- [ ] Clean: Any temporary files

---

### PHASE10-22: Create final checkpoint
**Effort**: S (5 min)
**Priority**: Critical
**Dependencies**: PHASE10-21

**Tasks**:
- [ ] Run: `git add -A` (if any final changes)
- [ ] Commit: `git commit -m "docs: update documentation post-restructure"` (if changes)
- [ ] Checkpoint: "Phase 10 complete - Restructure verified"
- [ ] Tag: `git tag v0.2.0-restructure` (optional)

---

### PHASE10-23: Document completion
**Effort**: M (30 min)
**Priority**: High
**Dependencies**: PHASE10-22

**Tasks**:
- [ ] Update: This TODO list - mark all complete!
- [ ] Calculate: Total actual time vs estimated
- [ ] Document: Lessons learned
- [ ] Note: Any issues encountered
- [ ] Create: Summary of changes for team/users
- [ ] Celebrate: Major restructure complete! 🎉

---

## 📈 Progress Tracking

**Overall Progress**: 0/87 tasks (0%)

| Phase | Tasks | Status | Est. Time | Actual Time |
|-------|-------|--------|-----------|-------------|
| Phase 1: Structure | 8 | Not Started | 2h | - |
| Phase 2: Core | 13 | Not Started | 1d | - |
| Phase 3: Identification | 15 | Not Started | 1d | - |
| Phase 4: Pipeline | 18 | Not Started | 1.5d | - |
| Phase 5: Batch Helpers | 7 | Not Started | 0.5d | - |
| Phase 6: GUI | 42 | Not Started | 3d | - |
| Phase 7: Tests | 12 | Not Started | 1d | - |
| Phase 8: Backwards Compat | 8 | Not Started | 0.5d | - |
| Phase 9: Config | 9 | Not Started | 0.5d | - |
| Phase 10: Verification | 23 | Not Started | 1d | - |
| **Total** | **155** | **0%** | **~12d** | **0d** |

---

## 📊 Velocity Tracking

**Day 1**: [ ] tasks completed
**Day 2**: [ ] tasks completed
**Day 3**: [ ] tasks completed
**Day 4**: [ ] tasks completed
**Day 5**: [ ] tasks completed

**Average velocity**: [Calculate after Day 1]
**Adjusted timeline**: [Update if needed]

---

## 🎯 Milestones

### Milestone 1: Foundation (End of Week 1)
**Target**: Day 5
**Tasks**: Phases 1-5 complete (61 tasks)
**Success**: Core, identification, pipeline, batch all moved and working

### Milestone 2: GUI Complete (End of Day 7)
**Target**: Day 7
**Tasks**: Phase 6 complete (42 tasks)
**Success**: All GUI code moved, no sys.path manipulation, explorer works

### Milestone 3: Restructure Complete (End of Week 2)
**Target**: Day 10
**Tasks**: Phases 7-10 complete (52 tasks)
**Success**: Tests pass, config updated, backwards compat works, all verified

---

## 🐛 Known Issues / Blockers

- [Track issues as they arise]
- [Reference GitHub issues if applicable]

---

## ⚠️ Risks & Concerns

**Technical Risks**:
1. **Import chains**: Complex import dependencies may cause circular imports
   - **Mitigation**: Test imports after each phase
2. **ResultsIntegration.py**: 1836 lines, may have many dependencies
   - **Mitigation**: Update imports carefully, test thoroughly
3. **GUI may break**: Many Streamlit components to update
   - **Mitigation**: Test in browser after Phase 6

**Timeline Risks**:
1. **Phase 6 is large**: 42 tasks, may take longer than estimated
   - **Mitigation**: Can split into sub-phases if needed
2. **Testing may reveal issues**: Phase 10 may uncover problems
   - **Mitigation**: Checkpoints allow rollback

---

## 💡 Tips for Implementation

- **Use checkpoints liberally**: Create checkpoint after each phase
- **Test incrementally**: Don't wait until end to test imports
- **Use `/rewind` if needed**: Rollback if something breaks
- **Update this TODO**: Mark tasks complete as you go
- **Take breaks**: Phase 6 is long, pace yourself
- **Ask if stuck**: Better to ask than guess

---

## 🔧 Quick Commands Reference

```bash
# Verify structure
find uorca -name "*.py" | head -20

# Check for sys.path
grep -r "sys\.path" uorca/

# Test imports
python -c "from uorca.core import TaskManager"

# Run tests
pytest tests/

# Type checking
pyright uorca/

# Build package
uv build

# Launch explorer
uv run uorca explore /tmp
```

---

## 📚 References

- **Original Plan**: `docs/plans/codebase_restructure.md`
- **Project Context**: `CLAUDE.md`
- **Test Patterns**: `tests/conftest.py`
- **No-CLI Plan**: `docs/plans/no_cli_uorca_implementation.md`

---

*This TODO list is a living document. Update it as work progresses. Good luck with the restructure!*
