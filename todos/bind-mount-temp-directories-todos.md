# TODO List: Bind-Mount Temp Directories for Docker

**Generated from**: `docs/plans/bind_mount_temp_directories_for_docker.md`
**Generated on**: 2025-10-02
**Total Tasks**: 15 core tasks + 4 test cases
**Estimated Effort**: 1-2 hours implementation + 30 min testing
**Target Completion**: Same day

---

## 📋 Quick Summary

- **Critical tasks**: 3 (Phase 1 & 2 implementation)
- **High priority**: 4 (Testing and verification)
- **Medium priority**: 2 (Documentation updates)
- **Total estimated time**: 1.5-2.5 hours

**Goal**: Fix Docker disk space errors by bind-mounting temp directories to host filesystem, allowing UORCA to process large datasets with default Docker Desktop settings (50GB).

---

## 🎯 Project Goals

Fix the "disk-limit exceeded" error that occurs when processing large RNA-seq datasets locally. Currently, fasterq-dump writes temp files to Docker's virtual disk (50GB limit), causing failures even when the host has 200GB+ free space. Solution: bind-mount temp directories to host filesystem, matching SLURM's existing behavior.

## ✅ Success Criteria

**Must achieve**:
- [ ] GSE103111 dataset (50GB peak) processes successfully with 50GB Docker allocation
- [ ] No "disk-limit exceeded" or "No space left on device" errors
- [ ] Temp directories automatically created and cleaned up
- [ ] SLURM execution completely unaffected (zero changes to SLURM code)

## 🚫 Out of Scope

- ❌ Adding `--temp_dir` CLI parameter (future enhancement)
- ❌ Modifying SLURM execution (already works correctly)
- ❌ Changing tmpfs mount configuration
- ❌ Adding Docker disk space checking logic
- ❌ Modifying base `BatchProcessor` class

---

## 📊 Progress Tracking

**Overall Progress**: 0/19 tasks (0%)

| Phase | Tasks | Status |
|-------|-------|--------|
| Phase 1: Bind Mount Setup | 0/7 | Not Started |
| Phase 2: Cleanup Logic | 0/4 | Not Started |
| Phase 3: Integration Testing | 0/4 | Not Started |
| Documentation Updates | 0/2 | Not Started |
| Regression Verification | 0/2 | Not Started |

---

## 🔴 Phase 1: Add Temp Directory Bind Mount (Critical)

### ✅ PHASE1-1: Create Git Checkpoint

**Priority**: Critical
**Effort**: XS (2 min)
**Dependencies**: None
**Skills**: Git
**Assigned**: [Unassigned]

**Description**:
Create initial checkpoint before making any code changes to enable easy rollback if needed.

**Tasks**:
- [ ] Commit current state with message: "Checkpoint: Before Phase 1 - Bind-mount temp directory"
- [ ] Note current commit hash for reference
- [ ] Verify working directory is clean: `git status`

**Acceptance Criteria**:
- [ ] Git checkpoint created successfully
- [ ] Working directory is clean
- [ ] Commit hash recorded

**Commands**:
```bash
git add -A
git commit -m "Checkpoint: Before Phase 1 - Bind-mount temp directory"
git rev-parse HEAD  # Note this hash
```

---

### ✅ PHASE1-2: Read and Understand Current Implementation

**Priority**: Critical
**Effort**: XS (5 min)
**Dependencies**: None
**Skills**: Python, Docker
**Assigned**: [Unassigned]

**Description**:
Review the current `run_single_dataset_local()` function to understand Docker command construction and identify exact insertion points for temp directory code.

**Tasks**:
- [ ] Read `uorca/batch/local.py:27-219` (full function)
- [ ] Identify where directories are created (around line 79)
- [ ] Identify where Docker volumes are mounted (around line 96-99)
- [ ] Identify where environment variables are set (around line 104-110)
- [ ] Note the function structure for try/finally wrapping

**Acceptance Criteria**:
- [ ] Understand current Docker command construction
- [ ] Identified all three insertion points
- [ ] Ready to make changes

**Files to Read**:
- `uorca/batch/local.py:27-219`

---

### ✅ PHASE1-3: Create Temp Directory on Host

**Priority**: Critical
**Effort**: S (5 min)
**Dependencies**: PHASE1-2
**Skills**: Python, Pathlib
**Assigned**: [Unassigned]

**Description**:
Add code to create a per-dataset temp directory on the host filesystem at `{output_dir}/{accession}/tmp/`.

**Location**: `uorca/batch/local.py:run_single_dataset_local()`
**Insert After**: Line 79 (after `resource_path.mkdir(parents=True, exist_ok=True)`)

**Code to Add**:
```python
# Create temp directory for large intermediate files (bind-mounted to bypass Docker disk limits)
temp_dir = output_path / accession / "tmp"
temp_dir.mkdir(parents=True, exist_ok=True)
```

**Acceptance Criteria**:
- [ ] Code inserted at correct location
- [ ] Uses `Path` object from existing code
- [ ] Creates directory with `parents=True, exist_ok=True`
- [ ] Type checking passes: `pyright uorca/batch/local.py`

**Files**:
- `uorca/batch/local.py:~79`

**Testing**:
```bash
# Verify syntax
python -m py_compile uorca/batch/local.py
pyright uorca/batch/local.py
```

---

### ✅ PHASE1-4: Add Temp Directory to Docker Mounts

**Priority**: Critical
**Effort**: S (5 min)
**Dependencies**: PHASE1-3
**Skills**: Python, Docker
**Assigned**: [Unassigned]

**Description**:
Add bind mount for temp directory to Docker volume mounts, mapping host temp dir to `/workspace/temp` in container.

**Location**: `uorca/batch/local.py:run_single_dataset_local()`
**Insert After**: Line 99 (in volume mounts section)

**Code to Add**:
```python
# Mount volumes
cmd += [
    '-v', f'{output_path}:{docker_output_dir}',
    '-v', f'{resource_path}:{docker_resource_dir}',
    '-v', f'{project_root}:/workspace/src',
    '-v', f'{temp_dir}:/workspace/temp',  # ← NEW: Bind-mount temp directory
    '--workdir', '/workspace/src',
]
```

**Acceptance Criteria**:
- [ ] Temp directory mount added to cmd list
- [ ] Mount path is `/workspace/temp` (not `/tmp`)
- [ ] Variable name `temp_dir` matches previous step
- [ ] Type checking passes

**Files**:
- `uorca/batch/local.py:~96-100`

**Testing**:
```bash
pyright uorca/batch/local.py
```

---

### ✅ PHASE1-5: Set Environment Variables for Temp Directory

**Priority**: Critical
**Effort**: S (5 min)
**Dependencies**: PHASE1-4
**Skills**: Python, Unix/Linux
**Assigned**: [Unassigned]

**Description**:
Set `TMPDIR`, `TEMP`, and `TMP` environment variables to guide bioinformatics tools to use the bind-mounted temp directory.

**Location**: `uorca/batch/local.py:run_single_dataset_local()`
**Insert After**: Line 110 (in environment variables section, after `VIRTUAL_ENV`)

**Code to Add**:
```python
# Pass environment variables
cmd += [
    '-e', 'UV_NO_SYNC=1',
    '-e', 'PATH=/workspace/.venv/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin',
    '-e', 'VIRTUAL_ENV=/workspace/.venv',
    '-e', 'TMPDIR=/workspace/temp',    # ← NEW: Standard Unix temp dir
    '-e', 'TEMP=/workspace/temp',      # ← NEW: Alternative name
    '-e', 'TMP=/workspace/temp',       # ← NEW: Windows compatibility
    '-e', f'ENTREZ_EMAIL={os.getenv("ENTREZ_EMAIL", "")}',
    '-e', f'OPENAI_API_KEY={os.getenv("OPENAI_API_KEY", "")}',
    '-e', f'ENTREZ_API_KEY={os.getenv("ENTREZ_API_KEY", "")}',
]
```

**Why These Variables**:
- `TMPDIR`: Standard Unix/Linux environment variable (fasterq-dump uses this)
- `TEMP`: Alternative checked by some tools
- `TMP`: Windows compatibility (some cross-platform tools check this)

**Acceptance Criteria**:
- [ ] All three environment variables added
- [ ] Variables point to `/workspace/temp` (container path)
- [ ] Existing environment variables unchanged
- [ ] Type checking passes

**Files**:
- `uorca/batch/local.py:~104-112`

**Testing**:
```bash
pyright uorca/batch/local.py
```

---

### ✅ PHASE1-6: Verify Phase 1 Implementation

**Priority**: Critical
**Effort**: M (10 min)
**Dependencies**: PHASE1-3, PHASE1-4, PHASE1-5
**Skills**: Python, Testing
**Assigned**: [Unassigned]

**Description**:
Run automated checks and a quick manual test to verify Phase 1 changes are syntactically correct and temp directory is being created.

**Automated Checks**:
- [ ] Type checking: `pyright uorca/batch/local.py`
- [ ] Syntax check: `python -m py_compile uorca/batch/local.py`
- [ ] Verify no regression in imports

**Manual Test** (optional, can defer to Phase 3):
```bash
# Quick test with a small dataset
uorca run local \
  --input test_datasets/tiny_test.csv \
  --output scratch/phase1_test/ \
  --max_workers 1

# Check if temp directory was created
ls -la scratch/phase1_test/GSE*/tmp/
# Should show temp directory exists (cleanup not added yet)
```

**Acceptance Criteria**:
- [ ] All automated checks pass
- [ ] No import errors
- [ ] No type errors
- [ ] (Optional) Temp directory created during test run

**Note**: At this phase, temp directories will NOT be cleaned up automatically. That's expected - cleanup is added in Phase 2.

---

### ✅ PHASE1-7: Create Phase 1 Checkpoint

**Priority**: Critical
**Effort**: XS (2 min)
**Dependencies**: PHASE1-6
**Skills**: Git
**Assigned**: [Unassigned]

**Description**:
Commit Phase 1 changes as a checkpoint before moving to Phase 2.

**Tasks**:
- [ ] Stage changes: `git add uorca/batch/local.py`
- [ ] Review diff: `git diff --staged`
- [ ] Commit with descriptive message
- [ ] Note commit hash

**Acceptance Criteria**:
- [ ] Commit created with clear message
- [ ] Only `local.py` modified
- [ ] Commit hash recorded for potential rollback

**Commands**:
```bash
git add uorca/batch/local.py
git diff --staged  # Review changes
git commit -m "feat: add bind-mounted temp directories for Docker local execution

- Create per-dataset temp directory on host filesystem
- Bind-mount to /workspace/temp in container
- Set TMPDIR, TEMP, TMP environment variables
- Bypasses Docker's virtual disk limits

Refs: docs/plans/bind_mount_temp_directories_for_docker.md (Phase 1)"
git rev-parse HEAD  # Note hash
```

---

## 🟡 Phase 2: Add Cleanup Logic (High Priority)

### ✅ PHASE2-1: Create Phase 2 Checkpoint

**Priority**: High
**Effort**: XS (2 min)
**Dependencies**: PHASE1-7
**Skills**: Git
**Assigned**: [Unassigned]

**Description**:
Create checkpoint before Phase 2 to enable rollback to "bind mount working, no cleanup" state.

**Commands**:
```bash
git commit --allow-empty -m "Checkpoint: Before Phase 2 - Cleanup logic"
git rev-parse HEAD
```

**Acceptance Criteria**:
- [ ] Checkpoint commit created
- [ ] Commit hash recorded

---

### ✅ PHASE2-2: Add Shutil Import

**Priority**: High
**Effort**: XS (1 min)
**Dependencies**: PHASE2-1
**Skills**: Python
**Assigned**: [Unassigned]

**Description**:
Add `import shutil` to function imports for `shutil.rmtree()` usage in cleanup.

**Location**: `uorca/batch/local.py:run_single_dataset_local()`
**Add to imports** (around line 46-50):

**Code to Add**:
```python
import subprocess
import json
import shutil  # ← NEW: For directory removal
from datetime import datetime
from pathlib import Path
import re
```

**Acceptance Criteria**:
- [ ] `import shutil` added
- [ ] Import placed logically with other imports
- [ ] Type checking passes

**Files**:
- `uorca/batch/local.py:~46`

---

### ✅ PHASE2-3: Wrap Function Body in Try/Finally

**Priority**: High
**Effort**: M (15 min)
**Dependencies**: PHASE2-2
**Skills**: Python, Error Handling
**Assigned**: [Unassigned]

**Description**:
Wrap the entire processing logic in a try/finally block to ensure temp directory cleanup runs even on errors.

**Location**: `uorca/batch/local.py:run_single_dataset_local()`
**Structure Change**: Lines 52-218

**Before**:
```python
def run_single_dataset_local(...):
    """docstring"""
    import subprocess
    ...

    # Setup paths and create directories
    ...

    try:
        # Docker execution
        result = subprocess.run(cmd, ...)
        ...
    except subprocess.TimeoutExpired:
        ...
    except Exception as e:
        ...
```

**After**:
```python
def run_single_dataset_local(...):
    """docstring"""
    import subprocess
    import shutil
    ...

    # Setup paths (before try block)
    output_path = Path(output_dir).absolute()
    resource_path = Path(resource_dir).absolute()
    current_file = Path(__file__).absolute()
    project_root = current_file.parent.parent.parent

    # Create directories
    output_path.mkdir(parents=True, exist_ok=True)
    resource_path.mkdir(parents=True, exist_ok=True)

    # Create temp directory
    temp_dir = output_path / accession / "tmp"
    temp_dir.mkdir(parents=True, exist_ok=True)

    try:
        # All existing processing logic (status updates, Docker execution, etc.)
        status_file = Path(status_dir) / f"{accession}_status.json"
        ...
        # (existing code unchanged)
        ...

    except subprocess.TimeoutExpired:
        # (existing timeout handling unchanged)
        ...

    except Exception as e:
        # (existing error handling unchanged)
        ...

    finally:
        # ← NEW: Clean up temp directory regardless of success/failure
        if temp_dir.exists():
            try:
                shutil.rmtree(temp_dir)
                # Note: Logger not available in this process, so no logging here
            except Exception as cleanup_error:
                # Cleanup failure is non-critical, just skip
                pass
```

**Key Points**:
- Move path setup and temp directory creation BEFORE the try block
- Wrap all processing logic in try/except/finally
- Add cleanup in finally block
- Silently ignore cleanup errors (non-critical)

**Acceptance Criteria**:
- [ ] Path setup moved before try block
- [ ] All processing logic inside try block
- [ ] Cleanup logic in finally block
- [ ] Cleanup errors caught and ignored
- [ ] Type checking passes
- [ ] Function still returns correct values

**Files**:
- `uorca/batch/local.py:27-219` (entire function)

**Testing**:
```bash
pyright uorca/batch/local.py
python -m py_compile uorca/batch/local.py
```

**Caution**: This is the trickiest change. Be careful with indentation and ensure all code paths still return appropriate values.

---

### ✅ PHASE2-4: Update Function Docstring

**Priority**: High
**Effort**: S (5 min)
**Dependencies**: PHASE2-3
**Skills**: Documentation
**Assigned**: [Unassigned]

**Description**:
Update the function docstring to document temp directory behavior for future developers.

**Location**: `uorca/batch/local.py:run_single_dataset_local()`
**Update docstring** (lines 30-44):

**New Docstring**:
```python
def run_single_dataset_local(accession: str, output_dir: str, resource_dir: str,
                           cleanup: bool, status_dir: str, docker_image: str = "kevingchen/uorca:0.1.0",
                           container_tmpfs_gb: int = 20) -> Dict[str, Any]:
    """
    Run a single dataset analysis locally using Docker.

    This function is designed to be run in a separate process. It creates a temporary
    directory on the host filesystem ({output_dir}/{accession}/tmp/) and bind-mounts it
    to the container to bypass Docker's virtual disk limits. This allows processing of
    large datasets regardless of Docker Desktop's disk allocation.

    Args:
        accession: Dataset accession ID
        output_dir: Output directory on host
        resource_dir: Resource directory for Kallisto indices on host
        cleanup: Whether to cleanup intermediate files (FASTQ/SRA after analysis)
        status_dir: Directory for status tracking on host
        docker_image: Docker image to use
        container_tmpfs_gb: Size in GB for in-memory tmpfs (for small temp files)

    Returns:
        Dictionary with job results including container_id

    Note:
        The temp directory ({output_dir}/{accession}/tmp/) is automatically created
        and cleaned up (even on failure). Temporary files from fasterq-dump, pigz, etc.
        write to this bind-mounted directory on the host, not to Docker's virtual disk.
    """
```

**Acceptance Criteria**:
- [ ] Docstring updated to describe temp directory behavior
- [ ] Note added explaining automatic cleanup
- [ ] Args and Returns sections preserved

**Files**:
- `uorca/batch/local.py:30-44`

---

### ✅ PHASE2-5: Verify Phase 2 Implementation

**Priority**: High
**Effort**: M (10 min)
**Dependencies**: PHASE2-2, PHASE2-3, PHASE2-4
**Skills**: Python, Testing
**Assigned**: [Unassigned]

**Description**:
Verify Phase 2 changes work correctly with automated checks and manual testing.

**Automated Checks**:
- [ ] Type checking: `pyright uorca/batch/local.py`
- [ ] Syntax check: `python -m py_compile uorca/batch/local.py`

**Manual Tests**:

**Test 1: Successful Cleanup**
```bash
# Run with small dataset
uorca run local \
  --input test_datasets/tiny_test.csv \
  --output scratch/phase2_success/ \
  --max_workers 1

# After completion, verify NO temp directory
ls scratch/phase2_success/GSE*/tmp/
# Should return "No such file or directory"
```

**Test 2: Cleanup on Failure**
```bash
# Run with invalid accession (will fail)
echo "Accession,DatasetSizeGB" > /tmp/invalid.csv
echo "INVALID_GSE123,5" >> /tmp/invalid.csv

uorca run local \
  --input /tmp/invalid.csv \
  --output scratch/phase2_failure/ \
  --max_workers 1

# After failure, verify NO temp directory
ls scratch/phase2_failure/INVALID_GSE123/tmp/
# Should return "No such file or directory"
```

**Acceptance Criteria**:
- [ ] All automated checks pass
- [ ] Temp directory cleaned up after success
- [ ] Temp directory cleaned up after failure
- [ ] No leaked temp directories

---

### ✅ PHASE2-6: Create Phase 2 Checkpoint

**Priority**: High
**Effort**: XS (2 min)
**Dependencies**: PHASE2-5
**Skills**: Git
**Assigned**: [Unassigned]

**Description**:
Commit Phase 2 changes as checkpoint before integration testing.

**Commands**:
```bash
git add uorca/batch/local.py
git diff --staged
git commit -m "feat: add automatic cleanup for temp directories

- Wrap processing in try/finally block
- Clean up temp directory in finally (runs even on error)
- Update docstring to document temp directory behavior
- Use shutil.rmtree() for cleanup

Refs: docs/plans/bind_mount_temp_directories_for_docker.md (Phase 2)"
git rev-parse HEAD
```

**Acceptance Criteria**:
- [ ] Commit created
- [ ] Only `local.py` modified
- [ ] Commit hash recorded

---

## 🟢 Phase 3: Integration Testing (High Priority)

### ✅ PHASE3-1: Test Case 1 - Large Dataset with Default Docker

**Priority**: High
**Effort**: L (20 min)
**Dependencies**: PHASE2-6
**Skills**: Testing, Bioinformatics
**Assigned**: [Unassigned]

**Description**:
Run GSE103111 (the dataset that originally failed) with default 50GB Docker allocation to verify fix works.

**Setup**:
- Ensure Docker Desktop is set to default 50GB allocation
- Clear any previous test outputs

**Test Command**:
```bash
# Create CSV for GSE103111
echo "Accession,DatasetSizeGB" > test_datasets/GSE103111.csv
echo "GSE103111,20" >> test_datasets/GSE103111.csv

# Run the pipeline
uorca run local \
  --input test_datasets/GSE103111.csv \
  --output scratch/docker_disk_test/ \
  --max_workers 1
```

**Verification Checklist**:
- [ ] Pipeline completes successfully without errors
- [ ] No "disk-limit exceeded" errors in logs
- [ ] No "No space left on device" errors
- [ ] Check logs: `cat scratch/docker_disk_test/GSE103111/logs/*.log | grep -i "disk\|space"`
- [ ] Temp directory was created during processing
- [ ] Temp directory cleaned up after completion: `ls scratch/docker_disk_test/GSE103111/tmp/` (should not exist)
- [ ] Final outputs present: `ls scratch/docker_disk_test/GSE103111/results/`

**Expected Duration**: 15-20 minutes (dataset processing time)

**Acceptance Criteria**:
- [ ] All verification items pass
- [ ] Dataset processes to completion
- [ ] No disk-related errors

---

### ✅ PHASE3-2: Test Case 2 - Concurrent Datasets

**Priority**: High
**Effort**: M (15 min)
**Dependencies**: PHASE2-6
**Skills**: Testing
**Assigned**: [Unassigned]

**Description**:
Test parallel processing of 2 small datasets to verify temp directories don't interfere with each other.

**Setup**:
```bash
# Create CSV with 2 small datasets
echo "Accession,DatasetSizeGB" > test_datasets/two_small.csv
echo "GSE12345,5" >> test_datasets/two_small.csv
echo "GSE67890,5" >> test_datasets/two_small.csv
```

**Test Command**:
```bash
uorca run local \
  --input test_datasets/two_small.csv \
  --output scratch/parallel_test/ \
  --max_workers 2
```

**Verification Checklist**:
- [ ] Both datasets start processing
- [ ] Each has its own temp directory during processing
- [ ] No file conflicts between datasets
- [ ] Both datasets complete successfully
- [ ] Both temp directories cleaned up: `ls scratch/parallel_test/GSE*/tmp/` (should not exist for either)
- [ ] Both have final outputs

**Acceptance Criteria**:
- [ ] Parallel processing works correctly
- [ ] Temp directories isolated per dataset
- [ ] All cleanup successful

---

### ✅ PHASE3-3: Test Case 3 - SLURM Unchanged (If HPC Available)

**Priority**: High
**Effort**: M (15 min) - Skip if no HPC access
**Dependencies**: None (SLURM code unchanged)
**Skills**: HPC, SLURM
**Assigned**: [Unassigned]

**Description**:
Verify SLURM execution is completely unaffected by local execution changes.

**Test Command** (on HPC cluster):
```bash
# On HPC with SLURM
uorca run slurm \
  --input test_datasets/small_dataset.csv \
  --output scratch/slurm_test/ \
  --config slurm_config.yaml
```

**Verification Checklist**:
- [ ] Job submits to SLURM successfully
- [ ] Apptainer execution works (no errors)
- [ ] Temp directory handling unchanged (SLURM uses its own temp dir logic)
- [ ] No unexpected errors in SLURM logs
- [ ] Output matches expected format

**Acceptance Criteria**:
- [ ] SLURM execution identical to before changes
- [ ] No regressions

**Note**: If HPC access not available, review `uorca/batch/slurm.py` to confirm no imports or dependencies on modified `local.py` code.

---

### ✅ PHASE3-4: Test Case 4 - Cleanup on Failure

**Priority**: High
**Effort**: S (10 min)
**Dependencies**: PHASE2-6
**Skills**: Testing
**Assigned**: [Unassigned]

**Description**:
Verify temp directory cleanup works even when pipeline fails.

**Setup - Create Invalid Dataset**:
```bash
echo "Accession,DatasetSizeGB" > test_datasets/invalid.csv
echo "INVALID_ACCESSION,5" >> test_datasets/invalid.csv
```

**Test Command**:
```bash
uorca run local \
  --input test_datasets/invalid.csv \
  --output scratch/failure_test/ \
  --max_workers 1
```

**Verification Checklist**:
- [ ] Pipeline fails as expected (invalid accession)
- [ ] Temp directory was created: Look for creation in logs
- [ ] Temp directory cleaned up despite failure: `ls scratch/failure_test/INVALID_ACCESSION/tmp/` (should not exist)
- [ ] No temp directory leakage

**Acceptance Criteria**:
- [ ] Cleanup works correctly on failure
- [ ] No leaked temp directories

---

## 📚 Documentation Updates (Medium Priority)

### ✅ DOC-1: Update README.md

**Priority**: Medium
**Effort**: S (10 min)
**Dependencies**: PHASE3-1 (after successful testing)
**Skills**: Documentation, Markdown
**Assigned**: [Unassigned]

**Description**:
Update README.md to document the new behavior and remove old Docker disk allocation warnings.

**File**: `README.md`
**Section to Update**: "Running Locally" or "Disk Space Requirements"

**Content to Add**:
```markdown
### Disk Space Requirements

UORCA automatically handles temporary files by writing them directly to your system's disk (not Docker's virtual disk). This means:

- **Local execution**: Requires free disk space on your machine (not Docker's allocation)
- **Storage calculation**: Based on your actual available disk space
- **No Docker configuration needed**: Default Docker Desktop settings (50GB) work fine

For large datasets, ensure you have adequate free space on your system:
- Small datasets (< 10 samples): ~50GB free space
- Medium datasets (10-50 samples): ~200GB free space
- Large datasets (> 50 samples): ~500GB+ free space

**Note**: Temporary files are automatically cleaned up after processing.
```

**Acceptance Criteria**:
- [ ] README.md updated with new disk space information
- [ ] Old warnings about Docker disk allocation removed
- [ ] Clear guidance for users

**Files**:
- `README.md`

---

### ✅ DOC-2: Update CLAUDE.md Known Issues

**Priority**: Medium
**Effort**: S (5 min)
**Dependencies**: PHASE3-1
**Skills**: Documentation
**Assigned**: [Unassigned]

**Description**:
Update CLAUDE.md to mark Docker disk limit issue as resolved.

**File**: `CLAUDE.md`
**Section to Update**: "Gotchas & Known Issues"

**Find Issue 3** (or create if not exists):
```markdown
### ~~Issue 3: Docker Desktop Disk Limits~~ ✅ RESOLVED

**Status**: Fixed as of 2025-10-02

**Previous Problem**: Docker Desktop's virtual disk allocation caused "disk-limit exceeded" errors when processing large RNA-seq datasets locally. fasterq-dump wrote temp files to Docker's virtual disk (50GB limit) rather than the host filesystem.

**Solution**: Temporary files now write to host filesystem via bind-mounted directories (`{output_dir}/{accession}/tmp/`), bypassing Docker's virtual disk entirely. No Docker configuration required. Temp directories are automatically created and cleaned up.

**Impact**:
- ✅ Works with default Docker Desktop settings (50GB)
- ✅ Can process datasets of any size (limited by host disk only)
- ✅ Matches SLURM behavior (both use bind-mounted temps)
```

**Acceptance Criteria**:
- [ ] Issue marked as resolved
- [ ] Solution documented
- [ ] Date of fix recorded

**Files**:
- `CLAUDE.md`

---

## 🔍 Regression Verification (Medium Priority)

### ✅ REGRESS-1: Verify SLURM Code Unchanged

**Priority**: Medium
**Effort**: XS (5 min)
**Dependencies**: PHASE2-6
**Skills**: Code Review
**Assigned**: [Unassigned]

**Description**:
Verify that SLURM-related code has not been modified at all.

**Verification Steps**:
```bash
# Check git diff for SLURM files
git diff origin/master -- uorca/batch/slurm.py
git diff origin/master -- main_workflow/run_helpers/run_single_dataset.sbatch.j2

# Should show no changes
```

**Files to Check**:
- [ ] `uorca/batch/slurm.py` - Unchanged
- [ ] `main_workflow/run_helpers/run_single_dataset.sbatch.j2` - Unchanged
- [ ] `uorca/batch/base.py` - Unchanged (if verified)

**Acceptance Criteria**:
- [ ] SLURM files show zero diff
- [ ] Only `uorca/batch/local.py` modified

---

### ✅ REGRESS-2: Local Functionality Regression Tests

**Priority**: Medium
**Effort**: M (15 min)
**Dependencies**: PHASE3-1
**Skills**: Testing
**Assigned**: [Unassigned]

**Description**:
Verify existing local functionality still works (no regressions).

**Test Cases**:

**Test 1: Small Dataset**
```bash
# Should work exactly as before
uorca run local \
  --input test_datasets/small.csv \
  --output scratch/regression_small/ \
  --max_workers 1
```

**Test 2: max_workers=1**
```bash
uorca run local \
  --input test_datasets/medium.csv \
  --output scratch/regression_workers1/ \
  --max_workers 1
```

**Test 3: max_workers=3**
```bash
uorca run local \
  --input test_datasets/medium.csv \
  --output scratch/regression_workers3/ \
  --max_workers 3
```

**Verification**:
- [ ] All test cases complete successfully
- [ ] Output format unchanged
- [ ] No unexpected errors
- [ ] Temp directories cleaned up in all cases

**Acceptance Criteria**:
- [ ] No regressions in existing functionality
- [ ] All test cases pass

---

## 📈 Milestones

### Milestone 1: Phase 1 Complete
**Target**: +30 minutes
**Tasks**: PHASE1-1 through PHASE1-7
**Success**: Temp directories created and bind-mounted, code passes type checking

### Milestone 2: Phase 2 Complete
**Target**: +1 hour
**Tasks**: PHASE2-1 through PHASE2-6
**Success**: Cleanup logic working, temp directories automatically removed

### Milestone 3: All Testing Complete
**Target**: +1.5 hours
**Tasks**: PHASE3-1 through PHASE3-4, DOC-1, DOC-2
**Success**: GSE103111 processes successfully with 50GB Docker, documentation updated

### Milestone 4: Production Ready
**Target**: +2 hours
**Tasks**: All tasks complete
**Success**: All tests pass, code committed, ready for merge

---

## 🐛 Known Issues / Blockers

- None currently identified
- If testing reveals issues, add here

---

## ⚠️ Risks & Mitigation

**Risk 1: Try/Finally Refactoring Breaks Function Logic**
- **Likelihood**: Medium
- **Impact**: High (would break local execution)
- **Mitigation**:
  - Create checkpoints before each phase
  - Test incrementally after each change
  - Use `/rewind` if issues occur

**Risk 2: Cleanup Doesn't Run on All Error Paths**
- **Likelihood**: Low
- **Impact**: Low (leaked temp directories, manual cleanup needed)
- **Mitigation**:
  - Test with intentional failures (Test Case 4)
  - Verify temp directories cleaned in all cases

**Risk 3: Permission Issues with Temp Directory Cleanup**
- **Likelihood**: Very Low
- **Impact**: Low (temp directories left behind)
- **Mitigation**:
  - Container runs as root, so cleanup should always work
  - If issues occur, document workaround in README

---

## 📚 References

- **Implementation Plan**: `docs/plans/bind_mount_temp_directories_for_docker.md`
- **Original Failure Log**: `scratch/local_analysis_test/GSE103111/logs/20251002_013314.log`
- **SLURM Template Reference**: `main_workflow/run_helpers/run_single_dataset.sbatch.j2:14-40`
- **Docker Bind Mounts**: https://docs.docker.com/storage/bind-mounts/

---

## 🔧 Development Setup

No special setup required beyond existing UORCA environment:

```bash
# Verify environment
which uorca
python --version  # Should be 3.13+
docker --version
pyright --version

# Verify Docker Desktop running
docker ps

# Check Docker disk allocation (optional)
docker system df
```

---

## 💡 Tips for Implementation

1. **Work incrementally**: Complete each phase fully before moving to next
2. **Create checkpoints**: Use git commits at each phase boundary
3. **Test after each change**: Run `pyright` after every code modification
4. **Read carefully**: The try/finally refactoring (PHASE2-3) requires careful attention to indentation
5. **Don't skip testing**: Test Case 1 (GSE103111) is critical - it reproduces the original bug
6. **Use /rewind**: If something breaks, rewind to previous checkpoint and try again

---

## 📊 Time Tracking

**Phase 1 (Bind Mount)**:
- Estimated: 30-40 minutes
- Actual: [To be filled]

**Phase 2 (Cleanup)**:
- Estimated: 30-40 minutes
- Actual: [To be filled]

**Phase 3 (Testing)**:
- Estimated: 30-45 minutes
- Actual: [To be filled]

**Documentation**:
- Estimated: 15 minutes
- Actual: [To be filled]

**Total**:
- Estimated: 1.5-2.5 hours
- Actual: [To be filled]

---

*This TODO list is a living document. Check off tasks as completed and update estimates if reality differs.*
