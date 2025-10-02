# TODO List: No-CLI UORCA Implementation

**Generated from**: `docs/plans/no_cli_uorca_implementation.md`
**Generated on**: 2025-10-01
**Total Tasks**: 47
**Estimated Effort**: 5-6 weeks (25-30 working days)
**Target Completion**: 2025-11-12 (6 weeks from now)

---

## ⭐ Key Design Principle: GUI as Thin Wrapper

**CRITICAL**: The Streamlit interface should call existing CLI commands programmatically, NOT reimplement logic.

**Implementation Pattern**:
```python
# ✅ CORRECT: Call existing CLI wrapper
from uorca.identify import main as identify_main
sys.argv = ['identify', '-q', query, '-o', output_dir, ...]
identify_main()

# ❌ WRONG: Reimport and call internal functions
from dataset_identification.DatasetIdentification import extract_terms, perform_search
# ... duplicating CLI logic
```

**Benefits**: No duplication, single source of truth, easier maintenance, CLI/GUI stay in sync

---

## 📋 Quick Summary

- **Critical tasks**: 15 (must do for MVP)
- **High priority**: 18 (core functionality)
- **Medium priority**: 10 (enhancements)
- **Low priority**: 4 (nice-to-haves)

---

## 🎯 Project Goals

Make UORCA accessible to wet lab researchers by eliminating command-line requirements and providing a fully-featured desktop application with:
1. GUI-based dataset identification (calls `uorca identify` in background)
2. GUI-based pipeline execution (calls `uorca run local` in background)
3. Real-time progress tracking and monitoring
4. Desktop launchers (no terminal required)
5. Persistent task history across app restarts

---

## ✅ Success Criteria

- [ ] Users can launch UORCA with double-click (no terminal)
- [x] Complete workflow (identify → run → explore) in GUI
- [x] Real-time progress bars for long operations
- [x] Task history persists across restarts
- [x] All existing CLI commands still work
- [ ] Cross-platform support (Windows, macOS, Linux)
- [ ] No crashes during multi-hour pipeline runs

**Current Status (2025-10-01)**: Phase 1B (Dataset Identification) complete with fixes! Critical lessons learned documented for Phase 1C.

---

## 📚 Lessons Learned from Phase 1B (Dataset Identification)

**Date**: 2025-10-01
**Status**: ✅ Completed with comprehensive fixes

### Critical Issues Resolved

1. **✅ BrokenPipeError from tqdm** - Fixed with thread detection and logging
2. **✅ Streamlit caching staleness** - Database-first status queries
3. **✅ NCBI API key configuration** - Explicit reconfiguration + UI verification
4. **✅ Log rotation** - Timestamped files with automatic cleanup
5. **✅ Duplicate imports** - Code cleanup

**Full Details**: See `docs/lessons_learned/dataset_identification_gui_implementation.md`

### Must-Apply Patterns for Phase 1C (Pipeline Execution)

**⚠️ CRITICAL**: The following patterns MUST be implemented in pipeline tab to avoid the same issues:

#### Pattern 1: Thread-Safe tqdm Handling
```python
# ALWAYS detect thread context and disable tqdm
is_main_thread = threading.current_thread() is threading.main_thread()
for item in tqdm(items, disable=not is_main_thread):
    ...

# REPLACE all tqdm.write() with logging
logging.info("message")  # Not tqdm.write("message")
```

#### Pattern 2: Stdout/Stderr Redirection
```python
# In background task wrapper
original_stdout, original_stderr = sys.stdout, sys.stderr
sys.stdout, sys.stderr = io.StringIO(), io.StringIO()
try:
    pipeline_main()
finally:
    sys.stdout, sys.stderr = original_stdout, original_stderr
```

#### Pattern 3: Database-First Status Checks
```python
# TaskManager.get_task_status() already fixed to query DB first
# Just use it - don't cache status in session_state
status = task_manager.get_task_status(task_id)  # Fresh from DB
```

#### Pattern 4: API Key Pre-Flight Checks
```python
# Add to pipeline submission form
errors = []
if not os.getenv("ENTREZ_EMAIL"):
    errors.append("ENTREZ_EMAIL not set")
if not os.getenv("ENTREZ_API_KEY"):
    st.warning("⚠️ ENTREZ_API_KEY not set - slower rate limits")
# ... check Docker, Kallisto indices, etc.
```

#### Pattern 5: Organized Logging with Rotation
```python
# Create logs/pipeline_logs/ directory
# Use pattern: pipeline_execution_{timestamp}.log
# Keep only 5 most recent
```

### Testing Checklist for Phase 1C

Before marking pipeline tab as complete:

- [ ] Long-running test (3 datasets, 1-2 hours) - must complete successfully
- [ ] Cancellation test - verify graceful Docker container cleanup
- [ ] Restart test - app killed mid-pipeline, restart, task marked failed
- [ ] Error handling - invalid CSV shows clear error
- [ ] Progress accuracy - bar matches actual dataset completion
- [ ] Log rotation - 6 pipelines run, only 5 logs remain
- [ ] Cross-platform - launchers work on Windows/macOS/Linux
- [ ] No orphaned Docker containers after any operation

### Anticipated Challenges for Phase 1C

1. **LocalBatchProcessor likely has tqdm** - will need same thread detection fixes
2. **Docker containers need env vars** - must verify ENTREZ_* passed correctly
3. **Longer execution times** - more opportunities for staleness bugs
4. **Container cleanup** - cancellation must kill Docker containers
5. **Per-dataset logs** - consider separate log files in addition to master log

---

## 🚫 Out of Scope

- SLURM support in GUI (remains CLI-only)
- Multi-user/cloud deployment (future phase)
- Mobile/responsive design (desktop-focused)
- Real-time log streaming (status updates only)
- Advanced monitoring graphs (basic metrics only)

---

## 📊 Progress Tracking

**Overall Progress**: 0/47 tasks (0%)

| Phase | Tasks | Estimated Effort | Status |
|-------|-------|------------------|--------|
| Phase 1A: Background Task System | 0/8 | 1 week | Not Started |
| Phase 1B: Dataset Identification Tab | 0/10 | 1 week | Not Started |
| Phase 1C: Pipeline Execution Tab | 0/12 | 2 weeks | Not Started |
| Phase 2A: Batch Launchers | 0/7 | 3 days | Not Started |
| Phase 2B: Native Desktop App | 0/6 | 4 days (Optional) | Not Started |
| Phase 3: Polish & Testing | 0/4 | 1 week | Not Started |

---

## 🔴 Phase 1A: Background Task System (Week 1)

### ✅ CRIT-1A.1: Create TaskManager Core Class

**Priority**: Critical
**Effort**: M (4-6 hours)
**Dependencies**: None
**Skills**: Backend, Threading, SQLite
**Assigned**: [Unassigned]

**Description**:
Create the core `TaskManager` class with threading and SQLite persistence for background task execution.

**Implementation**:
- [ ] Create `main_workflow/reporting/core/task_manager.py`
- [ ] Implement `TaskStatus` enum (PENDING, RUNNING, COMPLETED, FAILED, CANCELLED)
- [ ] Implement singleton `TaskManager` class with thread-safe initialization
- [ ] Set up SQLite database schema in `~/.uorca/tasks.db`
- [ ] Implement `_setup_database()` method with tasks table
- [ ] Implement `_load_existing_tasks()` to mark interrupted tasks as failed

**Acceptance Criteria**:
- [ ] TaskManager initializes as singleton
- [ ] SQLite database created at `~/.uorca/tasks.db`
- [ ] Tasks table has correct schema (task_id, task_type, status, timestamps, progress, etc.)
- [ ] Interrupted tasks marked as failed on restart
- [ ] Thread-safe concurrent access

**Files**:
- `main_workflow/reporting/core/task_manager.py` (new, lines 1-180)

**Testing**:
```python
# Test initialization
manager = TaskManager()
assert Path.home() / ".uorca" / "tasks.db" exists()
```

---

### ✅ CRIT-1A.2: Implement Task Submission

**Priority**: Critical
**Effort**: M (4-6 hours)
**Dependencies**: CRIT-1A.1
**Skills**: Backend, Threading
**Assigned**: [Unassigned]

**Description**:
Implement task submission with background thread execution and callback support.

**Implementation**:
- [ ] Implement `submit_task()` method with parameters: task_id, task_type, task_func, parameters, callbacks
- [ ] Create task wrapper function for thread execution
- [ ] Store task in SQLite database with PENDING status
- [ ] Initialize in-memory task state dictionary
- [ ] Create and start daemon thread for task execution
- [ ] Handle task exceptions with traceback capture

**Acceptance Criteria**:
- [ ] Tasks can be submitted with unique task_id
- [ ] Tasks execute in separate threads without blocking
- [ ] Task parameters passed correctly to task_func
- [ ] Callbacks (on_progress, on_complete, on_error) invoked appropriately
- [ ] Thread is daemon (doesn't prevent app exit)

**Files**:
- `main_workflow/reporting/core/task_manager.py` (lines 181-253)

**Testing**:
```python
def test_task(value, progress_callback=None):
    return value * 2

task_id = manager.submit_task("test_1", "test", test_task, {"value": 5})
time.sleep(0.2)
status = manager.get_task_status(task_id)
assert status["result"] == 10
```

---

### ✅ CRIT-1A.3: Implement Task Status Tracking

**Priority**: Critical
**Effort**: S (2-4 hours)
**Dependencies**: CRIT-1A.2
**Skills**: Backend, SQLite
**Assigned**: [Unassigned]

**Description**:
Implement real-time status updates for running tasks with database persistence.

**Implementation**:
- [ ] Implement `_update_task_status()` with in-memory and DB updates
- [ ] Handle status transitions (PENDING → RUNNING → COMPLETED/FAILED)
- [ ] Update progress (0.0 to 1.0) and progress_message
- [ ] Store result_path and error_message
- [ ] Record timestamps for started_at and completed_at
- [ ] Implement `_on_progress()` callback wrapper

**Acceptance Criteria**:
- [ ] Status updates reflected in both memory and SQLite
- [ ] Progress values between 0.0 and 1.0
- [ ] Timestamps recorded accurately
- [ ] Error messages include full tracebacks
- [ ] Results persisted for completed tasks

**Files**:
- `main_workflow/reporting/core/task_manager.py` (lines 254-322)

**Testing**:
```python
# Test progress tracking
progress_updates = []
def track(p, m):
    progress_updates.append((p, m))

manager.submit_task("test", "test", slow_task, {}, on_progress=track)
time.sleep(0.5)
assert len(progress_updates) > 0
```

---

### ✅ CRIT-1A.4: Implement Task Query Methods

**Priority**: Critical
**Effort**: S (2-3 hours)
**Dependencies**: CRIT-1A.3
**Skills**: Backend, SQLite
**Assigned**: [Unassigned]

**Description**:
Implement methods to query task status and list historical tasks.

**Implementation**:
- [ ] Implement `get_task_status(task_id)` - returns current task state
- [ ] Check in-memory cache first, fall back to database
- [ ] Implement `list_tasks(task_type, limit)` - list recent tasks
- [ ] Support filtering by task_type (e.g., "identify", "pipeline")
- [ ] Limit results (default 50)
- [ ] Return tasks ordered by created_at DESC

**Acceptance Criteria**:
- [ ] get_task_status returns dict with status, progress, message, result, error
- [ ] Returns None for non-existent task_id
- [ ] list_tasks returns list of task dicts with all metadata
- [ ] Filtering by task_type works correctly
- [ ] Limit parameter respected

**Files**:
- `main_workflow/reporting/core/task_manager.py` (lines 323-385)

**Testing**:
```python
# Test task listing
manager.submit_task("task1", "identify", func1, {})
manager.submit_task("task2", "pipeline", func2, {})
identify_tasks = manager.list_tasks(task_type="identify")
assert len(identify_tasks) == 1
assert identify_tasks[0]["task_id"] == "task1"
```

---

### ✅ CRIT-1A.5: Implement Task Cancellation

**Priority**: High
**Effort**: S (2-3 hours)
**Dependencies**: CRIT-1A.3
**Skills**: Backend, Threading
**Assigned**: [Unassigned]

**Description**:
Implement best-effort task cancellation (Python threading limitations apply).

**Implementation**:
- [ ] Implement `cancel_task(task_id)` method
- [ ] Check if task is currently RUNNING
- [ ] Mark task as CANCELLED in database
- [ ] Update in-memory state
- [ ] Return True if cancelled, False otherwise
- [ ] Document limitation: Python threads can't be force-killed

**Acceptance Criteria**:
- [ ] cancel_task returns True for running tasks
- [ ] cancel_task returns False for completed/non-existent tasks
- [ ] Task status updated to CANCELLED in database
- [ ] Task function should check progress_callback return for cancellation signal

**Files**:
- `main_workflow/reporting/core/task_manager.py` (lines 386-399)

**Testing**:
```python
task_id = manager.submit_task("long_task", "test", long_running, {})
time.sleep(0.1)
assert manager.cancel_task(task_id) == True
status = manager.get_task_status(task_id)
assert status["status"] == TaskStatus.CANCELLED
```

---

### ✅ CRIT-1A.6: Update Core Module Exports

**Priority**: Critical
**Effort**: XS (15 minutes)
**Dependencies**: CRIT-1A.1
**Skills**: Backend
**Assigned**: [Unassigned]

**Description**:
Update `main_workflow/reporting/core/__init__.py` to export TaskManager classes.

**Implementation**:
- [ ] Import TaskManager from task_manager module
- [ ] Import TaskStatus enum
- [ ] Add to __all__ list for public API

**Acceptance Criteria**:
- [ ] Can import: `from main_workflow.reporting.core import TaskManager, TaskStatus`
- [ ] No import errors

**Files**:
- `main_workflow/reporting/core/__init__.py` (update lines 1-10)

**Testing**:
```python
from main_workflow.reporting.core import TaskManager, TaskStatus
assert TaskManager is not None
assert TaskStatus.PENDING is not None
```

---

### ✅ CRIT-1A.7: Unit Tests for TaskManager

**Priority**: Critical
**Effort**: M (4-6 hours)
**Dependencies**: CRIT-1A.1 through CRIT-1A.6
**Skills**: Testing, pytest
**Assigned**: [Unassigned]

**Description**:
Create comprehensive unit tests for TaskManager covering all functionality.

**Implementation**:
- [ ] Create `tests/unit/reporting/test_task_manager.py`
- [ ] Test basic task submission and completion
- [ ] Test progress tracking callbacks
- [ ] Test error handling and failed tasks
- [ ] Test concurrent tasks (5+ simultaneous)
- [ ] Test SQLite persistence across manager instances
- [ ] Test task cancellation
- [ ] Use temporary database for testing

**Acceptance Criteria**:
- [ ] All tests pass: `pytest tests/unit/reporting/test_task_manager.py -v`
- [ ] Coverage >90% for task_manager.py
- [ ] Tests use temporary database (cleanup after)
- [ ] Tests complete in <5 seconds

**Files**:
- `tests/unit/reporting/test_task_manager.py` (new, ~200 lines)

**Testing**:
```bash
pytest tests/unit/reporting/test_task_manager.py -v
pytest tests/unit/reporting/test_task_manager.py --cov=main_workflow/reporting/core/task_manager
```

---

### ✅ CRIT-1A.8: Manual Integration Testing

**Priority**: High
**Effort**: S (2-3 hours)
**Dependencies**: CRIT-1A.7
**Skills**: Testing
**Assigned**: [Unassigned]

**Description**:
Manually test TaskManager in realistic scenarios to validate behavior.

**Test Scenarios**:
- [ ] Submit test task and verify progress updates appear
- [ ] Submit multiple tasks concurrently (10+)
- [ ] Restart Python process mid-task
- [ ] Verify interrupted tasks marked as failed in database
- [ ] Test task history retrieval
- [ ] Test with large result objects
- [ ] Test with long-running tasks (1+ minute)

**Acceptance Criteria**:
- [ ] No race conditions or deadlocks
- [ ] Database queries remain fast (<50ms)
- [ ] Memory usage stable with many tasks
- [ ] Crashed tasks properly recorded

**Files**:
- Manual testing script (optional: `tests/manual/test_task_manager_integration.py`)

---

## 🔴 Phase 1B: Dataset Identification Tab (Week 2)

### ✅ CRIT-1B.1: Create Dataset Identification Tab Structure

**Priority**: Critical
**Effort**: M (3-4 hours)
**Dependencies**: CRIT-1A (TaskManager complete)
**Skills**: Frontend, Streamlit
**Assigned**: [Unassigned]

**Description**:
Create the basic structure for the dataset identification tab with form UI.

**Implementation**:
- [ ] Create `main_workflow/reporting/tabs/dataset_identification.py`
- [ ] Implement `show_dataset_identification_tab()` main function
- [ ] Create header with description
- [ ] Initialize TaskManager
- [ ] Check for running tasks in session_state
- [ ] Create main input form with st.form()

**Acceptance Criteria**:
- [ ] Tab displays with clear header
- [ ] Form renders with all input fields
- [ ] TaskManager initializes without errors
- [ ] Session state properly checked

**Files**:
- `main_workflow/reporting/tabs/dataset_identification.py` (new, lines 1-100)

**Testing**:
```bash
# Manual: Open UORCA Explorer and verify new tab appears
uv run uorca explore
```

---

### ✅ CRIT-1B.2: Implement Identification Form Inputs

**Priority**: Critical
**Effort**: M (3-4 hours)
**Dependencies**: CRIT-1B.1
**Skills**: Frontend, Streamlit
**Assigned**: [Unassigned]

**Description**:
Implement all form inputs for dataset identification parameters.

**Implementation**:
- [ ] Research question text_area (placeholder example)
- [ ] max_results number_input (10-500, default 100)
- [ ] num_datasets number_input (1-50, default 10)
- [ ] model selectbox (gpt-4o-mini, gpt-4o, gpt-4-turbo)
- [ ] samples_required number_input (2-10, default 3)
- [ ] output_dir text_input (default: ~/UORCA_Identification/timestamp)
- [ ] Submit button (primary, full width)
- [ ] Form layout with columns for organization

**Acceptance Criteria**:
- [ ] All inputs have clear labels and help text
- [ ] Default values are sensible
- [ ] Output directory defaults to timestamped folder
- [ ] Form layout is clean and organized

**Files**:
- `main_workflow/reporting/tabs/dataset_identification.py` (lines 100-200)

**Testing**:
- [ ] All inputs accept valid values
- [ ] Help tooltips are informative
- [ ] Form submits successfully

---

### ✅ CRIT-1B.3: Implement Form Validation

**Priority**: Critical
**Effort**: S (2-3 hours)
**Dependencies**: CRIT-1B.2
**Skills**: Frontend, Python
**Assigned**: [Unassigned]

**Description**:
Add validation for form inputs before task submission.

**Implementation**:
- [ ] Validate research_query is not empty
- [ ] Check OPENAI_API_KEY environment variable
- [ ] Check ENTREZ_EMAIL environment variable
- [ ] Validate output_dir is writable
- [ ] Create output directory if doesn't exist
- [ ] Display clear error messages with st.error()

**Acceptance Criteria**:
- [ ] Empty query shows error: "Please enter a research question"
- [ ] Missing API key shows error with setup instructions
- [ ] Missing email shows error with setup instructions
- [ ] Invalid output path shows clear error
- [ ] All errors prevent task submission

**Files**:
- `main_workflow/reporting/tabs/dataset_identification.py` (lines 200-250)

**Testing**:
```python
# Test cases:
# 1. Submit with empty query → error
# 2. Submit without OPENAI_API_KEY → error
# 3. Submit without ENTREZ_EMAIL → error
# 4. Valid submission → no errors
```

---

### ✅ CRIT-1B.4: Implement Task Submission Logic

**Priority**: Critical
**Effort**: M (3-4 hours)
**Dependencies**: CRIT-1B.3
**Skills**: Backend, Frontend
**Assigned**: [Unassigned]

**Description**:
Submit identification task to TaskManager when form is valid.

**Implementation**:
- [ ] Generate unique task_id (identify_YYYYMMDD_HHMMSS)
- [ ] Call task_manager.submit_task() with correct parameters
- [ ] Pass research_query, max_results, num_datasets, model, samples_required, output_dir
- [ ] Pass _run_identification as task_func
- [ ] Store task_id in st.session_state.current_identify_task
- [ ] Trigger st.rerun() to show progress view

**Acceptance Criteria**:
- [ ] Task submits successfully
- [ ] Task_id stored in session state
- [ ] Page reruns to show progress view
- [ ] Task appears in TaskManager database

**Files**:
- `main_workflow/reporting/tabs/dataset_identification.py` (lines 250-300)

**Testing**:
```python
# Submit valid form
# Verify task appears in: manager.list_tasks(task_type="identify")
```

---

### ✅ CRIT-1B.5: Implement Running Task Progress View

**Priority**: Critical
**Effort**: L (6-8 hours)
**Dependencies**: CRIT-1B.4
**Skills**: Frontend, Streamlit
**Assigned**: [Unassigned]

**Description**:
Display real-time progress for running identification tasks.

**Implementation**:
- [ ] Implement `_show_running_task()` function
- [ ] Display progress bar with st.progress()
- [ ] Show progress message with st.info()
- [ ] Display Cancel button
- [ ] Auto-refresh while status is RUNNING
- [ ] Show success view when COMPLETED
- [ ] Show error view when FAILED
- [ ] Display results preview (CSV dataframe)
- [ ] Add download button for results CSV

**Acceptance Criteria**:
- [ ] Progress bar updates in real-time
- [ ] Progress messages are informative
- [ ] Page auto-refreshes during execution
- [ ] Completed tasks show results preview
- [ ] Failed tasks show error message
- [ ] CSV download works

**Files**:
- `main_workflow/reporting/tabs/dataset_identification.py` (lines 300-400)

**Testing**:
```python
# Manual testing:
# 1. Submit task
# 2. Verify progress bar moves
# 3. Verify messages update
# 4. Verify completion shows results
```

---

### ✅ CRIT-1B.6: Implement Identification Backend Function

**Priority**: Critical
**Effort**: L (8-10 hours)
**Dependencies**: CRIT-1B.4
**Skills**: Backend, AI, Bioinformatics
**Assigned**: [Unassigned]

**Description**:
Implement `_run_identification()` function that executes the identification workflow.

**Implementation**:
- [ ] Import DatasetIdentification module functions
- [ ] Step 1: extract_terms() from research query (progress: 1/8)
- [ ] Step 2: perform_search() for each term (progress: 2/8)
- [ ] Step 3: get_basic_dataset_info() (progress: 3/8)
- [ ] Step 4: validate_datasets() with SRA metadata (progress: 4/8)
- [ ] Step 5: embed_datasets() (progress: 5/8)
- [ ] Step 6: cluster_datasets() (progress: 6/8)
- [ ] Step 7: select_representative_datasets() (progress: 7/8)
- [ ] Step 8: repeated_relevance() AI assessment (progress: 8/8)
- [ ] Save results to CSV
- [ ] Call progress_callback at each step
- [ ] Handle errors with proper exception handling

**Acceptance Criteria**:
- [ ] All 8 steps execute in correct order
- [ ] Progress callback invoked with 0.125, 0.25, ..., 1.0
- [ ] Results saved to Dataset_identification_result.csv
- [ ] Selected datasets saved to selected_datasets.csv
- [ ] Errors caught and returned with traceback
- [ ] Function returns output directory path

**Files**:
- `main_workflow/reporting/tabs/dataset_identification.py` (lines 400-600)

**Testing**:
```python
# Integration test with mocked API calls
result = _run_identification(
    research_query="test query",
    max_results=10,
    num_datasets=3,
    model="gpt-4o-mini",
    samples_required=3,
    output_dir="/tmp/test"
)
assert Path(result).exists()
assert (Path(result) / "Dataset_identification_result.csv").exists()
```

---

### ✅ CRIT-1B.7: Implement Task History View

**Priority**: High
**Effort**: M (3-4 hours)
**Dependencies**: CRIT-1B.5
**Skills**: Frontend, Streamlit
**Assigned**: [Unassigned]

**Description**:
Display recent identification tasks with expandable details.

**Implementation**:
- [ ] Implement `_show_task_history()` function
- [ ] Call task_manager.list_tasks(task_type="identify", limit=10)
- [ ] Display each task in st.expander()
- [ ] Show status, timestamps, parameters
- [ ] Show result path for completed tasks
- [ ] Show error message for failed tasks
- [ ] Add "Start New Identification" button

**Acceptance Criteria**:
- [ ] Recent tasks listed newest first
- [ ] Each task shows status badge
- [ ] Timestamps formatted clearly
- [ ] Parameters preview (truncated if long)
- [ ] Failed tasks show error details

**Files**:
- `main_workflow/reporting/tabs/dataset_identification.py` (lines 600-700)

**Testing**:
```python
# Submit 3 tasks
# Verify all 3 appear in history
# Verify status badges correct
```

---

### ✅ CRIT-1B.8: Update Tab Module Exports

**Priority**: Critical
**Effort**: XS (15 minutes)
**Dependencies**: CRIT-1B.1
**Skills**: Backend
**Assigned**: [Unassigned]

**Description**:
Update tab module exports to include new dataset identification tab.

**Implementation**:
- [ ] Update `main_workflow/reporting/tabs/__init__.py`
- [ ] Import show_dataset_identification_tab
- [ ] Add to __all__ list

**Acceptance Criteria**:
- [ ] Can import: `from main_workflow.reporting.tabs import show_dataset_identification_tab`
- [ ] No import errors

**Files**:
- `main_workflow/reporting/tabs/__init__.py` (update)

---

### ✅ CRIT-1B.9: Integrate Tab into UORCA Explorer

**Priority**: Critical
**Effort**: S (1-2 hours)
**Dependencies**: CRIT-1B.8
**Skills**: Frontend, Streamlit
**Assigned**: [Unassigned]

**Description**:
Add dataset identification tab to main UORCA Explorer application.

**Implementation**:
- [ ] Update `main_workflow/reporting/uorca_explorer.py`
- [ ] Import show_dataset_identification_tab
- [ ] Add "🔍 Identify Datasets" to tabs list (first position)
- [ ] Call show_dataset_identification_tab() in tab context
- [ ] Shift existing tab indices by 1

**Acceptance Criteria**:
- [ ] New tab appears as first tab
- [ ] Existing tabs still work
- [ ] No import errors
- [ ] Tab renders correctly

**Files**:
- `main_workflow/reporting/uorca_explorer.py` (lines ~40-50)

**Testing**:
```bash
uv run uorca explore
# Verify "🔍 Identify Datasets" tab appears first
```

---

### ✅ CRIT-1B.10: Integration Tests for Identification Tab

**Priority**: High
**Effort**: M (4-6 hours)
**Dependencies**: CRIT-1B.9
**Skills**: Testing, pytest
**Assigned**: [Unassigned]

**Description**:
Create integration tests for complete identification workflow.

**Implementation**:
- [ ] Create `tests/integration/reporting/test_dataset_identification_tab.py`
- [ ] Mock API calls (extract_terms, perform_search, etc.)
- [ ] Test complete workflow with test query
- [ ] Test progress callback invocation
- [ ] Test result file creation
- [ ] Test error handling
- [ ] Test with invalid inputs

**Acceptance Criteria**:
- [ ] All tests pass: `pytest tests/integration/reporting/test_dataset_identification_tab.py -v`
- [ ] Tests use mocked API calls (no real OpenAI/NCBI)
- [ ] Tests verify CSV outputs created
- [ ] Tests complete in <10 seconds

**Files**:
- `tests/integration/reporting/test_dataset_identification_tab.py` (new, ~150 lines)

**Testing**:
```bash
pytest tests/integration/reporting/test_dataset_identification_tab.py -v
```

---

## 🔴 Phase 1C: Pipeline Execution Tab (Weeks 3-4)

### ✅ CRIT-1C.1: Create Pipeline Execution Tab Structure

**Priority**: Critical
**Effort**: M (3-4 hours)
**Dependencies**: CRIT-1A (TaskManager complete)
**Skills**: Frontend, Streamlit
**Assigned**: [Unassigned]

**Description**:
Create basic structure for pipeline execution tab with input selection.

**Implementation**:
- [ ] Create `main_workflow/reporting/tabs/pipeline_execution.py`
- [ ] Implement `show_pipeline_execution_tab()` main function
- [ ] Create header with description
- [ ] Initialize TaskManager
- [ ] Check for running pipeline in session_state
- [ ] Create input source selection (Upload CSV vs Use Identification Results)

**Acceptance Criteria**:
- [ ] Tab displays with clear header
- [ ] Input method radio buttons work
- [ ] TaskManager initializes
- [ ] Session state properly checked

**Files**:
- `main_workflow/reporting/tabs/pipeline_execution.py` (new, lines 1-100)

**Testing**:
```bash
uv run uorca explore
# Verify "🚀 Run Pipeline" tab appears
```

---

### ✅ CRIT-1C.2: Implement CSV Upload Input

**Priority**: Critical
**Effort**: M (3-4 hours)
**Dependencies**: CRIT-1C.1
**Skills**: Frontend, Streamlit
**Assigned**: [Unassigned]

**Description**:
Implement CSV file upload with preview and validation.

**Implementation**:
- [ ] Add st.file_uploader() for CSV
- [ ] Save uploaded file to temp location (~/.uorca/temp/)
- [ ] Display CSV preview with st.dataframe()
- [ ] Validate required columns: geo_accession, SafeStorageGB
- [ ] Show error if columns missing
- [ ] Store file path for pipeline submission

**Acceptance Criteria**:
- [ ] File upload works for .csv files
- [ ] Preview shows first rows
- [ ] Validation catches missing columns
- [ ] Clear error messages

**Files**:
- `main_workflow/reporting/tabs/pipeline_execution.py` (lines 100-200)

**Testing**:
```python
# Upload valid CSV → preview appears
# Upload invalid CSV → error shown
```

---

### ✅ CRIT-1C.3: Implement Identification Results Selection

**Priority**: Critical
**Effort**: M (3-4 hours)
**Dependencies**: CRIT-1C.1, CRIT-1B (Identification complete)
**Skills**: Frontend, Backend
**Assigned**: [Unassigned]

**Description**:
Allow selection from completed identification tasks.

**Implementation**:
- [ ] Query task_manager.list_tasks(task_type="identify")
- [ ] Filter for completed tasks
- [ ] Create selectbox with task_id + timestamp
- [ ] Load selected_datasets.csv from result directory
- [ ] Display preview of selected datasets
- [ ] Show dataset count
- [ ] Handle missing selected_datasets.csv

**Acceptance Criteria**:
- [ ] Dropdown shows completed identification tasks
- [ ] Selecting task loads CSV preview
- [ ] Shows count: "Found N datasets to process"
- [ ] Handles missing files gracefully

**Files**:
- `main_workflow/reporting/tabs/pipeline_execution.py` (lines 200-300)

**Testing**:
```python
# Complete identification task
# Open pipeline tab
# Verify task appears in dropdown
# Select task → CSV preview loads
```

---

### ✅ CRIT-1C.4: Implement Pipeline Configuration Form

**Priority**: Critical
**Effort**: M (3-4 hours)
**Dependencies**: CRIT-1C.2, CRIT-1C.3
**Skills**: Frontend, Streamlit
**Assigned**: [Unassigned]

**Description**:
Create form for pipeline execution parameters.

**Implementation**:
- [ ] output_dir text_input (default: ~/UORCA_Results/timestamp)
- [ ] max_workers number_input (1-16, default 4)
- [ ] max_storage_gb number_input (10-1000, default 100)
- [ ] container_tmpfs_gb number_input (5-100, default 20)
- [ ] Advanced expander with:
  - research_question text_area (optional)
  - docker_image text_input (default: kevingchen/uorca:0.1.0)
- [ ] Submit button (primary, full width)

**Acceptance Criteria**:
- [ ] All inputs have sensible defaults
- [ ] Help text explains each parameter
- [ ] Form layout is organized
- [ ] Advanced options hidden by default

**Files**:
- `main_workflow/reporting/tabs/pipeline_execution.py` (lines 300-400)

**Testing**:
- [ ] Form renders with all inputs
- [ ] Default values are correct
- [ ] Advanced expander toggles

---

### ✅ CRIT-1C.5: Implement Pipeline Validation

**Priority**: Critical
**Effort**: S (2-3 hours)
**Dependencies**: CRIT-1C.4
**Skills**: Backend, Python
**Assigned**: [Unassigned]

**Description**:
Validate pipeline parameters before submission.

**Implementation**:
- [ ] Validate input CSV exists and has required columns
- [ ] Check ENTREZ_EMAIL environment variable
- [ ] Validate output_dir is writable
- [ ] Create output directory
- [ ] Check Docker is running (optional check)
- [ ] Display clear error messages

**Acceptance Criteria**:
- [ ] Invalid CSV shows error with column requirements
- [ ] Missing email shows setup instructions
- [ ] Invalid output path shows error
- [ ] All validation errors prevent submission

**Files**:
- `main_workflow/reporting/tabs/pipeline_execution.py` (lines 400-450)

**Testing**:
```python
# Test cases:
# 1. Invalid CSV → error
# 2. Missing ENTREZ_EMAIL → error
# 3. Invalid output path → error
# 4. Valid inputs → no errors
```

---

### ✅ CRIT-1C.6: Implement Pipeline Task Submission

**Priority**: Critical
**Effort**: M (3-4 hours)
**Dependencies**: CRIT-1C.5
**Skills**: Backend, Frontend
**Assigned**: [Unassigned]

**Description**:
Submit pipeline execution task to TaskManager.

**Implementation**:
- [ ] Generate task_id (pipeline_YYYYMMDD_HHMMSS)
- [ ] Call task_manager.submit_task()
- [ ] Pass input_csv, output_dir, max_workers, max_storage_gb, container_tmpfs_gb, docker_image, research_question
- [ ] Pass _run_pipeline as task_func
- [ ] Store task_id in st.session_state.current_pipeline_task
- [ ] Trigger st.rerun()

**Acceptance Criteria**:
- [ ] Task submits successfully
- [ ] Task_id stored in session state
- [ ] Page reruns to show progress
- [ ] Task appears in database

**Files**:
- `main_workflow/reporting/tabs/pipeline_execution.py` (lines 450-500)

**Testing**:
```python
# Submit pipeline
# Verify task in: manager.list_tasks(task_type="pipeline")
```

---

### ✅ CRIT-1C.7: Implement Running Pipeline Progress View

**Priority**: Critical
**Effort**: L (8-10 hours)
**Dependencies**: CRIT-1C.6
**Skills**: Frontend, Streamlit
**Assigned**: [Unassigned]

**Description**:
Display real-time progress for running pipeline with detailed job status.

**Implementation**:
- [ ] Implement `_show_running_pipeline()` function
- [ ] Display overall progress bar
- [ ] Show progress message
- [ ] Display metrics: Completed, Running, Queued, Failed
- [ ] Show detailed job status table (dataframe)
- [ ] Add Cancel and Refresh buttons
- [ ] Auto-refresh every 5 seconds while RUNNING
- [ ] Show success view when COMPLETED
- [ ] Add "View Results" button to navigate to Results tab
- [ ] Show error view when FAILED

**Acceptance Criteria**:
- [ ] Progress bar updates in real-time
- [ ] Job metrics show current counts
- [ ] Job details table shows individual datasets
- [ ] Auto-refresh works (5 second interval)
- [ ] Cancel button works
- [ ] Navigation to Results tab works

**Files**:
- `main_workflow/reporting/tabs/pipeline_execution.py` (lines 500-700)

**Testing**:
```python
# Manual:
# 1. Submit pipeline with 3 datasets
# 2. Verify progress bar updates
# 3. Verify job metrics update
# 4. Verify job details table shows status
# 5. Test cancel button
```

---

### ✅ CRIT-1C.8: Implement Pipeline Backend Function

**Priority**: Critical
**Effort**: L (8-10 hours)
**Dependencies**: CRIT-1C.6
**Skills**: Backend, Docker, Batch Processing
**Assigned**: [Unassigned]

**Description**:
Implement `_run_pipeline()` function that executes the pipeline workflow.

**Implementation**:
- [ ] Import LocalBatchProcessor from uorca.batch.local
- [ ] Initialize processor
- [ ] Call processor.submit_datasets() with parameters
- [ ] Monitor job status in loop
- [ ] Calculate progress: completed / total_datasets
- [ ] Format progress message with job counts
- [ ] Call progress_callback with updates
- [ ] Poll every 5 seconds
- [ ] Return when status is 'completed' or 'failed'
- [ ] Return dict with output_dir, status, completed, failed, job_details

**Acceptance Criteria**:
- [ ] LocalBatchProcessor initializes
- [ ] Datasets submitted successfully
- [ ] Progress callback invoked regularly
- [ ] Loop exits on completion
- [ ] Returns complete status dict
- [ ] Errors caught and logged

**Files**:
- `main_workflow/reporting/tabs/pipeline_execution.py` (lines 700-900)

**Testing**:
```python
# Integration test with small dataset
result = _run_pipeline(
    input_csv="test.csv",
    output_dir="/tmp/test",
    max_workers=2,
    max_storage_gb=50,
    container_tmpfs_gb=10,
    docker_image="kevingchen/uorca:0.1.0",
    research_question=""
)
assert result["status"] in ["completed", "failed"]
```

---

### ✅ CRIT-1C.9: Add Progress Callback to LocalBatchProcessor

**Priority**: Critical
**Effort**: M (4-6 hours)
**Dependencies**: CRIT-1C.8
**Skills**: Backend, Python
**Assigned**: [Unassigned]

**Description**:
Modify `uorca/batch/local.py` to support progress callbacks.

**Implementation**:
- [ ] Add progress_callback parameter to submit_datasets() method
- [ ] Store callback in instance variable
- [ ] In monitor_jobs() method, call callback with progress updates
- [ ] Calculate progress: len(completed_jobs) / total_jobs
- [ ] Format message with running, completed, queued counts
- [ ] Call callback every monitoring cycle (5 seconds)
- [ ] Handle None callback gracefully

**Acceptance Criteria**:
- [ ] submit_datasets accepts progress_callback parameter
- [ ] Callback invoked during monitoring
- [ ] Progress values between 0.0 and 1.0
- [ ] Messages formatted clearly
- [ ] Works without callback (backwards compatible)

**Files**:
- `uorca/batch/local.py` (lines ~200, ~400)

**Testing**:
```python
progress_updates = []
def track(p, m):
    progress_updates.append((p, m))

processor = LocalBatchProcessor()
processor.submit_datasets(..., progress_callback=track)
# Wait for completion
assert len(progress_updates) > 0
```

---

### ✅ CRIT-1C.10: Implement Pipeline History View

**Priority**: High
**Effort**: M (3-4 hours)
**Dependencies**: CRIT-1C.7
**Skills**: Frontend, Streamlit
**Assigned**: [Unassigned]

**Description**:
Display recent pipeline executions with expandable details.

**Implementation**:
- [ ] Implement `_show_pipeline_history()` function
- [ ] Call task_manager.list_tasks(task_type="pipeline", limit=10)
- [ ] Display each task in st.expander()
- [ ] Show status, timestamps, parameters (workers, storage)
- [ ] Show output directory
- [ ] Add "View Results" button for completed tasks
- [ ] Show error message for failed tasks

**Acceptance Criteria**:
- [ ] Recent pipelines listed newest first
- [ ] Each pipeline shows status badge
- [ ] Parameters displayed clearly
- [ ] "View Results" button navigates to Results tab
- [ ] Failed pipelines show errors

**Files**:
- `main_workflow/reporting/tabs/pipeline_execution.py` (lines 900-1000)

**Testing**:
```python
# Submit 2 pipelines
# Verify both appear in history
# Click "View Results" → navigates to Results tab
```

---

### ✅ CRIT-1C.11: Integrate Pipeline Tab into UORCA Explorer

**Priority**: Critical
**Effort**: S (1-2 hours)
**Dependencies**: CRIT-1C.1
**Skills**: Frontend, Streamlit
**Assigned**: [Unassigned]

**Description**:
Add pipeline execution tab to main UORCA Explorer application.

**Implementation**:
- [ ] Update `main_workflow/reporting/tabs/__init__.py`
- [ ] Import show_pipeline_execution_tab
- [ ] Add to __all__ list
- [ ] Update `main_workflow/reporting/uorca_explorer.py`
- [ ] Add "🚀 Run Pipeline" to tabs list (second position)
- [ ] Call show_pipeline_execution_tab() in tab context

**Acceptance Criteria**:
- [ ] New tab appears as second tab
- [ ] All tabs render correctly
- [ ] No import errors

**Files**:
- `main_workflow/reporting/tabs/__init__.py` (update)
- `main_workflow/reporting/uorca_explorer.py` (update)

**Testing**:
```bash
uv run uorca explore
# Verify tab order: Identify, Run Pipeline, Results, ...
```

---

### ✅ CRIT-1C.12: Integration Tests for Pipeline Tab

**Priority**: High
**Effort**: M (4-6 hours)
**Dependencies**: CRIT-1C.11
**Skills**: Testing, pytest
**Assigned**: [Unassigned]

**Description**:
Create integration tests for pipeline execution workflow.

**Implementation**:
- [ ] Create `tests/integration/reporting/test_pipeline_execution_tab.py`
- [ ] Mock LocalBatchProcessor
- [ ] Test pipeline submission
- [ ] Test progress monitoring
- [ ] Test completion handling
- [ ] Test error handling
- [ ] Test CSV validation

**Acceptance Criteria**:
- [ ] All tests pass: `pytest tests/integration/reporting/test_pipeline_execution_tab.py -v`
- [ ] Tests use mocked batch processor
- [ ] Tests verify status updates
- [ ] Tests complete in <10 seconds

**Files**:
- `tests/integration/reporting/test_pipeline_execution_tab.py` (new, ~200 lines)

**Testing**:
```bash
pytest tests/integration/reporting/test_pipeline_execution_tab.py -v
```

---

## 🟡 Phase 2A: Batch Launcher Scripts (3 Days)

### ✅ HIGH-2A.1: Create Windows Launcher

**Priority**: High
**Effort**: S (2-3 hours)
**Dependencies**: Phase 1 complete
**Skills**: Batch scripting, Windows
**Assigned**: [Unassigned]

**Description**:
Create double-click launcher for Windows users.

**Implementation**:
- [ ] Create `launchers/launch_uorca.bat`
- [ ] Check if `uv` command exists
- [ ] Display error if uv not found
- [ ] Change to UORCA directory
- [ ] Display launch message
- [ ] Execute: `uv run uorca explore`
- [ ] Keep window open with `pause`

**Acceptance Criteria**:
- [ ] Double-clicking .bat launches UORCA
- [ ] Clear error if uv not installed
- [ ] Browser opens to UORCA Explorer
- [ ] Ctrl+C stops application
- [ ] Window remains open after exit

**Files**:
- `launchers/launch_uorca.bat` (new, ~40 lines)

**Testing**:
```bash
# On Windows:
# Double-click launch_uorca.bat
# Verify UORCA opens in browser
```

---

### ✅ HIGH-2A.2: Create macOS Launcher

**Priority**: High
**Effort**: S (2-3 hours)
**Dependencies**: Phase 1 complete
**Skills**: Bash scripting, macOS
**Assigned**: [Unassigned]

**Description**:
Create double-click launcher for macOS users.

**Implementation**:
- [ ] Create `launchers/launch_uorca.command`
- [ ] Add shebang: `#!/bin/bash`
- [ ] Get script directory
- [ ] Check if `uv` command exists
- [ ] Display error if uv not found
- [ ] Display launch message
- [ ] Execute: `uv run uorca explore`
- [ ] Keep terminal open with `read`

**Acceptance Criteria**:
- [ ] Double-clicking .command launches UORCA
- [ ] First time: User allows execution (Right-click → Open)
- [ ] Clear error if uv not installed
- [ ] Browser opens to UORCA Explorer
- [ ] Ctrl+C stops application

**Files**:
- `launchers/launch_uorca.command` (new, ~50 lines)

**Testing**:
```bash
# On macOS:
# chmod +x launchers/launch_uorca.command
# Double-click launch_uorca.command
# Verify UORCA opens in browser
```

---

### ✅ HIGH-2A.3: Create Linux Launcher

**Priority**: High
**Effort**: S (2-3 hours)
**Dependencies**: Phase 1 complete
**Skills**: Bash scripting, Linux
**Assigned**: [Unassigned]

**Description**:
Create launcher for Linux users with desktop environment integration.

**Implementation**:
- [ ] Create `launchers/launch_uorca.sh`
- [ ] Add shebang: `#!/bin/bash`
- [ ] Get script directory
- [ ] Check if `uv` command exists
- [ ] Display error if uv not found
- [ ] Check if running in terminal
- [ ] If not, open terminal (gnome-terminal or xterm)
- [ ] Display launch message
- [ ] Execute: `uv run uorca explore`

**Acceptance Criteria**:
- [ ] Right-click → Run launches UORCA
- [ ] Opens new terminal if not in one
- [ ] Clear error if uv not installed
- [ ] Browser opens to UORCA Explorer
- [ ] Works on Ubuntu, Fedora, Debian

**Files**:
- `launchers/launch_uorca.sh` (new, ~60 lines)

**Testing**:
```bash
# On Linux:
# chmod +x launchers/launch_uorca.sh
# Right-click → Run in Terminal
# Verify UORCA opens in browser
```

---

### ✅ HIGH-2A.4: Create Launcher README

**Priority**: High
**Effort**: S (1-2 hours)
**Dependencies**: HIGH-2A.1, HIGH-2A.2, HIGH-2A.3
**Skills**: Technical writing
**Assigned**: [Unassigned]

**Description**:
Document how to use launchers with troubleshooting guide.

**Implementation**:
- [ ] Create `launchers/README.md`
- [ ] Quick start section for each platform
- [ ] Requirements section (uv, Docker, API keys)
- [ ] Configuration instructions (.env file)
- [ ] Troubleshooting common issues
- [ ] Port configuration instructions

**Acceptance Criteria**:
- [ ] Clear instructions for Windows, macOS, Linux
- [ ] Prerequisites listed
- [ ] Common errors documented
- [ ] Solution steps provided

**Files**:
- `launchers/README.md` (new, ~100 lines)

**Testing**:
- [ ] Follow instructions on each platform
- [ ] Verify all troubleshooting steps work

---

### ✅ HIGH-2A.5: Set Execute Permissions

**Priority**: High
**Effort**: XS (15 minutes)
**Dependencies**: HIGH-2A.2, HIGH-2A.3
**Skills**: Shell, git
**Assigned**: [Unassigned]

**Description**:
Set executable permissions on launcher scripts for Unix systems.

**Implementation**:
- [ ] Run: `chmod +x launchers/launch_uorca.command`
- [ ] Run: `chmod +x launchers/launch_uorca.sh`
- [ ] Commit with executable permissions
- [ ] Verify permissions persist in git

**Acceptance Criteria**:
- [ ] Scripts executable after git clone
- [ ] Permissions: -rwxr-xr-x

**Files**:
- `launchers/launch_uorca.command`
- `launchers/launch_uorca.sh`

**Testing**:
```bash
ls -la launchers/
# Should show: -rwxr-xr-x for .command and .sh
```

---

### ✅ HIGH-2A.6: Update Main README with Launcher Instructions

**Priority**: High
**Effort**: S (1-2 hours)
**Dependencies**: HIGH-2A.4
**Skills**: Technical writing
**Assigned**: [Unassigned]

**Description**:
Update main README to prominently feature no-CLI launch method.

**Implementation**:
- [ ] Add "Quick Start (No Command Line!)" section at top
- [ ] Instructions for each platform
- [ ] Link to GUI_QUICKSTART.md
- [ ] Move CLI instructions lower
- [ ] Keep both options documented

**Acceptance Criteria**:
- [ ] GUI launch instructions prominent
- [ ] Each platform has clear steps
- [ ] CLI still documented
- [ ] Links work

**Files**:
- `README.md` (update, lines 1-50)

**Testing**:
- [ ] Follow README instructions
- [ ] Verify links work

---

### ✅ HIGH-2A.7: Cross-Platform Manual Testing

**Priority**: Critical
**Effort**: M (4-6 hours)
**Dependencies**: HIGH-2A.1, HIGH-2A.2, HIGH-2A.3
**Skills**: Testing, Multi-platform
**Assigned**: [Unassigned]

**Description**:
Test launchers on all three platforms with fresh installs.

**Test Scenarios**:
- [ ] Windows 10: Double-click .bat
- [ ] Windows 11: Double-click .bat
- [ ] macOS 12+: Double-click .command
- [ ] Ubuntu 20.04: Right-click → Run .sh
- [ ] Ubuntu 22.04: Right-click → Run .sh
- [ ] Test without uv installed → verify error
- [ ] Test without Docker running → verify graceful handling
- [ ] Test port 8501 busy → verify port conflict handling

**Acceptance Criteria**:
- [ ] All platforms launch successfully
- [ ] Error messages are helpful
- [ ] No platform-specific bugs
- [ ] Consistent behavior across platforms

**Files**:
- Testing checklist document (optional)

---

## 🟢 Phase 2B: Native Desktop App (4 Days) - OPTIONAL

### ✅ MED-2B.1: Create Desktop App Wrapper

**Priority**: Medium
**Effort**: M (4-6 hours)
**Dependencies**: Phase 1 complete
**Skills**: Python, PyWebView
**Assigned**: [Unassigned]

**Description**:
Create Python script to wrap Streamlit in native window.

**Implementation**:
- [ ] Create `desktop/uorca_desktop.py`
- [ ] Implement port availability check
- [ ] Implement find_available_port()
- [ ] Implement run_streamlit() subprocess launcher
- [ ] Wait for Streamlit to start (health check)
- [ ] Create PyWebView window
- [ ] Handle window close → terminate Streamlit

**Acceptance Criteria**:
- [ ] Finds available port automatically
- [ ] Starts Streamlit in background
- [ ] Opens native window (not browser)
- [ ] Closes cleanly (kills Streamlit)

**Files**:
- `desktop/uorca_desktop.py` (new, ~150 lines)

**Testing**:
```bash
python desktop/uorca_desktop.py
# Verify native window opens
```

---

### ✅ MED-2B.2: Create PyInstaller Build Script

**Priority**: Medium
**Effort**: M (4-6 hours)
**Dependencies**: MED-2B.1
**Skills**: PyInstaller, Packaging
**Assigned**: [Unassigned]

**Description**:
Create build script to package desktop app as executable.

**Implementation**:
- [ ] Create `desktop/build_desktop.py`
- [ ] Detect platform (Windows, macOS, Linux)
- [ ] Set platform-specific icon
- [ ] Configure PyInstaller with:
  - --onefile
  - --windowed
  - --icon
  - --add-data (include code)
  - --hidden-import (all dependencies)
  - --collect-all streamlit
- [ ] Build executable

**Acceptance Criteria**:
- [ ] Build script runs on all platforms
- [ ] Executable created in dist/
- [ ] All code and assets bundled
- [ ] Executable runs standalone

**Files**:
- `desktop/build_desktop.py` (new, ~80 lines)

**Testing**:
```bash
cd desktop
python build_desktop.py
./dist/UORCA  # or UORCA.exe on Windows
```

---

### ✅ MED-2B.3: Add Desktop Dependencies to pyproject.toml

**Priority**: Medium
**Effort**: XS (30 minutes)
**Dependencies**: None
**Skills**: Python, Poetry/uv
**Assigned**: [Unassigned]

**Description**:
Add optional desktop dependencies to project config.

**Implementation**:
- [ ] Add `[project.optional-dependencies]` section
- [ ] Add desktop group: pywebview>=4.0, pyinstaller>=6.0, requests>=2.31
- [ ] Document installation: `uv pip install -e ".[desktop]"`

**Acceptance Criteria**:
- [ ] Optional dependencies defined
- [ ] Installation works
- [ ] Doesn't affect normal installation

**Files**:
- `pyproject.toml` (update)

**Testing**:
```bash
uv pip install -e ".[desktop]"
python -c "import webview; import PyInstaller"
```

---

### ✅ MED-2B.4: Create Desktop Icons

**Priority**: Low
**Effort**: S (2-3 hours)
**Dependencies**: None
**Skills**: Design, Icon creation
**Assigned**: [Unassigned]

**Description**:
Create application icons for all platforms.

**Implementation**:
- [ ] Design UORCA icon (512x512 PNG)
- [ ] Convert to Windows .ico (multiple sizes)
- [ ] Convert to macOS .icns (multiple sizes)
- [ ] Save PNG for Linux
- [ ] Store in resources/ directory

**Acceptance Criteria**:
- [ ] Icon looks professional
- [ ] All formats created
- [ ] Icons work in built executables

**Files**:
- `resources/uorca_icon.png` (new)
- `resources/uorca_icon.ico` (new)
- `resources/uorca_icon.icns` (new)

**Testing**:
- [ ] Build executables with icons
- [ ] Verify icons display correctly

---

### ✅ MED-2B.5: Build and Test Windows Executable

**Priority**: Medium
**Effort**: M (3-4 hours)
**Dependencies**: MED-2B.2, MED-2B.3, MED-2B.4
**Skills**: Windows, Packaging
**Assigned**: [Unassigned]

**Description**:
Build and test standalone Windows executable.

**Implementation**:
- [ ] Install desktop dependencies on Windows
- [ ] Run build script
- [ ] Test executable on clean Windows VM
- [ ] Verify all features work
- [ ] Test installer (optional)

**Acceptance Criteria**:
- [ ] Executable builds without errors
- [ ] Runs on Windows 10/11
- [ ] No dependencies required
- [ ] Icon displays correctly

**Files**:
- `dist/UORCA.exe` (generated)

**Testing**:
- [ ] Build on Windows
- [ ] Run on fresh Windows install
- [ ] Complete full workflow

---

### ✅ MED-2B.6: Build and Test macOS/Linux Executables

**Priority**: Medium
**Effort**: M (3-4 hours)
**Dependencies**: MED-2B.2, MED-2B.3, MED-2B.4
**Skills**: macOS, Linux, Packaging
**Assigned**: [Unassigned]

**Description**:
Build and test standalone macOS and Linux executables.

**Implementation**:
- [ ] Install desktop dependencies on macOS
- [ ] Build macOS app bundle
- [ ] Sign app (optional)
- [ ] Install desktop dependencies on Linux
- [ ] Build Linux executable
- [ ] Create .desktop file (optional)

**Acceptance Criteria**:
- [ ] macOS app builds and runs
- [ ] Linux executable builds and runs
- [ ] Icons display correctly
- [ ] All features work

**Files**:
- `dist/UORCA.app` (macOS)
- `dist/uorca` (Linux)

**Testing**:
- [ ] Build on macOS and Linux
- [ ] Test on fresh installs
- [ ] Complete full workflow

---

## 🔵 Phase 3: Polish & Testing (Week 5)

### ✅ HIGH-3.1: Run Complete Automated Test Suite

**Priority**: High
**Effort**: M (4-6 hours)
**Dependencies**: All implementation complete
**Skills**: Testing, pytest
**Assigned**: [Unassigned]

**Description**:
Run all automated tests and achieve coverage targets.

**Implementation**:
- [ ] Run: `pytest tests/ -v`
- [ ] Fix any failing tests
- [ ] Run: `pytest --cov=main_workflow/reporting tests/`
- [ ] Verify coverage >80%
- [ ] Run: `pyright main_workflow/reporting/`
- [ ] Fix type errors
- [ ] Generate coverage report

**Acceptance Criteria**:
- [ ] All tests pass
- [ ] Coverage >80% overall
- [ ] No type errors
- [ ] Performance benchmarks pass

**Files**:
- All test files

**Testing**:
```bash
pytest tests/ -v
pytest --cov=main_workflow/reporting --cov-report=html tests/
pyright main_workflow/reporting/
```

---

### ✅ HIGH-3.2: Complete Manual User Testing

**Priority**: Critical
**Effort**: L (8-10 hours)
**Dependencies**: All implementation complete
**Skills**: QA, User testing
**Assigned**: [Unassigned]

**Description**:
Perform comprehensive manual testing with real-world scenarios.

**Test Scenarios**:
- [ ] Scenario 1: First-time user setup (fresh install)
- [ ] Scenario 2: Dataset identification workflow (small query)
- [ ] Scenario 3: Pipeline execution (2-3 datasets)
- [ ] Scenario 4: Cross-platform launcher testing
- [ ] Scenario 5: Error recovery (invalid inputs, crashes)
- [ ] Scenario 6: Long-running pipeline (5+ hours)
- [ ] Scenario 7: Task history persistence
- [ ] Scenario 8: Concurrent operations

**Acceptance Criteria**:
- [ ] All scenarios complete successfully
- [ ] No crashes or data loss
- [ ] Error messages are helpful
- [ ] Performance is acceptable
- [ ] UI remains responsive

**Files**:
- Testing checklist/report (optional)

---

### ✅ HIGH-3.3: Create User Documentation

**Priority**: High
**Effort**: M (6-8 hours)
**Dependencies**: All implementation complete
**Skills**: Technical writing
**Assigned**: [Unassigned]

**Description**:
Write comprehensive documentation for end users.

**Implementation**:
- [ ] Create `docs/GUI_QUICKSTART.md`
- [ ] Prerequisites section
- [ ] Setup instructions (.env configuration)
- [ ] Using UORCA section (3 steps: Identify → Run → Explore)
- [ ] Troubleshooting guide
- [ ] Create `docs/CLI_TO_GUI_MIGRATION.md`
- [ ] Command equivalents
- [ ] Advantages of GUI
- [ ] When to use CLI
- [ ] Migration strategy

**Acceptance Criteria**:
- [ ] Documentation is clear and complete
- [ ] New users can follow without help
- [ ] Common issues documented
- [ ] Screenshots included (optional)

**Files**:
- `docs/GUI_QUICKSTART.md` (new, ~200 lines)
- `docs/CLI_TO_GUI_MIGRATION.md` (new, ~150 lines)

**Testing**:
- [ ] Ask new user to follow quickstart
- [ ] Verify they succeed without assistance

---

### ✅ HIGH-3.4: Performance Optimization and Bug Fixes

**Priority**: High
**Effort**: M (6-8 hours)
**Dependencies**: HIGH-3.2 (testing complete)
**Skills**: Optimization, Debugging
**Assigned**: [Unassigned]

**Description**:
Fix any bugs found during testing and optimize performance.

**Tasks**:
- [ ] Fix all critical bugs from testing
- [ ] Optimize SQLite queries (<50ms)
- [ ] Optimize progress refresh rate (minimize reruns)
- [ ] Optimize memory usage for long-running tasks
- [ ] Fix UI responsiveness issues
- [ ] Profile slow operations
- [ ] Add loading indicators where needed

**Acceptance Criteria**:
- [ ] No critical bugs remain
- [ ] Database queries <50ms
- [ ] UI remains responsive during long operations
- [ ] Memory usage stable
- [ ] All performance benchmarks pass

**Files**:
- Various files as needed

**Testing**:
```bash
# Profile task manager
python -m cProfile -s cumtime test_task_manager.py

# Memory profiling
memory_profiler task_manager.py
```

---

## 📈 Milestones

### Milestone 1: Background Task System (End of Week 1)
**Target**: Day 5
**Tasks**: CRIT-1A.1 through CRIT-1A.8
**Success**: TaskManager working with unit tests passing

### Milestone 2: Dataset Identification GUI (End of Week 2)
**Target**: Day 10
**Tasks**: CRIT-1B.1 through CRIT-1B.10
**Success**: Can identify datasets entirely in GUI

### Milestone 3: Pipeline Execution GUI (End of Week 4)
**Target**: Day 20
**Tasks**: CRIT-1C.1 through CRIT-1C.12
**Success**: Can run pipelines entirely in GUI

### Milestone 4: Desktop Launchers (End of Week 4.5)
**Target**: Day 23
**Tasks**: HIGH-2A.1 through HIGH-2A.7
**Success**: Double-click launch works on all platforms

### Milestone 5: Production Ready (End of Week 6)
**Target**: Day 30
**Tasks**: HIGH-3.1 through HIGH-3.4
**Success**: Tests pass, docs complete, no critical bugs

---

## 🔄 Daily Standups Template

### Today's Focus: [Date]
- [ ] Task being worked on: [Task ID and name]
- [ ] Expected completion: [Today/Tomorrow/This week]
- [ ] Blockers: [None / List blockers]
- [ ] Completed yesterday: [Task IDs]

---

## 🐛 Known Issues / Blockers

*To be filled in during implementation*

- [ ] Issue #1: [Description]
  - **Impact**: [High/Medium/Low]
  - **Status**: [Open/In Progress/Resolved]
  - **Resolution**: [Approach]

---

## ⚠️ Risks & Concerns

### Technical Risks

**Risk 1: Python Threading Limitations**
- **Concern**: Cannot force-kill threads, only cooperative cancellation
- **Mitigation**: Document limitation, implement best-effort cancellation, add progress_callback checks
- **Impact**: Medium
- **Status**: Accepted

**Risk 2: Long-Running Pipeline Monitoring**
- **Concern**: Auto-refresh may cause browser performance issues on multi-hour runs
- **Mitigation**: Implement smart refresh (exponential backoff), allow manual refresh
- **Impact**: Medium
- **Status**: Mitigated

**Risk 3: Cross-Platform Launcher Compatibility**
- **Concern**: Different shell behaviors on Windows/macOS/Linux
- **Mitigation**: Thorough testing on all platforms, clear error messages
- **Impact**: High
- **Status**: Requires testing

### Timeline Risks

**Risk 4: Native Desktop App Complexity**
- **Concern**: PyInstaller packaging can be unpredictable
- **Mitigation**: Made Phase 2B optional, batch launchers sufficient for MVP
- **Impact**: Low (optional feature)
- **Status**: Deferred

**Risk 5: Integration Testing Time**
- **Concern**: Real pipeline runs take hours, hard to test
- **Mitigation**: Mock tests for speed, smaller datasets for integration tests
- **Impact**: Medium
- **Status**: Planned for

---

## 📚 References

- **Original Document**: `docs/plans/no_cli_uorca_implementation.md`
- **Current State**: 30% complete (existing UORCA Explorer, backend batch processing, Docker container)
- **Architecture**: Threading + SQLite (zero external dependencies)
- **Design Pattern**: Singleton TaskManager, background threads, polling for status

---

## 🔧 Development Setup

### Prerequisites
```bash
# Install uv
curl -LsSf https://astral.sh/uv/install.sh | sh  # macOS/Linux
# or download from https://docs.astral.sh/uv/  # Windows

# Install Docker Desktop
# Download from https://www.docker.com/products/docker-desktop/

# Clone UORCA
git clone https://github.com/uorca/uorca.git
cd uorca
```

### Setup
```bash
# Install dependencies
uv sync

# Create .env file
cat > .env <<EOF
ENTREZ_EMAIL=your.email@institution.edu
OPENAI_API_KEY=sk-proj-your-key-here
ENTREZ_API_KEY=your_ncbi_key  # Optional
EOF

# Run tests
pytest tests/

# Launch UORCA Explorer
uv run uorca explore
```

### Development Commands
```bash
# Run specific test file
pytest tests/unit/reporting/test_task_manager.py -v

# Run with coverage
pytest --cov=main_workflow/reporting tests/

# Type checking
pyright main_workflow/reporting/

# Format code
black main_workflow/reporting/

# Lint code
ruff check main_workflow/reporting/

# Run UORCA Explorer on custom port
uv run streamlit run main_workflow/reporting/uorca_explorer.py --server.port 8502
```

---

## 💡 Tips for Contributors

### Before Starting a Task
- [ ] Read the original implementation plan section
- [ ] Check dependencies are complete
- [ ] Review related code in codebase
- [ ] Create git branch: `feature/task-id-description`

### During Implementation
- [ ] Follow existing code style
- [ ] Add docstrings to all functions
- [ ] Add type hints where possible
- [ ] Write tests alongside code (TDD)
- [ ] Test manually before marking complete
- [ ] Update TODO status in this document

### Before Marking Complete
- [ ] All acceptance criteria met
- [ ] Unit tests written and passing
- [ ] Manual testing completed
- [ ] Code reviewed (self-review at minimum)
- [ ] Documentation updated (if needed)
- [ ] Commit with clear message

### Git Workflow
```bash
# Start new task
git checkout -b feature/crit-1a-1-task-manager-core

# Regular commits
git add main_workflow/reporting/core/task_manager.py
git commit -m "feat: implement TaskManager core class with SQLite persistence"

# Before PR
pytest tests/
pyright main_workflow/reporting/

# Submit PR
git push origin feature/crit-1a-1-task-manager-core
# Create PR on GitHub
```

---

## 📊 Velocity Tracking

*To be filled in during implementation*

| Week | Planned Tasks | Completed Tasks | Velocity | Notes |
|------|---------------|-----------------|----------|-------|
| Week 1 | 8 (Phase 1A) | 0 | 0% | Background Task System |
| Week 2 | 10 (Phase 1B) | 0 | 0% | Identification Tab |
| Week 3 | 6 (Phase 1C.1-6) | 0 | 0% | Pipeline Tab (Part 1) |
| Week 4 | 6 (Phase 1C.7-12) | 0 | 0% | Pipeline Tab (Part 2) |
| Week 5 | 11 (Phase 2A) | 0 | 0% | Launchers + Testing |
| Week 6 | 6 (Phase 3) | 0 | 0% | Polish & Documentation |

**Average velocity**: TBD after Week 1
**Adjusted timeline**: TBD after Week 2

---

## 🎯 Success Metrics

### Automated Metrics

**Test Coverage**:
- [ ] Unit tests: >80% coverage
- [ ] Integration tests: All critical paths covered
- [ ] Type checking: Zero pyright errors

**Performance**:
- [ ] TaskManager overhead: <100ms per operation
- [ ] UI responsiveness: <1s to update progress
- [ ] Database queries: <50ms average
- [ ] Memory usage: Stable over multi-hour runs

**Command to Check**:
```bash
pytest tests/ --cov=main_workflow/reporting --cov-report=term
pyright main_workflow/reporting/
```

### Manual Validation Checklist

**Usability** (Test with non-technical users):
- [ ] Can launch UORCA without command line help
- [ ] Can complete dataset identification without instructions
- [ ] Can run pipeline and understand progress
- [ ] Can navigate to results and export data
- [ ] Error messages are helpful (not cryptic)

**Reliability**:
- [ ] No crashes during 10-hour pipeline runs
- [ ] Task history persists across restarts
- [ ] Cancellation stops tasks within 10 seconds
- [ ] Failed datasets don't block queue
- [ ] Database never corrupts

**Cross-Platform**:
- [ ] Launchers work on Windows 10/11
- [ ] Launchers work on macOS 12+
- [ ] Launchers work on Ubuntu 20.04+
- [ ] Docker integration works on all platforms
- [ ] UI looks correct on all platforms

---

## 🚀 Future Enhancements (Phase 4+)

**Not in scope for this implementation, but documented for future work:**

1. **SLURM Support in GUI**: Add batch system selector in Pipeline tab
2. **Multi-User Deployment**: Upgrade to Redis for shared job queue
3. **Cloud Integration**: AWS Batch, Google Cloud integration
4. **Advanced Monitoring**: Real-time log streaming, resource graphs, container stats
5. **Mobile Access**: Responsive design for tablet/phone monitoring
6. **Collaboration Features**: Share results, comments, annotations
7. **Desktop App Improvements**: Auto-updater, system tray, notifications, offline mode
8. **Performance Optimization**: WebSocket for live updates, caching, lazy loading
9. **Advanced Analytics**: Usage tracking, error analytics, performance metrics
10. **Plugin System**: Allow custom analysis modules

---

## ❓ Questions & Support

- **Issues**: Report bugs at [GitHub Issues](https://github.com/uorca/uorca/issues)
- **Discussions**: Ask questions at [GitHub Discussions](https://github.com/uorca/uorca/discussions)
- **Documentation**: See `docs/` directory
- **This TODO**: `todos/no-cli-uorca-implementation-todos.md`

---

## 📝 Change Log

| Date | Change | By |
|------|--------|-----|
| 2025-10-01 | Initial TODO list generated | Claude |
|  |  |  |

---

*This TODO list is a living document. Update it as work progresses and requirements evolve.*

**Last Updated**: 2025-10-01
**Total Tasks**: 47
**Estimated Completion**: 2025-11-12 (6 weeks)
**Status**: 🔴 Not Started (0% complete)
