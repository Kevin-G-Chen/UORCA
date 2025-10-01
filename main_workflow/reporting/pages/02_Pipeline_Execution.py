"""
Pipeline Execution Page for UORCA Explorer.

Provides GUI interface for running UORCA RNA-seq analysis pipeline.
"""

import os
import sys
import streamlit as st
from pathlib import Path
from datetime import datetime
import pandas as pd
import time

# Add parent directories to path for imports
script_dir = Path(__file__).parent
reporting_dir = script_dir.parent
sys.path.insert(0, str(reporting_dir))
sys.path.insert(0, str(reporting_dir.parent.parent))  # Add project root

from core import TaskManager, TaskStatus

# Configure page
st.set_page_config(
    page_title="Pipeline Execution - UORCA",
    page_icon="🚀",
    layout="wide"
)


def main():
    """Main page function."""
    st.title("🚀 Run UORCA Pipeline")
    st.markdown("""
    Run RNA-seq analysis on selected datasets. Upload a CSV with dataset accessions
    or use results from Dataset Identification.
    """)

    # Initialize task manager
    task_manager = TaskManager()

    # Check for running pipeline
    if "current_pipeline_task" in st.session_state:
        task_id = st.session_state.current_pipeline_task
        status = task_manager.get_task_status(task_id)

        if status and status["status"] in (TaskStatus.RUNNING, TaskStatus.PENDING):
            show_running_pipeline(task_id, status, task_manager)
            return

    # Input source selection
    st.subheader("Input Dataset")

    input_method = st.radio(
        "Select input method:",
        ["Upload CSV", "Use Identification Results"],
        horizontal=True
    )

    input_csv_path = None

    if input_method == "Upload CSV":
        uploaded_file = st.file_uploader(
            "Upload dataset CSV",
            type=["csv"],
            help="CSV must contain 'Accession' and 'SafeStorageGB' columns"
        )

        if uploaded_file:
            # Save to temp location
            temp_dir = Path.home() / ".uorca" / "temp"
            temp_dir.mkdir(parents=True, exist_ok=True)

            temp_csv = temp_dir / f"upload_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
            temp_csv.write_bytes(uploaded_file.read())
            input_csv_path = str(temp_csv)

            # Preview
            df = pd.read_csv(temp_csv)
            st.dataframe(df.head(), use_container_width=True)

    else:  # Use Identification Results
        # List recent identification tasks
        identify_tasks = task_manager.list_tasks(task_type="identify", limit=20)
        completed_tasks = [t for t in identify_tasks if t['status'] == 'completed']

        if not completed_tasks:
            st.warning("No completed identification tasks found. Run 'Dataset Identification' first.")
            return

        task_options = {
            f"{t['task_id']} ({t['completed_at']})": t['result_path']
            for t in completed_tasks
        }

        selected_task = st.selectbox(
            "Select identification result:",
            options=list(task_options.keys())
        )

        if selected_task:
            result_dir = Path(task_options[selected_task])
            csv_path = result_dir / "selected_datasets.csv"

            if csv_path.exists():
                input_csv_path = str(csv_path)
                df = pd.read_csv(csv_path)
                st.dataframe(df, use_container_width=True)
                st.info(f"Found {len(df)} datasets to process")
            else:
                st.error(f"Could not find selected_datasets.csv in {result_dir}")
                return

    if not input_csv_path:
        st.info("Please provide input datasets to continue")
        return

    # Pipeline configuration form
    with st.form("pipeline_form"):
        st.subheader("Pipeline Configuration")

        col1, col2 = st.columns(2)

        with col1:
            output_dir = st.text_input(
                "Output Directory",
                value=str(Path.home() / "UORCA_Results" / datetime.now().strftime("%Y%m%d_%H%M%S")),
                help="Where to save analysis results"
            )

            max_workers = st.number_input(
                "Parallel Workers",
                min_value=1,
                max_value=16,
                value=4,
                help="Number of datasets to process in parallel"
            )

        with col2:
            max_storage_gb = st.number_input(
                "Storage Limit (GB)",
                min_value=10,
                max_value=1000,
                value=100,
                help="Maximum disk space to use at once"
            )

            container_tmpfs_gb = st.number_input(
                "Container Temp Storage (GB)",
                min_value=5,
                max_value=100,
                value=20,
                help="Temporary storage per container"
            )

        # Advanced options
        with st.expander("Advanced Options"):
            docker_image = st.text_input(
                "Docker Image",
                value="kevingchen/uorca:0.1.0",
                help="UORCA Docker container image"
            )

            cleanup = st.checkbox(
                "Cleanup Intermediate Files",
                value=True,
                help="Remove intermediate files after processing"
            )

        # Submit button
        submitted = st.form_submit_button("🚀 Start Pipeline", type="primary", use_container_width=True)

    if submitted:
        # Validate CSV
        try:
            df = pd.read_csv(input_csv_path)
            required_cols = ['Accession', 'SafeStorageGB']
            missing = [col for col in required_cols if col not in df.columns]
            if missing:
                st.error(f"CSV missing required columns: {missing}")
                return
        except Exception as e:
            st.error(f"Invalid CSV file: {e}")
            return

        # Validate API keys
        if not os.getenv("ENTREZ_EMAIL"):
            st.error("ENTREZ_EMAIL not found in environment")
            return

        # Create output directory
        output_path = Path(output_dir)
        try:
            output_path.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            st.error(f"Could not create output directory: {e}")
            return

        # Submit pipeline task
        task_id = f"pipeline_{datetime.now().strftime('%Y%m%d_%H%M%S')}"

        task_manager.submit_task(
            task_id=task_id,
            task_type="pipeline",
            task_func=run_pipeline,
            parameters={
                "input_csv": input_csv_path,
                "output_dir": str(output_path),
                "max_workers": max_workers,
                "max_storage_gb": max_storage_gb,
                "container_tmpfs_gb": container_tmpfs_gb,
                "docker_image": docker_image,
                "cleanup": cleanup
            }
        )

        st.session_state.current_pipeline_task = task_id
        st.rerun()

    # Show recent pipelines
    show_pipeline_history(task_manager)


def show_running_pipeline(task_id: str, status: dict, task_manager: TaskManager):
    """Display progress for running pipeline."""

    st.subheader("🔄 Pipeline Running")

    # Overall progress
    progress = status.get("progress", 0.0)
    st.progress(progress, text=f"{int(progress * 100)}% complete")

    # Status message
    message = status.get("progress_message", "Processing...")
    st.info(message)

    # Action buttons
    col1, col2 = st.columns([1, 4])

    with col1:
        if st.button("❌ Cancel Pipeline", type="secondary"):
            if task_manager.cancel_task(task_id):
                st.warning("Cancellation requested")
                del st.session_state.current_pipeline_task
                st.rerun()

    with col2:
        if st.button("🔄 Refresh Status"):
            st.rerun()

    # Auto-refresh while running
    if status["status"] == TaskStatus.RUNNING:
        time.sleep(5)  # Refresh every 5 seconds
        st.rerun()

    elif status["status"] == TaskStatus.COMPLETED:
        st.success("✅ Pipeline complete!")

        result_path = status.get("result")
        if result_path:
            st.markdown(f"**Results saved to:** `{result_path}`")

            # Quick navigation to main explorer
            if st.button("📊 View Results in Explorer", type="primary"):
                st.session_state.results_directory = result_path
                del st.session_state.current_pipeline_task
                # Switch to main app - user will need to navigate manually
                st.info("Navigate to the main UORCA Explorer page to view results")

        if st.button("Start New Pipeline"):
            del st.session_state.current_pipeline_task
            st.rerun()

    elif status["status"] == TaskStatus.FAILED:
        st.error(f"❌ Pipeline failed")

        error_msg = status.get("error", "Unknown error")
        with st.expander("Error Details"):
            st.code(error_msg)

        if st.button("Try Again"):
            del st.session_state.current_pipeline_task
            st.rerun()


def show_pipeline_history(task_manager: TaskManager):
    """Show recent pipeline executions."""

    st.divider()
    st.subheader("📜 Recent Pipeline Runs")

    tasks = task_manager.list_tasks(task_type="pipeline", limit=10)

    if not tasks:
        st.info("No previous pipeline runs found")
        return

    for task in tasks:
        with st.expander(f"{task['task_id']} - {task['status']} ({task['created_at']})"):
            col1, col2 = st.columns(2)

            with col1:
                st.markdown(f"**Status:** {task['status']}")
                st.markdown(f"**Started:** {task['started_at'] or 'N/A'}")
                st.markdown(f"**Completed:** {task['completed_at'] or 'N/A'}")

            with col2:
                params = task.get('parameters', {})
                st.markdown(f"**Workers:** {params.get('max_workers', 'N/A')}")
                st.markdown(f"**Storage Limit:** {params.get('max_storage_gb', 'N/A')} GB")
                st.markdown(f"**Output:** {params.get('output_dir', 'N/A')}")

            if task['status'] == 'failed':
                with st.expander("Error Details"):
                    st.code(task.get('error_message', 'Unknown'))


def run_pipeline(
    input_csv: str,
    output_dir: str,
    max_workers: int,
    max_storage_gb: float,
    container_tmpfs_gb: float,
    docker_image: str,
    cleanup: bool,
    progress_callback=None
) -> str:
    """
    Run UORCA pipeline in background thread.

    Returns:
        Output directory path
    """
    # Add project root to path
    project_root = Path(__file__).parent.parent.parent.parent
    sys.path.insert(0, str(project_root))

    try:
        from uorca.batch.local import LocalBatchProcessor

        # Initialize processor
        processor = LocalBatchProcessor()

        # Update progress
        if progress_callback:
            progress_callback(0.1, "Initializing pipeline...")

        # Submit datasets - this runs synchronously
        processor.submit_datasets(
            input_path=input_csv,
            output_dir=output_dir,
            max_workers=max_workers,
            max_storage_gb=max_storage_gb,
            container_tmpfs_gb=container_tmpfs_gb,
            docker_image=docker_image,
            cleanup=cleanup
        )

        # Pipeline completed
        if progress_callback:
            progress_callback(1.0, "Pipeline completed")

        return output_dir

    finally:
        if str(project_root) in sys.path:
            sys.path.remove(str(project_root))


if __name__ == "__main__":
    main()
