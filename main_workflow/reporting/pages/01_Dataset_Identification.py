"""
Dataset Identification Page for UORCA Explorer.

Provides GUI interface for finding relevant RNA-seq datasets from GEO.
"""

import os
import sys
import streamlit as st
from pathlib import Path
from datetime import datetime
import pandas as pd

# Add parent directories to path for imports
script_dir = Path(__file__).parent
reporting_dir = script_dir.parent
sys.path.insert(0, str(reporting_dir))

from core import TaskManager, TaskStatus

# Configure page
st.set_page_config(
    page_title="Dataset Identification - UORCA",
    page_icon="🔍",
    layout="wide"
)


def main():
    """Main page function."""
    st.title("🔍 Dataset Identification")
    st.markdown("""
    Find relevant RNA-seq datasets from GEO based on your research question.
    UORCA will search, validate, cluster, and assess datasets using AI.
    """)

    # Initialize task manager
    task_manager = TaskManager()

    # Check if there's a running identification task
    if "current_identify_task" in st.session_state:
        task_id = st.session_state.current_identify_task
        status = task_manager.get_task_status(task_id)

        if status and status["status"] in (TaskStatus.RUNNING, TaskStatus.PENDING):
            show_running_task(task_id, status, task_manager)
            return

    # Main input form
    with st.form("identify_form"):
        st.subheader("Search Parameters")

        # Research question
        research_query = st.text_area(
            "Research Question",
            placeholder="Example: What are the transcriptomic changes in response to hypoxia in cancer cells?",
            help="Describe your research question. UORCA will use AI to extract relevant search terms.",
            height=100
        )

        col1, col2 = st.columns(2)

        with col1:
            max_results = st.number_input(
                "Maximum Results per Search Term",
                min_value=10,
                max_value=500,
                value=100,
                step=10,
                help="How many datasets to retrieve for each search term"
            )

            num_datasets = st.number_input(
                "Number of Representative Datasets",
                min_value=1,
                max_value=50,
                value=10,
                step=1,
                help="How many datasets to select after clustering"
            )

        with col2:
            model = st.selectbox(
                "OpenAI Model",
                options=[
                    "gpt-5-mini",
                    "gpt-5",
                    "gpt-5-nano"
                ],
                index=0,
                help="Model for term extraction and relevance assessment (gpt-5-mini is faster and cheaper)"
            )

            samples_required = st.number_input(
                "Minimum Samples per Dataset",
                min_value=2,
                max_value=10,
                value=3,
                step=1,
                help="Minimum number of paired-end RNA-seq samples required"
            )

        # Output directory
        project_root = Path(__file__).parent.parent.parent.parent
        default_output = project_root / "scratch" / "UORCA_Identification" / datetime.now().strftime("%Y%m%d_%H%M%S")

        output_dir = st.text_input(
            "Output Directory",
            value=str(default_output),
            help="Where to save identification results"
        )

        # Submit button
        submitted = st.form_submit_button("🔍 Start Identification", type="primary", use_container_width=True)

    if submitted:
        if not research_query.strip():
            st.error("Please enter a research question")
            return

        # Validate API keys
        if not os.getenv("OPENAI_API_KEY"):
            st.error("OPENAI_API_KEY not found in environment. Please set it in your .env file.")
            return
        if not os.getenv("ENTREZ_EMAIL"):
            st.error("ENTREZ_EMAIL not found in environment. Please set it in your .env file.")
            return

        # Create output directory
        output_path = Path(output_dir)
        try:
            output_path.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            st.error(f"Could not create output directory: {e}")
            return

        # Submit task
        task_id = f"identify_{datetime.now().strftime('%Y%m%d_%H%M%S')}"

        task_manager.submit_task(
            task_id=task_id,
            task_type="identify",
            task_func=run_identification,
            parameters={
                "research_query": research_query,
                "max_results": max_results,
                "num_datasets": num_datasets,
                "model": model,
                "samples_required": samples_required,
                "output_dir": str(output_path)
            }
        )

        st.session_state.current_identify_task = task_id
        st.rerun()

    # Show recent tasks
    show_task_history(task_manager)


def show_running_task(task_id: str, status: dict, task_manager: TaskManager):
    """Display progress for a running identification task."""

    st.subheader("🔄 Identification in Progress")

    # Progress bar
    progress = status.get("progress", 0.0)
    st.progress(progress)

    # Status message
    message = status.get("progress_message", "Processing...")
    st.info(message)

    # Cancel button
    if st.button("❌ Cancel", type="secondary"):
        if task_manager.cancel_task(task_id):
            st.warning("Cancellation requested. Task may take a moment to stop.")
            del st.session_state.current_identify_task
            st.rerun()

    # Auto-refresh while running
    if status["status"] == TaskStatus.RUNNING:
        st.rerun()
    elif status["status"] == TaskStatus.COMPLETED:
        st.success("✅ Identification complete!")

        result_path = status.get("result")
        if result_path:
            st.markdown(f"**Results saved to:** `{result_path}`")

            # Show results preview
            csv_path = Path(result_path) / "Dataset_identification_result.csv"
            if csv_path.exists():
                df = pd.read_csv(csv_path)
                st.dataframe(df, use_container_width=True)

                # Download button
                st.download_button(
                    "📥 Download Results CSV",
                    data=df.to_csv(index=False),
                    file_name="dataset_identification_results.csv",
                    mime="text/csv"
                )

        if st.button("Start New Identification"):
            del st.session_state.current_identify_task
            st.rerun()

    elif status["status"] == TaskStatus.FAILED:
        st.error(f"❌ Identification failed")

        error_msg = status.get("error", "Unknown error")
        with st.expander("Error Details"):
            st.code(error_msg)

        if st.button("Try Again"):
            del st.session_state.current_identify_task
            st.rerun()


def show_task_history(task_manager: TaskManager):
    """Show recent identification tasks."""

    st.divider()
    st.subheader("📜 Recent Identifications")

    tasks = task_manager.list_tasks(task_type="identify", limit=10)

    if not tasks:
        st.info("No previous identifications found")
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
                query_preview = params.get('research_query', 'N/A')
                if len(query_preview) > 100:
                    query_preview = query_preview[:100] + "..."
                st.markdown(f"**Query:** {query_preview}")
                st.markdown(f"**Model:** {params.get('model', 'N/A')}")
                st.markdown(f"**Output:** {task.get('result_path', 'N/A')}")

            if task['status'] == 'failed':
                with st.expander("Error Details"):
                    st.code(task.get('error_message', 'Unknown'))


def run_identification(
    research_query: str,
    max_results: int,
    num_datasets: int,
    model: str,
    samples_required: int,
    output_dir: str,
    progress_callback=None
) -> str:
    """
    Run dataset identification in background thread.

    Args:
        research_query: Research question
        max_results: Max results per search term
        num_datasets: Number of representative datasets to select
        model: OpenAI model to use
        samples_required: Minimum samples per dataset (not used currently)
        output_dir: Output directory
        progress_callback: Function to call with (progress, message)

    Returns:
        Path to output directory
    """
    project_root = Path(__file__).parent.parent.parent.parent
    sys.path.insert(0, str(project_root))

    try:
        from uorca.identify import main as identify_main

        if progress_callback:
            progress_callback(0.1, "Starting dataset identification...")

        # Save original sys.argv
        original_argv = sys.argv.copy()

        # Build command line arguments for identify CLI
        sys.argv = ['identify']
        sys.argv.extend(['-q', research_query])
        sys.argv.extend(['-o', output_dir])
        sys.argv.extend(['-m', str(max_results)])
        sys.argv.extend(['-a', str(num_datasets)])
        sys.argv.extend(['--model', model])

        if progress_callback:
            progress_callback(0.2, "Extracting search terms and searching GEO...")

        # Call the identify main function
        identify_main()

        if progress_callback:
            progress_callback(1.0, "Identification complete")

        # Restore original sys.argv
        sys.argv = original_argv

        return output_dir

    finally:
        if str(project_root) in sys.path:
            sys.path.remove(str(project_root))


if __name__ == "__main__":
    main()
