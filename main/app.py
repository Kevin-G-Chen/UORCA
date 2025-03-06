"""
Flask application for the UORCA web interface.
"""
import logging
from flask import Flask, jsonify, request, abort
from pathlib import Path
import json
import os

from UORCA.config import get_settings
from UORCA.services.workflow_service import WorkflowService

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)

logger = logging.getLogger(__name__)

# Create and configure the app
app = Flask(__name__)
app.config.from_mapping(
    SECRET_KEY=os.environ.get("SECRET_KEY", "dev"),
)

# Initialize services
workflow_service = WorkflowService()
settings = get_settings()

@app.route("/api/health", methods=["GET"])
def health_check():
    """Health check endpoint."""
    return jsonify({"status": "ok"})

@app.route("/api/workflows", methods=["POST"])
def create_workflow():
    """Create a new workflow."""
    data = request.get_json()
    
    if not data or not data.get("research_question"):
        abort(400, description="Research question is required")
    
    workflow_state = workflow_service.create_workflow(data["research_question"])
    
    return jsonify({
        "workflow_id": workflow_state.workflow_id,
        "research_question": workflow_state.research_question,
        "state": workflow_state.state
    })

@app.route("/api/workflows/<workflow_id>/execute", methods=["POST"])
def execute_workflow(workflow_id):
    """Execute a workflow."""
    # Load the workflow
    workflow_state = workflow_service.load_workflow_state(workflow_id)
    if not workflow_state:
        abort(404, description=f"Workflow with ID {workflow_id} not found")
    
    # Execute the workflow asynchronously
    # In a production environment, this should be handled by a task queue
    # such as Celery, but for the initial implementation we'll execute synchronously
    updated_state = workflow_service.execute_workflow(workflow_state)
    
    return jsonify({
        "workflow_id": updated_state.workflow_id,
        "research_question": updated_state.research_question,
        "state": updated_state.state,
        "selected_datasets": updated_state.selected_datasets,
        "analyzed_datasets": list(updated_state.analysis_results.keys()),
        "integrated_insights_count": len(updated_state.integrated_insights),
        "errors": updated_state.errors,
        "warnings": updated_state.warnings
    })

@app.route("/api/workflows/<workflow_id>", methods=["GET"])
def get_workflow(workflow_id):
    """Get workflow status."""
    # Load the workflow
    workflow_state = workflow_service.load_workflow_state(workflow_id)
    if not workflow_state:
        abort(404, description=f"Workflow with ID {workflow_id} not found")
    
    return jsonify({
        "workflow_id": workflow_state.workflow_id,
        "research_question": workflow_state.research_question,
        "state": workflow_state.state,
        "search_queries_count": len(workflow_state.search_queries),
        "search_results_count": len(workflow_state.search_results),
        "dataset_relevance_scores_count": len(workflow_state.dataset_relevance_scores),
        "selected_datasets": workflow_state.selected_datasets,
        "fetched_datasets_count": len(workflow_state.fetched_datasets),
        "analysis_requests_count": len(workflow_state.analysis_requests),
        "analysis_results_count": len(workflow_state.analysis_results),
        "integrated_insights_count": len(workflow_state.integrated_insights),
        "errors": workflow_state.errors,
        "warnings": workflow_state.warnings
    })

@app.route("/api/workflows/<workflow_id>/summary", methods=["GET"])
def get_workflow_summary(workflow_id):
    """Get a summary of the workflow results."""
    # Load the workflow
    workflow_state = workflow_service.load_workflow_state(workflow_id)
    if not workflow_state:
        abort(404, description=f"Workflow with ID {workflow_id} not found")
    
    # Generate the summary
    summary = workflow_service.get_workflow_summary(workflow_state)
    
    return jsonify({
        "workflow_id": workflow_state.workflow_id,
        "research_question": workflow_state.research_question,
        "state": workflow_state.state,
        "summary": summary
    })

@app.route("/api/workflows/<workflow_id>/insights", methods=["GET"])
def get_workflow_insights(workflow_id):
    """Get insights from the workflow."""
    # Load the workflow
    workflow_state = workflow_service.load_workflow_state(workflow_id)
    if not workflow_state:
        abort(404, description=f"Workflow with ID {workflow_id} not found")
    
    # Check if the workflow has completed
    if workflow_state.state != "completed":
        return jsonify({
            "workflow_id": workflow_state.workflow_id,
            "state": workflow_state.state,
            "insights": [],
            "message": f"Workflow is not yet completed. Current state: {workflow_state.state}"
        })
    
    # Return the insights
    insights = [insight.dict() for insight in workflow_state.integrated_insights]
    
    return jsonify({
        "workflow_id": workflow_state.workflow_id,
        "research_question": workflow_state.research_question,
        "state": workflow_state.state,
        "insights": insights
    })

@app.route("/api/workflows/<workflow_id>/datasets", methods=["GET"])
def get_workflow_datasets(workflow_id):
    """Get datasets analyzed in the workflow."""
    # Load the workflow
    workflow_state = workflow_service.load_workflow_state(workflow_id)
    if not workflow_state:
        abort(404, description=f"Workflow with ID {workflow_id} not found")
    
    # Prepare dataset information
    datasets = []
    for dataset_id, relevance in workflow_state.dataset_relevance_scores.items():
        dataset_info = {
            "id": dataset_id,
            "relevance_score": relevance.relevance_score,
            "justification": relevance.justification,
            "selected": dataset_id in workflow_state.selected_datasets,
            "fetched": dataset_id in workflow_state.fetched_datasets,
            "analyzed": dataset_id in workflow_state.analysis_results
        }
        datasets.append(dataset_info)
    
    return jsonify({
        "workflow_id": workflow_state.workflow_id,
        "research_question": workflow_state.research_question,
        "state": workflow_state.state,
        "datasets": datasets
    })

@app.errorhandler(400)
def bad_request(error):
    """Handle bad requests."""
    return jsonify({
        "error": "Bad request",
        "message": error.description
    }), 400

@app.errorhandler(404)
def not_found(error):
    """Handle not found errors."""
    return jsonify({
        "error": "Not found",
        "message": error.description
    }), 404

@app.errorhandler(500)
def server_error(error):
    """Handle server errors."""
    logger.error(f"Server error: {str(error)}")
    return jsonify({
        "error": "Server error",
        "message": "An unexpected error occurred"
    }), 500

if __name__ == "__main__":
    # Create cache directory if it doesn't exist
    cache_dir = Path(settings.analysis_cache_dir)
    cache_dir.mkdir(exist_ok=True, parents=True)
    
    # Run the app
    app.run(host="0.0.0.0", port=5000, debug=settings.debug_mode)
