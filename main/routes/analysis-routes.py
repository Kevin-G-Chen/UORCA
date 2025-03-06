"""
Flask routes for analysis functionality.
"""
from flask import Blueprint, jsonify, request, abort
from typing import Dict, Any, List

from ..services.analysis_service import AnalysisService
from ..services.geo_service import GEOService
from ..models.analysis_models import AnalysisRequest, AnalysisType
from ..models.geo_models import OmicsType

# Create a blueprint for analysis routes
analysis_bp = Blueprint('analysis', __name__, url_prefix='/api/analysis')

# Initialize services
analysis_service = AnalysisService()
geo_service = GEOService()

@analysis_bp.route('/dataset/<dataset_id>', methods=['POST'])
def analyze_dataset(dataset_id):
    """
    Analyze a GEO dataset.
    
    Request body:
    {
        "research_question": "What genes are differentially expressed in...",
        "analysis_type": "differential_expression" (optional),
        "parameters": [
            {
                "name": "param_name",
                "value": "param_value",
                "description": "Parameter description"
            }
        ]
    }
    """
    data = request.get_json()
    
    if not data or not data.get('research_question'):
        abort(400, description="Research question is required")
    
    # Fetch the dataset
    try:
        if dataset_id.startswith('GDS'):
            dataset = geo_service.get_dataset_with_expression(dataset_id)
        elif dataset_id.startswith('GSE'):
            dataset = geo_service.get_series_with_expression(dataset_id)
        else:
            abort(400, description=f"Invalid dataset ID format: {dataset_id}")
            
        if not dataset.expression_matrix:
            abort(400, description=f"Dataset {dataset_id} has no expression data available")
    except Exception as e:
        abort(500, description=f"Error fetching dataset: {str(e)}")
    
    # Create analysis request
    analysis_type = data.get('analysis_type', 'custom')
    request_obj = AnalysisRequest(
        research_question=data['research_question'],
        analysis_type=AnalysisType(analysis_type) if analysis_type else AnalysisType.CUSTOM,
        dataset_ids=[dataset_id],
        omics_types=[dataset.metadata.omics_type],
        parameters=data.get('parameters', []),
        max_datasets=1
    )
    
    # Perform analysis
    try:
        result = analysis_service.analyze_dataset(dataset, request_obj)
        
        # Convert to dict for JSON response
        response = {
            'request_id': result.request_id,
            'status': result.status,
            'research_question': result.research_question,
            'datasets_analyzed': result.datasets_analyzed,
            'omics_types': [omics.value for omics in result.omics_types],
            'insights_count': len(result.insights)
        }
        
        return jsonify(response)
        
    except Exception as e:
        abort(500, description=f"Error analyzing dataset: {str(e)}")

@analysis_bp.route('/result/<request_id>', methods=['GET'])
def get_analysis_result(request_id):
    """Get the results of an analysis."""
    try:
        # Load result from cache (would normally be fetched from a database)
        result = analysis_service._load_result(request_id)
        
        if not result:
            abort(404, description=f"Analysis result with ID {request_id} not found")
            
        return jsonify(result.dict())
        
    except Exception as e:
        abort(500, description=f"Error fetching analysis result: {str(e)}")

@analysis_bp.route('/integrate', methods=['POST'])
def integrate_results():
    """
    Integrate results from multiple analyses.
    
    Request body:
    {
        "research_question": "What genes are differentially expressed in...",
        "analysis_ids": ["analysis_id1", "analysis_id2", ...]
    }
    """
    data = request.get_json()
    
    if not data or not data.get('research_question') or not data.get('analysis_ids'):
        abort(400, description="Research question and analysis IDs are required")
    
    # Load analysis results
    analysis_results = []
    for analysis_id in data['analysis_ids']:
        result = analysis_service._load_result(analysis_id)
        if result:
            analysis_results.append(result)
    
    if not analysis_results:
        abort(404, description="No valid analysis results found")
    
    # Integrate results
    try:
        insights = analysis_service.integrate_results(
            research_question=data['research_question'],
            analysis_results=analysis_results
        )
        
        # Generate summary
        summary = analysis_service.generate_summary(
            research_question=data['research_question'],
            datasets_analyzed=[id for result in analysis_results for id in result.datasets_analyzed],
            insights=insights
        )
        
        # Convert insights to dicts for JSON response
        insights_dicts = [insight.dict() for insight in insights]
        
        response = {
            'research_question': data['research_question'],
            'analysis_ids': data['analysis_ids'],
            'insights': insights_dicts,
            'summary': summary
        }
        
        return jsonify(response)
        
    except Exception as e:
        abort(500, description=f"Error integrating results: {str(e)}")

def register_blueprint(app):
    """Register the blueprint with the Flask app."""
    app.register_blueprint(analysis_bp)
