"""
Flask routes for GEO API functionality.
"""
from flask import Blueprint, jsonify, request, abort
from typing import Dict, Any, List

from ..services.geo_service import GEOService
from ..models.geo_models import GEOSearchQuery, OmicsType

# Create a blueprint for geo routes
geo_bp = Blueprint('geo', __name__, url_prefix='/api/geo')

# Initialize services
geo_service = GEOService()

@geo_bp.route('/search', methods=['POST'])
def search_datasets():
    """
    Search for datasets in GEO.
    
    Request body:
    {
        "term": "search term",
        "max_results": 10,
        "omics_type": "transcriptomics" (optional),
        "organism": "human" (optional)
    }
    """
    data = request.get_json()
    
    if not data or not data.get('term'):
        abort(400, description="Search term is required")
    
    # Create search query
    query = GEOSearchQuery(
        term=data['term'],
        max_results=int(data.get('max_results', 10)),
        omics_type=OmicsType(data['omics_type']) if data.get('omics_type') else None,
        organism=data.get('organism')
    )
    
    # Perform search
    try:
        result = geo_service.search_datasets(query)
        
        # Convert to dict for JSON response
        response = {
            'query': {
                'term': query.term,
                'max_results': query.max_results,
                'omics_type': query.omics_type.value if query.omics_type else None,
                'organism': query.organism
            },
            'total_count': result.total_count,
            'dataset_ids': result.dataset_ids,
            'series_ids': result.series_ids
        }
        
        return jsonify(response)
        
    except Exception as e:
        abort(500, description=f"Error searching datasets: {str(e)}")

@geo_bp.route('/dataset/<dataset_id>', methods=['GET'])
def get_dataset(dataset_id):
    """Get metadata for a GEO dataset."""
    try:
        if dataset_id.startswith('GDS'):
            dataset = geo_service.get_dataset(dataset_id)
            return jsonify(dataset.dict())
        elif dataset_id.startswith('GSE'):
            series = geo_service.get_series(dataset_id)
            return jsonify(series.dict())
        else:
            abort(400, description=f"Invalid dataset ID format: {dataset_id}")
    except Exception as e:
        abort(500, description=f"Error fetching dataset: {str(e)}")

@geo_bp.route('/dataset/<dataset_id>/expression', methods=['GET'])
def get_expression_data(dataset_id):
    """Get expression data for a GEO dataset."""
    try:
        if not dataset_id.startswith('GDS'):
            abort(400, description="Expression data retrieval is only supported for GDS datasets")
            
        expression = geo_service.get_expression_data(dataset_id)
        
        if not expression:
            return jsonify({
                'dataset_id': dataset_id,
                'error': 'No expression data available'
            })
        
        # For large datasets, we might want to limit the response size
        limit = request.args.get('limit', type=int)
        
        response = {
            'dataset_id': dataset_id,
            'gene_count': len(expression.gene_ids),
            'sample_count': len(expression.sample_ids),
            'sample_ids': expression.sample_ids
        }
        
        if limit:
            response['gene_ids'] = expression.gene_ids[:limit]
            response['values'] = expression.values[:limit]
        else:
            response['gene_ids'] = expression.gene_ids
            response['values'] = expression.values
            
        return jsonify(response)
        
    except Exception as e:
        abort(500, description=f"Error fetching expression data: {str(e)}")

@geo_bp.route('/sample/<sample_id>', methods=['GET'])
def get_sample(sample_id):
    """Get metadata for a GEO sample."""
    try:
        if not sample_id.startswith('GSM'):
            abort(400, description=f"Invalid sample ID format: {sample_id}")
            
        sample = geo_service.get_sample(sample_id)
        return jsonify(sample.dict())
        
    except Exception as e:
        abort(500, description=f"Error fetching sample: {str(e)}")

def register_blueprint(app):
    """Register the blueprint with the Flask app."""
    app.register_blueprint(geo_bp)
