"""
Service for performing analysis on omics datasets.
"""
import os
import logging
import tempfile
import pandas as pd
import numpy as np
from typing import Dict, List, Any, Optional, Tuple, Union
import json
import uuid
import importlib.util
import sys
from pathlib import Path

from ..config import get_settings
from ..models.geo_models import GEODatasetWithExpression
from ..models.analysis_models import (
    AnalysisRequest, 
    AnalysisResult, 
    AnalysisStatus,
    AnalysisParameter,
    AnalysisType,
    DifferentialExpressionResult,
    PathwayEnrichmentResult,
    Insight
)
from ..services.llm_service import LLMService

logger = logging.getLogger(__name__)

class AnalysisService:
    """Service for performing analysis on omics datasets."""
    
    def __init__(self):
        """Initialize the analysis service."""
        self.settings = get_settings()
        self.llm_service = LLMService()
        
        # Create cache directory if it doesn't exist
        cache_dir = Path(self.settings.analysis_cache_dir)
        cache_dir.mkdir(exist_ok=True, parents=True)
    
    def analyze_dataset(
        self, 
        dataset: GEODatasetWithExpression, 
        request: AnalysisRequest
    ) -> AnalysisResult:
        """
        Analyze a dataset based on the analysis request.
        
        Args:
            dataset: The dataset with expression data
            request: The analysis request
            
        Returns:
            The analysis result
        """
        # Generate a unique ID for this analysis
        request_id = f"analysis_{uuid.uuid4().hex[:8]}"
        
        # Log the start of analysis
        logger.info(f"Starting analysis {request_id} for dataset {dataset.metadata.id}")
        
        try:
            # Design the analysis using LLM
            analysis_design = self.llm_service.design_analysis(
                research_question=request.research_question,
                dataset=dataset
            )
            
            # Generate the analysis code
            analysis_code = self.llm_service.generate_analysis_code(
                research_question=request.research_question,
                dataset=dataset,
                analysis_design=analysis_design
            )
            
            # Execute the analysis code
            raw_results = self._execute_analysis_code(
                code=analysis_code,
                dataset=dataset
            )
            
            # Parse the raw results into structured results
            structured_results = self._parse_analysis_results(
                raw_results=raw_results,
                analysis_type=AnalysisType(analysis_design.get("analysis_type", "custom")),
                dataset_id=dataset.metadata.id
            )
            
            # Create the analysis result
            result = AnalysisResult(
                request_id=request_id,
                status=AnalysisStatus.COMPLETED,
                research_question=request.research_question,
                datasets_analyzed=[dataset.metadata.id],
                omics_types=[dataset.metadata.omics_type],
                code_generated={dataset.metadata.id: analysis_code},
                raw_results={dataset.metadata.id: raw_results}
            )
            
            # Add the structured results
            if analysis_design.get("analysis_type") == "differential_expression" and "differential_expression_result" in structured_results:
                result.differential_expression_results[dataset.metadata.id] = structured_results["differential_expression_result"]
            
            if analysis_design.get("analysis_type") == "pathway_enrichment" and "pathway_enrichment_result" in structured_results:
                result.pathway_enrichment_results[dataset.metadata.id] = structured_results["pathway_enrichment_result"]
            
            # Cache the results
            self._cache_results(request_id, result)
            
            logger.info(f"Completed analysis {request_id} for dataset {dataset.metadata.id}")
            return result
            
        except Exception as e:
            logger.error(f"Error analyzing dataset {dataset.metadata.id}: {str(e)}")
            
            # Return a failed result
            return AnalysisResult(
                request_id=request_id,
                status=AnalysisStatus.FAILED,
                research_question=request.research_question,
                datasets_analyzed=[dataset.metadata.id],
                omics_types=[dataset.metadata.omics_type],
                code_generated={},
                raw_results={"error": str(e)}
            )
    
    def integrate_results(
        self, 
        research_question: str, 
        analysis_results: List[AnalysisResult]
    ) -> List[Insight]:
        """
        Integrate results from multiple dataset analyses.
        
        Args:
            research_question: The research question
            analysis_results: List of analysis results
            
        Returns:
            List of integrated insights
        """
        # Prepare simplified results for the LLM
        simplified_results = []
        
        for result in analysis_results:
            # Only include completed analyses
            if result.status != AnalysisStatus.COMPLETED:
                continue
                
            dataset_id = result.datasets_analyzed[0] if result.datasets_analyzed else "unknown"
            
            # Extract significant findings
            significant_findings = []
            
            # Check for differential expression results
            if dataset_id in result.differential_expression_results:
                de_result = result.differential_expression_results[dataset_id]
                sig_genes = de_result.get_significant_genes()
                if sig_genes:
                    significant_findings.append({
                        "type": "differentially_expressed_genes",
                        "count": len(sig_genes),
                        "top_genes": sig_genes[:10]  # Limit to top 10
                    })
            
            # Check for pathway enrichment results
            if dataset_id in result.pathway_enrichment_results:
                pe_result = result.pathway_enrichment_results[dataset_id]
                sig_pathways = pe_result.get_significant_pathways()
                if sig_pathways:
                    significant_findings.append({
                        "type": "enriched_pathways",
                        "count": len(sig_pathways),
                        "top_pathways": sig_pathways[:10]  # Limit to top 10
                    })
            
            # Add to simplified results
            simplified_results.append({
                "dataset_id": dataset_id,
                "dataset_title": "Dataset " + dataset_id,  # Simplified for the example
                "analysis_type": next(iter(result.code_generated.keys()), "unknown"),
                "significant_findings": significant_findings
            })
        
        # Use LLM to integrate insights
        insights_data = self.llm_service.integrate_insights(
            research_question=research_question,
            analysis_results=simplified_results
        )
        
        # Convert to Insight objects
        insights = []
        for insight_data in insights_data:
            try:
                insight = Insight(
                    insight_type=insight_data.get("insight_type", "other"),
                    description=insight_data.get("description", ""),
                    supporting_datasets=insight_data.get("supporting_datasets", []),
                    confidence_score=float(insight_data.get("confidence_score", 0.5)),
                    related_entities=insight_data.get("related_entities", [])
                )
                insights.append(insight)
            except Exception as e:
                logger.error(f"Error creating insight: {str(e)}")
        
        return insights
    
    def generate_summary(
        self, 
        research_question: str, 
        datasets_analyzed: List[str], 
        insights: List[Insight]
    ) -> str:
        """
        Generate a human-readable summary of the analysis.
        
        Args:
            research_question: The research question
            datasets_analyzed: List of analyzed dataset IDs
            insights: List of integrated insights
            
        Returns:
            Formatted summary
        """
        # Convert insights to dictionaries
        insights_dict = [insight.dict() for insight in insights]
        
        # Use LLM to generate summary
        summary = self.llm_service.generate_summary_report(
            research_question=research_question,
            datasets_analyzed=datasets_analyzed,
            insights=insights_dict
        )
        
        return summary
    
    def _execute_analysis_code(
        self, 
        code: str, 
        dataset: GEODatasetWithExpression
    ) -> Dict[str, Any]:
        """
        Execute the generated analysis code.
        
        Args:
            code: The Python code to execute
            dataset: The dataset to analyze
            
        Returns:
            Dictionary with the analysis results
        """
        # For safety in a production environment, this should be executed in a sandbox
        # Here we use a simple approach with a temporary file and exec
        
        # Create a temporary directory for the analysis
        with tempfile.TemporaryDirectory() as temp_dir:
            # Save the dataset to a CSV file
            if dataset.expression_matrix:
                df = pd.DataFrame(
                    data=dataset.expression_matrix.values,
                    index=dataset.expression_matrix.gene_ids,
                    columns=dataset.expression_matrix.sample_ids
                )
                expression_file = os.path.join(temp_dir, "expression_matrix.csv")
                df.to_csv(expression_file)
                
                # Save sample metadata
                sample_meta = pd.DataFrame.from_dict(dataset.sample_metadata, orient='index')
                sample_meta_file = os.path.join(temp_dir, "sample_metadata.csv")
                sample_meta.to_csv(sample_meta_file)
                
                # Save the analysis code
                code_file = os.path.join(temp_dir, "analysis_code.py")
                with open(code_file, 'w') as f:
                    f.write(code)
                
                # Create a module to hold the results
                module_name = f"analysis_module_{uuid.uuid4().hex[:8]}"
                spec = importlib.util.spec_from_file_location(module_name, code_file)
                module = importlib.util.module_from_spec(spec)
                
                # Add dataset paths to the module
                module.expression_file = expression_file
                module.sample_meta_file = sample_meta_file
                module.results = {}
                
                try:
                    # Execute the module
                    spec.loader.exec_module(module)
                    
                    # Extract results
                    return getattr(module, "results", {})
                    
                except Exception as e:
                    logger.error(f"Error executing analysis code: {str(e)}")
                    return {"error": str(e)}
            else:
                return {"error": "No expression matrix available for analysis"}
    
    def _parse_analysis_results(
        self, 
        raw_results: Dict[str, Any], 
        analysis_type: AnalysisType, 
        dataset_id: str
    ) -> Dict[str, Any]:
        """
        Parse raw analysis results into structured results.
        
        Args:
            raw_results: The raw analysis results
            analysis_type: The type of analysis performed
            dataset_id: The dataset ID
            
        Returns:
            Dictionary with structured results
        """
        structured_results = {}
        
        # Check for errors
        if "error" in raw_results:
            return {"error": raw_results["error"]}
        
        # Parse differential expression results
        if analysis_type == AnalysisType.DIFFERENTIAL_EXPRESSION and "de_results" in raw_results:
            de_data = raw_results["de_results"]
            
            try:
                # Create a differential expression result object
                de_result = DifferentialExpressionResult(
                    dataset_id=dataset_id,
                    comparison=de_data.get("comparison", "Unknown comparison"),
                    gene_ids=de_data.get("gene_ids", []),
                    log2_fold_changes=de_data.get("log2_fold_changes", []),
                    p_values=de_data.get("p_values", []),
                    adjusted_p_values=de_data.get("adjusted_p_values", [])
                )
                structured_results["differential_expression_result"] = de_result
            except Exception as e:
                logger.error(f"Error parsing differential expression results: {str(e)}")
        
        # Parse pathway enrichment results
        if analysis_type == AnalysisType.PATHWAY_ENRICHMENT and "enrichment_results" in raw_results:
            pe_data = raw_results["enrichment_results"]
            
            try:
                # Create a pathway enrichment result object
                pe_result = PathwayEnrichmentResult(
                    dataset_id=dataset_id,
                    comparison=pe_data.get("comparison", "Unknown comparison"),
                    pathway_database=pe_data.get("pathway_database", "Unknown database"),
                    pathways=pe_data.get("pathways", []),
                    p_values=pe_data.get("p_values", []),
                    adjusted_p_values=pe_data.get("adjusted_p_values", []),
                    gene_counts=pe_data.get("gene_counts", []),
                    gene_lists=pe_data.get("gene_lists", [])
                )
                structured_results["pathway_enrichment_result"] = pe_result
            except Exception as e:
                logger.error(f"Error parsing pathway enrichment results: {str(e)}")
        
        return structured_results
    
    def _cache_results(self, request_id: str, result: AnalysisResult) -> None:
        """
        Cache analysis results for future use.
        
        Args:
            request_id: The analysis request ID
            result: The analysis result
        """
        try:
            # Convert to JSON-serializable format
            result_dict = result.dict()
            
            # Save to cache directory
            cache_file = Path(self.settings.analysis_cache_dir) / f"{request_id}.json"
            with open(cache_file, 'w') as f:
                json.dump(result_dict, f, indent=2)
                
            logger.debug(f"Cached analysis results to {cache_file}")
        except Exception as e:
            logger.error(f"Error caching analysis results: {str(e)}")
