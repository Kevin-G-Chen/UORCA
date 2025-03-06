"""
Pydantic models for analysis requests and results.
"""
from typing import Dict, List, Optional, Any, Union
from enum import Enum
from pydantic import BaseModel, Field, validator
from .geo_models import OmicsType, GEODatasetWithExpression


class AnalysisType(str, Enum):
    """Types of analysis that can be performed."""
    DIFFERENTIAL_EXPRESSION = "differential_expression"
    PATHWAY_ENRICHMENT = "pathway_enrichment"
    CLUSTERING = "clustering"
    CORRELATION = "correlation"
    CLASSIFICATION = "classification"
    CUSTOM = "custom"


class AnalysisStatus(str, Enum):
    """Status of an analysis job."""
    PENDING = "pending"
    IN_PROGRESS = "in_progress"
    COMPLETED = "completed"
    FAILED = "failed"


class AnalysisParameter(BaseModel):
    """Parameter for an analysis."""
    name: str = Field(..., description="Parameter name")
    value: Any = Field(..., description="Parameter value")
    description: Optional[str] = Field(None, description="Parameter description")
    
    class Config:
        """Pydantic configuration."""
        arbitrary_types_allowed = True


class AnalysisRequest(BaseModel):
    """Request for a dataset analysis."""
    research_question: str = Field(..., description="Research question to address")
    analysis_type: AnalysisType = Field(..., description="Type of analysis to perform")
    dataset_ids: Optional[List[str]] = Field(None, description="Specific dataset IDs to analyze")
    omics_types: List[OmicsType] = Field(default_factory=list, description="Types of omics data to consider")
    organism: Optional[str] = Field(None, description="Organism to filter datasets by")
    parameters: List[AnalysisParameter] = Field(default_factory=list, description="Analysis parameters")
    max_datasets: int = Field(3, description="Maximum number of datasets to analyze")
    
    def has_specific_datasets(self) -> bool:
        """Check if specific datasets were requested."""
        return self.dataset_ids is not None and len(self.dataset_ids) > 0


class DifferentialExpressionResult(BaseModel):
    """Results of differential expression analysis."""
    dataset_id: str = Field(..., description="Dataset ID")
    comparison: str = Field(..., description="Comparison description (e.g., 'Tumor vs Normal')")
    gene_ids: List[str] = Field(..., description="Differentially expressed gene IDs")
    log2_fold_changes: List[float] = Field(..., description="Log2 fold changes")
    p_values: List[float] = Field(..., description="Raw p-values")
    adjusted_p_values: List[float] = Field(..., description="Adjusted p-values (FDR)")
    
    def get_significant_genes(self, adj_p_threshold: float = 0.05, log2fc_threshold: float = 1.0) -> List[str]:
        """Get significantly differentially expressed genes."""
        significant_indices = [
            i for i, (adj_p, log2fc) in enumerate(zip(self.adjusted_p_values, self.log2_fold_changes))
            if adj_p < adj_p_threshold and abs(log2fc) >= log2fc_threshold
        ]
        return [self.gene_ids[i] for i in significant_indices]


class PathwayEnrichmentResult(BaseModel):
    """Results of pathway enrichment analysis."""
    dataset_id: str = Field(..., description="Dataset ID")
    comparison: str = Field(..., description="Comparison description")
    pathway_database: str = Field(..., description="Pathway database used (e.g., 'GO', 'KEGG')")
    pathways: List[str] = Field(..., description="Enriched pathway names")
    p_values: List[float] = Field(..., description="P-values")
    adjusted_p_values: List[float] = Field(..., description="Adjusted p-values (FDR)")
    gene_counts: List[int] = Field(..., description="Number of genes in each pathway")
    gene_lists: List[List[str]] = Field(..., description="Lists of genes in each pathway")
    
    def get_significant_pathways(self, adj_p_threshold: float = 0.05) -> List[str]:
        """Get significantly enriched pathways."""
        significant_indices = [
            i for i, adj_p in enumerate(self.adjusted_p_values)
            if adj_p < adj_p_threshold
        ]
        return [self.pathways[i] for i in significant_indices]


class InsightType(str, Enum):
    """Types of insights that can be derived from the analysis."""
    KEY_GENE = "key_gene"
    PATHWAY = "pathway"
    BIOMARKER = "biomarker"
    DRUG_TARGET = "drug_target"
    DISEASE_MECHANISM = "disease_mechanism"
    CORRELATION = "correlation"
    CONTRADICTION = "contradiction"
    OTHER = "other"


class Insight(BaseModel):
    """An insight derived from analysis of multiple datasets."""
    insight_type: InsightType = Field(..., description="Type of insight")
    description: str = Field(..., description="Description of the insight")
    supporting_datasets: List[str] = Field(..., description="IDs of datasets supporting this insight")
    confidence_score: float = Field(..., description="Confidence score (0-1)")
    evidence: Dict[str, Any] = Field(default_factory=dict, description="Evidence supporting the insight")
    related_entities: List[str] = Field(default_factory=list, description="Related genes, proteins, pathways, etc.")


class AnalysisResult(BaseModel):
    """Results of a complete analysis."""
    request_id: str = Field(..., description="Unique identifier for the analysis request")
    status: AnalysisStatus = Field(AnalysisStatus.COMPLETED, description="Status of the analysis")
    research_question: str = Field(..., description="Original research question")
    datasets_analyzed: List[str] = Field(..., description="IDs of analyzed datasets")
    omics_types: List[OmicsType] = Field(..., description="Types of omics data analyzed")
    
    # Analysis-specific results
    differential_expression_results: Dict[str, DifferentialExpressionResult] = Field(
        default_factory=dict, description="Differential expression results by dataset"
    )
    pathway_enrichment_results: Dict[str, PathwayEnrichmentResult] = Field(
        default_factory=dict, description="Pathway enrichment results by dataset"
    )
    
    # Integrated insights across datasets
    insights: List[Insight] = Field(default_factory=list, description="Insights derived from the analysis")
    
    # LLM-generated summary
    summary: str = Field("", description="Human-readable summary of the analysis")
    
    # Raw data and code
    code_generated: Dict[str, str] = Field(default_factory=dict, description="Code generated for analysis")
    raw_results: Dict[str, Any] = Field(default_factory=dict, description="Raw analysis results")
    
    def get_key_insights(self, min_confidence: float = 0.7, max_count: int = 5) -> List[Insight]:
        """Get the key insights with confidence above the threshold."""
        high_confidence_insights = [
            insight for insight in self.insights
            if insight.confidence_score >= min_confidence
        ]
        # Sort by confidence score in descending order
        sorted_insights = sorted(
            high_confidence_insights,
            key=lambda x: x.confidence_score,
            reverse=True
        )
        return sorted_insights[:max_count]
