"""
Pydantic models for GEO (Gene Expression Omnibus) data structures.
"""
from typing import Dict, List, Optional, Any, Union
from datetime import datetime
from enum import Enum
from pydantic import BaseModel, Field, validator


class OmicsType(str, Enum):
    """Types of omics data."""
    TRANSCRIPTOMICS = "transcriptomics"
    PROTEOMICS = "proteomics"
    GENOMICS = "genomics"
    METABOLOMICS = "metabolomics"
    EPIGENOMICS = "epigenomics"
    UNKNOWN = "unknown"


class GEOPlatform(BaseModel):
    """GEO Platform (GPL) information."""
    id: str = Field(..., description="GEO platform ID (GPL)")
    title: str = Field(..., description="Platform title")
    technology: str = Field(..., description="Technology type")
    organism: Optional[str] = Field(None, description="Organism")
    
    @validator('id')
    def validate_gpl_id(cls, v):
        """Ensure the platform ID starts with GPL."""
        if not v.startswith("GPL"):
            raise ValueError(f"Platform ID must start with GPL, got {v}")
        return v


class GEOSample(BaseModel):
    """GEO Sample (GSM) information."""
    id: str = Field(..., description="GEO sample ID (GSM)")
    title: str = Field(..., description="Sample title")
    attributes: Dict[str, str] = Field(default_factory=dict, description="Sample attributes")
    
    @validator('id')
    def validate_gsm_id(cls, v):
        """Ensure the sample ID starts with GSM."""
        if not v.startswith("GSM"):
            raise ValueError(f"Sample ID must start with GSM, got {v}")
        return v


class GEOSeries(BaseModel):
    """GEO Series (GSE) information."""
    id: str = Field(..., description="GEO series ID (GSE)")
    title: str = Field(..., description="Series title")
    summary: str = Field(..., description="Series summary")
    organism: Optional[str] = Field(None, description="Organism")
    platforms: List[GEOPlatform] = Field(default_factory=list, description="Platforms used")
    samples: List[str] = Field(default_factory=list, description="Sample IDs")
    publication_date: Optional[datetime] = Field(None, description="Publication date")
    submission_date: Optional[datetime] = Field(None, description="Submission date")
    last_update_date: Optional[datetime] = Field(None, description="Last update date")
    omics_type: OmicsType = Field(OmicsType.UNKNOWN, description="Type of omics data")
    
    @validator('id')
    def validate_gse_id(cls, v):
        """Ensure the series ID starts with GSE."""
        if not v.startswith("GSE"):
            raise ValueError(f"Series ID must start with GSE, got {v}")
        return v


class GEODataset(BaseModel):
    """GEO Dataset (GDS) information."""
    id: str = Field(..., description="GEO dataset ID (GDS)")
    title: str = Field(..., description="Dataset title")
    summary: str = Field(..., description="Dataset summary")
    organism: Optional[str] = Field(None, description="Organism")
    platform: GEOPlatform = Field(..., description="Platform used")
    samples: List[str] = Field(default_factory=list, description="Sample IDs")
    publication_date: Optional[datetime] = Field(None, description="Publication date")
    submission_date: Optional[datetime] = Field(None, description="Submission date")
    last_update_date: Optional[datetime] = Field(None, description="Last update date")
    omics_type: OmicsType = Field(OmicsType.UNKNOWN, description="Type of omics data")
    
    @validator('id')
    def validate_gds_id(cls, v):
        """Ensure the dataset ID starts with GDS."""
        if not v.startswith("GDS"):
            raise ValueError(f"Dataset ID must start with GDS, got {v}")
        return v


class GEOSearchQuery(BaseModel):
    """Structure for GEO search query parameters."""
    term: str = Field(..., description="Search terms")
    max_results: int = Field(10, description="Maximum number of results to return")
    omics_type: Optional[OmicsType] = Field(None, description="Filter by omics type")
    organism: Optional[str] = Field(None, description="Filter by organism")
    
    def to_query_params(self) -> Dict[str, str]:
        """Convert to query parameters for GEO API."""
        params = {
            "term": self.term,
            "retmax": str(self.max_results)
        }
        
        # Add filters if specified
        if self.omics_type:
            params["term"] += f" AND {self.omics_type.value}[DataType]"
        if self.organism:
            params["term"] += f" AND {self.organism}[Organism]"
            
        return params


class GEOSearchResult(BaseModel):
    """Results from a GEO search query."""
    query: GEOSearchQuery = Field(..., description="Original search query")
    total_count: int = Field(0, description="Total number of matching records")
    dataset_ids: List[str] = Field(default_factory=list, description="List of matching dataset IDs")
    series_ids: List[str] = Field(default_factory=list, description="List of matching series IDs")
    
    def is_empty(self) -> bool:
        """Check if the search returned any results."""
        return len(self.dataset_ids) == 0 and len(self.series_ids) == 0


# Composite model for storing a dataset with its metadata and expression matrix
class ExpressionMatrix(BaseModel):
    """Expression matrix with gene/probe IDs and sample values."""
    gene_ids: List[str] = Field(..., description="Gene or probe IDs")
    sample_ids: List[str] = Field(..., description="Sample IDs")
    values: List[List[float]] = Field(..., description="Expression values matrix")
    
    def validate_dimensions(self) -> bool:
        """Validate that the matrix dimensions match the gene and sample counts."""
        if not self.values:
            return False
        
        gene_count = len(self.gene_ids)
        sample_count = len(self.sample_ids)
        
        # Check if each row (gene) has the correct number of samples
        return all(len(row) == sample_count for row in self.values) and len(self.values) == gene_count


class GEODatasetWithExpression(BaseModel):
    """GEO dataset with its expression matrix."""
    metadata: Union[GEODataset, GEOSeries] = Field(..., description="Dataset metadata")
    expression_matrix: Optional[ExpressionMatrix] = Field(None, description="Expression matrix if available")
    sample_metadata: Dict[str, Dict[str, Any]] = Field(default_factory=dict, description="Metadata for each sample")
