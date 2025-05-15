# UORCA/main_workflow/shared/__init__.py

from pydantic import BaseModel, Field
from typing import Optional, List, Dict, Any

# --- Core context; minimal fields shared everywhere ---
class RNAseqCoreContext(BaseModel):
    accession: str = Field(..., description="NCBI GEO Accession ID for the run")
    output_dir: str = Field(..., description="Output directory for the run")
    fastq_dir: Optional[str] = Field(None, description="Directory for FASTQ files")
    metadata_path: Optional[str] = Field(None, description="Path to sample metadata file")
    organism: str = Field(..., description="Organism name")
    resource_dir: Optional[str] = Field(None, description="Directory containing resources such as Kallisto indices and transcript to gene (t2g) mapping files")
    files: Optional[List[str]] = Field(None, description="List of identified files in the system")
    dataset_information: Optional[str] = Field(None, description="Information about the dataset being identified")

    class Config:
        extra = "allow"  # allows enrichment by agents

# --- Extraction agent context (just uses core here) ---
class ExtractionContext(RNAseqCoreContext):
    pass

# --- Analysis agent context, adds abundance and more ---
class AnalysisContext(RNAseqCoreContext):
    analysis_history: Optional[List[Dict[str, str]]] = Field(default_factory=list,
           description="History of analysis prompts and outputs for reflection")
    abundance_files: Optional[List[str]] = Field(None, description="List of abundance files")
    merged_column: Optional[str] = Field(None, description="Column to merge on")
    unique_groups: Optional[List[str]] = Field(None, description="Unique groups for analysis")
    sample_mapping: Optional[Any] = Field(None, description="Sample mapping for analysis")
    contrast_matrix_df: Optional[Any] = Field(None, description="Contrast matrix DataFrame")
    contrasts: Optional[Any] = Field(None, description="Contrasts for analysis")
    deg_results_path: Optional[str] = Field(None, description="Path to DEGs results")


# --- Reporting context ---
class ReportingContext(RNAseqCoreContext):
    png_dir: Optional[str] = Field(None, description="Directory for PNG files")
    rst_folder: Optional[str] = Field(None, description="Directory for RST files")
    sphinx_output_folder: Optional[str] = Field(None, description="Directory for Sphinx output")
    log_path: Optional[str] = Field(None, description="Path to log file")
