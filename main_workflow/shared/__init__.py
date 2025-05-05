# UORCA/main_workflow/shared/__init__.py

from pydantic import BaseModel, Field
from typing import Optional, List, Dict, Any

# --- Core context; minimal fields shared everywhere ---
class RNAseqCoreContext(BaseModel):
    accession: str
    output_dir: str
    fastq_dir: Optional[str] = None
    metadata_path: Optional[str] = None
    kallisto_index_dir: Optional[str] = None
    organism: str = "human"
    tx2gene_path: Optional[str] = None

    class Config:
        extra = "allow"  # allows enrichment by agents

# --- Extraction agent context (just uses core here) ---
class ExtractionContext(RNAseqCoreContext):
    pass

# --- Analysis agent context, adds abundance and more ---
class AnalysisContext(RNAseqCoreContext):
    abundance_files: Optional[List[str]] = None
    merged_column: Optional[str] = None
    unique_groups: Optional[List[str]] = None
    sample_mapping: Optional[Any] = None
    contrast_matrix_df: Optional[Any] = None
    contrasts: Optional[Any] = None
    deg_results_path: Optional[str] = None

# --- Reporting context ---
class ReportingContext(RNAseqCoreContext):
    png_dir: Optional[str] = None
    rst_folder: Optional[str] = None
    sphinx_output_folder: Optional[str] = None
    log_path: Optional[str] = None
