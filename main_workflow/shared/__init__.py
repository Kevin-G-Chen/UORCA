from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any
import pandas as pd
from pydantic import BaseModel, Field, ConfigDict

from typing import Optional, List, Dict, Any
import pandas as pd
from pydantic import BaseModel, Field, ConfigDict

class RNAseqData(BaseModel):
    fastq_dir: Optional[str] = Field(
        None,
        description="Directory containing the raw FASTQ files for Kallisto"
    )
    metadata_path: Optional[str] = Field(
        None,
        description="Path to the sample-metadata table (CSV/TSV)"
    )
    kallisto_index_dir: Optional[str] = Field(
        None,
        description="Directory where Kallisto indices are stored or created"
    )
    output_dir: Optional[str] = Field(
        None,
        description="Root output directory for all analysis results"
    )
    organism: str = Field(
        "human",
        description="Organism keyword for reference selection; default is 'human'"
    )
    tx2gene_path: Optional[str] = Field(
        None,
        description="Path to a transcript-to-gene mapping TSV file"
    )
    kallisto_index_path: Optional[str] = Field(
        None,
        description="Optional path to an existing Kallisto index"
    )

    # Run-time and analysis artefacts
    metadata_df: Optional[Any] = Field(
        None,
        description="Loaded metadata (as pandas.DataFrame); populated after extraction"
    )
    abundance_files: Optional[List[str]] = Field(
        None,
        description="List of per-sample abundance TSVs produced by Kallisto"
    )
    merged_column: Optional[str] = Field(
        None,
        description="Column name for the analysis grouping variable (after merge)"
    )
    unique_groups: Optional[List[str]] = Field(
        None,
        description="Unique group values in the merged analysis column"
    )
    contrast_groups: Optional[Dict[str, Dict[str, str]]] = Field(
        None,
        description="Mapping of contrast names to group pairs"
    )
    sample_mapping: Optional[Any] = Field(
        None,
        description="Sample mapping DataFrame linking abundance files to metadata"
    )
    contrasts: Optional[Any] = Field(
        None,
        description="Contrast definitions (agent-generated; format can vary)"
    )
    deg_results_path: Optional[str] = Field(
        None,
        description="Path to the main differential expression results CSV"
    )
    contrast_path: Optional[str] = Field(
        None,
        description="Path to contrasts CSV file"
    )
    contrast_matrix_df: Optional[Any] = Field(
        None,
        description="DataFrame format of contrast matrix"
    )
    geo_summary: Optional[str] = Field(
        None,
        description="Plain text summary of the GEO/AE dataset"
    )
    model_config = ConfigDict(extra="allow")  # Enable pydantic extra fields if needed

class ReportContext(BaseModel):
    png_dir: str = Field(..., description="Directory where the PNG files are located.")
    rst_folder: str = Field(..., description="Directory to store the generated RST file and images.")
    sphinx_output_folder: str = Field(..., description="Directory for the Sphinx project and HTML output.")
    log_path: str = Field(..., description="Path to save the Sphinx build log file.")
    model_config = ConfigDict(extra="allow")
