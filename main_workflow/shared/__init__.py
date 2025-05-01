from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any
import pandas as pd

@dataclass
class RNAseqData:
    # Core paths
    fastq_dir: str | None = None
    metadata_path: str | None = None
    kallisto_index_dir: str | None = None
    output_dir: str | None = None
    organism: str = "human"

    # Runtime artefacts (populated by agents)
    metadata_df: Optional[pd.DataFrame] = None
    merged_column: Optional[str] = None
    unique_groups: Optional[List[str]] = None
    contrasts: Optional[Any] = None
    geo_summary: Optional[str] = None
