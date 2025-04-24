from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any
import pandas as pd
from pydantic import BaseModel, Field, ConfigDict

@dataclass
class RNAseqData:
    metadata_path: str
    metadata_df: Optional[pd.DataFrame] = None
    merged_column: Optional[str] = None
    unique_groups: Optional[List[str]] = None
    contrast_details: Optional[Dict[str, Any]] = None
    geo_summary: Optional[str] = None
    geo_accession: Optional[str] = None

class ContrastFormat(BaseModel):
    name: str
    expression: str
    description: Optional[str] = Field(default=None)
    justification: Optional[str] = Field(default=None)

class Contrasts(BaseModel):
    contrasts: List[ContrastFormat]
    summary: Optional[str] = None
    model_config = ConfigDict(extra="allow")
