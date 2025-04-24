from pydantic_ai import Agent
from pydantic_ai.mcp import MCPServerHTTP  # Change to HTTP transport
from dataclasses import dataclass
from typing import List, Optional

@dataclass

class MetadataContext:
    metadata_path: str
    metadata_df: Optional[pd.DataFrame] = None
    # This will hold the final analysis column name (or merged version)
    merged_column: Optional[str] = None
    # Store the unique groups found in the analysis column
    unique_groups: Optional[List[str]] = None
    # To store designed contrasts if needed
    contrast_details: Optional[Dict[str, Any]] = None
    # Track all iterations of contrasts for reflection
    all_contrasts_iterations: List[Dict[str, Any]] = field(default_factory=list)
    # Store GEO-related information
    geo_summary: Optional[str] = None
    geo_accession: Optional[str] = None

meta_server = MCPServerStdio(
    command="python",
    args=["-m", "mcp_servers.metadata_server"],   # module path
)

metadata_agent = Agent(
    "openai:gpt-4o",
    deps_type=MetadataContext,
    mcp_servers=[MCPServerHTTP(url="http://127.0.0.1:8000/sse")],
    system_prompt="""
    You are an RNA-seq metadata analysis expert. Use these tools to analyze metadata files:
    - process_metadata: Load and clean the metadata
    - merge_analysis_columns: Combine columns for analysis
    - extract_unique_values: Get unique values from a column
    - fetch_geo_summary: Retrieve GEO database information
    """
)
