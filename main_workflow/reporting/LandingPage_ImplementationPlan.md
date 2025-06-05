Comprehensive Implementation Plan for Agent-Based Reporting System

## 1. Architecture Overview

```
┌─────────────────┐      ┌───────────────────┐      ┌─────────────────┐
│                 │      │                   │      │                 │
│  Streamlit UI   │<─────│  Agent Controller │<─────│   MCP Servers   │
│                 │      │                   │      │                 │
└─────────────────┘      └───────────────────┘      └─────────────────┘
       ^                         ^                         ^
       │                         │                         │
       v                         v                         v
┌─────────────────┐      ┌───────────────────┐      ┌─────────────────┐
│  Auto-Startup   │      │    LLM Model      │      │  Data Sources   │
│    Handler      │      │  (GPT/Llama)      │      │  (CSV Files)    │
└─────────────────┘      └───────────────────┘      └─────────────────┘
```

## 2. Component Design

### 2.1. MCP Server Implementation

Create dedicated MCP servers for RNA-seq data analysis:

1. **Data Extraction Server (`mcp_data_extractor.py`)**
   - Focus on raw data extraction operations
   - Low-level access to CSV files and datasets
   - Simple data transformation operations

2. **Analysis Server (`mcp_analysis.py`)**
   - Higher-level analysis of extracted data
   - Statistical operations across datasets
   - Pattern recognition in gene expression data

### 2.2. Tool Design for MCP Servers

#### Data Extraction Tools:

```python
@mcp.tool()
async def get_deg_counts(contrast_id: str, threshold: float = 0.05, lfc_threshold: float = 1.0) -> str:
    """Get the number of differentially expressed genes for a specific contrast."""
    # Implementation...
    return json.dumps({"count": count, "up_regulated": up_count, "down_regulated": down_count})

@mcp.tool()
async def get_gene_info(gene_id: str, contrasts: List[str] = None) -> str:
    """Get expression information for a specific gene across selected contrasts."""
    # Implementation...
    return json.dumps({"contrasts": result_data})

@mcp.tool()
async def find_common_degs(contrasts: List[str], min_frequency: int = 2,
                          p_value_threshold: float = 0.05, lfc_threshold: float = 1.0) -> str:
    """Find genes that are differentially expressed across multiple contrasts."""
    # Implementation...
    return json.dumps({"genes": common_genes})

@mcp.tool()
async def get_contrast_metadata(contrast_id: str) -> str:
    """Get metadata information about a specific contrast."""
    # Implementation...
    return json.dumps({"description": description, "formula": formula})
```

#### Analysis Tools:

```python
@mcp.tool()
async def find_enriched_patterns(genes: List[str], method: str = "clustering") -> str:
    """Find expression patterns across a set of genes."""
    # Implementation...
    return json.dumps({"patterns": patterns})

@mcp.tool()
async def rank_contrasts_by_relevance(research_question: str, max_contrasts: int = 10) -> str:
    """Rank contrasts by relevance to a specified research question."""
    # Implementation...
    return json.dumps({"ranked_contrasts": ranked_contrasts})

@mcp.tool()
async def generate_hypothesis(genes: List[str], contrasts: List[str]) -> str:
    """Generate hypotheses based on gene expression patterns."""
    # Implementation...
    return json.dumps({"hypotheses": hypotheses})
```

### 2.3. Agent Controller

Develop an agent controller that integrates MCP servers with the LLM:

```python
def setup_reporting_agent(research_question: str, results_dir: str):
    """Set up the reporting agent with MCP servers and LLM."""
    # Initialize MCP servers with appropriate environment variables
    data_server = setup_mcp_server("data_extractor",
                                  env_vars={"RESULTS_DIR": results_dir})
    analysis_server = setup_mcp_server("analysis",
                                      env_vars={"RESULTS_DIR": results_dir})

    # Configure the agent with system prompt and tools
    agent = Agent(
        model="openai:gpt-4.1-mini",  # Or local model via Ollama
        model_settings={"temperature": 0.1},
        mcp_servers=[data_server, analysis_server],
        system_prompt=(
            f"You are an expert computational biologist analyzing RNA-seq results. "
            f"Your task is to identify key findings related to: {research_question}. "
            f"Use the available tools to extract data, perform analyses, and generate "
            f"a comprehensive report with visualizations."
        ),
    )

    return agent
```

### 2.4. Auto-Start Integration

```python
# Modified main Streamlit flow in uorca_explorer.py

@st.cache_resource
def get_integrator_and_agent(path):
    """Initialize both the integrator and reporting agent."""
    try:
        # Initialize the data integrator
        ri = ResultsIntegrator(results_dir=path)
        ri.load_data()

        # Initialize the reporting agent with generic research question
        agent = ReportingAgent(
            results_dir=path,
            model_name=os.getenv("OPENAI_MODEL", "openai:gpt-4.1-mini")
        )
        agent.setup_servers()

        # Return both objects
        return ri, agent, None
    except Exception as e:
        return None, None, str(e)

# Load data and initialize agent on startup
with st.sidebar.status("Loading data...", expanded=True) as status:
    ri, agent, error = get_integrator_and_agent(results_dir)

    if error:
        status.update(label=f"Error loading data: {error}", state="error")
    elif not ri or not ri.cpm_data:
        status.update(label="No data found. Please check the directory path.", state="error")
    else:
        # Success - data loaded successfully
        status.update(label=f"✅ Loaded {len(ri.cpm_data)} datasets", state="complete")

        # Set flag to indicate first-time load
        if 'initial_load_complete' not in st.session_state:
            st.session_state.initial_load_complete = True
            st.session_state.trigger_initial_analysis = True
```

## 3. Implementation Steps

### Phase 1: Core MCP Server Development (2-3 weeks)

1. **Set up basic MCP server structure**
   - Create base server scripts for data extraction and analysis
   - Define tool interfaces (input/output schemas)
   - Implement connection to data sources

2. **Implement core data extraction tools**
   - DEG count extraction
   - Gene information retrieval
   - Common DEG identification
   - Dataset/contrast metadata access

3. **Test individual tools**
   - Create unit tests for each tool function
   - Verify correct data extraction from test datasets
   - Test edge cases (empty datasets, missing files, etc.)

### Phase 2: Agent Integration (1-2 weeks)

1. **Develop agent controller**
   - Create wrapper for agent initialization
   - Implement system prompt templates
   - Set up tool routing and error handling

2. **LLM integration**
   - Test with different models (GPT-4, Llama)
   - Optimize prompts for RNA-seq analysis
   - Add guardrails for inference failures

3. **Build feedback mechanisms**
   - Error recovery strategies
   - Logging of tool usage patterns
   - Performance monitoring

### Phase 3: Auto-Start Integration (1 week)

1. **Implement automatic agent invocation**
   - Design startup flow to detect valid data loading
   - Create default analysis settings for initial run
   - Implement session state management to track analysis state

2. **Develop progressive disclosure UI**
   - Show initial analysis results immediately
   - Allow refinement with specific research questions
   - Add visual indicators of analysis status

3. **Add error recovery mechanisms**
   - Graceful handling of initialization failures
   - Allow manual analysis if automatic analysis fails
   - Cache partial results for faster recovery

### Phase 4: UI Integration (2 weeks)

1. **Replace existing AI landing page**
   - Ensure that the default view shows auto-generated results
   - Make the refinement controls clearly visible but secondary
   - Integrate with session state to persist analysis across tabs

2. **Add interactive elements**
   - Implement expandable/collapsible result sections
   - Add filters to focus on specific aspects of the analysis
   - Create options to regenerate or refine analyses

3. **Enhance visualization**
   - Auto-select the most relevant visualizations based on analysis
   - Create dynamic visualization generation based on findings
   - Add contextual explanations alongside visualizations

### Phase 5: Testing and Refinement (2 weeks)

1. **End-to-end testing**
   - Specifically test the automatic startup flow
   - Verify behavior with various data directories
   - Test recovery from interruptions and errors

2. **Performance optimization**
   - Optimize for speed during initial data loading
   - Implement tiered analysis (quick summary first, then detailed)
   - Improve caching for repeated analyses

3. **User experience testing**
   - Gather feedback on automatic analysis quality
   - Test with various research domains and data types
   - Refine prompting strategy for initial analysis

## 4. Code Samples for Implementation

### 4.1. Data Extraction MCP Server

```python
#!/usr/bin/env python3
# mcp_data_extractor.py

import os
import sys
import json
import pandas as pd
from typing import List, Dict, Optional, Tuple
from pathlib import Path
import logging

from mcp.server.fastmcp import FastMCP

# Initialize MCP server
mcp = FastMCP("data_extractor")
log = lambda m: print(f"MCP DATA: {m}", file=sys.stderr, flush=True)

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    stream=sys.stderr
)
logger = logging.getLogger("mcp-data-extractor")

# Global variables
RESULTS_DIR = os.environ.get("RESULTS_DIR", "")
_cache = {}  # Simple cache for repeated queries

def _get_deg_file_path(contrast_id: str) -> Optional[str]:
    """Find the DEG file path for a given contrast ID."""
    for root, dirs, files in os.walk(RESULTS_DIR):
        if contrast_id in os.path.basename(root) and "DEG.csv" in files:
            return os.path.join(root, "DEG.csv")
    return None

def _load_deg_data(contrast_id: str) -> Optional[pd.DataFrame]:
    """Load DEG data with caching."""
    cache_key = f"deg_{contrast_id}"
    if cache_key in _cache:
        return _cache[cache_key]

    file_path = _get_deg_file_path(contrast_id)
    if not file_path or not os.path.exists(file_path):
        return None

    try:
        df = pd.read_csv(file_path)
        _cache[cache_key] = df
        return df
    except Exception as e:
        logger.error(f"Error loading DEG file for {contrast_id}: {e}")
        return None

@mcp.tool()
async def get_deg_counts(contrast_id: str,
                        p_value_threshold: float = 0.05,
                        lfc_threshold: float = 1.0) -> str:
    """
    Get the number of differentially expressed genes for a specific contrast.

    Args:
        contrast_id: Identifier for the contrast.
        p_value_threshold: P-value threshold for significance.
        lfc_threshold: Log fold change threshold for significance.

    Returns:
        JSON string with count statistics.
    """
    try:
        df = _load_deg_data(contrast_id)
        if df is None:
            return json.dumps({"error": f"No DEG data found for contrast: {contrast_id}"})

        # Identify p-value and logFC columns
        p_value_col = next((c for c in df.columns if 'adj.P.Val' in c), None)
        if not p_value_col:
            p_value_col = next((c for c in df.columns if 'P.Value' in c), None)

        lfc_col = next((c for c in df.columns if 'logFC' in c), None)

        if not p_value_col or not lfc_col:
            return json.dumps({"error": "Missing required columns in DEG file"})

        # Count DEGs
        sig_mask = (df[p_value_col] < p_value_threshold) & (abs(df[lfc_col]) > lfc_threshold)
        up_mask = sig_mask & (df[lfc_col] > 0)
        down_mask = sig_mask & (df[lfc_col] < 0)

        total_count = sig_mask.sum()
        up_count = up_mask.sum()
        down_count = down_mask.sum()

        return json.dumps({
            "total_count": int(total_count),
            "up_regulated": int(up_count),
            "down_regulated": int(down_count),
            "thresholds": {
                "p_value": p_value_threshold,
                "lfc": lfc_threshold
            }
        })

    except Exception as e:
        logger.error(f"Error in get_deg_counts: {e}")
        return json.dumps({"error": str(e)})

@mcp.tool()
async def find_common_degs(contrasts: List[str],
                          min_frequency: int = 2,
                          p_value_threshold: float = 0.05,
                          lfc_threshold: float = 1.0,
                          max_genes: int = 50) -> str:
    """
    Find genes that are differentially expressed across multiple contrasts.

    Args:
        contrasts: List of contrast identifiers.
        min_frequency: Minimum number of contrasts a gene must appear in.
        p_value_threshold: P-value threshold for significance.
        lfc_threshold: Log fold change threshold for significance.
        max_genes: Maximum number of genes to return.

    Returns:
        JSON string with common DEGs and their properties.
    """
    try:
        # Track genes across contrasts
        gene_occurrences = {}
        gene_data = {}

        for contrast_id in contrasts:
            df = _load_deg_data(contrast_id)
            if df is None:
                continue

            # Identify columns
            p_value_col = next((c for c in df.columns if 'adj.P.Val' in c), None)
            if not p_value_col:
                p_value_col = next((c for c in df.columns if 'P.Value' in c), None)

            lfc_col = next((c for c in df.columns if 'logFC' in c), None)
            gene_col = next((c for c in df.columns if c.lower() in ['gene', 'symbol']), None)

            if not all([p_value_col, lfc_col, gene_col]):
                continue

            # Filter for significant genes
            sig_df = df[(df[p_value_col] < p_value_threshold) & (abs(df[lfc_col]) > lfc_threshold)]

            # Count occurrences
            for _, row in sig_df.iterrows():
                gene = row[gene_col]

                if gene not in gene_occurrences:
                    gene_occurrences[gene] = 0
                    gene_data[gene] = {
                        "contrasts": {},
                        "max_abs_lfc": 0
                    }

                gene_occurrences[gene] += 1
                gene_data[gene]["contrasts"][contrast_id] = {
                    "lfc": float(row[lfc_col]),
                    "p_value": float(row[p_value_col])
                }
                gene_data[gene]["max_abs_lfc"] = max(gene_data[gene]["max_abs_lfc"], abs(row[lfc_col]))

        # Filter by frequency
        common_genes = [
            {
                "gene": gene,
                "frequency": gene_occurrences[gene],
                "contrasts": gene_data[gene]["contrasts"],
                "max_abs_lfc": gene_data[gene]["max_abs_lfc"]
            }
            for gene, freq in gene_occurrences.items()
            if freq >= min_frequency
        ]

        # Sort by frequency (primary) and max_abs_lfc (secondary)
        common_genes.sort(key=lambda x: (x["frequency"], x["max_abs_lfc"]), reverse=True)

        # Limit number of returned genes
        common_genes = common_genes[:max_genes]

        return json.dumps({
            "genes": common_genes,
            "total_found": len(common_genes),
            "contrasts_analyzed": len(contrasts)
        })

    except Exception as e:
        logger.error(f"Error in find_common_degs: {e}")
        return json.dumps({"error": str(e)})

# Additional tools would be implemented here...

def main():
    """Main entry point for the MCP server."""
    logging.info(f"Starting data extractor MCP server with results_dir: {RESULTS_DIR}")
    if not RESULTS_DIR:
        logging.warning("RESULTS_DIR environment variable not set")

    # Run the MCP server
    mcp.run(transport="stdio")

if __name__ == "__main__":
    main()
```

### 4.2. Analysis MCP Server

```python
#!/usr/bin/env python3
# mcp_analysis.py

import os
import sys
import json
import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from typing import List, Dict, Optional, Tuple
import logging
from collections import Counter

from mcp.server.fastmcp import FastMCP

# Initialize MCP server
mcp = FastMCP("analysis")
log = lambda m: print(f"MCP ANALYSIS: {m}", file=sys.stderr, flush=True)

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    stream=sys.stderr
)
logger = logging.getLogger("mcp-analysis")

# Global variables
RESULTS_DIR = os.environ.get("RESULTS_DIR", "")
_cache = {}  # Simple cache for repeated queries

def _get_contrast_relevance(contrast_info, research_question):
    """Heuristic function to score contrast relevance based on research question."""
    # In practice, this would use a more sophisticated approach,
    # potentially calling the LLM for relevance scoring
    score = 0
    description = contrast_info.get("description", "").lower()
    question_terms = research_question.lower().split()

    # Count matching terms
    for term in question_terms:
        if len(term) > 3 and term in description:  # Avoid short words
            score += 1

    # Boost score for contrasts with more DEGs
    score *= (1 + 0.1 * min(10, contrast_info.get("deg_count", 0) / 100))

    return score

@mcp.tool()
async def rank_contrasts_by_relevance(research_question: str,
                                     max_contrasts: int = 10) -> str:
    """
    Rank contrasts by relevance to a specified research question.

    Args:
        research_question: The research question to analyze.
        max_contrasts: Maximum number of contrasts to return.

    Returns:
        JSON string with ranked contrasts and relevance scores.
    """
    try:
        all_contrasts = []

        # Find all contrasts and their metadata
        for root, dirs, files in os.walk(RESULTS_DIR):
            if "DEG.csv" in files:
                contrast_id = os.path.basename(root)
                analysis_id = os.path.basename(os.path.dirname(os.path.dirname(root)))

                # Get contrast description
                description = "No description available"
                metadata_file = os.path.join(RESULTS_DIR, analysis_id, "metadata", "contrasts.csv")
                if os.path.exists(metadata_file):
                    try:
                        contrasts_df = pd.read_csv(metadata_file)
                        match = contrasts_df[contrasts_df['name'] == contrast_id]
                        if not match.empty:
                            description = match.iloc[0].get('description', description)
                    except:
                        pass

                # Count DEGs
                deg_count = 0
                try:
                    deg_df = pd.read_csv(os.path.join(root, "DEG.csv"))
                    p_value_col = next((c for c in deg_df.columns if 'adj.P.Val' in c), None)
                    lfc_col = next((c for c in deg_df.columns if 'logFC' in c), None)

                    if p_value_col and lfc_col:
                        deg_count = ((deg_df[p_value_col] < 0.05) & (abs(deg_df[lfc_col]) > 1.0)).sum()
                except:
                    pass

                contrast_info = {
                    "contrast_id": contrast_id,
                    "analysis_id": analysis_id,
                    "description": description,
                    "deg_count": deg_count
                }

                all_contrasts.append(contrast_info)

        # Score contrasts by relevance
        for contrast in all_contrasts:
            contrast["relevance_score"] = _get_contrast_relevance(contrast, research_question)

        # Sort by relevance score
        all_contrasts.sort(key=lambda x: x["relevance_score"], reverse=True)

        # Take top N contrasts
        top_contrasts = all_contrasts[:max_contrasts]

        return json.dumps({
            "ranked_contrasts": top_contrasts,
            "research_question": research_question
        })

    except Exception as e:
        logger.error(f"Error in rank_contrasts_by_relevance: {e}")
        return json.dumps({"error": str(e)})

@mcp.tool()
async def analyze_gene_patterns(genes: List[str], contrasts: List[str]) -> str:
    """
    Analyze expression patterns across a set of genes.

    Args:
        genes: List of gene identifiers.
        contrasts: List of contrast identifiers.

    Returns:
        JSON string with pattern analysis.
    """
    try:
        # Build a gene expression matrix
        expr_matrix = []

        for contrast in contrasts:
            # Find the DEG file
            file_path = None
            for root, _, files in os.walk(RESULTS_DIR):
                if contrast in os.path.basename(root) and "DEG.csv" in files:
                    file_path = os.path.join(root, "DEG.csv")
                    break

            if not file_path:
                continue

            try:
                df = pd.read_csv(file_path)
                gene_col = next((c for c in df.columns if c.lower() in ['gene', 'symbol']), None)
                lfc_col = next((c for c in df.columns if 'logFC' in c), None)

                if not gene_col or not lfc_col:
                    continue

                # Filter to our genes of interest
                filtered_df = df[df[gene_col].isin(genes)]

                # Add data to matrix
                for gene in genes:
                    gene_row = filtered_df[filtered_df[gene_col] == gene]
                    if not gene_row.empty:
                        lfc = float(gene_row.iloc[0][lfc_col])
                    else:
                        lfc = 0.0

                    expr_matrix.append({
                        "gene": gene,
                        "contrast": contrast,
                        "lfc": lfc
                    })
            except Exception as e:
                logger.error(f"Error processing contrast {contrast}: {e}")

        # Convert to wide format for pattern analysis
        if not expr_matrix:
            return json.dumps({"error": "No expression data found"})

        expr_df = pd.DataFrame(expr_matrix)
        expr_wide = expr_df.pivot(index='gene', columns='contrast', values='lfc')
        expr_wide = expr_wide.fillna(0)

        # Calculate correlation between genes
        corr_matrix = expr_wide.T.corr().values
        np.fill_diagonal(corr_matrix, 0)  # Remove self-correlations

        # Group genes with similar patterns
        patterns = []
        for i, gene in enumerate(expr_wide.index):
            correlations = corr_matrix[i]
            positive_corr = [(expr_wide.index[j], float(correlations[j]))
                           for j in range(len(correlations))
                           if correlations[j] > 0.7]  # Strongly positively correlated
            negative_corr = [(expr_wide.index[j], float(correlations[j]))
                           for j in range(len(correlations))
                           if correlations[j] < -0.7]  # Strongly negatively correlated

            patterns.append({
                "gene": gene,
                "expression_values": expr_wide.loc[gene].to_dict(),
                "co_regulated_genes": positive_corr,
                "anti_regulated_genes": negative_corr
            })

        return json.dumps({
            "patterns": patterns,
            "contrasts_analyzed": len(contrasts)
        })

    except Exception as e:
        logger.error(f"Error in analyze_gene_patterns: {e}")
        return json.dumps({"error": str(e)})

# Additional tools would be implemented here...

def main():
    """Main entry point for the MCP server."""
    logging.info(f"Starting analysis MCP server with results_dir: {RESULTS_DIR}")
    if not RESULTS_DIR:
        logging.warning("RESULTS_DIR environment variable not set")

    # Run the MCP server
    mcp.run(transport="stdio")

if __name__ == "__main__":
    main()
```

### 4.3. Agent Controller Integration

```python
#!/usr/bin/env python3
# reporting_agent.py

import os
import sys
from typing import List, Dict, Any, Optional
from pathlib import Path
import json
import logging

# Import MCP server utilities
from mcp_utils import setup_mcp_server

# Import agent framework
from pydantic_ai.agent import Agent
from pydantic_ai.openai import OpenAIModel, OpenAIProvider

class ReportingAgent:
    """Controller for the RNA-seq reporting agent."""

    def __init__(self, results_dir: str, model_name: str = "gpt-4.1-mini"):
        """Initialize the reporting agent."""
        self.results_dir = os.path.abspath(results_dir)
        self.model_name = model_name
        self.servers = []
        self.agent = None

    def setup_servers(self):
        """Set up the MCP servers for data extraction and analysis."""
        # Data extraction server
        data_server = setup_mcp_server(
            "data_extractor",
            env_vars={"RESULTS_DIR": self.results_dir}
        )
        self.servers.append(data_server)

        # Analysis server
        analysis_server = setup_mcp_server(
            "analysis",
            env_vars={"RESULTS_DIR": self.results_dir}
        )
        self.servers.append(analysis_server)

        return self.servers

    def create_agent(self, research_question: str):
        """Create the agent with the specified research question."""
        # Set up servers if not already done
        if not self.servers:
            self.setup_servers()

        # Check if using local or remote model
        if self.model_name.startswith("openai:"):
            model_instance = self.model_name
        else:
            # For local models like Llama
            model_instance = OpenAIModel(
                model_name=self.model_name,
                provider=OpenAIProvider(base_url='http://localhost:11434/v1')
            )

        # Create the agent
        self.agent = Agent(
            model=model_instance,
            model_settings={"temperature": 0.1},
            mcp_servers=self.servers,
            system_prompt=self._get_system_prompt(research_question),
        )

        return self.agent

    def _get_system_prompt(self, research_question: str) -> str:
        """Generate the system prompt for the agent."""
        return f"""
        You are an expert computational biologist analyzing RNA-seq results.

        Your task is to identify key findings related to: '{research_question}'

        Follow this structured analysis approach:
        1. First, identify the most relevant contrasts using the rank_contrasts_by_relevance tool
        2. For these relevant contrasts, find common differentially expressed genes using find_common_degs
        3. Analyze expression patterns across these genes using analyze_gene_patterns
        4. Explore individual gene details as needed using get_gene_info

        Then synthesize the findings into a comprehensive report with:
        - A concise executive summary answering the research question
        - Relevance justifications for selected contrasts
        - Key gene pattern descriptions with biological interpretations
        - Suggested biological mechanisms and hypotheses

        Present your findings in a clear, scientific style with markdown formatting.
        """

    async def generate_report(self, research_question: str):
        """Generate a report for the specified research question."""
        # Create the agent if needed
        if not self.agent:
            self.create_agent(research_question)

        # Run the agent to generate the report
        response = await self.agent.run(
            "Please analyze the RNA-seq data and generate a comprehensive " +
            f"report addressing this research question: {research_question}"
        )

        return {
            "report": response,
            "research_question": research_question
        }

    def cleanup(self):
        """Clean up resources."""
        for server in self.servers:
            server.close()
```

### 4.4. Auto-Start Manager

```python
#!/usr/bin/env python3
# auto_start_manager.py

import asyncio
import os
import time
from typing import Dict, Any, Tuple, Optional
import logging
from pathlib import Path

from reporting_agent import ReportingAgent

class AutoStartManager:
    """Manages the automatic startup and initial analysis of RNA-seq data."""

    def __init__(self, results_dir: str, default_model: str = "openai:gpt-4.1-mini"):
        """Initialize the auto-start manager."""
        self.results_dir = os.path.abspath(results_dir)
        self.default_model = default_model
        self.agent = None
        self.logger = logging.getLogger("auto-start-manager")

    def validate_data_directory(self) -> Tuple[bool, Optional[str]]:
        """Validate that the results directory contains analyzable data."""
        try:
            # Check if directory exists
            if not os.path.exists(self.results_dir):
                return False, "Results directory does not exist"

            # Check for at least one RNAseqAnalysis directory
            found_rnaseq = False
            for root, dirs, _ in os.walk(self.results_dir):
                if "RNAseqAnalysis" in dirs:
                    found_rnaseq = True
                    break

            if not found_rnaseq:
                return False, "No RNAseqAnalysis directory found"

            # Success
            return True, None

        except Exception as e:
            self.logger.error(f"Error validating data directory: {e}")
            return False, f"Error: {str(e)}"

    async def run_initial_analysis(self,
                                 progress_callback=None) -> Dict[str, Any]:
        """
        Run initial analysis on the data without requiring a specific research question.

        Args:
            progress_callback: Optional function to call with progress updates (0-1)

        Returns:
            Dictionary with analysis results
        """
        try:
            # Create reporting agent if not already done
            if not self.agent:
                self.agent = ReportingAgent(
                    results_dir=self.results_dir,
                    model_name=self.default_model
                )
                self.agent.setup_servers()

            # Define default research question for initial overview
            default_question = (
                "Provide a concise overview of the key differential expression patterns "
                "and identify the most significant contrasts and genes in this dataset. "
                "Focus on general patterns rather than specific research questions."
            )

            # Call progress callback if provided
            if progress_callback:
                progress_callback(0.1, "Initializing agent")

            # Generate the report
            if progress_callback:
                # Hook up progress updates
                async def progress_updater():
                    for i in range(2, 10):
                        await asyncio.sleep(1)
                        progress_callback(i * 0.1, f"Analyzing data ({i*10}%)")

                # Start progress updates in background
                progress_task = asyncio.create_task(progress_updater())

                # Generate report
                result = await self.agent.generate_report(default_question)

                # Cancel progress task
                progress_task.cancel()
                progress_callback(1.0, "Analysis complete")
            else:
                # Generate report without progress updates
                result = await self.agent.generate_report(default_question)

            return {
                "success": True,
                "report": result,
                "research_question": default_question
            }

        except Exception as e:
            self.logger.error(f"Error during initial analysis: {e}")
            return {
                "success": False,
                "error": str(e)
            }

    def cleanup(self):
        """Clean up resources."""
        if self.agent:
            self.agent.cleanup()
```

### 4.5. Streamlit Integration with Auto-Start

```python
# In uorca_explorer.py

# Import the AutoStartManager
from auto_start_manager import AutoStartManager

# Modified results directory handling
results_dir = st.sidebar.text_input(
    "Results directory",
    value=default_dir
)

# Create a key for caching based on directory
cache_key = f"data_{results_dir}"

# Initialize and handle auto-start
if results_dir != st.session_state.get('previous_results_dir', None):
    # Path changed - reset state
    st.session_state.previous_results_dir = results_dir
    if 'initial_analysis_complete' in st.session_state:
        del st.session_state.initial_analysis_complete
    if 'agent_report' in st.session_state:
        del st.session_state.agent_report

# Initialize data and run auto-analysis
with st.sidebar.status("Working with data...", expanded=True) as status:
    # First load the data
    ri, error = get_integrator(results_dir)

    if error:
        status.update(label=f"Error loading data: {error}", state="error")
    elif not ri or not ri.cpm_data:
        status.update(label="No data found. Please check the directory path.", state="error")
    else:
        # Data loaded successfully
        status.update(label=f"✅ Loaded {len(ri.cpm_data)} datasets", state="complete")

        # If this is first time with this directory and auto-analysis not done yet
        if 'initial_analysis_complete' not in st.session_state:
            # Show analyzing message
            status.update(label="🧠 Analyzing data automatically...", state="running")

            # Create auto-start manager
            auto_start = AutoStartManager(results_dir)

            # Check if directory is valid for analysis
            is_valid, validation_msg = auto_start.validate_data_directory()

            if is_valid:
                # Valid directory - run analysis in a background task
                st.session_state.analysis_started = True

                # Will continue in landing page tab
                st.session_state.trigger_initial_analysis = True
            else:
                # Invalid directory - show error
                status.update(label=f"⚠️ {validation_msg}", state="error")
                st.session_state.initial_analysis_complete = False
```

### 4.6. Landing Page Implementation with Automatic Startup

```python
# In the landing page tab section of uorca_explorer.py

with tab_landing:
    st.header("🤖 AI Analysis & Interpretation")

    # Display initial loading message if this is the first load
    if st.session_state.get('trigger_initial_analysis', False):
        # Clear the trigger flag
        st.session_state.trigger_initial_analysis = False

        # Show initial loading message with spinner
        with st.spinner("🧠 Automatically analyzing your data..."):
            try:
                # Use a default research context for initial analysis
                default_question = "Provide an overview of key differential expression patterns"

                # Run the agent with this default question
                loop = asyncio.new_event_loop()
                asyncio.set_event_loop(loop)

                # Progress indicator
                progress_bar = st.progress(0)

                # Generate initial report
                for i in range(5):
                    # Simulate progress steps
                    progress_bar.progress((i+1)/5)
                    time.sleep(0.5)

                result = loop.run_until_complete(agent.generate_report(default_question))

                # Store result in session state
                st.session_state['agent_report'] = result
                st.session_state['research_question'] = default_question
                st.session_state['initial_analysis_complete'] = True

                # Force refresh to show the results
                st.rerun()

            except Exception as e:
                st.error(f"Error in automatic analysis: {str(e)}")
                st.session_state['initial_analysis_complete'] = False

    # Rest of the landing page UI code...
    # (This will show after the automatic analysis completes)
    render_agent_landing_page(results_dir, ri, agent)
```

### 4.7. Landing Page Rendering

```python
def render_agent_landing_page(results_dir, ri, agent, auto_run=False):
    """Render the agent-based landing page."""

    # Display initial report if available
    if 'agent_report' in st.session_state:
        report = st.session_state['agent_report']

        st.success("✅ Initial analysis complete!")
        st.info("This is an automatic analysis of your data. You can refine it with a specific research question below.")

        # Display the report...
        st.markdown("## 📊 Initial Analysis")
        st.markdown(report['report'])

        st.markdown("---")

    # Always show the option to ask a specific research question
    st.markdown("## 🔍 Ask a Specific Research Question")

    research_question = st.text_area(
        "🧬 Research Question",
        value="",
        height=100,
        placeholder="e.g., 'Identify key regulatory genes involved in inflammatory response'",
        help="Describe your biological research question for more targeted analysis"
    )

    # Controls for report generation
    col1, col2 = st.columns([3, 1])

    with col1:
        model_selection = st.selectbox(
            "🧠 LLM Model",
            ["openai:gpt-4.1-mini", "openai:gpt-4o", "llama3.3:latest"],
            index=0,
            help="Select which language model to use for analysis"
        )

    with col2:
        generate_button = st.button(
            "🚀 Generate Analysis",
            type="primary",
            use_container_width=True,
            help="Analyze data with your specific research question"
        )

    # Handle custom report generation
    if generate_button and research_question.strip():
        with st.spinner("🧠 AI agent analyzing RNA-seq data for your specific question..."):
            try:
                # Use the agent to generate the report
                loop = asyncio.new_event_loop()
                asyncio.set_event_loop(loop)

                # Set agent model if changed
                if agent.model_name != model_selection:
                    agent.model_name = model_selection
                    agent.agent = None  # Force recreation

                # Progress indicator
                progress_bar = st.progress(0)
                for i in range(10):
                    # Simulate progress steps
                    progress_bar.progress((i+1)/10)
                    time.sleep(0.3)

                # Generate the report
                result = loop.run_until_complete(agent.generate_report(research_question))

                # Store in session state
                st.session_state['agent_report'] = result
                st.session_state['research_question'] = research_question

                # Force refresh to show results
                st.rerun()

            except Exception as e:
                st.error(f"Error generating report: {str(e)}")
```

## 5. Setup Scripts

Create setup scripts to prepare the environment:

```bash
#!/bin/bash
# setup_agent_reporting.sh

# Create virtual environment
python -m venv venv
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt

# Copy MCP server scripts to appropriate locations
mkdir -p mcp_servers
cp mcp_data_extractor.py mcp_servers/
cp mcp_analysis.py mcp_servers/

# Update servers.toml
cat << EOF > servers.toml
[servers.data_extractor]
enabled = true
path = "mcp_servers/mcp_data_extractor.py"

[servers.analysis]
enabled = true
path = "mcp_servers/mcp_analysis.py"
EOF

echo "Agent reporting system setup complete!"
```

## 6. Evaluation Criteria

1. **Startup Performance**: Measure time from directory selection to initial results display
2. **Accuracy**: Compare agent-generated insights with known biological patterns
3. **Relevance**: Assess how well reported findings address the research question
4. **User Experience**: Evaluate the smoothness of the automatic startup flow
5. **Error Handling**: Test recovery from various failure scenarios

## 7. Risk Management

1. **Performance Risk**: Initial analysis may be slow for large datasets
   - **Mitigation**: Implement tiered analysis with quick summary first

2. **Error Handling Risk**: Automatic startup increases complexity of error states
   - **Mitigation**: Design comprehensive error handling with graceful fallbacks

3. **User Expectation Risk**: Default analysis may not match user expectations
   - **Mitigation**: Clearly explain what the automatic analysis provides and its limitations

4. **Resource Usage Risk**: Automatic analysis increases baseline resource requirements
   - **Mitigation**: Implement configurable analysis depth and resource limits

## 8. Future Extensions

1. **Advanced Pattern Recognition**: Implement more sophisticated statistical methods for identifying co-regulated gene sets
2. **Pathway Analysis**: Add pathway and gene ontology enrichment analysis tools
3. **Publication Support**: Generate draft methods and results sections for publications
4. **Multi-Omics Integration**: Extend to support integration with proteomics, metabolomics data
5. **Customizable Report Templates**: Allow users to define the structure and focus of generated reports

This comprehensive implementation plan provides a clear roadmap for creating an agent-based reporting system that automatically runs when the user specifies their results directory, with no additional steps required.
