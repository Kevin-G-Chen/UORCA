# ----------------------------
# Imports
# ----------------------------
import argparse
import os
import re
import json
import pandas as pd
from dotenv import load_dotenv
from rich.console import Console
from rich.panel import Panel
from dataclasses import dataclass
from typing import List, Optional
from unidecode import unidecode
from pydantic import BaseModel, Field, ConfigDict
from pydantic_ai import Agent, RunContext
from openai import OpenAI

# ----------------------------
# Configure logging
# ----------------------------


class LogLevel:
    MINIMAL = 0   # Only critical information
    NORMAL = 1    # Default level
    VERBOSE = 2   # Detailed information
    DEBUG = 3     # Maximum debugging information


CURRENT_LOG_LEVEL = LogLevel.NORMAL
console = Console()


def log(message, level=LogLevel.NORMAL, style=""):
    if CURRENT_LOG_LEVEL >= level:
        console.print(message, style=style if style else None)


def log_tool_header(tool_name, params=None):
    if CURRENT_LOG_LEVEL >= LogLevel.NORMAL:
        console.print(f"TOOL: {tool_name}", style="bold blue")
        if params:
            console.print(f"Parameters: {params}")


def log_tool_result(result):
    if CURRENT_LOG_LEVEL >= LogLevel.NORMAL:
        console.print(result)

# ----------------------------
# Dependency Class
# ----------------------------


@dataclass
class RNAseqData:
    """Container for metadata analysis data."""
    metadata_path: str
    metadata_df: Optional[pd.DataFrame] = None
    # This will hold the final analysis column name (or merged version)
    merged_column: Optional[str] = None
    # To store designed contrasts if needed
    contrast_details: Optional[dict] = None


# ----------------------------
# Load environment variables and initialize OpenAI client
# ----------------------------
load_dotenv()
openai_api_key = os.getenv("OPENAI_API_KEY")
client = OpenAI()

# ----------------------------
# Define output schema for candidate columns
# ----------------------------


class Contrast_format(BaseModel):
    """Schema for the output of the candidate contrasts."""
    name: str
    expression: str

class Contrasts(BaseModel):
    """Schema for the output of the candidate contrasts."""
    contrasts: List[Contrast_format]


# ----------------------------
# Create RNAseq metadata analysis agent
# ----------------------------
# Try reading your system prompt, otherwise use the fallback prompt
system_prompt_path = "systemprompt.txt"  # Update this path as needed
try:
    with open(system_prompt_path, 'r') as f:
        system_prompt = f.read()
except Exception as e:
    system_prompt = """
    #### Integrated Prompt for Metadata Processing and Grouping Variable Selection

    You are provided with RNAseq metadata from two different experiments. Your task is to identify the column(s) that contain biologically relevant information for differential expression analysis and to merge them into a single grouping variable if necessary. This grouping variable will be used in a single-factor analysis with edgeR/limma. In doing so, you must also evaluate each column to decide which ones provide informative biological variation and which ones should be excluded.

    General Guidelines:
	1.	Focus on Biologically Relevant Information:
	•	Include columns that capture sample-specific biological attributes such as tissue/disease type, genotype, or treatment conditions.
	•	Exclude technical columns (e.g., sample IDs, run/experiment numbers) and those with no variation (all values identical) or with unique values that do not group samples.
	2.	Merging Columns:
	•	If more than one column is informative (e.g., one column for tissue type and one for treatment), merge these into a single composite grouping variable (for example, merged_analysis_group).
	•	Ensure that the final grouping factor includes only information that is biologically significant for differential expression analysis.
	3.	Output:
	•	Return the name(s) of the final grouping column(s) and a brief explanation of your selection process and why the other columns were excluded.

    ⸻

    #### Examples

    Metadata Set 1

    Metadata Table:

    geo_accession	title	organism_ch1	characteristics_ch1	characteristics_ch1.1	characteristics_ch1.2	molecule_ch1	cell type:ch1	tissue:ch1	treatment:ch1	Run	Experiment
    GSM6443387	head & neck squamous cell carcinoma patient	Homo sapiens	tissue: head & neck squamous cell carcinoma	cell type: peripheral blood mononuclear cells (PBMC)	treatment: No treatment	total RNA	peripheral blood mononuclear cells (PBMC)	head & neck squamous cell carcinoma	No treatment	SRR21008833	SRX17025351
    GSM6443387	head & neck squamous cell carcinoma patient	Homo sapiens	tissue: head & neck squamous cell carcinoma	cell type: peripheral blood mononuclear cells (PBMC)	treatment: No treatment	total RNA	peripheral blood mononuclear cells (PBMC)	head & neck squamous cell carcinoma	No treatment	SRR21008834	SRX17025351
    GSM6443388	breast cancer patient	Homo sapiens	tissue: breast cancer	cell type: peripheral blood mononuclear cells (PBMC)	treatment: No treatment	total RNA	peripheral blood mononuclear cells (PBMC)	breast cancer	No treatment	SRR21008843	SRX17025346
    GSM6443388	breast cancer patient	Homo sapiens	tissue: breast cancer	cell type: peripheral blood mononuclear cells (PBMC)	treatment: No treatment	total RNA	peripheral blood mononuclear cells (PBMC)	breast cancer	No treatment	SRR21008844	SRX17025346
    GSM6443389	advanced melanoma patient 1, DMSO	Homo sapiens	tissue: advanced melanoma	cell type: peripheral blood mononuclear cells (PBMC)	treatment: DMSO	total RNA	peripheral blood mononuclear cells (PBMC)	advanced melanoma	DMSO	SRR21008841	SRX17025347
    GSM6443389	advanced melanoma patient 1, DMSO	Homo sapiens	tissue: advanced melanoma	cell type: peripheral blood mononuclear cells (PBMC)	treatment: DMSO	total RNA	peripheral blood mononuclear cells (PBMC)	advanced melanoma	DMSO	SRR21008842	SRX17025347
    GSM6443390	advanced melanoma patient 1, 5uM Ibrutinib	Homo sapiens	tissue: advanced melanoma	cell type: peripheral blood mononuclear cells (PBMC)	treatment: 5uM Ibrutinib	total RNA	peripheral blood mononuclear cells (PBMC)	advanced melanoma	5uM Ibrutinib	SRR21008839	SRX17025348
    GSM6443390	advanced melanoma patient 1, 5uM Ibrutinib	Homo sapiens	tissue: advanced melanoma	cell type: peripheral blood mononuclear cells (PBMC)	treatment: 5uM Ibrutinib	total RNA	peripheral blood mononuclear cells (PBMC)	advanced melanoma	5uM Ibrutinib	SRR21008840	SRX17025348
    GSM6443391	advanced melanoma patient 2, DMSO	Homo sapiens	tissue: advanced melanoma	cell type: peripheral blood mononuclear cells (PBMC)	treatment: DMSO	total RNA	peripheral blood mononuclear cells (PBMC)	advanced melanoma	DMSO	SRR21008837	SRX17025349
    GSM6443391	advanced melanoma patient 2, DMSO	Homo sapiens	tissue: advanced melanoma	cell type: peripheral blood mononuclear cells (PBMC)	treatment: DMSO	total RNA	peripheral blood mononuclear cells (PBMC)	advanced melanoma	DMSO	SRR21008838	SRX17025349
    GSM6443392	advanced melanoma patient 2, 5uM Ibrutinib	Homo sapiens	tissue: advanced melanoma	cell type: peripheral blood mononuclear cells (PBMC)	treatment: 5uM Ibrutinib	total RNA	peripheral blood mononuclear cells (PBMC)	advanced melanoma	5uM Ibrutinib	SRR21008835	SRX17025350
    GSM6443392	advanced melanoma patient 2, 5uM Ibrutinib	Homo sapiens	tissue: advanced melanoma	cell type: peripheral blood mononuclear cells (PBMC)	treatment: 5uM Ibrutinib	total RNA	peripheral blood mononuclear cells (PBMC)	advanced melanoma	5uM Ibrutinib	SRR21008836	SRX17025350

    Column Evaluation for Metadata Set 1:
	1.	geo_accession:
	•	Contains unique sample identifiers.
	•	Not included: Technical ID; no biological grouping information.
	2.	title:
	•	Provides a description (e.g., cancer type, treatment hint).
	•	Marginal utility: Redundant with structured columns; less reliable.
	3.	organism_ch1:
	•	Always “Homo sapiens”.
	•	Not included: No variation; does not aid grouping.
	4.	characteristics_ch1:
	•	Shows tissue/disease type (e.g., “tissue: head & neck squamous cell carcinoma”).
	•	Good candidate: Captures key biological context.
	5.	characteristics_ch1.1:
	•	Specifies cell type (all are PBMC).
	•	Not included: No variation across samples.
	6.	characteristics_ch1.2:
	•	Details treatment (e.g., “treatment: DMSO”, “treatment: 5uM Ibrutinib”).
	•	Good candidate: Provides important treatment differences.
	7.	molecule_ch1:
	•	Indicates molecule type (“total RNA”).
	•	Not included: Constant across samples.
	8.	cell type:ch1:
	•	Redundant to characteristics_ch1.1.
	•	Not included.
	9.	tissue:ch1:
	•	Repeats tissue/disease type.
	•	Good candidate (redundant with characteristics_ch1): Only one is needed.
	10.	treatment:ch1:
	•	Repeats treatment information.
	•	Good candidate (redundant with characteristics_ch1.2): Only one is needed.
	11.	Run:
	•	Sequencing run identifier.
	•	Not included: Technical detail.
	12.	Experiment:
	•	Sequencing experiment identifier.
	•	Not included: Technical detail.

    Final Assessment for Metadata Set 1:
    Merge the tissue/disease column (either characteristics_ch1 or tissue:ch1) with the treatment column (either characteristics_ch1.2 or treatment:ch1) into a composite grouping variable (e.g., merged_analysis_group). This captures the key biological differences across samples.

    ⸻

    Metadata Set 2

    Metadata Table:

    geo_accession	title	channel_count	source_name_ch1	organism_ch1	characteristics_ch1	characteristics_ch1.1	characteristics_ch1.2	characteristics_ch1.3	characteristics_ch1.4	characteristics_ch1.5	characteristics_ch1.6	characteristics_ch1.7	characteristics_ch1.8	supplementary_file_1	age:ch1	genotype:ch1	nomenclature for_ins1-cre_developed_by_the_thorens_group:ch1	nomenclature for_stim1_fl/fl:ch1	Sex:ch1	stock number_on_jackson_laboratory_(stim1fl/fl):ch1	nomenclature for_stim1_fl/fl:ch1	stock number_on_jackson_laboratory:ch1	tissue:ch1	treatment:ch1	Run	Experiment
    GSM6337959	Pancreatic islet, Control, 1	1	Pancreatic islet	Mus musculus	tissue: Pancreatic islet	genotype: Control (STIM1fl/fl, Cre-)	treatment: High-fat diet (ResearchDiets D12492) 8wk	Sex: Female	age: 16wk	nomenclature for_stim1_fl/fl: B6(Cg)-STIM1tm1Rao/J	stock number_on_jackson_laboratory_(stim1fl/fl): 23350	nomenclature for_ins1-cre_developed_by_the_thorens_group: B6(Cg)-Ins1tm1.1(cre)Thor/J	stock number_on_jackson_laboratory: 26801		NONE	16wk	Control (STIM1fl/fl, Cre-)	B6(Cg)-Ins1tm1.1(cre)Thor/J	B6(Cg)-STIM1tm1Rao/J	Female	23350	26801	Pancreatic islet	High-fat diet (ResearchDiets D12492) 8wk	SRR20166827	SRX16201852
    GSM6337960	Pancreatic islet, Control, 2	1	Pancreatic islet	Mus musculus	tissue: Pancreatic islet	genotype: Control (STIM1fl/fl, Cre-)	treatment: High-fat diet (ResearchDiets D12492) 8wk	Sex: Female	age: 16wk	nomenclature for_stim1_fl/fl: B6(Cg)-STIM1tm1Rao/J	stock number_on_jackson_laboratory_(stim1fl/fl): 23350	nomenclature for_ins1-cre_developed_by_the_thorens_group: B6(Cg)-Ins1tm1.1(cre)Thor/J	stock number_on_jackson_laboratory: 26801		NONE	16wk	Control (STIM1fl/fl, Cre-)	B6(Cg)-Ins1tm1.1(cre)Thor/J	B6(Cg)-STIM1tm1Rao/J	Female	23350	26801	Pancreatic islet	High-fat diet (ResearchDiets D12492) 8wk	SRR20166826	SRX16201853
    GSM6337961	Pancreatic islet, Control, 3	1	Pancreatic islet	Mus musculus	tissue: Pancreatic islet	genotype: Control (STIM1fl/fl, Cre-)	treatment: High-fat diet (ResearchDiets D12492) 8wk	Sex: Female	age: 16wk	nomenclature for_stim1_fl/fl: B6(Cg)-STIM1tm1Rao/J	stock number_on_jackson_laboratory_(stim1fl/fl): 23350	nomenclature for_ins1-cre_developed_by_the_thorens_group: B6(Cg)-Ins1tm1.1(cre)Thor/J	stock number_on_jackson_laboratory: 26801		NONE	16wk	Control (STIM1fl/fl, Cre-)	B6(Cg)-Ins1tm1.1(cre)Thor/J	B6(Cg)-STIM1tm1Rao/J	Female	23350	26801	Pancreatic islet	High-fat diet (ResearchDiets D12492) 8wk	SRR20166823	SRX16201854
    GSM6337962	Pancreatic islet, Control, 4	1	Pancreatic islet	Mus musculus	tissue: Pancreatic islet	genotype: Control (STIM1fl/fl, Cre-)	treatment: High-fat diet (ResearchDiets D12492) 8wk	Sex: Female	age: 16wk	nomenclature for_stim1_fl/fl: B6(Cg)-STIM1tm1Rao/J	stock number_on_jackson_laboratory_(stim1fl/fl): 23350	nomenclature for_ins1-cre_developed_by_the_thorens_group: B6(Cg)-Ins1tm1.1(cre)Thor/J	stock number_on_jackson_laboratory: 26801		NONE	16wk	Control (STIM1fl/fl, Cre-)	B6(Cg)-Ins1tm1.1(cre)Thor/J	B6(Cg)-STIM1tm1Rao/J	Female	23350	26801	Pancreatic islet	High-fat diet (ResearchDiets D12492) 8wk	SRR20166822	SRX16201855
    GSM6337963	Pancreatic islet, Control, 5	1	Pancreatic islet	Mus musculus	tissue: Pancreatic islet	genotype: Control (STIM1fl/fl, Cre-)	treatment: High-fat diet (ResearchDiets D12492) 8wk	Sex: Female	age: 16wk	nomenclature for_stim1_fl/fl: B6(Cg)-STIM1tm1Rao/J	stock number_on_jackson_laboratory_(stim1fl/fl): 23350	nomenclature for_ins1-cre_developed_by_the_thorens_group: B6(Cg)-Ins1tm1.1(cre)Thor/J	stock number_on_jackson_laboratory: 26801		NONE	16wk	Control (STIM1fl/fl, Cre-)	B6(Cg)-Ins1tm1.1(cre)Thor/J	B6(Cg)-STIM1tm1Rao/J	Female	23350	26801	Pancreatic islet	High-fat diet (ResearchDiets D12492) 8wk	SRR20166824	SRX16201856
    GSM6337964	Pancreatic islet, STIM1KO, 1	1	Pancreatic islet	Mus musculus	tissue: Pancreatic islet	genotype: [beta] cell STIM1-Knock out (STIM1fl/fl, Cre+)	treatment: High-fat diet (ResearchDiets D12492) 8wk	Sex: Female	age: 16wk	nomenclature for_stim1_fl/fl: B6(Cg)-STIM1tm1Rao/J	stock number_on_jackson_laboratory_(stim1fl/fl): 23350	nomenclature for_ins1-cre_developed_by_the_thorens_group: B6(Cg)-Ins1tm1.1(cre)Thor/J	stock number_on_jackson_laboratory: 26801		NONE	16wk	[beta] cell STIM1-Knock out (STIM1fl/fl, Cre+)	B6(Cg)-Ins1tm1.1(cre)Thor/J	B6(Cg)-STIM1tm1Rao/J	Female	23350	26801	Pancreatic islet	High-fat diet (ResearchDiets D12492) 8wk	SRR20166821	SRX16201857
    GSM6337965	Pancreatic islet, STIM1KO, 2	1	Pancreatic islet	Mus musculus	tissue: Pancreatic islet	genotype: [beta] cell STIM1-Knock out (STIM1fl/fl, Cre+)	treatment: High-fat diet (ResearchDiets D12492) 8wk	Sex: Female	age: 16wk	nomenclature for_stim1_fl/fl: B6(Cg)-STIM1tm1Rao/J	stock number_on_jackson_laboratory_(stim1fl/fl): 23350	nomenclature for_ins1-cre_developed_by_the_thorens_group: B6(Cg)-Ins1tm1.1(cre)Thor/J	stock number_on_jackson_laboratory: 26801		NONE	16wk	[beta] cell STIM1-Knock out (STIM1fl/fl, Cre+)	B6(Cg)-Ins1tm1.1(cre)Thor/J	B6(Cg)-STIM1tm1Rao/J	Female	23350	26801	Pancreatic islet	High-fat diet (ResearchDiets D12492) 8wk	SRR20166825	SRX16201858
    GSM6337966	Pancreatic islet, STIM1KO, 3	1	Pancreatic islet	Mus musculus	tissue: Pancreatic islet	genotype: [beta] cell STIM1-Knock out (STIM1fl/fl, Cre+)	treatment: High-fat diet (ResearchDiets D12492) 8wk	Sex: Female	age: 16wk	nomenclature for_stim1_fl/fl: B6(Cg)-STIM1tm1Rao/J	stock number_on_jackson_laboratory_(stim1fl/fl): 23350	nomenclature for_ins1-cre_developed_by_the_thorens_group: B6(Cg)-Ins1tm1.1(cre)Thor/J	stock number_on_jackson_laboratory: 26801		NONE	16wk	[beta] cell STIM1-Knock out (STIM1fl/fl, Cre+)	B6(Cg)-Ins1tm1.1(cre)Thor/J	B6(Cg)-STIM1tm1Rao/J	Female	23350	26801	Pancreatic islet	High-fat diet (ResearchDiets D12492) 8wk	SRR20166820	SRX16201859
    GSM6337967	Pancreatic islet, STIM1KO, 4	1	Pancreatic islet	Mus musculus	tissue: Pancreatic islet	genotype: [beta] cell STIM1-Knock out (STIM1fl/fl, Cre+)	treatment: High-fat diet (ResearchDiets D12492) 8wk	Sex: Female	age: 16wk	nomenclature for_stim1_fl/fl: B6(Cg)-STIM1tm1Rao/J	stock number_on_jackson_laboratory_(stim1fl/fl): 23350	nomenclature for_ins1-cre_developed_by_the_thorens_group: B6(Cg)-Ins1tm1.1(cre)Thor/J	stock number_on_jackson_laboratory: 26801		NONE	16wk	[beta] cell STIM1-Knock out (STIM1fl/fl, Cre+)	B6(Cg)-Ins1tm1.1(cre)Thor/J	B6(Cg)-STIM1tm1Rao/J	Female	23350	26801	Pancreatic islet	High-fat diet (ResearchDiets D12492) 8wk	SRR20166819	SRX16201860

    Column Evaluation for Metadata Set 2:
	1.	geo_accession:
	•	Contains unique sample identifiers.
	•	Not included: Only used as an identifier.
	2.	title:
	•	Descriptive title (e.g., “Pancreatic islet, Control, 1”).
	•	Marginal utility: Contains hints about group (e.g., Control vs. STIM1KO) but is less structured than dedicated genotype columns.
	3.	channel_count:
	•	Technical information (e.g., number of channels).
	•	Not included: Does not provide biologically relevant grouping.
	4.	source_name_ch1:
	•	States “Pancreatic islet” for all samples.
	•	Not included: Constant across samples; no grouping power.
	5.	organism_ch1:
	•	Always “Mus musculus”.
	•	Not included: No variation for grouping.
	6.	characteristics_ch1:
	•	Indicates tissue type (“tissue: Pancreatic islet”).
	•	Not included: Constant across all samples.
	7.	characteristics_ch1.1:
	•	Shows genotype information (e.g., “genotype: Control (STIM1fl/fl, Cre-)” vs. “[beta] cell STIM1-Knock out (STIM1fl/fl, Cre+)”).
	•	Good candidate: Provides key biological variation between control and knockout samples.
	8.	characteristics_ch1.2 to characteristics_ch1.8:
	•	Contain treatment details, sex, age, stock numbers, and nomenclature.
	•	Not included: Most of these are constant across samples or technical details; treatment (if present) is identical for all.
	9.	supplementary_file_1:
	•	Indicates additional file information (e.g., “NONE”).
	•	Not included: Not informative for grouping.
	10.	age:ch1:
	•	Age information (e.g., “16wk”).
	•	Not included: Constant for all samples.
	11.	genotype:ch1:
	•	Specifies genotype (e.g., “Control (STIM1fl/fl, Cre-)” vs. “[beta] cell STIM1-Knock out (STIM1fl/fl, Cre+)”).
	•	Good candidate: Captures the only biological variation in this dataset.
	12.	nomenclature for_ins1-cre_developed_by_the_thorens_group:ch1, nomenclature for_stim1_fl/fl:ch1, Sex:ch1, stock number_on_jackson_laboratory_(stim1fl/fl):ch1, stock number_on_jackson_laboratory:ch1:
	•	Contain technical or constant information.
	•	Not included: Do not contribute to sample grouping.
	13.	tissue:ch1:
	•	Indicates tissue type (“Pancreatic islet”).
	•	Not included: Constant across samples.
	14.	treatment:ch1:
	•	Shows treatment details (“High-fat diet (ResearchDiets D12492) 8wk”).
	•	Not included: Identical for all samples.
	15.	Run and Experiment:
	•	Technical sequencing identifiers.
	•	Not included: Used for QC and tracking only.

    Final Assessment for Metadata Set 2:
    The only column that exhibits meaningful biological variation is genotype:ch1. Use this column directly as the grouping variable for downstream analysis.

    """
    print("Using fallback system prompt instead.")

rnaseq_agent = Agent(
    'openai:gpt-4o',         # Use a powerful model
    deps_type=RNAseqData,
    system_prompt=system_prompt
)

# ----------------------------
# Utility function: Clean a string
# ----------------------------


@rnaseq_agent.tool
def clean_string(ctx: RunContext[RNAseqData], s: str) -> str:
    """
    Normalize and clean an input string by removing non-ASCII characters, extra whitespace, and unwanted symbols.
    """
    if pd.isna(s):
        return "NA"
    s = str(s).strip()
    s = unidecode(s)
    s = s.replace(" ", "_")
    s = re.sub(r'[^\w]', '', s)
    return s

# ----------------------------
# Tool: Process metadata
# ----------------------------


@rnaseq_agent.tool
async def process_metadata(ctx: RunContext[RNAseqData]) -> str:
    """
    Load metadata from ctx.deps.metadata_path, remove columns with identical values, clean all column names and
    cell values using clean_string, and store the cleaned DataFrame in the context.
    """
    try:
        log_tool_header("process_metadata", {
                        "metadata_path": ctx.deps.metadata_path})

        # Load metadata based on file extension
        if ctx.deps.metadata_path.endswith('.csv'):
            df = pd.read_csv(ctx.deps.metadata_path)
        elif ctx.deps.metadata_path.endswith(('.tsv', '.txt')):
            df = pd.read_csv(ctx.deps.metadata_path, sep='\t')
        else:
            df = pd.read_csv(ctx.deps.metadata_path, sep=None, engine='python')

        # Remove columns where all values are identical
        df = df.loc[:, df.nunique() > 1]
        # Remove columns where all values are different (i.e. all values are unique)
        df = df.loc[:, df.nunique() < df.shape[0]]

        # Clean column names
        new_columns = {col: clean_string(ctx, col) for col in df.columns}
        df.rename(columns=new_columns, inplace=True)

        # Clean each cell value
        for col in df.columns:
            df[col] = df[col].apply(
                lambda x: clean_string(ctx, x) if pd.notna(x) else x)

        # Store cleaned metadata in context
        ctx.deps.metadata_df = df

        summary = f"""Metadata processed: {df.shape[0]} rows and {df.shape[1]} columns.
Columns: {', '.join(df.columns)}
Preview: {df.head().to_dict(orient='records')}
"""
        log_tool_result(summary)
        return summary
    except Exception as e:
        error_msg = f"Error processing metadata: {str(e)}"
        log(error_msg, style="bold red")
        log_tool_result(error_msg)
        return error_msg

# ----------------------------
# Tool: Merge Analysis Columns
# ----------------------------


@rnaseq_agent.tool
async def merge_analysis_columns(ctx: RunContext[RNAseqData], candidate_cols_json: str) -> str:
    """
    Given a JSON string with candidate columns (output from the LLM prompt), this tool will:
      - If only one candidate column is present, set that column as the analysis column.
      - If multiple candidate columns are provided, merge them by joining the values with an underscore.

    The merged result is stored as a new column 'merged_analysis_group' in ctx.deps.metadata_df,
    and the chosen column name (or new merged column) is stored in ctx.deps.merged_column.

    Expected input JSON format:
      {"columns": ["col1", "col2", ...]}
    """
    try:
        candidates = json.loads(candidate_cols_json).get("columns", [])
        if not candidates:
            msg = "No candidate columns provided."
            log_tool_result(msg)
            return msg

        df = ctx.deps.metadata_df.copy()

        if len(candidates) == 1:
            # Single column; use it directly
            ctx.deps.merged_column = candidates[0]
            msg = f"Single candidate column '{candidates[0]}' selected as the analysis column."
        else:
            # Multiple columns; merge them
            merged_col = "merged_analysis_group"
            df[merged_col] = df[candidates].apply(
                lambda row: "_".join(row.values.astype(str)), axis=1)
            ctx.deps.metadata_df = df  # update the DataFrame
            ctx.deps.merged_column = merged_col
            msg = f"Multiple candidate columns {', '.join(candidates)} merged into column '{merged_col}'."

        log_tool_result(msg)
        return msg
    except Exception as e:
        error_msg = f"Error in merge_analysis_columns: {str(e)}"
        log(error_msg, style="bold red")
        log_tool_result(error_msg)
        return error_msg

# ----------------------------
# Tool: Design Contrasts (Simplified)
# ----------------------------


@rnaseq_agent.tool
async def design_contrasts(ctx: RunContext[RNAseqData]) -> str:
    """
    Identify the unique group values in the analysis column (either the merged column or the single candidate column)
    and return them as a simple comma-separated list. This tool does not compute contrasts but only reports
    the unique grouping values, leaving further contrast decisions to subsequent analysis.
    """
    try:
        if ctx.deps.metadata_df is None:
            return "Error: Metadata has not been processed."
        if not ctx.deps.merged_column:
            return "Error: Analysis column not defined. Please run merge_analysis_columns first."

        df = ctx.deps.metadata_df
        analysis_col = ctx.deps.merged_column

        # Extract unique values from the analysis column
        groups = df[analysis_col].dropna().unique().tolist()
        groups = [str(g) for g in groups]
        result = f"Unique groups in column '{analysis_col}': " + \
            ", ".join(groups)
        return result
    except Exception as e:
        error_msg = f"Error in design_contrasts: {str(e)}"
        return error_msg


# ----------------------------
# Main Execution
# ----------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="RNAseq metadata analysis pipeline with candidate selection, column merging, and unique group identification"
    )
    parser.add_argument("--metadata_path", type=str, default="metadata.csv",
                        help="Path to the metadata CSV file")
    args = parser.parse_args()

    analysis_data = RNAseqData(metadata_path=args.metadata_path)
    abs_path = os.path.abspath(args.metadata_path)
    print(f"Loading metadata from {abs_path}...")

    # Process metadata
    initial_prompt = "Please process the provided metadata. Identify the columns that will be most suitable for analysis, and subsequently determine the analysis contrasts that should be constructed."
    try:
        result = rnaseq_agent.run_sync(initial_prompt, deps=analysis_data, result_type = Contrasts)
        console.print(
            Panel("Metadata processing completed successfully!", style="bold green"))
        console.print("\n[bold yellow]Agent response:[/bold yellow]")
        console.print(result.data)
    except Exception as e:
        console.print(
            Panel(f"Error during metadata processing: {str(e)}", style="bold red"))
