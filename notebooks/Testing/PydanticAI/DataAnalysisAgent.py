# Create a global console instance
# console = Console()
# Imports
# ----------------------------
import os
import glob
import re
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import date
from rich.console import Console
from rich.panel import Panel
from dataclasses import dataclass
from typing import List, Dict, Optional, Union, Tuple, Any, Literal
from unidecode import unidecode
from pydantic import BaseModel, Field
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp
import nest_asyncio
nest_asyncio.apply()
console = Console()
from pydantic_ai import Agent, RunContext

# ----------------------------
# Define the dependency type
# ----------------------------
@dataclass
class RNAseqData:
    """Container for RNAseq analysis data and paths."""
    # Input data
    fastq_dir: str
    metadata_path: str
    kallisto_index_dir: str
    organism: str = "human"  # default to human
    output_dir: str = "output"
    tx2gene_path: Optional[str] = None

    # Runtime data that gets populated during analysis
    metadata_df: Optional[pd.DataFrame] = None
    abundance_files: List[str] = None
    merged_column: Optional[str] = None
    contrast_groups: Dict[str, Dict[str, str]] = None

# ----------------------------
# Create an RNAseq analysis agent
# ----------------------------

rnaseq_agent = Agent(
    'openai:gpt-4o-mini',
    deps_type=RNAseqData,
    system_prompt="""
    You are an expert RNAseq data analyst. Your task is to analyze RNAseq data using a series of bioinformatics tools.

    Follow these general principles:
    1. Work systematically through the RNA-seq analysis workflow
    2. Validate inputs at each step
    3. Provide clear explanations of what's happening
    4. Handle errors gracefully
    5. Generate appropriate visualizations when needed
    """
)

@rnaseq_agent.tool
async def run_gsea_analysis(ctx: RunContext[RNAseqData], contrast_name: str) -> str:
    """
    Run GSEA on expression data for a specific contrast.

    Args:
        contrast_name: Name of the contrast to analyze

    Returns:
        Summary of GSEA results
    """
    try:
        # Check if we have normalized expression data
        norm_counts_file = os.path.join(ctx.deps.output_dir, "DESeq2_normalized_counts.csv")
        if not os.path.exists(norm_counts_file):
            return "Error: Normalized counts file not found. Please run DESeq2 analysis first."

        # Load normalized counts
        expr_data = pd.read_csv(norm_counts_file, index_col=0)

        # Get contrast details
        contrast = ctx.deps.contrast_groups[contrast_name]
        group_a = contrast['numerator']
        group_b = contrast['denominator']

        # Create class vector
        class_vector = []
        for sample in expr_data.columns:
            if sample in ctx.deps.metadata_df[ctx.deps.metadata_df[ctx.deps.merged_column] == group_a].index:
                class_vector.append(group_a)
            else:
                class_vector.append(group_b)

        # Run GSEA
        gs_res = gp.gsea(
            data=expr_data,
            gene_sets='MSigDB_Hallmark_2020',
            cls=class_vector,
            permutation_type='phenotype',
            permutation_num=1000,
            outdir=os.path.join(ctx.deps.output_dir, f"gsea_{contrast_name}"),
            method='signal_to_noise',
            threads=4
        )

        # Generate plots
        terms = gs_res.res2d.Term
        gs_res.plot(terms=terms[:5],
                   show_ranking=True,
                   ofname=os.path.join(ctx.deps.output_dir, f"gsea_{contrast_name}_top5.png"))

        # Summarize results
        sig_pathways = gs_res.res2d[gs_res.res2d['FDR q-val'] < 0.25]
        return f"""
GSEA Analysis completed for contrast: {contrast_name}

Summary:
- Total pathways analyzed: {len(gs_res.res2d)}
- Significant pathways (FDR < 0.25): {len(sig_pathways)}
- Top enriched pathways:
{sig_pathways[['Term', 'NES', 'FDR q-val']].head().to_string()}

Generated files:
- GSEA results: gsea_{contrast_name}/
- Top pathways plot: gsea_{contrast_name}_top5.png
"""
    except Exception as e:
        return f"Error running GSEA analysis: {str(e)}"
    finally:
        plt.close('all')  # Ensure all plots are closed

@rnaseq_agent.tool
async def run_ssgsea_analysis(ctx: RunContext[RNAseqData]) -> str:
    """
    Run single sample GSEA analysis.

    Returns:
        Summary of ssGSEA results
    """
    try:
        # Check for normalized counts
        norm_counts_file = os.path.join(ctx.deps.output_dir, "DESeq2_normalized_counts.csv")
        if not os.path.exists(norm_counts_file):
            return "Error: Normalized counts file not found. Please run DESeq2 analysis first."

        # Run ssGSEA
        ss = gp.ssgsea(
            data=norm_counts_file,
            gene_sets='MSigDB_Hallmark_2020',
            outdir=os.path.join(ctx.deps.output_dir, "ssgsea"),
            sample_norm_method='rank',
            no_plot=False
        )

        # Create heatmap
        nes = ss.res2d.pivot(index='Term', columns='Name', values='NES')
        plt.figure(figsize=(12,8))
        sns.clustermap(nes,
                      cmap='RdBu_r',
                      center=0,
                      figsize=(12,8))
        plt.savefig(os.path.join(ctx.deps.output_dir, "ssgsea_heatmap.png"))
        plt.close()

        return f"""
Single Sample GSEA Analysis completed.

Generated files:
- ssGSEA results directory: ssgsea/
- Pathway score heatmap: ssgsea_heatmap.png

Number of pathways analyzed: {len(nes.index)}
Number of samples: {len(nes.columns)}
"""
    except Exception as e:
        return f"Error running ssGSEA analysis: {str(e)}"

@rnaseq_agent.tool
async def run_gsva_analysis(ctx: RunContext[RNAseqData]) -> str:
    """
    Run Gene Set Variation Analysis.

    Returns:
        Summary of GSVA results
    """
    try:
        # Check for normalized counts
        norm_counts_file = os.path.join(ctx.deps.output_dir, "DESeq2_normalized_counts.csv")
        if not os.path.exists(norm_counts_file):
            return "Error: Normalized counts file not found. Please run DESeq2 analysis first."

        # Run GSVA
        es = gp.gsva(
            data=norm_counts_file,
            gene_sets='MSigDB_Hallmark_2020',
            outdir=os.path.join(ctx.deps.output_dir, "gsva")
        )

        # Create visualization
        es_matrix = es.res2d.pivot(index='Term', columns='Name', values='ES')
        plt.figure(figsize=(12,8))
        sns.clustermap(es_matrix,
                      cmap='RdBu_r',
                      center=0,
                      figsize=(12,8))
        plt.savefig(os.path.join(ctx.deps.output_dir, "gsva_heatmap.png"))
        plt.close()

        return f"""
GSVA Analysis completed.

Generated files:
- GSVA results directory: gsva/
- Enrichment score heatmap: gsva_heatmap.png

Number of gene sets analyzed: {len(es_matrix.index)}
Number of samples: {len(es_matrix.columns)}
"""
    except Exception as e:
        return f"Error running GSVA analysis: {str(e)}"
# ----------------------------
# Utility Functions
# ----------------------------
@rnaseq_agent.tool
async def find_files(ctx: RunContext[RNAseqData], directory: str, suffix: str) -> List[str]:
    """
    Find files in a directory with a specific suffix.

    Args:
        directory: Directory to search in
        suffix: File suffix to look for (e.g., 'fastq.gz')

    Returns:
        List of absolute file paths matching the suffix
    """
    matched_files = []
    try:
        for root, _, files in os.walk(directory):
            for f in files:
                if f.endswith(suffix):
                    matched_files.append(os.path.join(root, f))
        return sorted(matched_files)
    except FileNotFoundError:
        return f"Error: Directory '{directory}' not found."
    except Exception as e:
        return f"Error: {str(e)}"

@rnaseq_agent.tool
async def load_metadata(ctx: RunContext[RNAseqData], metadata_path: str) -> str:
    """
    Load and validate the experiment metadata.

    Args:
        metadata_path: Path to the metadata file (CSV format expected)

    Returns:
        Description of the loaded metadata
    """
    try:
        # Try to load the metadata file
        if metadata_path.endswith('.csv'):
            metadata_df = pd.read_csv(metadata_path)
        elif metadata_path.endswith('.tsv') or metadata_path.endswith('.txt'):
            metadata_df = pd.read_csv(metadata_path, sep='\t')
        else:
            # Try to guess the format
            metadata_df = pd.read_csv(metadata_path, sep=None, engine='python')

        # Store in the context
        ctx.deps.metadata_df = metadata_df

        # Remove columns with all identical values as they're not useful for analysis
        useful_cols = metadata_df.loc[:, metadata_df.nunique() > 1].columns.tolist()

        # Provide information about the metadata
        return f"""
Metadata loaded successfully with {metadata_df.shape[0]} samples and {metadata_df.shape[1]} columns.
Columns with unique values: {', '.join(useful_cols)}
First 5 rows:
{metadata_df.head().to_string()}
        """
        return f"""
DESeq2 analysis completed successfully for contrast: {contrast_name}

Summary:
- Total genes analyzed: {len(results_df)}
- Significant genes (FDR < 0.05): {len(sig_results)}
- Upregulated in {numerator}: {len(sig_results[sig_results['log2FoldChange'] > 0])}
- Downregulated in {numerator}: {len(sig_results[sig_results['log2FoldChange'] < 0])}

Top 10 differentially expressed genes:
{sig_results.head(10)[['log2FoldChange', 'pvalue', 'padj']].to_string()}

Generated files:
- Results table: DESeq2_{contrast_name}_results.csv
- MA plot: DESeq2_{contrast_name}_MA_plot.png
- PCA plot: DESeq2_PCA_plot.png
- Heatmap: DESeq2_{contrast_name}_heatmap.png
- Volcano plot: DESeq2_{contrast_name}_volcano_plot.png
"""
    except Exception as e:
        return f"Error loading metadata: {str(e)}"

@rnaseq_agent.tool
async def clean_string(ctx: RunContext[RNAseqData], s: str) -> str:
    """
    Clean a string by normalizing special characters, replacing spaces with underscores,
    and removing non-word characters.

    Args:
        s: The string to clean

    Returns:
        The cleaned string
    """
    if pd.isna(s):
        return "NA"  # Handle missing values
    s = str(s).strip()  # Convert to string and remove leading/trailing whitespaces
    s = unidecode(s)  # Normalize special characters to ASCII
    s = s.replace(" ", "_")  # Replace spaces with underscores
    s = re.sub(r'[^\w]', '', s)  # Remove non-word characters
    return s

# ----------------------------
# Metadata Analysis Tools
# ----------------------------
@rnaseq_agent.tool
async def identify_analysis_columns(ctx: RunContext[RNAseqData]) -> str:
    """
    Identify which columns in the metadata should be used for differential expression analysis.

    Returns:
        Description of identified columns and recommendation
    """
    if ctx.deps.metadata_df is None:
        return "Error: Metadata not loaded. Please run load_metadata first."

    try:
        # Get columns with variability
        metadata_df = ctx.deps.metadata_df
        variable_cols = metadata_df.loc[:, metadata_df.nunique() > 1].columns.tolist()

        # Analyze potential biological factors
        biological_factors = []
        technical_factors = []

        # Common keywords for biological factors
        bio_keywords = ['treatment', 'condition', 'genotype', 'disease', 'cell', 'tissue',
                         'time', 'dose', 'age', 'sex', 'gender', 'strain', 'group']

        # Common keywords for technical factors
        tech_keywords = ['batch', 'run', 'lane', 'library', 'seq', 'date', 'id', 'rep']

        for col in variable_cols:
            col_lower = col.lower()
            is_bio = any(keyword in col_lower for keyword in bio_keywords)
            is_tech = any(keyword in col_lower for keyword in tech_keywords)

            if is_bio:
                biological_factors.append(col)
            elif is_tech:
                technical_factors.append(col)

        # Determine if columns should be merged
        merge_needed = len(biological_factors) > 1
        if merge_needed:
            recommendation = f"Multiple biological factors detected: {', '.join(biological_factors)}. Consider merging these columns for analysis."
            cols_to_merge = biological_factors
        else:
            if len(biological_factors) == 1:
                recommendation = f"One clear biological factor detected: {biological_factors[0]}. This can be used directly."
                cols_to_merge = [biological_factors[0]]
            else:
                # If no obvious biological factors, suggest the columns with the fewest unique values
                n_unique = {col: metadata_df[col].nunique() for col in variable_cols}
                sorted_cols = sorted(n_unique.items(), key=lambda x: x[1])
                cols_to_merge = [sorted_cols[0][0]]
                if len(sorted_cols) > 1 and sorted_cols[1][1] < 10:  # Only suggest merging if second column has few unique values
                    cols_to_merge.append(sorted_cols[1][0])
                    merge_needed = True
                    recommendation = f"No clear biological factors detected. Suggesting to merge columns with fewest unique values: {', '.join(cols_to_merge)}."
                else:
                    recommendation = f"No clear biological factors detected. Suggesting to use the column with fewest unique values: {cols_to_merge[0]}."

        # Store the recommendation in the context
        if merge_needed:
            ctx.deps.merged_column = "merged_analysis_group"
        else:
            ctx.deps.merged_column = cols_to_merge[0]

        return f"""
Analysis of metadata columns:
- Biological factors: {', '.join(biological_factors) if biological_factors else 'None clearly identified'}
- Technical factors: {', '.join(technical_factors) if technical_factors else 'None clearly identified'}
- Variable columns: {', '.join(variable_cols)}

Recommendation: {recommendation}
Columns to use: {', '.join(cols_to_merge)}
Merge needed: {merge_needed}
        """
    except Exception as e:
        return f"Error analyzing metadata columns: {str(e)}"

@rnaseq_agent.tool
async def merge_metadata_columns(ctx: RunContext[RNAseqData], columns: List[str], new_column_name: str = "merged_analysis_group") -> str:
    """
    Merge multiple metadata columns into a single column for analysis.

    Args:
        columns: List of column names to merge
        new_column_name: Name for the new merged column

    Returns:
        Description of the merged column
    """
    if ctx.deps.metadata_df is None:
        return "Error: Metadata not loaded. Please run load_metadata first."

    try:
        metadata_df = ctx.deps.metadata_df

        # Check that all columns exist
        missing_cols = [col for col in columns if col not in metadata_df.columns]
        if missing_cols:
            return f"Error: Columns not found in metadata: {', '.join(missing_cols)}"

        if len(columns) == 1:
            # Just rename the column if only one provided
            metadata_df[new_column_name] = metadata_df[columns[0]].apply(lambda x: clean_string(ctx, x))
            result_message = f"Renamed and cleaned column {columns[0]} to {new_column_name}."
        else:
            # Create the merged column by concatenating values with underscores
            metadata_df[new_column_name] = metadata_df[columns].apply(
                lambda row: '_'.join([clean_string(ctx, val) for val in row.values.astype(str)]),
                axis=1
            )
            result_message = f"Merged columns {', '.join(columns)} into new column {new_column_name}."

        # Update the metadata
        ctx.deps.metadata_df = metadata_df
        ctx.deps.merged_column = new_column_name

        # Get unique values in the merged column
        unique_values = metadata_df[new_column_name].unique().tolist()

        return f"""
{result_message}

The merged column '{new_column_name}' contains {len(unique_values)} unique values:
{', '.join(unique_values)}

Sample counts per group:
{metadata_df[new_column_name].value_counts().to_string()}
        """
    except Exception as e:
        return f"Error merging metadata columns: {str(e)}"

@rnaseq_agent.tool
async def design_contrasts(ctx: RunContext[RNAseqData]) -> str:
    """
    Design contrasts for differential expression analysis based on the metadata.

    Returns:
        Description of designed contrasts
    """
    if ctx.deps.metadata_df is None:
        return "Error: Metadata not loaded. Please run load_metadata first."

    if ctx.deps.merged_column is None:
        return "Error: Analysis column not identified. Please run identify_analysis_columns first."

    try:
        metadata_df = ctx.deps.metadata_df
        group_col = ctx.deps.merged_column

        # Get unique groups
        groups = metadata_df[group_col].unique().tolist()

        # Design contrasts based on the groups
        contrasts = []
        contrast_details = {}

        if len(groups) == 2:
            # Simple case with only two groups - one contrast
            contrasts.append(f"{groups[1]}_vs_{groups[0]}")
            contrast_details[f"{groups[1]}_vs_{groups[0]}"] = {
                'name': f"{groups[1]} vs {groups[0]}",
                'numerator': groups[1],
                'denominator': groups[0],
                'expression': f"{groups[1]} - {groups[0]}"
            }
        else:
            # For multiple groups, compare each to a reference
            # Try to identify a control/reference group
            potential_controls = [g for g in groups if any(kw in g.lower() for kw in
                                 ['control', 'ctrl', 'reference', 'ref', 'normal', 'wt', 'wild', 'mock', 'dmso', 'pbs', 'untreated'])]

            if potential_controls:
                # Use the identified control
                reference = potential_controls[0]
                other_groups = [g for g in groups if g != reference]

                for group in other_groups:
                    contrast_name = f"{group}_vs_{reference}"
                    contrasts.append(contrast_name)
                    contrast_details[contrast_name] = {
                        'name': f"{group} vs {reference}",
                        'numerator': group,
                        'denominator': reference,
                        'expression': f"{group} - {reference}"
                    }
            else:
                # If no obvious control, pick the first group as reference
                reference = groups[0]
                other_groups = groups[1:]

                for group in other_groups:
                    contrast_name = f"{group}_vs_{reference}"
                    contrasts.append(contrast_name)
                    contrast_details[contrast_name] = {
                        'name': f"{group} vs {reference}",
                        'numerator': group,
                        'denominator': reference,
                        'expression': f"{group} - {reference}"
                    }

                # Also create pairwise comparisons between non-reference groups if there aren't too many
                if len(other_groups) <= 4:  # Limit to avoid too many comparisons
                    for i, group1 in enumerate(other_groups):
                        for group2 in other_groups[i+1:]:
                            contrast_name = f"{group1}_vs_{group2}"
                            contrasts.append(contrast_name)
                            contrast_details[contrast_name] = {
                                'name': f"{group1} vs {group2}",
                                'numerator': group1,
                                'denominator': group2,
                                'expression': f"{group1} - {group2}"
                            }

        # Store the contrasts
        ctx.deps.contrast_groups = contrast_details

        return f"""
Designed contrasts based on {group_col} column with values: {', '.join(groups)}

Contrasts:
{pd.DataFrame([contrast_details[c] for c in contrasts]).to_string()}
        """
    except Exception as e:
        return f"Error designing contrasts: {str(e)}"

# ----------------------------
# Kallisto Quantification Tools
# ----------------------------
@rnaseq_agent.tool
async def find_kallisto_index(ctx: RunContext[RNAseqData]) -> str:
    """
    Find the appropriate Kallisto index based on the organism.

    Returns:
        Path to the Kallisto index file
    """
    try:
        organism = ctx.deps.organism.lower()
        index_dir = ctx.deps.kallisto_index_dir

        # Look for index files in the specified directory
        index_files = await find_files(ctx, index_dir, '.idx')

        if not index_files:
            return f"Error: No Kallisto index files found in {index_dir}"

        # Try to find an index matching the organism
        matching_indices = [idx for idx in index_files if organism in os.path.basename(os.path.dirname(idx)).lower()]

        if matching_indices:
            index_path = matching_indices[0]
            # Also find the transcript-to-gene mapping file if available
            tx2gene_files = await find_files(ctx, os.path.dirname(index_path), '.txt')
            if tx2gene_files:
                t2g_files = [f for f in tx2gene_files if any(x in os.path.basename(f).lower() for x in ['t2g', 'tx2gene'])]
                if t2g_files:
                    ctx.deps.tx2gene_path = t2g_files[0]

            return f"Found Kallisto index for {organism}: {index_path}"
        else:
            # If no organism-specific index found, return the first one
            return f"No index specific to {organism} found. Using the first available index: {index_files[0]}"

    except Exception as e:
        return f"Error finding Kallisto index: {str(e)}"

@rnaseq_agent.tool
async def run_kallisto_quantification(ctx: RunContext[RNAseqData]) -> str:
    """
    Run Kallisto quantification on the FASTQ files.

    Returns:
        Summary of the quantification process
    """
    try:
        # Find paired FASTQ files
        fastq_files = await find_files(ctx, ctx.deps.fastq_dir, 'fastq.gz')
        if not fastq_files:
            return f"Error: No FASTQ files found in {ctx.deps.fastq_dir}"

        # Find the Kallisto index
        index_result = await find_kallisto_index(ctx)
        if "Error" in index_result:
            return index_result

        # Extract the index path from the result
        index_path = None
        for line in index_result.splitlines():
            if '.idx' in line:
                # Extract the path, which should be after the colon
                if ':' in line:
                    index_path = line.split(':', 1)[1].strip()
                else:
                    # If no colon, look for a path with .idx
                    words = line.split()
                    for word in words:
                        if '.idx' in word:
                            index_path = word
                            break

        if not index_path:
            return "Error: Could not determine Kallisto index path"

        # Create output directory
        output_dir = ctx.deps.output_dir
        os.makedirs(output_dir, exist_ok=True)

        # Organize FASTQ files into pairs based on naming conventions
        paired_files = {}

        # Identify pairs using common naming patterns
        r1_pattern = re.compile(r'.*_(R1|1)\.fastq\.gz$')
        r2_pattern = re.compile(r'.*_(R2|2)\.fastq\.gz$')

        r1_files = [f for f in fastq_files if r1_pattern.match(f)]
        r2_files = [f for f in fastq_files if r2_pattern.match(f)]

        # Match R1 with R2 files
        for r1_file in r1_files:
            # Convert R1 to R2 in the filename
            expected_r2 = r1_file.replace('_R1', '_R2').replace('_1.fastq', '_2.fastq')
            if expected_r2 in r2_files:
                # Extract sample name from filename
                sample_name = os.path.basename(r1_file).split('_R1')[0].split('_1.fastq')[0]
                paired_files[sample_name] = (r1_file, expected_r2)

        if not paired_files:
            # If no pairs found, check if files are single-end
            single_end = all(not r1_pattern.match(f) and not r2_pattern.match(f) for f in fastq_files)
            if single_end:
                return "Error: Single-end reads detected. Kallisto requires paired-end reads or additional parameters for single-end analysis."
            else:
                return "Error: Could not identify paired FASTQ files"

        # Run Kallisto for each pair
        results = []
        for sample_name, (r1, r2) in paired_files.items():
            sample_output_dir = os.path.join(output_dir, sample_name)
            os.makedirs(sample_output_dir, exist_ok=True)

            # Build Kallisto command
            cmd = [
                "kallisto", "quant",
                "-i", index_path,
                "-o", sample_output_dir,
                "-t", "4",  # Use 4 threads
                "--plaintext",  # Output plaintext instead of HDF5
                "--bootstrap-samples=100",  # Number of bootstrap samples
                r1, r2
            ]

            # Run Kallisto
            process = subprocess.run(cmd, capture_output=True, text=True)

            if process.returncode == 0:
                results.append(f"Successfully processed {sample_name}")
            else:
                results.append(f"Error processing {sample_name}: {process.stderr}")

        # Collect paths to abundance files
        abundance_files = []
        for sample_name in paired_files.keys():
            abundance_file = os.path.join(output_dir, sample_name, "abundance.tsv")
            if os.path.exists(abundance_file):
                abundance_files.append(abundance_file)

        # Store the abundance file paths
        ctx.deps.abundance_files = abundance_files

        return f"""
Kallisto quantification completed for {len(results)} sample pairs.

Results:
{chr(10).join(results)}

Found {len(abundance_files)} abundance files for downstream analysis.
        """
    except Exception as e:
        return f"Error running Kallisto quantification: {str(e)}"

# ----------------------------
# Differential Expression Analysis Tools
# ----------------------------
@rnaseq_agent.tool
async def prepare_deseq2_analysis(ctx: RunContext[RNAseqData]) -> str:
    """
    Prepare data for DESeq2 analysis by integrating Kallisto quantification results with metadata.

    Returns:
        Description of the prepared data
    """
    try:
        # Check if we have abundance files
        if not ctx.deps.abundance_files:
            return "Error: No abundance files found. Please run Kallisto quantification first."

        # Check if we have metadata
        if ctx.deps.metadata_df is None:
            return "Error: Metadata not loaded. Please run load_metadata first."

        # Check if we have merged column
        if ctx.deps.merged_column is None:
            return "Error: Analysis column not identified. Please run identify_analysis_columns first."

        # Get sample names from abundance file paths
        sample_names = [os.path.basename(os.path.dirname(f)) for f in ctx.deps.abundance_files]

        # Create a mapping between abundance files and metadata
        # First, try to match based on exact sample names
        metadata_df = ctx.deps.metadata_df
        matched_samples = {}
        unmatched_samples = []

        for i, sample_name in enumerate(sample_names):
            # Try different ways to match samples
            # 1. Exact match in any column
            exact_matches = []
            for col in metadata_df.columns:
                matches = metadata_df[metadata_df[col].astype(str) == sample_name].index.tolist()
                exact_matches.extend(matches)

            if exact_matches:
                matched_samples[sample_name] = {
                    'abundance_file': ctx.deps.abundance_files[i],
                    'metadata_row': exact_matches[0]
                }
            else:
                # 2. Try to find a substring match
                substring_matches = []
                for col in metadata_df.columns:
                    # Check if the sample name is contained in any column value
                    matches = metadata_df[metadata_df[col].astype(str).str.contains(sample_name, case=False, na=False)].index.tolist()
                    substring_matches.extend(matches)

                if substring_matches:
                    matched_samples[sample_name] = {
                        'abundance_file': ctx.deps.abundance_files[i],
                        'metadata_row': substring_matches[0]
                    }
                else:
                    unmatched_samples.append(sample_name)

        if unmatched_samples:
            # Try reverse matching - look for metadata values in sample names
            for sample_name in unmatched_samples[:]:
                for col in metadata_df.columns:
                    for idx, value in metadata_df[col].items():
                        if str(value) in sample_name:
                            matched_samples[sample_name] = {
                                'abundance_file': ctx.deps.abundance_files[sample_names.index(sample_name)],
                                'metadata_row': idx
                            }
                            unmatched_samples.remove(sample_name)
                            break
                    if sample_name not in unmatched_samples:
                        break

        # If there are still unmatched samples, try a more aggressive approach
        if unmatched_samples and len(unmatched_samples) < len(sample_names) // 2:
            # If we have matched most samples, we can guess the rest based on order
            if len(metadata_df) == len(sample_names):
                sorted_samples = sorted(sample_names)
                sorted_metadata = metadata_df.sort_index()

                for i, sample_name in enumerate(sorted_samples):
                    if sample_name in unmatched_samples:
                        matched_samples[sample_name] = {
                            'abundance_file': ctx.deps.abundance_files[sample_names.index(sample_name)],
                            'metadata_row': sorted_metadata.index[i]
                        }
                unmatched_samples = []

        if unmatched_samples:
            return f"""
Warning: Could not match {len(unmatched_samples)} of {len(sample_names)} samples to metadata.
Unmatched samples: {', '.join(unmatched_samples)}

Please check that sample names in the FASTQ files correspond to identifiers in the metadata.
            """

        # Create a DataFrame for DESeq2 analysis
        analysis_df = pd.DataFrame(index=list(matched_samples.keys()))

        # Add the file paths
        analysis_df['abundance_file'] = [matched_samples[s]['abundance_file'] for s in analysis_df.index]

        # Add the metadata
        for col in metadata_df.columns:
            analysis_df[col] = [metadata_df.loc[matched_samples[s]['metadata_row'], col] for s in analysis_df.index]

        # Save the analysis dataframe for later use
        analysis_df.to_csv(os.path.join(ctx.deps.output_dir, "deseq2_analysis_samples.csv"))

        return f"""
Successfully prepared data for DESeq2 analysis with {len(analysis_df)} samples.

Sample mapping:
{analysis_df[['abundance_file', ctx.deps.merged_column]].head().to_string()}
{' ... ' if len(analysis_df) > 5 else ''}
{analysis_df[['abundance_file', ctx.deps.merged_column]].tail().to_string() if len(analysis_df) > 5 else ''}

Group counts:
{analysis_df[ctx.deps.merged_column].value_counts().to_string()}

Analysis is ready to proceed with the following groups: {', '.join(analysis_df[ctx.deps.merged_column].unique())}
        """
    except Exception as e:
        return f"Error preparing DESeq2 analysis: {str(e)}"

@rnaseq_agent.tool
async def run_deseq2_analysis(ctx: RunContext[RNAseqData], contrast_name: str) -> str:
    """
    Run DESeq2 differential expression analysis for a specific contrast.

    Args:
        contrast_name: Name of the contrast to analyze

    Returns:
        Summary of the differential expression results
    """
    try:
        # Check if we have the contrast definition
        if not ctx.deps.contrast_groups or contrast_name not in ctx.deps.contrast_groups:
            return f"Error: Contrast '{contrast_name}' not defined. Please run design_contrasts first."

        # Get contrast details
        contrast = ctx.deps.contrast_groups[contrast_name]
        numerator = contrast['numerator']
        denominator = contrast['denominator']

        # Check if we have the merged column
        if ctx.deps.merged_column is None:
            return "Error: Analysis column not identified. Please run identify_analysis_columns first."

        # Load the sample mapping file
        sample_file = os.path.join(ctx.deps.output_dir, "deseq2_analysis_samples.csv")
        if not os.path.exists(sample_file):
            return "Error: Sample mapping file not found. Please run prepare_deseq2_analysis first."

        sample_df = pd.read_csv(sample_file, index_col=0)

        # Create the R script for DESeq2 analysis
        r_script = f"""
        library(tximport)
        library(DESeq2)
        library(ggplot2)
        library(pheatmap)
        library(dplyr)
        library(tibble)
        library(tidyr)

        # Set working directory
        setwd("{ctx.deps.output_dir}")

        # Load sample information
        samples <- read.csv("{sample_file}", row.names=1)

        # Get file paths
        files <- samples$abundance_file

        # Get sample groups using the merged analysis column
        group <- factor(samples[["{ctx.deps.merged_column}"]])
        coldata <- data.frame(row.names=rownames(samples), group=group)

        # Check if tx2gene file exists
        tx2gene_exists <- {str(ctx.deps.tx2gene_path is not None).upper()}

        if (tx2gene_exists) {{
          # Import transcript to gene mapping
          tx2gene <- read.table("{ctx.deps.tx2gene_path}", sep="\\t", header=FALSE)
          tx2gene <- tx2gene[, c(1, 2)]
          txi <- tximport(files, type="kallisto", tx2gene=tx2gene, ignoreAfterBar=TRUE)
          dds <- DESeqDataSetFromTximport(txi, colData=coldata, design=~group)
        }} else {{
          txi <- tximport(files, type="kallisto", txOut=TRUE)
          dds <- DESeqDataSetFromTximport(txi, colData=coldata, design=~group)
        }}

        # Filter low count genes
        keep <- rowSums(counts(dds)) >= 10
        dds <- dds[keep,]

        # Set reference level
        dds$group <- relevel(dds$group, ref="{denominator}")
        dds <- DESeq(dds)
        res <- results(dds, contrast=c("group", "{numerator}", "{denominator}"))
        res <- res[order(res$padj),]
        write.csv(as.data.frame(res), file="DESeq2_{contrast_name}_results.csv")
        normalized_counts <- counts(dds, normalized=TRUE)
        write.csv(normalized_counts, file="DESeq2_normalized_counts.csv")

        png("DESeq2_{contrast_name}_MA_plot.png", width=800, height=600)
        plotMA(res, main="{contrast_name} MA Plot", ylim=c(-5,5))
        dev.off()

        png("DESeq2_PCA_plot.png", width=800, height=600)
        vsd <- vst(dds, blind=FALSE)
        plotPCA(vsd, intgroup=c("group")) + theme_bw() + ggtitle("PCA of samples")
        dev.off()

        png("DESeq2_{contrast_name}_heatmap.png", width=800, height=800)
        mat <- assay(vsd)
        topgenes <- head(rownames(res), 50)
        mat <- mat[topgenes,]
        mat <- mat - rowMeans(mat)
        pheatmap(mat, main="Top 50 DE genes")
        dev.off()

        png("DESeq2_{contrast_name}_volcano_plot.png", width=800, height=600)
        res_df <- as.data.frame(res)
        res_df$gene <- rownames(res_df)
        res_df$sig <- "Not Significant"
        res_df$sig[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Upregulated"
        res_df$sig[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Downregulated"
        ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=sig)) +
          geom_point(alpha=0.7) +
          scale_color_manual(values=c("blue", "gray", "red")) +
          theme_bw() +
          geom_vline(xintercept=c(-1,1), linetype="dashed") +
          geom_hline(yintercept=-log10(0.05), linetype="dashed") +
          labs(title="Volcano Plot", x="Log2 Fold Change", y="-Log10 Adjusted p-value") +
          theme(legend.title=element_blank())
        dev.off()

        writeLines(capture.output(sessionInfo()), "DESeq2_session_info.txt")
        """

        # Write the R script to a file
        r_script_path = os.path.join(ctx.deps.output_dir, f"run_deseq2_{contrast_name}.R")
        with open(r_script_path, "w") as f:
            f.write(r_script)

        # Run the R script
        process = subprocess.run(["Rscript", r_script_path], capture_output=True, text=True)

        if process.returncode != 0:
            return f"Error running DESeq2 analysis: {process.stderr}"

        # Check if the results file was created
        results_file = os.path.join(ctx.deps.output_dir, f"DESeq2_{contrast_name}_results.csv")
        if not os.path.exists(results_file):
            return "Error: DESeq2 analysis did not produce results file"

        # Load the results for summary
        results_df = pd.read_csv(results_file, index_col=0)

        # Filter for significant genes
        sig_results = results_df[results_df['padj'] < 0.05].sort_values('padj')

        return f"""
        DESeq2 analysis completed successfully for contrast: {contrast_name}

        Summary:
        - Total genes analyzed: {len(results_df)}
        - Significant genes (FDR < 0.05): {len(sig_results)}
        - Upregulated in {numerator}: {len(sig_results[sig_results['log2FoldChange'] > 0])}
        - Downregulated in {numerator}: {len(sig_results[sig_results['log2FoldChange'] < 0])}

        Top 10 differentially expressed genes:
        {sig_results.head(10)[['log2FoldChange', 'pvalue', 'padj']].to_string()}

        Generated files:
        - Results table: DESeq2_{contrast_name}_results.csv
        - MA plot: DESeq2_{contrast_name}_MA_plot.png
        - PCA plot: DESeq2_PCA_plot.png
        - Heatmap: DESeq2_{contrast_name}_heatmap.png
        - Volcano plot: DESeq2_{contrast_name}_volcano_plot.png
        """
    except Exception as e:
        return f"Error running DESeq2 analysis: {str(e)}"
# ----------------------------
# Main Execution
# ----------------------------
if __name__ == "__main__":
    # Create data instance for GSE262710
    test_data = RNAseqData(
        fastq_dir="./TestRNAseqData_SETBP1/GSE262710/fastq",
        metadata_path="./TestRNAseqData_SETBP1/GSE262710/GSE262710_metadata.csv",
        kallisto_index_dir="../../../data/kallisto_indices",
        organism="human",
        output_dir="analysis_output/GSE262710"
    )

    # Initialize conversation with analysis steps
    initial_prompt = """
    Please analyze the RNA-seq data with the following steps:
    1. Run Kallisto quantification on all paired-end FASTQ files
    2. Load and analyze the metadata to identify experimental groups
    3. Merge appropriate columns for analysis if needed
    4. Design appropriate contrasts for differential expression
    5. Run DESeq2 analysis for each contrast
    6. Perform pathway analysis using:
       - Standard GSEA
       - Single-sample GSEA
       - GSVA
    7. Generate appropriate visualizations

    Please provide updates at each step.
    """

    # Run the agent
    try:
        # Run the agent synchronously
        result = rnaseq_agent.run_sync(initial_prompt, deps=test_data)
        console.print(Panel("Analysis completed successfully!", style="bold green"))
        console.print("\n[bold yellow]Agent response:[/bold yellow]")
        console.print(result.data)
    except Exception as e:
        console.print(Panel(f"Error during analysis: {str(e)}", style="bold red"))
