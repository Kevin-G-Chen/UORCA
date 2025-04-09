#!/usr/bin/env python
"""
reporting_agent.py

A minimal reporting agent that takes in an analysis output directory,
locates the CSV file produced by differential expression analysis (DEGs),
and generates a Sphinx documentation report (in reStructuredText format).

Usage:
    python reporting_agent.py --output_dir /path/to/analysis_output
"""

import os
import glob
import argparse
import pandas as pd
from tabulate import tabulate
from dotenv import load_dotenv

from pydantic import BaseModel, Field
from pydantic_ai import Agent, RunContext

load_dotenv()
openai_api_key = os.getenv("OPENAI_API_KEY")

# ------------------------------------------------------------------------------
# Define a dependency model for the reporting agent.
# This model holds the output directory in which the analysis results (e.g. DEG CSV) are stored.
# ------------------------------------------------------------------------------


class ReportContext(BaseModel):
    output_dir: str = Field(...,
                            description="Directory where analysis results are stored.")


# ------------------------------------------------------------------------------
# Create the reporting agent using PydanticAI.
# The agent's system prompt is kept simple; you can expand this if needed.
# ------------------------------------------------------------------------------
reporting_agent = Agent(
    'openai:gpt-4o',
    deps_type=ReportContext,
    system_prompt="Reporting agent that generates Sphinx documentation from DEG CSV file."
)

# ------------------------------------------------------------------------------
# Define a tool to generate the DEG report as a Sphinx (reST) document.
# ------------------------------------------------------------------------------


@reporting_agent.tool
async def generate_deg_report(ctx: RunContext[ReportContext]) -> str:
    """
    Generate a Sphinx documentation report based on the DEG CSV file.

    Process:
      1. Locate the CSV file matching 'DEG_results*.csv' in the specified output directory.
      2. Load the CSV file into a pandas DataFrame.
      3. Convert the DataFrame to a reStructuredText formatted table using tabulate.
      4. Assemble a Sphinx-compliant reST document with a title, brief analysis summary, and the table.
      5. Write the report to 'deg_report.rst' in the output directory.

    Returns:
      A string message with the path to the generated report file.
    """
    # Locate CSV files that match the DEG results filename pattern
    csv_files = glob.glob(os.path.join(
        ctx.deps.output_dir, "DEG_results*.csv"))
    if not csv_files:
        return "No DEG CSV file found in the output directory."

    # For a minimal example, use the first found CSV file.
    deg_csv = csv_files[0]
    df = pd.read_csv(deg_csv)
    # Only keep the first 10 rows for the report
    df = df.head(10)

    # Convert the DataFrame to an reStructuredText formatted table using tabulate
    table_rst = tabulate(df, headers='keys', tablefmt='rst', showindex=False)

    # Assemble the Sphinx documentation text (in reST format)
    report_lines = [
        "Differentially Expressed Genes (DEG) Report",
        "===========================================",
        "",
        "This report documents the results of the differential expression analysis.",
        "",
        "The table below shows a summary of the top differentially expressed genes:",
        "",
        table_rst,
        "",
        "Generated on: " + pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"),
        ""
    ]
    report_text = "\n".join(report_lines)

    # Write the report to a file in the output directory
    report_file = os.path.join(ctx.deps.output_dir, "deg_report.rst")
    with open(report_file, "w") as f:
        f.write(report_text)

    return f"Report generated: {report_file}"

# ------------------------------------------------------------------------------
# Main function to run the reporting agent.
# This parses command-line arguments and runs the agent synchronously.
# ------------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Reporting agent that generates Sphinx documentation from DEG CSV."
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Directory where analysis results are stored."
    )
    args = parser.parse_args()

    # Create the context (dependencies) for the agent
    report_context = ReportContext(output_dir=args.output_dir)

    # Run the reporting agent synchronously; the agent will execute its tools.
    result = reporting_agent.run_sync(
        "Generate DEG report", deps=report_context)
    print(result.data)


# ------------------------------------------------------------------------------
# Script entry point
# ------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
