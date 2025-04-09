#!/usr/bin/env python3
"""
reporting_agent_sphinx_minimal.py

A minimal reporting agent that builds Sphinx HTML documentation from RST files.
It is designed to work with the RST file (e.g. deg_report.rst) produced by your DEG reporting step.

Usage:
  python reporting_agent_sphinx_minimal.py --rst_folder /path/to/rst_files \
    --output_folder /path/to/sphinx_project --log_path sphinx_build.log
"""

import os
import glob
import shutil
import subprocess
import argparse
import logging
from pydantic import BaseModel, Field
from pydantic_ai import Agent, RunContext
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()
OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")

# ------------------------------------
# Minimal Sphinx Builder Functions
# ------------------------------------


def append_line_to_file(file_path, line):
    """
    Append a line to the end of the specified file.

    Args:
        file_path (str): The path to the file.
        line (str): The line to append.
    """
    with open(file_path, 'a') as file:
        file.write(f"{line}\n")


def generate_index_rst(source_dir):
    """
    Generate an index.rst file for the Sphinx documentation.

    This function collects all .rst files in the source directory (excluding index.rst)
    and builds a simple toctree.

    Args:
        source_dir (str): The source directory of the Sphinx project.
    """
    logging.info("Generating index.rst from RST files...")
    # Collect all .rst files, excluding index.rst
    rst_files = sorted(glob.glob(os.path.join(source_dir, "*.rst")))
    rst_files = [f for f in rst_files if os.path.basename(f) != "index.rst"]
    # Get file names without the .rst extension.
    doc_files = [os.path.basename(f)[:-4] for f in rst_files]
    logging.info(f"Found document files: {doc_files}")

    index_content = [
        "DEG Report Documentation\n",
        "========================\n\n",
        ".. toctree::\n",
        "   :maxdepth: 1\n",
        "\n",
    ]
    for doc in doc_files:
        index_content.append(f"   {doc}\n")

    index_path = os.path.join(source_dir, "index.rst")
    with open(index_path, "w") as f:
        f.writelines(index_content)
    logging.info(f"Generated index.rst at {index_path}")


def build_sphinx_docs(rst_folder, output_folder, log_path):
    """
    Build Sphinx documentation from RST files.

    This function creates a new Sphinx project, copies over the RST files, generates an index.rst,
    and then builds the HTML documentation.

    Args:
        rst_folder (str): Folder containing the reStructuredText (.rst) files.
        output_folder (str): Folder to create the Sphinx project.
        log_path (str): Path to write the build log.

    Returns:
        bool: True if the build succeeded, False otherwise.
    """
    logging.info("Starting minimal Sphinx documentation build process...")
    try:
        # Ensure the output project folder exists.
        os.makedirs(output_folder, exist_ok=True)

        # Create a new Sphinx project (in quiet mode).
        subprocess.run([
            "sphinx-quickstart",
            output_folder,
            "--quiet",
            "--project", "DEG_Report",
            "--author", "",
            "--sep",
        ], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        logging.info(f"Created new Sphinx project in {output_folder}")

        # Define the source directory of the Sphinx project.
        source_dir = os.path.join(output_folder, "source")

        # Copy all RST files from rst_folder to the Sphinx source directory.
        for rst_file in glob.glob(os.path.join(rst_folder, "*.rst")):
            shutil.copy(rst_file, source_dir)
        logging.info(f"Copied RST files from {rst_folder} to {source_dir}")

        # Generate the index.rst file.
        generate_index_rst(source_dir)

        # Build the HTML documentation using sphinx-build.
        build_dir = os.path.join(output_folder, "build", "html")
        result = subprocess.run([
            "sphinx-build",
            "-b", "html",
            source_dir,
            build_dir,
        ], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Write build output and errors to the log file.
        with open(log_path, 'w') as f:
            f.write(result.stdout.decode())
            f.write("\n\n")
            f.write(result.stderr.decode())

        logging.info(f"Built HTML documentation at {build_dir}")
        logging.info(f"Build log saved to {log_path}")
        return True

    except subprocess.CalledProcessError as e:
        logging.error(f"Sphinx build error: {e}")
        with open(log_path, 'w') as f:
            f.write(f"Error: {str(e)}\n")
            if e.stdout:
                f.write("\nSTDOUT:\n")
                f.write(e.stdout.decode())
            if e.stderr:
                f.write("\nSTDERR:\n")
                f.write(e.stderr.decode())
        return False

    except Exception as e:
        logging.error(f"Unexpected error during Sphinx build: {e}")
        with open(log_path, 'w') as f:
            f.write(f"Unexpected error: {str(e)}\n")
        return False

# ------------------------------------
# Pydantic Model for the Reporting Agent
# ------------------------------------


class SphinxReportContext(BaseModel):
    rst_folder: str = Field(...,
                            description="Folder containing reStructuredText (.rst) files.")
    output_folder: str = Field(
        ..., description="Folder for the Sphinx project and HTML build output.")
    log_path: str = Field(...,
                          description="Path to the Sphinx build log file.")


# ------------------------------------
# Create the Reporting Agent using PydanticAI.
# ------------------------------------
sphinx_agent = Agent(
    'openai:gpt-4o',
    deps_type=SphinxReportContext,
    system_prompt="Minimal reporting agent that builds Sphinx documentation from RST files."
)

# ------------------------------------
# Define the Tool to Build the Report
# ------------------------------------


@sphinx_agent.tool
async def build_report(ctx: RunContext[SphinxReportContext]) -> str:
    """
    Build Sphinx documentation from the provided RST files.

    Process:
      1. Create a new Sphinx project.
      2. Copy RST files from the given folder to the project.
      3. Generate a minimal index.rst with a toctree.
      4. Build HTML documentation.

    Returns:
      A message indicating success or failure.
    """
    success = build_sphinx_docs(
        rst_folder=ctx.deps.rst_folder,
        output_folder=ctx.deps.output_folder,
        log_path=ctx.deps.log_path
    )
    if success:
        html_dir = os.path.join(ctx.deps.output_folder, "build", "html")
        return f"Sphinx HTML documentation built successfully in {html_dir} (log: {ctx.deps.log_path})"
    else:
        return "Error: Sphinx documentation build failed. Check the log for details."

# ------------------------------------
# Main function for standalone execution.
# ------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Minimal reporting agent that builds Sphinx HTML documentation from RST files."
    )
    parser.add_argument("--rst_folder", required=True,
                        help="Folder containing RST files (e.g., deg_report.rst)")
    parser.add_argument("--output_folder", required=True,
                        help="Output folder to create the Sphinx project")
    parser.add_argument("--log_path", default="sphinx_build.log",
                        help="Path to save the Sphinx build log")
    args = parser.parse_args()

    # Configure logging (prints timestamps to the console).
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.StreamHandler()]
    )

    # Create the dependency context for the agent.
    deps = SphinxReportContext(
        rst_folder=args.rst_folder,
        output_folder=args.output_folder,
        log_path=args.log_path
    )

    result = sphinx_agent.run_sync("Build Sphinx docs", deps=deps)
    print(result.data)


if __name__ == "__main__":
    main()
