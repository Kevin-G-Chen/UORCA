#!/usr/bin/env python3
"""
ReportingAgent_enhanced.py

An enhanced reporting agent that:
  1. Identifies PNG files in a user-specified PNG input folder.
  2. Analyzes and classifies these PNG files by type and purpose.
  3. Organizes them into a logical narrative flow for reporting.
  4. Generates a reStructuredText (RST) file that embeds these PNG images with meaningful titles and descriptions.
  5. Builds Sphinx HTML documentation from the generated RST file(s).

Output structure:
  - RST files are saved in <output_dir>/rst.
  - A timestamped subfolder is created inside <output_dir> to host the Sphinx project and build log.
  - The Sphinx build log (sphinx_build.log) is saved in the timestamped folder.
  - All images are copied to <output_dir>/rst/images, and the RST file references those copies.

Usage:
  python ReportingAgent_enhanced.py --png_dir path/to/png_folder --output_dir path/to/output_dir [--title "Report Title"] [--author "Author Name"]

Make sure your .env file (with OPENAI_API_KEY, etc.) is available.
"""

import os
import glob
import shutil
import subprocess
import argparse
import logging
import datetime
import re
from typing import List, Dict, Optional
from dataclasses import dataclass
from pydantic import BaseModel, Field, ConfigDict
from pydantic_ai import Agent, RunContext
from dotenv import load_dotenv

# Load environment variables
load_dotenv()
openai_api_key = os.getenv("OPENAI_API_KEY")

# -----------------------------------------------------------
# Image Classification Functionality
# -----------------------------------------------------------


@dataclass
class ImageMetadata:
    """Metadata extracted from an image file"""
    filepath: str
    filename: str
    category: str  # e.g., 'volcano_plot', 'heatmap', 'pca', etc.
    analysis_type: str  # e.g., 'quality_control', 'differential_expression'
    priority: int  # Lower number = higher priority in displaying
    # For plots associated with a specific contrast
    contrast_name: Optional[str] = None
    suggested_title: str = ""
    suggested_description: str = ""


def classify_image(image_path: str) -> ImageMetadata:
    """
    Analyze a PNG file and classify it based on filename patterns.
    This function uses regex patterns to identify common plot types in RNA-seq analyses.
    """
    filename = os.path.basename(image_path)
    lowercase_filename = filename.lower()

    # Default values
    category = "other"
    analysis_type = "unknown"
    priority = 100
    contrast_name = None
    suggested_title = filename  # Default to filename

    # Extract contrast name if present (common pattern in RNA-seq analysis)
    contrast_match = re.search(
        r'(?:contrast[_-]|^)([a-zA-Z0-9_-]+)(?:[_-]vs[_-]|[_-]vs$)([a-zA-Z0-9_-]+)', lowercase_filename)
    if contrast_match:
        contrast_name = f"{contrast_match.group(1)} vs {contrast_match.group(2)}"

    # Classify by plot type
    if any(x in lowercase_filename for x in ['pca', 'principal_component']):
        category = "pca_plot"
        analysis_type = "exploratory_analysis"
        priority = 10
        suggested_title = "Principal Component Analysis (PCA)"

    elif any(x in lowercase_filename for x in ['volcano']):
        category = "volcano_plot"
        analysis_type = "differential_expression"
        priority = 30
        suggested_title = "Volcano Plot of Differential Expression"

    elif any(x in lowercase_filename for x in ['heatmap', 'heat_map']):
        category = "heatmap"
        analysis_type = "differential_expression"
        priority = 40
        suggested_title = "Heatmap of Differential Expression"

    elif any(x in lowercase_filename for x in ['gsea', 'gene_set', 'pathway']):
        category = "pathway_analysis"
        analysis_type = "functional_analysis"
        priority = 50
        suggested_title = "Gene Set Enrichment Analysis"

    elif any(x in lowercase_filename for x in ['ma_plot', 'maplot']):
        category = "ma_plot"
        analysis_type = "differential_expression"
        priority = 35
        suggested_title = "MA Plot"

    elif any(x in lowercase_filename for x in ['boxplot', 'box_plot']):
        category = "boxplot"
        analysis_type = "exploratory_analysis"
        priority = 20
        suggested_title = "Expression Distribution Boxplot"

    elif any(x in lowercase_filename for x in ['qc', 'quality']):
        category = "quality_control"
        analysis_type = "quality_control"
        priority = 5
        suggested_title = "Quality Control Metrics"

    elif any(x in lowercase_filename for x in ['count', 'distribution']):
        category = "count_distribution"
        analysis_type = "exploratory_analysis"
        priority = 15
        suggested_title = "Read Count Distribution"

    # If contrast name is available, add it to the title
    if contrast_name and "vs" not in suggested_title:
        suggested_title += f" ({contrast_name})"

    return ImageMetadata(
        filepath=image_path,
        filename=filename,
        category=category,
        analysis_type=analysis_type,
        priority=priority,
        contrast_name=contrast_name,
        suggested_title=suggested_title,
        suggested_description=""  # Will be filled by generate_image_interpretations
    )

# -----------------------------------------------------------
# Dependency Class for Reporting Agent
# -----------------------------------------------------------


class ReportContext(BaseModel):
    """
    Enhanced configuration for the reporting agent.

    - png_dir: Input directory for PNG files.
    - rst_folder: Directory (inside output_dir/rst) where the RST file and copied images are saved.
    - sphinx_output_folder: Timestamped folder (inside output_dir) where the Sphinx project and HTML output are stored.
    - log_path: Full path to the Sphinx build log file.
    - title: Title for the generated report.
    - author: Author name for the generated report.
    - organized_images: Optional list of classified and organized image metadata.
    """
    png_dir: str = Field(...,
                         description="Directory where the PNG files are located.")
    rst_folder: str = Field(...,
                            description="Directory to store the generated RST file and images.")
    sphinx_output_folder: str = Field(
        ..., description="Directory for the Sphinx project and HTML output.")
    log_path: str = Field(...,
                          description="Path to save the Sphinx build log file.")
    title: str = Field("RNA-seq Analysis Report",
                       description="Title for the report")
    author: str = Field("Reporting Agent",
                        description="Author name for the report")
    organized_images: Optional[List[ImageMetadata]] = Field(
        None, description="Organized list of image metadata")
    model_config = ConfigDict(extra="allow")


# -----------------------------------------------------------
# Reporting Agent Tools
# -----------------------------------------------------------
system_prompt = """
You are an advanced reporting agent specialized in creating comprehensive RNA-seq analysis reports from PNG plots.
Your expertise includes:
1. Understanding different types of RNA-seq analysis plots (PCA, volcano plots, heatmaps, etc.)
2. Organizing plots into a logical narrative flow
3. Creating insightful titles and descriptions for each plot
4. Generating publication-quality reStructuredText documents

For each plot, provide:
1. A meaningful title based on the plot's content and purpose
2. An appropriate description explaining what the plot shows and its significance
3. Proper placement within the overall flow of the report

Organize the report following standard RNA-seq analysis practices:
1. Start with quality control and sample overview plots
2. Follow with exploratory data analysis
3. Continue with differential expression results
4. End with pathway/enrichment analyses

Always aim to create a coherent narrative that guides the reader through the analysis results.
"""

reporting_agent = Agent(
    'openai:gpt-4o',
    deps_type=ReportContext,
    system_prompt=system_prompt
)


@reporting_agent.tool
async def identify_and_classify_png_files(ctx: RunContext[ReportContext]) -> str:
    """
    Identify all PNG files in the folder and classify them by type, content, and purpose.
    Returns detailed information about the discovered images.
    """
    png_dir = ctx.deps.png_dir
    logging.info(
        f"[identify_and_classify_png_files] Scanning PNG input folder: {png_dir}")

    if not os.path.isdir(png_dir):
        err_msg = f"Error: PNG input folder '{png_dir}' does not exist."
        logging.error(err_msg)
        return err_msg

    png_files = []
    for root, _, files in os.walk(png_dir):
        for f in files:
            if f.lower().endswith('.png'):
                png_files.append(os.path.join(root, f))
    png_files = sorted(png_files)

    if not png_files:
        msg = f"No PNG files found in {png_dir}."
        logging.warning(msg)
        return msg

    # Classify each image and extract metadata
    classified_images = [classify_image(png_path) for png_path in png_files]

    # Store the classified images in the context for later use
    ctx.deps.organized_images = classified_images

    # Generate a summary of the classified images
    result = f"Found and classified {len(classified_images)} PNG file(s) in {png_dir}:\n"

    for img in classified_images:
        result += f"\n- File: {img.filename}"
        result += f"\n  Category: {img.category}"
        result += f"\n  Analysis Type: {img.analysis_type}"
        if img.contrast_name:
            result += f"\n  Contrast: {img.contrast_name}"

    logging.info(f"[identify_and_classify_png_files] {result}")
    return result


@reporting_agent.tool
async def organize_images_for_report(ctx: RunContext[ReportContext]) -> str:
    """
    Organize the images into a logical flow for the report based on their classification.
    This determines the order in which images will appear and which ones might be excluded.
    """
    if not ctx.deps.organized_images:
        error_msg = "Error: No classified images available. Please run identify_and_classify_png_files first."
        logging.error(error_msg)
        return error_msg

    # First, sort by analysis type and priority
    images = sorted(
        ctx.deps.organized_images,
        key=lambda x: (x.analysis_type, x.priority)
    )

    # Update the organized_images in the context
    ctx.deps.organized_images = images

    # Create a summary of the organization
    result = "Images have been organized in the following sequence for the report:\n\n"

    for i, img in enumerate(images, 1):
        result += f"{i}. {img.filename} ({img.category}, {img.analysis_type})\n"
        result += f"   Suggested title: {img.suggested_title}\n"

    logging.info(
        f"[organize_images_for_report] Organized {len(images)} images")
    return result


@reporting_agent.tool
async def generate_image_interpretations(ctx: RunContext[ReportContext]) -> str:
    """
    Generate meaningful titles and descriptions for each image based on its classification.
    This helps provide context and explanations in the final report.
    """
    if not ctx.deps.organized_images:
        error_msg = "Error: No organized images available. Please run organize_images_for_report first."
        logging.error(error_msg)
        return error_msg

    result = "Generated interpretations for images:\n\n"

    for i, img in enumerate(ctx.deps.organized_images):
        # Here, we would typically use more complex logic to generate interpretations
        # But for now, we'll use filename patterns to create basic interpretations

        filename = img.filename
        category = img.category

        # Generate title based on category and filename
        if "volcano" in filename.lower():
            title = f"Volcano Plot: Differential Expression Analysis"
            if img.contrast_name:
                title += f" for {img.contrast_name}"
            description = (
                "This volcano plot shows the distribution of differentially expressed genes. "
                "The x-axis represents the log2 fold change, while the y-axis shows the -log10 p-value. "
                "Points in the upper left and right corners represent genes with statistically significant "
                "differential expression."
            )
        elif "heatmap" in filename.lower():
            title = "Heatmap of Differentially Expressed Genes"
            if img.contrast_name:
                title += f" for {img.contrast_name}"
            description = (
                "This heatmap visualizes expression patterns across samples for the top differentially "
                "expressed genes. Rows represent genes, columns represent samples, and colors indicate "
                "expression levels (typically red for high expression, blue for low expression)."
            )
        elif "pca" in filename.lower():
            title = "Principal Component Analysis (PCA) of Samples"
            description = (
                "This PCA plot shows the relationship between samples based on their gene expression profiles. "
                "Samples clustering together have similar expression patterns, while samples far apart have "
                "more divergent profiles."
            )
        else:
            title = f"{category.replace('_', ' ').title()} Plot"
            description = f"This plot shows {category.replace('_', ' ')} data from the RNA-seq analysis."

        # Update the image metadata
        ctx.deps.organized_images[i].suggested_title = title
        ctx.deps.organized_images[i].suggested_description = description

        result += f"{i+1}. {filename}\n"
        result += f"   Title: {title}\n"
        result += f"   Description: {description[:100]}...\n\n"

    logging.info(
        f"[generate_image_interpretations] Generated interpretations for {len(ctx.deps.organized_images)} images")
    return result


@reporting_agent.tool
async def generate_enhanced_rst(ctx: RunContext[ReportContext]) -> str:
    """
    Generate an enhanced RST file that embeds the PNG images in a logical order with proper titles and descriptions.
    """
    rst_folder = ctx.deps.rst_folder
    images_dir = os.path.join(rst_folder, "images")
    logging.info(f"[generate_enhanced_rst] RST folder: {rst_folder}")

    if not ctx.deps.organized_images:
        error_msg = "Error: No organized images available. Please run generate_image_interpretations first."
        logging.error(error_msg)
        return error_msg

    # Create output directories
    os.makedirs(rst_folder, exist_ok=True)
    os.makedirs(images_dir, exist_ok=True)

    # Copy images into the local images directory
    logging.info(
        "[generate_enhanced_rst] Copying PNG files into local images directory...")
    copied_images = []

    for img_metadata in ctx.deps.organized_images:
        try:
            dest_path = os.path.join(
                images_dir, os.path.basename(img_metadata.filepath))
            shutil.copy(img_metadata.filepath, dest_path)
            copied_images.append((dest_path, img_metadata))
            logging.debug(f"Copied {img_metadata.filepath} to {dest_path}")
        except Exception as e:
            logging.error(f"Error copying {img_metadata.filepath}: {str(e)}")

    # Generate RST content with enhanced structure and organization
    rst_lines = [
        f"{ctx.deps.title}",
        "=" * len(ctx.deps.title),
        "",
        f":Author: {ctx.deps.author}",
        f":Date: {datetime.datetime.now().strftime('%Y-%m-%d')}",
        "",
        "Executive Summary",
        "-----------------",
        "",
        "This report presents the results of an RNA-seq data analysis, including quality control, ",
        "exploratory data analysis, differential expression, and pathway analysis.",
        "",
        "Contents",
        "--------",
        "",
        ".. contents:: Table of Contents",
        "   :depth: 2",
        "",
    ]

    # Group images by analysis type
    analysis_types = {}
    for img_path, img_metadata in copied_images:
        analysis_type = img_metadata.analysis_type
        if analysis_type not in analysis_types:
            analysis_types[analysis_type] = []
        analysis_types[analysis_type].append((img_path, img_metadata))

    # Add sections for each analysis type
    for analysis_type, images in sorted(analysis_types.items()):
        section_title = analysis_type.replace('_', ' ').title()
        rst_lines.extend([
            f"{section_title}",
            "-" * len(section_title),
            "",
        ])

        # Add images for this section
        for img_path, img_metadata in images:
            relative_path = os.path.relpath(img_path, rst_folder)
            rst_lines.extend([
                f".. _fig_{img_metadata.filename.replace('.', '_')}:",
                "",
                f".. figure:: {relative_path}",
                "   :align: center",
                "   :scale: 80%",
                "",
                f"   {img_metadata.suggested_title}",
                "",
                f"{img_metadata.suggested_description}",
                "",
                "----",
                "",
            ])

    # Write the RST file
    rst_text = "\n".join(rst_lines)
    rst_file = os.path.join(rst_folder, "analysis_report.rst")

    try:
        with open(rst_file, "w") as f:
            f.write(rst_text)
        msg = f"Enhanced RST file generated successfully: {rst_file}"
        logging.info(f"[generate_enhanced_rst] {msg}")
        return msg
    except Exception as e:
        err_msg = f"Error writing RST file: {str(e)}"
        logging.error(f"[generate_enhanced_rst] {err_msg}")
        return err_msg


# Keep the original identify_png_files function as a fallback
@reporting_agent.tool
async def identify_png_files(ctx: RunContext[ReportContext]) -> str:
    """
    Identify and list all PNG files in the folder specified by ctx.deps.png_dir.
    Logs the directory being scanned and the files found.
    """
    png_dir = ctx.deps.png_dir
    logging.info(f"[identify_png_files] Scanning PNG input folder: {png_dir}")

    if not os.path.isdir(png_dir):
        err_msg = f"Error: PNG input folder '{png_dir}' does not exist."
        logging.error(err_msg)
        return err_msg

    png_files = []
    for root, _, files in os.walk(png_dir):
        for f in files:
            if f.lower().endswith('.png'):
                png_files.append(os.path.join(root, f))
    png_files = sorted(png_files)

    if not png_files:
        msg = f"No PNG files found in {png_dir}."
        logging.warning(msg)
        return msg

    result = f"Found {len(png_files)} PNG file(s) in {png_dir}:\n" + \
        "\n".join(png_files)
    logging.info(f"[identify_png_files] " + result)
    return result


# Keep the original generate_rst_from_pngs as a fallback
@reporting_agent.tool
async def generate_rst_from_pngs(ctx: RunContext[ReportContext]) -> str:
    """
    Copy PNG files into an 'images' subfolder in the rst_folder,
    then generate an RST file that embeds these images using the .. figure:: directive
    with a static caption.
    """
    png_dir = ctx.deps.png_dir
    rst_folder = ctx.deps.rst_folder
    images_dir = os.path.join(rst_folder, "images")
    logging.info(f"[generate_rst_from_pngs] Using PNG folder: {png_dir}")
    logging.info(f"[generate_rst_from_pngs] RST folder: {rst_folder}")

    if not os.path.isdir(png_dir):
        err_msg = f"Error: PNG input folder '{png_dir}' does not exist."
        logging.error(err_msg)
        return err_msg

    os.makedirs(rst_folder, exist_ok=True)
    os.makedirs(images_dir, exist_ok=True)

    # Find PNG files
    png_files = glob.glob(os.path.join(png_dir, "**", "*.png"), recursive=True)
    png_files = sorted(png_files)

    if not png_files:
        err_msg = f"Error: No PNG files found in {png_dir}."
        logging.error(err_msg)
        return err_msg

    # Copy images into the local images directory
    logging.info(
        "[generate_rst_from_pngs] Copying PNG files into local images directory...")
    for png_path in png_files:
        try:
            dest_path = os.path.join(images_dir, os.path.basename(png_path))
            shutil.copy(png_path, dest_path)
            logging.debug(f"Copied {png_path} to {dest_path}")
        except Exception as e:
            logging.error(f"Error copying {png_path}: {str(e)}")

    # Generate RST content; reference images from the 'images' folder.
    rst_lines = [
        "Image Report",
        "============",
        "",
        "This document embeds the following PNG images:",
        ""
    ]

    static_caption = "Figure: Placeholder Caption"
    # Now iterate over the copied files in images_dir
    copied_images = glob.glob(os.path.join(images_dir, "*.png"))
    copied_images = sorted(copied_images)
    for img_path in copied_images:
        relative_path = os.path.relpath(img_path, rst_folder)
        rst_lines.append(f".. figure:: {relative_path}")
        rst_lines.append("   :scale: 50%")
        rst_lines.append("")
        rst_lines.append(f"   {static_caption}")
        rst_lines.append("")

    rst_text = "\n".join(rst_lines)
    rst_file = os.path.join(rst_folder, "report_images.rst")

    try:
        with open(rst_file, "w") as f:
            f.write(rst_text)
        msg = f"RST file generated successfully: {rst_file}"
        logging.info(f"[generate_rst_from_pngs] {msg}")
        return msg
    except Exception as e:
        err_msg = f"Error writing RST file: {str(e)}"
        logging.error(f"[generate_rst_from_pngs] {err_msg}")
        return err_msg


@reporting_agent.tool
async def build_report(ctx: RunContext[ReportContext]) -> str:
    """
    Build Sphinx HTML documentation from the RST files in rst_folder.
    Initialize a Sphinx project in sphinx_output_folder, copy the RST files (and associated images),
    generate an index.rst, and run sphinx-build. Logs each step.
    """
    rst_folder = ctx.deps.rst_folder
    sphinx_project = ctx.deps.sphinx_output_folder
    log_path = ctx.deps.log_path
    logging.info(
        f"[build_report] Initializing Sphinx project in: {sphinx_project}")

    try:
        os.makedirs(sphinx_project, exist_ok=True)
        quickstart_cmd = [
            "sphinx-quickstart",
            sphinx_project,
            "--quiet",
            "--project", ctx.deps.title.replace(" ", "_"),
            "--author", ctx.deps.author,
            "--sep"
        ]
        subprocess.run(quickstart_cmd, check=True,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        logging.info("[build_report] Sphinx project created successfully.")

        # The Sphinx source directory is inside the sphinx project folder.
        source_dir = os.path.join(sphinx_project, "source")
        # Copy RST files from rst_folder to source_dir.
        rst_files = glob.glob(os.path.join(rst_folder, "*.rst"))
        if not rst_files:
            err_msg = f"Error: No RST files found in {rst_folder}"
            logging.error(f"[build_report] {err_msg}")
            return err_msg

        for rst_file in rst_files:
            shutil.copy(rst_file, source_dir)
            logging.debug(
                f"[build_report] Copied RST file {rst_file} to source directory.")
        # Also copy the local images folder if it exists.
        images_src = os.path.join(rst_folder, "images")
        images_dest = os.path.join(source_dir, "images")
        if os.path.isdir(images_src):
            shutil.copytree(images_src, images_dest, dirs_exist_ok=True)
            logging.info(
                f"[build_report] Copied images from {images_src} to {images_dest}.")

        # Generate an index.rst in the source directory.
        index_content = [
            f"{ctx.deps.title}",
            "=" * len(ctx.deps.title),
            "",
            ".. toctree::",
            "   :maxdepth: 2",
            ""
        ]
        for rst_file in glob.glob(os.path.join(source_dir, "*.rst")):
            basename = os.path.basename(rst_file)
            if basename != "index.rst":
                doc_name = os.path.splitext(basename)[0]
                index_content.append(f"   {doc_name}")
        index_text = "\n".join(index_content)
        index_path = os.path.join(source_dir, "index.rst")
        with open(index_path, "w") as f:
            f.write(index_text)
        logging.info(f"[build_report] Generated index.rst at {index_path}")

        # Build the HTML documentation using sphinx-build.
        build_dir = os.path.join(sphinx_project, "build", "html")
        build_cmd = [
            "sphinx-build",
            "-b", "html",
            source_dir,
            build_dir
        ]
        result = subprocess.run(build_cmd, check=True,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        with open(log_path, "w") as log_file:
            log_file.write(result.stdout.decode())
            log_file.write("\n")
            log_file.write(result.stderr.decode())
        msg = f"Sphinx HTML documentation built successfully in {build_dir} (log: {log_path})"
        logging.info(f"[build_report] {msg}")
        return msg
    except subprocess.CalledProcessError as e:
        error_text = f"Sphinx build error: {str(e)}\nSTDOUT: {e.stdout.decode()}\nSTDERR: {e.stderr.decode()}"
        with open(log_path, "w") as log_file:
            log_file.write(error_text)
        logging.error(f"[build_report] {error_text}")
        return error_text
    except Exception as e:
        err_msg = f"Unexpected error during Sphinx build: {str(e)}"
        logging.error(f"[build_report] {err_msg}")
        return err_msg


# -----------------------------------------------------------
# Main Execution
# -----------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Enhanced Reporting Agent: Convert PNG files to an RST report with intelligent organization and descriptions."
    )
    parser.add_argument("--png_dir", required=True,
                        help="Input directory where PNG files are stored.")
    parser.add_argument("--output_dir", required=True,
                        help="Base output directory for all outputs (RST, report, and log).")
    parser.add_argument("--title", default="RNA-seq Analysis Report",
                        help="Title for the generated report.")
    parser.add_argument("--author", default="Reporting Agent",
                        help="Author name for the generated report.")
    args = parser.parse_args()

    # Setup directories and logging
    base_output_dir = os.path.abspath(args.output_dir)
    rst_dir = os.path.join(base_output_dir, "rst")
    os.makedirs(rst_dir, exist_ok=True)

    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    report_dir = os.path.join(base_output_dir, timestamp)
    os.makedirs(report_dir, exist_ok=True)

    sphinx_project_dir = report_dir
    log_file_path = os.path.join(report_dir, "sphinx_build.log")

    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.StreamHandler()]
    )

    logging.info("=== Enhanced Reporting Agent Started ===")
    logging.info(f"PNG input directory: {args.png_dir}")
    logging.info(f"Base output directory: {base_output_dir}")
    logging.info(f"RST output directory: {rst_dir}")
    logging.info(f"Report title: {args.title}")
    logging.info(f"Report author: {args.author}")

    deps = ReportContext(
        png_dir=os.path.abspath(args.png_dir),
        rst_folder=rst_dir,
        sphinx_output_folder=sphinx_project_dir,
        log_path=log_file_path,
        title=args.title,
        author=args.author
    )

    initial_prompt = """
    Please perform an enhanced reporting workflow:

    1. First, identify and classify all PNG files in the provided directory.
    2. Then organize these images into a logical sequence for the report.
    3. Generate appropriate titles and descriptions for each image.
    4. Create an enhanced RST document with proper organization and context.
    5. Build the Sphinx HTML documentation from the generated RST.

    Please ensure that the final report presents a coherent narrative flow with images
    organized by analysis type and with meaningful interpretations.
    """

    try:
        result = reporting_agent.run_sync(initial_prompt, deps=deps)
        logging.info("=== Enhanced Reporting Agent Completed ===")
        print("Agent Response:")
        print(result.data)
    except Exception as e:
        logging.error(f"Error during reporting: {str(e)}")
        print(f"Error during reporting: {str(e)}")


if __name__ == "__main__":
    main()
