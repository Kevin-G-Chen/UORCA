from __future__ import annotations
import os, glob, shutil, subprocess, logging, datetime
from typing import List, Optional, Dict, Any, Union

from pydantic_ai import Agent, RunContext
from dotenv import load_dotenv
from shared import ReportingContext
from shared.workflow_logging import log_tool

# Setup logging
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

# -----------------------------------------------------------
# Reporting Agent Definition
# -----------------------------------------------------------
system_prompt = """
You are a reporting agent tasked with converting PNG images into a reStructuredText document and then building a Sphinx HTML report.
For this proof of concept, every image should display the static caption 'Figure: Placeholder Caption'.
Your tasks are:
1. Identify the PNG files in the provided png_dir.
2. Copy these PNG files to a local 'images' subfolder within the RST folder.
3. Generate an RST file that embeds these PNG images (using the figure directive with the static caption).
4. Build Sphinx HTML documentation from the generated RST file(s).
Return a summary of each step.
"""

reporting_agent = Agent(
    'openai:o4-mini',
    deps_type=ReportingContext,
    system_prompt=system_prompt
)

# -----------------------------------------------------------
# Reporting Agent Tools
# -----------------------------------------------------------
@reporting_agent.tool
@log_tool
async def identify_png_files(ctx: RunContext[ReportingContext]) -> str:
    """
    Identify and list all PNG files in the folder specified by ctx.deps.png_dir.
    """
    png_dir = ctx.deps.png_dir
    logger.info("🔍 Scanning PNG input folder: %s", png_dir)

    if not os.path.isdir(png_dir):
        error_msg = f"Error: PNG input folder '{png_dir}' does not exist."
        logger.error("❌ %s", error_msg)
        return error_msg

    png_files = []
    for root, _, files in os.walk(png_dir):
        for f in files:
            if f.lower().endswith('.png'):
                png_files.append(os.path.join(root, f))
    png_files = sorted(png_files)

    if not png_files:
        msg = f"No PNG files found in {png_dir}."
        logger.warning("⚠️ %s", msg)
        return msg

    result = f"Found {len(png_files)} PNG file(s) in {png_dir}"
    if len(png_files) <= 5:
        logger.info("✅ %s: %s", result, png_files)
    else:
        logger.info("✅ %s (first 5): %s...", result, png_files[:5])

    return result


@reporting_agent.tool
@log_tool
async def generate_rst_from_pngs(ctx: RunContext[ReportingContext]) -> str:
    """
    Copy PNG files into an 'images' subfolder in the rst_folder,
    then generate an RST file that embeds these images using the .. figure:: directive
    with a static caption.
    """
    png_dir = ctx.deps.png_dir
    rst_folder = ctx.deps.rst_folder
    images_dir = os.path.join(rst_folder, "images")
    logger.info("🖼️ Preparing to generate RST from PNGs - source: %s, destination: %s", png_dir, rst_folder)

    if not os.path.isdir(png_dir):
        error_msg = f"Error: PNG input folder '{png_dir}' does not exist."
        logger.error("❌ %s", error_msg)
        return error_msg

    os.makedirs(rst_folder, exist_ok=True)
    os.makedirs(images_dir, exist_ok=True)
    logger.info("📁 Created output directories: %s, %s", rst_folder, images_dir)

    # Find PNG files
    png_files = glob.glob(os.path.join(png_dir, "**", "*.png"), recursive=True)
    png_files = sorted(png_files)

    if not png_files:
        error_msg = f"Error: No PNG files found in {png_dir}."
        logger.error("❌ %s", error_msg)
        return error_msg

    # Copy images into the local images directory
    logger.info("📋 Copying %d PNG files to local images directory", len(png_files))
    for png_path in png_files:
        try:
            dest_path = os.path.join(images_dir, os.path.basename(png_path))
            shutil.copy(png_path, dest_path)
            logger.debug("📄 Copied %s to %s", png_path, dest_path)
        except Exception as e:
            logger.error("❌ Error copying %s: %s", png_path, str(e), exc_info=True)

    # Generate RST content; reference images from the 'images' folder.
    logger.info("📝 Generating RST content")
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
        logger.info("✅ %s", msg)
        return msg
    except Exception as e:
        error_msg = f"Error writing RST file: {str(e)}"
        logger.error("❌ %s", error_msg, exc_info=True)
        return error_msg


@reporting_agent.tool
@log_tool
async def build_report(ctx: RunContext[ReportingContext]) -> str:
    """
    Build Sphinx HTML documentation from the RST file in rst_folder.
    Initialize a Sphinx project in sphinx_output_folder, copy the RST file (and associated images),
    generate an index.rst, and run sphinx-build.
    """
    rst_folder = ctx.deps.rst_folder
    sphinx_project = ctx.deps.sphinx_output_folder
    log_path = ctx.deps.log_path
    logger.info("🏗️ Initializing Sphinx project in: %s", sphinx_project)

    try:
        os.makedirs(sphinx_project, exist_ok=True)
        quickstart_cmd = [
            "sphinx-quickstart",
            sphinx_project,
            "--quiet",
            "--project", "Image_Report",
            "--author", "Reporting Agent",
            "--sep"
        ]
        logger.info("⚙️ Running sphinx-quickstart: %s", " ".join(quickstart_cmd))

        result = subprocess.run(quickstart_cmd, check=True,
                      stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        logger.info("✅ Sphinx project created successfully")

        # The Sphinx source directory is inside the sphinx project folder
        source_dir = os.path.join(sphinx_project, "source")

        # Copy RST files from rst_folder to source_dir
        rst_files = glob.glob(os.path.join(rst_folder, "*.rst"))
        if not rst_files:
            error_msg = f"Error: No RST files found in {rst_folder}"
            logger.error("❌ %s", error_msg)
            return error_msg

        for rst_file in rst_files:
            shutil.copy(rst_file, source_dir)
            logger.debug("📄 Copied RST file %s to source directory", rst_file)

        # Copy the local images folder if it exists
        images_src = os.path.join(rst_folder, "images")
        images_dest = os.path.join(source_dir, "images")
        if os.path.isdir(images_src):
            shutil.copytree(images_src, images_dest, dirs_exist_ok=True)
            logger.info("📁 Copied images from %s to %s", images_src, images_dest)

        # Generate an index.rst in the source directory
        logger.info("📝 Generating index.rst")
        index_content = [
            "Image Report Documentation",
            "============================",
            "",
            ".. toctree::",
            "   :maxdepth: 1",
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
        logger.info("✅ Generated index.rst at %s", index_path)

        # Build the HTML documentation using sphinx-build
        build_dir = os.path.join(sphinx_project, "build", "html")
        build_cmd = [
            "sphinx-build",
            "-b", "html",
            source_dir,
            build_dir
        ]

        logger.info("🔨 Building Sphinx documentation: %s", " ".join(build_cmd))
        result = subprocess.run(build_cmd, check=True,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        with open(log_path, "w") as log_file:
            log_file.write(result.stdout.decode())
            log_file.write("\n")
            log_file.write(result.stderr.decode())

        msg = f"Sphinx HTML documentation built successfully in {build_dir} (log: {log_path})"
        logger.info("✅ %s", msg)
        return msg

    except subprocess.CalledProcessError as e:
        error_text = f"Sphinx build error: {str(e)}\nSTDOUT: {e.stdout.decode()}\nSTDERR: {e.stderr.decode()}"
        with open(log_path, "w") as log_file:
            log_file.write(error_text)
        logger.error("❌ %s", error_text, exc_info=True)
        return error_text

    except Exception as e:
        error_msg = f"Unexpected error during Sphinx build: {str(e)}"
        logger.error("❌ %s", error_msg, exc_info=True)
        return error_msg


@log_tool
async def run_agent_async(prompt: str, deps: ReportingContext, usage=None):
    """Thin wrapper used by master.py to invoke the reporting agent asynchronously."""
    logger.info("📝 Reporting agent invoked by master – prompt: %s", prompt)
    return await reporting_agent.run(prompt, deps=deps, usage=usage)
