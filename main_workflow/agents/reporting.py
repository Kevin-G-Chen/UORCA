#!/usr/bin/env python3
"""
reporting_agent_two_parameters_logging.py

A reporting agent that:
  1. Identifies PNG files in a user-specified PNG input folder.
  2. Copies these PNG files to an "images" subfolder within the output RST folder.
  3. Generates a reStructuredText (RST) file that embeds these copied PNG images
     (using a static caption for every figure).
  4. Builds Sphinx HTML documentation from the generated RST file(s).

Output structure:
  - RST files are saved in <output_dir>/rst.
  - A timestamped subfolder is created inside <output_dir> to host the Sphinx project and build log.
  - The Sphinx build log (sphinx_build.log) is saved in the timestamped folder.
  - All images are copied to <output_dir>/rst/images, and the RST file references those copies.

Usage:
  python reporting_agent_two_parameters_logging.py --png_dir path/to/png_folder --output_dir path/to/output_dir

Make sure your .env file (with OPENAI_API_KEY, etc.) is available.
"""

import os
import glob
import shutil
import subprocess
import argparse
import logging
import datetime
from pydantic import BaseModel, Field, ConfigDict
from pydantic_ai import Agent, RunContext
from dotenv import load_dotenv
from shared import ReportContext

# Load environment variables
load_dotenv()
openai_api_key = os.getenv("OPENAI_API_KEY")

# -----------------------------------------------------------
# Reporting Agent Tools (with static captions and logging)
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
    'openai:o4-mini',  # Using a powerful model.
    deps_type=ReportContext,
    system_prompt=system_prompt
)


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
    Build Sphinx HTML documentation from the RST file in rst_folder.
    Initialize a Sphinx project in sphinx_output_folder, copy the RST file (and associated images),
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
            "--project", "Image_Report",
            "--author", "Reporting Agent",
            "--sep"
        ]
        subprocess.run(quickstart_cmd, check=True,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        logging.info("[build_report] Sphinx project created successfully.")

        # The Sphinx source directory is inside the sphinx project folder.
        source_dir = os.path.join(sphinx_project, "source")
        # Copy RST files from rst_folder (including report_images.rst) to source_dir.
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

async def run_agent_async(prompt: str, deps: ReportContext, usage=None):
    return await reporting_agent.run(prompt, deps=deps, usage=usage)
