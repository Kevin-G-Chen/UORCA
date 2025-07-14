"""
UORCA - Unified Omics Reference Corpus of Analyses

A fully containerized, AI-powered workflow for automated RNA-seq analysis
of public datasets from the Gene Expression Omnibus (GEO).
"""

__version__ = "0.1.0"
__author__ = "UORCA Development Team"
__description__ = "Automated RNA-seq analysis of public datasets"

# Main CLI entry point
from .cli import main

__all__ = ["main"]
