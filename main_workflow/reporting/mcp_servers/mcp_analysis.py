#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MCP Analysis Server for RNA-seq Analysis

This server provides higher-level analysis tools for RNA-seq data, building on
the raw data extraction capabilities. It focuses on pattern recognition,
statistical analysis across datasets, and relevance scoring.

Tools provided:
- rank_contrasts_by_relevance: Rank contrasts by relevance to research question
- analyze_gene_patterns: Analyze expression patterns across genes
- find_expression_clusters: Find clusters of co-expressed genes
- calculate_gene_statistics: Calculate summary statistics for genes across contrasts
- identify_consistent_degs: Find genes with consistent expression changes
"""

import os
import sys
import json
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist
from typing import List, Dict, Optional, Tuple, Any
from pathlib import Path
import logging
from collections import defaultdict
import re
import argparse

from mcp.server.fastmcp import FastMCP

# Initialize MCP server
mcp = FastMCP("analysis")
log = lambda m: print(f"[ANALYSIS] {m}", file=sys.stderr, flush=True)

# Set up logging
logging.basicConfig(
    level=logging.WARNING,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    stream=sys.stderr
)
logger = logging.getLogger("mcp-analysis")

# Add a decorator to log tool execution times
def log_tool_execution(func):
    async def wrapper(*args, **kwargs):
        import time
        logger.info(f"Starting tool: {func.__name__}")
        start_time = time.time()
        try:
            result = await func(*args, **kwargs)
            logger.info(f"Tool {func.__name__} completed in {time.time() - start_time:.2f}s")
            return result
        except Exception as e:
            logger.error(f"Tool {func.__name__} failed after {time.time() - start_time:.2f}s: {e}")
            raise
    return wrapper

# Global variables
RESULTS_DIR = os.environ.get("RESULTS_DIR", "")
_cache = {}  # Simple cache for repeated queries

def _calculate_text_similarity(text1: str, text2: str) -> float:
    """Calculate simple text similarity based on shared terms."""
    # Convert to lowercase and split into words
    words1 = set(re.findall(r'\w+', text1.lower()))
    words2 = set(re.findall(r'\w+', text2.lower()))
    
    # Remove common words
    stopwords = {'the', 'a', 'an', 'and', 'or', 'but', 'in', 'on', 'at', 'to', 'for',
                 'of', 'with', 'by', 'from', 'vs', 'versus', 'compared'}
    words1 = words1 - stopwords
    words2 = words2 - stopwords
    
    if not words1 or not words2:
        return 0.0
    
    # Calculate Jaccard similarity
    intersection = words1 & words2
    union = words1 | words2
    
    return len(intersection) / len(union) if union else 0.0

def _load_all_deg_data() -> Dict[str, pd.DataFrame]:
    """Load all DEG data from the results directory."""
    deg_data = {}
    
    for root, dirs, files in os.walk(RESULTS_DIR):
        if "DEG.csv" in files and "RNAseqAnalysis" in root:
            try:
                # Extract identifiers from path
                path_parts = root.split(os.sep)
                analysis_id = None
                contrast_id = None
                
                for i, part in enumerate(path_parts):
                    if part == "RNAseqAnalysis" and i > 0:
                        analysis_id = path_parts[i-1]
                        if i + 1 < len(path_parts):
                            contrast_id = path_parts[i+1]
                        break
                
                if analysis_id and contrast_id:
                    df = pd.read_csv(os.path.join(root, "DEG.csv"))
                    # Standardize gene column
                    for col in ['Gene', 'gene', 'SYMBOL', 'symbol']:
                        if col in df.columns:
                            df = df.rename(columns={col: 'Gene'})
                            break
                    
                    key = f"{analysis_id}:{contrast_id}"
                    deg_data[key] = df
                    
            except Exception as e:
                logger.warning(f"Error loading DEG file from {root}: {e}")
    
    return deg_data

def _load_contrast_descriptions() -> Dict[str, str]:
    """Load contrast descriptions from metadata files."""
    descriptions = {}
    
    # Look for contrasts.csv files
    for root, dirs, files in os.walk(RESULTS_DIR):
        if "contrasts.csv" in files:
            try:
                contrasts_df = pd.read_csv(os.path.join(root, "contrasts.csv"))
                for _, row in contrasts_df.iterrows():
                    if 'name' in row and 'description' in row:
                        descriptions[row['name']] = row['description']
            except Exception as e:
                logger.warning(f"Error loading contrasts from {root}: {e}")
    
    return descriptions

@mcp.tool()
@log_tool_execution
async def rank_contrasts_by_relevance(research_question: str,
                                     max_contrasts: int = 10,
                                     min_deg_count: int = 10) -> str:
    """
    Rank contrasts by relevance to a specified research question.
    
    Args:
        research_question: The research question to analyze.
        max_contrasts: Maximum number of contrasts to return.
        min_deg_count: Minimum number of DEGs required for a contrast to be considered.
    
    Returns:
        JSON string with ranked contrasts and relevance scores.
    """
    try:
        all_contrasts = []
        deg_data = _load_all_deg_data()
        descriptions = _load_contrast_descriptions()
        
        # Analyze each contrast
        for contrast_key, df in deg_data.items():
            analysis_id, contrast_id = contrast_key.split(":", 1)
            
            # Count DEGs
            deg_count = 0
            if 'adj.P.Val' in df.columns and 'logFC' in df.columns:
                deg_count = ((df['adj.P.Val'] < 0.05) & (abs(df['logFC']) > 1.0)).sum()
            elif 'P.Value' in df.columns and 'logFC' in df.columns:
                deg_count = ((df['P.Value'] < 0.05) & (abs(df['logFC']) > 1.0)).sum()
            
            if deg_count < min_deg_count:
                continue
            
            # Get description
            description = descriptions.get(contrast_id, f"Contrast: {contrast_id}")
            
            # Calculate relevance score based on text similarity
            relevance_score = _calculate_text_similarity(research_question, description)
            
            # Boost score based on DEG count (normalized)
            deg_boost = min(1.0, deg_count / 1000.0) * 0.3  # Up to 30% boost
            relevance_score = relevance_score * (1 + deg_boost)
            
            # Look for specific keywords that might indicate relevance
            keywords_in_question = set(re.findall(r'\w+', research_question.lower()))
            keywords_in_description = set(re.findall(r'\w+', description.lower()))
            
            # Check for exact matches of important terms (3+ characters)
            important_matches = [k for k in keywords_in_question & keywords_in_description 
                               if len(k) > 3]
            if important_matches:
                relevance_score += 0.2 * len(important_matches)
            
            all_contrasts.append({
                "analysis_id": analysis_id,
                "contrast_id": contrast_id,
                "description": description,
                "deg_count": int(deg_count),
                "relevance_score": round(relevance_score, 3),
                "matching_terms": important_matches
            })
        
        # Sort by relevance score
        all_contrasts.sort(key=lambda x: x["relevance_score"], reverse=True)
        
        # Take top contrasts
        top_contrasts = all_contrasts[:max_contrasts]
        
        return json.dumps({
            "ranked_contrasts": top_contrasts,
            "research_question": research_question,
            "total_contrasts_analyzed": len(deg_data),
            "contrasts_meeting_criteria": len(all_contrasts)
        })
    
    except Exception as e:
        logger.error(f"Error in rank_contrasts_by_relevance: {e}")
        return json.dumps({"error": str(e)})

@mcp.tool()
@log_tool_execution
async def analyze_gene_patterns(genes: List[str],
                               contrast_keys: List[str]) -> str:
    """
    Analyze expression patterns across a set of genes.
    
    Args:
        genes: List of gene identifiers.
        contrast_keys: List of contrast keys in format "analysis_id:contrast_id".
    
    Returns:
        JSON string with pattern analysis including correlations and clustering.
    """
    try:
        deg_data = _load_all_deg_data()
        
        # Build expression matrix
        expr_matrix = []
        gene_order = []
        contrast_order = []
        
        for gene in genes:
            gene_values = []
            for contrast_key in contrast_keys:
                if contrast_key not in deg_data:
                    gene_values.append(0.0)
                    continue
                
                df = deg_data[contrast_key]
                if 'Gene' in df.columns and 'logFC' in df.columns:
                    gene_row = df[df['Gene'] == gene]
                    if not gene_row.empty:
                        # Check if significant
                        p_col = 'adj.P.Val' if 'adj.P.Val' in df.columns else 'P.Value'
                        if p_col in df.columns:
                            p_val = gene_row.iloc[0][p_col]
                            lfc = gene_row.iloc[0]['logFC']
                            # Only include if significant
                            if p_val < 0.05 and abs(lfc) > 1.0:
                                gene_values.append(float(lfc))
                            else:
                                gene_values.append(0.0)
                        else:
                            gene_values.append(float(gene_row.iloc[0]['logFC']))
                    else:
                        gene_values.append(0.0)
                else:
                    gene_values.append(0.0)
            
            # Only add if gene has some non-zero values
            if any(v != 0 for v in gene_values):
                expr_matrix.append(gene_values)
                gene_order.append(gene)
        
        if not expr_matrix:
            return json.dumps({
                "error": "No expression data found for the specified genes and contrasts"
            })
        
        expr_array = np.array(expr_matrix)
        
        # Calculate correlation matrix between genes
        if len(gene_order) > 1:
            corr_matrix = np.corrcoef(expr_array)
        else:
            corr_matrix = np.array([[1.0]])
        
        # Perform hierarchical clustering if we have enough genes
        clusters = {}
        if len(gene_order) > 3:
            # Calculate linkage
            gene_linkage = linkage(pdist(expr_array), method='complete')
            # Get clusters (max 5 clusters)
            n_clusters = min(5, len(gene_order) // 2)
            cluster_labels = fcluster(gene_linkage, n_clusters, criterion='maxclust')
            
            for i, gene in enumerate(gene_order):
                cluster_id = int(cluster_labels[i])
                if cluster_id not in clusters:
                    clusters[cluster_id] = []
                clusters[cluster_id].append(gene)
        
        # Find genes with similar patterns
        patterns = []
        for i, gene in enumerate(gene_order):
            # Find highly correlated genes
            if len(gene_order) > 1:
                correlations = corr_matrix[i]
                similar_genes = []
                opposite_genes = []
                
                for j, other_gene in enumerate(gene_order):
                    if i != j:
                        if correlations[j] > 0.7:
                            similar_genes.append({
                                "gene": other_gene,
                                "correlation": round(float(correlations[j]), 3)
                            })
                        elif correlations[j] < -0.7:
                            opposite_genes.append({
                                "gene": other_gene,
                                "correlation": round(float(correlations[j]), 3)
                            })
            else:
                similar_genes = []
                opposite_genes = []
            
            # Calculate pattern statistics
            gene_expr = expr_array[i]
            non_zero_expr = gene_expr[gene_expr != 0]
            
            pattern_info = {
                "gene": gene,
                "expression_values": [round(float(v), 3) for v in gene_expr],
                "similar_genes": similar_genes,
                "opposite_genes": opposite_genes,
                "statistics": {
                    "mean_lfc": round(float(np.mean(non_zero_expr)), 3) if len(non_zero_expr) > 0 else 0,
                    "std_lfc": round(float(np.std(non_zero_expr)), 3) if len(non_zero_expr) > 0 else 0,
                    "percent_significant": round(len(non_zero_expr) / len(gene_expr) * 100, 1),
                    "consistent_direction": bool(np.all(non_zero_expr > 0) or np.all(non_zero_expr < 0))
                }
            }
            
            if clusters and i < len(cluster_labels):
                pattern_info["cluster"] = int(cluster_labels[i])
            
            patterns.append(pattern_info)
        
        # Perform PCA if we have enough dimensions
        pca_results = None
        if len(gene_order) > 2 and len(contrast_keys) > 2:
            try:
                scaler = StandardScaler()
                scaled_data = scaler.fit_transform(expr_array)
                pca = PCA(n_components=min(2, len(gene_order), len(contrast_keys)))
                pca_coords = pca.fit_transform(scaled_data)
                
                pca_results = {
                    "explained_variance": [round(float(v), 3) for v in pca.explained_variance_ratio_],
                    "coordinates": [
                        {
                            "gene": gene_order[i],
                            "pc1": round(float(pca_coords[i, 0]), 3),
                            "pc2": round(float(pca_coords[i, 1]), 3) if pca_coords.shape[1] > 1 else 0
                        }
                        for i in range(len(gene_order))
                    ]
                }
            except Exception as e:
                logger.warning(f"PCA failed: {e}")
        
        result = {
            "patterns": patterns,
            "contrasts_analyzed": len(contrast_keys),
            "genes_analyzed": len(gene_order),
            "clustering": clusters if clusters else None,
            "pca": pca_results
        }
        
        return json.dumps(result)
    
    except Exception as e:
        logger.error(f"Error in analyze_gene_patterns: {e}")
        return json.dumps({"error": str(e)})

@mcp.tool()
@log_tool_execution
async def find_expression_clusters(min_genes_per_cluster: int = 5,
                                  correlation_threshold: float = 0.7,
                                  p_value_threshold: float = 0.05,
                                  lfc_threshold: float = 1.0) -> str:
    """
    Find clusters of co-expressed genes across all contrasts.
    
    Args:
        min_genes_per_cluster: Minimum number of genes required to form a cluster.
        correlation_threshold: Minimum correlation required for genes to cluster together.
        p_value_threshold: P-value threshold for significance.
        lfc_threshold: Log fold change threshold.
    
    Returns:
        JSON string with identified gene clusters and their properties.
    """
    try:
        deg_data = _load_all_deg_data()
        
        # Build a matrix of all genes across all contrasts
        all_genes = set()
        for df in deg_data.values():
            if 'Gene' in df.columns:
                all_genes.update(df['Gene'].tolist())
        
        # Limit to reasonable number of genes
        if len(all_genes) > 1000:
            # Select top variable genes
            gene_variance = {}
            for gene in all_genes:
                values = []
                for df in deg_data.values():
                    if 'Gene' in df.columns and 'logFC' in df.columns:
                        gene_row = df[df['Gene'] == gene]
                        if not gene_row.empty:
                            p_col = 'adj.P.Val' if 'adj.P.Val' in df.columns else 'P.Value'
                            if p_col in df.columns:
                                if (gene_row.iloc[0][p_col] < p_value_threshold and 
                                    abs(gene_row.iloc[0]['logFC']) > lfc_threshold):
                                    values.append(gene_row.iloc[0]['logFC'])
                
                if len(values) > 1:
                    gene_variance[gene] = np.var(values)
            
            # Select top 500 most variable genes
            sorted_genes = sorted(gene_variance.items(), key=lambda x: x[1], reverse=True)
            all_genes = set([g[0] for g in sorted_genes[:500]])
        
        # Build expression matrix
        gene_list = list(all_genes)
        contrast_list = list(deg_data.keys())
        expr_matrix = np.zeros((len(gene_list), len(contrast_list)))
        
        for i, gene in enumerate(gene_list):
            for j, contrast_key in enumerate(contrast_list):
                df = deg_data[contrast_key]
                if 'Gene' in df.columns and 'logFC' in df.columns:
                    gene_row = df[df['Gene'] == gene]
                    if not gene_row.empty:
                        p_col = 'adj.P.Val' if 'adj.P.Val' in df.columns else 'P.Value'
                        if p_col in df.columns:
                            if (gene_row.iloc[0][p_col] < p_value_threshold and 
                                abs(gene_row.iloc[0]['logFC']) > lfc_threshold):
                                expr_matrix[i, j] = gene_row.iloc[0]['logFC']
        
        # Remove genes with too few significant values
        min_significant = 2
        significant_counts = np.sum(expr_matrix != 0, axis=1)
        valid_genes_mask = significant_counts >= min_significant
        expr_matrix = expr_matrix[valid_genes_mask]
        gene_list = [g for i, g in enumerate(gene_list) if valid_genes_mask[i]]
        
        if len(gene_list) < min_genes_per_cluster:
            return json.dumps({
                "clusters": [],
                "message": f"Not enough genes with significant expression ({len(gene_list)}) to form clusters"
            })
        
        # Calculate correlation matrix
        corr_matrix = np.corrcoef(expr_matrix)
        
        # Perform hierarchical clustering
        linkage_matrix = linkage(pdist(expr_matrix), method='average')
        
        # Find optimal number of clusters (between 2 and 10)
        from scipy.cluster.hierarchy import fcluster
        best_clusters = None
        best_score = -1
        
        for n_clusters in range(2, min(11, len(gene_list) // min_genes_per_cluster)):
            cluster_labels = fcluster(linkage_matrix, n_clusters, criterion='maxclust')
            
            # Calculate average within-cluster correlation
            cluster_scores = []
            for cluster_id in range(1, n_clusters + 1):
                cluster_genes_idx = np.where(cluster_labels == cluster_id)[0]
                if len(cluster_genes_idx) >= min_genes_per_cluster:
                    # Calculate average correlation within cluster
                    cluster_corr = corr_matrix[np.ix_(cluster_genes_idx, cluster_genes_idx)]
                    avg_corr = np.mean(cluster_corr[np.triu_indices_from(cluster_corr, k=1)])
                    if avg_corr >= correlation_threshold:
                        cluster_scores.append(avg_corr)
            
            if cluster_scores and np.mean(cluster_scores) > best_score:
                best_score = np.mean(cluster_scores)
                best_clusters = cluster_labels
        
        if best_clusters is None:
            return json.dumps({
                "clusters": [],
                "message": "No clusters found meeting the correlation threshold"
            })
        
        # Extract cluster information
        clusters = []
        for cluster_id in np.unique(best_clusters):
            cluster_genes_idx = np.where(best_clusters == cluster_id)[0]
            if len(cluster_genes_idx) >= min_genes_per_cluster:
                cluster_genes = [gene_list[i] for i in cluster_genes_idx]
                
                # Calculate cluster statistics
                cluster_expr = expr_matrix[cluster_genes_idx]
                mean_expr = np.mean(cluster_expr, axis=0)
                
                # Find contrasts where cluster is most active
                active_contrasts = []
                for j, contrast_key in enumerate(contrast_list):
                    if abs(mean_expr[j]) > lfc_threshold:
                        active_contrasts.append({
                            "contrast": contrast_key,
                            "mean_lfc": round(float(mean_expr[j]), 3)
                        })
                
                # Sort by absolute mean LFC
                active_contrasts.sort(key=lambda x: abs(x["mean_lfc"]), reverse=True)
                
                clusters.append({
                    "cluster_id": int(cluster_id),
                    "genes": cluster_genes,
                    "size": len(cluster_genes),
                    "mean_correlation": round(float(
                        np.mean(corr_matrix[np.ix_(cluster_genes_idx, cluster_genes_idx)])
                    ), 3),
                    "top_active_contrasts": active_contrasts[:5]
                })
        
        return json.dumps({
            "clusters": clusters,
            "total_genes_analyzed": len(gene_list),
            "total_contrasts": len(contrast_list)
        })
    
    except Exception as e:
        logger.error(f"Error in find_expression_clusters: {e}")
        return json.dumps({"error": str(e)})

@mcp.tool()
@log_tool_execution
async def calculate_gene_statistics(genes: List[str],
                                   min_contrasts: int = 1) -> str:
    """
    Calculate summary statistics for specified genes across all contrasts.
    
    Args:
        genes: List of gene identifiers.
        min_contrasts: Minimum number of contrasts a gene must appear in.
    
    Returns:
        JSON string with comprehensive statistics for each gene.
    """
    try:
        deg_data = _load_all_deg_data()
        gene_stats = []
        
        for gene in genes:
            # Collect all data for this gene
            lfc_values = []
            p_values = []
            significant_contrasts = []
            all_contrasts = []
            
            for contrast_key, df in deg_data.items():
                if 'Gene' not in df.columns:
                    continue
                
                gene_row = df[df['Gene'] == gene]
                if not gene_row.empty:
                    row = gene_row.iloc[0]
                    
                    # Get values
                    if 'logFC' in df.columns:
                        lfc = float(row['logFC'])
                        lfc_values.append(lfc)
                        
                        # Get p-value
                        p_val = None
                        if 'adj.P.Val' in df.columns:
                            p_val = float(row['adj.P.Val'])
                        elif 'P.Value' in df.columns:
                            p_val = float(row['P.Value'])
                        
                        if p_val is not None:
                            p_values.append(p_val)
                            
                            # Check if significant
                            if p_val < 0.05 and abs(lfc) > 1.0:
                                significant_contrasts.append({
                                    "contrast": contrast_key,
                                    "lfc": round(lfc, 3),
                                    "p_value": p_val
                                })
                        
                        all_contrasts.append(contrast_key)
            
            if len(all_contrasts) >= min_contrasts:
                # Calculate statistics
                stats = {
                    "gene": gene,
                    "total_contrasts": len(all_contrasts),
                    "significant_contrasts": len(significant_contrasts),
                    "percentage_significant": round(
                        len(significant_contrasts) / len(all_contrasts) * 100, 1
                    ) if all_contrasts else 0
                }
                
                if lfc_values:
                    stats["lfc_statistics"] = {
                        "mean": round(float(np.mean(lfc_values)), 3),
                        "median": round(float(np.median(lfc_values)), 3),
                        "std": round(float(np.std(lfc_values)), 3),
                        "min": round(float(np.min(lfc_values)), 3),
                        "max": round(float(np.max(lfc_values)), 3),
                        "consistent_direction": bool(
                            all(v > 0 for v in lfc_values) or all(v < 0 for v in lfc_values)
                        )
                    }
                
                if significant_contrasts:
                    sig_lfcs = [c["lfc"] for c in significant_contrasts]
                    stats["significant_lfc_statistics"] = {
                        "mean": round(float(np.mean(sig_lfcs)), 3),
                        "median": round(float(np.median(sig_lfcs)), 3),
                        "consistent_direction": bool(
                            all(v > 0 for v in sig_lfcs) or all(v < 0 for v in sig_lfcs)
                        )
                    }
                    
                    # Top significant contrasts
                    stats["top_significant_contrasts"] = sorted(
                        significant_contrasts,
                        key=lambda x: abs(x["lfc"]),
                        reverse=True
                    )[:5]
                
                if p_values:
                    stats["p_value_statistics"] = {
                        "min": float(np.min(p_values)),
                        "median": float(np.median(p_values)),
                        "geometric_mean": float(np.exp(np.mean(np.log(p_values))))
                    }
                
                gene_stats.append(stats)
        
        return json.dumps({
            "gene_statistics": gene_stats,
            "genes_analyzed": len(gene_stats),
            "total_contrasts_in_dataset": len(deg_data)
        })
    
    except Exception as e:
        logger.error(f"Error in calculate_gene_statistics: {e}")
        return json.dumps({"error": str(e)})

@mcp.tool()
@log_tool_execution
async def identify_consistent_degs(min_frequency: float = 0.8,
                                  consistency_threshold: float = 0.9,
                                  p_value_threshold: float = 0.05,
                                  lfc_threshold: float = 1.0,
                                  max_genes: int = 100) -> str:
    """
    Find genes with consistent differential expression across contrasts.
    
    Args:
        min_frequency: Minimum fraction of contrasts where gene must be significant.
        consistency_threshold: Fraction of significant contrasts that must have same direction.
        p_value_threshold: P-value threshold for significance.
        lfc_threshold: Log fold change threshold.
        max_genes: Maximum number of genes to return.
    
    Returns:
        JSON string with consistently differentially expressed genes.
    """
    try:
        deg_data = _load_all_deg_data()
        
        # Track gene behavior across contrasts
        gene_summary = defaultdict(lambda: {
            "significant_count": 0,
            "total_count": 0,
            "up_count": 0,
            "down_count": 0,
            "lfc_values": [],
            "contrasts": []
        })
        
        # Analyze each contrast
        for contrast_key, df in deg_data.items():
            if 'Gene' not in df.columns or 'logFC' not in df.columns:
                continue
            
            # Find p-value column
            p_col = None
            if 'adj.P.Val' in df.columns:
                p_col = 'adj.P.Val'
            elif 'P.Value' in df.columns:
                p_col = 'P.Value'
            
            if not p_col:
                continue
            
            # Process each gene
            for _, row in df.iterrows():
                gene = row['Gene']
                lfc = row['logFC']
                p_val = row[p_col]
                
                gene_summary[gene]["total_count"] += 1
                
                if p_val < p_value_threshold and abs(lfc) > lfc_threshold:
                    gene_summary[gene]["significant_count"] += 1
                    gene_summary[gene]["lfc_values"].append(float(lfc))
                    gene_summary[gene]["contrasts"].append(contrast_key)
                    
                    if lfc > 0:
                        gene_summary[gene]["up_count"] += 1
                    else:
                        gene_summary[gene]["down_count"] += 1
        
        # Filter for consistent genes
        consistent_genes = []
        
        for gene, data in gene_summary.items():
            if data["total_count"] == 0:
                continue
            
            # Check frequency
            frequency = data["significant_count"] / data["total_count"]
            if frequency < min_frequency:
                continue
            
            # Check consistency
            if data["significant_count"] > 0:
                direction_consistency = max(
                    data["up_count"] / data["significant_count"],
                    data["down_count"] / data["significant_count"]
                )
                
                if direction_consistency >= consistency_threshold:
                    # Calculate statistics
                    lfc_array = np.array(data["lfc_values"])
                    
                    consistent_genes.append({
                        "gene": gene,
                        "frequency": round(frequency, 3),
                        "consistency": round(direction_consistency, 3),
                        "direction": "up" if data["up_count"] > data["down_

def main():
    """Main entry point for the MCP analysis server."""
    parser = argparse.ArgumentParser(description="MCP Analysis Server")
    parser.add_argument("command", choices=["server"], help="Run as MCP server")
    parser.add_argument("--results-dir", help="Path to results directory")
    args = parser.parse_args()

    global RESULTS_DIR
    if args.results_dir:
        RESULTS_DIR = args.results_dir
        log(f"Using results directory from command line: {RESULTS_DIR}")

    logger.info(f"Starting analysis MCP server with results_dir: {RESULTS_DIR}")
    if not RESULTS_DIR:
        logger.warning("RESULTS_DIR environment variable not set")
    elif not os.path.exists(RESULTS_DIR):
        logger.warning(f"RESULTS_DIR does not exist: {RESULTS_DIR}")

    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()

