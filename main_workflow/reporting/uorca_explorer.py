#!/usr/bin/env python3
"""
UORCA Explorer - A Streamlit app for interactive exploration of UORCA RNA-seq results

This app allows you to:
- Select any subset of genes, datasets, or contrasts
- Instantly see updated heatmaps and expression plots
- Explore your RNA-seq results interactively

Usage:
  streamlit run uorca_explorer.py

Or with custom port:
  streamlit run uorca_explorer.py --server.port 8501
"""

import os
import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
import sys
import json
from datetime import datetime
from dataclasses import dataclass
from typing import List, Dict, Any, Optional, Tuple, Set
import statistics
import numpy as np
try:
    from scipy import stats
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    # Implement fallback normal distribution functions
    class NormalDistribution:
        @staticmethod
        def isf(p):
            """Inverse survival function (inverse of 1-cdf) - approximation."""
            if p <= 0:
                return float('inf')
            if p >= 1:
                return float('-inf')

            # Use Box-Muller approximation for inverse normal
            # This is a simple approximation - not as accurate as scipy
            if p > 0.5:
                return -NormalDistribution.isf(1 - p)

            # Rational approximation for inverse normal CDF
            # Based on Beasley-Springer-Moro algorithm
            a = [0, -3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
            b = [0, -5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01]

            if p < 1e-8:
                return -8.0  # Approximate for very small p

            y = math.sqrt(-2.0 * math.log(p))

            num = a[6]
            for i in range(5, -1, -1):
                num = num * y + a[i]

            den = b[5]
            for i in range(4, -1, -1):
                den = den * y + b[i]

            return y + num / den

        @staticmethod
        def sf(x):
            """Survival function (1-cdf) - approximation."""
            # Complementary error function approximation
            # For normal distribution: sf(x) = 0.5 * erfc(x/sqrt(2))
            if x > 6:
                return 0.0
            if x < -6:
                return 1.0

            # Abramowitz and Stegun approximation for erfc
            z = abs(x) / math.sqrt(2)
            t = 1.0 / (1.0 + 0.5 * z)

            erfc = t * math.exp(-z*z - 1.26551223 + t * (1.00002368 + t * (0.37409196 +
                   t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 +
                   t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 +
                   t * 0.17087277)))))))))

            if x >= 0:
                return 0.5 * erfc
            else:
                return 1.0 - 0.5 * erfc

    # Create a stats-like object for compatibility
    class StatsModule:
        norm = NormalDistribution()

    stats = StatsModule()

import re
from collections import Counter
import math

# Load environment variables from .env file
try:
    from dotenv import load_dotenv
    # Try loading from current directory first, then parent directories
    load_dotenv()
    # Also try loading from project root
    project_root = Path(__file__).parent.parent.parent
    env_file = project_root / ".env"
    if env_file.exists():
        load_dotenv(env_file)
except ImportError:
    # dotenv not available, continue without it
    pass

# Use current directory to import ResultsIntegration
script_dir = os.path.dirname(os.path.abspath(__file__))
# ResultsIntegration.py is in the same directory
from ResultsIntegration import ResultsIntegrator

# AI Landing Page functionality - integrated directly
try:
    import openai
    from openai import OpenAI
    LANDING_PAGE_AVAILABLE = True
except ImportError as e:
    LANDING_PAGE_AVAILABLE = False
    # Only show warning in sidebar, not main area

# Helper function for contrast labels
def short_label(full_label: str) -> str:
    """Create short labels for contrast multiselect display"""
    # keeps 'KO_vs_WT' out of 'GSE12345:KO_vs_WT – long sentence …'
    return full_label.split(":", 1)[-1].split(" –")[0].split(" - ")[0][:25]

# ========================================================================
# AI LANDING PAGE FUNCTIONALITY - INTEGRATED
# ========================================================================

@dataclass
class ContrastSelection:
    """Container for selected contrast with justification"""
    analysis_id: str
    contrast_id: str
    relevance_score: float
    justification: str
    deg_count: int

@dataclass
class ThresholdSelection:
    """Container for automatically selected statistical thresholds"""
    fdr_cutoff: float
    logfc_cutoff: float
    min_frequency: int
    justification: str

@dataclass
class LandingPageData:
    """Container for all landing page data"""
    selected_contrasts: List[ContrastSelection]
    thresholds: ThresholdSelection
    top_genes: List[str]
    heatmap_fig: Optional[go.Figure]
    gene_table: pd.DataFrame
    narrative: str

def generate_ai_landing_page(integrator, biological_prompt: str, max_contrasts: int = 8, max_genes: int = 50) -> Optional[LandingPageData]:
    """
    Generate AI-assisted landing page data using the provided integrator.
    """
    try:
        # Step 1: Select relevant contrasts
        selected_contrasts = select_relevant_contrasts(integrator, biological_prompt, max_contrasts)
        if not selected_contrasts:
            return None

        # Step 2: Select optimal thresholds
        thresholds = select_optimal_thresholds(selected_contrasts, biological_prompt)

        # Step 3: Aggregate and rank genes
        top_genes = aggregate_and_rank_genes(integrator, selected_contrasts, thresholds, max_genes)

        # Step 4: Create visualizations
        heatmap_fig = create_landing_heatmap(integrator, top_genes, selected_contrasts, thresholds)
        gene_table = create_gene_table(integrator, top_genes, selected_contrasts, thresholds)

        # Step 5: Generate narrative
        narrative = generate_narrative(selected_contrasts, thresholds, top_genes, biological_prompt)

        return LandingPageData(
            selected_contrasts=selected_contrasts,
            thresholds=thresholds,
            top_genes=top_genes,
            heatmap_fig=heatmap_fig,
            gene_table=gene_table,
            narrative=narrative
        )

    except Exception as e:
        st.error(f"Error generating AI landing page: {str(e)}")
        return None

class SimpleTfidfVectorizer:
    """Simple TF-IDF vectorizer as fallback for sklearn."""

    def __init__(self, max_features=384, stop_words='english'):
        self.max_features = max_features
        self.vocabulary_ = {}
        self.idf_ = {}
        self.stop_words = self._get_stop_words() if stop_words == 'english' else set()

    def _get_stop_words(self):
        # Basic English stop words
        return {'the', 'a', 'an', 'and', 'or', 'but', 'in', 'on', 'at', 'to', 'for', 'of', 'with', 'by', 'is', 'are', 'was', 'were', 'be', 'been', 'have', 'has', 'had', 'do', 'does', 'did', 'will', 'would', 'could', 'should', 'this', 'that', 'these', 'those', 'i', 'you', 'he', 'she', 'it', 'we', 'they', 'me', 'him', 'her', 'us', 'them'}

    def _tokenize(self, text):
        # Simple tokenization
        tokens = re.findall(r'\b\w+\b', text.lower())
        return [t for t in tokens if t not in self.stop_words and len(t) > 1]

    def fit(self, documents):
        # Count document frequencies
        df = Counter()
        all_tokens = set()

        for doc in documents:
            tokens = set(self._tokenize(doc))
            all_tokens.update(tokens)
            for token in tokens:
                df[token] += 1

        # Select top features
        n_docs = len(documents)
        if len(all_tokens) > self.max_features:
            # Sort by document frequency and take top features
            sorted_tokens = sorted(df.items(), key=lambda x: x[1], reverse=True)
            selected_tokens = [token for token, freq in sorted_tokens[:self.max_features]]
        else:
            selected_tokens = list(all_tokens)

        # Build vocabulary
        self.vocabulary_ = {token: i for i, token in enumerate(selected_tokens)}

        # Calculate IDF
        for token in selected_tokens:
            self.idf_[token] = math.log(n_docs / (df[token] + 1))

        return self

    def transform(self, documents):
        vectors = []
        for doc in documents:
            tokens = self._tokenize(doc)
            tf = Counter(tokens)

            # Create TF-IDF vector
            vector = np.zeros(len(self.vocabulary_))
            for token, count in tf.items():
                if token in self.vocabulary_:
                    idx = self.vocabulary_[token]
                    tf_score = count / len(tokens) if len(tokens) > 0 else 0
                    vector[idx] = tf_score * self.idf_.get(token, 0)

            # Normalize
            norm = np.linalg.norm(vector)
            if norm > 0:
                vector = vector / norm

            vectors.append(vector)

        return np.array(vectors)

def get_text_embedding(text: str, use_api: bool = True) -> np.ndarray:
    """Get embedding for text using OpenAI API or fallback to TF-IDF."""
    if use_api and LANDING_PAGE_AVAILABLE and os.getenv("OPENAI_API_KEY"):
        try:
            client = OpenAI()
            response = client.embeddings.create(
                model="text-embedding-3-small",
                input=text
            )
            return np.array(response.data[0].embedding)
        except Exception:
            # Fall back to TF-IDF if API fails
            pass

    # Fallback: Use TF-IDF with cached vectorizer
    if not hasattr(get_text_embedding, 'vectorizer'):
        get_text_embedding.vectorizer = SimpleTfidfVectorizer(max_features=384, stop_words='english')
        get_text_embedding.fitted = False

    # For fallback, we need to maintain a corpus to fit the vectorizer
    if not hasattr(get_text_embedding, 'corpus'):
        get_text_embedding.corpus = []

    if text not in get_text_embedding.corpus:
        get_text_embedding.corpus.append(text)
        get_text_embedding.fitted = False

    if not get_text_embedding.fitted:
        try:
            get_text_embedding.vectorizer.fit(get_text_embedding.corpus)
            get_text_embedding.fitted = True
        except:
            # Return zero vector if fitting fails
            return np.zeros(384)

    try:
        vector = get_text_embedding.vectorizer.transform([text])
        return vector[0]
    except:
        return np.zeros(384)

def shortlist_contrasts(query_text: str, contrasts: List[Dict], K: int = 30) -> List[Tuple[str, str, float]]:
    """Stage A: Semantic eligibility filter using embeddings."""
    if not contrasts:
        return []

    # Get query embedding
    query_embedding = get_text_embedding(query_text)

    # Get embeddings for all contrasts and compute similarities
    scored_contrasts = []
    for contrast in contrasts:
        # Combine title and description for embedding
        text = f"{contrast['contrast_id']} {contrast['description']}"
        contrast_embedding = get_text_embedding(text)

        # Compute cosine similarity
        if np.linalg.norm(query_embedding) > 0 and np.linalg.norm(contrast_embedding) > 0:
            similarity = np.dot(query_embedding, contrast_embedding) / (
                np.linalg.norm(query_embedding) * np.linalg.norm(contrast_embedding)
            )
        else:
            similarity = 0.0

        scored_contrasts.append((contrast['analysis_id'], contrast['contrast_id'], similarity, contrast))

    # Sort by similarity and return top K
    scored_contrasts.sort(key=lambda x: x[2], reverse=True)
    return scored_contrasts[:K]

def calculate_contrast_metrics(integrator, analysis_id: str, contrast_id: str) -> Dict:
    """Calculate quality metrics for a contrast."""
    metrics = {
        'sample_size': 0,
        'deg_count': 0,
        'deg_richness': 0.0,
        'data_quality': 0.0
    }

    # Get sample size from analysis info
    if analysis_id in integrator.analysis_info:
        metrics['sample_size'] = integrator.analysis_info[analysis_id].get('number_of_samples', 0)

    # Get DEG statistics
    if analysis_id in integrator.deg_data and contrast_id in integrator.deg_data[analysis_id]:
        df = integrator.deg_data[analysis_id][contrast_id]

        if 'adj.P.Val' in df.columns and 'logFC' in df.columns:
            # Count significant DEGs
            sig_genes = (df['adj.P.Val'] < 0.05) & (abs(df['logFC']) > 1.0)
            metrics['deg_count'] = sig_genes.sum()

            # Calculate DEG richness (proportion of tested genes that are significant)
            total_genes = len(df)
            metrics['deg_richness'] = metrics['deg_count'] / max(total_genes, 1)

            # Data quality score based on range of fold changes and p-value distribution
            if metrics['deg_count'] > 0:
                sig_df = df[sig_genes]
                if len(sig_df) > 0:
                    lfc_range = sig_df['logFC'].max() - sig_df['logFC'].min()
                    metrics['data_quality'] = min(lfc_range / 10.0, 1.0)  # Normalize to 0-1

    return metrics

def rank_contrasts(shortlisted: List[Tuple], integrator, max_final: int = 10,
                  weights: Tuple[float, float, float, float] = (0.4, 0.2, 0.2, 0.2)) -> List[ContrastSelection]:
    """Stage B: Evidence-aware ranking with diversity."""
    if not shortlisted:
        return []

    # Prepare contrast data with metrics
    candidates = []
    for analysis_id, contrast_id, similarity, original_data in shortlisted:
        metrics = calculate_contrast_metrics(integrator, analysis_id, contrast_id)

        candidates.append({
            'analysis_id': analysis_id,
            'contrast_id': contrast_id,
            'similarity': similarity,
            'sample_size': metrics['sample_size'],
            'deg_count': metrics['deg_count'],
            'deg_richness': metrics['deg_richness'],
            'data_quality': metrics['data_quality'],
            'original_data': original_data
        })

    # Normalize metrics to 0-1 range
    if candidates:
        max_similarity = max(c['similarity'] for c in candidates)
        max_sample_size = max(c['sample_size'] for c in candidates)
        max_deg_count = max(c['deg_count'] for c in candidates)
        max_deg_richness = max(c['deg_richness'] for c in candidates)

        for c in candidates:
            c['norm_similarity'] = c['similarity'] / max(max_similarity, 1e-6)
            c['norm_sample_size'] = c['sample_size'] / max(max_sample_size, 1)
            c['norm_deg_count'] = c['deg_count'] / max(max_deg_count, 1)
            c['norm_deg_richness'] = c['deg_richness'] / max(max_deg_richness, 1e-6)

    # Greedy selection with diversity
    selected = []
    selected_gene_sets = []

    while candidates and len(selected) < max_final:
        # Calculate diversity scores
        for candidate in candidates:
            # Get gene set for this contrast
            gene_set = set()
            if (candidate['analysis_id'] in integrator.deg_data and
                candidate['contrast_id'] in integrator.deg_data[candidate['analysis_id']]):
                df = integrator.deg_data[candidate['analysis_id']][candidate['contrast_id']]
                if 'adj.P.Val' in df.columns and 'logFC' in df.columns:
                    sig_genes = df[(df['adj.P.Val'] < 0.05) & (abs(df['logFC']) > 1.0)]
                    gene_set = set(sig_genes['Gene'].tolist())

            # Calculate diversity (1 - max jaccard with selected sets)
            diversity = 1.0
            if selected_gene_sets:
                max_jaccard = 0.0
                for selected_set in selected_gene_sets:
                    if len(gene_set) > 0 and len(selected_set) > 0:
                        intersection = len(gene_set & selected_set)
                        union = len(gene_set | selected_set)
                        jaccard = intersection / union if union > 0 else 0.0
                        max_jaccard = max(max_jaccard, jaccard)
                diversity = 1.0 - max_jaccard

            # Calculate composite score
            candidate['diversity'] = diversity
            candidate['composite_score'] = (
                weights[0] * candidate['norm_similarity'] +
                weights[1] * candidate['norm_sample_size'] +
                weights[2] * candidate['norm_deg_count'] +
                weights[3] * candidate['diversity']
            )

        # Select best candidate
        best_candidate = max(candidates, key=lambda x: x['composite_score'])

        # Add to selected
        selected.append(ContrastSelection(
            analysis_id=best_candidate['analysis_id'],
            contrast_id=best_candidate['contrast_id'],
            relevance_score=best_candidate['composite_score'] * 10,  # Scale to 0-10
            justification=f"Selected based on semantic similarity ({best_candidate['similarity']:.2f}), "
                         f"sample size ({best_candidate['sample_size']}), "
                         f"and {best_candidate['deg_count']} DEGs with diversity score {best_candidate['diversity']:.2f}",
            deg_count=best_candidate['deg_count']
        ))

        # Add gene set to selected sets
        gene_set = set()
        if (best_candidate['analysis_id'] in integrator.deg_data and
            best_candidate['contrast_id'] in integrator.deg_data[best_candidate['analysis_id']]):
            df = integrator.deg_data[best_candidate['analysis_id']][best_candidate['contrast_id']]
            if 'adj.P.Val' in df.columns and 'logFC' in df.columns:
                sig_genes = df[(df['adj.P.Val'] < 0.05) & (abs(df['logFC']) > 1.0)]
                gene_set = set(sig_genes['Gene'].tolist())
        selected_gene_sets.append(gene_set)

        # Remove from candidates
        candidates.remove(best_candidate)

    return selected

def select_relevant_contrasts(integrator, biological_prompt: str, max_contrasts: int) -> List[ContrastSelection]:
    """Select the most relevant contrasts using improved two-stage method."""

    # Get all available contrasts
    all_contrasts = []
    for analysis_id, contrasts in integrator.deg_data.items():
        for contrast_id in contrasts.keys():
            description = integrator._get_contrast_description(analysis_id, contrast_id)

            # Count potential DEGs
            df = contrasts[contrast_id]
            deg_count = 0
            if 'adj.P.Val' in df.columns and 'logFC' in df.columns:
                deg_count = ((df['adj.P.Val'] < 0.1) & (abs(df['logFC']) > 0.5)).sum()

            all_contrasts.append({
                'analysis_id': analysis_id,
                'contrast_id': contrast_id,
                'description': description,
                'deg_count': deg_count
            })

    if not all_contrasts:
        return []

    try:
        # Stage A: Semantic eligibility filter
        shortlisted = shortlist_contrasts(biological_prompt, all_contrasts, K=min(30, len(all_contrasts)))

        # Stage B: Evidence-aware ranking
        selected = rank_contrasts(shortlisted, integrator, max_final=max_contrasts)

        return selected

    except Exception as e:
        st.warning(f"Advanced selection failed, using fallback: {str(e)}")
        # Fallback to original heuristic method
        return score_contrasts_heuristic(all_contrasts, max_contrasts)

def score_contrasts_with_ai(contrasts: List[Dict], biological_prompt: str, max_contrasts: int) -> List[ContrastSelection]:
    """Score contrasts using OpenAI API."""

    client = OpenAI()

    # Prepare contrast descriptions
    contrast_descriptions = []
    for i, contrast in enumerate(contrasts):
        contrast_descriptions.append(
            f"{i+1}. {contrast['analysis_id']}_{contrast['contrast_id']}: "
            f"{contrast['description']} (Potential DEGs: {contrast['deg_count']})"
        )

    prompt = f"""
You are an expert RNA-seq analyst. Score the biological relevance of these differential expression contrasts for the research context: "{biological_prompt}"

Available contrasts:
{chr(10).join(contrast_descriptions)}

For each contrast, provide:
1. Relevance score (0-10, where 10 is most relevant)
2. Brief justification (1-2 sentences)

Consider:
- Biological significance for the research question
- Number of potential DEGs (more is generally better)
- Clarity of experimental design
- Scientific interest and interpretability

Respond in JSON format:
{{
  "scored_contrasts": [
    {{
      "contrast_number": 1,
      "relevance_score": 8.5,
      "justification": "This contrast directly addresses the research question by comparing..."
    }},
    ...
  ]
}}
"""

    try:
        response = client.chat.completions.create(
            model="o4-mini",
            messages=[{"role": "user", "content": prompt}],
            response_format={"type": "json_object"}
        )

        result = json.loads(response.choices[0].message.content)
        scored_contrasts = []

        for score_data in result.get("scored_contrasts", []):
            idx = score_data["contrast_number"] - 1
            if 0 <= idx < len(contrasts):
                contrast = contrasts[idx]
                scored_contrasts.append(ContrastSelection(
                    analysis_id=contrast['analysis_id'],
                    contrast_id=contrast['contrast_id'],
                    relevance_score=score_data["relevance_score"],
                    justification=score_data["justification"],
                    deg_count=contrast['deg_count']
                ))

        # Sort by relevance score and return top ones
        scored_contrasts.sort(key=lambda x: x.relevance_score, reverse=True)
        return scored_contrasts[:max_contrasts]

    except Exception as e:
        raise Exception(f"OpenAI API call failed: {str(e)}")

def score_contrasts_heuristic(contrasts: List[Dict], max_contrasts: int) -> List[ContrastSelection]:
    """Fallback scoring based on DEG counts and description informativeness."""

    scored_contrasts = []

    for contrast in contrasts:
        # Simple heuristic: score based on DEG count and description
        base_score = min(8, contrast['deg_count'] / 100 * 6)  # Up to 6 points for DEG count
        desc_score = min(2, len(contrast['description'].split()) / 20 * 2)  # Up to 2 points for description

        total_score = base_score + desc_score
        justification = f"Selected based on {contrast['deg_count']} potential DEGs and experimental design clarity."

        scored_contrasts.append(ContrastSelection(
            analysis_id=contrast['analysis_id'],
            contrast_id=contrast['contrast_id'],
            relevance_score=total_score,
            justification=justification,
            deg_count=contrast['deg_count']
        ))

    # Sort by score and return top ones
    scored_contrasts.sort(key=lambda x: x.relevance_score, reverse=True)
    return scored_contrasts[:max_contrasts]

def select_optimal_thresholds(selected_contrasts: List[ContrastSelection], biological_prompt: str) -> ThresholdSelection:
    """Select optimal statistical thresholds."""

    total_contrasts = len(selected_contrasts)
    median_degs = np.median([c.deg_count for c in selected_contrasts]) if selected_contrasts else 0

    # Try AI threshold selection first
    try:
        if LANDING_PAGE_AVAILABLE and os.getenv("OPENAI_API_KEY"):
            return select_thresholds_with_ai(total_contrasts, median_degs, biological_prompt)
        else:
            return select_thresholds_heuristic(total_contrasts, median_degs)
    except Exception:
        return select_thresholds_heuristic(total_contrasts, median_degs)

def select_thresholds_with_ai(n_contrasts: int, median_degs: float, biological_prompt: str) -> ThresholdSelection:
    """Use AI to select optimal thresholds for meta-analysis."""

    client = OpenAI()

    prompt = f"""
You are an expert statistician selecting optimal thresholds for RNA-seq meta-analysis across multiple contrasts.

Dataset characteristics:
- Number of contrasts: {n_contrasts}
- Median potential DEGs per contrast: {median_degs:.0f}
- Research context: {biological_prompt}

Select appropriate thresholds for meta-analysis considering:
1. Meta-FDR cutoff: Controls false discovery rate for meta-analysis results (0.01-0.1 range)
2. Minimum frequency: How many contrasts a gene must be significant in (1 to {max(1, n_contrasts//2)})
3. Effect size consideration: Implicit in meta-analysis ranking by median absolute logFC

Guidelines:
- With many contrasts ({n_contrasts}), can use stricter meta-FDR and higher min_frequency
- Meta-analysis provides more power, so can be more stringent than single-study analysis
- Min_frequency should balance sensitivity with reproducibility

Respond in JSON format:
{{
  "fdr_cutoff": 0.05,
  "logfc_cutoff": 0.0,
  "min_frequency": 2,
  "justification": "Selected meta-FDR=0.05 for balanced discovery with {n_contrasts} contrasts..."
}}
"""

    response = client.chat.completions.create(
        model="gpt-4o-mini",
        messages=[{"role": "user", "content": prompt}],
        response_format={"type": "json_object"},
        temperature=0
    )

    result = json.loads(response.choices[0].message.content)
    return ThresholdSelection(
        fdr_cutoff=result["fdr_cutoff"],
        logfc_cutoff=0.0,  # Not used in meta-analysis
        min_frequency=result["min_frequency"],
        justification=result["justification"]
    )

def select_thresholds_heuristic(n_contrasts: int, median_degs: float) -> ThresholdSelection:
    """Fallback threshold selection using heuristics for meta-analysis."""

    # For meta-analysis, focus on meta-FDR and minimum frequency
    if n_contrasts < 3:
        fdr_cutoff = 0.05
        min_frequency = 1  # Need at least 1 for very few contrasts
    elif n_contrasts < 6:
        fdr_cutoff = 0.05
        min_frequency = 2  # Require replication in at least 2 contrasts
    elif median_degs > 100:
        fdr_cutoff = 0.01  # Can be stricter with many DEGs
        min_frequency = max(2, n_contrasts // 4)  # Higher frequency requirement
    else:
        fdr_cutoff = 0.05
        min_frequency = max(2, n_contrasts // 3)

    # Cap minimum frequency to reasonable values
    min_frequency = min(min_frequency, max(2, n_contrasts // 2))

    justification = f"Meta-analysis heuristic: meta-FDR={fdr_cutoff}, min_frequency={min_frequency} based on {n_contrasts} contrasts and median {median_degs:.0f} DEGs per contrast. Genes must be significant in ≥{min_frequency} contrasts to be included."

    return ThresholdSelection(
        fdr_cutoff=fdr_cutoff,
        logfc_cutoff=0.0,  # Not used in meta-analysis approach
        min_frequency=min_frequency,
        justification=justification
    )

def benjamini_hochberg_correction(p_values: np.ndarray) -> np.ndarray:
    """Apply Benjamini-Hochberg FDR correction."""
    if len(p_values) == 0:
        return np.array([])

    # Sort p-values and keep track of original indices
    sorted_indices = np.argsort(p_values)
    sorted_p = p_values[sorted_indices]

    # Calculate adjusted p-values
    n = len(p_values)
    adjusted_p = np.zeros(n)

    for i in range(n-1, -1, -1):
        if i == n-1:
            adjusted_p[i] = sorted_p[i]
        else:
            adjusted_p[i] = min(sorted_p[i] * n / (i + 1), adjusted_p[i + 1])

    # Restore original order
    result = np.zeros(n)
    result[sorted_indices] = adjusted_p

    return result

def meta_analysis_stouffer(integrator, selected_contrasts: List[ContrastSelection],
                          gene_list: List[str]) -> Dict[str, Dict]:
    """Perform sign-blind meta-analysis using weighted Stouffer Z-method."""

    gene_meta_stats = {}

    for gene in gene_list:
        p_values = []
        weights = []
        log_fcs = []
        significant_counts = {'up': 0, 'down': 0, 'total': 0}

        for contrast in selected_contrasts:
            if (contrast.analysis_id in integrator.deg_data and
                contrast.contrast_id in integrator.deg_data[contrast.analysis_id]):

                df = integrator.deg_data[contrast.analysis_id][contrast.contrast_id]

                if 'Gene' in df.columns and 'adj.P.Val' in df.columns and 'logFC' in df.columns:
                    gene_row = df[df['Gene'] == gene]

                    if not gene_row.empty:
                        row = gene_row.iloc[0]
                        p_val = row['adj.P.Val']
                        logfc = row['logFC']

                        # Get sample size for weighting (fallback to 10 if not available)
                        sample_size = integrator.analysis_info.get(
                            contrast.analysis_id, {}
                        ).get('number_of_samples', 10)

                        p_values.append(p_val)
                        weights.append(np.sqrt(sample_size))
                        log_fcs.append(logfc)

                        # Count significant directions
                        if p_val < 0.05 and abs(logfc) > 1.0:
                            significant_counts['total'] += 1
                            if logfc > 0:
                                significant_counts['up'] += 1
                            else:
                                significant_counts['down'] += 1

        if len(p_values) >= 2:  # Need at least 2 contrasts for meta-analysis
            # Convert p-values to Z-scores (two-tailed)
            z_scores = [stats.norm.isf(p/2) if p > 0 else stats.norm.isf(1e-16/2) for p in p_values]

            # Weighted Stouffer combination
            weights = np.array(weights)
            z_scores = np.array(z_scores)

            # Combined Z-score
            combined_z = np.sum(weights * z_scores) / np.sqrt(np.sum(weights**2))

            # Two-tailed meta p-value
            meta_p = 2 * stats.norm.sf(abs(combined_z))

            # Calculate median absolute logFC for ranking
            median_abs_logfc = np.median([abs(lfc) for lfc in log_fcs])

            gene_meta_stats[gene] = {
                'meta_p': meta_p,
                'median_abs_logfc': median_abs_logfc,
                'n_contrasts': len(p_values),
                'up_count': significant_counts['up'],
                'down_count': significant_counts['down'],
                'total_sig': significant_counts['total'],
                'log_fcs': log_fcs
            }

    return gene_meta_stats

def aggregate_and_rank_genes(integrator, selected_contrasts: List[ContrastSelection],
                           thresholds: ThresholdSelection, max_genes: int) -> List[str]:
    """Aggregate DEGs using meta-analysis and rank by effect size and significance."""

    # Get all genes that appear in at least one contrast
    all_genes = set()
    for contrast in selected_contrasts:
        if (contrast.analysis_id in integrator.deg_data and
            contrast.contrast_id in integrator.deg_data[contrast.analysis_id]):

            df = integrator.deg_data[contrast.analysis_id][contrast.contrast_id]
            if 'Gene' in df.columns:
                all_genes.update(df['Gene'].tolist())

    all_genes = list(all_genes)

    if not all_genes:
        return []

    # Perform meta-analysis
    gene_meta_stats = meta_analysis_stouffer(integrator, selected_contrasts, all_genes)

    if not gene_meta_stats:
        return []

    # Apply FDR correction
    genes = list(gene_meta_stats.keys())
    p_values = np.array([gene_meta_stats[g]['meta_p'] for g in genes])
    fdr_corrected = benjamini_hochberg_correction(p_values)

    # Update gene stats with FDR
    for i, gene in enumerate(genes):
        gene_meta_stats[gene]['meta_fdr'] = fdr_corrected[i]

    # Filter by FDR threshold and minimum frequency
    significant_genes = []
    for gene, stats in gene_meta_stats.items():
        if (stats['meta_fdr'] < thresholds.fdr_cutoff and
            stats['total_sig'] >= thresholds.min_frequency):
            significant_genes.append(gene)

    # Rank by combined score: median_abs_logFC × -log10(meta_fdr)
    def ranking_score(gene):
        stats = gene_meta_stats[gene]
        if stats['meta_fdr'] > 0:
            return stats['median_abs_logfc'] * (-np.log10(stats['meta_fdr']))
        else:
            return stats['median_abs_logfc'] * 16  # Cap at -log10(1e-16)

    ranked_genes = sorted(significant_genes, key=ranking_score, reverse=True)

    return ranked_genes[:max_genes]

def create_landing_heatmap(integrator, top_genes: List[str], selected_contrasts: List[ContrastSelection], thresholds: ThresholdSelection) -> Optional[go.Figure]:
    """Create heatmap for landing page."""

    try:
        contrast_pairs = [(c.analysis_id, c.contrast_id) for c in selected_contrasts]

        # Use standard visualization thresholds for heatmap coloring
        # Genes are already selected via meta-analysis, so use reasonable display thresholds
        viz_p_threshold = 0.05
        viz_lfc_threshold = 1.0

        fig = integrator.create_lfc_heatmap(
            genes=top_genes,
            contrasts=contrast_pairs,
            output_file=None,
            p_value_threshold=viz_p_threshold,
            lfc_threshold=viz_lfc_threshold,
            hide_empty_rows_cols=True,
            font_size=11
        )

        if fig:
            fig.update_layout(
                title="Meta-Analysis Selected Genes Across Contrasts",
                title_font_size=16,
                height=max(400, min(800, len(top_genes) * 20))
            )

        return fig
    except Exception as e:
        st.error(f"Error creating heatmap: {str(e)}")
        return None

def create_gene_table(integrator, top_genes: List[str], selected_contrasts: List[ContrastSelection], thresholds: ThresholdSelection) -> pd.DataFrame:
    """Create gene summary table with meta-analysis results."""

    # Perform meta-analysis for the top genes
    gene_meta_stats = meta_analysis_stouffer(integrator, selected_contrasts, top_genes)

    # Apply FDR correction
    if gene_meta_stats:
        genes = list(gene_meta_stats.keys())
        p_values = np.array([gene_meta_stats[g]['meta_p'] for g in genes])
        fdr_corrected = benjamini_hochberg_correction(p_values)

        for i, gene in enumerate(genes):
            gene_meta_stats[gene]['meta_fdr'] = fdr_corrected[i]

    gene_data = []

    for gene in top_genes:
        if gene in gene_meta_stats:
            stats = gene_meta_stats[gene]

            # Create direction summary
            up_count = stats['up_count']
            down_count = stats['down_count']
            total_tested = stats['n_contrasts']
            total_sig = stats['total_sig']

            direction_summary = f"↑{up_count}/↓{down_count} (of {total_tested})"

            # Calculate median logFC (signed, not absolute)
            median_logfc = np.median(stats['log_fcs']) if stats['log_fcs'] else 0

            gene_data.append({
                'Gene': gene,
                'Meta-LogFC': round(median_logfc, 2),
                'Meta-FDR': f"{stats['meta_fdr']:.2e}" if stats['meta_fdr'] < 0.01 else f"{stats['meta_fdr']:.3f}",
                'Direction': direction_summary,
                'Tested': total_tested
            })

    df = pd.DataFrame(gene_data)
    if not df.empty:
        # Sort by ranking score: median_abs_logFC × -log10(meta_fdr)
        def ranking_score(row):
            fdr = float(row['Meta-FDR']) if isinstance(row['Meta-FDR'], str) else row['Meta-FDR']
            if fdr > 0:
                return abs(row['Meta-LogFC']) * (-np.log10(fdr))
            else:
                return abs(row['Meta-LogFC']) * 16

        df['_ranking_score'] = df.apply(ranking_score, axis=1)
        df = df.sort_values('_ranking_score', ascending=False)
        df = df.drop('_ranking_score', axis=1)

    return df

def generate_narrative(selected_contrasts: List[ContrastSelection], thresholds: ThresholdSelection, top_genes: List[str], biological_prompt: str) -> str:
    """Generate interpretive narrative."""

    try:
        if LANDING_PAGE_AVAILABLE and os.getenv("OPENAI_API_KEY"):
            return generate_narrative_with_ai(selected_contrasts, thresholds, top_genes, biological_prompt)
        else:
            return generate_narrative_fallback(selected_contrasts, thresholds, top_genes, biological_prompt)
    except Exception:
        return generate_narrative_fallback(selected_contrasts, thresholds, top_genes, biological_prompt)

def generate_narrative_with_ai(selected_contrasts: List[ContrastSelection], thresholds: ThresholdSelection, top_genes: List[str], biological_prompt: str) -> str:
    """Generate narrative using AI."""

    client = OpenAI()

    contrast_summaries = []
    for contrast in selected_contrasts:
        contrast_summaries.append(f"- {contrast.contrast_id}: (Score: {contrast.relevance_score:.1f})")

    gene_list_str = ', '.join(top_genes[:10])
    if len(top_genes) > 10:
        gene_list_str += f' and {len(top_genes)-10} others'

    prompt = f"""
You are an expert computational biologist writing a clear, accessible summary of RNA-seq differential expression results.

Research Context: {biological_prompt}

Analysis Summary:
- {len(selected_contrasts)} biologically relevant contrasts were automatically selected using semantic similarity and evidence-aware ranking
- Meta-analysis performed using weighted Stouffer Z-method with FDR correction
- Statistical thresholds: Meta-FDR < {thresholds.fdr_cutoff}, minimum frequency ≥ {thresholds.min_frequency}
- {len(top_genes)} key differentially expressed genes identified through meta-analysis

Selected Contrasts:
{chr(10).join(contrast_summaries)}

Top Genes: {gene_list_str}

Write a 2-3 paragraph narrative that:
1. Summarizes the key findings in accessible language
2. Highlights the most important genes and patterns discovered through meta-analysis
3. Explains the strength of evidence from combining multiple contrasts
4. Provides biological context and potential significance
5. Maintains scientific accuracy while being readable

Focus on biological insights and the robustness of findings across multiple experiments.
"""

    response = client.chat.completions.create(
        model="gpt-4o-mini",
        messages=[{"role": "user", "content": prompt}],
        max_tokens=1000,
        temperature=0  # Deterministic output
    )

    return response.choices[0].message.content.strip()

def generate_narrative_fallback(selected_contrasts: List[ContrastSelection], thresholds: ThresholdSelection, top_genes: List[str], biological_prompt: str) -> str:
    """Fallback narrative generation."""

    narrative = f"""
This analysis automatically identified {len(selected_contrasts)} biologically relevant experimental contrasts related to {biological_prompt} using a two-stage selection process combining semantic similarity and evidence-aware ranking. A meta-analysis was performed using the weighted Stouffer Z-method with sample size weighting, followed by FDR correction to control for multiple testing (meta-FDR < {thresholds.fdr_cutoff}).

The {len(top_genes)} most consistently dysregulated genes include {', '.join(top_genes[:5])}{'and others' if len(top_genes) > 5 else ''}. These genes were identified through meta-analysis across multiple experimental conditions, providing robust evidence for their involvement in the biological processes under investigation. The meta-analytic approach captures both strong effects in individual studies and weaker but consistent patterns across experiments.

The results provide a statistically rigorous foundation for understanding the molecular mechanisms underlying the observed phenotypes. The combination of evidence across multiple contrasts strengthens confidence in these findings and suggests high-priority targets for further experimental validation.
"""

    return narrative

# Check if fragment is available (Streamlit >=1.33.0)
# If not, fallback to experimental_fragment
try:
    from streamlit import fragment
    st.fragment = fragment
except ImportError:
    try:
        from streamlit import experimental_fragment
        st.fragment = experimental_fragment
    except ImportError:
        # Fallback for very old Streamlit versions
        def fragment(func):
            """Fallback decorator when st.fragment is not available."""
            def wrapper(*args, **kwargs):
                return func(*args, **kwargs)
            return wrapper
        st.fragment = fragment

    # Add a cache for figure objects to improve performance in older versions
    @st.cache_data(hash_funcs={go.Figure: lambda _: None})
    def cached_figure_creation(func_name, *args, **kwargs):
        """Cache figure objects to avoid recreating them."""
        if func_name == "create_lfc_heatmap":
            return ri.create_lfc_heatmap(*args, **kwargs)
        elif func_name == "create_expression_plots":
            return ri.create_expression_plots(*args, **kwargs)
        return None

# Set page configuration
st.set_page_config(
    page_title="UORCA Explorer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Add CSS for dataframe text wrapping and compact multiselect tags
st.markdown("""
<style>
  /* Fix for text wrapping in dataframes */
  .stDataFrame tbody tr td {
    white-space: normal !important;
    word-wrap: break-word !important;
    max-width: 300px;
  }
  .stDataFrame th {
    white-space: normal !important;
    word-wrap: break-word !important;
    max-width: 300px;
  }
  /* Ensure Description column has more space */
  .stDataFrame td:nth-child(3) {
    min-width: 250px;
    max-width: 500px;
  }
  /* Compact multiselect tags */
  span[data-baseweb="tag"] {
    font-size: 11px !important;
    padding: 0.25rem 0.5rem !important;
    height: 1.2rem !important;
  }
</style>
""", unsafe_allow_html=True)

# ---------- 1. sidebar: data & user parameters -----------------
st.sidebar.title("🧬 UORCA Explorer")

# Get the default results directory
default_dir = os.getenv("UORCA_DEFAULT_RESULTS_DIR")
if not default_dir:
    default_dir = os.path.join(os.path.dirname(os.path.dirname(script_dir)), "UORCA_results")
    if not os.path.exists(default_dir):
        default_dir = os.path.dirname(os.path.dirname(script_dir))

# absolute path on the server
results_dir = st.sidebar.text_input(
    "Results directory",
    value=default_dir
)

# expensive → cache resource (keeps the object alive for the whole session)
@st.cache_resource
def get_integrator(path):
    try:
        ri = ResultsIntegrator(results_dir=path)
        ri.load_data()
        return ri, None
    except Exception as e:
        return None, str(e)

# Load the integrator (cached)
with st.sidebar.status("Loading data...", expanded=True) as status:
    ri, error = get_integrator(results_dir)

    if error:
        status.update(label=f"Error loading data: {error}", state="error")
    elif not ri or not ri.cpm_data:
        status.update(label="No data found. Please check the directory path.", state="error")
    else:
        status.update(label=f"✅ Loaded {len(ri.cpm_data)} datasets", state="complete")

# Show AI landing page availability status
api_key_available = bool(os.getenv("OPENAI_API_KEY"))

if LANDING_PAGE_AVAILABLE and api_key_available:
    st.sidebar.success("🤖 AI Landing Page: Fully Available")
elif LANDING_PAGE_AVAILABLE and not api_key_available:
    st.sidebar.warning("🤖 AI Landing Page: Partial (Heuristic Mode)")
    with st.sidebar.expander("ℹ️ Enable Full AI Features", expanded=False):
        st.markdown("""
        **Current Status:** OpenAI package installed ✅, but API key missing ⚠️

        **Available:** Heuristic contrast selection and gene ranking
        **Missing:** AI-powered relevance scoring and narrative generation

        **To enable full AI features:**

        Option 1 - Environment variable:
        ```bash
        export OPENAI_API_KEY=your_key_here
        ```

        Option 2 - .env file (recommended):
        Create a `.env` file in the project root:
        ```
        OPENAI_API_KEY=your_key_here
        ```

        The landing page works without API key but uses simpler algorithms.
        """)
else:
    st.sidebar.error("🤖 AI Landing Page: Not Available")
    with st.sidebar.expander("ℹ️ Enable AI Features", expanded=False):
        st.markdown("""
        **Current Status:** OpenAI package not available ❌

        **To enable AI landing page:**
        ```bash
        uv add openai
        ```

        Then set your API key:
        - Environment: `export OPENAI_API_KEY=your_key`
        - .env file: Add `OPENAI_API_KEY=your_key` to `.env`

        AI features provide:
        - Automatic contrast relevance scoring
        - Optimal threshold selection
        - Biological narrative generation
        - Interpretive summaries
        """)

# Only show the rest of the UI if we successfully loaded data
if ri and ri.cpm_data:
    # Get all genes across all datasets
    all_genes = set()
    for cpm_df in ri.cpm_data.values():
        if 'Gene' in cpm_df.columns:
            all_genes.update(cpm_df['Gene'].tolist())
    all_genes = sorted(all_genes)

    # Advanced options (moved up so variables are available throughout)
    st.sidebar.subheader("Advanced Options")
    show_advanced = st.sidebar.checkbox("Show advanced options", value=False)

    # Initialize default values
    pvalue_thresh = 0.05
    lfc_thresh = 1.0
    min_unique = 1
    max_contrasts_unique = 2
    top_frequent_genes = 20
    top_unique_genes = 10
    hide_x_labels = True
    debug_mode = False

    if show_advanced:
        # Add slider and number input side by side for precision
        col1, col2 = st.sidebar.columns([3, 1])
        with col1:
            pvalue_thresh = st.slider("Adjusted P-value threshold", 0.0, 0.1, 0.05, 0.001, help="Genes must have adj.P.Val below this threshold to be considered significant")
        with col2:
            pvalue_text = st.text_input("", value=f"{pvalue_thresh:.3f}", key="pvalue_text")
            try:
                pvalue_thresh = float(pvalue_text)
            except ValueError:
                pass

        col1, col2 = st.sidebar.columns([3, 1])
        with col1:
            lfc_thresh = st.slider("Log2FC threshold", 0.0, 5.0, 1.0, 0.1, help="Absolute log2 fold change threshold")
        with col2:
            lfc_text = st.text_input("", value=f"{lfc_thresh:.1f}", key="lfc_text")
            try:
                lfc_thresh = float(lfc_text)
            except ValueError:
                pass

        st.sidebar.subheader("Auto-selection Parameters")
        st.sidebar.markdown("**Auto-selection uses two strategies:**")
        st.sidebar.markdown("1. **Frequent DEGs**: Genes significant across multiple contrasts")
        st.sidebar.markdown("2. **Contrast-specific DEGs**: High fold-change genes unique to few contrasts")

        col1, col2 = st.sidebar.columns([3, 1])
        with col1:
            top_frequent_genes = st.slider("Frequent DEGs to include", 5, 50, 20, help="Number of genes that are differentially expressed across the most contrasts")
        with col2:
            freq_text = st.text_input("", value=str(top_frequent_genes), key="freq_text")
            try:
                top_frequent_genes = int(freq_text)
            except ValueError:
                pass

        col1, col2 = st.sidebar.columns([3, 1])
        with col1:
            top_unique_genes = st.slider("Contrast-specific DEGs per contrast", 1, 20, 10, help="Number of high fold-change genes to select from each contrast that appear in few other contrasts")
        with col2:
            unique_text = st.text_input("", value=str(top_unique_genes), key="unique_text_simple")
            try:
                top_unique_genes = int(unique_text)
            except ValueError:
                pass

        col1, col2 = st.sidebar.columns([3, 1])
        with col1:
            max_contrasts_unique = st.slider("Max contrasts for 'contrast-specific'", 1, 10, 2, help="A gene is considered 'contrast-specific' if it appears as significant in this many contrasts or fewer")
        with col2:
            max_text = st.text_input("", value=str(max_contrasts_unique), key="max_text")
            try:
                max_contrasts_unique = int(max_text)
            except ValueError:
                pass

        col1, col2 = st.sidebar.columns([3, 1])
        with col1:
            min_unique = st.slider("Min contrast-specific DEGs per contrast", 0, 10, 1, help="Minimum number of contrast-specific genes that must be selected from each contrast (quality control)")
        with col2:
            min_text = st.text_input("", value=str(min_unique), key="min_text")
            try:
                min_unique = int(min_text)
            except ValueError:
                pass

        st.sidebar.subheader("Visualization Options")
        hide_x_labels = st.sidebar.checkbox("Hide x-axis labels in expression plots", value=True)

        st.sidebar.subheader("Heatmap Options")
        use_separate_heatmap_filters = st.sidebar.checkbox("Use different significance filters for heatmap", value=False, help="By default, heatmap uses same filters as Auto-selected DEGs. Enable this to set different thresholds.")

        if use_separate_heatmap_filters:
            st.sidebar.markdown("**Heatmap-specific filters:**")
            col1, col2 = st.sidebar.columns([3, 1])
            with col1:
                heatmap_pvalue_thresh = st.slider("Heatmap P-value threshold", 0.0, 0.1, 0.05, 0.001, help="P-value threshold for heatmap coloring")
            with col2:
                heatmap_pvalue_text = st.text_input("", value=f"{heatmap_pvalue_thresh:.3f}", key="heatmap_pvalue_text")
                try:
                    heatmap_pvalue_thresh = float(heatmap_pvalue_text)
                except ValueError:
                    pass

            col1, col2 = st.sidebar.columns([3, 1])
            with col1:
                heatmap_lfc_thresh = st.slider("Heatmap LFC threshold", 0.0, 5.0, 1.0, 0.1, help="LFC threshold for heatmap coloring")
            with col2:
                heatmap_lfc_text = st.text_input("", value=f"{heatmap_lfc_thresh:.1f}", key="heatmap_lfc_text")
                try:
                    heatmap_lfc_thresh = float(heatmap_lfc_text)
                except ValueError:
                    pass

            use_dynamic_filtering = True
            effective_pvalue_thresh = heatmap_pvalue_thresh
            effective_lfc_thresh = heatmap_lfc_thresh
        else:
            use_dynamic_filtering = True
            effective_pvalue_thresh = pvalue_thresh
            effective_lfc_thresh = lfc_thresh
            st.sidebar.info("🎯 Heatmap uses same significance filters as Auto-selected DEGs")

        hide_empty_rows_cols = st.sidebar.checkbox("Hide genes/contrasts with no significant values", value=False, help="Remove rows/columns where no values meet significance criteria")
        if hide_empty_rows_cols:
            st.sidebar.info("🧹 Genes/contrasts with no significant values will be completely removed from the heatmap")

        debug_mode = st.sidebar.checkbox("Debug mode", value=False)

    # Initialize session state for selections
    if 'datasets_selected' not in st.session_state:
        st.session_state['datasets_selected'] = None

    if 'contrasts_selected' not in st.session_state:
        st.session_state['contrasts_selected'] = None

    # Gene selection
    st.sidebar.subheader("Gene Selection")
    gene_select_method = st.sidebar.radio(
        "Selection method:",
        ["Auto-selected DEGs", "Custom"],
        index=0  # Default to "Auto-selected DEGs"
    )

    if gene_select_method == "Custom":
        st.sidebar.write("Enter genes (one per line):")
        gene_input = st.sidebar.text_area(
            "Gene list:",
            height=150,
            placeholder="MYCN\nALK\nPHOX2B\nATRX\nTP53\nTERT\nARID1A\nARID1B\nNF1\nBARD1",
            label_visibility="collapsed"
        )

        if gene_input.strip():
            # Parse genes from textarea (one per line)
            input_genes = [gene.strip() for gene in gene_input.strip().split('\n') if gene.strip()]
            # Filter to only genes that exist in the data
            gene_sel = [gene for gene in input_genes if gene in all_genes]

            # Show validation info
            if len(input_genes) != len(gene_sel):
                missing_genes = [gene for gene in input_genes if gene not in all_genes]
                st.sidebar.warning(f"⚠️ {len(missing_genes)} genes not found in data: {', '.join(missing_genes[:5])}{' ...' if len(missing_genes) > 5 else ''}")

            if gene_sel:
                st.sidebar.success(f"✅ {len(gene_sel)} genes selected")
                st.sidebar.caption(f"🧬 {len(gene_sel)} genes selected")
            else:
                st.sidebar.error("❌ No valid genes found")
        else:
            gene_sel = []
    elif gene_select_method == "Auto-selected DEGs":
        # Use p-value and lfc thresholds (these will be defined either here or in advanced options)
        if not show_advanced:
            # Add slider and number input side by side for precision
            col1, col2 = st.sidebar.columns([3, 1])
            with col1:
                pvalue_thresh = st.slider("Adjusted P-value threshold", 0.0, 0.1, 0.05, 0.001, help="Use adjusted P-value (adj.P.Val) for filtering DEGs")
            with col2:
                pvalue_text = st.text_input("", value=f"{pvalue_thresh:.3f}", key="pvalue_text_simple")
                try:
                    pvalue_thresh = float(pvalue_text)
                except ValueError:
                    pass

            col1, col2 = st.sidebar.columns([3, 1])
            with col1:
                lfc_thresh = st.slider("abs(Log2FC) threshold", 0.0, 5.0, 1.0, 0.1, help="Genes must have absolute log2 fold change above this threshold to be considered significant")
            with col2:
                lfc_text = st.text_input("", value=f"{lfc_thresh:.1f}", key="lfc_text_simple")
                try:
                    lfc_thresh = float(lfc_text)
                except ValueError:
                    pass

            # In simple mode, always use same filters for heatmap as important genes
            use_dynamic_filtering = True
            effective_pvalue_thresh = pvalue_thresh
            effective_lfc_thresh = lfc_thresh
            hide_empty_rows_cols = st.sidebar.checkbox("Hide genes/contrasts with no significant values", value=False, help="Remove rows/columns where no values meet significance criteria")
            if hide_empty_rows_cols:
                st.sidebar.info("🧹 Genes/contrasts with no significant values will be completely removed from the heatmap")
            st.sidebar.info("🎯 Heatmap uses same significance filters as Auto-selected DEGs")

            # Auto-selection parameters
            st.sidebar.markdown("**Auto-selection uses two strategies:**")
            st.sidebar.markdown("1. **Frequent DEGs**: Genes significant across multiple contrasts")
            st.sidebar.markdown("2. **Contrast-specific DEGs**: High fold-change genes unique to few contrasts")

            col1, col2 = st.sidebar.columns([3, 1])
            with col1:
                top_frequent_genes = st.slider("Frequent DEGs to include", 5, 50, 20, help="Number of genes that are differentially expressed across the most contrasts")
            with col2:
                freq_text = st.text_input("", value=str(top_frequent_genes), key="freq_text_simple")
                try:
                    top_frequent_genes = int(freq_text)
                except ValueError:
                    pass

            col1, col2 = st.sidebar.columns([3, 1])
            with col1:
                top_unique_genes = st.slider("Contrast-specific DEGs per contrast", 1, 20, 10, help="Number of high fold-change genes to select from each contrast that appear in few other contrasts")
            with col2:
                unique_text = st.text_input("", value=str(top_unique_genes), key="unique_text_simple")
                try:
                    top_unique_genes = int(unique_text)
                except ValueError:
                    pass

            col1, col2 = st.sidebar.columns([3, 1])
            with col1:
                max_contrasts_unique = st.slider("Max contrasts for 'contrast-specific'", 1, 10, 2, help="A gene is considered 'contrast-specific' if it appears as significant in this many contrasts or fewer")
            with col2:
                max_text = st.text_input("", value=str(max_contrasts_unique), key="max_text_simple")
                try:
                    max_contrasts_unique = int(max_text)
                except ValueError:
                    pass

        # Get all significant genes
        top_genes = ri.identify_important_genes(
            top_frequent=top_frequent_genes,
            top_unique=top_unique_genes,
            max_contrasts_for_unique=max_contrasts_unique,
            min_unique_per_contrast=min_unique if show_advanced else 1,
            p_value_threshold=pvalue_thresh,
            lfc_threshold=lfc_thresh
        )
        top_genes = sorted(top_genes)
        col1, col2 = st.sidebar.columns([3, 1])
        with col2:
            st.write("")  # Add some spacing
            if st.button("Select All", key="select_all_genes_important") and limited_genes:
                st.session_state['gene_sel_important'] = limited_genes  # Use limited and ranked genes
                try:
                    st.rerun()
                except AttributeError:
                    try:
                        st.experimental_rerun()
                    except AttributeError:
                        pass

        # Limit genes to 200 maximum, selecting by highest LFC if needed
        if len(top_genes) > 200:
            # Get LFC data for ranking
            gene_lfc_map = {}
            for analysis_id, contrasts in ri.deg_data.items():
                for contrast_id, df in contrasts.items():
                    if 'Gene' in df.columns and 'logFC' in df.columns:
                        for _, row in df.iterrows():
                            gene = row['Gene']
                            if gene in top_genes:
                                current_max = gene_lfc_map.get(gene, 0)
                                gene_lfc_map[gene] = max(current_max, abs(row['logFC']))

            # Sort by highest LFC and take top 200
            sorted_genes = sorted(top_genes, key=lambda g: gene_lfc_map.get(g, 0), reverse=True)
            limited_genes = sorted_genes[:200]
            st.sidebar.warning(f"⚠️ Showing top 200 of {len(top_genes)} important genes (ranked by highest LFC). {len(top_genes)-200} genes excluded.")
        else:
            limited_genes = top_genes

        # Use session state if "Select All" was clicked or use default selection (all limited genes)
        default_genes = st.session_state.get('gene_sel_important', limited_genes if limited_genes else [])
        gene_sel = st.sidebar.multiselect(
            "Auto-selected DEGs:",
            limited_genes,
            default=default_genes,
            max_selections=200
        )

        # Add count badge
        if gene_sel:
            st.sidebar.caption(f"🧬 {len(gene_sel)} genes selected")

        # Clear the session state after use
        if 'gene_sel_important' in st.session_state:
            del st.session_state['gene_sel_important']



    # Visualization options
    st.sidebar.subheader("Plot Options")
    plot_type = "violin"  # Always use violin plots

    # Add pagination controls if needed
    if gene_sel:
        genes_per_page = 30
        total_pages = (len(gene_sel) + genes_per_page - 1) // genes_per_page

        if total_pages > 1:
            st.sidebar.markdown("---")
            st.sidebar.subheader("Gene Pagination")

            st.sidebar.info(f"Your selection contains {len(gene_sel)} genes, showing {genes_per_page} per page.")

            # Create three columns for pagination controls
            col1, col2, col3 = st.sidebar.columns([1, 2, 1])

            with col1:
                prev_disabled = 'page_num' in st.session_state and st.session_state.page_num <= 1
                if st.button("◀", disabled=prev_disabled, key="prev_page"):
                    st.session_state.page_num = max(1, st.session_state.get('page_num', 1) - 1)
                    try:
                        st.rerun()
                    except AttributeError:
                        # Fallback for older Streamlit versions
                        try:
                            st.experimental_rerun()
                        except AttributeError:
                            import streamlit.runtime.scriptrunner.magic as _m
                            _m._set_stop_thread(False)  # hidden API for very old versions

            with col2:
                # Initialize page_num in session state if not present
                if 'page_num' not in st.session_state:
                    st.session_state.page_num = 1

                page_num = st.number_input(
                    "Page",
                    min_value=1,
                    max_value=total_pages,
                    value=st.session_state.page_num,
                    step=1,
                    key="page_input"
                )
                st.session_state.page_num = page_num

            with col3:
                next_disabled = 'page_num' in st.session_state and st.session_state.page_num >= total_pages
                if st.button("▶", disabled=next_disabled, key="next_page"):
                    st.session_state.page_num = min(total_pages, st.session_state.get('page_num', 1) + 1)
                    try:
                        st.rerun()
                    except AttributeError:
                        # Fallback for older Streamlit versions
                        try:
                            st.experimental_rerun()
                        except AttributeError:
                            import streamlit.runtime.scriptrunner.magic as _m
                            _m._set_stop_thread(False)  # hidden API for very old versions

            # Show page indicator
            st.sidebar.caption(f"Page {st.session_state.page_num} of {total_pages}")
        else:
            if 'page_num' not in st.session_state:
                st.session_state.page_num = 1
    else:
        if 'page_num' not in st.session_state:
            st.session_state.page_num = 1

    # These options have been moved up in the UI

    # ---------- 2. main tabs ---------------------------------------
    tab_landing, tab_sel, tab1, tab2, tab3, tab4, tab5 = st.tabs(["🎨 Landing Page", "📊 Selections", "🌡️ Heat-map", "🎻 Expression", "📈 Analysis Plots", "📋 Dataset Info", "🔍 Contrast Info"])

    # Initialize session state for selections if not exists
    if 'selected_datasets' not in st.session_state:
        # Default to first 5 datasets
        all_dataset_ids = list(ri.cpm_data.keys())
        st.session_state['selected_datasets'] = set(all_dataset_ids[:5])

    if 'selected_contrasts' not in st.session_state:
        # Default to all contrasts for selected datasets
        selected_contrasts = set()
        for analysis_id, contrasts in ri.deg_data.items():
            if analysis_id in st.session_state['selected_datasets']:
                for contrast_id in contrasts.keys():
                    selected_contrasts.add((analysis_id, contrast_id))
        st.session_state['selected_contrasts'] = selected_contrasts



    with tab_landing:
        st.header("🎨 AI-Assisted Landing Page")
        st.markdown("""
        **Generate intelligent summaries of your RNA-seq results with minimal interaction.**

        This AI-powered landing page automatically selects the most relevant contrasts, determines optimal statistical thresholds,
        and creates interpretive narratives for your differential expression analysis.
        """)

        # Landing page controls
        col1, col2 = st.columns([2, 1])

        with col1:
            biological_prompt = st.text_area(
                "🧬 Research Context",
                value="General differential expression analysis",
                height=100,
                help="Describe your biological research question. This guides the AI in selecting relevant contrasts and generating interpretations.",
                placeholder="e.g., 'MYCN amplification effects in neuroblastoma tumors' or 'immune checkpoint response in melanoma'"
            )

        with col2:
            st.markdown("**⚙️ Options**")
            max_contrasts_landing = st.slider("Max contrasts", 3, 15, 8, key="landing_max_contrasts")
            max_genes_landing = st.slider("Max genes", 20, 100, 50, key="landing_max_genes")

            generate_landing = st.button(
                "🚀 Generate AI Landing Page",
                type="primary",
                use_container_width=True,
                help="Uses AI to automatically select relevant contrasts and create interpretive summaries"
            )

        # Check if landing page functionality is available
        if not LANDING_PAGE_AVAILABLE:
            st.error("❌ AI landing page functionality requires OpenAI package and landing page module.")
            st.info("💡 To enable: `pip install openai` and ensure landing_page_generator.py is available.")
        elif not ri or not ri.cpm_data:
            st.info("📁 Please load data first using the sidebar controls.")
        elif generate_landing:
            # Generate AI-assisted landing page
            with st.spinner("🤖 AI is analyzing your data and generating the landing page..."):
                try:
                    # Generate landing page data using integrated functionality
                    landing_data = generate_ai_landing_page(
                        integrator=ri,
                        biological_prompt=biological_prompt,
                        max_contrasts=max_contrasts_landing,
                        max_genes=max_genes_landing
                    )

                    if landing_data:
                        # Store in session state
                        st.session_state.landing_data = landing_data
                    else:
                        st.error("❌ Failed to generate landing page - no suitable data found.")

                except Exception as e:
                    st.error(f"❌ Error generating landing page: {str(e)}")
                    if "openai" in str(e).lower():
                        st.info("💡 Tip: Make sure your OpenAI API key is set as an environment variable: `export OPENAI_API_KEY=your_key`")

        # Display landing page results if available
        if hasattr(st.session_state, 'landing_data') and st.session_state.landing_data:
            landing_data = st.session_state.landing_data

            st.success("✅ AI landing page generated successfully!")

            # Summary metrics
            col1, col2, col3, col4 = st.columns(4)
            with col1:
                st.metric("🔍 Selected Contrasts", len(landing_data.selected_contrasts))
            with col2:
                st.metric("🧬 Meta-Analysis Genes", len(landing_data.top_genes))
            with col3:
                st.metric("📊 Meta-FDR Cutoff", f"{landing_data.thresholds.fdr_cutoff}")
            with col4:
                st.metric("🔄 Min Frequency", f"{landing_data.thresholds.min_frequency}")

            # Key findings narrative
            st.markdown("### 📝 AI-Generated Key Findings")
            st.markdown(f"""
            <div style="background: #e8f4fd; padding: 1.5rem; border-radius: 0.5rem; border-left: 4px solid #3498db; margin: 1rem 0;">
            {landing_data.narrative}
            </div>
            """, unsafe_allow_html=True)

            # Analysis parameters
            st.markdown("### 🎯 Automatic Parameter Selection")
            st.markdown(f"""
            <div style="background: #f0f8e8; padding: 1rem; border-radius: 0.5rem; border-left: 4px solid #27ae60; margin: 1rem 0;">
            <strong>Meta-Analysis Parameters:</strong> Meta-FDR < {landing_data.thresholds.fdr_cutoff}, Min. frequency ≥ {landing_data.thresholds.min_frequency} contrasts<br>
            <strong>Method:</strong> Weighted Stouffer Z-method with sample size weighting and FDR correction<br>
            <strong>AI Justification:</strong> {landing_data.thresholds.justification}
            </div>
            """, unsafe_allow_html=True)

            # Tabbed results view
            landing_tab1, landing_tab2, landing_tab3 = st.tabs(["🔍 AI-Selected Contrasts", "🧬 Top Genes", "🌡️ Expression Heatmap"])

            with landing_tab1:
                st.markdown("#### Contrasts Selected by AI with Justifications")

                # Create contrast DataFrame
                contrast_df = pd.DataFrame([
                    {
                        "Dataset": c.analysis_id,
                        "Contrast": c.contrast_id,
                        "Relevance Score": f"{c.relevance_score:.1f}",
                        "DEG Count": c.deg_count,
                        "AI Justification": c.justification
                    }
                    for c in landing_data.selected_contrasts
                ])

                st.dataframe(
                    contrast_df,
                    use_container_width=True,
                    column_config={
                        "AI Justification": st.column_config.TextColumn(
                            "AI Justification",
                            width="large"
                        ),
                        "Relevance Score": st.column_config.TextColumn("Relevance Score", width="small"),
                        "DEG Count": st.column_config.NumberColumn("DEG Count", width="small")
                    }
                )

                # Auto-apply selections to other tabs
                if st.button("📋 Apply These Selections to Other Tabs", type="secondary"):
                    # Update session state for other tabs
                    selected_datasets = set(c.analysis_id for c in landing_data.selected_contrasts)
                    selected_contrasts = set((c.analysis_id, c.contrast_id) for c in landing_data.selected_contrasts)

                    st.session_state['selected_datasets'] = selected_datasets
                    st.session_state['selected_contrasts'] = selected_contrasts

                    st.success("✅ Selections applied! Switch to Heat-map or Expression tabs to see visualizations with these selections.")

            with landing_tab2:
                st.markdown("#### Top Differentially Expressed Genes")

                if not landing_data.gene_table.empty:
                    st.dataframe(
                        landing_data.gene_table,
                        use_container_width=True,
                        column_config={
                            "Gene": st.column_config.TextColumn("Gene Symbol", width="medium"),
                            "Meta-LogFC": st.column_config.NumberColumn("Meta Log2FC", format="%.2f", help="Median log2 fold change across contrasts"),
                            "Meta-FDR": st.column_config.TextColumn("Meta FDR", help="FDR-corrected meta-analysis p-value"),
                            "Direction": st.column_config.TextColumn("↑/↓ Pattern", help="Up/down regulation pattern across contrasts"),
                            "Tested": st.column_config.NumberColumn("Tested", help="Number of contrasts where gene was tested")
                        }
                    )

                    # Gene list for easy copying
                    with st.expander("📋 Gene List (for external tools)", expanded=False):
                        gene_list_text = '\n'.join(landing_data.top_genes)
                        st.text_area(
                            "Copy gene symbols:",
                            value=gene_list_text,
                            height=150,
                            help="Gene symbols, one per line. Copy and paste into other tools like DAVID, Enrichr, etc."
                        )
                else:
                    st.info("No genes met the AI-selected significance criteria.")

            with landing_tab3:
                st.markdown("#### AI-Generated Expression Heatmap")

                if landing_data.heatmap_fig:
                    st.plotly_chart(landing_data.heatmap_fig, use_container_width=True)

                    st.markdown("""
                    💡 **How to interpret this meta-analysis heatmap:**
                    - Each row represents a gene selected through meta-analysis, each column represents a biological contrast
                    - Colors indicate log2 fold change: red (upregulated), blue (downregulated), white (not significant in individual study)
                    - Hover over cells to see detailed information including contrast descriptions
                    - Genes and contrasts are automatically clustered to reveal biological patterns
                    - Genes were selected based on meta-analysis significance across multiple contrasts
                    - Individual cell coloring uses standard thresholds (p<0.05, |logFC|>1.0) for visualization
                    """)
                else:
                    st.warning("Could not generate heatmap. This may occur if no genes meet the AI-selected significance criteria.")

            # Export options
            st.markdown("---")
            st.markdown("### 💾 Export AI Results")

            col1, col2 = st.columns(2)

            with col1:
                if not landing_data.gene_table.empty:
                    csv_data = landing_data.gene_table.to_csv(index=False)
                    st.download_button(
                        label="📊 Download Gene Table CSV",
                        data=csv_data,
                        file_name=f"ai_selected_genes_{pd.Timestamp.now().strftime('%Y%m%d_%H%M%S')}.csv",
                        mime="text/csv"
                    )

            with col2:
                # Create a summary report
                summary_data = {
                    'research_context': biological_prompt,
                    'selected_contrasts': len(landing_data.selected_contrasts),
                    'key_genes': len(landing_data.top_genes),
                    'thresholds': {
                        'fdr': landing_data.thresholds.fdr_cutoff,
                        'logfc': landing_data.thresholds.logfc_cutoff,
                        'min_frequency': landing_data.thresholds.min_frequency
                    },
                    'ai_narrative': landing_data.narrative,
                    'generated_at': pd.Timestamp.now().isoformat()
                }

                st.download_button(
                    label="📄 Download AI Summary JSON",
                    data=json.dumps(summary_data, indent=2),
                    file_name=f"ai_landing_summary_{pd.Timestamp.now().strftime('%Y%m%d_%H%M%S')}.json",
                    mime="application/json"
                )

        elif LANDING_PAGE_AVAILABLE:
            # Show helpful information about the landing page
            st.markdown("---")
            st.markdown("### 🚀 Get Started with AI-Assisted Analysis")

            col1, col2 = st.columns(2)

            with col1:
                st.markdown("""
                **🤖 What the AI Landing Page Does:**

                1. **Smart Contrast Selection**: AI analyzes all your contrasts and picks the most biologically relevant ones
                2. **Optimal Thresholds**: Automatically determines the best statistical cutoffs for your data
                3. **Gene Prioritization**: Ranks genes by importance across multiple contrasts
                4. **Biological Interpretation**: Generates accessible summaries of your findings
                5. **Interactive Visualizations**: Creates publication-ready plots with detailed annotations
                """)

            with col2:
                st.markdown("""
                **💡 Tips for Best Results:**

                - **Be Specific**: Describe your research focus (e.g., "MYCN amplification in neuroblastoma")
                - **Include Context**: Mention key pathways, genes, or processes of interest
                - **State Comparisons**: Specify what you're comparing (tumor vs normal, treated vs control)
                - **Add Goals**: Include what you hope to discover or validate

                **Example Prompts:**
                - "Immune checkpoint response in melanoma patients"
                - "Developmental gene expression in neural differentiation"
                - "Drug resistance mechanisms in cancer cell lines"
                """)

    with tab_sel:
        st.header("📊 Selections")
        st.markdown("**Select datasets and contrasts for analysis.** Use the tables below to choose your data, then click 'Regenerate Plots' to update visualizations.")

        # Dataset selection table
        st.subheader("Choose Datasets")
        dataset_rows = []
        for analysis_id, info in ri.analysis_info.items():
            title = getattr(ri, "dataset_info", {}).get(analysis_id, {}).get("title", "")
            title_display = title[:100] + ("..." if len(title) > 100 else "") if title else "No title available"

            dataset_rows.append({
                "✔": analysis_id in st.session_state.get('selected_datasets', set()),
                "Accession": info.get("accession", "Unknown"),
                "Dataset ID": analysis_id,
                "Title": title_display,
                "Organism": info.get("organism", "Unknown"),
                "Samples": info.get("number_of_samples", 0),
                "Contrasts": info.get("number_of_contrasts", 0)
            })

        if dataset_rows:
            ds_df = pd.DataFrame(dataset_rows)
            edited_ds = st.data_editor(
                ds_df,
                hide_index=True,
                use_container_width=True,
                column_config={
                    "✔": st.column_config.CheckboxColumn("Select", default=False),
                    "Title": st.column_config.TextColumn("Title", width="large"),
                    "Accession": st.column_config.TextColumn("Accession", width="medium"),
                    "Dataset ID": st.column_config.TextColumn("Dataset ID", width="medium")
                },
                key="selections_dataset_editor"
            )

            # Update selected datasets based on edited data
            if not edited_ds.empty:
                selected_datasets = set(edited_ds.loc[edited_ds["✔"], "Dataset ID"].tolist())
                st.session_state['selected_datasets'] = selected_datasets

                st.caption(f"📊 {len(selected_datasets)} datasets selected")

        # Contrast selection table
        st.subheader("Choose Contrasts")
        if st.session_state['selected_datasets']:
            contrast_rows = []
            for analysis_id, contrasts in ri.deg_data.items():
                if analysis_id in st.session_state['selected_datasets']:
                    for contrast_id in contrasts.keys():
                        # Get description
                        description = ri._get_contrast_description(analysis_id, contrast_id)
                        if description.startswith("Contrast: "):
                            description = description[10:]  # Remove "Contrast: " prefix

                        # Count DEGs
                        deg_count = 0
                        if analysis_id in ri.deg_data and contrast_id in ri.deg_data[analysis_id]:
                            df = ri.deg_data[analysis_id][contrast_id]
                            if 'adj.P.Val' in df.columns and 'logFC' in df.columns:
                                deg_count = ((df['adj.P.Val'] < pvalue_thresh) & (abs(df['logFC']) > lfc_thresh)).sum()

                        contrast_rows.append({
                            "✔": (analysis_id, contrast_id) in st.session_state.get('selected_contrasts', set()),
                            "Dataset": analysis_id,
                            "Accession": ri.analysis_info.get(analysis_id, {}).get("accession", "Unknown"),
                            "Contrast": contrast_id,
                            "Description": description[:150] + ("..." if len(description) > 150 else ""),
                            "DEGs": deg_count
                        })

            if contrast_rows:
                ctr_df = pd.DataFrame(contrast_rows)
                edited_ctr = st.data_editor(
                    ctr_df,
                    hide_index=True,
                    use_container_width=True,
                    column_config={
                        "✔": st.column_config.CheckboxColumn("Select", default=False),
                        "Description": st.column_config.TextColumn("Description", width="large"),
                        "Dataset": st.column_config.TextColumn("Dataset", width="medium"),
                        "Contrast": st.column_config.TextColumn("Contrast", width="medium"),
                        "DEGs": st.column_config.NumberColumn("DEGs", format="%d")
                    },
                    key="selections_contrast_editor"
                )

                # Update selected contrasts based on edited data
                if not edited_ctr.empty:
                    selected_contrasts = set()
                    for _, row in edited_ctr.iterrows():
                        if row["✔"]:
                            selected_contrasts.add((row["Dataset"], row["Contrast"]))
                    st.session_state['selected_contrasts'] = selected_contrasts

                    st.caption(f"🔍 {len(selected_contrasts)} contrasts selected")
            else:
                st.info("No contrasts available for selected datasets.")
        else:
            st.info("Please select at least one dataset to see available contrasts.")

        # Selection summary
        st.markdown("---")
        selected_datasets_count = len(st.session_state.get('selected_datasets', set()))
        selected_contrasts_count = len(st.session_state.get('selected_contrasts', set()))
        genes_count = len(gene_sel) if gene_sel else 0

        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Datasets", selected_datasets_count)
        with col2:
            st.metric("Contrasts", selected_contrasts_count)
        with col3:
            st.metric("Genes", genes_count)

        if selected_datasets_count > 0 and selected_contrasts_count > 0 and genes_count > 0:
            st.success("✅ Ready for analysis! Switch to Heat-map or Expression tabs to view results.")
        else:
            missing = []
            if selected_datasets_count == 0:
                missing.append("datasets")
            if selected_contrasts_count == 0:
                missing.append("contrasts")
            if genes_count == 0:
                missing.append("genes")
            st.info(f"ℹ️ Please select {', '.join(missing)} to enable plot generation.")

    with tab1:
        st.header("Differential Expression Heatmap")
        st.markdown("**📊 Interactive heatmap showing log2 fold changes for selected genes across contrasts.** Hover over cells for details. Use the sidebar to adjust significance thresholds and filtering options.")

        # Display settings for heatmap
        with st.sidebar.expander("🎨 Display Settings", expanded=False):
            heatmap_font_size = st.slider("Font size", 8, 16, 12, key="heatmap_font")
            show_grid_lines = st.checkbox("Show grid lines", value=True, key="heatmap_grid")
            grid_opacity = st.slider("Grid opacity", 0.1, 1.0, 0.3, key="heatmap_grid_opacity")

        selected_contrasts = list(st.session_state.get('selected_contrasts', set()))

        if not gene_sel:
            st.info("Please select genes from the sidebar.")
        elif not selected_contrasts:
            st.info("Please select contrasts in the 'Selections' tab.")
        else:
                # Create and display the heatmap using isolated fragment
                @st.fragment
                def draw_heatmap(gene_selection, contrast_pairs, show_debug=False, show_adv=False):
                    with st.spinner("Generating heatmap..."):
                        try:
                            if show_adv and show_debug:
                                st.info("Debug mode: Showing detailed heatmap generation info")
                                import logging
                                import io

                                # Set up a string IO to capture log messages
                                log_stream = io.StringIO()
                                handler = logging.StreamHandler(log_stream)
                                handler.setLevel(logging.DEBUG)
                                formatter = logging.Formatter('%(levelname)s - %(message)s')
                                handler.setFormatter(formatter)

                                # Add the handler to the logger
                                logger = logging.getLogger("ResultsIntegration")
                                logger.setLevel(logging.DEBUG)
                                logger.addHandler(handler)

                            # Create heatmap with possibly modified parameters
                            # Use cached version if available for older Streamlit versions
                            if 'cached_figure_creation' in globals():
                                fig = cached_figure_creation("create_lfc_heatmap",
                                                            gene_selection,
                                                            contrast_pairs,
                                                            None)
                            else:
                                # Apply dynamic filtering if enabled
                                p_thresh = effective_pvalue_thresh if use_dynamic_filtering else None
                                lfc_thresh_val = effective_lfc_thresh if use_dynamic_filtering else None

                                fig = ri.create_lfc_heatmap(
                                    genes=gene_selection,
                                    contrasts=contrast_pairs,
                                    output_file=None,
                                    p_value_threshold=p_thresh,
                                    lfc_threshold=lfc_thresh_val,
                                    hide_empty_rows_cols=hide_empty_rows_cols,
                                    font_size=heatmap_font_size,
                                    show_grid_lines=show_grid_lines,
                                    grid_opacity=grid_opacity
                                )

                            # Display settings are now handled in the create_lfc_heatmap function

                            if show_adv and show_debug:
                                # Remove the handler to avoid duplicates
                                logger.removeHandler(handler)

                                # Display the log
                                st.expander("Debug Log", expanded=True).code(log_stream.getvalue())

                            if fig:
                                # Add notes about functionality
                                info_messages = ["💡 Hover over cells in the heatmap to see contrast descriptions and gene information."]
                                # if use_dynamic_filtering:
                                #    if show_advanced and 'use_separate_heatmap_filters' in locals() and use_separate_heatmap_filters:
                                #        info_messages.append(f"🎯 Only genes meeting heatmap-specific thresholds (p<{effective_pvalue_thresh:.3f}, #|LFC|>{effective_lfc_thresh:.1f}) are colored.")
                                #    else:
                                #        info_messages.append(f"🎯 Only genes meeting Auto-selected DEGs thresholds (p<{effective_pvalue_thresh:.3f}, |LFC|>{effective_lfc_thresh:.1f}) are colored.")
                                #if hide_empty_rows_cols:
                                #    info_messages.append("🧹 Genes/contrasts with no significant values have been removed.")

                                for msg in info_messages:
                                    st.info(msg)

                                # Display the plot
                                st.plotly_chart(fig, use_container_width=True)
                            else:
                                st.error("Could not generate heatmap. Please check your selections.")
                        except Exception as e:
                            st.error(f"Error generating heatmap: {str(e)}")

                # Call the fragment with just the input parameters
                draw_heatmap(gene_sel, selected_contrasts, debug_mode, show_advanced)

    with tab2:
        st.header("Gene Expression Plots")
        st.markdown("**🎻 Violin plots showing gene expression distributions across sample groups.** Each panel represents one gene, with samples grouped by experimental conditions.")

        # Display settings for expression plots
        with st.sidebar.expander("🎨 Display Settings", expanded=False):
            facet_font_size = st.slider("Facet title size", 8, 16, 10, key="violin_font")
            lock_y_axis = st.checkbox("Lock y-axis across genes", value=False, key="violin_lock_y")
            show_raw_points = st.checkbox("Show raw points", value=True, key="violin_points")
            legend_position = st.selectbox("Legend position", ["Bottom", "Right", "Top"], index=0, key="violin_legend")

        selected_datasets = list(st.session_state.get('selected_datasets', set()))

        if not gene_sel:
            st.info("Please select genes from the sidebar.")
        elif not selected_datasets:
            st.info("Please select datasets in the 'Selections' tab.")
        else:
                # Create and display the expression plots using isolated fragment
                @st.fragment
                def draw_expression_plots(gene_selection, dataset_selection, plot_style, hide_labels, page_num, total_pgs, facet_font_size, lock_y_axis, show_raw_points, legend_position):
                    with st.spinner("Generating expression plots..."):
                        try:
                            # Calculate gene slice for the current page
                            genes_per_page = 30
                            current_page = page_num
                            start_idx = (current_page - 1) * genes_per_page
                            end_idx = min(start_idx + genes_per_page, len(gene_selection))
                            current_genes = gene_selection[start_idx:end_idx]

                            fig2 = ri.create_expression_plots(
                                genes=current_genes,
                                analyses=dataset_selection,
                                plot_type="violin",
                                output_file=None,
                                hide_x_labels=hide_labels,
                                page_number=page_num,
                                facet_font_size=facet_font_size,
                                lock_y_axis=lock_y_axis,
                                show_raw_points=show_raw_points,
                                legend_position=legend_position,
                                show_grid_lines=True,
                                grid_opacity=0.3
                            )

                            if fig2:
                                # Add page navigation info if we have multiple pages
                                if total_pgs > 1:
                                    # Create pagination controls at the top of the plot for convenience
                                    cols = st.columns([2, 1, 1, 1, 2])
                                    with cols[1]:
                                        prev_disabled = page_num <= 1
                                        if st.button("◀ Previous", disabled=prev_disabled, key="prev_main"):
                                            st.session_state.page_num = max(1, page_num - 1)
                                            try:
                                                st.rerun()
                                            except AttributeError:
                                                # Fallback for older Streamlit versions
                                                try:
                                                    st.experimental_rerun()
                                                except AttributeError:
                                                    import streamlit.runtime.scriptrunner.magic as _m
                                                    _m._set_stop_thread(False)  # hidden API for very old versions
                                    with cols[2]:
                                        st.markdown(f"**Page {page_num}/{total_pgs}**")
                                    with cols[3]:
                                        next_disabled = page_num >= total_pgs
                                        if st.button("Next ▶", disabled=next_disabled, key="next_main"):
                                            st.session_state.page_num = min(total_pgs, page_num + 1)
                                            try:
                                                st.rerun()
                                            except AttributeError:
                                                # Fallback for older Streamlit versions
                                                try:
                                                    st.experimental_rerun()
                                                except AttributeError:
                                                    import streamlit.runtime.scriptrunner.magic as _m
                                                    _m._set_stop_thread(False)  # hidden API for very old versions

                                st.plotly_chart(fig2, use_container_width=True)
                            else:
                                st.error("Could not generate expression plots. Please check your selections.")
                        except Exception as e:
                            st.error(f"Error generating expression plots: {str(e)}")

                # Call the fragment with all input parameters including display settings
                draw_expression_plots(gene_sel, selected_datasets, "violin", hide_x_labels, st.session_state.page_num, total_pages, facet_font_size, lock_y_axis, show_raw_points, legend_position)

    with tab3:
        st.header("RNA-seq Analysis Plots")
        st.markdown("**📈 Quality control and differential expression plots from individual datasets.** View MDS plots, normalization diagnostics, volcano plots, and MA plots for detailed analysis.")

        # Display settings for analysis plots
        with st.sidebar.expander("🎨 Display Settings", expanded=False):
            plot_font_size = st.slider("Font size", 8, 16, 12, key="analysis_font")
            plot_grid = st.checkbox("Show grid lines", value=True, key="analysis_grid")

        # Select a dataset first - isolate this in a fragment to prevent full rerun
        @st.fragment
        def analysis_plots_tab():
            dataset_options = list(ri.cpm_data.keys()) if ri else []
            selected_dataset = st.selectbox("Select a dataset to view analysis plots:", dataset_options)

            if not selected_dataset:
                st.info("Please select a dataset to view RNA-seq analysis plots.")
                return

            # Now handle the plot display for the selected dataset
            # Find QC plots for this dataset
            base_path = os.path.join(results_dir, selected_dataset, "RNAseqAnalysis")

            # Check if the path exists
            if not os.path.exists(base_path):
                st.error(f"Could not find RNAseqAnalysis folder for {selected_dataset}")
            else:
                # General QC plots
                st.subheader("Quality Control and Normalization")

                qc_plot_files = {
                    "MDS Plot": {"file": os.path.join(base_path, "MDS.png"),
                                "description": "**Multi-dimensional scaling plot** showing sample relationships based on their gene expression profiles. This is analogous to a PCA plot.\n\n- **What to look for**: Samples from the same experimental group should cluster together. Outliers may indicate technical issues.\n\n- **Technical details**: Distances on the plot represent leading log2-fold changes between samples (the average of the largest absolute log2-fold changes)."},
                    "Filtering Density": {"file": os.path.join(base_path, "filtering_density.png"),
                                         "description": "**Density plot** showing distribution of gene expression levels (log-CPM) before and after filtering low-expression genes.\n\n- **What to look for**: After filtering (red curve), the left peak of very low-expressed genes should be reduced or eliminated.\n\n- **Technical details**: Low-expression genes are filtered out because they provide little statistical power for differential expression analysis and can increase the multiple testing burden."},
                    "Normalization Boxplots": {"file": os.path.join(base_path, "normalization_boxplots.png"),
                                    "description": "**Boxplots** showing expression distributions before (left) and after (right) normalization.\n\n- **What to look for**: After normalization, all samples should have similar median expression and comparable distributions.\n\n- **Technical details**: Normalization adjusts for technical differences in sequencing depth and composition between libraries to make samples comparable."},
                    "Mean-Variance Relationship (voom)": {"file": os.path.join(base_path, "voom_mean_variance.png"),
                                                 "description": "**Voom transformation plot** showing how the variance of gene expression depends on the mean expression level.\n\n- **What to look for**: A smooth decreasing trend where variance is higher for low-expressed genes and stabilizes for highly-expressed genes.\n\n- **Technical details**: The voom method transforms count data to log-CPM values and estimates the mean-variance relationship to assign appropriate weights for linear modeling."},
                    "SA Plot": {"file": os.path.join(base_path, "sa_plot.png"),
                              "description": "**Sigma vs Average plot** showing the standard deviation of normalized expression values against their average expression.\n\n- **What to look for**: A smooth trend without unusual patterns or outliers.\n\n- **Technical details**: This diagnostic plot from limma helps visualize how gene-wise variances change with expression level after fitting the linear model."}
                }

                # Display QC plots in a structured layout with tabs
                plot_tabs = st.tabs(["Quality Overview", "Expression Filtering", "Normalization", "Variance Modeling"])

                # Tab 1: Quality Overview (MDS plot)
                with plot_tabs[0]:
                    if os.path.exists(qc_plot_files["MDS Plot"]["file"]):
                        st.image(qc_plot_files["MDS Plot"]["file"])
                        with st.expander("What does this plot mean?", expanded=True):
                            st.markdown(qc_plot_files["MDS Plot"]["description"])
                    else:
                        st.info("MDS plot not available for this dataset.")

                # Tab 2: Expression Filtering
                with plot_tabs[1]:
                    if os.path.exists(qc_plot_files["Filtering Density"]["file"]):
                        st.image(qc_plot_files["Filtering Density"]["file"])
                        with st.expander("What does this plot mean?", expanded=True):
                            st.markdown(qc_plot_files["Filtering Density"]["description"])
                    else:
                        st.info("Filtering density plot not available for this dataset.")

                # Tab 3: Normalization
                with plot_tabs[2]:
                    if os.path.exists(qc_plot_files["Normalization Boxplots"]["file"]):
                        st.image(qc_plot_files["Normalization Boxplots"]["file"])
                        with st.expander("What does this plot mean?", expanded=True):
                            st.markdown(qc_plot_files["Normalization Boxplots"]["description"])
                    else:
                        st.info("Normalization boxplots not available for this dataset.")

                # Tab 4: Variance Modeling
                with plot_tabs[3]:
                    col1, col2 = st.columns(2)

                    with col1:
                        if os.path.exists(qc_plot_files["Mean-Variance Relationship (voom)"]["file"]):
                            st.subheader("Voom Transformation")
                            st.image(qc_plot_files["Mean-Variance Relationship (voom)"]["file"])
                            with st.expander("What does this plot mean?"):
                                st.markdown(qc_plot_files["Mean-Variance Relationship (voom)"]["description"])

                    with col2:
                        if os.path.exists(qc_plot_files["SA Plot"]["file"]):
                            st.subheader("Variance Trend")
                            st.image(qc_plot_files["SA Plot"]["file"])
                            with st.expander("What does this plot mean?"):
                                st.markdown(qc_plot_files["SA Plot"]["description"])

                # Contrast-specific plots
                st.markdown("---")
                st.subheader("Differential Expression Results")

                # Find all contrast directories
                contrast_dirs = [d for d in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, d)) and d != "logs"]

                if contrast_dirs:
                    # Find contrast descriptions if available
                    contrast_descriptions = {}
                    for c_id in contrast_dirs:
                        if hasattr(ri, "contrast_info") and c_id in ri.contrast_info:
                            contrast_descriptions[c_id] = ri.contrast_info[c_id].get('description', '')

                    # Create a list of contrast options with descriptions
                    contrast_options = []
                    for c_id in contrast_dirs:
                        label = f"{c_id}"
                        if c_id in contrast_descriptions and contrast_descriptions[c_id]:
                            # Truncate description if too long
                            desc = contrast_descriptions[c_id]
                            if len(desc) > 80:
                                desc = desc[:77] + "..."
                            label += f" - {desc}"
                        contrast_options.append(label)

                    # Create a mapping from display string to actual contrast ID
                    contrast_display_to_id = {}
                    for i, opt in enumerate(contrast_options):
                        contrast_display_to_id[opt] = contrast_dirs[i]

                    # Debug: Log available contrasts
                    if debug_mode:
                        st.write(f"Available contrasts: {', '.join(contrast_dirs)}")

                    selected_contrast_display = st.selectbox("Select a contrast:", contrast_options)
                    selected_contrast = contrast_display_to_id[selected_contrast_display]

                    if selected_contrast:
                        contrast_path = os.path.join(base_path, selected_contrast)

                        # Show the full contrast description if available
                        if selected_contrast in contrast_descriptions and contrast_descriptions[selected_contrast]:
                            st.markdown(f"**Description:** {contrast_descriptions[selected_contrast]}")

                        # Define plot information with improved descriptions
                        contrast_plot_files = {
                            "Heatmap of Top DEGs": {"file": os.path.join(contrast_path, "heatmap_top50.png"),
                                      "description": "**Hierarchical clustering heatmap** showing the top 50 differentially expressed genes for this contrast.\n\n- **What to look for**: Distinct expression patterns that separate the experimental groups. Groups of co-regulated genes (clusters) may represent functional modules.\n\n- **Technical details**: Genes are clustered based on their expression similarity. Color intensity represents expression level (red = high, blue = low)."},
                            "MA Plot": {"file": os.path.join(contrast_path, "ma_plot.png"),
                                      "description": "**MA plot** showing the relationship between average expression level (x-axis) and log2 fold change (y-axis).\n\n- **What to look for**: Symmetrical distribution around y=0 with significant genes (red points) distributed across expression levels. Asymmetry might indicate normalization issues.\n\n- **Technical details**: Each point represents a gene. M (y-axis) is the log-ratio of expression (log2 fold change), and A (x-axis) is the average expression. Red points are statistically significant DEGs."},
                            "Volcano Plot": {"file": os.path.join(contrast_path, "volcano_plot.png"),
                                          "description": "**Volcano plot** showing the relationship between statistical significance (-log10 p-value on y-axis) and biological significance (log2 fold change on x-axis).\n\n- **What to look for**: Genes in the upper left and right corners are both statistically significant and have large fold changes, making them the most interesting candidates.\n\n- **Technical details**: Each point represents a gene. Red points are statistically significant DEGs after multiple testing correction."}
                        }

                        # Create tabs for the contrast-specific plots
                        de_plot_tabs = st.tabs(["Volcano Plot", "MA Plot", "Heatmap"])

                        # Tab 1: Volcano Plot
                        with de_plot_tabs[0]:
                            if os.path.exists(contrast_plot_files["Volcano Plot"]["file"]):
                                st.image(contrast_plot_files["Volcano Plot"]["file"])
                                with st.expander("What does this plot mean?", expanded=True):
                                    st.markdown(contrast_plot_files["Volcano Plot"]["description"])
                            else:
                                st.info("Volcano plot not available for this contrast.")

                        # Tab 2: MA Plot
                        with de_plot_tabs[1]:
                            if os.path.exists(contrast_plot_files["MA Plot"]["file"]):
                                st.image(contrast_plot_files["MA Plot"]["file"])
                                with st.expander("What does this plot mean?", expanded=True):
                                    st.markdown(contrast_plot_files["MA Plot"]["description"])
                            else:
                                st.info("MA plot not available for this contrast.")

                        # Tab 3: Heatmap
                        with de_plot_tabs[2]:
                            if os.path.exists(contrast_plot_files["Heatmap of Top DEGs"]["file"]):
                                st.image(contrast_plot_files["Heatmap of Top DEGs"]["file"])
                                with st.expander("What does this plot mean?", expanded=True):
                                    st.markdown(contrast_plot_files["Heatmap of Top DEGs"]["description"])
                            else:
                                st.info("Heatmap not available for this contrast.")

                        # Add link to view DEG table
                        deg_file = os.path.join(contrast_path, "DEG.csv")
                        if os.path.exists(deg_file):
                            st.markdown("---")
                            with st.expander("View Differentially Expressed Genes Table", expanded=False):
                                try:
                                    deg_df = pd.read_csv(deg_file)

                                    # Create a copy for display, preserving original numeric types
                                    display_df = deg_df.copy()

                                    # Format p-values for display while keeping them numeric for sorting
                                    column_config = {}
                                    if 'adj.P.Val' in display_df.columns:
                                        # Ensure p-values are numeric for proper sorting
                                        try:
                                            display_df['adj.P.Val'] = pd.to_numeric(display_df['adj.P.Val'])
                                            column_config['adj.P.Val'] = st.column_config.NumberColumn(
                                                "adj.P.Val",
                                                format="%.2e",
                                                help="Adjusted p-value (corrected for multiple testing)"
                                            )
                                        except:
                                            st.warning("Could not convert adj.P.Val to numeric format for sorting.")
                                    if 'P.Value' in display_df.columns:
                                        # Ensure p-values are numeric for proper sorting
                                        try:
                                            display_df['P.Value'] = pd.to_numeric(display_df['P.Value'])
                                            column_config['P.Value'] = st.column_config.NumberColumn(
                                                "P.Value",
                                                format="%.2e",
                                                help="Unadjusted p-value"
                                            )
                                        except:
                                            st.warning("Could not convert P.Value to numeric format for sorting.")
                                    if 'logFC' in display_df.columns:
                                        column_config['logFC'] = st.column_config.NumberColumn(
                                            "logFC",
                                            format="%.2f",
                                            help="Log2 fold change"
                                        )
                                    if 'Gene' in display_df.columns:
                                        column_config['Gene'] = st.column_config.TextColumn(
                                            "Gene",
                                            width="auto",
                                            help="Gene symbol or identifier"
                                        )
                                    if 'Description' in display_df.columns:
                                        column_config['Description'] = st.column_config.TextColumn(
                                            "Description",
                                            width="medium",
                                            help="Contrast description"
                                        )

                                    # Show the dataframe with proper formatting
                                    st.dataframe(
                                        display_df,
                                        use_container_width=True,
                                        column_config=column_config
                                    )

                                    # Add download button for DEG file
                                    csv = pd.read_csv(deg_file).to_csv(index=False)  # Use original values for download
                                    st.download_button(
                                        label="Download DEG Table as CSV",
                                        data=csv,
                                        file_name=f"{selected_dataset}_{selected_contrast}_DEG.csv",
                                        mime="text/csv"
                                    )
                                except Exception as e:
                                    st.error(f"Error loading DEG file: {str(e)}")
                else:
                    st.info("No contrast-specific plots found for this dataset.")

        # Call the isolated fragment
        analysis_plots_tab()

    with tab4:
        st.header("Dataset Information")
        st.markdown("**📋 Browse and filter dataset metadata.** View study details, organism information, sample counts, and experimental descriptions. Use filters to find specific datasets of interest.")

        # Isolate datasets tab with fragment to prevent recomputation
        @st.fragment
        def datasets_tab():
            # Create a DataFrame with dataset information
            dataset_info = []
            for analysis_id, info in ri.analysis_info.items():
                # Build a dataset info dictionary
                dataset_dict = {
                    "Accession": info.get("accession", "Unknown"),
                    "Organism": info.get("organism", "Unknown"),
                    "Number of Samples": info.get("number_of_samples", 0),
                    "Number of Contrasts": info.get("number_of_contrasts", 0),
                    "Dataset ID": analysis_id  # Keep dataset ID but place it last
                }

                # Add dataset_info.txt content split into three columns if available
                if hasattr(ri, "dataset_info") and analysis_id in getattr(ri, "dataset_info", {}):
                    # Remove any "Title:" prefix from the title field
                    title = ri.dataset_info[analysis_id].get("title", "")
                    if isinstance(title, str) and title.startswith("Title:"):
                        title = title[6:].strip()
                    dataset_dict["Title"] = title

                    dataset_dict["Summary"] = ri.dataset_info[analysis_id].get("summary", "")
                    dataset_dict["Design"] = ri.dataset_info[analysis_id].get("design", "")

                dataset_info.append(dataset_dict)

            # Display the dataset information
            if dataset_info:
                df = pd.DataFrame(dataset_info)

                # Add filtering options
                st.subheader("Filter Datasets")
                col1, col2 = st.columns(2)

                with col1:
                    organism_filter = st.multiselect(
                        "Filter by Organism",
                        options=sorted(df["Organism"].unique()),
                        default=[]
                    )

                with col2:
                    search_filter = st.text_input("Search Datasets", "")

                # Apply filters
                filtered_df = df
                if organism_filter:
                    filtered_df = filtered_df[filtered_df["Organism"].isin(organism_filter)]

                if search_filter:
                    search_mask = filtered_df.apply(
                        lambda row: any(search_filter.lower() in str(val).lower() for val in row),
                        axis=1
                    )
                    filtered_df = filtered_df[search_mask]

                # Display the filtered dataset information
                if not filtered_df.empty:
                    # Add checkbox column for selection
                    display_df = filtered_df.copy()
                    display_df["✔"] = display_df["Dataset ID"].isin(st.session_state.get('selected_datasets', set()))

                    # Display dataset information with dataframe for interactivity
                    edited_df = st.data_editor(
                        display_df,
                        hide_index=True,
                        use_container_width=True,
                        column_config={
                            "✔": st.column_config.CheckboxColumn("Select", default=False),
                            "Dataset ID": st.column_config.TextColumn("Dataset ID", width="medium")
                        },
                        key="dataset_info_editor"
                    )

                    # Update selections based on checkboxes
                    if not edited_df.empty:
                        selected_from_info = set(edited_df.loc[edited_df["✔"], "Dataset ID"].tolist())
                        st.session_state['selected_datasets'] = selected_from_info

                    # Show dataset details if available (always show by default)
                    if any(col in filtered_df.columns for col in ["Title", "Summary", "Design"]):
                        for _, row in filtered_df.iterrows():
                            dataset_id = row.get("Dataset ID", row.get("Accession", "Unknown"))
                            accession = row.get("Accession", "")
                            with st.expander(f"Details for {dataset_id} {f'({accession})' if accession else ''}", expanded=True):
                                if "Title" in filtered_df.columns and pd.notna(row.get("Title")):
                                    st.subheader("Title")
                                    st.markdown(str(row["Title"]))

                                if "Summary" in filtered_df.columns and pd.notna(row.get("Summary")):
                                    st.subheader("Summary")
                                    st.markdown(str(row["Summary"]))

                                if "Design" in filtered_df.columns and pd.notna(row.get("Design")):
                                    st.subheader("Overall Design")
                                    st.markdown(str(row["Design"]))

                    # Add quick selection button
                    if st.button("Select all visible datasets", key="select_all_visible_datasets"):
                        visible_datasets = set(display_df["Dataset ID"].tolist())
                        st.session_state['selected_datasets'] = visible_datasets
                        # Reset page number when changing datasets
                        st.session_state.page_num = 1
                        st.success(f"Selected {len(visible_datasets)} datasets for analysis!")
                        st.info("Switch to the Heat-map or Expression tab to view updated visualizations.")
                else:
                    st.info("No datasets match the current filters.")
            else:
                st.info("No dataset information available.")

        # Call the isolated fragment
        datasets_tab()

    with tab5:
        st.header("Contrast Information")
        st.markdown("**🔍 Browse and filter contrast details.** View contrast descriptions, DEG counts, and filter by dataset or significance. Use this to understand what each comparison represents.")

        # Isolate contrasts tab with fragment to prevent recomputation
        @st.fragment
        def contrasts_tab():
            # Create a DataFrame with contrast information
            # Display contrast information
            contrast_info = []
            for aid, contrasts in ri.deg_data.items():
                for contrast_id in contrasts.keys():
                    # Get original contrast name and description
                    original_name = contrast_id
                    description = "No description available"

                    # Try to get information from contrast_info
                    if hasattr(ri, "contrast_info") and contrast_id in ri.contrast_info:
                        description = ri.contrast_info[contrast_id].get('description', "No description available")
                        # Use original name if available
                        if 'name' in ri.contrast_info[contrast_id]:
                            original_name = ri.contrast_info[contrast_id]['name']

                    # Count DEGs for this contrast
                    deg_count = 0
                    if aid in ri.deg_data and contrast_id in ri.deg_data[aid]:
                        df = ri.deg_data[aid][contrast_id]
                        # Use exact column names from DEG.csv file - prefer adjusted p-value
                        p_value_col = None
                        if 'adj.P.Val' in df.columns:
                            p_value_col = 'adj.P.Val'  # Adjusted p-value (preferred)
                        elif 'P.Value' in df.columns:
                            p_value_col = 'P.Value'  # Unadjusted p-value (fallback)

                        lfc_col = 'logFC' if 'logFC' in df.columns else None

                        if p_value_col and lfc_col:
                            deg_count = ((df[p_value_col] < pvalue_thresh) &
                                         (abs(df[lfc_col]) > lfc_thresh)).sum()

                    contrast_info.append({
                        "Dataset": aid,
                        "Accession": ri.analysis_info.get(aid, {}).get("accession", "Unknown"),
                        "Contrast": original_name,
                        "Original ID": contrast_id,
                        "Description": description,
                        "DEGs": deg_count
                    })

            # Display the contrast information
            if contrast_info:
                df = pd.DataFrame(contrast_info)
                # Format DEG count column as integer
                if 'DEGs' in df.columns:
                    df['DEGs'] = df['DEGs'].astype(int)

                # Add filtering options
                st.subheader("Filter Contrasts")
                col1, col2, col3 = st.columns(3)

                with col1:
                    dataset_column = "Dataset" if "Dataset" in df.columns else "Accession"
                    dataset_filter = st.multiselect(
                        "Filter by Dataset",
                        options=sorted(df[dataset_column].unique()),
                        default=[]
                    )

                with col2:
                    min_degs = st.number_input("Minimum DEGs", min_value=0, value=0)

                with col3:
                    search_filter = st.text_input("Search Contrasts", "")

                # Apply filters
                filtered_df = df
                if dataset_filter:
                    # Get the column to filter by
                    dataset_column = "Dataset" if "Dataset" in filtered_df.columns else "Accession"
                    filtered_df = filtered_df[filtered_df[dataset_column].isin(dataset_filter)]

                if min_degs > 0:
                    filtered_df = filtered_df[filtered_df["DEGs"] >= min_degs]

                if search_filter:
                    search_mask = filtered_df.apply(
                        lambda row: any(search_filter.lower() in str(val).lower() for val in row),
                        axis=1
                    )
                    filtered_df = filtered_df[search_mask]

                # Display the filtered contrast information
                if not filtered_df.empty:
                    # Sort by DEG count by default
                    filtered_df = filtered_df.sort_values("DEGs", ascending=False)

                    # Add checkbox column for selection
                    display_df = filtered_df.copy()
                    display_df["✔"] = display_df.apply(
                        lambda row: (row['Dataset'], row.get('Original ID', row['Contrast'])) in st.session_state.get('selected_contrasts', set()),
                        axis=1
                    )

                    # Display contrast information with dataframe for interactivity
                    edited_df = st.data_editor(
                        display_df,
                        hide_index=True,
                        use_container_width=True,
                        column_config={
                            "✔": st.column_config.CheckboxColumn("Select", default=False),
                            "Dataset": st.column_config.TextColumn("Dataset", width="medium"),
                            "Contrast": st.column_config.TextColumn("Contrast", width="medium"),
                            "Original ID": st.column_config.TextColumn("Original ID", width="medium")
                        },
                        key="contrast_info_editor"
                    )

                    # Update selections based on checkboxes
                    if not edited_df.empty:
                        selected_from_info = set()
                        for _, row in edited_df.iterrows():
                            if row["✔"]:
                                dataset = row['Dataset']
                                contrast = row.get('Original ID', row['Contrast'])
                                selected_from_info.add((dataset, contrast))
                        st.session_state['selected_contrasts'] = selected_from_info

                    # Add quick selection button
                    if st.button("Select all visible contrasts", key="select_all_visible_contrasts"):
                        visible_contrasts = set()
                        for _, row in display_df.iterrows():
                            dataset = row['Dataset']
                            contrast = row.get('Original ID', row['Contrast'])
                            visible_contrasts.add((dataset, contrast))
                        st.session_state['selected_contrasts'] = visible_contrasts
                        # Reset page number when changing contrasts
                        st.session_state.page_num = 1
                        st.success(f"Selected {len(visible_contrasts)} contrasts for analysis!")
                        st.info("Switch to the Heat-map tab to view updated visualizations.")
                else:
                    st.info("No contrasts match the current filters.")
            else:
                st.info("No contrast information available.")

        # Call the isolated fragment
        contrasts_tab()

    # ---------- 3. additional features -------------------------------
    st.sidebar.divider()
    st.sidebar.subheader("Export Options")
    export_format = st.sidebar.selectbox("Export format:", ["HTML", "CSV"])

    if st.sidebar.button("Export Current View"):
            if export_format == "HTML":
                try:
                    # Create a temporary directory for the export
                    import tempfile
                    import shutil

                    with tempfile.TemporaryDirectory() as tmpdir:
                        # Generate the interactive report
                        output_dir = ri.create_integrated_report(
                            top_frequent=top_frequent_genes,
                            top_unique=top_unique_genes,
                            plot_type="violin",
                            gene_list=gene_sel,
                            max_genes=100,
                            min_unique_per_contrast=min_unique,
                            p_value_threshold=pvalue_thresh,
                            lfc_threshold=lfc_thresh,
                            max_contrasts_for_unique=max_contrasts_unique,
                            hide_x_labels=hide_x_labels,
                            output_dir=tmpdir
                        )

                        # Create a zip file
                        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                        zip_filename = f"uorca_explorer_export_{timestamp}"
                        zip_path = os.path.join(tmpdir, f"{zip_filename}.zip")
                        shutil.make_archive(os.path.join(tmpdir, zip_filename), 'zip', output_dir)

                        # Provide download link
                        with open(zip_path, "rb") as f:
                            st.sidebar.download_button(
                                label="Download HTML Report",
                                data=f,
                                file_name="uorca_explorer_export.zip",
                                mime="application/zip"
                            )
                except Exception as e:
                    st.sidebar.error(f"Error exporting HTML: {str(e)}")
            else:
                try:
                    # Create CSV exports
                    import io

                    # Create a buffer for the zip file
                    buffer = io.BytesIO()
                    import zipfile

                    with zipfile.ZipFile(buffer, 'w') as zf:
                        # Export gene expression data
                        selected_datasets = st.session_state.get('selected_datasets', set())
                        for dataset_id in selected_datasets:
                            if dataset_id in ri.cpm_data:
                                df = ri.cpm_data[dataset_id]
                                if 'Gene' in df.columns:
                                    # Filter to selected genes if any
                                    if gene_sel:
                                        df = df[df['Gene'].isin(gene_sel)]
                                    csv_data = df.to_csv(index=False)
                                    zf.writestr(f"{dataset_id}_expression.csv", csv_data)

                        # Export DEG data
                        selected_contrasts = st.session_state.get('selected_contrasts', set())
                        for aid, cid in selected_contrasts:
                            if aid in ri.deg_data and cid in ri.deg_data[aid]:
                                df = ri.deg_data[aid][cid]
                                if 'Gene' in df.columns:
                                    # Filter to selected genes if any
                                    if gene_sel:
                                        df = df[df['Gene'].isin(gene_sel)]
                                    csv_data = df.to_csv(index=False)
                                    zf.writestr(f"{aid}_{cid}_DEG.csv", csv_data)

                    # Reset buffer position
                    buffer.seek(0)

                    # Provide download button
                    st.sidebar.download_button(
                        label="Download CSV Data",
                        data=buffer,
                        file_name="uorca_explorer_data.zip",
                        mime="application/zip"
                    )
                except Exception as e:
                    st.sidebar.error(f"Error exporting CSV: {str(e)}")

    # ---------- 4. housekeeping ------------------------------------
    st.sidebar.divider()
    st.sidebar.markdown(
        """
        *Powered by Streamlit & Plotly*

        **How to access this app remotely:**
        1. Run this app on your server: `streamlit run uorca_explorer.py`
        2. Create an SSH tunnel from your local machine:
           ```
           ssh -L 8000:127.0.0.1:8501 user@server
           ```
        3. Open `http://127.0.0.1:8000` in your browser
        """
    )
else:
    # Show help if no data is loaded
    st.info(
        """
        ## Welcome to UORCA Explorer

        This app allows you to interactively explore RNA-seq results from UORCA analyses.

        ### Getting Started
        1. Enter the path to your UORCA results directory in the sidebar
        2. The app will load your data and display interactive visualizations
        3. Use the sidebar controls to filter genes, datasets, and contrasts

        ### Troubleshooting
        - Make sure the path contains valid UORCA analysis results
        - Each analysis should have a directory structure with:
          - RNAseqAnalysis/ directory with CPM.csv file
          - metadata/ directory with contrasts.csv and edger_analysis_samples.csv files
        """
    )

if __name__ == "__main__":
    # This is just here for documentation
    pass
