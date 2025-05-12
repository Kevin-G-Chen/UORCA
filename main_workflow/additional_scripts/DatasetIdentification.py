#!/usr/bin/env python3
from __future__ import annotations

import argparse
import asyncio
import json
import os
import time
from collections import Counter
from typing import List, Dict, Any

import pandas as pd
from Bio import Entrez
from dotenv import load_dotenv
from openai import OpenAI
from pydantic import BaseModel, Field, ConfigDict
from tqdm import tqdm
from pathlib import Path

# -----------------------------------------------------------------------------#
# Environment
# -----------------------------------------------------------------------------#
load_dotenv()
Entrez.email = os.getenv("ENTREZ_EMAIL")
Entrez.api_key = os.getenv("ENTREZ_API_KEY")

openai_api_key = os.getenv("OPENAI_API_KEY")
if not openai_api_key:
    raise RuntimeError("OPENAI_API_KEY not set")

client = OpenAI(api_key=openai_api_key)

# -----------------------------------------------------------------------------#
# Pydantic models → JSON Schema
# -----------------------------------------------------------------------------#


class ExtractedTerms(BaseModel):
    model_config = ConfigDict(extra="forbid")
    extracted_terms: List[str] = Field(
        description="Biologically specific term(s) pulled directly from the query"
    )
    expanded_terms: List[str] = Field(
        description="Related or synonymous terms helpful for search expansion"
    )


class Assessment(BaseModel):
    model_config = ConfigDict(extra="forbid")
    ID: str = Field(description="Dataset identifier (GDS/GSE/GSM)")
    RelevanceScore: int = Field(description="0 = not relevant, 10 = highly relevant"
    )
    Justification: str = Field(
        description="Reason for the score; ≤ 30 words"
    )


class Assessments(BaseModel):
    model_config = ConfigDict(extra="forbid")
    assessments: List[Assessment]


TERMS_SCHEMA = ExtractedTerms.model_json_schema()
ASSESSMENTS_SCHEMA = Assessments.model_json_schema()

# -----------------------------------------------------------------------------#
# OpenAI helpers
# -----------------------------------------------------------------------------#


def call_openai_json(prompt: str, schema: Dict[str, Any], name: str) -> dict:
    response = client.responses.create(
        model="gpt-4o-mini",
        input=prompt,
        text={
            "format": {
                "type": "json_schema",
                "name": name,
                "schema": schema,
                "strict": True,
            }
        },
    )
    return json.loads(response.output_text)


# -----------------------------------------------------------------------------#
# Core steps
# -----------------------------------------------------------------------------#


def extract_terms(user_query: str) -> ExtractedTerms:
    system_prompt = (
        "You are an expert in biological text mining. "
        "Extract biologically specific terms from the query and generate useful synonyms."
    )
    prompt = f"{system_prompt}\nQuery: {user_query}"
    data = call_openai_json(prompt, TERMS_SCHEMA, "extracted_terms")
    return ExtractedTerms.model_validate(data)


def perform_search(term: str, retmax: int = 50) -> Dict[str, Any]:
    with Entrez.esearch(db="gds", term=term, retmode="xml", retmax=retmax) as handle:
        return Entrez.read(handle)


def extract_geo_info_batch(geo_ids: List[str]) -> List[Dict[str, Any]]:
    ids_str = ",".join(geo_ids)
    with Entrez.esummary(db="gds", id=ids_str, retmode="xml") as handle:
        summaries = Entrez.read(handle)

    data = []
    for geo_id, s in zip(geo_ids, summaries):
        if isinstance(s, dict):
            data.append(
                {
                    "ID": geo_id,
                    "Title": s.get("title", ""),
                    "Summary": s.get("summary", ""),
                    "Accession": s.get("Accession", ""),
                    "Species": s.get("taxon", ""),
                    "Date": s.get("PDAT", ""),
                }
            )
        else:
            data.append(
                {
                    "ID": geo_id,
                    "Title": "Error",
                    "Summary": "Unable to fetch",
                    "Accession": "Error",
                    "Species": "Error",
                    "Date": "Error",
                }
            )
    return data


async def assess_relevance_async(
    df: pd.DataFrame, query: str, batch_size: int = 20
) -> List[Assessment]:
    async def _single_batch(batch_df: pd.DataFrame) -> List[Assessment]:
        dataset_json = batch_df.to_json(orient="records")
        prompt = (
            "You are a highly knowledgeable biologist.\n"
            f"Research query: \"{query}\"\n\n"
            "Assess the following datasets (as JSON list) and output a JSON object "
            'with key "assessments" that satisfies the provided schema, one entry '
            "per dataset."
        )
        full_input = f"{prompt}\n\nDatasets:\n{dataset_json}"
        data = await asyncio.to_thread(
            call_openai_json, full_input, ASSESSMENTS_SCHEMA, "assessments"
        )
        return [Assessment.model_validate(a) for a in data["assessments"]]

    tasks = []
    for i in range(0, len(df), batch_size):
        tasks.append(_single_batch(df.iloc[i : i + batch_size]))

    results_nested = await asyncio.gather(*tasks)
    # Flatten
    return [ass for sub in results_nested for ass in sub]


async def repeated_relevance(
    df: pd.DataFrame, query: str, repeats: int = 3, batch_size: int = 20
) -> List[Dict[str, Any]]:
    all_runs = []
    for _ in range(repeats):
        all_runs.append(await assess_relevance_async(df, query, batch_size=batch_size))

    collated: Dict[str, Dict[str, Any]] = {}
    for run in all_runs:
        for ass in run:
            d = collated.setdefault(ass.ID, {"scores": [], "justifications": []})
            d["scores"].append(ass.RelevanceScore)
            d["justifications"].append(ass.Justification)

    majority = []
    for id_, d in collated.items():
        score_counts = Counter(d["scores"])
        majority_score, _ = score_counts.most_common(1)[0]
        record = {"ID": id_, "RelevanceScore": majority_score}
        for i, (score, just) in enumerate(zip(d["scores"], d["justifications"])):
            record[f"Run{i+1}Score"] = score
            record[f"Run{i+1}Justification"] = just
        majority.append(record)
    return majority


# -----------------------------------------------------------------------------#
# CLI entry‑point
# -----------------------------------------------------------------------------#


def main() -> None:
    p = argparse.ArgumentParser(
        description="Identify relevant GEO datasets for a biological research query."
    )
    p.add_argument("-q", "--query", required=True, help="Research query")
    p.add_argument("--num-queries", "-n", type=int, default=3, help="Relevance repeats")
    p.add_argument("--retmax", type=int, default=50, help="Max GEO records per term")
    p.add_argument("-o", "--output", default="results.csv", help="Output CSV")
    args = p.parse_args()

    # 1. Extract biological terms
    print("🧠 Extracting terms …")
    terms = extract_terms(args.query)
    search_terms = sorted(set(terms.extracted_terms + terms.expanded_terms))
    print("Search terms:", ", ".join(search_terms))

    # 2. GEO search
    print("🔍 Searching GEO …")
    geo_ids: List[str] = []
    for term in tqdm(search_terms, desc="Entrez search"):
        out = perform_search(term, retmax=args.retmax)
        geo_ids.extend(out.get("IdList", []))

    if not geo_ids:
        print("No GEO datasets found.")
        return

    # Deduplicate preserving order
    geo_ids = list(dict.fromkeys(geo_ids))

    # 3. Fetch summaries
    print("⬇️ Fetching summaries …")
    rows = []
    for i in tqdm(range(0, len(geo_ids), 10), desc="Summaries"):
        rows.extend(extract_geo_info_batch(geo_ids[i : i + 10]))
        time.sleep(0.3)  # NCBI etiquette

    df = pd.DataFrame(rows)

    # 4. Relevance scoring
    print("📊 Scoring relevance …")
    relevance_records = asyncio.run(
        repeated_relevance(df, args.query, repeats=args.num_queries)
    )
    rel_df = pd.DataFrame(relevance_records)

    # 5. Merge & save
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    final_df = df.merge(rel_df, on="ID", how="left")
    final_df.to_csv(args.output, index=False)
    print(f"✅ Saved {len(final_df)} rows to {args.output}")


if __name__ == "__main__":
    main()
