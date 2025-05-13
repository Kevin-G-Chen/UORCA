#!/usr/bin/env python3
"""
Integrated GEO-to-SRA metadata extractor.

This script combines two workflows:
 1. Identify relevant GEO datasets for a biological research query (uses OpenAI + NCBI Entrez).
 2. Fetch linked SRA runinfo metadata for each GEO accession

The final output joins GEO dataset info with selected SRA run metadata fields.
"""
from __future__ import annotations
import argparse
import asyncio
import json
import os
import statistics
import time
import sys
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd
from Bio import Entrez
from dotenv import load_dotenv
from openai import OpenAI
from pydantic import BaseModel, ConfigDict, Field
from tqdm import tqdm
import xml.etree.ElementTree as ET

# ----------------------
# Environment & prompts
# ----------------------
load_dotenv()
Entrez.email = os.getenv("ENTREZ_EMAIL")
Entrez.api_key = os.getenv("ENTREZ_API_KEY")
openai_api_key = os.getenv("OPENAI_API_KEY")
if not openai_api_key:
    raise RuntimeError("OPENAI_API_KEY not set")
client = OpenAI(api_key=openai_api_key)
PROMPT_DIR = os.getenv("PROMPT_DIR", "./main_workflow/prompts/dataset_identification")

def load_prompt(fname: str) -> str:
    path = Path(PROMPT_DIR) / fname
    return path.read_text().strip()

# ----------------------
# Pydantic models
# ----------------------
class ExtractedTerms(BaseModel):
    model_config = ConfigDict(extra="forbid")
    extracted_terms: List[str]
    expanded_terms: List[str]

class Assessment(BaseModel):
    model_config = ConfigDict(extra="forbid")
    ID: str
    RelevanceScore: float
    Justification: str

class Assessments(BaseModel):
    model_config = ConfigDict(extra="forbid")
    assessments: List[Assessment]

TERMS_SCHEMA = ExtractedTerms.model_json_schema()
ASSESS_SCHEMA = Assessments.model_json_schema()

# ----------------------
# OpenAI JSON helper
# ----------------------
def call_openai_json(prompt: str, schema: Dict[str, Any], name: str) -> dict:
    resp = client.responses.create(
        model="gpt-4.1-mini",
        input=prompt,
        text={"format": {"type": "json_schema", "name": name, "schema": schema, "strict": True}},
    )
    return json.loads(resp.output_text)

# ----------------------
# Step 1 & 2: GEO search
# ----------------------
def extract_terms(query: str) -> ExtractedTerms:
    sys_prompt = load_prompt("extract_terms.txt")
    prompt = f"{sys_prompt}\nQuery: {query}"
    data = call_openai_json(prompt, TERMS_SCHEMA, "extracted_terms")
    return ExtractedTerms.model_validate(data)



def perform_search(term: str, retmax: int = 50) -> List[str]:
    # Apply specific search filter
    search = f"{term} AND (\"Expression profiling by high throughput sequencing\"[Filter])"
    with Entrez.esearch(db="gds", term=search, retmode="xml", retmax=retmax) as h:
        out = Entrez.read(h)
    return out.get("IdList", [])


def fetch_geo_summaries(ids: List[str]) -> pd.DataFrame:
    rows = []
    print("⬇️ Fetching GEO summaries …")
    for batch in tqdm([ids[i:i+10] for i in range(0, len(ids), 10)], desc="GEO summaries"):
        with Entrez.esummary(db="gds", id=','.join(batch), retmode="xml") as h:
            summaries = Entrez.read(h)
        for rec in summaries:
            rows.append({
                'ID': rec.get('Accession', ''),
                'Title': rec.get('title', ''),
                'Summary': rec.get('summary', ''),
                'Accession': rec.get('Accession', ''),
                'Species': rec.get('taxon', ''),
                'Date': rec.get('PDAT', ''),
            })
        time.sleep(0.3)
    return pd.DataFrame(rows)

# ----------------------
# Step 3: Relevance scoring
# ----------------------
async def assess_batch(df: pd.DataFrame, query: str, batch_size: int=20) -> List[Assessment]:
    sys_prompt = load_prompt("assess_relevance.txt")
    out: List[Assessment] = []
    for i in range(0, len(df), batch_size):
        sub = df.iloc[i:i+batch_size]
        prompt = f"{sys_prompt}\nResearch query: \"{query}\"\nDatasets:\n" + sub.to_json(orient='records')
        data = await asyncio.to_thread(call_openai_json, prompt, ASSESS_SCHEMA, "assessments")
        out += [Assessment.model_validate(a) for a in data['assessments']]
    return out

async def repeated_relevance(df: pd.DataFrame, query: str, repeats: int) -> pd.DataFrame:
    print("📊 Scoring relevance …")
    all_runs: List[List[Assessment]] = []
    for _ in range(repeats):
        all_runs.append(await assess_batch(df, query))
    coll: Dict[str, Dict[str, Any]] = {}
    for run in all_runs:
        for a in run:
            entry = coll.setdefault(a.ID, {'scores': [], 'justifications': []})
            entry['scores'].append(a.RelevanceScore)
            entry['justifications'].append(a.Justification)
    records: List[Dict[str, Any]] = []
    for id_, v in coll.items():
        rec: Dict[str, Any] = {'ID': id_, 'RelevanceScore': round(statistics.mean(v['scores']), 2)}
        for i, (score, just) in enumerate(zip(v['scores'], v['justifications'])):
            rec[f'Run{i+1}Score'] = score
            rec[f'Run{i+1}Justification'] = just
        records.append(rec)
    return pd.DataFrame(records)

# ----------------------
# SRA extraction utils
# ----------------------
def strip_ns(root):
    for e in root.iter():
        if isinstance(e.tag, str) and '}' in e.tag:
            e.tag = e.tag.split('}',1)[1]

def get_geo_uid(acc: str) -> str:
    with Entrez.esearch(db="gds", term=acc) as h:
        rec = Entrez.read(h)
    return rec['IdList'][0]

def get_sra_uid(geo_uid: str) -> str:
    with Entrez.elink(dbfrom="gds", db="sra", id=geo_uid) as h:
        rec = Entrez.read(h)
    return rec[0]['LinkSetDb'][0]['Link'][0]['Id']

def fetch_sra_xml(uid: str) -> str:
    with Entrez.efetch(db="sra", id=uid, rettype="xml", retmode="text") as h:
        return h.read()

def parse_bioproject(xml_str: str) -> str:
    root = ET.fromstring(xml_str)
    strip_ns(root)
    stud = root.find('.//STUDY')
    if stud is not None and 'accession' in stud.attrib:
        return stud.attrib['accession']
    proj = root.find('.//PROJECT')
    if proj is not None and 'accession' in proj.attrib:
        return proj.attrib['accession']
    raise ValueError('No BioProject accession found')

def fetch_runinfo(proj_acc: str) -> pd.DataFrame:
    with Entrez.esearch(db="sra", term=proj_acc) as h:
        rec = Entrez.read(h)
    uid = rec['IdList'][0]
    with Entrez.efetch(db="sra", id=uid, rettype="runinfo", retmode="text") as h:
        df = pd.read_csv(h, low_memory=False)
    return df

# ----------------------
# Main
# ----------------------
def main():
    p = argparse.ArgumentParser()
    p.add_argument('-q','--query', required=True)
    p.add_argument('-n','--num-queries', type=int, default=3)
    p.add_argument('--retmax', type=int, default=50)
    p.add_argument('--max-assess', type=int, default=None, help="Maximum number of datasets to assess relevance")
    p.add_argument('-o','--output', default='final_combined.csv')
    args = p.parse_args()

    # GEO identification
    print("🧠 Extracting terms …")
    terms = extract_terms(args.query)
    print("🔍 Searching GEO …")
    search_terms = set(terms.extracted_terms + terms.expanded_terms)
    geo_ids: List[str] = []
    for t in search_terms:
        print(f"Searching for: {t}")
        geo_ids += perform_search(t, args.retmax)
    geo_ids = list(dict.fromkeys(geo_ids))
    if not geo_ids:
        print('No GEO IDs found')
        sys.exit(0)

    geo_df = fetch_geo_summaries(geo_ids)

    # Optionally limit number of assessments
    if args.max_assess is not None:
        geo_df = geo_df.head(args.max_assess)

    # Relevance scoring
    rel_df = asyncio.run(repeated_relevance(geo_df, args.query, args.num_queries))
    geo_full = geo_df.merge(rel_df, on='ID', how='left')

    # SRA metadata extraction
    print("🔗 Fetching SRA metadata …")
    runs = []
    for acc in tqdm(geo_full['Accession'].unique(), desc='Processing SRA'):
        try:
            g_uid = get_geo_uid(acc)
            s_uid = get_sra_uid(g_uid)
            xml = fetch_sra_xml(s_uid)
            bp = parse_bioproject(xml)
            df_run = fetch_runinfo(bp)
            df_run.insert(0, 'GEO_Accession', acc)
            runs.append(df_run)
        except Exception as e:
            print(f'Warning, skipping {acc}: {e}')
    if not runs:
        print('No SRA runs retrieved')
        sys.exit(0)
    sra_df = pd.concat(runs, ignore_index=True)

    # Join and subset
    final = geo_full.merge(sra_df, left_on='Accession', right_on='GEO_Accession', how='inner')
    # Keep only original GEO-identification columns + chosen SRA fields
    keep_cols = list(geo_full.columns) + ['LibrarySource', 'LibraryLayout', 'LibraryStrategy']
    final = final[keep_cols]

    # Save
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    final.to_csv(args.output, index=False)
    print(f'Saved combined table to {args.output}')

if __name__=='__main__':
    main()
