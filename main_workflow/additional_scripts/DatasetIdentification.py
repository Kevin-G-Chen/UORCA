#!/usr/bin/env python3
"""
Integrated GEO-to-SRA metadata extractor with relevance filtering and validity flag.

This script combines two workflows:
 1. Identify relevant GEO datasets for a biological research query (uses OpenAI + NCBI Entrez).
 2. Fetch linked SRA runinfo metadata for each GEO accession above a relevance threshold.

The final output joins GEO dataset info with selected SRA run metadata fields, includes all GEO datasets,
 and marks runs as valid based on library properties.
"""
from __future__ import annotations
import argparse
import asyncio
import json
import os
import statistics
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd
from Bio import Entrez
from dotenv import load_dotenv
from openai import OpenAI
from pydantic import BaseModel, ConfigDict
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
    search = f"{term} AND (\"Expression profiling by high throughput sequencing\"[Filter])"
    with Entrez.esearch(db="gds", term=search, retmode="xml", retmax=retmax) as h:
        out = Entrez.read(h)
    return out.get("IdList", [])


def fetch_geo_summaries(ids: List[str]) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    print("⬇️ Fetching GEO summaries …")
    for batch in tqdm([ids[i:i+10] for i in range(0, len(ids), 10)], desc="GEO summaries"):
        with Entrez.esummary(db="gds", id=','.join(batch), retmode="xml") as h:
            summaries = Entrez.read(h)
        for rec in summaries:
            # Safely extract PubMed IDs if available
            pmids = []
            try:
                if 'PubMedIds' in rec and rec['PubMedIds']:
                    pmids = [int(pmid) for pmid in rec['PubMedIds']]
            except (ValueError, TypeError) as e:
                logger.warning(f"Error extracting PubMed IDs for {rec.get('Accession', '')}: {e}")

            # Build row with optional PubMed data
            row_data = {
                'ID': rec.get('Accession', ''),
                'Title': rec.get('title', ''),
                'Summary': rec.get('summary', ''),
                'Accession': rec.get('Accession', ''),
                'Species': rec.get('taxon', ''),
                'Date': rec.get('PDAT', '')
            }

            # Add PubMed information when available
            if pmids:
                row_data['PrimaryPubMedID'] = pmids[0]
                row_data['AllPubMedIDs'] = ','.join(str(pid) for pid in pmids)
            else:
                row_data['PrimaryPubMedID'] = None
                row_data['AllPubMedIDs'] = None

            rows.append(row_data)
        time.sleep(0.3)
    return pd.DataFrame(rows)

# ----------------------
# Step 3: Relevance scoring
# ----------------------
async def assess_subbatch(sub_df: pd.DataFrame, query: str, schema: Dict[str, Any], name: str, rep: int, batch_idx: int, sem: asyncio.Semaphore) -> List[Assessment]:
    async with sem:
        print(f"⏳ Relevance rep {rep+1}, batch {batch_idx+1}/{(len(sub_df.index) // args.batch_size) + 1}")
        sys_prompt = load_prompt("assess_relevance.txt")
        prompt = f"{sys_prompt}\nResearch query: \"{query}\"\nDatasets to assess (batch {batch_idx+1}):\n" + sub_df.to_json(orient='records')
        data = await asyncio.to_thread(call_openai_json, prompt, schema, name)
        return [Assessment.model_validate(a) for a in data['assessments']]

async def repeated_relevance(df: pd.DataFrame, query: str, repeats: int, batch_size: int, concurrency: int) -> pd.DataFrame:
    print(f"📊 Starting relevance scoring: {repeats} repetitions, batch size {batch_size}, concurrency {concurrency}")
    sem = asyncio.Semaphore(concurrency)
    tasks = []
    # create sub-batches
    batches = [df.iloc[i:i+batch_size] for i in range(0, len(df), batch_size)]
    for rep in range(repeats):
        for idx, sub in enumerate(batches):
            tasks.append(assess_subbatch(sub, query, ASSESS_SCHEMA, "assessments", rep, idx, sem))
    all_results = await asyncio.gather(*tasks)
    # flatten and aggregate
    coll: Dict[str, Dict[str, Any]] = {}
    for result in all_results:
        for a in result:
            entry = coll.setdefault(a.ID, {'scores': [], 'justifications': []})
            entry['scores'].append(a.RelevanceScore)
            entry['justifications'].append(a.Justification)
    records: List[Dict[str, Any]] = []
    for id_, v in coll.items():
        rec = {'ID': id_, 'RelevanceScore': round(statistics.mean(v['scores']), 2)}
        for i, (score, just) in enumerate(zip(v['scores'], v['justifications']), 1):
            rec[f'Run{i}Score'] = score
            rec[f'Run{i}Justification'] = just
        records.append(rec)
    return pd.DataFrame(records)

# ----------------------
# SRA extraction utils
# ----------------------
def strip_ns(root: ET.Element):
    for e in root.iter():
        if isinstance(e.tag, str) and '}' in e.tag:
            e.tag = e.tag.split('}', 1)[1]


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
    parser = argparse.ArgumentParser()
    parser.add_argument('-q', '--query', required=True)
    parser.add_argument('-n', '--num-queries', type=int, default=3,
                        help='Number of relevance scoring repeats')
    parser.add_argument('--retmax', type=int, default=50,
                        help='Max GEO results per term')
    parser.add_argument('--max-assess', type=int, default=None,
                        help='Max datasets to assess relevance')
    parser.add_argument('-t', '--relevance-threshold', type=float, default=7.0,
                        help='Minimum relevance score to fetch SRA metadata')
    parser.add_argument('--batch-size', type=int, default=20,
                        help='Number of datasets per relevance batch')
    parser.add_argument('--concurrency', type=int, default=1,
                        help='Max parallel relevance API calls')
    parser.add_argument('-o', '--output', default='final_combined.csv',
                        help='Output CSV path')
    parser.add_argument('--generate-multi-csv', action='store_true',
                       help='Generate an additional CSV file formatted for multi-dataset analysis')

    global args
    args = parser.parse_args()

    print("🧠 Extracting terms …")
    terms = extract_terms(args.query)
    search_terms = set(terms.extracted_terms + terms.expanded_terms)

    print("🔍 Searching GEO …")
    geo_ids: List[str] = []
    for term in search_terms:
        print(f"Searching for: {term}")
        geo_ids.extend(perform_search(term, args.retmax))
    geo_ids = list(dict.fromkeys(geo_ids))
    if not geo_ids:
        print('No GEO IDs found')
        sys.exit(0)

    geo_df = fetch_geo_summaries(geo_ids)
    if args.max_assess is not None:
        geo_df = geo_df.head(args.max_assess)

    rel_df = asyncio.run(repeated_relevance(geo_df, args.query, args.num_queries, args.batch_size, args.concurrency))
    geo_full = geo_df.merge(rel_df, on='ID', how='left')

    # Filter for SRA fetch
    to_fetch = geo_full[geo_full['RelevanceScore'] >= args.relevance_threshold]
    print(f"🔗 Fetching SRA metadata for {len(to_fetch)} datasets (threshold >= {args.relevance_threshold}) …")

    runs: List[pd.DataFrame] = []
    for acc in tqdm(to_fetch['Accession'].unique(), desc='Processing SRA'):
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

    if runs:
        sra_df = pd.concat(runs, ignore_index=True)
    else:
        sra_df = pd.DataFrame(columns=['GEO_Accession'])

    # Merge all GEO datasets with fetched SRA info
    merged = geo_full.merge(
        sra_df,
        left_on='Accession',
        right_on='GEO_Accession',
        how='left'
    )

    # Select and flag validity
    cols_to_keep = list(geo_full.columns) + ['LibrarySource', 'LibraryLayout', 'LibraryStrategy']
    final = merged[cols_to_keep].copy()
    final['Valid'] = final.apply(
        lambda row: 'Yes' if (row.get('LibraryLayout') == 'PAIRED' \
                              and row.get('LibrarySource') == 'TRANSCRIPTOMIC') else 'No',
        axis=1
    )

    # For convenience, generate URLs to PubMed and GEO
    def create_pubmed_url(pmid):
        if pd.notna(pmid):
            return f"https://pubmed.ncbi.nlm.nih.gov/{int(pmid)}/"
        return None

    def create_geo_url(accession):
        if pd.notna(accession):
            return f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}"
        return None

    final['PubMedURL'] = final['PrimaryPubMedID'].apply(create_pubmed_url)
    final['GEOURL'] = final['Accession'].apply(create_geo_url)

    final = final.drop_duplicates(subset=['Accession'])


    # Save
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    final.to_csv(args.output, index=False)
    print(f'Saved combined table to {args.output}')

    if args.generate_multi_csv:
        multi_df = final[['Accession', 'Species', 'PrimaryPubMedID', 'RelevanceScore', 'Valid']].copy()
        multi_df = multi_df[multi_df['Valid'] == 'Yes']
        multi_df = multi_df.sort_values('RelevanceScore', ascending=False)
        multi_df = multi_df.rename(columns={'Species': 'organism'})
        multi_csv_path = os.path.join(os.path.dirname(args.output), 'multi_dataset_input.csv')
        multi_df.to_csv(multi_csv_path, index=False)
        print(f"Generated multi-dataset input CSV at {multi_csv_path}")

if __name__ == '__main__':
    main()
