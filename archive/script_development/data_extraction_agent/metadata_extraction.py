#!/usr/bin/env python3
"""
geo_to_sra_metadata_debug.py
↳  GEO Series ➜ GSM ➜ SRX ➜ SRR  (long-format CSV) with verbose prints
"""

import argparse, subprocess, re, pandas as pd, GEOparse as gp, sys
from textwrap import indent

# ---------- helpers ----------------------------------------------------------
def sh(cmd):
    return subprocess.run(cmd, shell=True, check=True,
                          stdout=subprocess.PIPE, text=True).stdout.strip()

re_srx  = re.compile(r"(SR[XP]\d+)")
re_srr  = re.compile(r"(SRR\d+)")

def srx_from_rel(rel):
    if not rel:
        return None
    for line in rel:                 # <-- NEW: scan the whole list
        hit = re_srx.search(line)
        if hit:
            return hit.group(1)
    return None

def srrs_from_srx(srx):
    cmd = ("esearch -db sra -query " + srx +
           " | efetch -format runinfo | cut -d',' -f1 | grep ^SRR")
    try:
        return [s for s in sh(cmd).splitlines() if re_srr.match(s)]
    except subprocess.CalledProcessError:
        return []

# ---------- CLI --------------------------------------------------------------
ap = argparse.ArgumentParser()
ap.add_argument("gse")
ap.add_argument("--out", default="meta_long.csv")
args = ap.parse_args()

print(f"⏬  Fetching {args.gse}")
gse  = gp.get_GEO(args.gse, destdir=".", silent=True)             # docs  [oai_citation_attribution:2‡Gist](https://gist.github.com/tubuliferous/7aea86b008a84e28b07379e1545f83a1?utm_source=chatgpt.com)
gsms = gse.gsms
print(f"   ↳ {len(gsms)} samples")

# ---------- wide sample table -------------------------------------------------
meta_wide = (pd.DataFrame.from_dict(
                {n: {k: (v[0] if v else None)
                     for k, v in g.metadata.items()}
                 for n, g in gsms.items()},
                orient="index")
             .reset_index()
             .rename(columns={"index": "GSM"}))                   # force GSM col

# ---------- long run table ----------------------------------------------------
rows = []
for gsm, g in gsms.items():
    srx = srx_from_rel(g.metadata.get("relation"))
    if not srx:
        continue
    for srr in srrs_from_srx(srx):                                # Entrez Direct  [oai_citation_attribution:3‡NCBI Trace Archive](https://trace.ncbi.nlm.nih.gov/Traces/study1/?go=help&utm_source=chatgpt.com)
        rows.append({"GSM": gsm, "SRX": srx, "SRR": srr})

long_df = pd.DataFrame(rows)

# ---------- diagnostics -------------------------------------------------------
def head(df, n=3):
    return indent(df.head(n).to_string(index=False), "   ")

print("\nlong_df columns:", long_df.columns.tolist(), "rows:", len(long_df))
print(head(long_df))
print("\nmeta_wide columns:", meta_wide.columns.tolist(), "rows:", len(meta_wide))
print(head(meta_wide[['GSM','title']]))  # print only two cols to stay narrow

# ---------- merge & save ------------------------------------------------------
out_df = long_df.merge(meta_wide, on="GSM", how="left", validate="many_to_one")
out_df.to_csv(args.out, index=False)
print(f"\n✅  Wrote {len(out_df)} rows to {args.out}")
