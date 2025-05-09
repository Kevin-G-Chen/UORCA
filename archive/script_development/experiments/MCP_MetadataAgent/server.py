from mcp.server.fastmcp import FastMCP                 # 🛠 server toolkit
from pydantic_ai import Agent, RunContext              # 🔗 local LLM agent
from pydantic import BaseModel
from dotenv import load_dotenv
import pandas as pd, json, os, re
from unidecode import unidecode
from Bio import Entrez

# --- shared models ---------------------------------------------------------
from metadata.models import RNAseqData, Contrasts, ContrastFormat

load_dotenv()
OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")

metadata_system_prompt_path = "../../prompt_development/metadata_processing_systemprompt.txt"

try:
    system_prompt_text = metadata_system_prompt_path.read_text(encoding="utf-8")
except Exception:
    system_prompt_text = (
        "You are provided with RNA-seq metadata.  Identify biologically relevant "
        "columns, build a grouping factor, list unique groups and propose "
        "contrasts."
    )

# llm_agent = Agent("openai:o4-mini", system_prompt=system_prompt_text,)

server = FastMCP("RNA-seq Metadata MCP")

# ────────────────────────────── utility ────────────────────────────────────
def clean_string(s: str) -> str:
    if pd.isna(s):
        return "NA"
    s = unidecode(str(s).strip()).replace(" ", "_")
    return re.sub(r"[^\w]", "", s)

# ────────────────────────────── Tool 1 ─────────────────────────────────────
class MetadataFile(BaseModel):
    metadata_path: str

@server.tool(description="Process and clean an RNA-seq metadata file")
async def process_metadata(req: MetadataFile) -> dict:
    """Loads, cleans and summarises a metadata file."""
    p = req.metadata_path
    df = (pd.read_csv(p, sep="\t" if p.endswith((".tsv", ".txt")) else ",")
          .drop_duplicates())
    df = df.loc[:, df.nunique() < df.shape[0]]          # drop all-unique cols
    df.rename(columns={c: clean_string(c) for c in df.columns}, inplace=True)
    for col in df.columns:
        df[col] = df[col].apply(clean_string)
    stats = {c: {"unique_count": df[c].nunique()} for c in df.columns}
    df.to_pickle(p + ".clean.pkl")   # ← persist for later tool calls
    return {"shape": df.shape, "column_stats": stats}

# ────────────────────────────── Tool 2 ─────────────────────────────────────
class MergeArgs(BaseModel):
    metadata_path: str
    columns: list[str]

@server.tool(description="Merge columns into a single analysis factor")
async def merge_analysis_columns(args: MergeArgs) -> dict:
    df = pd.read_pickle(args.metadata_path + ".clean.pkl")
    if len(args.columns) == 1:
        merged = args.columns[0]
    else:
        merged = "merged_analysis_group"
        df[merged] = df[args.columns].astype(str).agg("_".join, axis=1)
    df.to_pickle(args.metadata_path + ".clean.pkl")
    return {"merged_column": merged, "preview": df[merged].unique()[:5].tolist()}

# ────────────────────────────── Tool 3 ─────────────────────────────────────
class UniqueReq(BaseModel):
    metadata_path: str
    analysis_column: str

@server.tool(description="Extract unique groups from the analysis column")
async def extract_unique_values(req: UniqueReq) -> dict:
    df = pd.read_pickle(req.metadata_path + ".clean.pkl")
    uni = sorted(df[req.analysis_column].unique().tolist())
    return {"unique_values": uni, "count": len(uni)}

# ────────────────────────────── Tool 4 ─────────────────────────────────────
@server.tool(description="Fetch a GEO summary for an accession found in the metadata")
async def fetch_geo_summary(accession: str) -> dict:
    Entrez.email = "you@example.com"
    search = Entrez.esearch(db="gds", term=f"{accession}[ACCN]", retmode="xml")
    uids = Entrez.read(search)["IdList"]; search.close()
    if not uids:
        return {"success": False, "message": "Accession not found"}
    summ = Entrez.read(Entrez.esummary(db="gds", id=uids[0], retmode="xml"))[0]
    return {"success": True, "summary": summ.get("summary", "No summary")}

# ────────────────────────────── Tool 5 ─────────────────────────────────────
class FullRun(BaseModel):
    metadata_path: str

#@server.tool(description="Run the complete five-step metadata workflow")
#async def analyse_metadata(req: FullRun) -> dict:
#    deps = RNAseqData(metadata_path=req.metadata_path)
#    prompt = (
#        "Please:\n"
#        "1. Process the metadata\n2. Pick relevant columns\n"
#        "3. Build a final grouping column\n4. List its unique values\n"
#        "5. Propose contrasts")
#    result = llm_agent.run_sync(prompt, deps=deps, output_type=Contrasts)
#    return json.loads(result.json())

# ───────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    server.run()
