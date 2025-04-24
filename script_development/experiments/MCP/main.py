from MetadataMCP_agent import metadata_agent, MetadataContext
from pydantic import BaseModel
from typing import List

# Define result type
class Contrasts(BaseModel):
    contrasts: List[dict]
    summary: str = None

async def run_analysis(metadata_path):
    # Create context
    analysis_data = MetadataContext(metadata_path=metadata_path)

    # Use context manager to start MCP servers
    async with metadata_agent.run_mcp_servers():
        # Run initial metadata processing
        initial = await metadata_agent.run(
            "Analyze this metadata file and suggest contrasts for differential expression",
            deps=analysis_data,
            result_type=Contrasts
        )
        return initial

if __name__ == "__main__":
    import asyncio
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata_path", required=True)
    args = parser.parse_args()

    result = asyncio.run(run_analysis(args.metadata_path))
    print(result.data)
