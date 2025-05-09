import asyncio, argparse, os
from dotenv import load_dotenv
from pydantic_ai import Agent
from pydantic_ai.mcp import MCPServerStdio          # stdio transport
from metadata.models import Contrasts

load_dotenv()
OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")



metadata_server = MCPServerStdio(
    "uv", args=["run", "server.py", "server"]
)

agent = Agent(
    "openai:o4-mini",
    system_prompt=(
        "You are an RNA-seq metadata expert.  Use whatever tools are "
        "necessary to satisfy the user."
    ),
    mcp_servers=[metadata_server],
)

async def main(path: str):
    async with agent.run_mcp_servers():
        reply = await agent.run(f"Analyse the metadata file at {path}",
            output_type =Contrasts)
    print(reply.output)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata_path", required=True)
    asyncio.run(main(parser.parse_args().metadata_path))
