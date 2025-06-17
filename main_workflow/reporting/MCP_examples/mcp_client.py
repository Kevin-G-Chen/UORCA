import asyncio, json, os

from pydantic_ai.mcp import MCPServerStdio
from dotenv import load_dotenv
from pathlib import Path
from pydantic_ai import Agent
from mcp import ClientSession, StdioServerParameters
from rich.console import Console

load_dotenv()
OPENAI_API_KEY = os.getenv('OPENAI_API_KEY')
console = Console()

async def main():
    console.rule("[bold]MCP agentic test")
    results_dir = "/data/tki_agpdev/kevin/phd/aim1/UORCA_results/2025-06-10_TCellFull"
    # 1️⃣  spin up the memory server as a subprocess --------------------------
    env = os.environ
    server = MCPServerStdio(
        command="uv",
        args=["run", "mcp_server.py", "server"],
        env=env,
        timeout = 30
    )

    # 3️⃣  create the agent, attach the MCP connection -----------------------
    agent = Agent(
        model="openai:gpt-4o-mini",
        model_settings={"temperature": 0.1},
        mcp_servers=[server],
        system_prompt=(
            "You are an assistant that helps analyze RNA-seq datasets from UORCA. "
            "You can describe datasets, list available datasets"
        ),
    )

    # 4️⃣  run a single turn --------------------------------------------------
    async with agent.run_mcp_servers():
        result = await agent.run(f"Tell me about the first dataset in {results_dir}")
    console.print(result.output)

if __name__ == "__main__":
    asyncio.run(main())
