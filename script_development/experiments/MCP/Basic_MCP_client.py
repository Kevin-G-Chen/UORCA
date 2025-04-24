import asyncio
import os

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client
from dotenv import load_dotenv
from pydantic_ai import Agent
from pydantic_ai.mcp import MCPServerStdio

load_dotenv()
OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")

server = MCPServerStdio(
    'uv', args=['run', 'Basic_MCP_server.py', 'server']
)

agent = Agent(
    'openai:gpt-4o-mini',
    system_prompt="You are free to decide which tools to use",
    mcp_servers=[server],          # ① expose server tools
)

async def main():
    async with agent.run_mcp_servers():   # ② starts and stops subprocess
        result = await agent.run("Please write a short poem about jewels")
    print(result.output)                  # ③ may or may not call `poet`

if __name__ == "__main__":
    asyncio.run(main())
