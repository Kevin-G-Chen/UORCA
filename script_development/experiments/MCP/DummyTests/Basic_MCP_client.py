# client_autonomous.py
import asyncio, os
from dotenv import load_dotenv
from pydantic_ai import Agent
from pydantic_ai.mcp import MCPServerStdio          # runs the server as a subprocess

load_dotenv()
OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")

# Point to the local server script we just wrote
jewel_server = MCPServerStdio(
    "uv",
    args=["run", "Basic_MCP_server.py", "server"]
)

agent = Agent(
    "openai:gpt-4o-mini",
    system_prompt=(
        "You are free to decide which tools to use in order to satisfy the user."
    ),
    mcp_servers=[jewel_server]                      # <- manifest injected here
)

async def main():
    async with agent.run_mcp_servers():             # spawns + cleans up subprocess
        reply = await agent.run(
            "I have 23 words.  \n"
            "1. Tell me how many words that is in total.\n"
            "2. Then write a rhyming poem about jewels.\n"
            "3. Finally, reverse the word 'brilliance'."
        )
    print(reply.output)                             # agent may chain several tools

if __name__ == "__main__":
    asyncio.run(main())
