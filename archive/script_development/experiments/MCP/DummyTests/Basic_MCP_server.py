from mcp.server.fastmcp import FastMCP
from dotenv import load_dotenv
from pydantic_ai import Agent
from pydantic import BaseModel
import os

load_dotenv()
OPENAI_API_KEY = os.getenv('OPENAI_API_KEY')

# A local agent the tools can delegate to (optional)
llm_agent = Agent(
    "openai:gpt-4o-mini",
    system_prompt="Always answer in rhyme"
)

server = FastMCP("Jewels-and-Text Demo")

# ── Tool 1: poem generator (unchanged) ────────────────────────────────────────
@server.tool()                                   # 2️⃣ signature → JSON schema
async def poet(theme: str) -> str:
    """Generate a short rhyming poem about *theme*."""
    r = await llm_agent.run(f"Write a four-line poem about {theme}")
    return r.output


# ── Tool 2: reverse any string ───────────────────────────────────────────────
class ReverseInput(BaseModel):
    text: str

@server.tool()
async def reverse_text(args: ReverseInput) -> str:
    """Return the input string reversed."""
    return args.text[::-1]


# ── Tool 3: word counter ─────────────────────────────────────────────────────
@server.tool(description="Count the words in a piece of text")
async def word_count(text: str) -> int:
    return len(text.split())


# ── Tool 4: gem price “estimator” (toy example with two fields) ─────────────
class Gem(BaseModel):
    name: str
    carat: float

@server.tool()
async def gem_value(gem: Gem) -> str:
    """Very rough value estimate of a gem."""
    est = gem.carat * 1200          # made-up math
    return f"A {gem.carat} ct {gem.name} might retail for about ${est:,.0f}."


# Start the server when run as a script
if __name__ == "__main__":
    server.run()
