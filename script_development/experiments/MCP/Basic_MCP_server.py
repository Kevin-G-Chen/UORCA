from mcp.server.fastmcp import FastMCP
from dotenv import load_dotenv
from pydantic_ai import Agent
import os

load_dotenv()
OPENAI_API_KEY = os.getenv('OPENAI_API_KEY')

server = FastMCP('PydanticAI Server')
server_agent = Agent(
    'openai:gpt-4o-mini', system_prompt='always reply in rhyme'
)


@server.tool()
async def poet(theme: str) -> str:
    """Poem generator"""
    r = await server_agent.run(f'write a poem about {theme}')
    return r.output

@server.tool()
async def chat(user_message: str) -> any:
    # Let the agent decide which of its tools to invoke
    result = await server_agent.run(user_message)
    return result

if __name__ == '__main__':
    server.run()
