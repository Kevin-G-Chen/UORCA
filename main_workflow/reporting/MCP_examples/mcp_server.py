from mcp.server.fastmcp import FastMCP

from pydantic_ai import Agent
import shutil, os, json

server = FastMCP("Utility-Tools")

@server.tool()
async def disk_usage(path: str) -> str:
    """Return `du -sh` for a folder."""
    out = shutil.disk_usage(path)
    return json.dumps({"total": out.total, "used": out.used, "free": out.free})

if __name__ == "__main__":
    server.run()
