import os
import shutil
from typing import List

from dotenv import load_dotenv

from mcp.server.fastmcp import FastMCP
from pydantic import BaseModel
from pydantic_ai import Agent

server = FastMCP('PydanticAI Server')

@server.tool()
async def list_files(path: str = ".") -> List[str]:
    """Return a list of files and directories at *path* (defaults to CWD)."""
    # NB: no recursion or globbing – keep it simple for the demo.
    return os.listdir(path)


class DiskUsage(BaseModel):
    total: int
    used: int
    free: int


@server.tool()
async def disk_usage(path: str = ".") -> DiskUsage:
    """Return disk-usage statistics for the filesystem containing *path*."""
    total, used, free = shutil.disk_usage(path)
    return DiskUsage(total=total, used=used, free=free)


if __name__ == "__main__":
    server.run()
