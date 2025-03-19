from dataclasses import dataclass
from typing import List
from pydantic_ai import Agent, RunContext
from rich.console import Console

console = Console()

# Define a minimal dependency container with a list of file paths.
@dataclass
class TestDependency:
    file_paths: List[str]

# Create an agent that uses TestDependency. The system prompt instructs the LLM
# to list all file paths present in the dependency.
agent = Agent(
    'openai:gpt-4o-mini',
    deps_type=TestDependency,
    system_prompt="""
You are a helpful assistant. When given a dependency that contains a list of file paths, please enumerate each file path in your answer.
"""
)

# Define a single tool (an LLM call) that simply returns the dependency's file paths.
@agent.tool
async def list_file_paths(ctx: RunContext[TestDependency]) -> str:
    # Read the file paths from the dependency and return them as a comma-separated list.
    return "I see the following file paths: " + ", ".join(ctx.deps.file_paths)

if __name__ == "__main__":
    # Create a dependency instance with some example file paths.
    test_dep = TestDependency(
        file_paths=[
            "/path/to/file1.txt",
            "/path/to/file2.txt",
            "/another/path/to/file3.log"
        ]
    )
    # Run the agent synchronously with the provided dependency.
    result = agent.run_sync("Please list the file paths.", deps=test_dep)
    console.print(result.data)
