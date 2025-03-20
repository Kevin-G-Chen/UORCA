from dataclasses import dataclass
from typing import List, Union
import os
from pydantic_ai import Agent, RunContext
from rich.console import Console

console = Console()

# Define a minimal dependency container with a test directory.
@dataclass
class TestDependency:
    test_dir: str

# Create an agent that uses TestDependency.
agent = Agent(
    'openai:gpt-4o-mini',
    deps_type=TestDependency,
    system_prompt="""
You are a helpful assistant. When given a dependency that contains a test directory path,
please use the provided tool to list all matching file paths.
"""
)

# Replicate the custom find_files function
@agent.tool
async def find_files(ctx: RunContext[TestDependency], directory: str, suffix: Union[str, List[str]]) -> List[str]:
    """
    Recursively search for and return a sorted list of files within the provided directory
    that have the given suffix (or any in the suffix list).
    """
    try:
        if not hasattr(ctx.deps, '_logged_context'):
            console.log(f"[bold blue]Initial Context.deps details:[/]\n{vars(ctx.deps)}")
            setattr(ctx.deps, '_logged_context', True)
        console.log(f"[bold blue]Tool Called:[/] find_files with directory: {directory}, suffix: {suffix}")
        console.log(f"[bold blue]Context.deps details:[/]\n{vars(ctx.deps)}")
        matched_files = []
        for root, _, files in os.walk(directory):
            for f in files:
                if isinstance(suffix, str):
                    condition = f.endswith(suffix)
                else:
                    condition = any(f.endswith(s) for s in suffix)
                if condition:
                    matched_files.append(os.path.join(root, f))
        console.log(f"[bold yellow]Progress:[/] Found {len(matched_files)} files matching suffix '{suffix}' in directory: {directory}")
        console.log(f"[bold green]Tool Completed:[/] find_files found {len(matched_files)} files")
        return sorted(matched_files)
    except FileNotFoundError:
        msg = f"Error: Directory '{directory}' not found."
        console.log(f"[bold red]Tool Error:[/] {msg}")
        return [msg]
    except Exception as e:
        error_msg = f"Error: {str(e)}"
        console.log(f"[bold red]Tool Exception:[/] {error_msg}")
        return [error_msg]

# Define a tool which uses find_files, always obtaining its directory from the dependency container.
@agent.tool
async def list_files(ctx: RunContext[TestDependency]) -> str:
    """
    Use the find_files tool to list files from the dependency's test_dir that match a given suffix.
    Here, we use the suffix '.txt' as an example.
    """
    directory = ctx.deps.test_dir  # get the directory from the dependency container
    suffix = ".txt"  # dynamic parameter, for example purposes
    files = await find_files(ctx, directory, suffix)
    return "The files found are: " + ", ".join(files)

if __name__ == "__main__":
    # Create a dependency instance with an example test directory.
    # (Adjust the directory value to one that exists in your testing environment.)
    test_dep = TestDependency(test_dir="../TestRNAseqData_SETBP1")

    # Create a dependency instance with an example test directory.
    # (Adjust the directory value to one that exists in your testing environment.)
    test_dep = TestDependency(test_dir="../TestRNAseqData_SETBP1")
    result = agent.run_sync("Please list the files with suffix .txt", deps=test_dep)
    console.print(result.data)
