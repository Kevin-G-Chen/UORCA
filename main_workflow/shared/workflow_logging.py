from __future__ import annotations
import logging, datetime, pathlib, inspect, functools, asyncio, os, json

__all__ = [
    "setup_logging",
    "log_tool",
    "log_tool_for_reflection",
]

DEFAULT_FMT = "%(asctime)s  %(levelname)-8s  %(name)s ▶  %(message)s"

################################################################################
# 1. Root configuration – called **once** at the start of the run              #
################################################################################

def setup_logging(log_dir: str | os.PathLike = "logs", *, level: int = logging.WARNING,
                  run_id: str | None = None) -> pathlib.Path:
    """Create *logs/YYYYMMDD_HHMMSS.log* and hook both file & console handlers.

    Returns the full path to the log file that was created so callers can show
    it (e.g. `print(f\"Log: {log_path}\")`).  Subsequent calls become NO‑OPs
    so agents can import this module safely without re‑initialising handlers.
    """
    if logging.getLogger().handlers:      # already configured → skip
        return pathlib.Path()

    ts = run_id or datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    log_dir = pathlib.Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{ts}.log"

    logging.basicConfig(
        level=level,
        format=DEFAULT_FMT,                 # uses %(asctime)s above
        datefmt="%Y-%m-%d %H:%M:%S",        # full timestamp
        handlers=[
            logging.FileHandler(log_file, encoding="utf‑8"),
            logging.StreamHandler(),
        ],
        force=True,                         # overrides any prior handlers
    )
    # Silence noisy deps (optional)
    logging.getLogger("urllib3").setLevel(logging.WARNING)
    logging.getLogger("openai").setLevel(logging.WARNING)

    return log_file

################################################################################
# 2.  Decorator – drop‑in per‑tool instrumentation                             #
################################################################################

def _as_dict(bound):
    return {k: v for k, v in bound.arguments.items() if k != "ctx"}


def log_tool(func):
    """Decorate any *agent.tool* function to get uniform **start / end / crash**
    messages.  Works for async + sync functions.
    """
    logger = logging.getLogger(func.__module__)

    if asyncio.iscoroutinefunction(func):

        @functools.wraps(func)
        async def wrapper(*args, **kwargs):
            bound = inspect.signature(func).bind_partial(*args, **kwargs)
            logger.info("🛠️  %s called – params=%s", func.__name__, _as_dict(bound))
            try:
                res = await func(*args, **kwargs)
                logger.info("✅  %s finished", func.__name__)
                return res
            except Exception:
                logger.exception("❌  %s crashed", func.__name__)
                raise

    else:

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            bound = inspect.signature(func).bind_partial(*args, **kwargs)
            logger.info("🛠️  %s called – params=%s", func.__name__, _as_dict(bound))
            try:
                res = func(*args, **kwargs)
                logger.info("✅  %s finished", func.__name__)
                return res
            except Exception:
                logger.exception("❌  %s crashed", func.__name__)
                raise

    return wrapper

################################################################################
# 3.  Agent tool decorator – like log_tool but omits parameters                #
################################################################################


def log_agent_tool(func):
    """Decorate agent-calling functions to get simplified logging without parameters.

    Similar to log_tool but doesn't log the parameters, which can be extremely large
    for agent calls (deps, prompt, etc).
    """
    logger = logging.getLogger(func.__module__)

    if asyncio.iscoroutinefunction(func):

        @functools.wraps(func)
        async def wrapper(*args, **kwargs):
            logger.info("🤖 %s called", func.__name__)
            try:
                res = await func(*args, **kwargs)
                logger.info("✅ %s finished", func.__name__)
                return res
            except Exception:
                logger.exception("❌ %s crashed", func.__name__)
                raise

    else:

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            logger.info("🤖 %s called", func.__name__)
            try:
                res = func(*args, **kwargs)
                logger.info("✅ %s finished", func.__name__)
                return res
            except Exception:
                logger.exception("❌ %s crashed", func.__name__)
                raise

    return wrapper

################################################################################
# 4.  Tool decorator that logs detailed information for reflection             #
################################################################################

def log_tool_for_reflection(func):
    """Decorate tools to log detailed information for reflection.
    
    Similar to log_tool but also records information in ctx.deps.tool_logs for
    later reflection. Particularly tracks:
    - Tool name
    - Parameters (excluding ctx)
    - Output
    - Success/failure status
    
    This information helps the reflection process understand what was attempted
    and with what results.
    """
    logger = logging.getLogger(func.__module__)

    if asyncio.iscoroutinefunction(func):

        @functools.wraps(func)
        async def wrapper(*args, **kwargs):
            # Get tool name and parameters (excluding ctx)
            tool_name = func.__name__
            bound = inspect.signature(func).bind_partial(*args, **kwargs)
            params = _as_dict(bound)
            
            # Extract ctx from args
            ctx = args[0] if args else None
            
            logger.info("🛠️ %s called for reflection – params=%s", tool_name, params)
            
            tool_log_entry = {
                "tool_name": tool_name,
                "parameters": params,
                "timestamp": datetime.datetime.now().isoformat(),
                "success": None,
                "output": None,
                "error": None
            }
            
            try:
                res = await func(*args, **kwargs)
                
                # Record successful output
                tool_log_entry["success"] = True
                tool_log_entry["output"] = res
                
                # Special tracking for important files
                if tool_name == "run_kallisto_quantification" and ctx and hasattr(ctx.deps, "kallisto_index_used"):
                    # Extract the index path from parameters
                    if "kallisto_index" in params:
                        ctx.deps.kallisto_index_used = params["kallisto_index"]
                
                if tool_name == "run_edger_limma_analysis" and ctx and hasattr(ctx.deps, "tx2gene_file_used"):
                    # Extract the tx2gene path from parameters
                    if "tx2gene_path" in params:
                        ctx.deps.tx2gene_file_used = params["tx2gene_path"]
                
                if tool_name == "list_files" and ctx:
                    # For list_files, simplify the output to just the count to avoid very large logs
                    if isinstance(res, list):
                        tool_log_entry["output"] = f"Found {len(res)} files"
                        if len(res) > 0 and len(res) <= 5:
                            tool_log_entry["output"] += f": {res}"
                        elif len(res) > 5:
                            tool_log_entry["output"] += f": {res[:5]} and {len(res)-5} more"
                
                logger.info("✅ %s finished for reflection", tool_name)
                
                # Add the log entry to context if available
                if ctx and hasattr(ctx.deps, "tool_logs"):
                    ctx.deps.tool_logs.append(tool_log_entry)
                
                return res
            except Exception as e:
                # Record error information
                tool_log_entry["success"] = False
                tool_log_entry["error"] = str(e)
                
                # Add the log entry to context even on failure
                if ctx and hasattr(ctx.deps, "tool_logs"):
                    ctx.deps.tool_logs.append(tool_log_entry)
                
                logger.exception("❌ %s crashed during reflection", tool_name)
                raise

    else:

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            # Get tool name and parameters (excluding ctx)
            tool_name = func.__name__
            bound = inspect.signature(func).bind_partial(*args, **kwargs)
            params = _as_dict(bound)
            
            # Extract ctx from args
            ctx = args[0] if args else None
            
            logger.info("🛠️ %s called for reflection – params=%s", tool_name, params)
            
            tool_log_entry = {
                "tool_name": tool_name,
                "parameters": params,
                "timestamp": datetime.datetime.now().isoformat(),
                "success": None,
                "output": None,
                "error": None
            }
            
            try:
                res = func(*args, **kwargs)
                
                # Record successful output
                tool_log_entry["success"] = True
                tool_log_entry["output"] = res
                
                # Special tracking for important files
                if tool_name == "run_kallisto_quantification" and ctx and hasattr(ctx.deps, "kallisto_index_used"):
                    # Extract the index path from parameters
                    if "kallisto_index" in params:
                        ctx.deps.kallisto_index_used = params["kallisto_index"]
                
                if tool_name == "run_edger_limma_analysis" and ctx and hasattr(ctx.deps, "tx2gene_file_used"):
                    # Extract the tx2gene path from parameters
                    if "tx2gene_path" in params:
                        ctx.deps.tx2gene_file_used = params["tx2gene_path"]
                
                if tool_name == "list_files" and ctx:
                    # For list_files, simplify the output to just the count to avoid very large logs
                    if isinstance(res, list):
                        tool_log_entry["output"] = f"Found {len(res)} files"
                        if len(res) > 0 and len(res) <= 5:
                            tool_log_entry["output"] += f": {res}"
                        elif len(res) > 5:
                            tool_log_entry["output"] += f": {res[:5]} and {len(res)-5} more"
                
                logger.info("✅ %s finished for reflection", tool_name)
                
                # Add the log entry to context if available
                if ctx and hasattr(ctx.deps, "tool_logs"):
                    ctx.deps.tool_logs.append(tool_log_entry)
                
                return res
            except Exception as e:
                # Record error information
                tool_log_entry["success"] = False
                tool_log_entry["error"] = str(e)
                
                # Add the log entry to context even on failure
                if ctx and hasattr(ctx.deps, "tool_logs"):
                    ctx.deps.tool_logs.append(tool_log_entry)
                
                logger.exception("❌ %s crashed during reflection", tool_name)
                raise

    return wrapper
