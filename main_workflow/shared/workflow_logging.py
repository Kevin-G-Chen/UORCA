from __future__ import annotations
import logging, datetime, pathlib, inspect, functools, asyncio, os

__all__ = [
    "setup_logging",
    "log_tool",
]

DEFAULT_FMT = "%(asctime)s  %(levelname)-8s  %(name)s ▶  %(message)s"

################################################################################
# 1. Root configuration – called **once** at the start of the run              #
################################################################################

def setup_logging(log_dir: str | os.PathLike = "logs", *, level: int = logging.INFO,
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
