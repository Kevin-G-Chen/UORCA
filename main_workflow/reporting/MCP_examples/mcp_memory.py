#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
MCP (Model Context Protocol) memory server — **clean three‑table design**
-----------------------------------------------------------------------
A lightweight vector store backed by SQLite + sqlite‑vec exposing three
simple operations for Pydantic‑AI agents via MCP tools:

* **add_document**   (content + optional metadata → stores embedding)
* **search**         (k‑NN with optional metadata filter)
* **delete_document**

Everything else from the older server (logging, CLI wrapper, graceful
shutdown, etc.) is preserved. Only the storage backend and MCP‑tool
surface changed.
"""
import os, sys, json, signal, threading, sqlite3, hashlib
import argparse
import logging
#from contextlib import contextmanager
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import sqlite_vec                           # pip install sqlite-vec
from ollama import embed                    # or any embedding backend

from mcp.server.fastmcp import FastMCP
from autonify.config import ServerConfig

# -----------------------------------------------------------------------------
# MCP server setup & thin logging helpers
# -----------------------------------------------------------------------------
mcp = FastMCP("memory")
log  = lambda m: print(f"MCP MEMORY: {m}", file=sys.stderr, flush=True)


# -----------------------------------------------------------------------------
# Logging (stderr only)
# -----------------------------------------------------------------------------
os.environ["ANONYMIZED_TELEMETRY"] = "false"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s – %(name)s – %(levelname)s – %(message)s",
    stream=sys.stderr
)
logger = logging.getLogger("mcp-memory-server")

# -----------------------------------------------------------------------------
# Utility
# -----------------------------------------------------------------------------
def _resolve_db_path(raw_db: Optional[str] = None) -> str:
    raw = raw_db or os.environ.get("MEMORY_DB_PATH", "memory.db")
    if raw == ":memory:":
        return raw
    p = Path(raw).expanduser()
    if not p.is_absolute():
        p = Path.cwd() / p
    return str(p)

def _f32(v):                                   # → little‑endian float32 BLOB
    if isinstance(v, list):
        v = np.asarray(v, dtype=np.float32)
    return v.astype(np.float32).tobytes()

# -----------------------------------------------------------------------------
# Memory layer (three‑table)
# -----------------------------------------------------------------------------
class MemoryDB:
    """
    A lean, append‑only “memory” that stores **only** the final output
    documents of a graph run.  Each distinct document body appears once;
    its md5 hash is the primary key.  Embeddings are stored with
    sqlite‑vec so you can run k‑NN queries later.
    """

    # --------------------------------------------------------------------- #
    #  Construction / schema
    # --------------------------------------------------------------------- #
    def __init__(
        self,
        db_path: Optional[str] = None,
        model: str = os.getenv("EMBED_MODEL", "nomic-embed-text:latest"),
    ) -> None:
        self.db_path = _resolve_db_path(db_path)
        self.model = model
        self._tls = threading.local()

        # probe the embedding dimensionality once
        self.dim = len(self._embed("dimension_probe"))

        with self._conn() as c:
            # main documents table (one row per unique body)
            c.execute(
                """
                CREATE TABLE IF NOT EXISTS documents (
                    id      TEXT PRIMARY KEY,
                    content TEXT NOT NULL UNIQUE,
                    created TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                )
            """
            )

            # vector index (sqlite‑vec virtual table)
            c.execute(
                f"""
                CREATE VIRTUAL TABLE IF NOT EXISTS doc_embeddings
                USING vec0(embedding FLOAT[{self.dim}])
            """
            )

            # simple key‑value metadata table (many‑to‑one)
            c.execute(
                """
                CREATE TABLE IF NOT EXISTS document_metadata (
                    document_id  TEXT NOT NULL,
                    meta_key     TEXT NOT NULL,
                    meta_value   TEXT NOT NULL,
                    PRIMARY KEY (document_id, meta_key),
                    FOREIGN KEY (document_id)
                        REFERENCES documents(id)
                        ON DELETE CASCADE
                )
            """
            )
            c.execute(
                """
                CREATE INDEX IF NOT EXISTS idx_meta_key_value
                ON document_metadata(meta_key, meta_value);
                """
            )

            c.execute("PRAGMA foreign_keys = ON")
            c.commit()

    def _conn(self) -> sqlite3.Connection:
        """Thread‑local cached connection with sqlite‑vec loaded."""
        if not hasattr(self._tls, "conn"):
            conn = sqlite3.connect(self.db_path, timeout=30)
            conn.execute("PRAGMA journal_mode=WAL")
            conn.execute("PRAGMA busy_timeout = 5000;")      # wait up to 5 s on lock
            conn.execute("PRAGMA foreign_keys = ON;")
            conn.enable_load_extension(True)
            sqlite_vec.load(conn)
            conn.enable_load_extension(False)
            self._tls.conn = conn
        return self._tls.conn

    # --------------------------------------------------------------------- #
    #  Helpers
    # --------------------------------------------------------------------- #

    def _embed(self, text: str) -> np.ndarray:
        v = embed(model=self.model, input=text)["embeddings"]
        vec = np.asarray(v, dtype=np.float32)
        return vec[0] if vec.ndim == 2 else vec

    def _rowid(self, doc_id: str) -> int:
        row = self._conn().execute(
            "SELECT rowid FROM documents WHERE id = ?", (doc_id,)
        ).fetchone()
        if row is None:
            raise KeyError(f"Document id not found: {doc_id}")
        return row[0]

    # --------------------------------------------------------------------- #
    #  Public API
    # --------------------------------------------------------------------- #
    def add_document(
        self,
        content: str,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> str:
        """
        Persist an immutable output document.

        • Primary key is md5(content).
        • If the same body was stored before, we re‑use that row and just
          upsert metadata.
        • Embedding is written once (or repaired if it was missing).
        """
        if not content.strip():
            raise ValueError("Output document must have non‑empty content")

        db_id = hashlib.md5(content.encode()).hexdigest()
        conn = self._conn()

        with conn:
            # 1 · insert main row (ignored when body already seen)
            conn.execute(
                "INSERT OR IGNORE INTO documents(id, content) VALUES (?, ?)",
                (db_id, content),
            )

            # 2 · ensure an embedding exists
            rowid = self._rowid(db_id)
            vec_exists = conn.execute(
                "SELECT 1 FROM doc_embeddings WHERE rowid = ?", (rowid,)
            ).fetchone()
            if not vec_exists:
                vec_blob = _f32(self._embed(content))
                conn.execute(
                    "INSERT INTO doc_embeddings(rowid, embedding) VALUES (?, ?)",
                    (rowid, vec_blob),
                )

            # 3 · upsert metadata (last‑writer‑wins per key)
            for key, val in (metadata or {}).items():
                json_val = json.dumps(val, sort_keys=True, separators=(",",":"))
                conn.execute(
                  """
                  INSERT INTO document_metadata(document_id, meta_key, meta_value)
                  VALUES (?, ?, ?)
                  ON CONFLICT(document_id, meta_key)
                  DO UPDATE SET meta_value = excluded.meta_value
                  """,
                  (db_id, key, json_val),
                )

        return db_id


    # ------------------------------------------------------------- search ----
    # ------------------------------------------------------------------------- #
    #
    def _format_hit(self, row: Any) -> Dict[str, Any]:
        """
        Normalize a DB row into the dict with keys id, content, distance, metadata.
        Expects row = (id, content) or (id, content, distance).
        """
        doc_id, content, *rest = row
        dist = rest[0] if rest else None

        conn = self._conn()
        meta = {}
        for k, raw in conn.execute(
            "SELECT meta_key, meta_value FROM document_metadata WHERE document_id = ?",
            (doc_id,),
        ).fetchall():
            meta[k] = json.loads(raw)
        return {"id":doc_id, "content":content, "distance":dist, "metadata":meta}


    def _ids_for_metadata(self, metadata_filter: Dict[str,Any]) -> List[str]:
        if not metadata_filter:
            return []
        clauses, params = [], []
        for key, val in metadata_filter.items():
            clauses.append(
              "SELECT document_id FROM document_metadata WHERE meta_key=? AND meta_value=?"
            )
            params.extend([key,
              json.dumps(val, sort_keys=True, separators=(",",":"))
            ])
        sql = " INTERSECT ".join(clauses)
        rows = self._conn().execute(sql, params).fetchall()
        return [r[0] for r in rows]


    def calc_novelty(self, query: str, k: int = 15) -> float:
        """
        Perform a vector search and compute novelty (mean distance to the k neighbors).
        Returns:
          {
            "query":    <str>,
            "novelty":  <float>,    # average distance
            "hits":     [<hit dicts>]
          }
        """
        # 1) get the top-k hits with distances
        hits = self._vector_search(query, k)

        # 2) compute novelty as the mean distance
        if not hits:
            novelty = float("inf")
        else:
            total_dist = sum(hit["distance"] for hit in hits)
            novelty    = total_dist / len(hits)

        # 3) return everything
        return novelty

    def search(
            self,
            query: str,
            k: int = 5,
            metadata_filter: Optional[Dict[str, str]] = None,
        ) -> List[Dict[str, Any]]:
        """
        Route to one of three search modes:

        1) Pure vector search (query only).
        2) Metadata-only search (metadata_filter only).
        3) Hybrid search (both query and metadata_filter).
        """
        has_query = bool(query and query.strip())
        has_filter = bool(metadata_filter)

        if has_query and has_filter:
            return self._hybrid_search(query, metadata_filter, k)
        elif has_query:
            return self._vector_search(query, k)
        elif has_filter:
            return self._metadata_search(metadata_filter, k)
        else:
            raise ValueError("search() requires at least a query or a metadata_filter")

    # ------------------------------------------------------------------------- #
    # 1) Pure vector search
    # ------------------------------------------------------------------------- #
    def _vector_search(self, query: str, k: int) -> List[Dict[str, Any]]:
        """
        Return the top-k documents by cosine distance on their embeddings.
        """
        vec = json.dumps(self._embed(query).tolist())
        sql = """
            WITH knn AS (
                SELECT rowid, distance
                FROM doc_embeddings
                WHERE embedding MATCH json(?)
                ORDER BY distance
                LIMIT ?
            )
            SELECT d.id, d.content, knn.distance
                FROM knn
                JOIN documents d ON d.rowid = knn.rowid
                ORDER BY knn.distance
                LIMIT ?
        """
        conn = self._conn()
        rows = conn.execute(sql, (vec, k, k)).fetchall()
        return [self._format_hit(r) for r in rows]

    # ------------------------------------------------------------------------- #
    # 2) Metadata-only search
    # ------------------------------------------------------------------------- #
    def _metadata_search(
        self, metadata_filter: Dict[str, str], k: int
    ) -> List[Dict[str, Any]]:
        """
        Return up to k documents whose metadata match **all** key/value pairs.
        """
        # Helper to fetch all IDs matching the filter:
        candidate_ids = self._ids_for_metadata(metadata_filter)
        if not candidate_ids:
            return []

        # Limit to the first k
        ids = candidate_ids[:k]
        placeholders = ",".join("?" * len(ids))
        sql = f"SELECT id, content FROM documents WHERE id IN ({placeholders}) LIMIT ?"

        conn = self._conn()
        rows = conn.execute(sql, (*ids, k)).fetchall()
        return [self._format_hit((row[0], row[1], None)) for row in rows]

    # ------------------------------------------------------------------------- #
    # 3) Hybrid vector + metadata search
    # ------------------------------------------------------------------------- #
    def _hybrid_search(
        self,
        query: str,
        metadata_filter: Dict[str, str],
        k: int,
    ) -> List[Dict[str, Any]]:
        """
        First filter by metadata, then run a vector k-NN over the survivors.
        """
        # 1) metadata filter
        candidate_ids = self._ids_for_metadata(metadata_filter)
        if not candidate_ids:
            return []

        # 2) vector search over that set
        vec = json.dumps(self._embed(query).tolist())
        overshoot = k * 5
        placeholders = ",".join("?" * len(candidate_ids))

        sql = f"""
            WITH knn AS (
                SELECT rowid, distance
                FROM doc_embeddings
                WHERE embedding MATCH json(?)
                ORDER BY distance
                LIMIT ?
            )
            SELECT d.id, d.content, knn.distance
                FROM knn
                JOIN documents d ON d.rowid = knn.rowid
                WHERE d.id IN ({placeholders})
                ORDER BY knn.distance
                LIMIT ?
        """

        conn = self._conn()
        rows = conn.execute(
            sql, (vec, overshoot, *candidate_ids, k)
        ).fetchall()

        return [self._format_hit(r) for r in rows]


    # ----------------------------------------------------------- delete ----
    def delete_document(self, doc_id: str) -> bool:
        """Remove a document (plus its vector & metadata)."""
        with self._conn() as c:
            row = c.execute("SELECT rowid FROM documents WHERE id=?", (doc_id,)).fetchone()
            if not row:
                return False

            rowid = row[0]
            c.execute("DELETE FROM doc_embeddings WHERE rowid=?", (rowid,))
            cur = c.execute("DELETE FROM documents WHERE id=?", (doc_id,))  # <- keep cursor
            c.commit()
            return cur.rowcount > 0


    # (optional utils kept for completeness)
    def list_tags(self):
        conn = self._conn()

        rows = conn.execute(
            "SELECT DISTINCT meta_value FROM document_metadata WHERE meta_key='tag' AND meta_value!=''"
        ).fetchall()
        return [r[0] for r in rows]

# Shared singleton
_db: Optional[MemoryDB] = None
def get_db():
    global _db
    if _db is None:
        _db = MemoryDB()
    return _db

# -----------------------------------------------------------------------------
# MCP tools — **only add_document, search, delete_document**
# -----------------------------------------------------------------------------
@mcp.tool()
async def add_document(content: str,
                       metadata: Optional[Dict[str, str]] = None,
                       ) -> str:
    """
        Insert or update a document and its vector+metadata.

        Args:
            content: Non-empty text to index.
            metadata: Arbitrary JSON-serializable key/value pairs.

        Returns:
            The document ID (md5 hash).
    """
    try:
        return get_db().add_document(content, metadata or {})
    except Exception as exc:
        log(f"ERROR add_document: {exc}")
        return f"error: {exc}"

@mcp.tool()
async def search(query: str,
                 k: int = 5,
                 metadata_filter: Optional[Dict[str, str]] = None) -> str:
    """
    k-NN search over document embeddings with optional metadata filtering.

    Args:
    query: Free-text query.
    k:      Max number of hits.
    metadata_filter: Exact-match filters on metadata (keys → JSON values).

    Returns:
    JSON string:
    {
        "results": [
        {
            "id": "…",
            "content": "…",
            "distance": 0.123,
            "metadata": { … }
        },
        …
        ]
    }
    """
    try:
        hits = get_db().search(query, k, metadata_filter)
        return json.dumps({"results": hits}, ensure_ascii=False)
    except Exception as exc:
        log(f"ERROR search: {exc}")
        return json.dumps({"results": [], "error": str(exc)}, ensure_ascii=False)

@mcp.tool()
async def delete_document(doc_id: str) -> str:
    """
    Permanently remove a document (and its embedding+metadata).

    Args:
        doc_id: The ID returned by add_document.

    Returns:
        "ok"        – document was deleted
        "not_found" – no such ID existed
        "error:…"   – on failure
    """
    try:
        ok = get_db().delete_document(doc_id)
        return "ok" if ok else "not_found"
    except Exception as exc:
        log(f"ERROR delete_document: {exc}")
        return f"error: {exc}"

# -----------------------------------------------------------------------------
# Graceful shutdown + CLI unchanged
# -----------------------------------------------------------------------------

def _signal_handler(sig, frame):
    log(f"received signal {sig} – exiting")
    sys.exit(0)

signal.signal(signal.SIGINT, _signal_handler)
signal.signal(signal.SIGTERM, _signal_handler)

# -----------------------------------------------------------------------------
# CLI entry point
# -----------------------------------------------------------------------------
def load_configuration():
    try:
        cfg = ServerConfig().get_server_config("memory") or {}
        return {"verbose": cfg.get("config", {}).get("verbose", False)}
    except Exception as e:
        logger.warning(f"config load error: {e}")
        return {"verbose": False}

def main():
    parser = argparse.ArgumentParser(description="MCP Memory server")
    parser.add_argument(
        "--transport", default="stdio", choices=["stdio", "sse"],
        help="Transport protocol (stdio or server-sent events)"
    )
    parser.add_argument(
        "--debug", action="store_true", help="Enable debug logging"
    )
    args = parser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
        cfg = load_configuration()
        log(f"using configuration: {cfg}")

    mcp.run(transport=args.transport)

if __name__ == "__main__":
    main()
