"""
Configuration settings for the UORCA project.
"""
import os
from typing import Dict, Any, Optional
from pydantic import BaseSettings, Field


class UORCASettings(BaseSettings):
    """Configuration settings for UORCA application."""
    
    # API Keys and External Services
    llm_api_key: Optional[str] = Field(None, env="UORCA_LLM_API_KEY")
    llm_api_url: str = Field("https://api.openai.com/v1/chat/completions", env="UORCA_LLM_API_URL")
    llm_model: str = Field("gpt-4o", env="UORCA_LLM_MODEL")
    
    # NCBI GEO API Settings
    geo_api_base_url: str = Field("https://eutils.ncbi.nlm.nih.gov/entrez/eutils", env="UORCA_GEO_API_BASE_URL")
    geo_api_tool: str = Field("esearch,efetch,esummary", env="UORCA_GEO_API_TOOL")
    geo_api_db: str = Field("gds", env="UORCA_GEO_API_DB")  # GEO DataSets database
    geo_api_email: Optional[str] = Field(None, env="UORCA_GEO_API_EMAIL")  # Optional, but recommended by NCBI
    geo_api_rate_limit: float = Field(0.34, env="UORCA_GEO_API_RATE_LIMIT")  # Requests per second (3 per second)
    
    # Analysis Settings
    max_datasets_to_analyze: int = Field(3, env="UORCA_MAX_DATASETS")  # For test cases, limit to 3 datasets
    analysis_cache_dir: str = Field("./cache", env="UORCA_ANALYSIS_CACHE_DIR")
    
    # Development Settings
    debug_mode: bool = Field(False, env="UORCA_DEBUG_MODE")
    log_level: str = Field("INFO", env="UORCA_LOG_LEVEL")
    
    class Config:
        env_file = ".env"
        env_file_encoding = "utf-8"


# Create a global settings object
settings = UORCASettings()


def get_settings() -> UORCASettings:
    """Return the global settings object."""
    return settings
