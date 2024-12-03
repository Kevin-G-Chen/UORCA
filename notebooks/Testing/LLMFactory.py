# This will contain experimentation about using LLMs to generate code for data extraction. Before, I pretty had no AI input (I did not see a reason to need to implement AI), however I want to explore the possibility of using LLMs to generate this code.
# Separately, I am also interested in using LLMFactory. This will allow me to produce Pydantic output regardless of LLM used, which will be extremely useful for assessing the performance of different LLMs.
# %% Cell 1
# Load modules
import sys
import os
import instructor
import httpx
from openai import OpenAI
import pandas as pd
from typing import Any, Dict, List, Type, Optional, Union
from pydantic import BaseModel, Field
from pydantic_settings import BaseSettings
from functools import lru_cache
from dotenv import load_dotenv
import time
from datetime import datetime
# %% Cell 2
# Define the LLM Factory class

class LLMFactory:
    def __init__(self, provider: str):
        self.provider = provider
        self.settings = getattr(get_settings(), provider)
        self.client = self._initialize_client()

    def _initialize_client(self) -> Any:
        client_initializers = {
            "openai": lambda s: instructor.from_openai(OpenAI(api_key=s.api_key)),
            "anthropic": lambda s: instructor.from_anthropic(
                Anthropic(api_key=s.api_key)
            ),
            "ollama": lambda s: instructor.patch(OpenAI(
                base_url=f"{s.base_url}/v1",
                api_key="not-needed"
            ), mode=instructor.Mode.JSON),
            "llama": lambda s: instructor.from_openai(
                OpenAI(base_url=s.base_url, api_key=s.api_key),
                mode=instructor.Mode.JSON,
            )
        }

        initializer = client_initializers.get(self.provider)
        if initializer:
            return initializer(self.settings)
        raise ValueError(f"Unsupported LLM provider: {self.provider}")

    def create_completion(
        self, response_model: Type[BaseModel], messages: List[Dict[str, str]], **kwargs
    ) -> Any:
        completion_params = {
            "model": kwargs.get("model", self.settings.default_model),
            "temperature": kwargs.get("temperature", self.settings.temperature),
            "max_retries": kwargs.get("max_retries", self.settings.max_retries),
            "max_tokens": kwargs.get("max_tokens", self.settings.max_tokens),
            "response_model": response_model,
            "messages": messages,
        }
        return self.client.chat.completions.create(**completion_params)


if __name__ == "__main__":

    class CompletionModel(BaseModel):
        response: str = Field(description="Your response to the user.")
        reasoning: str = Field(description="Explain your reasoning for the response.")

    messages = [
        {"role": "system", "content": "You are a helpful assistant."},
        {
            "role": "user",
            "content": "If it takes 2 hours to dry 1 shirt out in the sun, how long will it take to dry 5 shirts?",
        },
    ]
# %% Cell 3
# Define the provider settings
load_dotenv()

class LLMProviderSettings(BaseSettings):
    temperature: float = 0.0
    max_tokens: Optional[int] = None
    max_retries: int = 3
    base_url: Optional[str] = None

class OpenAISettings(LLMProviderSettings):
    api_key: str = os.getenv("OPENAI_API_KEY")
    default_model: str = "gpt-4o-mini"

class OllamaSettings(LLMProviderSettings):
    api_key: str = "not-needed"  # Ollama doesn't need an API key
    default_model: str = "llama2"  # Can be any model you have pulled via Ollama
    base_url: str = "http://localhost:11434"


class Settings(BaseSettings):
    app_name: str = "GenAI Project Template"
    openai: OpenAISettings = OpenAISettings()
    llama: LlamaSettings = LlamaSettings()
    ollama: OllamaSettings = OllamaSettings()


@lru_cache
def get_settings():
    return Settings()

# %% Cell 4
# Make an LLM call via LLMFactory
# I will also need to be making sure that I can call LLMs via Ollama - more generally, that I can call local LLMs
if __name__ == "__main__":

    class CompletionModel(BaseModel):
        response: str = Field(description="Your response to the user.")
        reasoning: str = Field(description="Explain your reasoning for the response.")

    messages = [
        {"role": "system", "content": "You are a helpful assistant."},
        {
            "role": "user",
            "content": "If it takes 2 hours to dry 1 shirt out in the sun, how long will it take to dry 5 shirts?",
        },
    ]

    print("\n=== Testing OpenAI ===")
    llm = LLMFactory("openai")
    start_time = time.time()

    try:
        completion = llm.create_completion(
            response_model=CompletionModel,
            messages=messages,
        )
        end_time = time.time()

        print(f"Model: {llm.settings.default_model}")
        print(f"Time taken: {end_time - start_time:.2f} seconds")
        if hasattr(completion, 'usage'):
            print(f"Input tokens: {completion.usage.prompt_tokens}")
            print(f"Output tokens: {completion.usage.completion_tokens}")
            print(f"Total tokens: {completion.usage.total_tokens}")
        print("\nResponse:", completion.response)
        print("\nReasoning:", completion.reasoning)

    except Exception as e:
        print(f"Error with OpenAI: {str(e)}")

    print("\n=== Testing Ollama ===")
    try:
        ollama_llm = LLMFactory("ollama")
        start_time = time.time()

        ollama_completion = ollama_llm.create_completion(
            response_model=CompletionModel,
            messages=messages,
            model="phi3"  # or any other model you've pulled via Ollama
        )
        end_time = time.time()

        print(f"Model: phi3")
        print(f"Time taken: {end_time - start_time:.2f} seconds")
        if hasattr(ollama_completion, 'usage'):
            print(f"Input tokens: {ollama_completion.usage.prompt_tokens}")
            print(f"Output tokens: {ollama_completion.usage.completion_tokens}")
            print(f"Total tokens: {ollama_completion.usage.total_tokens}")
        print("\nResponse:", ollama_completion.response)
        print("\nReasoning:", ollama_completion.reasoning)

    except Exception as e:
        print(f"Error with Ollama: {str(e)}")
        print("Make sure you have Ollama installed and running, and have pulled the required model")
# Well... being a smart hat aside... this works VERY well. I will be looking to implement this in my workflows.
