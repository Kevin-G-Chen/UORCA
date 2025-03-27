import os
import asyncio
from lightrag import LightRAG, QueryParam
from lightrag.llm.openai import gpt_4o_mini_complete, openai_embed
from lightrag.kg.shared_storage import initialize_pipeline_status

WORKING_DIR = "./seals"

if not os.path.exists(WORKING_DIR):
    os.mkdir(WORKING_DIR)


async def initialize_rag():
    rag = LightRAG(
        working_dir=WORKING_DIR,
        embedding_func=openai_embed,
        llm_model_func=gpt_4o_mini_complete,
    )

    await rag.initialize_storages()
    await initialize_pipeline_status()
    return rag


def main():
    # Initialize RAG instance
    rag = asyncio.run(initialize_rag())

    # Read the plain text from the file created in the previous step
    with open("./seal.txt", "r", encoding="utf-8") as f:
        rag.insert(f.read())

    # Now you can perform your queries with LightRAG
    print(
        rag.query(
            "Can you compare and contrast facts about seals?", param=QueryParam(mode="naive")
        )
    )
    print(
        rag.query(
            "Can you compare and contrast facts about seals?", param=QueryParam(mode="local")
        )
    )
    print(
        rag.query(
            "Can you compare and contrast facts about seals?", param=QueryParam(mode="global")
        )
    )
    print(
        rag.query(
            "Can you compare and contrast facts about seals?", param=QueryParam(mode="hybrid")
        )
    )


if __name__ == "__main__":
    main()
