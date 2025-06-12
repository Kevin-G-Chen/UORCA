#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AI Landing Page for UORCA Explorer
---------------------------------
Interactive AI-powered interface for exploring UORCA RNA-seq analysis results.
"""

import asyncio
import json
import logging
from typing import Dict, List, Optional, Any

import streamlit as st
import pandas as pd

from ai_agent_factory import get_uorca_agent, validate_agent_setup

logger = logging.getLogger(__name__)

def render_ai_landing_page(results_dir: str):
    """
    Render the AI landing page for UORCA Explorer.

    Args:
        results_dir: Path to UORCA results directory
    """
    st.title("🤖 UORCA AI Assistant")
    st.markdown(
        """
        Welcome to the UORCA AI Assistant! Ask questions about your RNA-seq analysis results,
        explore gene expression patterns, and get insights into your data.
        """
    )

    # Initialize the AI agent
    try:
        with st.spinner("Initializing AI assistant..."):
            agent = get_uorca_agent(results_dir)

        # Validate agent setup
        if not validate_agent_setup(agent):
            st.error("AI assistant validation failed. Some features may not work properly.")
            return

        st.success("✅ AI assistant ready!")

    except Exception as e:
        st.error(f"Failed to initialize AI assistant: {e}")
        st.info("Please check your configuration and try again.")
        return

    # Create tabs for different interaction modes
    tab1, tab2, tab3 = st.tabs(["🎯 Quick Queries", "🔍 Data Explorer", "💬 Chat"])

    with tab1:
        render_quick_queries_tab(agent)

    with tab2:
        render_data_explorer_tab(agent)

    with tab3:
        render_chat_tab(agent)

def render_quick_queries_tab(agent):
    """Render the quick queries tab with pre-defined query templates."""
    st.header("Quick Data Queries")
    st.markdown("Get instant answers to common questions about your analysis.")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("📊 Overview Queries")

        if st.button("📈 How many contrasts do I have?", use_container_width=True):
            query_agent(agent, "How many contrasts are available in total? Just give me the number and a brief summary.")

        if st.button("🧬 How many genes were analyzed?", use_container_width=True):
            query_agent(agent, "How many unique genes are available in the analysis? Provide just the count.")

        if st.button("📋 List all my contrasts", use_container_width=True):
            query_agent(agent, "List all available contrasts in a well-formatted list. Group by analysis if multiple analyses are present.")

    with col2:
        st.subheader("🎯 Specific Queries")

        # Gene-specific query
        st.markdown("**Query specific gene:**")
        gene_input = st.text_input("Enter gene symbol:", placeholder="e.g., TP53, MYC, GAPDH", key="quick_gene")
        if st.button("🔍 Find this gene", disabled=not gene_input, use_container_width=True):
            query_agent(agent, f"Tell me about gene {gene_input}. In which contrasts is it present and what are its fold changes?")

        # Top genes query
        if st.button("⭐ Show me top upregulated genes", use_container_width=True):
            query_agent(agent, "What are some of the most highly upregulated genes across all contrasts? Show me the top examples with their fold changes.")

def render_data_explorer_tab(agent):
    """Render the interactive data explorer tab."""
    st.header("Interactive Data Explorer")
    st.markdown("Explore specific gene-contrast combinations with detailed information.")

    # Get available data for dropdowns
    with st.spinner("Loading available data..."):
        contrasts_data = query_agent_silent(agent, "List all contrasts as a simple JSON array")
        genes_data = query_agent_silent(agent, "List the first 100 genes as a simple JSON array")

    # Parse the data for dropdowns
    contrasts = parse_list_response(contrasts_data, "contrasts")
    genes = parse_list_response(genes_data, "genes")

    if not contrasts:
        st.warning("Could not load contrasts. Please check your data.")
        return

    col1, col2 = st.columns(2)

    with col1:
        selected_contrast = st.selectbox(
            "Select a contrast:",
            options=contrasts,
            help="Choose an experimental contrast to explore"
        )

    with col2:
        gene_option = st.radio(
            "Gene selection:",
            ["Choose from list", "Enter manually"],
            horizontal=True
        )

        if gene_option == "Choose from list":
            selected_gene = st.selectbox(
                "Select a gene:",
                options=genes if genes else ["No genes loaded"],
                disabled=not genes
            )
        else:
            selected_gene = st.text_input(
                "Enter gene symbol:",
                placeholder="e.g., TP53, MYC, GAPDH"
            )

    # Query button
    if st.button("🔬 Analyze Gene-Contrast Pair", use_container_width=True):
        if selected_contrast and selected_gene:
            query = f"Get the log fold change for gene {selected_gene} in contrast {selected_contrast}. Also explain what this value means biologically."
            query_agent(agent, query)
        else:
            st.warning("Please select both a contrast and a gene.")

    # Batch analysis section
    st.subheader("📊 Batch Gene Analysis")
    st.markdown("Analyze multiple genes at once in the selected contrast.")

    genes_text = st.text_area(
        "Enter gene symbols (one per line):",
        placeholder="TP53\nMYC\nGAPDH\nACTB",
        height=100
    )

    if st.button("🔍 Analyze Gene List", disabled=not genes_text or not selected_contrast):
        genes_list = [gene.strip() for gene in genes_text.split('\n') if gene.strip()]
        if genes_list:
            query = f"For contrast {selected_contrast}, get the log fold changes for these genes: {', '.join(genes_list)}. Present the results in a clear table format with biological interpretation."
            query_agent(agent, query)

def render_chat_tab(agent):
    """Render the free-form chat interface."""
    st.header("Chat with UORCA AI")
    st.markdown("Ask any question about your RNA-seq analysis in natural language.")

    # Initialize chat history
    if "chat_history" not in st.session_state:
        st.session_state.chat_history = []

    # Display chat history
    chat_container = st.container()
    with chat_container:
        for i, (query, response) in enumerate(st.session_state.chat_history):
            with st.chat_message("user"):
                st.write(query)
            with st.chat_message("assistant"):
                st.write(response)

    # Chat input
    user_input = st.chat_input("Ask me anything about your UORCA analysis...")

    if user_input:
        # Add user message to history
        with chat_container:
            with st.chat_message("user"):
                st.write(user_input)

        # Get AI response
        with chat_container:
            with st.chat_message("assistant"):
                response_placeholder = st.empty()
                with st.spinner("Thinking..."):
                    response = query_agent_silent(agent, user_input)
                    if response:
                        response_placeholder.write(response)
                        # Add to history
                        st.session_state.chat_history.append((user_input, response))
                    else:
                        response_placeholder.error("Sorry, I couldn't process that query.")

    # Clear chat button
    if st.button("🗑️ Clear Chat History"):
        st.session_state.chat_history = []
        st.rerun()

def query_agent(agent, query: str) -> Optional[str]:
    """
    Query the agent and display the response in Streamlit.

    Args:
        agent: The AI agent instance
        query: The query string

    Returns:
        The agent's response or None if failed
    """
    try:
        with st.status("🤖 AI is thinking...", expanded=True) as status:
            st.write(f"**Query:** {query}")

            # Run the async query
            response = asyncio.run(agent.run(query))

            st.write("**Response:**")
            st.write(response)

            status.update(label="✅ Complete!", state="complete")
            return response

    except Exception as e:
        logger.error(f"Query failed: {e}")
        st.error(f"❌ Query failed: {e}")
        return None

def query_agent_silent(agent, query: str) -> Optional[str]:
    """
    Query the agent without Streamlit UI elements (for background operations).

    Args:
        agent: The AI agent instance
        query: The query string

    Returns:
        The agent's response or None if failed
    """
    try:
        response = asyncio.run(agent.run(query))
        return response
    except Exception as e:
        logger.error(f"Silent query failed: {e}")
        return None

def parse_list_response(response: Optional[str], key: str) -> List[str]:
    """
    Parse a JSON response containing a list.

    Args:
        response: The agent's response
        key: The JSON key containing the list

    Returns:
        Parsed list or empty list if parsing fails
    """
    if not response:
        return []

    try:
        # Try to extract JSON from the response
        import re
        json_match = re.search(r'\{.*\}', response, re.DOTALL)
        if json_match:
            data = json.loads(json_match.group())
            return data.get(key, [])
        else:
            # If no JSON found, try to parse as direct list
            if '[' in response and ']' in response:
                list_match = re.search(r'\[.*\]', response, re.DOTALL)
                if list_match:
                    return json.loads(list_match.group())
    except Exception as e:
        logger.warning(f"Failed to parse list response: {e}")

    return []

def show_ai_landing_page_info():
    """Show information about the AI landing page feature."""
    st.info(
        """
        **🤖 AI Assistant Features:**

        - **Quick Queries**: Get instant answers to common questions
        - **Data Explorer**: Interactive exploration of gene-contrast combinations
        - **Chat Interface**: Natural language conversations about your data

        The AI assistant can help you:
        - Understand your analysis results
        - Find specific genes and their expression patterns
        - Compare conditions and contrasts
        - Interpret biological significance of findings
        """
    )
