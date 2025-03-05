# %% [markdown]
# # UORCA - Unified Omics Research and Computational Analysis
# 
# This notebook implements an agentic workflow for automated analysis of NCBI GEO RNA-seq datasets. The workflow:
# 
# 1. Processes a research question
# 2. Identifies relevant GEO datasets
# 3. Retrieves and parses dataset metadata
# 4. Designs and executes analysis
# 5. Synthesizes findings across datasets
# 
# ## Setup and Dependencies

# %% [imports]
import os
import re
import json
import time
import pandas as pd
import numpy as np
import requests
import xmltodict
import seaborn as sns
import matplotlib.pyplot as plt
from typing import List, Dict, Optional, Union, Any, Set, Tuple
from pydantic import BaseModel, Field, validator
from enum import Enum
import logging
from tenacity import retry, stop_after_attempt, wait_exponential
from concurrent.futures import ThreadPoolExecutor
from io import StringIO
import GEOparse

# For PydanticAI (AI agents framework)
from langchain.chat_models import ChatOpenAI
from langchain.schema import HumanMessage, SystemMessage
from langchain.chains import LLMChain
from langchain.prompts import ChatPromptTemplate, HumanMessagePromptTemplate, SystemMessagePromptTemplate
from langchain.output_parsers import PydanticOutputParser

# %% [logging_setup]
# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("UORCA")

# %% [config]
# Configuration settings
class UORCAConfig(BaseModel):
    """Configuration for the UORCA system"""
    # API settings
    openai_model: str = "gpt-4-0125-preview"  # You can replace with your model choice
    geo_base_url: str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    
    # Runtime settings
    max_datasets: int = 5  # Maximum number of datasets to analyze
    max_retry_attempts: int = 3  # Maximum number of retry attempts
    
    # Analysis settings
    p_value_threshold: float = 0.05
    log2fc_threshold: float = 1.0  # Log2 fold change threshold
    
    # Output settings
    output_dir: str = "./uorca_results"
    
    @validator('output_dir')
    def create_output_dir(cls, v):
        os.makedirs(v, exist_ok=True)
        return v

# Initialize config
config = UORCAConfig()

# Make sure the output directory exists
os.makedirs(config.output_dir, exist_ok=True)

# %% [markdown]
# ## Core Data Models
# 
# These Pydantic models define the structured data used throughout the workflow.

# %% [data_models]
class OmicsDataType(str, Enum):
    """Types of omics data"""
    BULK_RNA_SEQ = "bulk_rna_seq"
    SINGLE_CELL_RNA_SEQ = "single_cell_rna_seq"
    PROTEOMICS = "proteomics"
    GENOMICS = "genomics"
    METHYLATION = "methylation"
    OTHER = "other"

class DiseaseContext(BaseModel):
    """Disease context for a research question"""
    name: str
    synonyms: List[str] = []
    subtypes: List[str] = []
    related_conditions: List[str] = []

class ResearchQuestion(BaseModel):
    """Structured representation of a biomedical research question"""
    question: str
    disease_contexts: List[DiseaseContext] = []
    target_genes: List[str] = []
    biological_processes: List[str] = []
    omics_types: List[OmicsDataType] = [OmicsDataType.BULK_RNA_SEQ]
    
    # Search terms derived from the research question
    search_terms: List[str] = []
    
    class Config:
        schema_extra = {
            "example": {
                "question": "What genes are differentially expressed in lung cancer compared to normal tissue?",
                "disease_contexts": [{"name": "lung cancer", "synonyms": ["lung carcinoma", "pulmonary cancer"]}],
                "target_genes": [],
                "biological_processes": ["oncogenesis", "tumor progression"],
                "omics_types": ["bulk_rna_seq"],
                "search_terms": ["lung cancer", "RNA-seq", "differential expression"]
            }
        }

class DatasetMetadata(BaseModel):
    """Metadata for a GEO dataset"""
    geo_accession: str
    title: str
    summary: str
    organism: List[str]
    overall_design: Optional[str] = None
    platform_id: Optional[str] = None
    sample_count: int = 0
    supplementary_files: List[str] = []
    condition_groups: Dict[str, List[str]] = {}  # Group name -> list of sample IDs
    raw_metadata: Dict = {}  # Store the full raw metadata
    
    # Relevance score and reasoning provided by an agent
    relevance_score: float = 0.0
    relevance_reasoning: str = ""

class SampleMetadata(BaseModel):
    """Metadata for a GEO sample"""
    gsm_id: str
    title: str
    source_name: Optional[str] = None
    organism: Optional[str] = None
    characteristics: Dict[str, str] = {}  # Key-value pairs of characteristics
    condition: Optional[str] = None  # The condition this sample belongs to
    
    @validator('characteristics', pre=True)
    def parse_characteristics(cls, v):
        if isinstance(v, list):
            result = {}
            for item in v:
                if isinstance(item, dict) and 'tag' in item and 'value' in item:
                    key = item['tag'].lower().strip()
                    value = item['value'].strip()
                    result[key] = value
            return result
        return v

class AnalysisDesign(BaseModel):
    """Design matrix and parameters for differential expression analysis"""
    dataset_id: str
    control_group: str
    experimental_group: str
    sample_mapping: Dict[str, str] = {}  # Sample ID -> Group
    covariates: List[str] = []
    code_snippet: str = ""  # Generated analysis code

class AnalysisResult(BaseModel):
    """Results from analyzing a single dataset"""
    dataset_id: str
    analysis_design: AnalysisDesign
    deg_count: int = 0
    top_degs: List[Dict[str, Any]] = []  # List of top differentially expressed genes
    enriched_pathways: List[Dict[str, Any]] = []  # List of enriched pathways
    execution_log: str = ""  # Log of the analysis execution
    
    # Metrics for eval
    success: bool = False
    error_message: Optional[str] = None

class Insight(BaseModel):
    """An insight derived from analyzing multiple datasets"""
    title: str
    description: str
    supporting_datasets: List[str] = []
    supporting_evidence: List[str] = []
    genes_involved: List[str] = []
    confidence_score: float = 0.0  # 0.0 to 1.0
    
    # For evaluation
    is_validated: bool = False
    validation_method: Optional[str] = None

class InsightCollection(BaseModel):
    """A collection of insights derived from multiple datasets"""
    research_question: ResearchQuestion
    insights: List[Insight] = []
    common_degs: List[str] = []
    dataset_summary: Dict[str, Any] = {}
    
    def add_insight(self, insight: Insight):
        self.insights.append(insight)
    
    def get_insights_by_confidence(self, min_confidence: float = 0.7) -> List[Insight]:
        return [i for i in self.insights if i.confidence_score >= min_confidence]

# %% [markdown]
# ## Agent Definitions
# 
# The following classes define the agentic components of the UORCA workflow.

# %% [agent_base]
class BaseAgent(BaseModel):
    """Base class for all agents in the UORCA system"""
    name: str
    description: str
    llm: Any = None
    
    class Config:
        arbitrary_types_allowed = True
    
    def initialize_llm(self, model_name: str = None):
        """Initialize the language model for this agent"""
        model = model_name or config.openai_model
        self.llm = ChatOpenAI(
            model=model,
            temperature=0.1
        )
        return self.llm
    
    @retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=2, max=10))
    def execute_with_retry(self, func, *args, **kwargs):
        """Execute a function with retry logic"""
        try:
            return func(*args, **kwargs)
        except Exception as e:
            logger.error(f"Error in {self.name}: {str(e)}")
            raise

# %% [research_question_agent]
class ResearchQuestionAgent(BaseAgent):
    """Agent responsible for understanding and structuring the research question"""
    name: str = "ResearchQuestionAnalyzer"
    description: str = "Analyzes and structures biomedical research questions"
    
    def __init__(self, **data):
        super().__init__(**data)
        self.initialize_llm()
    
    def analyze_question(self, question_text: str) -> ResearchQuestion:
        """Analyze a research question and convert it to a structured format"""
        logger.info(f"Analyzing research question: {question_text}")
        
        # Create a parser for the ResearchQuestion model
        parser = PydanticOutputParser(pydantic_object=ResearchQuestion)
        
        # Define the prompt
        system_message = SystemMessagePromptTemplate.from_template(
            """You are an expert in biomedical research. Your task is to analyze a research question and 
            extract key components including diseases, genes, biological processes, and relevant search terms.
            
            For RNA-seq studies, think about:
            - What disease or condition is being studied?
            - What specific genes or pathways might be relevant?
            - What biological processes are involved?
            - What terms would be useful for searching GEO datasets?
            
            Format your response according to this JSON schema: {format_instructions}
            """
        )
        
        human_message = HumanMessagePromptTemplate.from_template(
            "Research Question: {question}\n\nPlease analyze this question and provide a structured response."
        )
        
        chat_prompt = ChatPromptTemplate.from_messages([system_message, human_message])
        
        # Format the prompt with the parser instructions and question
        formatted_prompt = chat_prompt.format_prompt(
            format_instructions=parser.get_format_instructions(),
            question=question_text
        ).to_messages()
        
        # Get the response from the LLM
        response = self.llm(formatted_prompt)
        
        # Parse the response into a ResearchQuestion object
        try:
            research_question = parser.parse(response.content)
            
            # If search terms are empty, generate them from other fields
            if not research_question.search_terms:
                search_terms = []
                for disease in research_question.disease_contexts:
                    search_terms.append(disease.name)
                    search_terms.extend(disease.synonyms[:2])  # Add first two synonyms if available
                
                search_terms.extend(research_question.target_genes[:5])  # Add up to 5 genes
                search_terms.extend(research_question.biological_processes[:3])  # Add up to 3 processes
                search_terms.append("RNA-seq")  # Always include RNA-seq since that's our focus
                
                # Remove duplicates and empty strings
                research_question.search_terms = list(set(filter(None, search_terms)))
            
            return research_question
        
        except Exception as e:
            logger.error(f"Failed to parse research question: {e}")
            logger.debug(f"LLM Response: {response.content}")
            
            # Create a basic research question with just the text
            return ResearchQuestion(
                question=question_text,
                disease_contexts=[],
                omics_types=[OmicsDataType.BULK_RNA_SEQ],
                search_terms=[question_text, "RNA-seq"]
            )

# %% [dataset_finder_agent]
class DatasetFinderAgent(BaseAgent):
    """Agent responsible for finding relevant GEO datasets"""
    name: str = "DatasetFinder"
    description: str = "Searches for relevant GEO datasets based on research questions"
    
    def __init__(self, **data):
        super().__init__(**data)
        self.initialize_llm()
    
    @retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=4, max=20))
    def search_geo_datasets(self, search_term: str, retmax: int = 20) -> Dict:
        """Search GEO for datasets matching the search term"""
        logger.info(f"Searching GEO for: {search_term}")
        
        # Construct the search URL
        base_url = f"{config.geo_base_url}/esearch.fcgi"
        params = {
            "db": "gds",  # GEO DataSets
            "term": f"{search_term} AND gse[ETYP] AND \"expression profiling by high throughput sequencing\"[DataSet Type]",
            "retmax": retmax,
            "retmode": "json"
        }
        
        response = requests.get(base_url, params=params)
        
        if response.status_code != 200:
            logger.error(f"GEO API error: {response.status_code}")
            return {"count": 0, "ids": []}
        
        try:
            data = response.json()
            return {
                "count": int(data["esearchresult"]["count"]),
                "ids": data["esearchresult"].get("idlist", [])
            }
        except (json.JSONDecodeError, KeyError) as e:
            logger.error(f"Failed to parse GEO search results: {e}")
            return {"count": 0, "ids": []}
    
    @retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=4, max=20))
    def fetch_dataset_summary(self, dataset_id: str) -> Dict:
        """Fetch summary information for a specific GEO dataset"""
        logger.info(f"Fetching summary for dataset: {dataset_id}")
        
        # Construct the URL
        base_url = f"{config.geo_base_url}/esummary.fcgi"
        params = {
            "db": "gds",
            "id": dataset_id,
            "retmode": "json"
        }
        
        response = requests.get(base_url, params=params)
        
        if response.status_code != 200:
            logger.error(f"GEO API error: {response.status_code}")
            return {}
        
        try:
            data = response.json()
            return data["result"][dataset_id]
        except (json.JSONDecodeError, KeyError) as e:
            logger.error(f"Failed to parse GEO summary: {e}")
            return {}
    
    def find_datasets(self, research_question: ResearchQuestion) -> List[str]:
        """Find GEO datasets relevant to the research question"""
        all_dataset_ids = set()
        
        # Search for each term in the research question
        for term in research_question.search_terms:
            result = self.search_geo_datasets(term)
            all_dataset_ids.update(result["ids"])
        
        # Convert set back to list for consistent ordering
        return list(all_dataset_ids)
    
    def evaluate_dataset_relevance(self, dataset_summary: Dict, research_question: ResearchQuestion) -> Tuple[float, str]:
        """
        Evaluate the relevance of a dataset to the research question
        Returns a tuple of (score, reasoning)
        """
        # Define the prompt
        system_message = SystemMessagePromptTemplate.from_template(
            """You are an expert in bioinformatics and genomics. 
            Evaluate how relevant this GEO dataset is to the provided research question.
            
            Consider these factors:
            1. Does the dataset study the same disease/condition?
            2. Does it use the appropriate omics method (RNA-seq)?
            3. Does it have sufficient samples for statistical analysis?
            4. Does the experimental design match what we're looking for?
            
            Rate relevance from 0.0 (not relevant) to 1.0 (highly relevant).
            Format your response as JSON with two fields:
            {
                "relevance_score": [float between 0.0-1.0],
                "reasoning": [brief explanation of your rating]
            }
            """
        )
        
        human_message = HumanMessagePromptTemplate.from_template(
            """Research Question: {question}
            
            Dataset Information:
            Title: {title}
            Summary: {summary}
            Organism: {organism}
            Samples: {samples}
            
            Evaluate the relevance of this dataset to the research question.
            """
        )
        
        chat_prompt = ChatPromptTemplate.from_messages([system_message, human_message])
        
        # Format the question and dataset info
        disease_contexts = ", ".join([d.name for d in research_question.disease_contexts])
        organism = dataset_summary.get("organism", "Not specified")
        samples = dataset_summary.get("n_samples", "Unknown")
        
        formatted_prompt = chat_prompt.format_prompt(
            question=f"{research_question.question} (Focuses on: {disease_contexts})",
            title=dataset_summary.get("title", ""),
            summary=dataset_summary.get("summary", ""),
            organism=organism,
            samples=samples
        ).to_messages()
        
        # Get the response
        response = self.llm(formatted_prompt)
        
        # Parse the response
        try:
            result = json.loads(response.content)
            score = float(result["relevance_score"])
            reasoning = result["reasoning"]
            return (score, reasoning)
        except (json.JSONDecodeError, KeyError) as e:
            logger.error(f"Failed to parse relevance evaluation: {e}")
            return (0.0, "Failed to evaluate relevance")
    
    def get_top_relevant_datasets(self, research_question: ResearchQuestion, max_datasets: int = None) -> List[DatasetMetadata]:
        """Find and rank datasets by relevance to the research question"""
        if max_datasets is None:
            max_datasets = config.max_datasets
            
        logger.info(f"Searching for datasets relevant to: {research_question.question}")
        
        # Find candidate datasets
        dataset_ids = self.find_datasets(research_question)
        logger.info(f"Found {len(dataset_ids)} candidate datasets")
        
        if not dataset_ids:
            logger.warning("No datasets found")
            return []
        
        # Limit to prevent too many API calls during development
        dataset_ids = dataset_ids[:min(30, len(dataset_ids))]
        
        # Fetch and evaluate each dataset
        datasets = []
        for dataset_id in dataset_ids:
            summary = self.fetch_dataset_summary(dataset_id)
            if not summary:
                continue
                
            # Evaluate relevance
            relevance_score, relevance_reasoning = self.evaluate_dataset_relevance(summary, research_question)
            
            # Create dataset metadata
            geo_accession = summary.get("accession", "")
            
            # Skip if no accession (shouldn't happen)
            if not geo_accession:
                continue
                
            dataset = DatasetMetadata(
                geo_accession=geo_accession,
                title=summary.get("title", ""),
                summary=summary.get("summary", ""),
                organism=[summary.get("organism", "")],
                sample_count=int(summary.get("n_samples", 0)),
                raw_metadata=summary,
                relevance_score=relevance_score,
                relevance_reasoning=relevance_reasoning
            )
            
            datasets.append(dataset)
        
        # Sort by relevance score
        datasets.sort(key=lambda x: x.relevance_score, reverse=True)
        
        # Return top datasets
        return datasets[:max_datasets]

# %% [dataset_analyzer_agent]
class DatasetAnalyzerAgent(BaseAgent):
    """Agent responsible for analyzing GEO datasets"""
    name: str = "DatasetAnalyzer"
    description: str = "Analyzes GEO datasets and extracts metadata"
    
    def __init__(self, **data):
        super().__init__(**data)
        self.initialize_llm()
    
    @retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=4, max=20))
    def fetch_geo_dataset(self, geo_accession: str) -> Any:
        """Fetch a GEO dataset using GEOparse"""
        logger.info(f"Fetching GEO dataset: {geo_accession}")
        try:
            return GEOparse.get_GEO(geo_accession)
        except Exception as e:
            logger.error(f"Failed to fetch GEO dataset {geo_accession}: {e}")
            return None
    
    def parse_dataset_metadata(self, gse: Any) -> DatasetMetadata:
        """Parse a GSE object into DatasetMetadata"""
        logger.info(f"Parsing metadata for dataset: {gse.name}")
        
        # Get basic metadata
        metadata = DatasetMetadata(
            geo_accession=gse.name,
            title=gse.metadata.get("title", ["Unknown"])[0],
            summary=gse.metadata.get("summary", ["No summary available"])[0],
            organism=gse.metadata.get("organism_ch1", ["Unknown"]),
            overall_design=gse.metadata.get("overall_design", [""])[0],
            platform_id=list(gse.gpls.keys())[0] if gse.gpls else None,
            sample_count=len(gse.gsms),
            supplementary_files=gse.metadata.get("supplementary_file", []),
            raw_metadata=gse.metadata
        )
        
        return metadata
    
    def parse_sample_metadata(self, gse: Any) -> List[SampleMetadata]:
        """Parse sample metadata from a GSE object"""
        samples = []
        
        for gsm_id, gsm in gse.gsms.items():
            sample = SampleMetadata(
                gsm_id=gsm_id,
                title=gsm.metadata.get("title", ["Unknown"])[0],
                source_name=gsm.metadata.get("source_name_ch1", ["Unknown"])[0] if "source_name_ch1" in gsm.metadata else None,
                organism=gsm.metadata.get("organism_ch1", ["Unknown"])[0] if "organism_ch1" in gsm.metadata else None,
                characteristics={}
            )
            
            # Parse characteristics
            characteristics = {}
            for key in gsm.metadata:
                if key.startswith("characteristics_ch1"):
                    for char in gsm.metadata[key]:
                        parts = char.split(":", 1)
                        if len(parts) == 2:
                            char_name = parts[0].strip().lower()
                            char_value = parts[1].strip()
                            characteristics[char_name] = char_value
            
            sample.characteristics = characteristics
            samples.append(sample)
        
        return samples
    
    def identify_condition_groups(self, samples: List[SampleMetadata]) -> Dict[str, List[str]]:
        """
        Identify experimental condition groups from sample metadata
        Returns a dict mapping condition names to lists of sample IDs
        """
        # Create a prompt to identify conditions
        system_message = SystemMessagePromptTemplate.from_template(
            """You are an expert in bioinformatics and experimental design.
            Analyze these RNA-seq sample metadata to identify distinct experimental conditions or groups.
            
            Your task:
            1. Look for characteristics that define experimental conditions (e.g., disease/healthy, treated/untreated)
            2. Identify the key characteristic that separates samples into meaningful groups
            3. Create logical groupings based on this characteristic
            
            Format your response as JSON:
            {
                "condition_field": "the characteristic field name that defines conditions",
                "groups": {
                    "group_name_1": ["sample_id1", "sample_id2", ...],
                    "group_name_2": ["sample_id3", "sample_id4", ...],
                    ...
                },
                "reasoning": "Brief explanation of why you chose this grouping"
            }
            
            Focus on creating clean control vs. experimental groups. Aim for 2-4 groups maximum.
            """
        )
        
        human_message = HumanMessagePromptTemplate.from_template(
            """Here are the RNA-seq samples from a GEO dataset:
            
            {samples_json}
            
            Identify the meaningful experimental conditions/groups in this dataset.
            """
        )
        
        chat_prompt = ChatPromptTemplate.from_messages([system_message, human_message])
        
        # Format the samples as JSON
        samples_json = json.dumps([{
            "sample_id": s.gsm_id,
            "title": s.title,
            "source": s.source_name,
            "characteristics": s.characteristics
        } for s in samples], indent=2)
        
        formatted_prompt = chat_prompt.format_prompt(
            samples_json=samples_json
        ).to_messages()
        
        # Get the response
        response = self.llm(formatted_prompt)
        
        # Parse the response
        try:
            result = json.loads(response.content)
            groups = result.get("groups", {})
            
            # Update samples with their assigned condition
            for group_name, sample_ids in groups.items():
                for sample_id in sample_ids:
                    for sample in samples:
                        if sample.gsm_id == sample_id:
                            sample.condition = group_name
            
            return groups
        
        except (json.JSONDecodeError, KeyError) as e:
            logger.error(f"Failed to parse condition groups: {e}")
            
            # Fallback: create a basic grouping based on title patterns
            groups = {"unclassified": []}
            
            # Simple heuristic matching
            for sample in samples:
                title_lower = sample.title.lower()
                
                if "control" in title_lower or "normal" in title_lower or "healthy" in title_lower:
                    if "control" not in groups:
                        groups["control"] = []
                    groups["control"].append(sample.gsm_id)
                    sample.condition = "control"
                elif "patient" in title_lower or "disease" in title_lower or "tumor" in title_lower or "cancer" in title_lower:
                    if "disease" not in groups:
                        groups["disease"] = []
                    groups["disease"].append(sample.gsm_id)
                    sample.condition = "disease"
                else:
                    groups["unclassified"].append(sample.gsm_id)
                    sample.condition = "unclassified"
            
            # Remove empty groups
            return {k: v for k, v in groups.items() if v}
    
    def create_analysis_design(self, dataset_id: str, condition_groups: Dict[str, List[str]]) -> List[AnalysisDesign]:
        """
        Create analysis design matrices for differential expression analysis
        Returns potential designs for control vs. experimental comparisons
        """
        # If we have exactly two groups, it's straightforward
        if len(condition_groups) == 2:
            groups = list(condition_groups.keys())
            
            # Find which one is likely the control
            control_candidates = ["control", "normal", "healthy", "wt", "wildtype", "wild_type", "untreated", "mock"]
            control_group = groups[0]
            exp_group = groups[1]
            
            for cand in control_candidates:
                for group in groups:
                    if cand in group.lower():
                        control_group = group
                        exp_group = [g for g in groups if g != control_group][0]
                        break
            
            # Create a sample mapping
            sample_mapping = {}
            for group, samples in condition_groups.items():
                for sample in samples:
                    sample_mapping[sample] = group
            
            return [AnalysisDesign(
                dataset_id=dataset_id,
                control_group=control_group,
                experimental_group=exp_group,
                sample_mapping=sample_mapping
            )]
            
        # If we have more than two groups, create pairwise comparisons with the control
        elif len(condition_groups) > 2:
            designs = []
            
            # Try to identify the control group
            control_candidates = ["control", "normal", "healthy", "wt", "wildtype", "wild_type", "untreated", "mock"]
            control_group = None
            
            for cand in control_candidates:
                for group in condition_groups.keys():
                    if cand in group.lower():
                        control_group = group
                        break
                if control_group:
                    break
            
            # If no control identified, use the first group
            if control_group is None:
                control_group = list(condition_groups.keys())[0]
            
            # Create pairwise designs
            for exp_group in condition_groups.keys():
                if exp_group != control_group:
                    # Create a sample mapping
                    sample_mapping = {}
                    for sample in condition_groups[control_group]:
                        sample_mapping[sample] = control_group
                    for sample in condition_groups[exp_group]:
                        sample_mapping[sample] = exp_group
                    
                    designs.append(AnalysisDesign(
                        dataset_id=dataset_id,
                        control_group=control_group,
                        experimental_group=exp_group,
                        sample_mapping=sample_mapping
                    ))
            
            return designs
            
        # If we only have one group, we can't do DE analysis
        else:
            logger.warning(f"Dataset {dataset_id} has only one condition group, can't create analysis design")
            return []
    
    def generate_analysis_code(self, analysis_design: AnalysisDesign) -> str:
        """Generate analysis code for differential expression analysis"""
        # Create a prompt for code generation
        system_message = SystemMessagePromptTemplate.from_template(
            """You are an expert bioinformatician who specializes in RNA-seq analysis. 
            Generate a Python code snippet that performs differential expression analysis using DESeq2 in R 
            (via rpy2) for the provided experimental design.
            
            The code should:
            1. Set up the DESeq2 design matrix based on the provided sample mapping
            2. Run differential expression analysis
            3. Extract the results with adjusted p-value and log2 fold change thresholds
            4. Return the list of DEGs sorted by significance
            
            Use these thresholds:
            - Adjusted p-value < 0.05
            - |log2FoldChange| > 1.0
            
            Assume the count matrix is available as 'counts_df' with genes as rows and samples as columns.
            """
        )
        
        human_message = HumanMessagePromptTemplate.from_template(
            """Please generate a code snippet for differential expression analysis with this design:
            
            Dataset: {dataset_id}
            Control group: {control_group}
            Experimental group: {experimental_group}
            Sample mapping:
            {sample_mapping}
            
            Generate readable, well-commented R code for DESeq2 analysis via rpy2.
            """
        )
        
        chat_prompt = ChatPromptTemplate.from_messages([system_message, human_message])
        
        # Format the analysis design
        formatted_prompt = chat_prompt.format_prompt(
            dataset_id=analysis_design.dataset_id,
            control_group=analysis_design.control_group,
            experimental_group=analysis_design.experimental_group,
            sample_mapping=json.dumps(analysis_design.sample_mapping, indent=2)
        ).to_messages()
        
        # Get the response
        response = self.llm(formatted_prompt)
        
        # Extract code from the response (handling markdown code blocks)
        code = response.content
        if "```" in code:
            code_blocks = re.findall(r"```(?:python|r)?\n(.*?)```", code, re.DOTALL)
            if code_blocks:
                code = "\n".join(code_blocks)
        
        # Update the analysis design with the code
        analysis_design.code_snippet = code
        
        return code

    def analyze_dataset(self, geo_accession: str) -> Tuple[DatasetMetadata, List[SampleMetadata], List[AnalysisDesign]]:
        """
        Full analysis of a GEO dataset
        Returns metadata, samples, and analysis designs
        """
        # Fetch the dataset
        gse = self.fetch_geo_dataset(geo_accession)
        if not gse:
            logger.error(f"Failed to fetch dataset {geo_accession}")
            return None, [], []
        
        # Parse metadata
        metadata = self.parse_dataset_metadata(gse)
        
        # Parse samples
        samples = self.parse_sample_metadata(gse)
        
        # Identify condition groups
        condition_groups = self.identify_condition_groups(samples)
        metadata.condition_groups = condition_groups
        
        # Create analysis designs
        analysis_designs = self.create_analysis_design(geo_accession, condition_groups)
        
        # Generate analysis code for each design
        for design in analysis_designs:
            self.generate_analysis_code(design)
        
        return metadata, samples, analysis_designs

# %% [insights_generator_agent]
class InsightsGeneratorAgent(BaseAgent):
    """Agent responsible for generating insights from analysis results"""
    name: str = "InsightsGenerator"
    description: str = "Generates insights from multiple dataset analysis results"
    
    def __init__(self, **data):
        super().__init__(**data)
        self.initialize_llm()
    
    def find_common_degs(self, analysis_results: List[AnalysisResult]) -> List[str]:
        """Find genes that are differentially expressed across multiple datasets"""
        if not analysis_results:
            return []
        
        # Extract gene lists from each result
        gene_sets = []
        for result in analysis_results:
            genes = [gene["gene_id"] for gene in result.top_degs if "gene_id" in gene]
            if genes:
                gene_sets.append(set(genes))
        
        # Find intersection if we have at least two sets
        if len(gene_sets) >= 2:
            common_genes = gene_sets[0]
            for gene_set in gene_sets[1:]:
                common_genes = common_genes.intersection(gene_set)
            return list(common_genes)
        elif len(gene_sets) == 1:
            return list(gene_sets[0])
        else:
            return []
    
    def generate_insights(
        self, 
        research_question: ResearchQuestion,
        datasets: List[DatasetMetadata],
        analysis_results: List[AnalysisResult]
    ) -> InsightCollection:
        """Generate insights from analysis results across multiple datasets"""
        logger.info(f"Generating insights for {len(analysis_results)} datasets")
        
        # Create a basic insight collection
        collection = InsightCollection(
            research_question=research_question,
            insights=[],
            common_degs=self.find_common_degs(analysis_results),
            dataset_summary={
                "total_datasets": len(datasets),
                "datasets_analyzed": len(analysis_results),
                "dataset_ids": [d.geo_accession for d in datasets]
            }
        )
        
        # If no analysis results, return empty collection
        if not analysis_results:
            return collection
        
        # Create a prompt for insight generation
        system_message = SystemMessagePromptTemplate.from_template(
            """You are an expert in bioinformatics and systems biology.
            Generate meaningful scientific insights by analyzing results from multiple RNA-seq datasets.
            
            Focus on:
            1. Common differentially expressed genes across datasets
            2. Shared biological pathways and processes
            3. Potential biomarkers or therapeutic targets
            4. Disease mechanisms suggested by the data
            
            Format each insight as:
            {
                "title": "Brief descriptive title",
                "description": "Detailed explanation of the insight",
                "supporting_datasets": ["dataset_ids"],
                "supporting_evidence": ["specific findings that support this insight"],
                "genes_involved": ["list of genes"],
                "confidence_score": float (0.0-1.0)
            }
            
            Return a JSON array of insights, with each insight following the above format.
            Generate 3-5 high-quality insights.
            """
        )
        
        human_message = HumanMessagePromptTemplate.from_template(
            """Research Question: {question}
            
            Datasets:
            {datasets_summary}
            
            Analysis Results:
            {results_summary}
            
            Common DEGs across datasets: {common_degs}
            
            Generate scientific insights based on these results.
            """
        )
        
        chat_prompt = ChatPromptTemplate.from_messages([system_message, human_message])
        
        # Format the datasets and results summaries
        datasets_summary = "\n".join([
            f"- {d.geo_accession}: {d.title} ({len(d.condition_groups)} condition groups)"
            for d in datasets
        ])
        
        results_summary = "\n".join([
            f"- {r.dataset_id} ({r.analysis_design.control_group} vs {r.analysis_design.experimental_group}): "
            f"{r.deg_count} DEGs, top genes: "
            f"{', '.join([g.get('gene_id', 'Unknown') for g in r.top_degs[:5]])}"
            for r in analysis_results
        ])
        
        formatted_prompt = chat_prompt.format_prompt(
            question=research_question.question,
            datasets_summary=datasets_summary,
            results_summary=results_summary,
            common_degs=", ".join(collection.common_degs[:20])  # Show up to 20 common DEGs
        ).to_messages()
        
        # Get the response
        response = self.llm(formatted_prompt)
        
        # Parse the insights
        try:
            # Extract JSON from the response
            json_content = response.content
            if "```json" in json_content:
                json_content = re.search(r"```json\n(.*?)\n```", json_content, re.DOTALL).group(1)
            elif "```" in json_content:
                json_content = re.search(r"```\n(.*?)\n```", json_content, re.DOTALL).group(1)
            
            insights_data = json.loads(json_content)
            
            # Convert to Insight objects
            for insight_data in insights_data:
                insight = Insight(**insight_data)
                collection.add_insight(insight)
            
        except (json.JSONDecodeError, AttributeError) as e:
            logger.error(f"Failed to parse insights: {e}")
            # Create a single basic insight
            fallback_insight = Insight(
                title="Potential disease biomarkers",
                description=f"Analysis of {len(analysis_results)} GEO datasets identified {len(collection.common_degs)} "
                            f"genes that are consistently differentially expressed in the context of "
                            f"{research_question.disease_contexts[0].name if research_question.disease_contexts else 'the studied condition'}.",
                supporting_datasets=[r.dataset_id for r in analysis_results],
                genes_involved=collection.common_degs[:10],
                confidence_score=0.6
            )
            collection.add_insight(fallback_insight)
        
        return collection

# %% [markdown]
# ## Workflow Orchestration
# 
# The workflow coordinator manages the execution of the end-to-end analysis process.

# %% [workflow_coordinator]
class UORCAWorkflowCoordinator:
    """Coordinates the end-to-end workflow for the UORCA system"""
    
    def __init__(self):
        """Initialize the workflow coordinator with all required agents"""
        self.research_question_agent = ResearchQuestionAgent()
        self.dataset_finder_agent = DatasetFinderAgent()
        self.dataset_analyzer_agent = DatasetAnalyzerAgent()
        self.insights_generator_agent = InsightsGeneratorAgent()
        
        self.config = config
        
        # Create output directory
        os.makedirs(self.config.output_dir, exist_ok=True)
    
    def run_workflow(self, question_text: str) -> InsightCollection:
        """Run the full workflow from question to insights"""
        # Step 1: Analyze the research question
        research_question = self.research_question_agent.analyze_question(question_text)
        
        # Step 2: Find relevant datasets
        datasets = self.dataset_finder_agent.get_top_relevant_datasets(
            research_question, 
            max_datasets=self.config.max_datasets
        )
        
        # Step 3: Analyze datasets
        dataset_results = []
        for dataset in datasets:
            metadata, samples, analysis_designs = self.dataset_analyzer_agent.analyze_dataset(dataset.geo_accession)
            
            # Store updated metadata
            if metadata:
                # Merge relevance information
                metadata.relevance_score = dataset.relevance_score
                metadata.relevance_reasoning = dataset.relevance_reasoning
                
                dataset_results.append({
                    "metadata": metadata,
                    "samples": samples,
                    "analysis_designs": analysis_designs
                })
        
        # Here we'd normally execute the analyses, but we'll simulate results for this prototype
        simulated_results = self.simulate_analysis_results(dataset_results)
        
        # Step 4: Generate insights
        insights = self.insights_generator_agent.generate_insights(
            research_question,
            [r["metadata"] for r in dataset_results],
            simulated_results
        )
        
        # Save results
        self.save_results(research_question, dataset_results, simulated_results, insights)
        
        return insights
    
    def simulate_analysis_results(self, dataset_results: List[Dict]) -> List[AnalysisResult]:
        """Simulate analysis results for demonstration purposes"""
        analysis_results = []
        
        for result in dataset_results:
            metadata = result["metadata"]
            analysis_designs = result["analysis_designs"]
            
            for design in analysis_designs:
                # Create a simulated result
                deg_count = np.random.randint(100, 1000)
                
                # Generate some fake top DEGs
                top_degs = []
                gene_symbols = ["BRCA1", "TP53", "EGFR", "KRAS", "PTEN", "AKT1", "PIK3CA", 
                               "MAPK1", "TNF", "IL6", "STAT3", "MYC", "VEGFA", "CDH1"]
                
                for i in range(min(10, deg_count)):
                    # Randomly sample a gene if we need more than our list
                    gene_symbol = gene_symbols[i] if i < len(gene_symbols) else f"GENE{i+1}"
                    
                    top_degs.append({
                        "gene_id": gene_symbol,
                        "log2FoldChange": np.random.uniform(-5, 5),
                        "padj": np.random.uniform(0, 0.05)
                    })
                
                analysis_result = AnalysisResult(
                    dataset_id=metadata.geo_accession,
                    analysis_design=design,
                    deg_count=deg_count,
                    top_degs=top_degs,
                    enriched_pathways=[],
                    execution_log="Simulated analysis execution",
                    success=True
                )
                
                analysis_results.append(analysis_result)
        
        return analysis_results
    
    def save_results(
        self, 
        research_question: ResearchQuestion,
        dataset_results: List[Dict],
        analysis_results: List[AnalysisResult],
        insights: InsightCollection
    ):
        """Save all results to output files"""
        output_dir = self.config.output_dir
        
        # Save research question
        with open(f"{output_dir}/research_question.json", "w") as f:
            f.write(research_question.json(indent=2))
        
        # Save datasets
        datasets = [r["metadata"] for r in dataset_results]
        with open(f"{output_dir}/datasets.json", "w") as f:
            datasets_json = [d.dict() for d in datasets]
            json.dump(datasets_json, f, indent=2)
        
        # Save analysis results
        with open(f"{output_dir}/analysis_results.json", "w") as f:
            results_json = [r.dict() for r in analysis_results]
            json.dump(results_json, f, indent=2)
        
        # Save insights
        with open(f"{output_dir}/insights.json", "w") as f:
            f.write(insights.json(indent=2))
        
        # Create a summary report
        summary = {
            "research_question": research_question.question,
            "datasets_analyzed": len(datasets),
            "total_degs": sum(r.deg_count for r in analysis_results),
            "common_degs": insights.common_degs,
            "insight_count": len(insights.insights),
            "top_insights": [
                {"title": i.title, "confidence": i.confidence_score}
                for i in sorted(insights.insights, key=lambda x: x.confidence_score, reverse=True)[:3]
            ]
        }
        
        with open(f"{output_dir}/summary.json", "w") as f:
            json.dump(summary, f, indent=2)

# %% [markdown]
# ## Test and Demonstration
# 
# The following cells demonstrate the workflow with test cases.

# %% [test_workflow]
def run_test_workflow():
    """Run a test workflow with a sample research question"""
    coordinator = UORCAWorkflowCoordinator()
    
    # Sample research questions
    test_questions = [
        "What genes are differentially expressed in lung cancer compared to normal tissue?",
        "What are the key pathways affected in Alzheimer's disease based on RNA-seq data?",
        "Identify potential biomarkers for breast cancer from public GEO datasets."
    ]
    
    # Run the workflow for the first question
    question = test_questions[0]
    print(f"Running test workflow for question: {question}")
    
    insights = coordinator.run_workflow(question)
    
    # Print a summary of results
    print("\nWorkflow completed. Summary:")
    print(f"- Analyzed {insights.dataset_summary['total_datasets']} datasets")
    print(f"- Found {len(insights.common_degs)} common differentially expressed genes")
    print(f"- Generated {len(insights.insights)} insights")
    
    print("\nTop insights:")
    for i, insight in enumerate(sorted(insights.insights, key=lambda x: x.confidence_score, reverse=True)):
        if i >= 3:
            break
        print(f"- {insight.title} (confidence: {insight.confidence_score:.2f})")
        print(f"  Genes: {', '.join(insight.genes_involved[:5])}")
    
    return insights

# %% [main]
if __name__ == "__main__":
    # Run the test workflow
    insights = run_test_workflow()
    
    # Display the results
    print("\nWorkflow execution complete. Results saved to:", config.output_dir)
