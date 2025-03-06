"""
Service for interacting with LLM APIs.
"""
import logging
import json
import requests
from typing import Dict, List, Any, Optional, Union
from pydantic import BaseModel

from ..config import get_settings
from ..models.geo_models import GEODatasetWithExpression, OmicsType
from ..models.workflow_models import DatasetRelevanceScore
from ..models.analysis_models import AnalysisType, AnalysisParameter, Insight, InsightType

logger = logging.getLogger(__name__)

class Message(BaseModel):
    """Message in a conversation with an LLM."""
    role: str
    content: str


class LLMService:
    """Service for interacting with LLM APIs."""
    
    def __init__(self):
        """Initialize the LLM service."""
        self.settings = get_settings()
        self.api_key = self.settings.llm_api_key
        self.api_url = self.settings.llm_api_url
        self.model = self.settings.llm_model
    
    def generate_search_queries(self, research_question: str, max_queries: int = 3) -> List[Dict[str, str]]:
        """
        Generate search queries for a research question.
        
        Args:
            research_question: The research question to generate queries for
            max_queries: Maximum number of queries to generate
            
        Returns:
            List of dictionaries with 'term' and 'description' for each query
        """
        prompt = self._create_prompt(
            system_message="""You are a bioinformatics expert helping to generate search queries for the NCBI GEO database.
            Given a research question, generate specific search queries that will help find relevant datasets.
            Focus on key biological terms, conditions, diseases, tissue types, and omics technologies mentioned in the question.
            Each search query should be targeted and specific to maximize the chance of finding relevant datasets.""",
            
            user_message=f"""Research Question: {research_question}

            Please generate up to {max_queries} search queries for the NCBI GEO database that would help find datasets 
            relevant to this research question. For each query, provide a brief description of what aspect of the 
            research question it addresses.
            
            Format your response as a JSON array with each query as an object containing 'term' and 'description' fields.
            Example:
            [
                {{
                    "term": "your search term here",
                    "description": "explanation of this search strategy"
                }}
            ]
            """
        )
        
        response = self._call_llm_api(prompt)
        
        # Extract JSON from the response
        try:
            # Find JSON-like content in the response
            start_idx = response.find('[')
            end_idx = response.rfind(']') + 1
            
            if start_idx >= 0 and end_idx > start_idx:
                json_content = response[start_idx:end_idx]
                queries = json.loads(json_content)
                
                # Ensure we don't exceed the max number of queries
                return queries[:max_queries]
            else:
                logger.error("Failed to extract JSON from LLM response")
                return []
        except Exception as e:
            logger.error(f"Error parsing search queries from LLM response: {str(e)}")
            return []
    
    def evaluate_dataset_relevance(
        self, 
        research_question: str, 
        dataset_metadata: Dict[str, Any]
    ) -> DatasetRelevanceScore:
        """
        Evaluate the relevance of a dataset to a research question.
        
        Args:
            research_question: The research question
            dataset_metadata: Metadata for the dataset
            
        Returns:
            DatasetRelevanceScore with relevance score and justification
        """
        # Create a simplified version of the metadata for the prompt
        metadata_str = json.dumps({
            "id": dataset_metadata.get("id", ""),
            "title": dataset_metadata.get("title", ""),
            "summary": dataset_metadata.get("summary", ""),
            "organism": dataset_metadata.get("organism", ""),
            "omics_type": dataset_metadata.get("omics_type", ""),
            # Add other relevant fields
        }, indent=2)
        
        prompt = self._create_prompt(
            system_message="""You are a bioinformatics expert evaluating the relevance of genomic datasets to a specific research question.
            Assess how well the dataset matches the research question, considering factors like:
            1. Organism relevance
            2. Disease/condition match
            3. Tissue type appropriateness
            4. Experimental design
            5. Omics technology suitability
            6. Sample size and quality
            
            Provide a relevance score from 0.0 to 1.0 where:
            - 0.0-0.2: Not relevant
            - 0.3-0.5: Somewhat relevant
            - 0.6-0.8: Relevant
            - 0.9-1.0: Highly relevant
            
            Also provide a detailed justification explaining your reasoning.""",
            
            user_message=f"""Research Question: {research_question}

            Dataset Metadata:
            {metadata_str}
            
            Evaluate how relevant this dataset is to the research question.
            Provide your response as a JSON object with the following fields:
            - relevance_score: A float between 0.0 and 1.0
            - justification: A string explaining your reasoning
            
            Example:
            {{
                "relevance_score": 0.75,
                "justification": "This dataset is relevant because..."
            }}
            """
        )
        
        response = self._call_llm_api(prompt)
        
        # Extract JSON from the response
        try:
            # Find JSON-like content in the response
            start_idx = response.find('{')
            end_idx = response.rfind('}') + 1
            
            if start_idx >= 0 and end_idx > start_idx:
                json_content = response[start_idx:end_idx]
                evaluation = json.loads(json_content)
                
                # Create the relevance score object
                return DatasetRelevanceScore(
                    dataset_id=dataset_metadata.get("id", ""),
                    relevance_score=float(evaluation.get("relevance_score", 0.0)),
                    justification=evaluation.get("justification", "")
                )
            else:
                logger.error("Failed to extract JSON from LLM response")
                return DatasetRelevanceScore(
                    dataset_id=dataset_metadata.get("id", ""),
                    relevance_score=0.0,
                    justification="Failed to evaluate relevance"
                )
        except Exception as e:
            logger.error(f"Error parsing relevance score from LLM response: {str(e)}")
            return DatasetRelevanceScore(
                dataset_id=dataset_metadata.get("id", ""),
                relevance_score=0.0,
                justification=f"Error evaluating relevance: {str(e)}"
            )
    
    def design_analysis(
        self,
        research_question: str,
        dataset: GEODatasetWithExpression
    ) -> Dict[str, Any]:
        """
        Design an appropriate analysis for a dataset based on the research question.
        
        Args:
            research_question: The research question
            dataset: The dataset with expression data
            
        Returns:
            Dictionary with analysis type and parameters
        """
        # Simplify dataset for the prompt
        metadata = dataset.metadata
        metadata_dict = {
            "id": metadata.id,
            "title": metadata.title,
            "summary": metadata.summary[:500] + "..." if len(metadata.summary) > 500 else metadata.summary,
            "organism": metadata.organism,
            "omics_type": str(metadata.omics_type),
            "sample_count": len(dataset.sample_metadata),
        }
        
        # Add sample groups information if available
        sample_groups = {}
        for sample_id, metadata in dataset.sample_metadata.items():
            for key, value in metadata.items():
                if key.lower() in ["tissue", "disease", "condition", "treatment", "group"]:
                    if key not in sample_groups:
                        sample_groups[key] = set()
                    sample_groups[key].add(value)
        
        # Convert sets to lists for JSON serialization
        for key in sample_groups:
            sample_groups[key] = list(sample_groups[key])
        
        metadata_dict["sample_groups"] = sample_groups
        
        # Information about expression matrix if available
        if dataset.expression_matrix:
            metadata_dict["expression_matrix"] = {
                "gene_count": len(dataset.expression_matrix.gene_ids),
                "sample_count": len(dataset.expression_matrix.sample_ids),
                "has_values": len(dataset.expression_matrix.values) > 0
            }
        
        metadata_str = json.dumps(metadata_dict, indent=2)
        
        prompt = self._create_prompt(
            system_message="""You are a bioinformatics expert designing analyses for omics datasets.
            Given a research question and dataset metadata, determine the most appropriate analysis approach.
            Consider the dataset type, available samples, and the specific biological question being asked.
            Design an analysis that will extract valuable insights relevant to the research question.""",
            
            user_message=f"""Research Question: {research_question}

            Dataset Metadata:
            {metadata_str}
            
            Design an appropriate analysis for this dataset that would help answer the research question.
            Consider:
            1. What type of analysis is most suitable (e.g., differential expression, clustering, pathway enrichment)?
            2. What parameters or settings should be used?
            3. How should samples be grouped or compared?
            
            Provide your response as a JSON object with the following fields:
            - analysis_type: One of ["differential_expression", "pathway_enrichment", "clustering", "correlation", "classification", "custom"]
            - description: A brief description of the analysis approach
            - parameters: An array of parameter objects, each with "name", "value", and "description" fields
            - code_template: A Python code snippet template for performing this analysis
            
            Example:
            {{
                "analysis_type": "differential_expression",
                "description": "Compare gene expression between tumor and normal samples",
                "parameters": [
                    {{
                        "name": "group_var",
                        "value": "condition",
                        "description": "Variable to use for grouping samples"
                    }},
                    {{
                        "name": "test_method",
                        "value": "DESeq2",
                        "description": "Statistical method for differential testing"
                    }}
                ],
                "code_template": "# Python code template for the analysis..."
            }}
            """
        )
        
        response = self._call_llm_api(prompt)
        
        # Extract JSON from the response
        try:
            # Find JSON-like content in the response
            start_idx = response.find('{')
            end_idx = response.rfind('}') + 1
            
            if start_idx >= 0 and end_idx > start_idx:
                json_content = response[start_idx:end_idx]
                analysis_design = json.loads(json_content)
                return analysis_design
            else:
                logger.error("Failed to extract JSON from LLM response")
                return {
                    "analysis_type": "custom",
                    "description": "Failed to design analysis",
                    "parameters": [],
                    "code_template": "# Failed to generate code template"
                }
        except Exception as e:
            logger.error(f"Error parsing analysis design from LLM response: {str(e)}")
            return {
                "analysis_type": "custom",
                "description": f"Error designing analysis: {str(e)}",
                "parameters": [],
                "code_template": "# Error generating code template"
            }
    
    def generate_analysis_code(
        self,
        research_question: str,
        dataset: GEODatasetWithExpression,
        analysis_design: Dict[str, Any]
    ) -> str:
        """
        Generate Python code to perform the designed analysis on the dataset.
        
        Args:
            research_question: The research question
            dataset: The dataset with expression data
            analysis_design: The analysis design from design_analysis()
            
        Returns:
            Python code string
        """
        # Extract relevant info
        dataset_id = dataset.metadata.id
        analysis_type = analysis_design.get("analysis_type", "custom")
        parameters = analysis_design.get("parameters", [])
        code_template = analysis_design.get("code_template", "# No code template available")
        
        # Build a prompt with more context
        prompt = self._create_prompt(
            system_message="""You are a bioinformatics expert generating Python code for analyzing omics datasets.
            Given a research question, dataset details, and analysis design, generate complete, executable Python code.
            Your code should be well-documented, robust, and include appropriate error handling.
            Use standard bioinformatics libraries like pandas, numpy, scipy, sklearn, and specialized libraries as needed.""",
            
            user_message=f"""Research Question: {research_question}

            Dataset ID: {dataset_id}
            Analysis Type: {analysis_type}
            
            Analysis Parameters:
            {json.dumps(parameters, indent=2)}
            
            Code Template:
            {code_template}
            
            Generate complete, well-documented Python code to perform this analysis on the dataset.
            Include functions for loading the data, preprocessing, performing the analysis, and visualizing the results.
            Make sure the code is properly structured, includes error handling, and is ready to execute.
            """
        )
        
        response = self._call_llm_api(prompt)
        
        # Extract code from the response
        try:
            # Find code blocks in markdown format
            code_start = response.find("```python")
            code_end = response.rfind("```")
            
            if code_start >= 0 and code_end > code_start:
                # Extract the code without the markdown indicators
                code_start += 9  # Length of "```python" + newline
                code = response[code_start:code_end].strip()
                return code
            else:
                logger.warning("No code block found in response, returning the entire response")
                return response
        except Exception as e:
            logger.error(f"Error extracting code from LLM response: {str(e)}")
            return f"# Error generating code: {str(e)}\n\n{response}"
    
    def integrate_insights(
        self,
        research_question: str,
        analysis_results: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        """
        Integrate insights from multiple dataset analyses.
        
        Args:
            research_question: The research question
            analysis_results: Results from multiple dataset analyses
            
        Returns:
            List of integrated insights
        """
        # Simplify the analysis results for the prompt
        simplified_results = []
        
        for result in analysis_results:
            simplified_result = {
                "dataset_id": result.get("dataset_id", "unknown"),
                "dataset_title": result.get("dataset_title", "unknown"),
                "analysis_type": result.get("analysis_type", "unknown"),
                "significant_findings": result.get("significant_findings", [])
            }
            simplified_results.append(simplified_result)
        
        results_str = json.dumps(simplified_results, indent=2)
        
        prompt = self._create_prompt(
            system_message="""You are a bioinformatics expert integrating insights across multiple omics analyses.
            Given a research question and results from analyses of different datasets, identify common patterns, 
            contradictions, and novel biological trends.
            Focus on extracting meaningful insights that address the research question and could lead to testable hypotheses.
            Consider how findings from different datasets complement or contradict each other.""",
            
            user_message=f"""Research Question: {research_question}

            Analysis Results:
            {results_str}
            
            Integrate these findings to derive meaningful insights that address the research question.
            Identify:
            1. Common biological patterns or trends across datasets
            2. Key genes/pathways that appear in multiple analyses
            3. Potential contradictions or dataset-specific findings
            4. Novel hypotheses that could be tested experimentally
            
            Provide your insights as a JSON array of insight objects, each with the following fields:
            - insight_type: One of ["key_gene", "pathway", "biomarker", "drug_target", "disease_mechanism", "correlation", "contradiction", "other"]
            - description: A detailed description of the insight
            - supporting_datasets: Array of dataset IDs supporting this insight
            - confidence_score: A float between 0.0 and 1.0 indicating confidence
            - related_entities: Array of related genes, proteins, pathways, etc.
            
            Example:
            [
                {{
                    "insight_type": "key_gene",
                    "description": "Gene XYZ is consistently upregulated in all tumor samples across datasets",
                    "supporting_datasets": ["GDS1234", "GSE5678"],
                    "confidence_score": 0.85,
                    "related_entities": ["XYZ", "PATHWAY_ABC"]
                }}
            ]
            """
        )
        
        response = self._call_llm_api(prompt)
        
        # Extract JSON from the response
        try:
            # Find JSON-like content in the response
            start_idx = response.find('[')
            end_idx = response.rfind(']') + 1
            
            if start_idx >= 0 and end_idx > start_idx:
                json_content = response[start_idx:end_idx]
                insights = json.loads(json_content)
                return insights
            else:
                logger.error("Failed to extract JSON from LLM response")
                return []
        except Exception as e:
            logger.error(f"Error parsing insights from LLM response: {str(e)}")
            return []
    
    def generate_summary_report(
        self,
        research_question: str,
        datasets_analyzed: List[str],
        insights: List[Dict[str, Any]]
    ) -> str:
        """
        Generate a human-readable summary report of the analysis.
        
        Args:
            research_question: The research question
            datasets_analyzed: List of dataset IDs that were analyzed
            insights: List of integrated insights
            
        Returns:
            Formatted summary report
        """
        insights_str = json.dumps(insights, indent=2)
        
        prompt = self._create_prompt(
            system_message="""You are a bioinformatics expert creating summary reports of multi-omics analyses.
            Given a research question, analyzed datasets, and integrated insights, create a clear, well-structured summary.
            The summary should be understandable by researchers without detailed bioinformatics knowledge.
            Focus on the biological significance of the findings and their relevance to the research question.""",
            
            user_message=f"""Research Question: {research_question}

            Datasets Analyzed: {', '.join(datasets_analyzed)}
            
            Integrated Insights:
            {insights_str}
            
            Generate a comprehensive yet concise summary report of the analysis.
            The report should:
            1. Restate the research question and its context
            2. Summarize the datasets that were analyzed
            3. Present the key findings and insights, with emphasis on their biological significance
            4. Discuss potential implications for further research or clinical applications
            5. Note any limitations or areas requiring further investigation
            
            Format the report with clear section headings and well-structured paragraphs.
            """
        )
        
        response = self._call_llm_api(prompt)
        return response
    
    def _create_prompt(self, system_message: str, user_message: str) -> List[Message]:
        """
        Create a prompt with system and user messages.
        
        Args:
            system_message: The system message
            user_message: The user message
            
        Returns:
            List of Message objects
        """
        return [
            Message(role="system", content=system_message),
            Message(role="user", content=user_message)
        ]
    
    def _call_llm_api(self, messages: List[Message]) -> str:
        """
        Call the LLM API with the given messages.
        
        Args:
            messages: List of Message objects
            
        Returns:
            The model's response text
        """
        if not self.api_key:
            logger.warning("No API key provided for LLM service, using mock response")
            return "This is a mock response since no API key was provided."
        
        # Format the messages as expected by the API
        formatted_messages = [{"role": msg.role, "content": msg.content} for msg in messages]
        
        # Prepare the API request
        headers = {
            "Content-Type": "application/json",
            "Authorization": f"Bearer {self.api_key}"
        }
        
        payload = {
            "model": self.model,
            "messages": formatted_messages,
            "temperature": 0.3,  # Lower temperature for more deterministic responses
            "max_tokens": 2000
        }
        
        try:
            logger.debug(f"Calling LLM API with model: {self.model}")
            response = requests.post(self.api_url, headers=headers, json=payload)
            response.raise_for_status()
            
            # Parse the response
            result = response.json()
            return result.get("choices", [{}])[0].get("message", {}).get("content", "")
            
        except requests.RequestException as e:
            logger.error(f"Error calling LLM API: {str(e)}")
            return f"Error calling LLM API: {str(e)}"