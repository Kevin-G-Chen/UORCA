"""
Main entry point for local testing of the UORCA workflow.
"""
import argparse
import logging
import sys
import json
from pathlib import Path

from UORCA.services.workflow_service import WorkflowService
from UORCA.config import get_settings

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler("uorca.log")
    ]
)

logger = logging.getLogger(__name__)

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="UORCA: Unified -omics reference corpus of analyses")
    
    # Main arguments
    parser.add_argument(
        "--research-question", 
        type=str,
        help="Research question to address"
    )
    parser.add_argument(
        "--workflow-id", 
        type=str,
        help="ID of an existing workflow to resume or analyze"
    )
    
    # Action arguments
    action_group = parser.add_mutually_exclusive_group()
    action_group.add_argument(
        "--create", 
        action="store_true",
        help="Create a new workflow"
    )
    action_group.add_argument(
        "--execute", 
        action="store_true",
        help="Execute the workflow"
    )
    action_group.add_argument(
        "--summarize", 
        action="store_true",
        help="Generate a summary of the workflow results"
    )
    action_group.add_argument(
        "--status", 
        action="store_true",
        help="Check the status of a workflow"
    )
    
    # Optional arguments
    parser.add_argument(
        "--output", 
        type=str,
        help="Output file for the results"
    )
    parser.add_argument(
        "--verbose", 
        action="store_true",
        help="Enable verbose logging"
    )
    
    return parser.parse_args()

def main():
    """Execute the UORCA workflow."""
    args = parse_args()
    
    # Set log level if verbose
    if args.verbose:
        logger.setLevel(logging.DEBUG)
        for handler in logger.handlers:
            handler.setLevel(logging.DEBUG)
    
    # Create the workflow service
    workflow_service = WorkflowService()
    
    # Handle the requested action
    try:
        if args.create:
            # Create a new workflow
            if not args.research_question:
                logger.error("Research question is required to create a workflow")
                sys.exit(1)
                
            workflow_state = workflow_service.create_workflow(args.research_question)
            logger.info(f"Created workflow with ID: {workflow_state.workflow_id}")
            
            # Save the workflow ID to a file for reference
            with open("workflow_id.txt", "w") as f:
                f.write(workflow_state.workflow_id)
            
            # Print the workflow state
            print(f"Workflow created with ID: {workflow_state.workflow_id}")
            print(f"Research Question: {workflow_state.research_question}")
            
        elif args.execute:
            # Execute a workflow
            workflow_state = None
            
            if args.workflow_id:
                # Load the workflow
                workflow_state = workflow_service.load_workflow_state(args.workflow_id)
                if not workflow_state:
                    logger.error(f"Workflow with ID {args.workflow_id} not found")
                    sys.exit(1)
            elif args.research_question:
                # Create a new workflow
                workflow_state = workflow_service.create_workflow(args.research_question)
                logger.info(f"Created workflow with ID: {workflow_state.workflow_id}")
            else:
                # Try to load from file
                try:
                    with open("workflow_id.txt", "r") as f:
                        workflow_id = f.read().strip()
                        workflow_state = workflow_service.load_workflow_state(workflow_id)
                        if not workflow_state:
                            logger.error(f"Workflow with ID {workflow_id} not found")
                            sys.exit(1)
                except FileNotFoundError:
                    logger.error("No workflow ID provided and no workflow_id.txt file found")
                    sys.exit(1)
            
            # Execute the workflow
            updated_state = workflow_service.execute_workflow(workflow_state)
            
            # Save the output
            if args.output:
                with open(args.output, "w") as f:
                    json.dump(updated_state.dict(), f, indent=2)
            
            # Print the workflow state
            print(f"Workflow {updated_state.workflow_id} executed")
            print(f"Final state: {updated_state.state}")
            print(f"Selected datasets: {', '.join(updated_state.selected_datasets)}")
            print(f"Analyzed datasets: {', '.join(updated_state.analysis_results.keys())}")
            print(f"Integrated insights: {len(updated_state.integrated_insights)}")
            
            if updated_state.errors:
                print(f"Errors: {len(updated_state.errors)}")
                for error in updated_state.errors:
                    print(f"  - {error}")
            
            if updated_state.warnings:
                print(f"Warnings: {len(updated_state.warnings)}")
                for warning in updated_state.warnings:
                    print(f"  - {warning}")
            
        elif args.summarize:
            # Generate a summary
            workflow_state = None
            
            if args.workflow_id:
                # Load the workflow
                workflow_state = workflow_service.load_workflow_state(args.workflow_id)
                if not workflow_state:
                    logger.error(f"Workflow with ID {args.workflow_id} not found")
                    sys.exit(1)
            else:
                # Try to load from file
                try:
                    with open("workflow_id.txt", "r") as f:
                        workflow_id = f.read().strip()
                        workflow_state = workflow_service.load_workflow_state(workflow_id)
                        if not workflow_state:
                            logger.error(f"Workflow with ID {workflow_id} not found")
                            sys.exit(1)
                except FileNotFoundError:
                    logger.error("No workflow ID provided and no workflow_id.txt file found")
                    sys.exit(1)
            
            # Generate the summary
            summary = workflow_service.get_workflow_summary(workflow_state)
            
            # Save the output
            if args.output:
                with open(args.output, "w") as f:
                    f.write(summary)
            else:
                print("\n" + "=" * 80)
                print("WORKFLOW SUMMARY")
                print("=" * 80)
                print(summary)
                print("=" * 80)
            
        elif args.status:
            # Check the status of a workflow
            workflow_state = None
            
            if args.workflow_id:
                # Load the workflow
                workflow_state = workflow_service.load_workflow_state(args.workflow_id)
                if not workflow_state:
                    logger.error(f"Workflow with ID {args.workflow_id} not found")
                    sys.exit(1)
            else:
                # Try to load from file
                try:
                    with open("workflow_id.txt", "r") as f:
                        workflow_id = f.read().strip()
                        workflow_state = workflow_service.load_workflow_state(workflow_id)
                        if not workflow_state:
                            logger.error(f"Workflow with ID {workflow_id} not found")
                            sys.exit(1)
                except FileNotFoundError:
                    logger.error("No workflow ID provided and no workflow_id.txt file found")
                    sys.exit(1)
            
            # Print the workflow state
            print(f"Workflow ID: {workflow_state.workflow_id}")
            print(f"Research Question: {workflow_state.research_question}")
            print(f"Current State: {workflow_state.state}")
            print(f"Search Queries: {len(workflow_state.search_queries)}")
            print(f"Search Results: {len(workflow_state.search_results)}")
            print(f"Dataset Relevance Scores: {len(workflow_state.dataset_relevance_scores)}")
            print(f"Selected Datasets: {len(workflow_state.selected_datasets)}")
            print(f"Fetched Datasets: {len(workflow_state.fetched_datasets)}")
            print(f"Analysis Requests: {len(workflow_state.analysis_requests)}")
            print(f"Analysis Results: {len(workflow_state.analysis_results)}")
            print(f"Integrated Insights: {len(workflow_state.integrated_insights)}")
            
            if workflow_state.errors:
                print(f"Errors: {len(workflow_state.errors)}")
                for error in workflow_state.errors:
                    print(f"  - {error}")
            
            if workflow_state.warnings:
                print(f"Warnings: {len(workflow_state.warnings)}")
                for warning in workflow_state.warnings:
                    print(f"  - {warning}")
            
        else:
            # No action specified
            logger.error("No action specified. Use --create, --execute, --summarize, or --status")
            sys.exit(1)
            
    except Exception as e:
        logger.error(f"Error executing UORCA: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
