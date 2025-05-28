#!/usr/bin/env python3
"""
UORCA Explorer Launcher with Integrated AI Landing Page

Simple script to launch UORCA Explorer with proper configuration and dependency checking.
The AI Landing Page provides automatic contrast selection and interpretive narratives.

Usage:
    python launch_uorca_explorer.py [--port PORT] [--help]
"""

import os
import sys
import subprocess
import argparse
from pathlib import Path

# Load environment variables from .env file
try:
    from dotenv import load_dotenv
    # Try loading from current directory first, then parent directories
    load_dotenv()
    # Also try loading from project root
    project_root = Path(__file__).parent.parent.parent
    env_file = project_root / ".env"
    if env_file.exists():
        load_dotenv(env_file)
except ImportError:
    # dotenv not available, continue without it
    pass

def print_banner():
    """Print welcome banner."""
    print("🧬 UORCA Explorer with AI Landing Page")
    print("=" * 50)
    print("Interactive RNA-seq Results Analysis with AI-Assisted Interpretation")
    print()

def check_uv_available():
    """Check if uv is available."""
    try:
        result = subprocess.run(["uv", "--version"], capture_output=True, text=True)
        if result.returncode == 0:
            print(f"✅ UV package manager: {result.stdout.strip()}")
            return True
        else:
            print("❌ UV not found or not working")
            return False
    except FileNotFoundError:
        print("❌ UV package manager not found")
        print("   Install UV: https://docs.astral.sh/uv/getting-started/installation/")
        return False

def check_ai_status():
    """Check AI features availability."""
    print("\n🤖 AI Features Status:")
    
    # Check if running with uv (packages available)
    try:
        result = subprocess.run([
            "uv", "run", "python", "-c", 
            "import openai; print('OpenAI available'); import os; print('API key:', 'SET' if os.getenv('OPENAI_API_KEY') else 'NOT SET')"
        ], capture_output=True, text=True, cwd=Path(__file__).parent.parent.parent)
        
        if result.returncode == 0:
            lines = result.stdout.strip().split('\n')
            if "OpenAI available" in lines[0]:
                print("  ✅ OpenAI package: Available")
                if "API key: SET" in lines[1]:
                    print("  ✅ API Key: Configured")
                    return "full"
                else:
                    print("  ⚠️ API Key: Not set (heuristic mode available)")
                    return "heuristic"
            else:
                print("  ❌ OpenAI package: Error")
                return "none"
        else:
            print("  ❌ Package check failed")
            return "none"
            
    except Exception as e:
        print(f"  ❌ Error checking AI status: {e}")
        return "none"

def get_default_results_dir():
    """Find a reasonable default results directory."""
    candidates = [
        "UORCA_results",
        "../UORCA_results", 
        "../../UORCA_results",
        os.path.expanduser("~/UORCA_results")
    ]
    
    for candidate in candidates:
        if os.path.exists(candidate) and os.path.isdir(candidate):
            return os.path.abspath(candidate)
    
    return os.getcwd()

def show_ai_features(ai_status):
    """Show what AI features are available."""
    print(f"\n✨ Available Features:")
    
    if ai_status == "full":
        print("  🎯 Full AI Landing Page")
        print("    • Smart contrast selection with biological relevance scoring")
        print("    • AI-optimized statistical thresholds")
        print("    • Intelligent narrative generation") 
        print("    • Biological interpretation and insights")
    elif ai_status == "heuristic":
        print("  🔧 Basic AI Landing Page (Heuristic Mode)")
        print("    • Automatic contrast selection based on DEG counts")
        print("    • Statistical threshold optimization")
        print("    • Template-based narrative generation")
        print("  • To enable full AI:")
        print("    - Environment: export OPENAI_API_KEY=your_key")
        print("    - .env file: Add OPENAI_API_KEY=your_key to .env")
    else:
        print("  📊 Manual Analysis Tools")
        print("    • Interactive dataset and contrast selection")
        print("    • Expression heatmaps and plots") 
        print("    • Quality control visualizations")
        print("    • To enable AI: uv add openai && export OPENAI_API_KEY=your_key")
    
    print("  📈 Standard Features (Always Available)")
    print("    • Interactive expression heatmaps with clustering")
    print("    • Gene expression violin plots across samples")
    print("    • RNA-seq QC plots (MDS, normalization, etc.)")
    print("    • Dataset and contrast browsing")
    print("    • Export capabilities for further analysis")

def launch_streamlit(port=8501):
    """Launch the Streamlit app."""
    print(f"\n🚀 Launching UORCA Explorer on port {port}...")
    
    # Set environment variable for default directory
    default_dir = get_default_results_dir()
    env = os.environ.copy()
    env["UORCA_DEFAULT_RESULTS_DIR"] = default_dir
    
    print(f"📁 Default results directory: {default_dir}")
    print(f"🌐 App will be available at: http://localhost:{port}")
    print(f"💡 The AI Landing Page is the first tab in the interface")
    print()
    print("=" * 50)
    print("🎉 Starting UORCA Explorer...")
    print("   Press Ctrl+C to stop the server")
    print("=" * 50)
    
    # Launch streamlit with uv run
    script_path = Path(__file__).parent / "uorca_explorer.py"
    
    try:
        subprocess.run([
            "uv", "run", "streamlit", "run", str(script_path),
            "--server.port", str(port),
            "--server.headless", "false"
        ], env=env, cwd=Path(__file__).parent.parent.parent)
    except KeyboardInterrupt:
        print("\n👋 UORCA Explorer stopped")
    except FileNotFoundError:
        print("❌ Could not find uorca_explorer.py")
        print("   Make sure you're running this script from the correct directory")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error launching Streamlit: {e}")
        sys.exit(1)

def main():
    """Main launcher function."""
    parser = argparse.ArgumentParser(
        description="Launch UORCA Explorer with AI Landing Page",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--port", "-p",
        type=int, 
        default=8501,
        help="Port for Streamlit server"
    )
    parser.add_argument(
        "--check-only", "-c",
        action="store_true",
        help="Only check dependencies, don't launch"
    )
    
    args = parser.parse_args()
    
    print_banner()
    
    # Check UV availability
    if not check_uv_available():
        print("\n❌ Cannot proceed without UV package manager")
        sys.exit(1)
    
    # Check AI status
    ai_status = check_ai_status()
    
    # Show available features
    show_ai_features(ai_status)
    
    if args.check_only:
        print(f"\n✅ Dependency check complete")
        if ai_status == "full":
            print("🎉 Ready to launch with full AI features!")
        elif ai_status == "heuristic":
            print("⚠️ Ready to launch with basic AI features")
        else:
            print("📊 Ready to launch with manual analysis tools")
        return
    
    # Launch the app
    try:
        launch_streamlit(args.port)
    except KeyboardInterrupt:
        print("\n👋 Launch cancelled")
        sys.exit(0)

if __name__ == "__main__":
    main()