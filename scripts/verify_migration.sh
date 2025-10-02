#!/bin/bash
set -e

echo "=== UORCA Migration Verification ==="
echo

# 1. Check for sys.path manipulation
echo "1. Checking for sys.path manipulation..."
if grep -r "sys.path.insert.*main_workflow" uorca/gui/ 2>/dev/null; then
    echo "❌ FAIL: Found sys.path manipulation in uorca/gui/"
    exit 1
else
    echo "✅ PASS: No sys.path manipulation found"
fi

# 2. Check legacy files are removed
echo
echo "2. Checking legacy files are removed..."
LEGACY_FILES=(
    "main_workflow/reporting/ResultsIntegration.py"
    "main_workflow/reporting/single_analysis_plots.py"
    "main_workflow/reporting/ortholog_mapper.py"
    "main_workflow/reporting/ai_agent_factory.py"
    "main_workflow/reporting/config_loader.py"
    "main_workflow/reporting/mcp_server_core.py"
)

for file in "${LEGACY_FILES[@]}"; do
    if [ -f "$file" ]; then
        echo "❌ FAIL: Legacy file still exists: $file"
        exit 1
    fi
done
echo "✅ PASS: All legacy files removed"

# 3. Test imports
echo
echo "3. Testing imports..."
uv run python -c "from uorca.gui.results_integration import ResultsIntegrator; print('  ✅ ResultsIntegrator')"
uv run python -c "from uorca.gui.single_analysis_plots import create_pca_plot; print('  ✅ single_analysis_plots')"
uv run python -c "from uorca.gui.ortholog_mapper import expand_genes_all_vs_all; print('  ✅ ortholog_mapper')"
uv run python -c "from uorca.gui.ai import create_uorca_agent; print('  ✅ AI components')"
uv run python -c "from uorca.gui.mcp import get_filtered_dataframe; print('  ✅ MCP server')"

# 4. Type checking
echo
echo "4. Running type checks..."
uv run pyright uorca/gui/ --level warning 2>&1 | head -20

# 5. Test Streamlit app startup (just import, don't run)
echo
echo "5. Testing Streamlit app startup..."
uv run python -c "import sys; sys.argv = ['test']; import uorca.gui.uorca_explorer; print('  ✅ Streamlit app imports successfully')"

echo
echo "=== ✅ All Verification Checks Passed ==="
