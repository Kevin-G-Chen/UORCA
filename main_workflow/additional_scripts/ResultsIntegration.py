# Dependencies:
# - Base: os, glob, argparse, pandas, logging, pathlib
# - Visualization: panel, bokeh - install with: pip install panel bokeh
#
# New dependencies for improved visualization:
# Panel provides a high-level API for building interactive web applications
# Bokeh provides the underlying visualization capabilities
import os
import glob
import argparse
import pandas as pd
import logging
from pathlib import Path
import panel as pn
import bokeh
from bokeh.models import TableColumn, DataTable, ColumnDataSource, CustomJS

# Set up basic logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger(__name__)

def write_datatables_html(df: pd.DataFrame, path: Path, title="Integration Matrix"):
    # 1  ------- constants for CDN assets -------
    DATATABLES_CSS  = "https://cdn.datatables.net/2.0.7/css/dataTables.dataTables.min.css"
    BUTTONS_CSS     = "https://cdn.datatables.net/buttons/2.4.1/css/buttons.dataTables.min.css"

    JQUERY_JS       = "https://code.jquery.com/jquery-3.7.1.min.js"
    DATATABLES_JS   = "https://cdn.datatables.net/2.0.7/js/dataTables.min.js"
    BUTTONS_JS      = "https://cdn.datatables.net/buttons/2.4.1/js/dataTables.buttons.min.js"
    COLVIS_JS       = "https://cdn.datatables.net/buttons/2.4.1/js/buttons.colVis.min.js"

    # 2  ------- give the index a label so it has a header cell -------
    df = df.copy()
    df.index.name = "Gene"

    # 3  ------- build the HTML document -------
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{title}</title>

  <!-- DataTables+Buttons CSS -->
  <link rel="stylesheet" href="{DATATABLES_CSS}">
  <link rel="stylesheet" href="{BUTTONS_CSS}">

  <!-- jQuery and DataTables JS -->
  <script src="{JQUERY_JS}"></script>
  <script src="{DATATABLES_JS}"></script>
  <script src="{BUTTONS_JS}"></script>
  <script src="{COLVIS_JS}"></script>

  <!-- Tiny page styling -->
  <style>
    body {{ font-family: sans-serif; margin: 2rem; }}
    .toolbar {{ margin-bottom: 1rem; }}
    .toolbar select, .toolbar input {{ margin-right: 0.5rem; }}
  </style>
</head>
<body>

<h1>{title}</h1>

<div class="toolbar">
  <!-- (a) contrast-column selector -->
  <label for="contrastFilter">Show rows with 1 in&nbsp;</label>
  <select id="contrastFilter">
    <option value="">— any contrast —</option>
    {"".join(f'<option value="{c}">{c}</option>' for c in df.columns if c != "Num_DE_Contrasts")}
  </select>

  <!-- (b) Num_DE_Contrasts threshold -->
  <label for="numComparator">and&nbsp;Num_DE_Contrasts&nbsp;</label>
  <select id="numComparator">
    <option value=">=">≥</option>
    <option value="<=">≤</option>
  </select>
  <input type="number" id="numThreshold" value="1" min="0" style="width:6rem;">
</div>

{df.to_html(index=True, classes="display", table_id="matrix")}

<script>
$(function () {{

  // ------------- initialise the DataTable -------------
  var table = $('#matrix').DataTable({{
      dom: 'Bfrtip',                           // Buttons + filter + table
      buttons: [ 'colvis' ],                   // column visibility menu
      pageLength: 25,
      scrollX: true,
      deferRender: true,                       // build rows on demand
      scrollCollapse: true,                    // trim unused space
      scroller: true                           // adds virtual scrolling plugin
  }});

  // ------------- custom row-filter function -----------
  $.fn.dataTable.ext.search.push(function(settings, data, dataIndex) {{
      // (1) contrast filter
      var chosenCol = $('#contrastFilter').val();     // e.g. "PHOX2B_Knockdown_vs_Control"
      if (chosenCol) {{
          var colIdx = table.column(chosenCol + ':name').index(); // find its index
          if (data[colIdx] !== '1') return false;     // keep only rows with 1
      }}

      // (2) Num_DE_Contrasts threshold
      var cmp = $('#numComparator').val();            // ">=" or "<="
      var thresh = parseInt($('#numThreshold').val(), 10) || 0;
      var numIdx = table.column('Num_DE_Contrasts:name').index();
      var num = parseInt(data[numIdx], 10) || 0;
      return (cmp === '>=' ? num >= thresh : num <= thresh);
  }});

  // ------------- re-draw when controls change ---------
  $('#contrastFilter, #numComparator, #numThreshold').on('change keyup', function() {{
      table.draw();
  }});

}});
</script>

</body>
</html>"""
    path.write_text(html, encoding="utf-8")

def write_panel_bokeh_html(df: pd.DataFrame, path: Path, title="Integration Matrix"):
    """
    Generate interactive HTML visualization using Panel and Bokeh.
    
    This function creates a Tabulator widget with custom filtering capabilities
    and exports it as a standalone HTML file.
    
    Parameters:
    -----------
    df : pd.DataFrame
        The DataFrame to visualize
    path : Path
        Where to save the HTML file
    title : str
        Title for the visualization
    """
    # Initialize Panel extensions
    pn.extension('tabulator')
    
    # Set index name
    df = df.copy()
    df.index.name = "Gene"
    
    # Reset index to make Gene a column for easier interaction
    df_reset = df.reset_index()
    
    # Create widgets for filtering
    contrast_columns = [col for col in df.columns if col != "Num_DE_Contrasts"]
    contrast_select = pn.widgets.Select(
        name='Show rows with 1 in',
        options=['— any contrast —'] + contrast_columns,
        value='— any contrast —'
    )
    
    comparator_select = pn.widgets.Select(
        name='Num_DE_Contrasts',
        options=['≥', '≤'],
        value='≥'
    )
    
    threshold_input = pn.widgets.IntInput(
        name='Threshold',
        value=1,
        step=1,
        start=0
    )
    
    # Function to filter the data based on widget values
    def filter_data(contrast, comparator, threshold):
        filtered_df = df_reset.copy()
        
        # Filter by contrast if selected
        if contrast != '— any contrast —':
            filtered_df = filtered_df[filtered_df[contrast] == 1]
        
        # Filter by Num_DE_Contrasts
        if comparator == '≥':
            filtered_df = filtered_df[filtered_df['Num_DE_Contrasts'] >= threshold]
        else:
            filtered_df = filtered_df[filtered_df['Num_DE_Contrasts'] <= threshold]
            
        return filtered_df
    
    # Create interactive table with filtered data
    table = pn.widgets.Tabulator(
        pn.bind(filter_data, contrast_select, comparator_select, threshold_input),
        pagination='remote',
        page_size=25,
        height=600,
        layout='fit_data_fill',
        disabled=True  # Makes the table read-only
    )
    
    # Create the layout with title, filters, and table
    layout = pn.Column(
        f"# {title}",
        pn.Row(contrast_select, comparator_select, threshold_input),
        table
    )
    
    # Save to standalone HTML
    layout.save(str(path))
    
    return f"Panel+Bokeh visualization saved to {path}"

def find_contrast_csvs(root_dir: str) -> list[str]:
    """Recursively find all contrasts.csv under root_dir."""
    pattern = os.path.join(root_dir, "**", "contrasts.csv")
    matches = glob.glob(pattern, recursive=True)
    logger.info(f"Found {len(matches)} contrast CSV files")
    return matches

def read_contrast_names(contrast_csv: str) -> list[str]:
    """Read the 'name' column from a contrasts.csv."""
    df = pd.read_csv(contrast_csv)
    names = df["name"].dropna().astype(str).tolist()
    logger.debug(f"Read {len(names)} contrast names from {contrast_csv}")
    return names

def find_deg_file(root_dir: str, contrast_name: str) -> str | None:
    """Find a file named deg_{contrast_name}.csv anywhere under root_dir."""
    pattern = os.path.join(root_dir, "**", f"deg_{contrast_name}.csv")
    matches = glob.glob(pattern, recursive=True)
    return matches[0] if matches else None

def build_binary_matrix(root_dir: str, lfc_threshold: float = 0, padj_threshold: float = 1.0) -> pd.DataFrame:
    logger.info(f"Building binary matrix with |logFC| >= {lfc_threshold} and adj.P.Val <= {padj_threshold}")

    # 1) Collect all contrast names
    contrast_csvs = find_contrast_csvs(root_dir)
    contrast_names: list[str] = []
    for cfile in contrast_csvs:
        contrast_names += read_contrast_names(cfile)
    contrast_names = sorted(set(contrast_names))
    logger.info(f"Found {len(contrast_names)} unique contrast names")

    # 2) For each contrast, load the DEG file and record its genes that pass thresholds
    gene_sets: dict[str, set[str]] = {}
    all_genes: set[str] = set()
    total_genes_before_filtering = 0
    found_deg_files = 0

    for name in contrast_names:
        deg_path = find_deg_file(root_dir, name)
        if not deg_path:
            logger.warning(f"No DEG file found for contrast: {name}")
            continue

        found_deg_files += 1
        df = pd.read_csv(deg_path)
        if "Gene" not in df.columns:
            logger.error(f"'Gene' column not found in {deg_path}")
            raise KeyError(f"'Gene' column not found in {deg_path}")

        all_genes_in_file = set(df["Gene"].astype(str))
        total_genes_before_filtering += len(all_genes_in_file)

        # Filter genes by logFC and adjusted p-value thresholds
        if "logFC" in df.columns and "adj.P.Val" in df.columns:
            filtered_df = df[(abs(df["logFC"]) >= lfc_threshold) &
                            (df["adj.P.Val"] <= padj_threshold)]
            genes = set(filtered_df["Gene"].astype(str))

            logger.info(f"Contrast '{name}': {len(genes)}/{len(df)} genes pass thresholds "
                      f"(|logFC| >= {lfc_threshold}, adj.P.Val <= {padj_threshold})")
        else:
            # Fall back to old behavior if expected columns not found
            logger.warning(f"Missing logFC or adj.P.Val columns in {deg_path}. Using all genes.")
            genes = all_genes_in_file

        gene_sets[name] = genes
        all_genes |= genes

    logger.info(f"Processed {found_deg_files}/{len(contrast_names)} DEG files")
    logger.info(f"Total genes before filtering: {total_genes_before_filtering}")
    logger.info(f"Total unique genes after filtering: {len(all_genes)}")

    # 3) Build the binary presence/absence DataFrame
    matrix = pd.DataFrame(
        0,
        index=sorted(all_genes),
        columns=sorted(gene_sets.keys()),
        dtype=int,
    )
    for name, genes in gene_sets.items():
        matrix.loc[list(genes), name] = 1

    logger.info(f"Created binary matrix with {matrix.shape[0]} rows × {matrix.shape[1]} columns")
    return matrix

def main():
    p = argparse.ArgumentParser(description="Build gene×contrast binary matrix")
    p.add_argument("root_dir",
                   help="Root folder under which analysis results (contrasts.csv & deg_*.csv) live")
    p.add_argument("--outdir", "-o",
                   default="integration_matrix.csv",
                   help="Path to write matrix CSV")
    p.add_argument("--lfc", "-l", type=float, default=1.0,
                   help="Log fold change threshold (absolute value, default=1.0)")
    p.add_argument("--padj", "-p", type=float, default=0.05,
                   help="Adjusted p-value threshold (default=0.05)")
    p.add_argument("--verbose", "-v", action="store_true",
                   help="Enable verbose logging")
    p.add_argument("--html", action="store_true",
                   help="Also write an interactive HTML file")
    p.add_argument("--html-type", choices=["datatables", "panel"], default="panel",
                   help="Type of HTML to generate: datatables (old) or panel (new, default)")
    args = p.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    logger.info(f"Starting analysis on directory: {args.root_dir}")
    logger.info(f"Using thresholds: |logFC| >= {args.lfc}, adj.P.Val <= {args.padj}")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)  # safe even if it already exists

    logger.info(f"Writing all outputs under: {outdir.resolve()}")

    matrix_path = outdir / "integrated_results.csv"
    filtered_path = outdir / "integrated_results_filtered.csv"
    html_path = outdir / "integrated_results_filtered.html"

    # --- original logic ---
    matrix = build_binary_matrix(args.root_dir, args.lfc, args.padj)
    
    # Add explicit Gene header
    matrix.index.name = "Gene"
    matrix.to_csv(matrix_path)
    logger.info(f"Wrote integration matrix to: {matrix_path}")

    filtered = matrix[matrix.sum(axis=1) > 0].copy()
    filtered["Num_DE_Contrasts"] = filtered.sum(axis=1)
    
    # Add explicit Gene header
    filtered.index.name = "Gene"
    filtered.to_csv(filtered_path)

    logger.info(f"Wrote filtered integration matrix to: {filtered_path}")

    # --- write HTML if requested ---------------------------
    if args.html:
        if args.html_type == "datatables":
            write_datatables_html(filtered, html_path)
            logger.info(f"Wrote interactive DataTables HTML table to: {html_path}")
        else:  # panel
            write_panel_bokeh_html(filtered, html_path)
            logger.info(f"Wrote interactive Panel+Bokeh HTML visualization to: {html_path}")

    print(
        f"✓ Results:\n"
        f"  • {matrix_path.name} ({matrix.shape[0]}×{matrix.shape[1]})\n"
        f"  • {filtered_path.name} ({filtered.shape[0]}×{filtered.shape[1]})\n"
        f"  • {html_path.name if args.html else '(HTML skipped)'}"
    )


if __name__ == "__main__":
    main()
