"""Script and documentation generation utilities for reproducible heatmaps."""


def build_repro_script() -> str:
    """
    Return contents of a standalone script that reproduces the heatmap using Plotly.

    This script mirrors the Streamlit app's heatmap (px.imshow + layout), ensuring identical visuals
    when run with default settings. It supports an optional JSON file to remap x-axis labels.

    Returns:
        Complete Python script as a string
    """
    script_lines = [
        '#!/usr/bin/env python3',
        '"""',
        'Reproduce the UORCA heatmap from `heatmap_data.csv` using Plotly (matches the app).',
        '',
        'Usage (recommended inside UORCA repo with uv):',
        '  uv run python reproduce_heatmap.py --input heatmap_data.csv --output heatmap.pdf',
        '',
        'Optionally remap x-axis labels with a JSON file:',
        '  uv run python reproduce_heatmap.py --input heatmap_data.csv --output heatmap.pdf --xlabel-mapping labels.json',
        '',
        'The optional `labels.json` should map existing column names to new labels,',
        'e.g. {"analysis1_contrastA": "Control vs Treated", "analysis2_contrastB": "Condition X"}',
        '"""',
        'import argparse',
        'import json',
        'import pandas as pd',
        'import numpy as np',
        'import plotly.express as px',
        'import plotly.graph_objects as go',
        '',
        'def main():',
        '    meta_help = (',
        "        'Metadata tweaks (edit metadata.json):\\n'",
        "        '  - output_format: pdf|png|svg (default pdf)\\n'",
        "        '  - dpi: integer DPI for export (default 300)\\n'",
        "        '  - font_size: integer/float font size (default 12)\\n'",
        "        '  - show_grid_lines: true|false (default true)\\n'",
        "        '  - grid_opacity: 0.0-1.0 (default 0.3)\\n'",
        "        '  - x_axis_labels: {column_name: display_label, ...}\\n'",
        "        '  - clustered_contrasts: [column_name, ...] (order applied)\\n'",
        '    )',
        '    parser = argparse.ArgumentParser(description=meta_help, formatter_class=argparse.RawDescriptionHelpFormatter)',
        "    parser.add_argument('--input', default='heatmap_data.csv')",
        "    parser.add_argument('--output', default=None, help='Output file path (default heatmap.<format> based on metadata)')",
        "    parser.add_argument('--xlabel-mapping', default=None, help='JSON file mapping original column names to new labels (overrides metadata)')",
        '    args = parser.parse_args()',
        '',
        '    print("[1/6] Loading metadata.json (if present)...", flush=True)',
        '    # Load metadata.json if present and use it as defaults',
        '    meta = {}',
        '    try:',
        "        with open('metadata.json', 'r') as mfh:",
        '            meta = json.load(mfh)',
        '    except Exception:',
        '        meta = {}',
        '',
        '    # Resolve parameters primarily from metadata.json',
        "    out_fmt = (meta.get('output_format') or 'pdf').lower()",
        "    if out_fmt not in ('pdf','png','svg'): out_fmt = 'pdf'",
        "    dpi = meta.get('dpi', 300)",
        "    font_size = meta.get('font_size', 12)",
        "    show_grid = True if meta.get('show_grid_lines', True) else False",
        "    grid_opacity = meta.get('grid_opacity', 0.3)",
        '',
        '    # Determine output path',
        "    default_output = f'heatmap.{out_fmt}'",
        "    output_path = args.output or default_output",
        '',
        '    print(f"[2/6] Reading input heatmap data: {args.input}")',
        '    # Read input data',
        '    df = pd.read_csv(args.input)',
        "    if 'Gene' in df.columns:",
        "        df = df.set_index('Gene')",
        '    df = df.fillna(0)',
        '',
        '    print("[3/6] Applying label mapping and column order...")',
        '    # Apply xlabel mapping: CLI file if provided, otherwise metadata x_axis_labels',
        '    xlabel_map = {}',
        '    if args.xlabel_mapping:',
        "        with open(args.xlabel_mapping, 'r') as fh:",
        '            xlabel_map = json.load(fh)',
        '    else:',
        "        xlabel_map = meta.get('x_axis_labels', {})",
        '',
        '    # Reorder columns if metadata provides clustered_contrasts',
        "    clustered = meta.get('clustered_contrasts')",
        '    if clustered and set(clustered).issubset(set(df.columns)):',
        '        df = df[clustered]',
        '',
        '    # Replace column labels according to mapping, preserving order',
        '    simplified_cols = [xlabel_map.get(c, c) for c in df.columns]',
        '',
        '    print("[4/6] Computing layout and color scaling...")',
        '    # Compute symmetric color scale range around 0',
        '    max_abs = float(max(abs(df.values.min()), abs(df.values.max())))',
        '    if max_abs == 0:\n            max_abs = 1.0',
        '',
        '    # Dynamic sizing to match app logic',
        '    try:',
        '        max_gene_label_len = max((len(str(g)) for g in df.index), default=10)',
        '    except Exception:',
        '        max_gene_label_len = 10',
        '    row_height = max(18, int(font_size * 1.8))',
        '    dynamic_height = int(140 + df.shape[0] * row_height)',
        '    height_px = max(650, min(4500, dynamic_height))',
        '    left_margin = max(250, min(520, int(40 + max_gene_label_len * font_size * 0.7)))',
        '    width_px = max(700, min(1700, df.shape[1] * 90))',
        '',
        '    print("[5/6] Building figure with Plotly...")',
        '    # Build figure',
        '    fig = px.imshow(',
        '        df,',
        "        color_continuous_scale='RdBu_r',",
        "        labels=dict(x='Contrast', y='Gene', color='Log2FC'),",
        "        title='Differential Expression Heatmap (Log2 Fold Change)',",
        "        aspect='auto',",
        '        height=height_px,',
        '        width=width_px,',
        '        zmin=-max_abs,',
        '        zmax=max_abs,',
        '        color_continuous_midpoint=0',
        '    )',
        '',
        '    fig.update_layout(',
        '        margin=dict(l=left_margin, r=20, t=60, b=120),',
        "        coloraxis_colorbar=dict(title='Log2FC'),",
        "        font_family='Inter, sans-serif',",
        '        font_size=font_size',
        '    )',
        '',
        '    # Update x-axis tick labels with simplified names, respecting clustered order',
        '    tick_vals = list(range(len(simplified_cols)))',
        '    fig.update_xaxes(tickmode="array", tickvals=tick_vals, ticktext=simplified_cols, tickangle=45, title="Biological contrast")',
        '',
        '    # Grid lines',
        '    if show_grid:',
        "        rgba = f'rgba(200,200,200,{grid_opacity})'",
        '        fig.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor=rgba)',
        '        fig.update_xaxes(showgrid=True, gridwidth=0.5, gridcolor=rgba)',
        '    else:',
        '        fig.update_yaxes(showgrid=False)',
        '        fig.update_xaxes(showgrid=False)',
        '',
        '    print(f"[6/6] Saving figure ({out_fmt.upper()}) to: {output_path}")',
        '    # Save',
        '    try:',
        '        if out_fmt == "pdf":',
        '            fig.write_image(output_path, format="pdf", scale=2, width=width_px, height=height_px)',
        '        elif out_fmt == "svg":',
        '            fig.write_image(output_path, format="svg", scale=2, width=width_px, height=height_px)',
        '        else:',
        '            fig.write_image(output_path, format="png", scale=2, width=width_px, height=height_px)',
        '    except Exception as exc:',
        '        import traceback',
        '        report = {',
        "            'error': str(exc),",
        "            'traceback': traceback.format_exc(),",
        "            'note': 'Ensure plotly and kaleido are available or run via the UORCA uv environment (e.g. `uv run python reproduce_heatmap.py ...`).'",
        '        }',
        "        with open('error_report.txt', 'w') as ef:",
        "            ef.write('Error reproducing heatmap\\n')",
        "            ef.write('Error: ' + report['error'] + '\\n\\n')",
        "            ef.write(report['traceback'])",
        "            ef.write('\\n\\n')",
        "            ef.write(report['note'])",
        '        raise',
        '',
        "if __name__ == '__main__':",
        '    main()',
        ''
    ]
    return '\n'.join(script_lines) + '\n'


def build_readme_text() -> str:
    """
    Generate README text for reproducible heatmap packages.

    Returns:
        README text as a string
    """
    return (
        "UORCA Heatmap Reproducible Package\n"
        "--------------------------------\n\n"
        "Included files:\n"
        "- heatmap_data.csv: CSV table of genes (rows) and contrasts (columns) containing log2 fold change values (non-significant values set to 0).\n"
        "- metadata.json: JSON metadata including thresholds, clustered contrast order, and default x-axis label mapping.\n"
        "- reproduce_heatmap.py: Script to generate the same heatmap as in the app.\n"
        "- README.txt: This file.\n\n"
        "Quick start (recommended):\n"
        "1. Unzip the package.\n"
        "2. Optional but helpful: move this folder into your UORCA repository root so you can use the uv environment (e.g., UORCA/your_heatmap_package/).\n"
        "3. Run via uv to produce a high-resolution PDF (default):\n"
        "   uv run python reproduce_heatmap.py\n\n"
        "Tips:\n"
        "- Provide an --xlabel-mapping JSON to rename contrast axis labels for publication.\n"
        "- PDF is the default format for publication-quality vector output; PNG/SVG are also supported.\n"
        "- The script reproduces the Streamlit heatmap exactly under default settings.\n\n"
        "If you need the original DEG/CPM source files included for provenance, ask to include them in the package."
    )