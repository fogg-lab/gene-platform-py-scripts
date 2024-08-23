import argparse
import os
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.cluster import hierarchy


def read_data(counts_file: str, coldata_file: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Read the counts and coldata files.

    Args:
        counts_file (str): Path to the VST-transformed counts CSV file.
        coldata_file (str): Path to the coldata CSV file.

    Returns:
        tuple[pd.DataFrame, pd.DataFrame]: A tuple containing the counts and coldata DataFrames.
    """
    counts = pd.read_csv(counts_file, index_col=0)
    counts.index.name = "Ensembl_gene_id"

    coldata = pd.read_csv(coldata_file, index_col=0)
    coldata.index.name = "sample_id"

    return counts, coldata


def calculate_correlation(counts: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate the sample-to-sample correlation matrix using np.corrcoef.

    Args:
        counts (pd.DataFrame): The VST-transformed counts DataFrame.

    Returns:
        pd.DataFrame: The correlation matrix.
    """
    corr = np.corrcoef(counts.T)
    corr = pd.DataFrame(corr, index=counts.columns, columns=counts.columns)
    # return correlation matrix with samples sorted based on hierarchical clustering
    linkage = hierarchy.linkage(corr, method="average")
    samp_order = hierarchy.dendrogram(linkage, no_plot=True)["leaves"]
    return corr.iloc[samp_order, samp_order]


def create_heatmap(
    corr_matrix: pd.DataFrame, coldata: pd.DataFrame, output_dir: str
) -> None:
    """
    Create and save an interactive heatmap of sample-to-sample correlation using Plotly.

    Args:
        corr_matrix (pd.DataFrame): The correlation matrix.
        coldata (pd.DataFrame): The coldata DataFrame.
        output_dir (str): Directory to save the output HTML file.

    """

    def get_sample_info(sample):
        sample_info = [f"{sample}"]
        for trait, value in list(coldata.loc[sample].items())[:3]:
            sample_info.append(f"{trait}: {value}")
        return "<br>".join(sample_info)

    hover_text = [
        [
            f"X: {x}<br>" f"Y:{y}<br>" f"r = {corr_matrix.loc[y, x]:.4f}"
            for x in corr_matrix.columns
        ]
        for y in corr_matrix.index
    ]

    fig = go.Figure(
        data=go.Heatmap(
            z=corr_matrix.values,
            x=corr_matrix.columns,
            y=corr_matrix.index,
            hoverongaps=False,
            text=hover_text,
            hoverinfo="text",
            colorscale="Viridis",
            colorbar=dict(
                title="<b>r</b>",
            ),
            hovertemplate="%{text}<extra></extra>",
        )
    )

    fig.update_layout(
        title="Sample-to-Sample Pearson Correlation Heatmap<br><sup>Hover over a cell to view sample IDs.</sup>",
        xaxis=dict(showticklabels=False, title=None),
        yaxis=dict(showticklabels=False, title=None),
        width=800,
        height=800,
        hoverlabel=dict(
            bgcolor="rgb(255, 255, 255)",
            font=dict(color="black"),
        ),
    )

    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, "sample_correlation_heatmap.html")
    fig.write_html(output_file)

    # Add CSS to make the hover box 60% transparent
    # This is a workaround for the issue with the hover box covering the heatmap
    # First thought would be to change the hoverlabel bgcolor to rgba(255, 255, 255, 0.4)
    # but this has no effect unless hovermode="x unified" is also set. what the heck plotly
    with open(output_file, "r") as f:
        html = f.read()
    css = "<style>.hovertext { fill-opacity: 0.4; stroke-opacity: 1; }</style>"
    html = html.replace("</head>", f"{css}</head>")
    with open(output_file, "w") as f:
        f.write(html)  # write the modified html back to the file

    print(f"Saved sample correlation heatmap to {output_file}")


def main(counts_file: str, coldata_file: str, output_dir: str) -> None:
    """
    Main function to create a sample-to-sample correlation heatmap.

    Args:
        counts_file (str): Path to the VST-transformed counts CSV file.
        coldata_file (str): Path to the coldata CSV file.
        output_dir (str): Directory to save the output HTML file.

    """
    counts, coldata = read_data(counts_file, coldata_file)
    corr_matrix = calculate_correlation(counts)
    create_heatmap(corr_matrix, coldata, output_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create a sample-to-sample correlation heatmap using Plotly"
    )
    parser.add_argument(
        "counts_file", help="Path to the VST-transformed counts CSV file"
    )
    parser.add_argument("coldata_file", help="Path to the coldata CSV file")
    parser.add_argument("output_dir", help="Directory to save the output HTML file")

    args = parser.parse_args()

    main(args.counts_file, args.coldata_file, args.output_dir)
