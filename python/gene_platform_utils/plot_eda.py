import numpy as np
import plotly.graph_objects as go
from scipy.cluster import hierarchy


def calculate_correlation(counts: np.ndarray) -> np.ndarray:
    """
    Calculate the sample-to-sample correlation matrix using np.corrcoef.

    Args:
        counts (np.ndarray): The VST-transformed counts matrix.

    Returns:
        np.ndarray: The correlation matrix.
    """
    corr = np.corrcoef(counts.T)
    # return correlation matrix with samples sorted based on hierarchical clustering
    linkage = hierarchy.linkage(corr, method="average")
    samp_order = hierarchy.dendrogram(linkage, no_plot=True)["leaves"]
    corr = corr[samp_order, :][:, samp_order]
    return corr, samp_order


def create_correlation_heatmap(counts: np.ndarray, sample_ids: list[str]) -> str:
    """
    Create an interactive heatmap of sample-to-sample correlation using Plotly.

    Args:
        counts (np.ndarray): The counts matrix.
        sample_ids (list[str]): The sample IDs.

    Returns:
        str: The HTML string of the heatmap.
    """
    corr_matrix, samp_order = calculate_correlation(counts)
    sample_ids = [sample_ids[i] for i in samp_order]

    fig = go.Figure(
        data=go.Heatmap(
            z=corr_matrix,
            x=sample_ids,
            y=sample_ids,
            hoverongaps=False,
            colorscale="Viridis",
            colorbar=dict(
                title="<b>r</b>",
            ),
            hovertemplate="X: %{x}<br>Y: %{y}<br>r = %{z:.4f}<extra></extra>",
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

    html = fig.to_html()
    plotly_cdn = '<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>'
    css = "<style>.hovertext { fill-opacity: 0.4; stroke-opacity: 1; }</style>"
    html = html.replace("</head>", f"{plotly_cdn}{css}</head>")

    return html
