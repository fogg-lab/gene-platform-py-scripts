import plotly.graph_objects as go
import numpy as np


def create_volcano_plot(
    data: np.ndarray,
    row_names: list[str],
    column_names: list[str],
    lfc_thresh: float,
    pval_thresh: float,
    cohort_name: str,
) -> str:
    """
    Create an interactive volcano plot using Plotly.
    """
    log2fc_col = column_names.index("logFC")
    pvalue_col = column_names.index("P.Value")

    diff_expr = np.full(len(row_names), "Not sig.")
    diff_expr[
        (data[:, log2fc_col] > lfc_thresh) & (data[:, pvalue_col] < pval_thresh)
    ] = "Up"
    diff_expr[
        (data[:, log2fc_col] < -lfc_thresh) & (data[:, pvalue_col] < pval_thresh)
    ] = "Down"

    color_map = {"Up": "#B31B21", "Down": "#1465AC", "Not sig.": "darkgray"}

    fig = go.Figure()

    for category in ["Up", "Down", "Not sig."]:
        mask = diff_expr == category
        fig.add_trace(
            go.Scatter(
                x=data[mask, log2fc_col],
                y=-np.log10(data[mask, pvalue_col]),
                mode="markers",
                name=category,
                marker=dict(color=color_map[category], size=5),
                text=[row_names[i] for i in np.where(mask)[0]],
                hoverinfo="text+x+y",
            )
        )

    fig.add_shape(
        type="line",
        x0=-lfc_thresh,
        x1=-lfc_thresh,
        y0=0,
        y1=1,
        yref="paper",
        line=dict(color="#B31B21", width=1, dash="dash"),
    )
    fig.add_shape(
        type="line",
        x0=lfc_thresh,
        x1=lfc_thresh,
        y0=0,
        y1=1,
        yref="paper",
        line=dict(color="#B31B21", width=1, dash="dash"),
    )
    fig.add_shape(
        type="line",
        x0=-8,
        x1=8,
        y0=-np.log10(pval_thresh),
        y1=-np.log10(pval_thresh),
        line=dict(color="#B31B21", width=1, dash="dash"),
    )

    fig.update_layout(
        title=f"Volcano Plot - {cohort_name}",
        xaxis_title="Log2 Fold Change",
        yaxis_title="-Log10 P-value",
        legend_title="Differential Expression",
        template="plotly_white",
    )

    return fig.to_html()


def create_mean_difference_plot(
    data: np.ndarray,
    row_names: list[str],
    column_names: list[str],
    fdr: float,
    cohort_name: str,
) -> str:
    """
    Create an interactive mean difference plot using Plotly.
    """
    basemean_col = column_names.index("AveExpr")
    log2fc_col = column_names.index("logFC")
    padj_col = column_names.index("adj.P.Val")

    mean_expression = np.log2((data[:, basemean_col] + 0.5) / 2)
    significant = data[:, padj_col] < fdr

    fig = go.Figure()

    # Add non-significant points
    fig.add_trace(
        go.Scatter(
            x=mean_expression[~significant],
            y=data[~significant, log2fc_col],
            mode="markers",
            marker=dict(color="darkgray", size=5),
            name="Not significant",
            text=[row_names[i] for i in np.where(~significant)[0]],
            hoverinfo="text+x+y",
        )
    )

    # Add significant points
    fig.add_trace(
        go.Scatter(
            x=mean_expression[significant],
            y=data[significant, log2fc_col],
            mode="markers",
            marker=dict(
                color=[
                    "#B31B21" if x > 0 else "#1465AC"
                    for x in data[significant, log2fc_col]
                ],
                size=5,
            ),
            name="Significant",
            text=[row_names[i] for i in np.where(significant)[0]],
            hoverinfo="text+x+y",
        )
    )

    fig.update_layout(
        title=f"Mean Difference Plot - {cohort_name}",
        xaxis_title="Log2 mean expression",
        yaxis_title="Log2 fold change",
        legend_title="Significance",
        template="plotly_white",
    )

    return fig.to_html()
