import os
import pandas as pd
import plotly.graph_objects as go
import numpy as np


def create_volcano_plot(
    df: pd.DataFrame,
    lfc_thresh: float,
    pval_thresh: float,
    cohort_name: str,
) -> str:
    """
    Create an interactive volcano plot using Plotly.

    Args:
        df (pd.DataFrame): DataFrame containing differential expression results.
        lfc_thresh (float): Log fold change threshold for significance.
        pval_thresh (float): P-value threshold for significance.
        cohort_name (str): Name of the cohort for plot title.

    """
    df["diff_expr"] = "Not sig."
    df.loc[
        (df["log2FoldChange"] > lfc_thresh) & (df["pvalue"] < pval_thresh), "diff_expr"
    ] = "Up"
    df.loc[
        (df["log2FoldChange"] < -lfc_thresh) & (df["pvalue"] < pval_thresh), "diff_expr"
    ] = "Down"

    color_map = {"Up": "#B31B21", "Down": "#1465AC", "Not sig.": "darkgray"}

    fig = go.Figure()

    for category in ["Up", "Down", "Not sig."]:
        subset = df[df["diff_expr"] == category]
        fig.add_trace(
            go.Scatter(
                x=subset["log2FoldChange"],
                y=-np.log10(subset["pvalue"]),
                mode="markers",
                name=category,
                marker=dict(color=color_map[category], size=5),
                text=subset["gene"],
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


def create_mean_difference_plot(df: pd.DataFrame, fdr: float, cohort_name: str) -> str:
    """
    Create an interactive mean difference plot using Plotly.

    Args:
        df (pd.DataFrame): DataFrame containing differential expression results.
        fdr (float): False discovery rate threshold for significance.
        cohort_name (str): Name of the cohort for plot title.

    """
    df["mean_expression"] = np.log2((df["baseMean"] + 0.5) / 2)
    df["significant"] = df["padj"] < fdr

    fig = go.Figure()

    # Add non-significant points
    fig.add_trace(
        go.Scatter(
            x=df[~df["significant"]]["mean_expression"],
            y=df[~df["significant"]]["log2FoldChange"],
            mode="markers",
            marker=dict(color="darkgray", size=5),
            name="Not significant",
            text=df[~df["significant"]]["gene"],
            hoverinfo="text+x+y",
        )
    )

    # Add significant points
    fig.add_trace(
        go.Scatter(
            x=df[df["significant"]]["mean_expression"],
            y=df[df["significant"]]["log2FoldChange"],
            mode="markers",
            marker=dict(
                color=df[df["significant"]]["log2FoldChange"].apply(
                    lambda x: "#B31B21" if x > 0 else "#1465AC"
                ),
                size=5,
            ),
            name="Significant",
            text=df[df["significant"]]["gene"],
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


def main(
    input_file: str,
    output_dir: str,
    lfc_thresh: float = 1,
    pval_thresh: float = 0.05,
    fdr: float = 0.05,
) -> None:
    """
    Main function to create interactive Plotly plots for differential expression analysis results.

    Args:
        input_file (str): Path to the input CSV file containing differential expression results.
        output_dir (str): Directory to save the output HTML files.
        lfc_thresh (float, optional): Log fold change threshold for significance. Defaults to 1.
        pval_thresh (float, optional): P-value threshold for significance. Defaults to 0.05.
        fdr (float, optional): False discovery rate threshold for significance. Defaults to 0.05.

    """
    df = pd.read_csv(input_file)
    cohort_name = os.path.splitext(os.path.basename(input_file))[0].replace(
        "DE_results_", ""
    )

    create_volcano_plot(df, lfc_thresh, pval_thresh, output_dir, cohort_name)
    create_mean_difference_plot(df, fdr, output_dir, cohort_name)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Create interactive Plotly plots for differential expression analysis results"
    )
    parser.add_argument("input_file", help="Path to the input CSV file")
    parser.add_argument("output_dir", help="Directory to save the output HTML files")
    parser.add_argument(
        "--lfc_thresh", type=float, default=1, help="Log fold change threshold"
    )
    parser.add_argument(
        "--pval_thresh", type=float, default=0.05, help="P-value threshold"
    )
    parser.add_argument("--fdr", type=float, default=0.05, help="False discovery rate")

    args = parser.parse_args()

    main(args.input_file, args.output_dir, args.lfc_thresh, args.pval_thresh, args.fdr)
