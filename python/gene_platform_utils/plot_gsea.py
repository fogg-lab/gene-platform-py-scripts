from typing import TYPE_CHECKING
import networkx as nx
import plotly.graph_objects as go
import numpy as np

if TYPE_CHECKING:
    import pandas as pd


def gene_concept_network_plot(
    gsea_res: pd.DataFrame,
    de_res: pd.DataFrame,
    genes_df: pd.DataFrame,
    color_metric: str = "log2FoldChange",
    pvalue_threshold: float = 0.05,
    layout_seed: int = 0,
    color_seed: int = 0,
) -> str:
    # Filter DE results by p-value
    de_res = de_res[de_res["pvalue"] < pvalue_threshold]

    # Create a new graph
    G = nx.Graph()

    # Assign unique colors for pathways
    rng = np.random.default_rng(color_seed)
    pathway_colors = {
        pathway_id: f"rgb({rng.integers(0, 256)}, {rng.integers(0, 256)}, {rng.integers(0, 256)})"
        for pathway_id in gsea_res.index
    }

    # Add pathway nodes and their sizes
    pathway_sizes = {}
    for pathway_id in gsea_res.index:
        pathway_genes = set(
            gsea_res.loc[pathway_id, "core_enrichment_ensembl"].split("/")
        )
        pathway_sizes[pathway_id] = len(pathway_genes)
        G.add_node(
            pathway_id,
            size=len(pathway_genes),
            label=gsea_res.loc[pathway_id, "Description"],
        )

    # Add gene nodes and edges
    for gene in de_res.index:
        for pathway_id in gsea_res.index:
            pathway_genes = set(
                gsea_res.loc[pathway_id, "core_enrichment_ensembl"].split("/")
            )
            if gene in pathway_genes:
                gene_sym = genes_df.loc[
                    genes_df.index.str.startswith(gene), "gene_name"
                ].values[0]
                G.add_node(gene_sym, color_metric=de_res.loc[gene, color_metric])
                G.add_edge(gene_sym, pathway_id, color=pathway_colors[pathway_id])

    # Layout with seed for reproducibility
    pos = nx.fruchterman_reingold_layout(G, k=1, iterations=400, seed=layout_seed)

    # Prepare the data for Plotly
    edge_traces = []
    for edge in G.edges(data=True):
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_trace = go.Scatter(
            x=[x0, x1, None],
            y=[y0, y1, None],
            line=dict(width=0.5, color=edge[2]["color"]),
            hoverinfo="none",
            mode="lines",
            showlegend=False,
        )
        edge_traces.append(edge_trace)

    pathway_traces = []
    gene_x = []
    gene_y = []
    gene_text = []
    gene_hovertext = []
    gene_color = []
    gene_size = []

    for node in G.nodes():
        x, y = pos[node]
        if "label" in G.nodes[node]:
            # Pathway nodes
            pathway_node_x = [x]
            pathway_node_y = [y]
            pathway_node_size = [G.nodes[node]["size"] / 5]
            pathway_node_text = [""]
            pathway_node_hovertext = [G.nodes[node]["label"]]

            pathway_trace = go.Scatter(
                x=pathway_node_x,
                y=pathway_node_y,
                mode="markers",
                text=pathway_node_text,
                hovertext=pathway_node_hovertext,
                hoverinfo="text",
                marker=dict(
                    color=pathway_colors[node], size=pathway_node_size, line_width=2
                ),
                name=G.nodes[node]["label"],
                showlegend=True,
            )
            pathway_traces.append(pathway_trace)
        else:
            # Gene nodes
            gene_x.append(x)
            gene_y.append(y)
            gene_text.append("")
            gene_hovertext.append(node)
            gene_color.append(G.nodes[node]["color_metric"])
            gene_size.append(10)

    gene_trace = go.Scatter(
        x=gene_x,
        y=gene_y,
        mode="markers",
        text=gene_text,
        hovertext=gene_hovertext,
        hoverinfo="text",
        marker=dict(
            showscale=True,
            colorscale="Viridis",
            color=gene_color,
            size=gene_size,
            colorbar=dict(
                thickness=15,
                title=color_metric,
                xanchor="left",
                titleside="right",
                y=0.3,
                len=0.5,
            ),
            line_width=2,
        ),
        showlegend=False,
    )

    fig = go.Figure(
        data=edge_traces + [gene_trace] + pathway_traces,
        layout=go.Layout(
            width=1000,  # Set the width of the figure
            height=800,  # Set the height of the figure
            showlegend=True,
            hovermode="closest",
            margin=dict(b=0, l=0, r=0, t=20),
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        ),
    )

    return fig.to_html()
