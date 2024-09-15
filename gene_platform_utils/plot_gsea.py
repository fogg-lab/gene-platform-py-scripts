from typing import Union
import networkx as nx
import plotly.graph_objects as go
import numpy as np


def gene_concept_network_plot(
    gsea_res: dict[str, dict[str, Union[str, float]]],
    de_res: dict[str, dict[str, Union[str, float]]],
    color_metric: str = "logFC",
    pvalue_threshold: float = 0.05,
    layout_seed: int = 0,
    color_seed: int = 0,
) -> str:
    # Filter DE results by p-value
    de_res = {k: v for k, v in de_res.items() if v["P.Value"] < pvalue_threshold}

    # Create a new graph
    G = nx.Graph()

    # Assign unique colors for pathways
    rng = np.random.default_rng(color_seed)
    pathway_colors = {
        pathway_id: f"rgb({rng.integers(0, 256)}, {rng.integers(0, 256)}, {rng.integers(0, 256)})"
        for pathway_id in gsea_res.keys()
    }

    # Add gene nodes and edges, as well as pathway nodes and their sizes
    pathway_sizes = {}
    all_genes = set(de_res)
    genes_not_added = all_genes.copy()
    for pathway_id, pathway_data in gsea_res.items():
        pathway_genes = set(pathway_data["core_enrichment"].split("/"))
        pathway_genes.intersection_update(all_genes)
        pathway_sizes[pathway_id] = len(pathway_genes)
        G.add_node(
            pathway_id,
            size=len(pathway_genes),
            label=pathway_data["Description"],
        )
        new_node_genes = pathway_genes.intersection(genes_not_added)
        genes_not_added -= new_node_genes
        for gene in new_node_genes:
            G.add_node(gene, color_metric=de_res[gene][color_metric])
        for gene in pathway_genes:
            G.add_edge(gene, pathway_id, color=pathway_colors[pathway_id])

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
            width=1000,
            height=800,
            showlegend=True,
            hovermode="closest",
            margin=dict(b=0, l=0, r=0, t=20),
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        ),
    )

    return fig.to_html()
