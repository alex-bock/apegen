
import sys

import networkx as nx
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.spatial.distance import pdist, squareform


if __name__ == "__main__":

    # load template

    template_df = pd.read_csv("./helper_files/Proper_files/Template_DB.csv")
    if len(sys.argv) > 1:
        template_df = template_df[template_df.MHC == sys.argv[1]]
    template_df = template_df[template_df.peptide_length == 9]
    print(template_df)

    # compute pairwise template distances

    peptide_arrs = template_df.peptide.apply(lambda x: [ord(c) for c in x])
    n_peptides = len(peptide_arrs)
    X = np.array(peptide_arrs.values.tolist())
    dist_mat = squareform(pdist(X, metric="hamming"))

    # add node for each template

    graph = nx.MultiGraph()
    graph.add_nodes_from(range(len(peptide_arrs)))

    # add edge layer(s) for distance threshold(s)

    thresholds = [(i / 9) for i in range(9)]

    for i in range(1, len(thresholds)):

        print(thresholds[i - 1], thresholds[i])
        print(((dist_mat >= thresholds[i - 1]) & (dist_mat <= thresholds[i])).nonzero())
        
        graph.add_edges_from(np.transpose(((dist_mat >= thresholds[i - 1]) & (dist_mat <= thresholds[i])).nonzero()).tolist(), dist=i)

    # compute node layout

    layout = nx.spring_layout(graph)

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=[layout[i][0] for i in range(len(graph.nodes))],
            y=[layout[i][1] for i in range(len(graph.nodes))],
            mode="markers",
            marker=dict(symbol="circle", size=10),
            hoverinfo="text",
            text=template_df.peptide
        )
    )

    # draw edge layer(s)

    for i in range(1, len(thresholds)):

        edges = [(u, v) for (u, v, *_, data) in graph.edges(data=True) if data["dist"] == i]
        print(len(edges))

        fig.add_trace(
            go.Scatter(
                x=np.array([[float(layout[u][0]), float(layout[v][0]), None] for (u, v) in edges]).flatten(),
                y=np.array([[float(layout[u][1]), float(layout[v][1]), None] for (u, v) in edges]).flatten(),
                mode="lines",
                opacity=0.5
            )
        )

    fig.show()