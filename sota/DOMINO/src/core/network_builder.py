import networkx as nx
import pandas as pd


def build_network(network_file):
    """"""
    edges_dataset = pd.read_csv(network_file, sep='\t', header=0, dtype=str)
    edges = []
    for ind, row in edges_dataset.iterrows():
        if ind == 0:
            continue
        edges.append((row.iloc[0], row.iloc[1]))
    
    G = nx.Graph()
    G.add_edges_from(edges)
    nx.set_node_attributes(G, 0, 'score')

    return G