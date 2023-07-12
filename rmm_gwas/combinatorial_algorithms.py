import networkx as nx
import random as rn
import math
import pprint as pp
import time
import datetime

current_milli_time = lambda: int(round(time.time() * 1000))


##################################################################
##################################################################
##################################################################

def compute_density_of_the_induced_subgraph(Vs, g):
    induced_subgraph = g.subgraph(Vs)
    return get_sum_of_all_edge_weights(induced_subgraph) / float(len(Vs))


###
def get_sum_of_all_edge_weights(graph):
    sum_of_all_edge_weights = sum([weight for u, v, weight in graph.edges.data('weight', default=1.)])
    return float(sum_of_all_edge_weights)


###
def create_Goldberg_graph(G):
    best_Vs = set(G.nodes())
    #
    m = get_sum_of_all_edge_weights(G)
    n = G.number_of_nodes()
    #
    graph_Goldberg = nx.DiGraph(G)
    for di_edge in graph_Goldberg.edges(data=False):
        node_u = di_edge[0]
        node_v = di_edge[1]
        weight = G.get_edge_data(node_u, node_v)['weight']
        graph_Goldberg[node_u][node_v]['weight'] = weight
        graph_Goldberg[node_v][node_u]['weight'] = weight
    #
    source_Goldberg = -1  # source node
    sink_Goldberg = -2  # sink node
    graph_Goldberg.add_node(source_Goldberg)
    graph_Goldberg.add_node(sink_Goldberg)
    for node in G:
        graph_Goldberg.add_edge(source_Goldberg, node, weight=float(m))
    #
    g = 0.
    for node in G:
        # degree_of_node = G.degree(node)
        degree_of_node = G.degree(node, weight="weight")
        weight_to_sink = m + 2. * g - degree_of_node
        graph_Goldberg.add_edge(node, sink_Goldberg, weight=weight_to_sink)
    #
    return graph_Goldberg, source_Goldberg, sink_Goldberg


###
def update_weights_on_Goldberg_graph(original_graph, graph_Goldberg, density_of_subgraph, sink_Goldberg):
    # update the weights on the edges to the sink node.
    g = float(density_of_subgraph)
    # m = original_graph.number_of_edges()
    m = get_sum_of_all_edge_weights(original_graph)
    for node in original_graph:
        degree_of_node = original_graph.degree(node, weight="weight")
        weight_to_sink = m + 2. * g - degree_of_node
        graph_Goldberg[node][sink_Goldberg]['weight'] = weight_to_sink
    return


###
def compute_Min_Cut_on_Goldberg_graph(graph_Goldberg, source_Goldberg, sink_Goldberg):
    cut_value, partition = nx.minimum_cut(graph_Goldberg, source_Goldberg, sink_Goldberg, capacity='weight')
    return cut_value, partition


###
def Goldberg_Algorithm(G, starting_min_density=None, starting_max_density=None):
    t_0 = current_milli_time()
    best_Vs = set()
    #
    graph_Goldberg, source_Goldberg, sink_Goldberg = create_Goldberg_graph(G)
    t_1 = current_milli_time()
    #
    min_node_degree = min(G.degree(weight='weight'), key=lambda x: x[-1])[-1]
    max_node_degree = max(G.degree(weight='weight'), key=lambda x: x[-1])[-1]
    # print("dsfg min_node_degree", min_node_degree)
    # print("dsfg max_node_degree", max_node_degree)
    #
    n = float(G.number_of_nodes())
    density_granularity = min_node_degree / (n * (n - 1.))
    #
    max_density = starting_max_density if starting_max_density is not None else max_node_degree / 2.
    min_density = starting_min_density if starting_min_density is not None else min_node_degree / 2.
    #
    max_number_of_iterations = 1 + math.ceil(math.log((max_density - min_density) / density_granularity) / math.log(2))
    #
    current_density = -1.
    num_iterations = 0
    while (num_iterations <= max_number_of_iterations) and ((max_density - min_density) >= density_granularity):
        #
        # print("graph_Goldberg ", len(graph_Goldberg), len(graph_Goldberg.edges()))
        #
        num_iterations += 1
        #
        current_density = (max_density + min_density) / 2.
        if (max_density == current_density) or (min_density == current_density):
            # print(" max_density==current_density", max_density == current_density)
            # print(" min_density==current_density", min_density == current_density)
            break
        #
        current_time = datetime.datetime.now()
        # print()
        # print(" Goldberg_Algorithm:", "num_iterations", num_iterations, "current_density", current_density,
        #       "current_time", current_time)
        # print(" max_density", max_density)
        # print(" min_density", min_density)
        # print(" density_granularity", density_granularity)
        # print(" max_number_of_iterations", max_number_of_iterations)
        # print(" num_iterations", num_iterations)
        # print()
        # print(" max_density==current_density", max_density == current_density)
        # print(" min_density==current_density", min_density == current_density)
        # print(" max_density - min_density", max_density - min_density)
        # print(" ((max_density - min_density) >= density_granularity)",
        #       ((max_density - min_density) >= density_granularity))
        #
        update_weights_on_Goldberg_graph(G, graph_Goldberg, current_density, sink_Goldberg)
        #
        cut_value, partition = compute_Min_Cut_on_Goldberg_graph(graph_Goldberg, source_Goldberg, sink_Goldberg)
        #
        S = partition[0]
        T = partition[1]
        if source_Goldberg in partition[1]:
            S = partition[1]
            T = partition[0]
        #
        if len(S) == 1:
            max_density = current_density
            # print("   max_density = current_density")
        else:
            min_density = current_density
            best_Vs = set(S)
            # print("   min_density = current_density")
        #
        # print()
        # print("cut_value: " + str(cut_value))
        # print("|best_Vs|: ", len(best_Vs) - 1)
        ### print("partition: " + str(partition))
        # print()
    #
    best_Vs.discard(source_Goldberg)
    t_2 = current_milli_time()
    #
    total_exec_time_msec = t_2 - t_0
    # total_exec_time_WITHOUT_GRAPH_CONSTRUCTION_msec = t_1 - t_0
    total_exec_time_WITHOUT_GRAPH_CONSTRUCTION_msec = t_2 - t_1
    #
    return best_Vs, total_exec_time_msec, total_exec_time_WITHOUT_GRAPH_CONSTRUCTION_msec

##################################################################
##################################################################
##################################################################
