import networkx as nx
import time
import copy
import pprint as pp

from spectral_algorithms import make_the_weighed_graph_connected
from spectral_algorithms import compute_main_eigenvector__scipy
from spectral_algorithms import compute_main_eigenvector_fairified_matrix__scipy__MORE_THAN_TWO_COLORS
from rounding import perform_paired_sweep_on_sorted_main_eigen_vector
from combinatorial_algorithms import Goldberg_Algorithm, compute_density_of_the_induced_subgraph

current_milli_time = lambda: int(round(time.time() * 1000))


##################################################################
##################################################################
##################################################################


def spectral_peeling(input_graph, map__node_id__color_id, using_fair_projection=True,
                     perform_peeling=True,
                     removing_policy="log",
                     smooth_landing=False,
                     at_most_one_node_for_each_color_in_the_final_solution=False):
    #
    best_solution = None
    #
    w_graph = copy.deepcopy(input_graph)
    #
    # add color attribute to each node in the graph.
    for c_node_id in w_graph:
        w_graph.nodes[c_node_id]["color"] = map__node_id__color_id[c_node_id]
    #
    set__all_color_ids = set(map__node_id__color_id.values())
    #
    iteration_number = 0
    must_continue_with_peeling = True
    while must_continue_with_peeling:
        iteration_number += 1
        #
        is_the_graph_connected = nx.is_connected(w_graph)
        #
        if not is_the_graph_connected:
            make_the_weighed_graph_connected(w_graph)
        #
        # create bijection node_id matrix index
        list_of_nodes_representing_the_mapping_index_node = list(w_graph.nodes())
        list_of_nodes_representing_the_mapping_index_node.sort()
        map__index__node_id = {}
        map__node_id__index = {}
        for c_index, c_node in enumerate(list_of_nodes_representing_the_mapping_index_node):
            map__index__node_id[c_index] = c_node
            map__node_id__index[c_node] = c_index
        #
        # create the adjacency matrix
        # scipy_sparse_adj_matrix = nx.convert_matrix.to_scipy_sparse_matrix(w_graph,
        #                           nodelist=list_of_nodes_representing_the_mapping_index_node,
        #                           dtype=None, weight='weight', format='csc')

        scipy_sparse_adj_matrix = \
            nx.convert_matrix.to_scipy_sparse_array(w_graph, nodelist=list_of_nodes_representing_the_mapping_index_node,
                                                    dtype=None, weight='weight', format='csc')

        ########################################################
        # Compute EigenValues and EigenVectors
        # of the projected matrix M.
        ########################################################
        map__index__color_id = {}
        for c_index in map__index__node_id:
            c_node_id = map__index__node_id[c_index]
            c_color_id = w_graph.nodes[c_node_id]["color"]
            map__index__color_id[c_index] = c_color_id
        t_0 = current_milli_time()
        #
        main_eve_M__w_power__scipy = None
        num_iterations_meve_M_sp = None
        absence_of_convergence = False
        try:
            if using_fair_projection:
                main_eve_M__w_power__scipy, num_iterations_meve_M_sp = \
                    compute_main_eigenvector_fairified_matrix__scipy__MORE_THAN_TWO_COLORS(
                        scipy_sparse_adj_matrix, map__index__color_id, max_number_of_iterations=100000,
                        tollerance=1e-08)
            else:
                main_eve_M__w_power__scipy, num_iterations_meve_M_sp = compute_main_eigenvector__scipy(
                    scipy_sparse_adj_matrix, max_number_of_iterations=100000, tollerance=1e-08)
        except:
            absence_of_convergence = True
        #
        t_1 = current_milli_time()
        total_time = t_1 - t_0
        #
        ########################################################

        ###########################################################
        # Conversion of the format of the MainEigenVector
        # of the projected matrix M.
        ###########################################################
        #
        main_eve_M = []
        if not absence_of_convergence:
            main_eve_M = main_eve_M__w_power__scipy.tolist()[0]
        else:
            for c_index in map__index__node_id:
                c_node_id = map__index__node_id[c_index]
                c_degree = input_graph.degree[c_node_id]
                main_eve_M = [c_degree] * len(map__index__color_id)
        #
        #
        ###########################################################
        # Sort indexes in the MainEigeinVector
        # of the projected matrix M
        # according to their value.
        ###########################################################
        #
        sorted_main_eigen_vector_M = [[index, value] for index, value in enumerate(main_eve_M)]
        sorted_main_eigen_vector_M.sort(key=lambda x: -x[1])
        #
        #
        # balance_level_M, map__color__sum_M = check_fairness(sorted_main_eigen_vector_M, map__index__color_id)

        ###########
        # SwEeP ###
        ###########
        list_of_solutions_from_Paired_Sweep_rounding_on_M = perform_paired_sweep_on_sorted_main_eigen_vector(
            w_graph,
            sorted_main_eigen_vector_M,
            map__index__node_id,
            map__index__color_id,
            set__all_color_ids,
            only_densest_fair_solutions=True,
            exactly_one_node_for_each_color_in_the_final_solution=True,
            at_most_one_node_for_each_color_in_the_final_solution=False)
        #
        if len(list_of_solutions_from_Paired_Sweep_rounding_on_M) == 0:
            must_continue_with_peeling = False
            break
        if (len(list_of_solutions_from_Paired_Sweep_rounding_on_M) == 1) and (
                list_of_solutions_from_Paired_Sweep_rounding_on_M[0] is None):
            must_continue_with_peeling = False
            break
        #
        best_solutions_from_Paired_Sweep_rounding_on_M = max(list_of_solutions_from_Paired_Sweep_rounding_on_M,
                                                             key=lambda x: (x[3], x[2]))
        #
        if at_most_one_node_for_each_color_in_the_final_solution:
            #
            subgraph_of_w_graph = w_graph.subgraph(best_solutions_from_Paired_Sweep_rounding_on_M[0])
            #print("subgraph_of_w_graph", subgraph_of_w_graph)
            best_Vs, total_exec_time_msec, total_exec_time_WITHOUT_GRAPH_CONSTRUCTION_msec = \
                Goldberg_Algorithm(subgraph_of_w_graph)
            c_best_solution = [frozenset(best_Vs),
                               best_solutions_from_Paired_Sweep_rounding_on_M[1],
                               compute_density_of_the_induced_subgraph(best_Vs, subgraph_of_w_graph),
                               None,
                               len(best_Vs)]
            #
            if (best_solution is None) or (c_best_solution[2] > best_solution[2]):
                best_solution = c_best_solution
        else:
            if (best_solution is None) or (best_solutions_from_Paired_Sweep_rounding_on_M[2] > best_solution[2]):
                best_solution = best_solutions_from_Paired_Sweep_rounding_on_M
        #
        print()
        print("Iteration number:", iteration_number)
        print("Density of the best solution found so far:", best_solution[-3])
        print("best_solution", best_solution)
        #
        order = best_solutions_from_Paired_Sweep_rounding_on_M[1]
        if order == 0:
            sorted_main_eigen_vector_M.sort(key=lambda x: (+1. * x[1], x[0]))
        if order == 1:
            sorted_main_eigen_vector_M.sort(key=lambda x: (-1. * x[1], x[0]))
        if order == 2:
            sorted_main_eigen_vector_M.sort(key=lambda x: (+1. * abs(x[1]), x[0]))
        if order == 3:
            sorted_main_eigen_vector_M.sort(key=lambda x: (-1. * abs(x[1]), x[0]))
        #
        #
        c_set__all_color_ids = set()
        map__color_id__num_indexes = {}
        map__index__color_id = {}
        map__color_id__top_index = {}
        map__color_id__bottom_index = {}
        for c_index, c_Meve_index_id__Meve_score in enumerate(sorted_main_eigen_vector_M):
            #
            c_Meve_index_id = c_Meve_index_id__Meve_score[0]
            c_node_id = map__index__node_id[c_Meve_index_id]
            #
            c_color = w_graph.nodes[c_node_id]["color"]
            #
            map__index__color_id[c_index] = c_color
            map__color_id__num_indexes[c_color] = map__color_id__num_indexes.get(c_color, 0) + 1
            #
            c_set__all_color_ids.add(c_color)
            #
        if not at_most_one_node_for_each_color_in_the_final_solution:
            if len(c_set__all_color_ids) < len(set__all_color_ids):
                must_continue_with_peeling = False
                break
        #
        #
        nodes_to_remove_in_this_iteration = 1
        if removing_policy == "log":
            nodes_to_remove_in_this_iteration = int(1 + len(sorted_main_eigen_vector_M) / 2)
            #
        if removing_policy == "sqrt":
            nodes_to_remove_in_this_iteration = int(len(sorted_main_eigen_vector_M) ** 0.5)
            #
        if removing_policy.isnumeric():
            nodes_to_remove_in_this_iteration = int(removing_policy)
            #
        #
        if smooth_landing:
            if len(sorted_main_eigen_vector_M) <= 2 * len(set__all_color_ids):
                nodes_to_remove_in_this_iteration = 1
        #
        list__indexes_to_remove = list()
        for j in range(len(sorted_main_eigen_vector_M)):
            #
            if len(list__indexes_to_remove) == nodes_to_remove_in_this_iteration:
                break
            #
            c_index = len(sorted_main_eigen_vector_M) - j - 1
            c_color = map__index__color_id[c_index]
            #
            if map__color_id__num_indexes[c_color] == 1:
                continue
            #
            list__indexes_to_remove.append(c_index)
            map__color_id__num_indexes[c_color] -= 1
        #
        #
        #########################
        # remove node from graph
        #########################
        #
        if not perform_peeling:
            must_continue_with_peeling = False
        #
        if len(list__indexes_to_remove) == 0:
            must_continue_with_peeling = False
        for c_index in list__indexes_to_remove:
            node_to_peel_out_id = map__index__node_id[sorted_main_eigen_vector_M[c_index][0]]
            w_graph.remove_node(node_to_peel_out_id)
        #
    #
    return [best_solution]
