from utilities import compute_balance_level


###
def compute_Lower_Bound_for_DSG_with_single_sweep_on_sorted_main_eigen_vector(graph,
                                                                              sorted_main_eigen_vector,
                                                                              map__index__node):
    #
    best_density = -1.
    #
    for order in [0, 1, 2, 3]:
        if order == 0:
            sorted_main_eigen_vector.sort(key=lambda x: (+1. * x[1], x[0]))
        if order == 1:
            sorted_main_eigen_vector.sort(key=lambda x: (-1. * x[1], x[0]))
        if order == 2:
            sorted_main_eigen_vector.sort(key=lambda x: (+1. * abs(x[1]), x[0]))
        if order == 3:
            sorted_main_eigen_vector.sort(key=lambda x: (-1. * abs(x[1]), x[0]))
        #
        #
        current_set_of_nodes = set()
        current_num_edges = 0
        current_sum_weights_on_edges = 0.
        for row_id, value in sorted_main_eigen_vector:
            #
            # print(" order_type", order, "#nodes", len(current_set_of_nodes))
            #
            c_node_id = map__index__node[row_id]
            #
            for c_neig_node_id in graph[c_node_id].keys():
                if (c_neig_node_id in current_set_of_nodes):
                    current_num_edges += 1
                    current_sum_weights_on_edges += graph[c_node_id][c_neig_node_id]["weight"]
                    #
            #
            current_set_of_nodes.add(c_node_id)
            #
            # compute density of induced subgraph
            num_edges_current_sub_graph = current_num_edges
            sum_weights_on_edges_current_sub_graph = current_sum_weights_on_edges
            num_nodes_current_sub_graph = len(current_set_of_nodes)
            # unweighted_density_current_sub_graph = float(num_edges_current_sub_graph) / num_nodes_current_sub_graph
            density_current_sub_graph = sum_weights_on_edges_current_sub_graph / num_nodes_current_sub_graph
            #
            best_density = density_current_sub_graph if density_current_sub_graph > best_density else best_density
            #
    return best_density


###
def perform_single_sweep_on_sorted_main_eigen_vector(graph, sorted_main_eigen_vector,
                                                     map__index__node,
                                                     map__index__color,
                                                     all_colors,
                                                     only_fair_solutions=False):
    integrated_sorted_main_eigen_vector = []
    integrated_sorted_main_eigen_vector__as__set = set()
    #
    for order in [0, 1, 2, 3]:
        if order == 0:
            sorted_main_eigen_vector.sort(key=lambda x: (+1. * x[1], x[0]))
        if order == 1:
            sorted_main_eigen_vector.sort(key=lambda x: (-1. * x[1], x[0]))
        if order == 2:
            sorted_main_eigen_vector.sort(key=lambda x: (+1. * abs(x[1]), x[0]))
        if order == 3:
            sorted_main_eigen_vector.sort(key=lambda x: (-1. * abs(x[1]), x[0]))
        #
        #
        map__color__num_nodes = {}
        for color in all_colors:
            map__color__num_nodes[color] = 0
        #
        current_set_of_nodes = set()
        current_num_edges = 0
        current_sum_weights_on_edges = 0.
        for row_id, value in sorted_main_eigen_vector:
            #
            # print(" order_type", order, "#nodes", len(current_set_of_nodes))
            #
            c_node_id = map__index__node[row_id]
            #
            for c_neig_node_id in graph[c_node_id].keys():
                if (c_neig_node_id in current_set_of_nodes):
                    current_num_edges += 1
                    current_sum_weights_on_edges += graph[c_node_id][c_neig_node_id]["weight"]
                    #
            #
            current_set_of_nodes.add(c_node_id)
            # current_list_of_nodes.append(c_node_id)
            #
            #
            # compute density of induced subgraph
            num_edges_current_sub_graph = current_num_edges
            sum_weights_on_edges_current_sub_graph = current_sum_weights_on_edges
            num_nodes_current_sub_graph = len(current_set_of_nodes)
            unweighted_density_current_sub_graph = float(num_edges_current_sub_graph) / num_nodes_current_sub_graph
            density_current_sub_graph = sum_weights_on_edges_current_sub_graph / num_nodes_current_sub_graph
            #
            # compute level of fairness of induced subgraph
            color_of_row_id = map__index__color[row_id]
            map__color__num_nodes[color_of_row_id] += 1
            level_of_ballance_current_sub_graph = compute_balance_level(map__color__num_nodes)
            # add density and level of balance of current_induced_subgraph to integrated_sorted_main_eigen_vector
            #
            if not only_fair_solutions or (level_of_ballance_current_sub_graph >= 0.9999):
                integrated_sorted_main_eigen_vector__as__set.add(
                    (frozenset(current_set_of_nodes), value, density_current_sub_graph,
                     level_of_ballance_current_sub_graph, num_nodes_current_sub_graph))
            #
    integrated_sorted_main_eigen_vector = list(integrated_sorted_main_eigen_vector__as__set)
    return integrated_sorted_main_eigen_vector


###
def perform_paired_sweep_on_sorted_main_eigen_vector(graph,
                                                     sorted_main_eigen_vector,
                                                     map__index__node,
                                                     map__index__color,
                                                     all_colors,
                                                     only_densest_fair_solutions=True,
                                                     exactly_one_node_for_each_color_in_the_final_solution=False,
                                                     at_most_one_node_for_each_color_in_the_final_solution=False):
    #
    integrated_sorted_main_eigen_vector = []
    max_DENSITY_encountered_so_far = -1.
    best_SOLUTION_as_output_record_encountered_so_far = None
    #
    for order in [0, 1, 2, 3]:
        if order == 0:
            sorted_main_eigen_vector.sort(key=lambda x: (+1. * x[1], x[0]))
        if order == 1:
            sorted_main_eigen_vector.sort(key=lambda x: (-1. * x[1], x[0]))
        if order == 2:
            sorted_main_eigen_vector.sort(key=lambda x: (+1. * abs(x[1]), x[0]))
        if order == 3:
            sorted_main_eigen_vector.sort(key=lambda x: (-1. * abs(x[1]), x[0]))
        #
        #
        if at_most_one_node_for_each_color_in_the_final_solution:
            c_best_SOLUTION_as_output_record_encountered_so_far = perform__sweep_on_COLORS(sorted_main_eigen_vector,
                                                                                           order, graph,
                                                                                           map__index__node,
                                                                                           map__index__color,
                                                                                           all_colors)
            #
            density_current_sub_graph = c_best_SOLUTION_as_output_record_encountered_so_far[2]
            #
            if (density_current_sub_graph > max_DENSITY_encountered_so_far):
                max_DENSITY_encountered_so_far = density_current_sub_graph
                best_SOLUTION_as_output_record_encountered_so_far = c_best_SOLUTION_as_output_record_encountered_so_far
            continue
        #
        #
        map__color__MEve = {}
        for color in all_colors:
            map__color__MEve[color] = []
        #
        for row_id, value in sorted_main_eigen_vector:
            color_of_row_id = map__index__color[row_id]
            #
            map__color__MEve[color_of_row_id].append((row_id, value))
        #
        #
        map__color__num_nodes = {}
        for color in all_colors:
            map__color__num_nodes[color] = 0
        #
        #
        current_set_of_nodes = set()
        current_num_edges = 0
        current_sum_weights_on_edges = 0.
        index = -1
        to_break_for_only_C_nodes_in_solution = False
        while index < len(sorted_main_eigen_vector):
            index += 1
            #
            for c_color in map__color__MEve:
                if index >= len(map__color__MEve[c_color]):
                    index = len(sorted_main_eigen_vector)
                    break
                #
                row_id = map__color__MEve[c_color][index][0]
                value = map__color__MEve[c_color][index][1]
                #
                #
                c_node_id = map__index__node[row_id]
                #
                for c_neig_node_id in graph[c_node_id].keys():
                    if (c_neig_node_id in current_set_of_nodes):
                        current_num_edges += 1
                        current_sum_weights_on_edges += graph[c_node_id][c_neig_node_id]["weight"]
                #
                current_set_of_nodes.add(c_node_id)
                #
                # print("paired order_type", order, "#nodes", len(current_set_of_nodes))
                #
                #
                # compute density of induced subgraph
                num_edges_current_sub_graph = current_num_edges
                sum_weights_on_edges_current_sub_graph = current_sum_weights_on_edges
                num_nodes_current_sub_graph = len(current_set_of_nodes)
                unweighted_density_current_sub_graph = float(num_edges_current_sub_graph) / num_nodes_current_sub_graph
                density_current_sub_graph = sum_weights_on_edges_current_sub_graph / num_nodes_current_sub_graph
                #
                # compute level of fairness of induced subgraph
                color_of_row_id = map__index__color[row_id]
                map__color__num_nodes[color_of_row_id] += 1
                level_of_fairness_current_sub_graph = compute_balance_level(map__color__num_nodes)
                if (level_of_fairness_current_sub_graph < 1.):
                    continue
                #
                # add density and level of fairness of current_induced_subgraph to integrated_sorted_main_eigen_vector
                if (density_current_sub_graph > max_DENSITY_encountered_so_far):
                    max_DENSITY_encountered_so_far = density_current_sub_graph
                    best_SOLUTION_as_output_record_encountered_so_far = [frozenset(current_set_of_nodes), order,
                                                                         density_current_sub_graph,
                                                                         level_of_fairness_current_sub_graph,
                                                                         num_nodes_current_sub_graph]
                #
                #
                if exactly_one_node_for_each_color_in_the_final_solution:
                    to_break_for_only_C_nodes_in_solution = True
                    break
                #
                # check if exactly one rep for each color
                # to_break_ = False
                # num_colss___ = 0
                # for col_, num_nodes_ in map__color__num_nodes.items():
                #    if num_nodes_ == 1:
                #        num_colss___ += 1
                # if len(all_colors) == num_colss___:
                #    break
                #
                #
            if exactly_one_node_for_each_color_in_the_final_solution:
                if to_break_for_only_C_nodes_in_solution:
                    break
    #
    integrated_sorted_main_eigen_vector.append(best_SOLUTION_as_output_record_encountered_so_far)
    final_integrated_sorted_main_eigen_vector = integrated_sorted_main_eigen_vector
    #
    #
    #
    # final_integrated_sorted_main_eigen_vector = []
    # for c_set_of_nodes_sub_graph, c_value, c_density_current_sub_graph, c_level_of_fairness_current_sub_graph, c_num_nodes_current_sub_graph in integrated_sorted_main_eigen_vector:
    #    if (c_density_current_sub_graph >= max_DENSITY_encountered_so_far):
    #        final_integrated_sorted_main_eigen_vector.append(
    #            [c_set_of_nodes_sub_graph, c_value, c_density_current_sub_graph, c_level_of_fairness_current_sub_graph,
    #             c_num_nodes_current_sub_graph])
    #
    return final_integrated_sorted_main_eigen_vector


def perform__sweep_on_COLORS(sorted_main_eigen_vector, order, graph,
                             map__index__node,
                             map__index__color,
                             all_colors):
    #
    current_set_of_nodes = set()
    density_current_sub_graph = None
    level_of_fairness_current_sub_graph = None
    num_nodes_current_sub_graph = None
    best_SOLUTION_as_output_record_encountered_so_far = [frozenset(current_set_of_nodes), order,
                                                         density_current_sub_graph,
                                                         level_of_fairness_current_sub_graph,
                                                         num_nodes_current_sub_graph]
    #
    set__colors_encountered_so_far = set()
    list__top_row_id__color = []
    for row_id, value in sorted_main_eigen_vector:
        c_color = map__index__color[row_id]
        if c_color in set__colors_encountered_so_far:
            continue
        #
        list__top_row_id__color.append((row_id, c_color))
        #
        set__colors_encountered_so_far.add(c_color)
        #
        if len(set__colors_encountered_so_far) == len(all_colors):
            break
        #
    #
    # Let's SWEEP
    current_set_of_nodes = set()
    current_num_edges = 0
    current_sum_weights_on_edges = 0.
    num_edges_current_sub_graph = 0
    sum_weights_on_edges_current_sub_graph = 0
    num_nodes_current_sub_graph = 0
    unweighted_density_current_sub_graph = 0
    density_current_sub_graph = 0
    max_DENSITY_encountered_so_far = -1.
    #
    map__color__num_nodes = {}
    for color in all_colors:
        map__color__num_nodes[color] = 0
    #
    for row_id, c_color in list__top_row_id__color:
        #
        c_node_id = map__index__node[row_id]
        #
        for c_neig_node_id in graph[c_node_id].keys():
            if (c_neig_node_id in current_set_of_nodes):
                current_num_edges += 1
                current_sum_weights_on_edges += graph[c_node_id][c_neig_node_id]["weight"]
        #
        current_set_of_nodes.add(c_node_id)
        #
        # compute density of induced subgraph
        num_edges_current_sub_graph = current_num_edges
        sum_weights_on_edges_current_sub_graph = current_sum_weights_on_edges
        num_nodes_current_sub_graph = len(current_set_of_nodes)
        unweighted_density_current_sub_graph = float(num_edges_current_sub_graph) / num_nodes_current_sub_graph
        density_current_sub_graph = sum_weights_on_edges_current_sub_graph / num_nodes_current_sub_graph
        #
        # compute level of fairness of induced subgraph
        color_of_row_id = map__index__color[row_id]
        map__color__num_nodes[color_of_row_id] += 1
        level_of_fairness_current_sub_graph = compute_balance_level(map__color__num_nodes)
        #
        # add density and level of fairness of current_induced_subgraph to integrated_sorted_main_eigen_vector
        if (density_current_sub_graph > max_DENSITY_encountered_so_far):
            max_DENSITY_encountered_so_far = density_current_sub_graph
            best_SOLUTION_as_output_record_encountered_so_far = [frozenset(current_set_of_nodes), order,
                                                                 density_current_sub_graph,
                                                                 level_of_fairness_current_sub_graph,
                                                                 num_nodes_current_sub_graph]
        #
    #
    return best_SOLUTION_as_output_record_encountered_so_far
