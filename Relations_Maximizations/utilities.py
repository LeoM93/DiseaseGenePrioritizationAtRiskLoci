import numpy as np
import math

from numpy.linalg import pinv
from scipy.sparse.linalg import eigsh


####################################################################################
####################################################################################
####################################################################################

###
def build_F(number_of_colors, map__index__color):
    f = []
    #
    sorted_list_of_indexes = list(map__index__color.keys())
    sorted_list_of_indexes.sort()
    num_indexes = len(sorted_list_of_indexes)
    #
    norm = 0.
    for c_color_id in range(0, number_of_colors - 1):
        c_f_row = list()
        for c_index in sorted_list_of_indexes:
            if map__index__color[c_index] == 0:
                c_f_row.append(1.)
                norm += 1.
            elif map__index__color[c_index] == c_color_id + 1:
                c_f_row.append(-1.)
                norm += 1.
            else:
                c_f_row.append(0)
        f.append(c_f_row)
    #
    norm = math.sqrt(norm)
    #
    for c_row_index in range(len(f)):
        for c_col_index in range(len(f[c_row_index])):
            f[c_row_index][c_col_index] /= norm
    #
    F = np.array(f).T
    #
    return F


###
def build_F_in_the_binary_color_scenario(map__index__color, color_factor_1=1.):
    f = []
    sorted_list_of_indexes = list(map__index__color.keys())
    sorted_list_of_indexes.sort()
    norm = 0.
    for index in sorted_list_of_indexes:
        if map__index__color[index] == 0:
            f.append(1.)
            norm += 1.
        else:
            f.append(-1. * color_factor_1)
            norm += (1. * color_factor_1) ** 2
    norm = math.sqrt(norm)
    # print("norm: " + str(norm))
    # print("norm**2: " + str(norm**2))
    for i in range(len(f)):
        f[i] /= norm
    F = np.array(f).reshape(len(f), 1)
    return F


def create_M(F, adj_matrix):
    F_x_Fpinv = F @ pinv(F)
    I = np.identity(len(F))
    I_minus_F_x_Fpinv = (I - F_x_Fpinv)
    B = I_minus_F_x_Fpinv @ adj_matrix
    M = B @ I_minus_F_x_Fpinv
    return M


####################################################################################
####################################################################################
####################################################################################


###
def check_fairness(sorted_main_eigen_vector, map__index__color):
    map__color__sum = {}
    #
    fairness_level = 0.
    for index, value in sorted_main_eigen_vector:
        c_color = map__index__color[index]
        #
        if c_color == 0:
            fairness_level += value
        if c_color == 1:
            fairness_level -= value
        #
        if c_color not in map__color__sum:
            map__color__sum[c_color] = 0.
        map__color__sum[c_color] += value
    #
    return fairness_level, map__color__sum


###
def compute_balance_level(map__color__num_nodes):
    #
    max_num_nodes = 0.
    min_num_nodes = float("+inf")
    for c_color, c_num_nodes in map__color__num_nodes.items():
        if c_num_nodes < min_num_nodes:
            min_num_nodes = c_num_nodes
        if c_num_nodes > max_num_nodes:
            max_num_nodes = c_num_nodes
        if c_num_nodes == 0:
            return 0.
    #
    balance_level = min_num_nodes / max_num_nodes
    #
    return balance_level


#
def compute_balance_level__ONLY_FOR_TWO_COLORS(map__color__num_nodes):
    if len(map__color__num_nodes) != 2:
        return -1.
    if map__color__num_nodes[0] == 0 or map__color__num_nodes[1] == 0:
        return 0.
    balance_level = min([float(map__color__num_nodes[0]) / map__color__num_nodes[1],
                         float(map__color__num_nodes[1]) / map__color__num_nodes[0]])
    return balance_level


###
def compute_balance_of_the_graph(graph, map__node__color, all_colors):
    map__color__num_nodes = {}
    for color in all_colors:
        map__color__num_nodes[color] = 0
    #
    for c_node in graph.nodes:
        c_color = map__node__color[c_node]
        # print(str(c_color))
        map__color__num_nodes[c_color] += 1
    #
    balance_level = compute_balance_level(map__color__num_nodes)
    #
    return balance_level


####################################################################################
####################################################################################
####################################################################################


###
def get_set_of_nodes_of_best_densest_subgraph(integrated_sorted_main_eigen_vector=[],
                                              min_acceptable_level_of_fairness=1.):
    best_set_of_nodes = set()
    best_density = 0.
    level_of_fairness_of_best_set_of_nodes = 0.
    num_nodes_of_best_set_of_nodes = 0
    #
    for list_of_nodes, value, density, level_of_fairness, num_nodes in integrated_sorted_main_eigen_vector:
        set_of_nodes = set(list_of_nodes)
        if level_of_fairness < min_acceptable_level_of_fairness:
            continue
        if density >= best_density:
            best_set_of_nodes = set_of_nodes
            best_density = density
            level_of_fairness_of_best_set_of_nodes = level_of_fairness
            num_nodes_of_best_set_of_nodes = num_nodes
    #
    return best_set_of_nodes, best_density, level_of_fairness_of_best_set_of_nodes, num_nodes_of_best_set_of_nodes

####################################################################################
####################################################################################
####################################################################################
