import sys
import time
import numpy as np

import scipy as sp
from scipy.sparse.linalg import eigsh

from utilities import build_F
from utilities import build_F_in_the_binary_color_scenario

##################################################################
##################################################################
##################################################################

current_milli_time = lambda: int(round(time.time() * 1000))


def compute_UPPER_BOUND_on_Density_of_DSG(main_eve_A__w_power__scipy, csc_A):
    UPPER_BOUND_on_Density_of_DSG = 0.
    #
    numerator = (main_eve_A__w_power__scipy @ csc_A) @ main_eve_A__w_power__scipy.transpose()
    denomerator = main_eve_A__w_power__scipy @ main_eve_A__w_power__scipy.transpose()
    print("numerator  =", numerator)
    print("denomerator=", denomerator)
    UPPER_BOUND_on_Density_of_DSG = numerator / denomerator
    UPPER_BOUND_on_Density_of_DSG = UPPER_BOUND_on_Density_of_DSG[0, 0]
    #
    return UPPER_BOUND_on_Density_of_DSG


#
##
###
def compute_main_eigenvector_fairified_matrix__scipy__MORE_THAN_TWO_COLORS(csc_A, map__index__color,
                                                                           max_number_of_iterations=10000,
                                                                           tollerance=1e-08):
    main_eve = None
    #
    x = np.random.rand(csc_A.shape[1])
    # x = [1.] * csc_A.shape[1]
    x = sp.asarray(np.array(x).reshape(1, csc_A.shape[1]))
    x = x / sp.linalg.norm(x)
    #
    t_0 = current_milli_time()
    number_of_colors = len({color_id for index, color_id in map__index__color.items()})
    # F = sp.asarray(build_F(number_of_colors, map__index__color))
    F = build_F(number_of_colors, map__index__color)
    t_1 = current_milli_time()
    total_time = t_1 - t_0
    # print("F creation (to_scipy_sparse_matrix)", total_time, "msec")
    #
    t_0 = current_milli_time()
    # F_pinv = sp.asarray(np.linalg.pinv(F))
    F_pinv = np.linalg.pinv(F)
    t_1 = current_milli_time()
    total_time = t_1 - t_0
    # print("F_pinv creation (to_scipy_sparse_matrix)", total_time, "msec")
    #
    #
    previous_dist = float("+inf")
    dist = 0.
    num_iterations = 0
    epsilon = sys.float_info.min
    next_x = None
    for _ in range(max_number_of_iterations):
        num_iterations += 1
        #
        # print("num_iterations", num_iterations)
        #
        y = x @ F
        # y = y / sp.linalg.norm(y)
        #
        z = y @ F_pinv
        # z = z / sp.linalg.norm(z)
        #
        w = x - z
        # w = w / sp.linalg.norm(w)
        #
        t = w @ csc_A
        # t = t / sp.linalg.norm(t)
        #
        y = t @ F
        # y = y / sp.linalg.norm(y)
        #
        z = y @ F_pinv
        # z = z / sp.linalg.norm(z)
        next_x = t - z
        #
        #
        next_x = next_x / sp.linalg.norm(next_x)
        #
        dist = sp.linalg.norm(next_x - x)
        #
        # print()
        # print("dist                   =", dist)
        # print("previous_dist          =", previous_dist)
        # print("previous_dist + epsilon=", previous_dist + epsilon)
        # print("(dist - previous_dist) =", (dist - previous_dist))
        # print("|dist - previous_dist| =", abs(dist - previous_dist))
        # print("epsilon                =", epsilon)
        # print("((dist - previous_dist) <= epsilon) = ", ((dist - previous_dist) <= epsilon))
        # print("(|dist - previous_dist| <= epsilon) = ", (abs(dist - previous_dist) <= epsilon))
        #
        if dist <= tollerance:
            break
        if (abs(dist - previous_dist) <= epsilon):
            break
        # if (abs(dist - previous_dist) == 0.):
        # if (dist - previous_dist == 0.):
        #    break
        #
        x = next_x
        previous_dist = dist
        #
    #
    if (dist > tollerance) and not (abs(dist - previous_dist) <= epsilon):
        raise Exception('The provided maximum number of iterations has been reached.')
    #
    main_eve = next_x
    return main_eve, num_iterations


#
##
###
def compute_main_eigenvector__scipy(csc_A, max_number_of_iterations=100000, tollerance=1e-08):
    main_eve = None
    #
    x = np.random.rand(csc_A.shape[1])
    x = sp.asarray(np.array(x).reshape(1, csc_A.shape[1]))
    #
    dist = 0.
    num_iterations = 0
    next_x = None
    for _ in range(max_number_of_iterations):
        num_iterations += 1
        #
        # print("num_iterations", num_iterations)
        #
        next_x = x @ csc_A
        #
        next_x = next_x / sp.linalg.norm(next_x)
        #
        dist = sp.linalg.norm(next_x - x)
        #
        # print("csc_A", csc_A)
        # print("next_x", next_x)
        # print("num_iterations", num_iterations, "dist", dist)
        if dist <= tollerance:
            break
        #
        x = next_x
        #
        # print(x)
        #
    #
    if dist > tollerance:
        raise Exception('The provided maximum number of iterations has been reached.')
    #
    main_eve = next_x
    return main_eve, num_iterations


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################


def compute_main_eigenvector_fairified_matrix__scipy(csc_A, map__index__color,
                                                     max_number_of_iterations=10000,
                                                     tollerance=1e-06):
    main_eve = None
    #
    x = np.random.rand(csc_A.shape[1])
    x = sp.asarray(np.array(x).reshape(1, csc_A.shape[1]))
    #
    t_0 = current_milli_time()
    F = sp.asarray(build_F_in_the_binary_color_scenario(map__index__color))
    t_1 = current_milli_time()
    total_time = t_1 - t_0
    print("F creation (to_scipy_sparse_matrix)", total_time, "msec")
    #
    t_0 = current_milli_time()
    F_pinv = sp.asarray(np.linalg.pinv(F))
    # F_pinv = sp.linalg.pinv(F)
    t_1 = current_milli_time()
    total_time = t_1 - t_0
    print("F_pinv creation (to_scipy_sparse_matrix)", total_time, "msec")
    #
    #
    dist = 0.
    num_iterations = 0
    next_x = None
    for _ in range(max_number_of_iterations):
        num_iterations += 1
        #
        # print("num_iterations", num_iterations)
        #
        y = x @ F
        z = y @ F_pinv
        w = x - z
        #
        t = w @ csc_A
        #
        y = t @ F
        z = y @ F_pinv
        next_x = t - z
        #
        #
        next_x = next_x / sp.linalg.norm(next_x)
        #
        dist = sp.linalg.norm(next_x - x)
        if dist <= tollerance:
            break
        #
        x = next_x
        #
    #
    if dist > tollerance:
        raise Exception('The provided maximum number of iterations has been reached.')
    #
    main_eve = next_x
    return main_eve, num_iterations


def compute_main_eigenvector__numpy(A, max_number_of_iterations=10000, tollerance=1e-06):
    main_eve = None
    #
    b_k = np.random.rand(A.shape[1])
    next_unnormalized__b_k = None
    norm_of_next_unnormalized__b_k = 0.
    next__b_k = None
    dist = 0.
    #
    num_iterations = 0
    for _ in range(max_number_of_iterations):
        num_iterations += 1
        #
        # next_unnormalized__b_k = np.dot(A, b_k)
        next_unnormalized__b_k = np.dot(b_k, A)
        norm_of_next_unnormalized__b_k = np.linalg.norm(next_unnormalized__b_k)
        next__b_k = next_unnormalized__b_k / norm_of_next_unnormalized__b_k
        #
        dist = np.linalg.norm(next__b_k - b_k)
        if dist <= tollerance:
            break
        #
        b_k = next__b_k
        #
    #
    if dist > tollerance:
        raise Exception('The provided maximum number of iterations has been reached.')
    #
    main_eve = next__b_k
    return main_eve, num_iterations


def make_the_weighed_graph_connected(w_graph):
    #
    # print("Add epsilon-weighted-edges to make the graph complete and connected.")
    num_nodes_in_the_original_graph = w_graph.number_of_nodes()
    min_edge_weight = 10 ** -10
    epsilon_weight = min_edge_weight * (1. / (num_nodes_in_the_original_graph * (num_nodes_in_the_original_graph)))
    #
    set__node_id_in_input_graph = set(w_graph.nodes())
    for c_node_id in w_graph:
        #
        for c_new_neigh_node_id in w_graph:
            if c_new_neigh_node_id in w_graph[c_node_id]:
                continue
            #
            w_graph.add_edge(c_node_id, c_new_neigh_node_id, weight=epsilon_weight)
            #
    print()
    # is_NOW_the_input_graph_connected = nx.is_connected(w_graph)
    # print(" is_NOW_the_input_graph_connected? " + str(is_NOW_the_input_graph_connected))
    #
    return
