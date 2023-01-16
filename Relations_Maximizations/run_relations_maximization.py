import sys

if len(sys.argv) != 2:
    print()
    print("Wrong number of arguments: provided " + str(len(sys.argv) - 1) + " required 1.")
    print()
    print("Correct usage:")
    print(" python3 run_relations_maximization.py graph_input_file_name")
    print()
    exit(-1)

graph_file_name = sys.argv[1]

import networkx as nx
from networkx import *
import pprint as pp
import csv
import time
from datetime import datetime

from spectral_peeling_algorithms import spectral_peeling

current_milli_time = lambda: int(round(time.time() * 1000))

c_time = str(datetime.now()).replace("-", "_").replace(" ", "__").replace(":", "_").replace(".", "_")
output_file_handler = open("output__" + str(c_time) + ".tsv", "w")
csv_writer = csv.writer(output_file_handler, delimiter='\t', quoting=csv.QUOTE_NONE)
csv_writer.writerow(
    ["input graph name",
     "total number of genes", "total number of loci",
     "algorithm",
     "density of the solution",
     "number of genes in the solution",
     "set of genes in the solution"])
output_file_handler.flush()
algo_name = "Relations Maximization Algorithm"
#
#
########################
# Load graph from disk #
########################
#
t_0 = current_milli_time()
graph = nx.Graph()
map__product_id__node_id = {}
max_node_id_so_far = -1
#
map__category__set_node_ids = {}
#
input_file = open(graph_file_name, "r")
csv_reader = csv.reader(input_file, delimiter='\t', quotechar='"', quoting=csv.QUOTE_NONE)
header = csv_reader.__next__()
is_the_graph_a_weighted_graph = False
num_records = 0
for c_record in csv_reader:
    num_records += 1
    #
    c_product_id_1 = c_record[0]
    c_product_id_2 = c_record[1]
    c_category_1 = c_record[2]
    c_category_2 = c_record[3]
    #
    c_weight = 1.
    if len(c_record) > 4:
        c_weight = float(c_record[4])
        is_the_graph_a_weighted_graph = True
    #
    c_node_id_1 = map__product_id__node_id.get(c_product_id_1, -1)
    if (c_node_id_1 == -1):
        max_node_id_so_far += 1
        c_node_id_1 = max_node_id_so_far
        map__product_id__node_id[c_product_id_1] = c_node_id_1
    #
    c_node_id_2 = map__product_id__node_id.get(c_product_id_2, -1)
    if (c_node_id_2 == -1):
        max_node_id_so_far += 1
        c_node_id_2 = max_node_id_so_far
        map__product_id__node_id[c_product_id_2] = c_node_id_2
    #
    if (c_category_1 not in map__category__set_node_ids):
        map__category__set_node_ids[c_category_1] = set()
    map__category__set_node_ids[c_category_1].add(c_node_id_1)
    #
    if (c_category_2 not in map__category__set_node_ids):
        map__category__set_node_ids[c_category_2] = set()
    map__category__set_node_ids[c_category_2].add(c_node_id_2)
    #
    #
    graph.add_edge(c_node_id_1, c_node_id_2, weight=c_weight)
    #
    #
#
map__inner_node_id__outer_node_id = {}
for product_id, node_id in map__product_id__node_id.items():
    map__inner_node_id__outer_node_id[node_id] = product_id
#
map__category__color = {}
map__color__category = {}
for c_color, c_category in enumerate(map__category__set_node_ids):
    map__category__color[c_category] = c_color
    map__color__category[c_color] = c_category
#
map__node__color = {}
for c_category, c_set_node_ids in map__category__set_node_ids.items():
    c_color = map__category__color[c_category]
    for c_node_id in c_set_node_ids:
        map__node__color[c_node_id] = c_color
#
all_colors = set(map__color__category.keys())
input_file.close()
t_1 = current_milli_time()
total_time = t_1 - t_0
#
print()
print("Graph file name:", graph_file_name)
print("Graph with " + str(len(graph)) + " nodes and " + str(nx.number_of_edges(graph)) + " edges")
print("Graph creation:", total_time, "msec")
print()

########################################
# Run Relations-Maximization Algorithm #
########################################
#
t_0 = current_milli_time()
#
list_of_solutions_from_Single_Sweep_rounding_on_M = spectral_peeling(graph, map__node__color,
                                                                     using_fair_projection=True)
#
t_1 = current_milli_time()
total_time = t_1 - t_0
print()
print(" " + algo_name + " total_time in msec: ", total_time, "msec")
print()
#########################
#########################
#########################
density_of_solution = list_of_solutions_from_Single_Sweep_rounding_on_M[0][-3]
set__original_node_id_in_solution = set()
for inner_node_id in list_of_solutions_from_Single_Sweep_rounding_on_M[0][0]:
    set__original_node_id_in_solution.add(map__inner_node_id__outer_node_id[inner_node_id])
#
row = [graph_file_name, len(graph.nodes), len(all_colors), algo_name, density_of_solution,
       len(set__original_node_id_in_solution), set__original_node_id_in_solution]
csv_writer.writerow(row)
output_file_handler.flush()
output_file_handler.close()

print()
