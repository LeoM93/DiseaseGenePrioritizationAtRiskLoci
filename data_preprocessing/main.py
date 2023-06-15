import sys
import os
import csv
import networkx as nx

from algorithm.protein_ontology_feature_extractor import ProteinOntologyFeatureExtractor

class Main():
	
	def __init__(
		self,
		network_file_path,
		gene_to_locus_file_path,
		output_file_path,

		):

		self.network_file_path = network_file_path
		self.network_dir_path = [item for index,item in enumerate(network_file_path.split("/")) if index != len(network_file_path.split("/")) -1 ]
		
		
		self.gene_to_locus_file_path = gene_to_locus_file_path
		self.weighted_network_resnik_sim = "/".join(self.network_dir_path) + "/biologically_weighted_"+ network_file_path.split("/")[-1]
		self.output_file_path = output_file_path

		print(self.weighted_network_resnik_sim)

	def __load_map_gene_locus__(self,delimiter = "\t"):

		map__gene_name__locus = {}

		with open(self.gene_to_locus_file_path,'r') as fp:
			csv_reader = csv.reader(fp,delimiter = delimiter)
			for row in csv_reader:
				gene_name = row[0]
				locus = row[1]

				map__gene_name__locus[gene_name] = locus

		return map__gene_name__locus



	def __get_network_weight__(
		self,

		map__gene_name__locus,
		regulation_costant = 1.0,
		only_reliable_info = True,
		database = "GO",
		interested_score = "resnik_similarity_biological_process",
		):
						
		fe = ProteinOntologyFeatureExtractor(
				ontology_file_path = "/Users/leonardomartini/Documents/network_medicine/Data/db/annotated_db/GO_10_02_2021.txt",
				ontology_graph_file_path = "/Users/leonardomartini/Documents/network_medicine/Data/db/annotated_db/ontology/go_graph.txt",
				network_file_path = self.network_file_path,
				output_file_path = self.weighted_network_resnik_sim,
				similarity_score_def = "resnik_similarity"
				
		)
		
		fe.run()

		

		csv_reader = csv.reader(open(self.weighted_network_resnik_sim,"r"),delimiter = "\t")
		weighted_network = [["gene_name_1","gene_name_2","snp_gene_1", "snp_gene_2", "weight"]]
		column_index = -1
			
		for index, row in enumerate(csv_reader):

			if index == 0:
				for col_id, col_name in enumerate(row):
					if col_name == interested_score:
						column_index = col_id		
				continue

			node_1 = row[0]
			node_2 = row[1]

			if node_1 not in map__gene_name__locus or node_2 not in map__gene_name__locus:
				continue

			locus_1 = map__gene_name__locus[node_1]
			locus_2 = map__gene_name__locus[node_2]

			if locus_1 == locus_2:
				continue

			if column_index != -1:
				interedsted_score = float(row[column_index])
					
				if interedsted_score > 0.0 and only_reliable_info:
			
					weighted_network.append([node_1,node_2,locus_1,locus_2,interedsted_score])


				elif not only_reliable_info:

					weighted_network.append([node_1,node_2,locus_1,locus_2,interedsted_score + regulation_costant])

			else:
				interedsted_score = 1.0
				weighted_network.append([node_1,node_2,locus_1,locus_2,interedsted_score])

		csv_writer = csv.writer(open(self.output_file_path, "w"),delimiter = "\t")
		csv_writer.writerows(weighted_network)

	def run(self,):

		print()
		print("Loading Unbiased DB...")
		map__gene_name__locus = self.__load_map_gene_locus__()

		print()
		print("Loading and Weighting Co-Regulation Network...")
		self.__get_network_weight__(map__gene_name__locus)


if __name__ == '__main__':


	if len(sys.argv) != 4:
	    print()
	    print("Wrong number of arguments: provided " + str(len(sys.argv) - 1) + " required 1.")
	    print()
	    print("Correct usage:")
	    print(" python3 main.py graph_input_file_path output_file_path gene_to_locus_file_path")
	    print()
	    exit(-1)
	
	m = Main(

		network_file_path = sys.argv[1],
		gene_to_locus_file_path = sys.argv[3],
		output_file_path = sys.argv[2],
			
	)

	m.run()

