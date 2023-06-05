import os
import subprocess
import csv
import scipy.stats as sps

import pandas as pd
import numpy as np

from data_preprocessing.snp_gene_creation.create_bipartite_graph import SNPGeneGraph
class AlgorithmWrapper():
	
	def __init__(self, 
		disease_dir_path,

	):
		self.disease_genes_dir_path = disease_dir_path + "output_vegas_2/"
		self.GWAS_dir = disease_dir_path + "GWAS/"
		self.disease_dir_path = disease_dir_path
		
		self.input_dir_path = disease_dir_path + "seed/"
		self.ensembl_db = "datasets/curated_db/mart_export.txt"
		
		self.algorithm_dir_path = self.disease_dir_path + "algorithms/"

		if not os.path.exists(self.algorithm_dir_path):
			os.makedirs(self.algorithm_dir_path)
		
		self.disease_seed_set_file_paths = [self.disease_genes_dir_path + file for file in os.listdir(self.disease_genes_dir_path) if file[0] != "."]


	def __load_gene_in_2MB_windows__(self,file_path):
		
		bg = SNPGeneGraph(
				rs_id_position_file_path = file_path + ".tsv",
				gene_overlap_2mb_windows_file_path = "datasets/GWAS/refSeqGenes_hg19.txt",
				output_file_path = ""
				)
			
		self.gene_set_2_MB = set()
		
		for g in bg.run():
			if g in self.map__gene__ensembl_id:
				if self.map__gene__ensembl_id[g] in self.V:
					self.gene_set_2_MB.add(self.map__gene__ensembl_id[g])

	
	
	def __load_node_PPI_network__(self, network_file_path = "/Users/leonardomartini/Documents/network_medicine/DiseaseGenePrioritizationAtRiskLoci/datasets/networks/STRING_PPI.tsv"):
		csv_reader = csv.reader(open(network_file_path,"r"),delimiter = "\t")
		self.V = set()
		
		for index, row in enumerate(csv_reader):
			if index == 0:
				continue
			self.V.add(row[0])
			self.V.add(row[1])
	


	
	def __load_ensembl_db__(self,):
		self.map__gene__ensembl_id = {}
		self. map__ensembl_id__gene = {}
		csv_reader = csv.reader(open(self.ensembl_db , "r"),delimiter = "\t")
		
		for index, row in enumerate(csv_reader):

			if index == 0:
				continue

			gene_name = row[1]
			ensembl_id = row[0]

			self.map__gene__ensembl_id[gene_name] = ensembl_id
			self.map__ensembl_id__gene[ensembl_id] = gene_name
	

	def __compute_input_data_for_network_based_approach__(self,):
		file_paths = []
		for file in self.disease_seed_set_file_paths:
			csv_reader = csv.reader(open(file,"r"),delimiter = " ")
			scored_gene_set = {}
			seed_set = []
			
			self.__load_gene_in_2MB_windows__(self.GWAS_dir + file.split("/")[-1])

			for index,  row in enumerate(csv_reader):
				if row[1] in self.map__gene__ensembl_id:	
					scored_gene_set[self.map__gene__ensembl_id[row[1]]] =  float(row[7])
			
			for v in self.gene_set_2_MB:
				if v in scored_gene_set:
					seed_set.append([v, scored_gene_set[v]])
				else:
					seed_set.append([v, 0.05])


			
			csv_writer = csv.writer(open(self.input_dir_path + file.split("/")[-1],"w"),delimiter = "\t")
			csv_writer.writerows(seed_set)
			file_paths.append(self.input_dir_path + file.split("/")[-1])
		return file_paths


	def __run_DmGWAS__(self,
		PPI_network_file_path =  "/Users/leonardomartini/Documents/network_medicine/DiseaseGenePrioritizationAtRiskLoci/datasets/networks/STRING_PPI.tsv",):
		pass
	
	def __run_SigMod__(self,
	gene_score_file_path, 
	network_file_path = "/Users/leonardomartini/Documents/network_medicine/DiseaseGenePrioritizationAtRiskLoci/datasets/networks/STRING_PPI.tsv"):
		
		algorithm_path = self.algorithm_dir_path + "SIGMOID"
		if not os.path.exists(algorithm_path):
			os.makedirs(algorithm_path)
		output_file_path = algorithm_path +"/"+ gene_score_file_path.split("/")[-1]
		SIGMOD = f'cd sota/SigMod/; Rscript src.R {network_file_path} {gene_score_file_path} {output_file_path}'
		print(SIGMOD)
		subprocess.call(SIGMOD, shell=True,env = os.environ,stdout=subprocess.PIPE)
		
	
	def __run_domino__(self,
		seed_set_file_path,
		PPI_network_file_path =  "/Users/leonardomartini/Documents/network_medicine/DiseaseGenePrioritizationAtRiskLoci/datasets/networks/STRING_PPI.tsv",
		slicer_file_path = "/Users/leonardomartini/Documents/network_medicine/DiseaseGenePrioritizationAtRiskLoci/datasets/networks/STRING.sif"
		):

		algorithm_name = "DOMINO"
		output_file_name = seed_set_file_path.split("/")[-1]
		current_algorithm_dir = self.algorithm_dir_path + algorithm_name + "/"

		if not os.path.exists(current_algorithm_dir):
			os.makedirs(current_algorithm_dir)

		DOMINO = 'cd sota/DOMINO/src/; python3 main.py'
		slicer = 'cd sota/DOMINO/src/; python3 slicer_creator.py'
		
		if not os.path.exists(slicer_file_path):
			command = f'{slicer} -n {PPI_network_file_path} -o {slicer_file_path}'
			subprocess.call(command, shell=True,env = os.environ,stdout=subprocess.PIPE)
		
		if not os.path.exists(current_algorithm_dir + output_file_name):
			command = f'{DOMINO} -n {PPI_network_file_path} -a {seed_set_file_path} -s {slicer_file_path} -o {current_algorithm_dir} -f {output_file_name}'
			print(command)
		
			subprocess.call(command, shell=True,env = os.environ,stdout=subprocess.PIPE)
		
		solution = []
		
		with open(current_algorithm_dir + output_file_name,"r") as fp:
			lines = fp.readlines()
			for line in lines:
				
				current_line = line.replace("[","")
				current_line = current_line.replace("]","")
				current_line = current_line.replace(" ","")
				current_line = current_line.replace('"','')

				vector = current_line.split(",")
				for item in vector:
					solution.append([item])

		csv_writer = csv.writer(open(current_algorithm_dir + output_file_name,"w"),delimiter = "\t")
		csv_writer.writerows(solution)



	def __run_lean__(self,
	gene_score_file_path,
	network_file_path = "/Users/leonardomartini/Documents/network_medicine/DiseaseGenePrioritizationAtRiskLoci/datasets/networks/STRING_PPI.sif"):
		algorithm_path = self.algorithm_dir_path + "LEAN"
		if not os.path.exists(algorithm_path):
			os.makedirs(algorithm_path)
		output_file_path = algorithm_path +"/"+ gene_score_file_path.split("/")[-1]
		LEAN = f'cd sota/LEAN/; Rscript src.R {gene_score_file_path} {network_file_path} {output_file_path}'
		print(LEAN)
	
	def __run_depict__(self, configuration_file):
		
		algorithm_name = "DEPICT"
		current_algorithm_dir = self.algorithm_dir_path + algorithm_name + "/"
		
		if not os.path.exists(current_algorithm_dir):
			os.makedirs(current_algorithm_dir)
		
		DEPICT = 'cd sota/DEPICT/src/python/; python3 depict.py'
		command = f'{DEPICT} {configuration_file}'

		subprocess.call(command, shell=True,env = os.environ,stdout=subprocess.PIPE)


	def __run_DmGWAS_on_datasets__(self,):
		pass

	def run(self,):
		
		self.__load_node_PPI_network__()
		self.__load_ensembl_db__()

		seed_set_files = self.__compute_input_data_for_network_based_approach__()

		for file in seed_set_files:
			self.__run_lean__(file)
			exit(1)
			self.__run_domino__(file)
			self.__run_SigMod__(file)




aw = AlgorithmWrapper(
	disease_dir_path = "/Users/leonardomartini/Documents/network_medicine/Data/exps/RMM-GWAS/algorithm_comparison/",
	)

aw.run()