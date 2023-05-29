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
		self.disease_genes_dir_path = disease_dir_path + "seed_set/"
		self.disease_dir_path = disease_dir_path

		self.ensembl_db = "datasets/curated_db/mart_export.txt"
		
		self.algorithm_dir_path = self.disease_dir_path + "algorithms/"

		if not os.path.exists(self.algorithm_dir_path):
			os.makedirs(self.algorithm_dir_path)
		
		self.disease_seed_set_file_paths = [self.disease_genes_dir_path + file for file in os.listdir(self.disease_genes_dir_path) if file[0] != "."]



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
			self. map__ensembl_id__gene[ensembl_id] = gene_name
	

	def __compute_input_data_for_domino__(self,gene_pool_file_path):
		csv_reader = csv.reader(open(gene_pool_file_path,"r"),delimiter = "\t")
		
		domino_seed_file_path = gene_pool_file_path.replace(".tsv","__domino_seed.tsv")
		domino_input_data = []
		
		for index,  row in enumerate(csv_reader):
			if row[0] in self.map__gene__ensembl_id:	
				domino_input_data.append([self.map__gene__ensembl_id[row[0]],1])
		
		csv_writer = csv.writer(open(domino_seed_file_path,"w"),delimiter = "\t")
		csv_writer.writerows(domino_input_data)
		return domino_seed_file_path


	def __run_DmGWAS__(self,
		PPI_network_file_path =  "/Users/leonardomartini/Documents/network_medicine/DiseaseGenePrioritizationAtRiskLoci/datasets/networks/STRING_PPI.tsv",):
		pass
	
	def __run_SigMod__(self,gene_score_file_path, network_file_path):
		
		SIGMOD = f'cd sota/SigMod/; Rscript src.R {network_file_path} {gene_score_file_path}'
		subprocess.call(SIGMOD, shell=True,env = os.environ,stdout=subprocess.PIPE)
		
	
	def __run_domino__(self,
		seed_set_file_path,
		PPI_network_file_path =  "/Users/leonardomartini/Documents/network_medicine/DiseaseGenePrioritizationAtRiskLoci/datasets/networks/STRING_PPI.tsv",
		slicer_file_path = "/Users/leonardomartini/Documents/network_medicine/DiseaseGenePrioritizationAtRiskLoci/datasets/networks/STRING.sif"
		):

		algorithm_name = "DOMINO"
		output_file_name = seed_set_file_path.split("__")[1]
		current_algorithm_dir = self.algorithm_dir_path + algorithm_name + "/"

		if not os.path.exists(current_algorithm_dir):
			os.makedirs(current_algorithm_dir)

		DOMINO = 'cd sota/DOMINO/src/; python3 main.py'
		slicer = 'cd sota/DOMINO/src/; python3 slicer_creator.py'
		
		if not os.path.exists(slicer_file_path):
			command = f'{slicer} -n {PPI_network_file_path} -o {slicer_file_path}'
			subprocess.call(command, shell=True,env = os.environ,stdout=subprocess.PIPE)
		
		command = f'{DOMINO} -n {PPI_network_file_path} -a {seed_set_file_path} -s {slicer_file_path} -o {current_algorithm_dir} -f {output_file_name}'
		print(command)
		
		subprocess.call(command, shell=True,env = os.environ,stdout=subprocess.PIPE)

	def __run_depict__(self, configuration_file):
		
		algorithm_name = "DEPICT"
		current_algorithm_dir = self.algorithm_dir_path + algorithm_name + "/"
		
		if not os.path.exists(current_algorithm_dir):
			os.makedirs(current_algorithm_dir)
		
		DEPICT = 'cd sota/DEPICT/src/python/; python3 depict.py'
		command = f'{DEPICT} {configuration_file}'

		subprocess.call(command, shell=True,env = os.environ,stdout=subprocess.PIPE)
		
	def run(self,):
		
		for seed_set_file_path in self.disease_seed_set_file_paths:
			
			if len(seed_set_file_path.split("__")) > 2:
				continue

			gene_pool_file_path = seed_set_file_path.replace(".tsv","__gene_pool.tsv")
			
			bg = SNPGeneGraph(
				rs_id_position_file_path = seed_set_file_path,
				gene_overlap_2mb_windows_file_path = "datasets/GWAS/refSeqGenes_hg19.txt",
				output_file_path = gene_pool_file_path
				)
			
			bg.run()

			domino_seed_file_path = self.__compute_input_data_for_domino__(gene_pool_file_path)
			
			#Wrap for domino
			self.__run_domino__(domino_seed_file_path)
			
		configuration_file = [
			'/Users/leonardomartini/Documents/network_medicine/DiseaseGenePrioritizationAtRiskLoci/sota/DEPICT/example/configuration.cfg'
			]
			
		for configuration_file in configuration_files:
			#Wrap for depict:
			self.__run_depict__(configuration_file)
		
		gene_score_file_paths = []
		for gene_score_file_path in gene_score_file_paths:
			self.__run_SigMod__(gene_score_file_path = gene_score_file_path)


aw = AlgorithmWrapper(
	disease_dir_path = "/Users/leonardomartini/Documents/network_medicine/Data/exps/RMM-GWAS/algorithm_comparison/",
	)
aw.run()

