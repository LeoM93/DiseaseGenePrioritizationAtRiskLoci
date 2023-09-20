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
		GWAS_dir_path,

	):
		self.GWAS_dir = GWAS_dir_path + "studies/"
		self.GWAS_config_file = "/Users/leonardomartini/Documents/network_medicine/DiseaseGenePrioritizationAtRiskLoci/sota/DEPICT/GWAS/"
		
		self.disease_dir_path = disease_dir_path

		if not os.path.exists(self.disease_dir_path):
			os.makedirs(self.disease_dir_path)
		
		self.input_dir_path = GWAS_dir_path + "seed/"
		self.vagas_2_dir = GWAS_dir_path + "vegas_2/output/"
		
		self.input_dir_path_for_RMM_GWAS = GWAS_dir_path + "seed_RMM-GWAS/"
		self.input_network_dir_path_for_RMM_GWAS = GWAS_dir_path + "newtork_RMM-GWAS/"
		self.rmm_gwas_bipartite_input_graph = []
		
		self.ensembl_db = "datasets/curated_db/mart_export.txt"
		
		self.algorithm_dir_path = self.disease_dir_path + "algorithms/"

		if not os.path.exists(self.input_dir_path_for_RMM_GWAS):
			os.makedirs(self.input_dir_path_for_RMM_GWAS)
		
		if not os.path.exists(self.input_network_dir_path_for_RMM_GWAS):
			os.makedirs(self.input_network_dir_path_for_RMM_GWAS)
		
		if not os.path.exists(self.input_dir_path):
			os.makedirs(self.input_dir_path)
		
		if not os.path.exists(self.algorithm_dir_path):
			os.makedirs(self.algorithm_dir_path)
		
		self.GWAS_file_paths = [self.GWAS_dir + file for file in os.listdir(self.GWAS_dir) if file[0] != "."]
		self.GWAS_config_files = [self.GWAS_config_file + file for file in os.listdir(self.GWAS_config_file) if file[0] != "."]
	
	def __load_vegas_2__(self,file_name):
		map__gene__vegas_pval = {}
		try:
			csv_reader = csv.reader(open(self.vagas_2_dir + file_name, "r"), delimiter = " ")
		except:
			return {}

		for index, row in enumerate(csv_reader):
			if index == 0:
				continue
			row[1]
			float(row[-1])

			map__gene__vegas_pval[row[1]] = float(row[-1])
		return map__gene__vegas_pval
	
	def __load_gene_in_2MB_windows__(self,file_path):
		
		bg = SNPGeneGraph(
				rs_id_position_file_path = file_path,
				gene_overlap_2mb_windows_file_path = "datasets/GWAS/refSeqGenes_hg19.txt",
				output_file_path = self.input_dir_path_for_RMM_GWAS + file_path.split("/")[-1]
				)
		self.rmm_gwas_bipartite_input_graph.append(self.input_dir_path_for_RMM_GWAS + file_path.split("/")[-1])
		map__gene__pval = bg.run()
		table = []
		for k,v in map__gene__pval.items():
			table.append([k,v])
		
		return table


		
	
	
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
		window_file_paths = []
		vegas_2_file_paths =[] 
		for file in self.GWAS_file_paths:
			
			filtered_map__gene__pval = self.__load_gene_in_2MB_windows__(file)
			csv_writer = csv.writer(open(self.input_dir_path + file.split("/")[-1].replace(".tsv","") + "__closest_genes.tsv","w"),delimiter = "\t")
			csv_writer.writerows(filtered_map__gene__pval)
			
			file_paths.append(self.input_dir_path + file.split("/")[-1].replace(".tsv","") + "__closest_genes.tsv")

		return file_paths


	def __run_DmGWAS__(self, file,
		PPI_network_file_path =  "/Users/leonardomartini/Documents/network_medicine/DiseaseGenePrioritizationAtRiskLoci/datasets/networks/STRING_PPI.tsv",):
		algorithm_path = self.algorithm_dir_path + "DmGWAS"
		if not os.path.exists(algorithm_path):
			os.makedirs(algorithm_path)
		output_file_path = algorithm_path +"/"+ file.split("/")[-1]
		dmGWAS = f'cd sota/DmGWAS/; Rscript src.R {PPI_network_file_path} {file} {output_file_path}'
		print(dmGWAS)
		subprocess.call(dmGWAS, shell=True,env = os.environ,stdout=subprocess.PIPE)
	
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
					solution.append([item.strip()])

		csv_writer = csv.writer(open(current_algorithm_dir +output_file_name,"w"),delimiter = "\t")
		csv_writer.writerows(solution)



	def __run_lean__(self,
	gene_score_file_path,
	
	network_file_path = "/Users/leonardomartini/Documents/network_medicine/DiseaseGenePrioritizationAtRiskLoci/datasets/networks/STRING_PPI.sif"):
		algorithm_path = self.algorithm_dir_path + "LEAN"
		
		if not os.path.exists(algorithm_path):
			os.makedirs(algorithm_path)
		output_file_path = algorithm_path +"/"+ gene_score_file_path.split("/")[-1]
		LEAN = f'cd sota/LEAN/; Rscript src.R {gene_score_file_path} {network_file_path} {output_file_path}'
		
		subprocess.call(LEAN, shell=True,env = os.environ,stdout=subprocess.PIPE)
		csv_reader = csv.reader(open(output_file_path),delimiter = ",")
		solution = []
		for index, row in enumerate(csv_reader):
			if index == 0:
				continue
			solution.append([row[0]])
		
		csv_writer = csv.writer(open(output_file_path,"w"),delimiter = "\t")
		csv_writer.writerows(solution)
	
	
	def __run_depict__(self, configuration_file):
		
		algorithm_name = "DEPICT"
		current_algorithm_dir = self.algorithm_dir_path + algorithm_name + "/"
		
		if not os.path.exists(current_algorithm_dir):
			os.makedirs(current_algorithm_dir)
		
		DEPICT = 'cd sota/DEPICT/src/python/; python3 depict.py'
		command = f'{DEPICT} {configuration_file}'
		print(command)
		subprocess.call(command, shell=True,env = os.environ,stdout=subprocess.PIPE)


	def __run__RMM_GWAS__(self,
	file_path, 
	network_file_path = "/Users/leonardomartini/Documents/network_medicine/DiseaseGenePrioritizationAtRiskLoci/datasets/networks/co_regulated_network_1e-5_thresholded.txt"):
		algorithm_name = "RMM-GWAS"
		current_algorithm_dir = self.algorithm_dir_path + algorithm_name + "/"
		disease_name = file_path.split("/")[-1]
		if not os.path.exists(current_algorithm_dir):
			os.makedirs(current_algorithm_dir)
		
		preprocessing = 'cd data_preprocessing; python3 main.py'
		command = f'{preprocessing} {network_file_path} {self.input_network_dir_path_for_RMM_GWAS + disease_name} {file_path}'
		
		print(command)
		subprocess.call(command, shell=True,env = os.environ,stdout=subprocess.PIPE)

		run_rmm_gwas =  f'cd rmm_gwas; python3 run_relations_maximization.py {self.input_network_dir_path_for_RMM_GWAS + disease_name} {current_algorithm_dir + disease_name}'
		print(run_rmm_gwas)
		subprocess.call(run_rmm_gwas, shell=True,env = os.environ,stdout=subprocess.PIPE)

	def run(self, only_our_framework = False):
				
		self.__load_node_PPI_network__()
		self.__load_ensembl_db__()

		file_paths = self.__compute_input_data_for_network_based_approach__()
		
		
		for file_path in self.rmm_gwas_bipartite_input_graph:
			self.__run__RMM_GWAS__(file_path)
		
		if only_our_framework:
			return
		
		for file in file_paths:
			self.__run_SigMod__(file)
			self.__run_domino__(file)
	
				



aw = AlgorithmWrapper(
	disease_dir_path = "/Users/leonardomartini/Documents/network_medicine/DiseaseGenePrioritizationAtRiskLoci/experiments/algorithm_comparison/",
	GWAS_dir_path = "/Users/leonardomartini/Documents/network_medicine/DiseaseGenePrioritizationAtRiskLoci/experiments/GWAS/"
	)

aw.run()
