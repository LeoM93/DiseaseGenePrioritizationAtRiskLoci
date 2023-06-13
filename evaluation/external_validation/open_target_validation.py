import csv
import os
import numpy as np
import pandas as pd
import random

import seaborn as sns
sns.set_theme(style="whitegrid")
import matplotlib.pyplot as plt
import matplotlib as mpl
class OpenTarget():

	def __init__(self, 
		
		disease_experiment_dir_path,
		mouse_phenotype_db,
		random_solution_dir,
		ensembl_db_file_path,
		filter_
		):
		
		self.random_solution_dir = random_solution_dir
		self.disease_experiment_dir_path = disease_experiment_dir_path
		self.algorithms_file_path = disease_experiment_dir_path + "algorithms/"
		self.validation_dir_path = disease_experiment_dir_path + "validation/"
		self.seed_dir_path =  "../../experiments/input/seed/"
		if not os.path.exists(self.validation_dir_path):
			os.makedirs(self.validation_dir_path)

		self.filter_ = filter_
		self.ensembl_db = ensembl_db_file_path
		self.__load_ensembl_db__()
		self.map__algorithm__solutions = {dir_: self.__load_solutions__(self.algorithms_file_path + dir_ +"/",dir_) for dir_ in os.listdir(self.algorithms_file_path) if dir_[0] != "." and dir_ not in self.filter_}
		

		self.mouse_phenotype_db = mouse_phenotype_db

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
	
	def __load_solutions__(self,directory_path, algorithm_name):
		
		file_paths = [directory_path + file for file in os.listdir(directory_path) if file[0] != "." ]
		map__disease_name__disease_module = {}
		
		for file_path in file_paths:
			csv_reader = csv.reader(open(file_path , "r"),delimiter = "\t")
			set_ = set()
		
			for row in csv_reader:
				if algorithm_name != "RMM-GWAS":
					if row[0] in self.map__ensembl_id__gene:
						set_.add(self.map__ensembl_id__gene[row[0]])
				else:
					set_.add(row[0])
			
			map__disease_name__disease_module[file_path.split("/")[-1].split("__")[0].replace(".tsv","")] = set_

		return map__disease_name__disease_module


	def __load_random_solution__(self,):
		self.random_solution_file_paths = [self.random_solution_dir + file for file in os.listdir(self.random_solution_dir) if file[0] != "."]
		solutions = []
		for file in self.random_solution_file_paths:
			csv_reader = csv.reader(open(file, "r"),delimiter = "\t")

			for index, row in enumerate(csv_reader):
				if index == 0:
					continue

				solution = row[-1]
				solution = solution.replace("{","")
				solution = solution.replace("}","")
				solution = solution.replace("'","")
				solution = solution.replace(" ","")
				solution = set(solution.split(","))


				solutions.append(solution)
		return solutions



	def __load_disease_seed__(self):

		file_paths = [file for file in os.listdir(self.seed_dir_path) if file[0] != "."]
		map__disease__seed = {}
		for file in file_paths:
			csv_reader = csv.reader(open(self.seed_dir_path +file, "r"),delimiter = "\t")
			set_ = set()
			for row in csv_reader:
				if row[0] in self.map__ensembl_id__gene:
					set_.add(self.map__ensembl_id__gene[row[0]])

			map__disease__seed[file.split("/")[-1].replace(".tsv","")] = set_
		return map__disease__seed

	def __load_mouse_phenotype_db__(self,):
		
		csv_reader = csv.reader(open(self.mouse_phenotype_db , "r"),delimiter = "\t")
		
		self.map__gene_name__phenotype = {}
		for index, row in enumerate(csv_reader):

			if index == 0:
				continue

			gene_name  = row[0]
			label = row[1]
			class_ = row[2]
			phenotype_id = row[2]

			if gene_name not in self.map__gene_name__phenotype:
				self.map__gene_name__phenotype[gene_name] = set()

			self.map__gene_name__phenotype[gene_name].add((label,class_,phenotype_id))

	
	def __compute_number_of_gene_in_phenotype_class__(self, phenotype_class = "respiratory system phenotype"):

		target_nodes = set()
		for gene in self.map__gene_name__phenotype:
				
			for phenotype in self.map__gene_name__phenotype[gene]:
				label,class_,phenotype_id = phenotype


				if class_ == phenotype_class:
					target_nodes.add(gene)
					break


		return target_nodes

	def run(self,trial = 1000, validation_dir = "mouse_phenotypes/"):
		
		open_target_validation_dir =self.validation_dir_path + validation_dir
		
		if not os.path.exists(open_target_validation_dir):
			os.makedirs(open_target_validation_dir)
		
		self.__load_mouse_phenotype_db__()
		
		map__disease__seed = self.__load_disease_seed__()

		precision_table = []
		drug_indication = []
		algorithm_pval = []
		random_distribution = []
		print(map__disease__seed)
		
		for algorithm, solutions in self.map__algorithm__solutions.items():
			for gwas, solution in solutions.items():
				if gwas in ["GCST009841", "GCST007692","GCST004748"]:
					target_nodes = self.__compute_number_of_gene_in_phenotype_class__()
					target_nodes = target_nodes.intersection(map__disease__seed[gwas])
					phi =  len(solution.intersection(target_nodes))
					precision_table.append([algorithm,gwas,phi/len(solution)])

					counter = 0

					for i in range(trial):
				
						random_solution = set(random.sample(map__disease__seed[gwas],len(solution)))
						phi_random = len(random_solution.intersection(target_nodes))
						random_distribution.append([algorithm, gwas,phi_random/len(random_solution)])

						if phi <= phi_random:
							counter += 1
					if counter == 0:
						counter_str = "p_{val} < 10^{-6}"
					else:
						counter_str = "p_{val} ~ " + str(counter/trial) 


					algorithm_pval.append([algorithm, gwas,counter_str])
		
		pd.DataFrame(precision_table, columns = ["Algorithm","GWAS", "Score"]).to_csv(open_target_validation_dir + "metrics.tsv", sep = "\t")
		pd.DataFrame(random_distribution,columns = ["Algorithm","GWAS","Score"]).to_csv(open_target_validation_dir + "metrics_random_distribution.tsv", sep = "\t")
		pd.DataFrame(algorithm_pval,columns = ["Algorithm","GWAS","p_val"]).to_csv(open_target_validation_dir + "metrics_p_val.tsv", sep = "\t")


	def run_on_randomize_bipartite_graph(self,):
		self.__load_mouse_phenotype_db__()
		
		solutions = self.__load_random_solution__()
		target_nodes = self.__compute_number_of_gene_in_phenotype_class__()

		random_phis = []
		counter = 0
		phi = len(self.map__algorithm__solution["RMA"].intersection(target_nodes))/len(self.map__algorithm__solution["RMA"])
		for solution in solutions:
			phi_random = len(solution.intersection(target_nodes))/len(solution)
			random_phis.append(phi_random)

			if phi < phi_random:
				counter += 1


		data_frame = pd.DataFrame(random_phis, columns = ["Random Distribution"])
		data_frame.to_csv(self.validation_dir_path + "open_target_on_randomize_bipartite_graph.tsv", sep = "\t")

		p_val = round(float(counter/len(solutions)), 4)
		pd.DataFrame([["RMA",p_val]],columns = ["Algorithm","p_val"]).to_csv(self.validation_dir_path + "open_target_p_value_on_randomize_bipartite_graph.tsv", sep = "\t")

		sns_plot = sns.histplot(data = data_frame, x = "Random Distribution")
		sns_plot.annotate("$P_{val}$: " + str(p_val), xy=(phi, 0.5), xytext=(phi, 10),
            arrowprops=dict(arrowstyle="->",color='blue'))
		sns_plot.set_title("Mouse Phenotypes Prediction")

		plt.show()

			


op = OpenTarget(
	disease_experiment_dir_path = "../../experiments/algorithm_comparison_GWAS_2Mb/", 
	random_solution_dir = "../../experiments/Robustness_Experiment/",
	mouse_phenotype_db = "../../datasets/curated_db/mouse_phenotype_COPD.tsv",
	ensembl_db_file_path = "../../datasets/curated_db/mart_export.txt",

	filter_ = []
)
op.run()


