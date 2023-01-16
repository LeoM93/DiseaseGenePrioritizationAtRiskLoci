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
		gene_pool_file_path,
		mouse_phenotype_db,
		random_solution_dir,
		filter_
		):
		
		self.random_solution_dir = random_solution_dir
		self.disease_experiment_dir_path = disease_experiment_dir_path
		self.algorithms_file_path = disease_experiment_dir_path + "algorithms/"
		self.validation_dir_path = disease_experiment_dir_path + "validation/"
		self.gene_pool_file_path = gene_pool_file_path

		if not os.path.exists(self.validation_dir_path):
			os.makedirs(self.validation_dir_path)

		self.filter_ = filter_
		self.map__algorithm__solution = {dir_: self.__load_solution__(self.algorithms_file_path + dir_ + "/solution.txt") for dir_ in os.listdir(self.algorithms_file_path) if dir_[0] != "." and dir_ not in self.filter_}


		self.mouse_phenotype_db = mouse_phenotype_db

	def __load_solution__(self,file_path):

		csv_reader = csv.reader(open(file_path , "r"),delimiter = "\t")
		set_ = set()
		
		for row in csv_reader:
			set_.add(row[0])

		return set_


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



	def __load_gene_pool__(self):

		csv_reader = csv.reader(open(self.gene_pool_file_path , "r"),delimiter = "\t")
		self.candidate_set = set()
		
		for row in csv_reader:
			self.candidate_set.add(row[0])

		return self.candidate_set

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

	def run(self,trial = 1000):
		
		self.__load_mouse_phenotype_db__()
		target_nodes = self.__compute_number_of_gene_in_phenotype_class__()
		self.__load_gene_pool__()

		table_1 = []
		table_2 = []
		table_3 = []


		for algorithm, solution in self.map__algorithm__solution.items():
			
			counter = 0

			print(algorithm)
			print(solution.intersection(target_nodes))

			phi =  len(solution.intersection(target_nodes))
			table_1.append([algorithm,phi/len(solution)])

			for i in range(trial):
				
				random_solution = set(random.sample(self.candidate_set,len(solution)))
				phi_random = len(random_solution.intersection(target_nodes))
				table_2.append([algorithm,phi_random/len(random_solution)])

				if phi < phi_random:
					counter += 1

			if counter == 0:
				counter_str = "p_{val} < 10^{-6}"
			else:
				counter_str = "p_{val} ~ " + str(counter/trial) 


			table_3.append([algorithm, counter_str])


		pd.DataFrame(table_1,columns = ["Algorithm","Score"]).to_csv(self.validation_dir_path + "open_target_intersection.tsv", sep = "\t")
		pd.DataFrame(table_2,columns = ["Algorithm","Score"]).to_csv(self.validation_dir_path + "open_target_intersection_random.tsv", sep = "\t")
		pd.DataFrame(table_3,columns = ["Algorithm","p_val"]).to_csv(self.validation_dir_path + "open_target_intersection_p_val.tsv", sep = "\t")


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
	disease_experiment_dir_path = "../../experiments/algorithm_comparison/", 
	gene_pool_file_path = "../../datasets/GWAS/unbiased_copd_snp_database.txt", 
	random_solution_dir = "../../experiments/Robustness_Experiment/",
	mouse_phenotype_db = "../../datasets/curated_db/mouse_phenotype_COPD.tsv",
	filter_ = []
)
op.run()
op.run_on_randomize_bipartite_graph()


