import pandas as pd
import csv
import os
import random

import seaborn as sns
sns.set_theme(style="whitegrid")
import matplotlib.pyplot as plt
import matplotlib as mpl
import json
  
class DrugHubValidation():
	
	def __init__(
		
		self,
		disease_experiment_dir_path,
		drug_table_file_path,
		random_solution_dir,
		config_file,
		ensembl_db_file_path,
		filter_
  		
		):
		self.config_file = config_file
		self.map__trait__disease_info = json.load(open(self.config_file))

		self.ensembl_db = ensembl_db_file_path

		self.random_solution_dir = random_solution_dir
		self.validation_dir_path = disease_experiment_dir_path + "validation/"
		
		self.disease_experiment_dir_path = disease_experiment_dir_path
		self.algorithms_file_path = disease_experiment_dir_path + "algorithms/"
		
		self.seed_dir_path = "../../experiments/GWAS/seed/"
		self.seed_dir_rmm_gwas_path = "../../experiments/GWAS/seed_RMM-GWAS/"

		if not os.path.exists(self.validation_dir_path):
			os.makedirs(self.validation_dir_path)

		self.filter_ = filter_
		self.__load_ensembl_db__()
		self.map__algorithm__solutions = {dir_: self.__load_solutions__(self.algorithms_file_path + dir_ +"/", algorithm_name = dir_) for dir_ in os.listdir(self.algorithms_file_path) if dir_[0] != "." and dir_ not in self.filter_}

		self.drug_table_file_path = drug_table_file_path

	def __load_random_solution__(self,):
		
		self.random_solution_file_paths = [ self.random_solution_dir+ file for file in os.listdir(self.random_solution_dir) if file[0] != "."]
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
			
			map__disease_name__disease_module[file_path.split("/")[-1].replace(".tsv","")] = set_

		return map__disease_name__disease_module

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

	def __load_disease_seed_for_RMM_GWAS__(self):

		file_paths = [file for file in os.listdir(self.seed_dir_rmm_gwas_path) if file[0] != "."]
		map__disease__seed = {}
		map__disease__gene__locus = {}
		for file in file_paths:
			map__disease__gene__locus[file.split("/")[-1].replace(".tsv","")] = {}
			csv_reader = csv.reader(open(self.seed_dir_rmm_gwas_path +file, "r"),delimiter = "\t")
			set_ = set()
			for row in csv_reader:
				set_.add(row[0])
				map__disease__gene__locus[file.split("/")[-1].replace(".tsv","")][row[0]] = row[1]

			map__disease__seed[file.split("/")[-1].replace(".tsv","")] = set_
		return map__disease__seed,map__disease__gene__locus

	def __load_drug_repurposing_hub__(self,):
	
		self.map__target__drugs = {}
		self.map__clinical_phase_disease_area_indication__targets = {}
		
		with open(self.drug_table_file_path, "r") as fp:
			csv_reader = csv.reader(fp,delimiter = "\t")

			for index,row in enumerate(csv_reader):

				if index == 0:
					continue

				targets = row[3].split("|")
				disease_area = row[4]
				clinical_phase = row[1]
				mou = row[2]
				drug_name = row[0]
				indication = row[5]

				if (clinical_phase,disease_area,indication) not in self.map__clinical_phase_disease_area_indication__targets:
					self.map__clinical_phase_disease_area_indication__targets[(clinical_phase,disease_area,indication)] = set()
				

				for target in targets:
					
					if target in self.map__target__drugs:
						self.map__target__drugs[target].append((drug_name, disease_area, mou,clinical_phase,indication))
					else:
						self.map__target__drugs[target] = [(drug_name, disease_area, mou,clinical_phase,indication)]

					self.map__clinical_phase_disease_area_indication__targets[(clinical_phase,disease_area,indication)].add(target)



	def __compute_drug_targets__(self, query_disease_area = None, query_indications = None): 	
		target_nodes = set()
		

		for node in self.map__target__drugs:
			if node == "":
				continue

			target_drugs = self.map__target__drugs[node]

			for drug_name, disease_area, mou,clinical_phase,indication in target_drugs:


				if query_disease_area is not None:
					
					if query_disease_area in disease_area:
						

						if query_indications is not None:
							query_indications_vector = query_indications.split("|")
							for query_indication in query_indications_vector:
								if query_indication in indication:
									target_nodes.add(node)
						else:
							target_nodes.add(node)

				else:
					target_nodes.add(node)

							
		return target_nodes



	def run(self, trial = 100, validation_dir = "drug_hub/"):
		drug_hub_validation_dir =self.validation_dir_path + validation_dir
		
		if not os.path.exists(drug_hub_validation_dir):
			os.makedirs(drug_hub_validation_dir)

		self.__load_drug_repurposing_hub__()
		
		map__disease__seed = self.__load_disease_seed__()
		map__disease__seed_RMM_GWAS, map__disease__gene__locus = self.__load_disease_seed_for_RMM_GWAS__()

		precision_table = []
		drug_indication = []
		algorithm_pval = []
		random_distribution = []
		
		for algorithm, solutions in self.map__algorithm__solutions.items():
			
			
			for considered_gwas in self.map__trait__disease_info.keys():
				print(algorithm,considered_gwas)
				if considered_gwas not in solutions:
					continue

				drug_targets = self.__compute_drug_targets__(query_disease_area = self.map__trait__disease_info[considered_gwas][0], query_indications = None)
				
				
				drug_targets = drug_targets.intersection(map__disease__seed_RMM_GWAS[considered_gwas.split('__')[0]])

				if len(drug_targets) == 0:
					print(considered_gwas)
					print(map__disease__seed[considered_gwas])
					continue

				solution = solutions[considered_gwas]
				phi =  len(solution.intersection(drug_targets))
				precision_table.append([algorithm,considered_gwas.split('__')[0],phi/len(drug_targets)])
				intersected_targets = solution.intersection(drug_targets)

				for target in intersected_targets:
					drugs = self.map__target__drugs[target]
					for drug in drugs:
						
						drug_name, disease_area, mou,clinical_phase,indication = drug
						record = [target,drug_name,mou,clinical_phase, algorithm, considered_gwas]
						drug_indication.append(record)
			
				counter = 0
				
				for i in range(trial):

					
					if algorithm == 'RMM-GWAS':
						random_solution = set(random.sample(map__disease__seed_RMM_GWAS[considered_gwas],len(solution)))
					else:

						random_solution = set(random.sample(map__disease__seed[considered_gwas],len(solution)))
					
					phi_random = len(random_solution.intersection(drug_targets))
					
					random_distribution.append([algorithm, considered_gwas.split('__')[0],phi_random/len(drug_targets)])
					
					if phi <= phi_random:
						counter += 1

				counter_str = ""
				if counter == 0:
					counter_str = "p_{val} < 10^{-3}"
				else:
					counter_str = "p_{val} ~ " + str(counter/trial) 



				algorithm_pval.append([algorithm, considered_gwas.split('__')[0],counter_str])

		


		pd.DataFrame(precision_table, columns = ["Algorithm","GWAS", "Recall"]).to_csv(drug_hub_validation_dir + "metrics.tsv", sep = "\t")
		pd.DataFrame(random_distribution,columns = ["Algorithm","GWAS","Recall"]).to_csv(drug_hub_validation_dir + "metrics_random_distribution.tsv", sep = "\t")
		pd.DataFrame(algorithm_pval,columns = ["Algorithm","GWAS","p_val"]).to_csv(drug_hub_validation_dir + "metrics_p_val.tsv", sep = "\t")
		pd.DataFrame(drug_indication,columns = ["Gene","Drug Name", "Mechanism of Action", "Drug development phase", "Algorithm", "GWAS"]).to_csv(drug_hub_validation_dir  + "discovered.tsv", sep = "\t")



	def run_on_randomize_bipartite_graph(self,):
		self.__load_drug_repurposing_hub__()
		
		drug_targets = self.__compute_drug_targets__()
		self.__load_gene_pool__()
		solutions = self.__load_random_solution__()


		table_1 = []
		table_2 = []
		
		random_phis = []
		counter = 0
		phi = len(self.map__algorithm__solution["RMA"].intersection(drug_targets))/len(self.map__algorithm__solution["RMA"])
		
		for solution in solutions:
			phi_random = len(solution.intersection(drug_targets))/len(solution)

			if phi < phi_random:
				counter += 1

			random_phis.append(phi_random)
		
		data_frame = pd.DataFrame(random_phis, columns = ["Random Distribution"])
		p_val = round(float(counter/len(solutions)), 4)

		data_frame.to_csv(self.validation_dir_path + "drug_target_on_randomize_bipartite_graph.tsv", sep = "\t")
		pd.DataFrame([["RMA",p_val]],columns = ["Algorithm","p_val"]).to_csv(self.validation_dir_path + "drug_target_p_value_on_randomize_bipartite_graph.tsv", sep = "\t")
		

		sns_plot = sns.histplot(data = data_frame, x = "Random Distribution")
		sns_plot.annotate("$P_{val}$: " + str(p_val), xy=(phi, 4), xytext=(phi, 10),
            arrowprops=dict(arrowstyle="->",color='blue'))
		sns_plot.set_title("Drug Target Prediction")

		plt.show()

if __name__ == '__main__':
	
	
	d = DrugHubValidation(
		disease_experiment_dir_path = "../../experiments/algorithm_comparison/",
		drug_table_file_path = "../../datasets/curated_db/drug_repurposing_hub.tsv",
		random_solution_dir = "../../experiments/Robustness_Experiment/",
		config_file = "config_files/drug_hub.json",
		ensembl_db_file_path = "../../datasets/curated_db/mart_export.txt",
		filter_ = []
		)
	
	d.run( trial = 1000)
	





