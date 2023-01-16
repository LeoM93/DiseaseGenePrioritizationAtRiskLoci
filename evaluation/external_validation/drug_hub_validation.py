import pandas as pd
import csv
import os
import random

import seaborn as sns
sns.set_theme(style="whitegrid")
import matplotlib.pyplot as plt
import matplotlib as mpl

class DrugHubValidation():
	
	def __init__(
		
		self,
		disease_experiment_dir_path,
		gene_pool_file_path,
		drug_table_file_path,
		random_solution_dir,
		filter_
  		
		):
		self.random_solution_dir = random_solution_dir
		self.validation_dir_path = disease_experiment_dir_path + "validation/"
		
		self.disease_experiment_dir_path = disease_experiment_dir_path
		self.algorithms_file_path = disease_experiment_dir_path + "algorithms/"
		self.gene_pool_file_path = gene_pool_file_path

		if not os.path.exists(self.validation_dir_path):
			os.makedirs(self.validation_dir_path)

		self.filter_ = filter_
		self.map__algorithm__solution = {dir_: self.__load_solution__(self.algorithms_file_path + dir_ + "/solution.txt") for dir_ in os.listdir(self.algorithms_file_path) if dir_[0] != "." and dir_ not in self.filter_}


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


	def __load_solution__(self,file_path):

		csv_reader = csv.reader(open(file_path , "r"),delimiter = "\t")
		set_ = set()
		
		for row in csv_reader:
			set_.add(row[0])

		return set_

	def __load_gene_pool__(self):

		csv_reader = csv.reader(open(self.gene_pool_file_path , "r"),delimiter = "\t")
		self.candidate_set = set()
		
		for row in csv_reader:
			self.candidate_set.add(row[0])

		return self.candidate_set


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

			target_drugs = self.map__target__drugs[node]

			for drug_name, disease_area, mou,clinical_phase,indication in target_drugs:


				if query_disease_area is not None and query_indications is not None:
					
					if query_disease_area in disease_area:
						query_indications_vector = query_indications.split("|")

						for query_indication in query_indications_vector:
							if query_indication in indication:
								target_nodes.add(node)
				else:
					target_nodes.add(node)

							
		return target_nodes



	def run(self, trial = 1000000, pulmonary = False):

		self.__load_drug_repurposing_hub__()
		
		if pulmonary:
			drug_targets = self.__compute_drug_targets__(query_disease_area = "pulmonary", query_indications = "chronic obstructive pulmonary disease (COPD)|asthma|bronchospasm")
		else:
			drug_targets = self.__compute_drug_targets__()

		self.__load_gene_pool__()


		table_1 = []
		table_2 = []
		table_3 = []
		for algorithm, solution in self.map__algorithm__solution.items():
			phi =  len(solution.intersection(drug_targets))
			table_1.append([algorithm,phi/len(solution)])
			print(algorithm,solution.intersection(drug_targets))

			counter = 0
			
			
			for i in range(trial):
				
				random_solution = set(random.sample(self.candidate_set,len(solution)))
				phi_random = len(random_solution.intersection(drug_targets))
				table_2.append([algorithm,phi_random/len(random_solution)])

				if phi < phi_random:
					counter += 1

			counter_str = ""
			if counter == 0:
				counter_str = "p_{val} < 10^{-6}"
			else:
				counter_str = "p_{val} ~ " + str(counter/trial) 


			table_3.append([algorithm, counter_str])



		if pulmonary:
			pd.DataFrame(table_1,columns = ["Algorithm","Score"]).to_csv(self.validation_dir_path + "drug_hub_drug_target_intersection_pulmonary.tsv", sep = "\t")
			pd.DataFrame(table_2,columns = ["Algorithm","Score"]).to_csv(self.validation_dir_path + "drug_hub_drug_target_intersection_pulmonary_random.tsv", sep = "\t")
			pd.DataFrame(table_3,columns = ["Algorithm","p_val"]).to_csv(self.validation_dir_path + "drug_hub_drug_target_intersection_pulmonary_p_val.tsv", sep = "\t")
		else:
			pd.DataFrame(table_1,columns = ["Algorithm","Score"]).to_csv(self.validation_dir_path + "drug_hub_drug_target_intersection.tsv", sep = "\t")
			pd.DataFrame(table_2,columns = ["Algorithm","Score"]).to_csv(self.validation_dir_path + "drug_hub_drug_target_intersection_random.tsv", sep = "\t")
			pd.DataFrame(table_3,columns = ["Algorithm","p_val"]).to_csv(self.validation_dir_path + "drug_hub_drug_target_intersection_p_val.tsv", sep = "\t")
		


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
		gene_pool_file_path = "../../datasets/GWAS/unbiased_copd_snp_database.txt",
		drug_table_file_path = "../../datasets/curated_db/drug_repurposing_hub.tsv",
		random_solution_dir = "../../experiments/Robustness_Experiment/",
		filter_ = []
		)
	
	d.run( trial = 1000, pulmonary = False)
	d.run( trial = 1000, pulmonary = True)
	d.run_on_randomize_bipartite_graph()	
	





