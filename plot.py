import os
import math
import csv
import pandas as pd
import seaborn as sns
sns.set_theme(style="whitegrid")
import string

import matplotlib.pyplot as plt
import matplotlib as mpl

class Plot():
	def __init__(self,
		disease_experiment_dir_path,
 	):
		
		self.reactome_enrichment_dir = disease_experiment_dir_path + "validation/gsea_reactome/"
		try:
			self.reactome_enrichment_file_paths = [self.reactome_enrichment_dir + file for file in os.listdir(self.reactome_enrichment_dir) if file[0] != "." and ".tsv" in file]
			self.reactome_enrichment_alg_file_paths = {}
			
			for alg in os.listdir(self.reactome_enrichment_dir):
				if alg[0] != "." and ".tsv" not in alg:

					if alg not in self.reactome_enrichment_alg_file_paths:

						self.reactome_enrichment_alg_file_paths[alg] = [self.reactome_enrichment_dir +alg+"/" + file for file in os.listdir(self.reactome_enrichment_dir +alg+"/") if file[0] != "."]

		except:
			pass

		self.validation_dir_path = disease_experiment_dir_path + "validation/"

		self.mouse_phenotype_validation = self.validation_dir_path + "open_target_intersection.tsv"
		self.mouse_phenotype_random = self.validation_dir_path + "open_target_intersection_random.tsv"
		self.mouse_phenotype_p_val = self.validation_dir_path + "open_target_intersection_p_val.tsv"

		self.drug_hub_validation = self.validation_dir_path + "drug_hub_drug_target_intersection.tsv"
		self.drug_hub_random = self.validation_dir_path + "drug_hub_drug_target_intersection_random.tsv"
		self.drug_hub_p_val = self.validation_dir_path + "drug_hub_drug_target_intersection_p_val.tsv"

		self.drug_hub_pulmonary_validation = self.validation_dir_path + "drug_hub_drug_target_intersection_pulmonary.tsv"
		self.drug_hub_pulmonary_random = self.validation_dir_path + "drug_hub_drug_target_intersection_pulmonary_random.tsv"
		self.drug_hub_pulmonary_p_val = self.validation_dir_path + "drug_hub_drug_target_intersection_pulmonary_p_val.tsv"


		self.mouse_phenotype_random_bg_distribution = self.validation_dir_path + "open_target_on_randomize_bipartite_graph.tsv"
		self.mouse_phenotype_random_bg_p_value = self.validation_dir_path + "open_target_p_value_on_randomize_bipartite_graph.tsv"

		self.drug_hub_random_bg_distribution = self.validation_dir_path + "drug_target_on_randomize_bipartite_graph.tsv"
		self.drug_hub_random_bg_p_value = self.validation_dir_path + "drug_target_p_value_on_randomize_bipartite_graph.tsv"

		self.odd_vs_even_bma_validation = self.validation_dir_path + "bma_values.tsv"
		self.odd_vs_even_bma_random = self.validation_dir_path + "bma_distribution.tsv"
		self.odd_vs_even_bma_p_values = self.validation_dir_path + "bma_p_values.tsv"
		
		self.odd_vs_even_np_validation = self.validation_dir_path + "network_proximity_values.tsv"
		self.odd_vs_even_np_random = self.validation_dir_path + "network_proximity_distribution.tsv"
		self.odd_vs_even_np_p_values = self.validation_dir_path + "network_proximity_p_values.tsv"
		
		self.np_missing_data_validation = self.validation_dir_path + "network_distance_distribution_missing_data.tsv"
		self.bma_missing_data_validation = self.validation_dir_path + "bma_distribution_missing_data.tsv"

	

	def __load_reactome_enrichment_data_frame__(self, file_paths, str_ = " Homo sapiens"):
		self.map__pathways__algorithm__pval = {}
		
		for file_path in file_paths:
			algorithm_file_path = file_path.split("/")[-1]
			
			csv_reader = csv.reader(open(file_path, "r"), delimiter = "\t")

			for index, row in enumerate(csv_reader):
				if index == 0:
					continue

				term = row[1]
				p_value  = float(row[2])


				if str_ in term:
					pathways_name = term.split(str_)[0]

					if pathways_name not in self.map__pathways__algorithm__pval:
						self.map__pathways__algorithm__pval[pathways_name] = {}
						
					algorithm_name = algorithm_file_path.split(".")[0]
					
					if algorithm_name not in self.map__pathways__algorithm__pval[pathways_name]:
						self.map__pathways__algorithm__pval[pathways_name][algorithm_name] = -math.log10(p_value)


	def __load_sampled_reactome_enrichment_data_frame__(self,str_ = " Homo sapiens"):
		
		self.map__algorithm__pathways__score = {}
		
		for algorithm_name, file_paths in self.reactome_enrichment_alg_file_paths.items():

			counter = len(file_paths)
			map__pathways__p_val = {}
			self.map__algorithm__pathways__score[algorithm_name] = {}
			
			for file_path in file_paths:
				
				csv_reader = csv.reader(open(file_path, "r"), delimiter = "\t")

				for index, row in enumerate(csv_reader):
					if index == 0:
						continue

					term = row[1]
					p_value  = float(row[2])


					if str_ in term:
					
						pathways_name = term.split(str_)[0]

						if pathways_name not in map__pathways__p_val:
							map__pathways__p_val[pathways_name] = 0

						map__pathways__p_val[pathways_name] += p_value/float(counter)

			
			
			for k,v in map__pathways__p_val.items():
				self.map__algorithm__pathways__score[algorithm_name][k] = -math.log10(v)


	

	def plot_sampled_reactome_enrichment(self,threshold = -math.log10(0.05)):

		self.__load_sampled_reactome_enrichment_data_frame__()

		data_frame = pd.DataFrame(self.map__algorithm__pathways__score)

		data_frame = data_frame.fillna(0.0)
		data_frame = data_frame[(data_frame > threshold).any(1)]
		data_frame = data_frame[data_frame['RMA'] > threshold] 
		data_frame.to_csv('./imgs/sampled_Kegg_enrichment_MRA.tsv',sep = "\t")


	def __load_df__(self,file_path, is_float = False):
		csv_reader = csv.reader(open(file_path,"r"),delimiter = "\t")
		map__algorithm__value = {}

		for index, row in enumerate(csv_reader):
			if index == 0:
				continue
			if is_float:
				map__algorithm__value[row[1]] = float(row[2])
			else:
				map__algorithm__value[row[1]] = "$" + row[2] + "$"

		return map__algorithm__value




	def plot_mouse_phenotypes_and_drug_hub(self,):

		
		f,axes = plt.subplots(1,3,figsize = (20,5))
		
		drug_hub_random_data_frame  = pd.read_csv(self.drug_hub_pulmonary_random, sep = "\t", index_col = 0, header = 0)
		drug_hub_data_frame  = pd.read_csv(self.drug_hub_pulmonary_validation, sep = "\t", index_col = 0, header = 0)

		map__algorithm__precision = self.__load_df__(self.drug_hub_pulmonary_validation,is_float = True)
		map__algorithm__p_value = self.__load_df__(self.drug_hub_pulmonary_p_val, is_float = False)

		print(drug_hub_random_data_frame)
		sns_plot = sns.boxplot(data = drug_hub_random_data_frame,x="Algorithm", y="Score",showfliers = False,color='white',ax = axes[0])
		sns_plot.set_ylabel("Portion of target predicted")

		sns_plot = sns.scatterplot(data = drug_hub_data_frame,x = "Algorithm", y = "Score",s=200, color=".2", hue= "Algorithm", style = "Algorithm",ax = axes[0])
		axes[0].legend([],[], frameon=False)
		sns_plot.annotate(map__algorithm__p_value["DOMINO"], xy=("DOMINO", map__algorithm__precision["DOMINO"]), xytext=("DOMINO", map__algorithm__precision["DOMINO"] + 0.008),
            arrowprops=dict(arrowstyle="->",color='blue'))
		sns_plot.annotate(map__algorithm__p_value["FUMA"], xy=("FUMA", map__algorithm__precision["FUMA"]), xytext=("FUMA", map__algorithm__precision["FUMA"] + 0.005),
            arrowprops=dict(arrowstyle="->",color='blue'))
		sns_plot.annotate(map__algorithm__p_value["DEPICT"], xy=("DEPICT", map__algorithm__precision["DEPICT"]), xytext=("DEPICT", map__algorithm__precision["DEPICT"] + 0.005),
            arrowprops=dict(arrowstyle="->",color='blue'))
		sns_plot.annotate(map__algorithm__p_value["MENDELVAR"], xy=("MENDELVAR", map__algorithm__precision["MENDELVAR"]), xytext=("MENDELVAR", map__algorithm__precision["MENDELVAR"] + 0.005),
            arrowprops=dict(arrowstyle="->",color='blue'))
		sns_plot.annotate(map__algorithm__p_value["RMA"], xy=("RMA", map__algorithm__precision["RMA"]), xytext=("RMA", map__algorithm__precision["RMA"] + 0.005),
            arrowprops=dict(arrowstyle="->",color='blue'))
		

		sns_plot.set_title("Pulmonary Drug Target Prediction")


		drug_hub_random_data_frame  = pd.read_csv(self.drug_hub_random, sep = "\t", index_col = 0, header = 0)
		drug_hub_data_frame  = pd.read_csv(self.drug_hub_validation, sep = "\t", index_col = 0, header = 0)

		map__algorithm__precision = self.__load_df__(self.drug_hub_validation,is_float = True)
		map__algorithm__p_value = self.__load_df__(self.drug_hub_p_val, is_float = False)

		print(drug_hub_random_data_frame)
		sns_plot = sns.boxplot(data = drug_hub_random_data_frame,x="Algorithm", y="Score",showfliers = False,color='white',ax=axes[1])
		sns_plot.set_ylabel("Portion of target predicted")

		sns_plot = sns.scatterplot(data = drug_hub_data_frame,x = "Algorithm", y = "Score",s=200, color=".2", hue= "Algorithm", style = "Algorithm",ax=axes[1])
		
		sns_plot.set_title("Drug Target Prediction")

		axes[1].legend([],[], frameon=False)
		sns_plot.annotate(map__algorithm__p_value["DOMINO"], xy=("DOMINO", map__algorithm__precision["DOMINO"]), xytext=("DOMINO", map__algorithm__precision["DOMINO"] - 0.025),
            arrowprops=dict(arrowstyle="->",color='blue'))
		sns_plot.annotate(map__algorithm__p_value["FUMA"], xy=("FUMA", map__algorithm__precision["FUMA"]), xytext=("FUMA", map__algorithm__precision["FUMA"] - 0.025),
            arrowprops=dict(arrowstyle="->",color='blue'))
		sns_plot.annotate(map__algorithm__p_value["DEPICT"], xy=("DEPICT", map__algorithm__precision["DEPICT"]), xytext=("DEPICT", map__algorithm__precision["DEPICT"] + 0.025),
            arrowprops=dict(arrowstyle="->",color='blue'))
		sns_plot.annotate(map__algorithm__p_value["MENDELVAR"], xy=("MENDELVAR", map__algorithm__precision["MENDELVAR"]), xytext=("MENDELVAR", map__algorithm__precision["MENDELVAR"] + 0.025),
            arrowprops=dict(arrowstyle="->",color='blue'))
		sns_plot.annotate(map__algorithm__p_value["RMA"], xy=("RMA", map__algorithm__precision["RMA"]), xytext=("RMA", map__algorithm__precision["RMA"] - 0.02),
            arrowprops=dict(arrowstyle="->",color='blue'))
		
		

		random_df  = pd.read_csv(self.mouse_phenotype_random, sep = "\t", index_col = 0, header = 0)
		df  = pd.read_csv(self.mouse_phenotype_validation, sep = "\t", index_col = 0, header = 0)

		map__algorithm__precision = self.__load_df__(self.mouse_phenotype_validation,is_float = True)
		map__algorithm__p_value = self.__load_df__(self.mouse_phenotype_p_val, is_float = False)

		sns_plot = sns.boxplot(data = random_df,x="Algorithm", y="Score",showfliers = False,color='white',ax=axes[2])
		sns_plot = sns.scatterplot(data = df,x = "Algorithm", y = "Score", s=200, color=".2", hue= "Algorithm", style = "Algorithm",ax=axes[2])
		
		sns_plot.annotate(map__algorithm__p_value["DOMINO"], xy=("DOMINO", map__algorithm__precision["DOMINO"]), xytext=("DOMINO", map__algorithm__precision["DOMINO"] + 0.025),
            arrowprops=dict(arrowstyle="->",color='blue'))
		sns_plot.annotate(map__algorithm__p_value["FUMA"], xy=("FUMA", map__algorithm__precision["FUMA"]), xytext=("FUMA", map__algorithm__precision["FUMA"] + 0.025),
            arrowprops=dict(arrowstyle="->",color='blue'))
		sns_plot.annotate(map__algorithm__p_value["DEPICT"], xy=("DEPICT", map__algorithm__precision["DEPICT"]), xytext=("DEPICT", map__algorithm__precision["DEPICT"] + 0.025),
            arrowprops=dict(arrowstyle="->",color='blue'))
		sns_plot.annotate(map__algorithm__p_value["MENDELVAR"], xy=("MENDELVAR", map__algorithm__precision["MENDELVAR"]), xytext=("MENDELVAR", map__algorithm__precision["MENDELVAR"] + 0.03),
            arrowprops=dict(arrowstyle="->",color='blue'))
		sns_plot.annotate(map__algorithm__p_value["RMA"], xy=("RMA", map__algorithm__precision["RMA"]), xytext=("RMA", map__algorithm__precision["RMA"] + 0.025),
            arrowprops=dict(arrowstyle="->",color='blue'))
		
		sns_plot.set_title("Mouse Phenotypes Prediction")
		for n, ax in enumerate(axes):

			ax.text(-0.1, 1.1, string.ascii_uppercase[n], transform=ax.transAxes, 
            size=20, weight='bold')

		f.savefig('./imgs/external_validation.pdf',bbox_inches='tight')


	def plot_mouse_phenotypes_and_drug_hub_on_randomized_graph(self,):
		
		f,axes = plt.subplots(1,2,figsize = (20,5))

		random_df  = pd.read_csv(self.mouse_phenotype_random_bg_distribution, sep = "\t", index_col = 0, header = 0)
		
		map__algorithm__precision = self.__load_df__(self.mouse_phenotype_validation,is_float = True)
		map__algorithm__p_value = self.__load_df__(self.mouse_phenotype_random_bg_p_value, is_float = True)

		sns_plot = sns.histplot(data = random_df, x = "Random Distribution", ax = axes[0])
		sns_plot.annotate("$P_{val}$: " + str(map__algorithm__p_value["RMA"]), xy=(map__algorithm__precision["RMA"], 0.5), xytext=(map__algorithm__precision["RMA"], 10),
            arrowprops=dict(arrowstyle="->",color='blue'))
		sns_plot.set_title("Mouse Phenotypes Prediction")

		random_df  = pd.read_csv(self.drug_hub_random_bg_distribution, sep = "\t", index_col = 0, header = 0)
		
		map__algorithm__precision = self.__load_df__(self.drug_hub_validation,is_float = True)
		map__algorithm__p_value = self.__load_df__(self.drug_hub_random_bg_p_value, is_float = True)

		sns_plot = sns.histplot(data = random_df, x = "Random Distribution",ax = axes[1])
		sns_plot.annotate("$P_{val}$: " + str(map__algorithm__p_value["RMA"]), xy=(map__algorithm__precision["RMA"], 0.5), xytext=(map__algorithm__precision["RMA"], 10),
            arrowprops=dict(arrowstyle="->",color='blue'))
		sns_plot.set_title("Mouse Phenotypes Prediction")

		for n, ax in enumerate(axes):

			ax.text(-0.1, 1.1, string.ascii_uppercase[n], transform=ax.transAxes, 
            size=20, weight='bold')

		f.savefig('./imgs/external_validation_random_biipartite_graph.pdf',bbox_inches='tight')


	def plot_odd_vs_even(self, option = ['55%','65%','75%','85%','95%']):
		
		f,axes = plt.subplots(2,2,figsize = (16,9))


		bma_random_data_frame  = pd.read_csv(self.odd_vs_even_bma_random, sep = "\t", index_col = 0, header = 0)
		bma_data_frame  = pd.read_csv(self.odd_vs_even_bma_validation, sep = "\t", index_col = 0, header = 0)

		map__algorithm__precision = self.__load_df__(self.odd_vs_even_bma_validation,is_float = True)
		map__algorithm__p_value = self.__load_df__(self.odd_vs_even_bma_p_values, is_float = True)

		sns_plot = sns.boxplot(data = bma_random_data_frame,x="Algorithm", y="Score",showfliers = False,color='white',ax = axes[0][0])
		sns_plot = sns.scatterplot(data = bma_data_frame,x = "Algorithm", y = "Score",s=200, color=".2", hue= "Algorithm", style = "Algorithm",ax = axes[0][0])		
		
		sns_plot.annotate("$P_{val} \sim$" + str(map__algorithm__p_value["DOMINO"]), xy=("DOMINO", map__algorithm__precision["DOMINO"]), xytext=("DOMINO", map__algorithm__precision["DOMINO"] + 0.05),
            arrowprops=dict(arrowstyle="->",color='blue'))
		
		sns_plot.annotate("$P_{val} \sim$" +str(map__algorithm__p_value["FUMA"]), xy=("FUMA", map__algorithm__precision["FUMA"]), xytext=("FUMA", map__algorithm__precision["FUMA"] + 0.05),
            arrowprops=dict(arrowstyle="->",color='blue'))
		
		sns_plot.annotate("$P_{val} \sim$" +str(map__algorithm__p_value["MENDELVAR"]), xy=("MENDELVAR", map__algorithm__precision["MENDELVAR"]), xytext=("MENDELVAR", map__algorithm__precision["MENDELVAR"] - 0.05),
            arrowprops=dict(arrowstyle="->",color='blue'))
		
		sns_plot.annotate("$P_{val} \sim$" +str(map__algorithm__p_value["RMA"]), xy=("RMA", map__algorithm__precision["RMA"]), xytext=("RMA", map__algorithm__precision["RMA"] + 0.05),
            arrowprops=dict(arrowstyle="->",color='blue'))
		

		sns_plot.set_title("Odd versus Even Predicted gene sets")
		sns_plot.set_ylabel("Biological Similarity")
		sns_plot.set_xlabel(" ")

		np_random_data_frame  = pd.read_csv(self.odd_vs_even_np_random, sep = "\t", index_col = 0, header = 0)
		np_data_frame  = pd.read_csv(self.odd_vs_even_np_validation, sep = "\t", index_col = 0, header = 0)
		map__algorithm__precision = self.__load_df__(self.odd_vs_even_np_validation,is_float = True)
		map__algorithm__p_value = self.__load_df__(self.odd_vs_even_np_p_values, is_float = True)

		
		sns_plot = sns.boxplot(data = np_random_data_frame,x="Algorithm", y="Score",showfliers = False,color='white',ax= axes[1][0])
		sns_plot = sns.scatterplot(data = np_data_frame,x = "Algorithm", y = "Score",s=200, color=".2", hue= "Algorithm", style = "Algorithm",ax= axes[1][0])
		
		sns_plot.annotate("$P_{val} \sim$" + str(map__algorithm__p_value["DOMINO"]), xy=("DOMINO", map__algorithm__precision["DOMINO"]), xytext=("DOMINO", map__algorithm__precision["DOMINO"] + 0.05),
            arrowprops=dict(arrowstyle="->",color='blue'))
		
		sns_plot.annotate("$P_{val} \sim$" +str(map__algorithm__p_value["FUMA"]), xy=("FUMA", map__algorithm__precision["FUMA"]), xytext=("FUMA", map__algorithm__precision["FUMA"] + 0.05),
            arrowprops=dict(arrowstyle="->",color='blue'))
		
		sns_plot.annotate("$P_{val} \sim$" +str(map__algorithm__p_value["MENDELVAR"]), xy=("MENDELVAR", map__algorithm__precision["MENDELVAR"]), xytext=("MENDELVAR", map__algorithm__precision["MENDELVAR"] - 0.05),
            arrowprops=dict(arrowstyle="->",color='blue'))
		
		sns_plot.annotate("$P_{val} \sim$" +str(map__algorithm__p_value["RMA"]), xy=("RMA", map__algorithm__precision["RMA"]), xytext=("RMA", map__algorithm__precision["RMA"] + 0.05),
            arrowprops=dict(arrowstyle="->",color='blue'))
		

		sns_plot.set_title("Odd versus Even Predicted gene sets")
		sns_plot.set_ylabel("Network Distance")
		


		bma_missing_data_df = pd.read_csv(self.bma_missing_data_validation, sep = "\t", index_col = 0, header = 0)
		bma_missing_data_df = bma_missing_data_df.sort_values(by=['Percentage of remaining SNPs'])
		bma_missing_data_df = bma_missing_data_df[bma_missing_data_df['Percentage of remaining SNPs'].isin(option)] 

		sns_plot = sns.boxplot(data = bma_missing_data_df,x="Percentage of remaining SNPs", y="Biological Similarity",showfliers = True,color='white',ax= axes[0][1])
		sns_plot.set_xlabel("")
		sns_plot.set_title("Robustness analysis to missing data")

		np_missing_data_df = pd.read_csv(self.np_missing_data_validation, sep = "\t", index_col = 0, header = 0)
		np_missing_data_df = np_missing_data_df.sort_values(by=['Percentage of remaining SNPs'])
		np_missing_data_df = np_missing_data_df[np_missing_data_df['Percentage of remaining SNPs'].isin(option)] 

		sns_plot = sns.boxplot(data = np_missing_data_df,x="Percentage of remaining SNPs", y="Network Distance",showfliers = True,color='white',ax= axes[1][1])
	
		sns_plot.set_title("Robustness analysis to missing data")

		for n, ax in enumerate(axes.flat):

			ax.text(-0.1, 1.1, string.ascii_uppercase[n], transform=ax.transAxes, 
            size=20, weight='bold')

		f.savefig('./imgs/internal_validation.pdf',bbox_inches='tight')


p = Plot(disease_experiment_dir_path = "./experiments/algorithm_comparison/")

p.plot_mouse_phenotypes_and_drug_hub()
p.plot_mouse_phenotypes_and_drug_hub_on_randomized_graph()

p = Plot(disease_experiment_dir_path = "./experiments/missing_data/")
p.plot_odd_vs_even()