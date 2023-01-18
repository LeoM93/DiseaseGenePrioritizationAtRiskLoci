import csv
import networkx as nx
import pandas as pd
import seaborn as sns
import random
import os
import seaborn as sns
import matplotlib.pyplot as plt
import string
class CrossStudyValidation():

	def __init__(self,
		ontology_file_path,
		network_file_path,
		gwas_catalog_file_path,
		solution_file_path,
		disease_dir		
		):

		self.ontology_file_path =  ontology_file_path
		self.network_file_path = network_file_path
		self.gwas_catalog_file_path = gwas_catalog_file_path
	

		self.solution_file_path = solution_file_path
		self.disease_dir = disease_dir

		
	def __load_solution__(self,):
		csv_reader = csv.reader(open(self.solution_file_path,"r"),delimiter = "\t")
		self.map__study__solution = {}
		self.universe = set()
		
		for index, row in enumerate(csv_reader):

			if index == 0:
				continue

			name = row[0].split("/")[-2]

			solution = set(row[-2].replace("{","").replace("}","").replace(" ","").replace("'","").split(","))
			self.universe = self.universe.union(solution)
			self.map__study__solution[name] = solution


	def __load_study_information__(self,):
		
		dirs = [self.disease_dir + dir_ + "/" for dir_ in os.listdir(self.disease_dir) if os.path.isdir(self.disease_dir + dir_ + "/")]
		self.map__study__record = {}
		for dir_ in dirs:
			
			study_name = dir_.split("/")[-2]

			network_file_path = dir_ + "biological_process_co_regulation_network_1e-2_thresholded.txt"

			V,E,risk_loci = self.__load_co_regulation_network__(network_file_path)
			
			self.map__study__record[study_name] = (V, E, risk_loci)
	

	def __load_gwas_catalog__(self):
		csv_reader = csv.reader(open(self.gwas_catalog_file_path , "r"),delimiter = "\t")
		self.map__study_id__disease = {}
		for index, row in enumerate(csv_reader):
			if index ==0:
				continue

			study_id = row[-2]
			self.map__study_id__disease[study_id] =  row[-8]


	def __load_gene_ontology__(self,
		current_domain = "biological_process"
	):

		self.map__ensembl__id__ontologies = {}

		with open(self.ontology_file_path,'r') as fp:
			csv_reader = csv.reader(fp, delimiter = "\t")
			
			for index,row in enumerate(csv_reader):

				if index == 0:
					continue
					
				gene_name = row[0]
				term_id = row[1]
				domain = row[3]


				if domain != current_domain:
					continue

				if gene_name not in self.map__ensembl__id__ontologies:
					self.map__ensembl__id__ontologies[gene_name] = set()

				self.map__ensembl__id__ontologies[gene_name].add(term_id)



	def __load_network__(self,has_header = True):

		self.G = nx.Graph()

		with open(self.network_file_path, "r") as fp:
			csv_reader = csv.reader(fp, delimiter = "\t")
			
			for index,row in enumerate(csv_reader):
				
				if has_header and index == 0:
					continue

				source = row[0]
				target = row[1]

				self.G.add_edge(source, target)

		return self.G

	def __load_co_regulation_network__(self,network_file_path,has_header = True):

		G = nx.Graph()
		risk_loci = set()
		with open(network_file_path, "r") as fp:
			csv_reader = csv.reader(fp, delimiter = "\t")
			
			for index,row in enumerate(csv_reader):
				
				if has_header and index == 0:
					continue

				source = row[0]
				target = row[1]
				snp_1 = row[2]
				snp_2 = row[3]

				risk_loci.add(snp_1)
				risk_loci.add(snp_2)

				G.add_edge(source, target)

		return set(G.nodes()),set(G.edges()),risk_loci


	def __precompute_shortest_path_distance__(self):
		
		
		self.map__node__pair__sp_l = {}

		for index,i in enumerate(self.universe):
			for j in self.universe:
				
				try:
					sp_l = nx.shortest_path_length(self.G, source=i, target=j)
					self.map__node__pair__sp_l[(i,j)] = sp_l
					
				except:
					pass

	def __get_sp_dist_between_set__(self, set_1, set_2):

		size = 0
		mean = 0.0
		total_size = 0
		for node_1 in set_1:
			
			for node_2 in set_2:
				
				if (node_1, node_2) in self.map__node__pair__sp_l:
					
					mean += self.map__node__pair__sp_l[node_1, node_2]
					total_size += 1
				

		if total_size != 0:
			mean = mean/total_size
		else:
			mean = 6
		
		return mean

	def __get_network_proximity__(self,set_1, set_2):

	
		s_a_b = self.__get_sp_dist_between_set__( set_1, set_2)
		s_a_a = self.__get_sp_dist_between_set__( set_1, set_1)
		s_b_b = self.__get_sp_dist_between_set__( set_2, set_2)


		return s_a_b - (s_a_a + s_b_b)/2

	def __BMA__(self, set_1, set_2):
	
		set_1_scores = []
		set_2_scores = []

		#
		map__biological_processes = {}

		for item_1 in set_1:
			max_score = 0.0
			
			for item_2 in set_2:

				if item_1 in self.map__ensembl__id__ontologies and item_2 in self.map__ensembl__id__ontologies:

					intersection = len(self.map__ensembl__id__ontologies[item_1].intersection(self.map__ensembl__id__ontologies[item_2]))
					union = len(self.map__ensembl__id__ontologies[item_1].union(self.map__ensembl__id__ontologies[item_2]))

					score = intersection/union

					if score > max_score:
						
						max_score = score


			set_1_scores.append(max_score)

		for item_1 in set_2:
			max_score = 0.0

			for item_2 in set_1:

				if item_1 in self.map__ensembl__id__ontologies and item_2 in self.map__ensembl__id__ontologies:

					intersection = len(self.map__ensembl__id__ontologies[item_1].intersection(self.map__ensembl__id__ontologies[item_2]))
					union = len(self.map__ensembl__id__ontologies[item_1].union(self.map__ensembl__id__ontologies[item_2]))

					score = intersection/union

					if score > max_score:
						max_score = score

			set_2_scores.append(max_score)


		sum_1 = sum(set_1_scores)
		sum_2 = sum(set_2_scores)

		return 0.5*(sum_1/len(set_1) + sum_2/len(set_2))


	def visualize_bma_np(self,ax_1, ax_2, ax_3):

		self.__load_solution__()
		self.__load_network__()
		self.__load_gene_ontology__()
		self.__precompute_shortest_path_distance__()


		studies = list(self.map__study__solution.keys())
		map__study_1__study_2__bma_score = {}
		map__study_1__study_2__np_score = {}
		map__study_1__study_2__jaccard_solution = {}


		for study_1, solution_1 in self.map__study__solution.items():
			
			map__study_1__study_2__bma_score[study_1] = {}
			map__study_1__study_2__np_score[study_1] = {}
			map__study_1__study_2__jaccard_solution[study_1] = {}

			for study_2, solution_2 in self.map__study__solution.items():

				bma = self.__BMA__(solution_1,solution_2)
				np = self.__get_network_proximity__(solution_1,solution_2)

				intersection = len(solution_1.intersection(solution_2))
				union = len(solution_1.union(solution_2))

				map__study_1__study_2__bma_score[study_1][study_2] = bma
				map__study_1__study_2__np_score[study_1][study_2]= np
				map__study_1__study_2__jaccard_solution[study_1][study_2] = intersection/union

		data_frame = pd.DataFrame(map__study_1__study_2__jaccard_solution)
		
		#data_frame = data_frame.reindex(self.filter_)
		#data_frame = data_frame.reindex(self.filter_, axis="columns")
		
		sns_plot = sns.heatmap(data_frame,annot = True, ax= ax_1 )
		sns_plot.set_title("Jaccard Index of Solution")

		data_frame = pd.DataFrame(map__study_1__study_2__bma_score)
		
		#data_frame = data_frame.reindex(self.filter_)
		#data_frame = data_frame.reindex(self.filter_, axis="columns")

		sns_plot = sns.heatmap(data_frame,annot = True,ax= ax_2 )
		sns_plot.set_title("Biological Similarity")


		data_frame = pd.DataFrame(map__study_1__study_2__np_score)
		
		#data_frame = data_frame.reindex(self.filter_)
		#data_frame = data_frame.reindex(self.filter_, axis="columns")

		sns_plot = sns.heatmap(data_frame,annot = True, ax= ax_3)
		sns_plot.set_title("Network Distance")


	def compare_dataset(self,ax):
		
		self.__load_study_information__()
		self.__load_gwas_catalog__()

		
		summary_table = []
		total_risk_loci = set()
		map__study_1__study_2__jaccard = {}
		risk_locus_distribution = []
		
		for k in self.map__study__record.keys():
			V, E, risk_loci = self.map__study__record[k]
			total_risk_loci = total_risk_loci.union(risk_loci)
			risk_locus_distribution.append(len(risk_loci))
			
			summary_table.append([k,self.map__study_id__disease[k],len(risk_loci),len(V),len(E)])
		
		summary_table.sort(key = lambda x: x[2])
		summary_data_frame = pd.DataFrame(summary_table, columns = ["Study ID","Population","N. of Risk Loci","N. of Nodes", "N. of Edges"])
		risk_locus_distribution_data_frame = pd.DataFrame(risk_locus_distribution, columns = ["Risk Locus Distribution"])
		
		new_index = [i[0] for i in summary_table]
				
		for study_1, record_1 in self.map__study__record.items():
			map__study_1__study_2__jaccard[study_1] = {}
			V_1, E_1, risk_loci_1 = record_1
			
			for study_2, record_2 in self.map__study__record.items():
				V_2, E_2, risk_loci_2 = record_2
				intersection_V = len(V_1.intersection(V_2))
				intersection_E = len(E.intersection(E))
				intersection_rl = len(risk_loci_1.intersection(risk_loci_2))

				union_V = len(V_1.union(V_2))
				union_E = len(E.union(E))
				union_rl = len(risk_loci_1.union(risk_loci_2))

				map__study_1__study_2__jaccard[study_1][study_2] = intersection_rl/union_rl

		data_frame = pd.DataFrame(map__study_1__study_2__jaccard)
		data_frame = data_frame.reindex(new_index)
		data_frame = data_frame.reindex(new_index, axis="columns")

		sns_plot = sns.heatmap(data_frame,annot = True, ax= ax)
		sns_plot.set_title("Jaccard Index of SNPs")

			

	def run(self,):
		f,axes = plt.subplots(2,2,figsize=(20,20))
		self.compare_dataset(axes[0][0],)
		self.visualize_bma_np(axes[0][1],axes[1][0],axes[1][1])

		new_axes = axes.flat
		for n, ax in enumerate(new_axes):

			ax.text(-0.1, 1.1, string.ascii_uppercase[n], transform=ax.transAxes, 
            size=20, weight='bold')
		
		f.savefig('../../imgs/cross_validation.pdf',bbox_inches='tight')

c = CrossStudyValidation(
	ontology_file_path =  "../../datasets/curated_db/GO_10_02_2021.txt",
	network_file_path = "../../datasets/PPI_network/network.tsv",
	gwas_catalog_file_path= "../../datasets/GWAS/gwas_catalog_v1.0.2-studies_r2022-06-15.tsv",
	solution_file_path = "../../experiments/breast_cancer/algorithms/RMA/RMA.tsv",
	disease_dir = "../../datasets/GWAS/breast_cancer/",
	)
c.run()
