import csv
import networkx as nx
import pandas as pd
import seaborn as sns
import random
import matplotlib.pyplot as plt
import os

class OddVsEven():
	
	def __init__(self, 
		disease_experiment_dir_path,
		solution_without_missing_information_file_path, 
		network_file_path,
		ontology_file_path,
		locus_gene_file_path,
		filter_):

		self.ontology_file_path =  ontology_file_path

		self.solution_without_missing_information_file_path = solution_without_missing_information_file_path
		self.disease_experiment_dir_path = disease_experiment_dir_path

		self.algorithms_file_path = disease_experiment_dir_path + "algorithms/"
		self.validation_dir_path = disease_experiment_dir_path + "validation/"

		if not os.path.exists(self.validation_dir_path):
			os.makedirs(self.validation_dir_path)

		self.filter_ = filter_
		self.locus_gene_file_path = locus_gene_file_path
		
		self.chromosome_splitting_dir_path = disease_experiment_dir_path + "chromosome_splitting/"
		
		self.__load_locus_genes__()
		self.map__algorithm__solution = {dir_: (self.__load_solution__(self.algorithms_file_path + dir_ + "/even.txt"),self.__load_solution__(self.algorithms_file_path + dir_ + "/odd.txt")) for dir_ in os.listdir(self.algorithms_file_path) if dir_[0] != "." and dir_ not in self.filter_}

		self.network_file_path = network_file_path

		self.odd_loci = {
			'rs7650602','rs11118406','rs12185268','rs10760580','rs798565','rs11655567','rs2040732','rs72673419',
			'rs10866659','rs156394','rs2096468','rs2897075','rs34651','rs7866939','rs34727469','rs55676755',
			'rs4757118','rs76841360','rs3009947','rs4093840','rs629619','rs1529672','rs13073544','rs7642001',
			'rs9435731','rs1441358','rs9525927','rs10152300','rs979453','rs1551943','rs10114763','rs12519165',
			'rs72626215','rs62259026','rs8080772','rs62065216','rs117261012','rs10037493','rs4660861','rs17759204',
			'rs2442776','rs803923','rs62375246','rs72731149','rs2955083','rs11579382','rs153916'
			
		}

		self.even_loci = {
			'rs34712979','rs1570221','rs7671261','rs7307510','rs8044657','rs9617650','rs7068966','rs4888379',
			'rs1334576','rs73158393','rs2571445','rs7958945','rs10929386','rs647097','rs646695','rs16825267',
			'rs13198656','rs1631199','rs13140176','rs3095329','rs2070600','rs56134392','rs9399401','rs4585380',
			'rs62191105','rs955277','rs2806356','rs11049386','rs721917','rs9329170','rs9350191','rs2579762',
			'rs12466981','rs72699855','rs72902175'
		}
		self.__load_chromosome_splitting__()
			
		self.__load_network__()
		self.__load_gene_ontology__()
		self.__precompute_shortest_path_distance__()


	def __load_chromosome_splitting__(self,threshold = 100):
		
		percentuage_files = [self.chromosome_splitting_dir_path + file for file in os.listdir(self.chromosome_splitting_dir_path) if file[0] != "." and ".tsv" in file]
		self.map__percentuage__solutions = {}
		
		for file in percentuage_files:
			
			csv_reader = csv.reader(open(file,"r"),delimiter = "\t")
			
			for index, row in enumerate(csv_reader):
				
				if index == 0:
					continue

				if index > threshold:
					break
				percenage = row[0].split("__")[-3]
				solution = set(row[-1].replace("{","").replace("}","").replace(" ","").replace("'","").split(","))
				
				if percenage not in self.map__percentuage__solutions:
					self.map__percentuage__solutions[percenage] = []

				self.map__percentuage__solutions[percenage].append(solution)

		return self.map__percentuage__solutions
	

	def __load_solution__(self,file_path):
		solution = set()

		csv_reader = csv.reader(open(file_path, "r"), delimiter = "\t")

		for row in csv_reader:
			if row[0] in self.universe:
				solution.add(row[0])
	
		return solution


	def __load_gene_ontology__(self,
		current_domain = "cellular_component"
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

	def __load_locus_genes__(self,):

		self.map__locus__genes = {}
		self.map__gene__locus = {}
		
		with open(self.locus_gene_file_path, "r") as fp:
			csv_reader = csv.reader(fp, delimiter = "\t")
			for row in csv_reader:
				gene_name = row[0]
				locus = row[1]

				if locus not in self.map__locus__genes:
					self.map__locus__genes[locus] = set()

				
				self.map__locus__genes[locus].add(gene_name)
				self.map__gene__locus[gene_name] = locus


		self.universe = set(self.map__gene__locus.keys())

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


	def __precompute_shortest_path_distance__(self):
		
		universe = set(self.map__gene__locus.keys())
		self.map__node__pair__sp_l = {}

		for index,i in enumerate(universe):
			for j in universe:
				
				try:
					sp_l = nx.shortest_path_length(self.G, source=i, target=j)
					self.map__node__pair__sp_l[(i,j)] = sp_l
					
				except:
					pass
	
	def __create_random_solution__(self,):
		
		random_solution_even = set()
		random_solution_odd = set()

		for locus,genes in self.map__locus__genes.items():
			random_gene = random.sample(genes,1)[0]

			if locus in self.odd_loci:
				random_solution_odd.add(random_gene)

			if locus in self.even_loci:
				random_solution_even.add(random_gene)

		return random_solution_odd, random_solution_even


	def __create_random_solution_all_algorithm__(self, odd_size, even_size):

		odd_genes = set()
		even_genes = set()

		random_solution_odd = set()
		random_solution_even = set()
		
		for locus in self.odd_loci:
			odd_genes = odd_genes.union(self.map__locus__genes[locus])

		for locus in self.even_loci:
			even_genes = even_genes.union(self.map__locus__genes[locus])
			
		random_solution_odd = set(random.sample(odd_genes,odd_size))

		random_solution_even = set(random.sample(even_genes,even_size))

				
		return random_solution_odd, random_solution_even



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

	def __compute_pval_biological_similarity__(self,alg,even,odd,trial = 100):
		
		bma = self.__BMA__( even, odd)
		counter = 0
		table = []
		
		for i in range(trial):

			if alg != "RMA":
				random_solution_odd, random_solution_even = self.__create_random_solution_all_algorithm__(len(odd),len(even))
			else:
				random_solution_odd, random_solution_even = self.__create_random_solution__()

			bma_i = self.__BMA__( random_solution_even, random_solution_odd)

			table.append([alg,bma_i])

			if bma_i > bma:
				counter += 1

		print("p_val distance between odd and even", alg,counter / trial)
		p_val = counter / trial
		return table, bma, p_val


	def __compute_pval_network_proximity__(self,alg,even,odd,trial = 100):
		
		phi = self.__get_network_proximity__( even, odd)
		
		counter = 0
		table = []

		for i in range(trial):
			if alg != "RMA":
				random_solution_odd, random_solution_even = self.__create_random_solution_all_algorithm__(len(odd),len(even))
			else:
				random_solution_odd, random_solution_even = self.__create_random_solution__()
			
			phi_i = self.__get_network_proximity__( random_solution_odd, random_solution_even)
			table.append([alg,phi_i])

			if phi_i < phi:
				counter += 1

		print("p_val distance between odd and even", alg,counter / trial)

		p_val = counter / trial
		return table, phi, p_val

			

	def run(self,):
		
		bmas = []
		bmas_r = []
		bmas_pvals = []

		net_prox = []
		net_prox_r = []
		net_prox_pvals = []
		for alg, solutions in self.map__algorithm__solution.items():
			even,odd = solutions
			
			df, phi,p_val = self.__compute_pval_network_proximity__(alg,even,odd)
			net_prox.extend(df)
			net_prox_r.append([alg,phi])
			net_prox_pvals.append([alg,p_val])

			df, phi,p_val = self.__compute_pval_biological_similarity__(alg,even,odd)
			bmas.extend(df)
			bmas_r.append([alg,phi])
			bmas_pvals.append([alg, p_val])


		pd.DataFrame(bmas,columns = ["Algorithm","Score"]).to_csv(self.validation_dir_path + "bma_distribution.tsv", sep = "\t")
		pd.DataFrame(bmas_r,columns = ["Algorithm","Score"]).to_csv(self.validation_dir_path + "bma_values.tsv", sep = "\t")
		pd.DataFrame(bmas_pvals,columns = ["Algorithm","p_val"]).to_csv(self.validation_dir_path + "bma_p_values.tsv", sep = "\t")

		pd.DataFrame(net_prox,columns = ["Algorithm","Score"]).to_csv(self.validation_dir_path + "network_proximity_distribution.tsv", sep = "\t")
		pd.DataFrame(net_prox_r,columns = ["Algorithm","Score"]).to_csv(self.validation_dir_path + "network_proximity_values.tsv", sep = "\t")
		pd.DataFrame(net_prox_pvals,columns = ["Algorithm","p_val"]).to_csv(self.validation_dir_path + "network_proximity_p_values.tsv", sep = "\t")

	def run_chromosome_splitting(self,):
		solution = self.__load_solution__(self.solution_without_missing_information_file_path)
		table_1 = []
		table_2 = []
		
		for parcentage, solutions in self.map__percentuage__solutions.items():

			for index, current_solution in enumerate(solutions):
			
				JC = len(solution.intersection(current_solution))/len(solution.union(current_solution))
				
				phi_network_proximity = self.__get_network_proximity__( solution, current_solution)
				phi_biological_similarity = self.__BMA__(solution, current_solution)
				
				print(parcentage,JC)

				table_1.append([parcentage + "%", phi_network_proximity])
				table_2.append([parcentage + "%", phi_biological_similarity])
				
		pd.DataFrame(table_1,columns = ["Percentage of remaining SNPs","Network Distance"]).to_csv(self.validation_dir_path + "network_distance_distribution_missing_data.tsv", sep = "\t")

		pd.DataFrame(table_2,columns = ["Percentage of remaining SNPs","Biological Similarity"]).to_csv(self.validation_dir_path +"bma_distribution_missing_data.tsv" , sep = "\t")
		
ovse = OddVsEven(

	disease_experiment_dir_path = "../../experiments/missing_data/",
	solution_without_missing_information_file_path = "../../experiments/algorithm_comparison/algorithms/RMA/solution.txt", 
	network_file_path = "../../datasets/PPI_network/network.tsv",
	ontology_file_path = "../../datasets/curated_db/GO_10_02_2021.txt",
	locus_gene_file_path = "../../datasets/GWAS/unbiased_copd_snp_database.txt",
	filter_ = []
)
ovse.run()
ovse.run_chromosome_splitting()