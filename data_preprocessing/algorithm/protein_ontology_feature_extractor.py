import csv
import networkx as nx
import os
from algorithm.ontology_similarity_measures.resnik_sim import ResnikSimilarity

class ProteinOntologyFeatureExtractor():
	
	def __init__(self,
		
		ontology_file_path,
		ontology_graph_file_path,
		network_file_path,
		output_file_path,

		similarity_score_def = "resnik_similarity"
		
		):
		
		self.ontology_file_path = ontology_file_path
		self.ontology_graph_file_path = ontology_graph_file_path
		self.network_file_path = network_file_path


		self.output_file_path = output_file_path
		self.similarity_score_def = similarity_score_def

		
	
	def __load_V__(self,):
		
		V = set()
		E = set()
		
		with open(self.network_file_path,'r') as fp:
			csv_reader = csv.reader(fp, delimiter = "\t")
			for index,row in enumerate(csv_reader):

				if index == 0:
					continue
				
				u = row[0]
				v = row[1]

				V.add(u)
				V.add(v)

				E.add((u,v))
		
		return E,V


	def __load_acyclic_network__(self,):
		
		G = nx.DiGraph()
		
		with open(self.ontology_graph_file_path,'r') as fp:
			csv_reader = csv.reader(fp, delimiter = "\t")
			for index,row in enumerate(csv_reader):
				
				u = row[0]
				v = row[1]


				G.add_edge(u,v,cagegory = row[2])
		
		V = set(G.nodes())

		return G,V


	def __load_map__ontology_domain__ontologies__(self,V,ES = "IEA"):
		
		map__ontology_domain__ontologies = {}	
		map__ontology_domain__ensembl__id__ontologies = {}
		with open(self.ontology_file_path,'r') as fp:
			
			csv_reader = csv.reader(fp, delimiter = "\t")
			for index,row in enumerate(csv_reader):

				if index == 0:
					continue
				gene_stable_id = row[0]
				accession_number = row[4]
				evidence_score = row[6]
				domain = row[7]
				
				if gene_stable_id in V:
					if domain != '':

						if evidence_score != ES:

							if domain not in map__ontology_domain__ontologies:
								map__ontology_domain__ontologies[domain] = set()

							map__ontology_domain__ontologies[domain].add(accession_number)

							if domain not in map__ontology_domain__ensembl__id__ontologies:
								map__ontology_domain__ensembl__id__ontologies[domain] = {}

							if gene_stable_id not in map__ontology_domain__ensembl__id__ontologies[domain]:
								map__ontology_domain__ensembl__id__ontologies[domain][gene_stable_id] = set()

							map__ontology_domain__ensembl__id__ontologies[domain][gene_stable_id].add(accession_number)

		return map__ontology_domain__ontologies,map__ontology_domain__ensembl__id__ontologies


	def __load_map__ontology_domain__ontologies_for_a_better_version__(self,V):
		
		map__ontology_domain__ontologies = {}	
		map__ontology_domain__ensembl__id__ontologies = {}
		with open(self.ontology_file_path,'r') as fp:
			
			csv_reader = csv.reader(fp, delimiter = "\t")
			for index,row in enumerate(csv_reader):

				if index == 0:
					continue
				gene_stable_id = row[0]
				accession_number = row[1]
				domain = row[3]
				
				if gene_stable_id in V:
					if domain != '':

						if domain not in map__ontology_domain__ontologies:
							map__ontology_domain__ontologies[domain] = set()

						map__ontology_domain__ontologies[domain].add(accession_number)

						if domain not in map__ontology_domain__ensembl__id__ontologies:
							map__ontology_domain__ensembl__id__ontologies[domain] = {}

						if gene_stable_id not in map__ontology_domain__ensembl__id__ontologies[domain]:
							map__ontology_domain__ensembl__id__ontologies[domain][gene_stable_id] = set()

						map__ontology_domain__ensembl__id__ontologies[domain][gene_stable_id].add(accession_number)

		return map__ontology_domain__ontologies,map__ontology_domain__ensembl__id__ontologies




	def __save_acyclic_GO_graph__(self,ontology_obo_file_path):
		csv_writer = csv.writer(open(self.ontology_graph_file_path,"w"),delimiter = "\t")

		with open(ontology_obo_file_path, "r") as fp:
			csv_reader = csv.reader(fp, delimiter = "\n")
			
			loading_flag = False
			map__term__description ={}
			
			for row in csv_reader:
				
				if len(row) == 0:


					if loading_flag:
						
						if "is_a" in map__term__description:
							parents = map__term__description["is_a"]

							for parent in parents:
								csv_writer.writerow([map__term__description["id"],parent, map__term__description["namespace"]])

						if "part_of" in map__term__description:
							parents = map__term__description["part_of"]

							for parent in parents:
								csv_writer.writerow([map__term__description["id"],parent, map__term__description["namespace"]])


					loading_flag = False
					map__term__description ={}
					continue

				line = row[0]

				if "[Term]" == line:
					loading_flag = True

				if loading_flag:
					if "id: " in line and "alt_id: " not in line:

						
						map__term__description["id"] = line.split(" ")[-1]
					if "namespace: " in line:
						map__term__description["namespace"] = line.split(" ")[-1]

					if "is_a: " in line:
						if "is_a" not in map__term__description:
							map__term__description["is_a"] = set()
							
						map__term__description["is_a"].add(line.split(" ")[1])

					if "relationship: part_of" in line:
						if "part_of" not in map__term__description:
							map__term__description["part_of"] = set()
							
						map__term__description["part_of"].add(line.split(" ")[2])


	def __load_term_similarity_dictionary__(self,file_paths,ontology_domains):
		
		map__ontology_domain__term_pair__score = {}

		for file_path,ontology_domain in zip(file_paths,ontology_domains):

			assert ontology_domain in file_path
			
			map__ontology_domain__term_pair__score[ontology_domain] = {}
			with open(file_path,"r") as fp:
				csv_reader = csv.reader(fp,delimiter = "\t")

				for index, row in enumerate(csv_reader):
					if index == 0:
						continue

					term_1 = row[0]
					term_2 = row[1]
					score = float(row[2])
					map__ontology_domain__term_pair__score[ontology_domain][(term_1,term_2)] = score

		return map__ontology_domain__term_pair__score



	def __get_best_match_average__(self, p_T, q_T,ontology_domain):

		score = 0
		n = float(len(p_T))
		m = float(len(q_T))

		score_1 = 0
		score_2 = 0

		for t in p_T:
			max_ = 0
			
			for q in q_T:
				pair_score = 0.0

				if (t,q) in self.map__ontology_domain__term_pair__score[ontology_domain]:
					pair_score = self.map__ontology_domain__term_pair__score[ontology_domain][(t,q)]

				elif (q,t) in self.map__ontology_domain__term_pair__score[ontology_domain]:
					pair_score = self.map__ontology_domain__term_pair__score[ontology_domain][(q,t)]

				if pair_score > max_:
						max_ = pair_score

			score_1 += max_

		for q in q_T:
			
			max_ = 0
			
			for t in p_T:
				pair_score = 0.0

				if (t,q) in self.map__ontology_domain__term_pair__score[ontology_domain]:
					pair_score = self.map__ontology_domain__term_pair__score[ontology_domain][(t,q)]

				elif (q,t) in self.map__ontology_domain__term_pair__score[ontology_domain]:
					pair_score = self.map__ontology_domain__term_pair__score[ontology_domain][(q,t)]

				if pair_score > max_:
						max_ = pair_score

			score_2 += max_


		score = (score_1/n) + (score_2/m)


		return score

	def __iterate_over_protein_pairs__(self,E,ontology_domains):

		n_of_pairs = len(E)
		iteration = 0

		header = ["u", "v"]

		csv_writer = csv.writer(open(self.output_file_path,"w"),delimiter = "\t")
		for ontology_domain in ontology_domains:
			header.append(self.similarity_score_def + "_" + ontology_domain)

		csv_writer.writerow(header)
		
		for p, q in E:
			
			iteration += 1
			
			record = [p,q]
			
			for ontology_domain in ontology_domains:
				
				if p in self.map__ontology_domain__ensembl__id__ontologies[ontology_domain]:
					if q in self.map__ontology_domain__ensembl__id__ontologies[ontology_domain]:
							
						p_T = self.map__ontology_domain__ensembl__id__ontologies[ontology_domain][p]
						q_T = self.map__ontology_domain__ensembl__id__ontologies[ontology_domain][q]

						score = self.__get_best_match_average__(p_T,q_T, ontology_domain)
						record.append(score)
					else:
						record.append(0.0)
				else:
					record.append(0.0)

			csv_writer.writerow(record)


	def run(self,):
		print(self.output_file_path)
		if os.path.exists(self.output_file_path):
			return
		print(self.output_file_path)
		G_ontology, V_ontology = self.__load_acyclic_network__()
		E,V = self.__load_V__()

		map__ontology_domain__ontologies, self.map__ontology_domain__ensembl__id__ontologies = self.__load_map__ontology_domain__ontologies_for_a_better_version__(V)

		
		############# Resnik Similarity #############################
		if self.similarity_score_def == "resnik_similarity":
			
			rs = ResnikSimilarity(
				map__ontology_domain__ontologies = map__ontology_domain__ontologies,
				G = G_ontology,
				output_directory_path = "../datasets/ontology/"
				)
			
			resnik_sim_file_paths,ontology_domains = rs.run()

		############# Resnik Similarity #############################


		self.map__ontology_domain__term_pair__score = self.__load_term_similarity_dictionary__(resnik_sim_file_paths,ontology_domains)	
		self.__iterate_over_protein_pairs__(E,ontology_domains)

