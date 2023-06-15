import os
import math
import csv


class ResnikSimilarity():
	
	def __init__(self,map__ontology_domain__ontologies,G,output_directory_path,):
		
		self.map__ontology_domain__ontologies = map__ontology_domain__ontologies

		self.G = G

		self.output_directory_path = output_directory_path

		if not os.path.exists(self.output_directory_path):
			os.mkdir(self.output_directory_path)



	
	def __get_sub_concept__(self,node,H):
		
		count = set()

		for source,target in H.in_edges(node):
			count.add(source)
			
			for item in self.__get_sub_concept__(source,H):
				count.add(item)

		return count

	def __get_ancestor_concepts__(self,node,H):
		
		count = set()

		for source,target in H.out_edges(node):
			count.add(target)
			
			for item in self.__get_ancestor_concepts__(target,H):
				count.add(item)

		return count


	def __create_term_probability_by_domain_Graph__(self,H,N):

		map__term__probability = {}

		iteration = 0
		total = N

		for node in H:
			iteration += 1
			
			freq = len(self.__get_sub_concept__(node,H)) + 1
			
			if freq > 1:
				map__term__probability[node] = - math.log10(freq/N)

		return map__term__probability

	def __create_map__term__ancestors__(self,H,N):

		iteration = 0
		total = N
		map__term__ancestors = {}
		for node in H:
			iteration += 1
			
			map__term__ancestors[node] = self.__get_ancestor_concepts__(node,H)

		return map__term__ancestors



	def __create_resnik_similarity__(self,):

		output_file_paths = []
		ontology_domains = []
		for domain, ontologies in self.map__ontology_domain__ontologies.items():

			self.output_similarity_file_path = self.output_directory_path + domain +"__resnik_similarity.txt"
			
			output_file_paths.append(self.output_similarity_file_path)
			ontology_domains.append(domain)
			
			if os.path.exists(self.output_similarity_file_path):
				continue
			
			csv_writer = csv.writer(open(self.output_similarity_file_path, "w"),delimiter = "\t")
			csv_writer.writerow(["term_ID", "term_ID", "score"])
	
			selected_edges = [(u,v) for u,v,e in self.G.edges(data=True) if e['cagegory'] == domain]
			H = self.G.edge_subgraph(selected_edges)

			N = len(H.nodes())
			V = list(H.nodes())

			map__term__probability = self.__create_term_probability_by_domain_Graph__(H,N)
			map__term__ancestors = self.__create_map__term__ancestors__(H,N)
			

			iteration = 0
			total = N
			
			for i in range(N):
				iteration += 1
				
				for j in range(i+1,N):
						
					term_i = V[i]
					term_j = V[j]

					if term_i in ontologies and term_j in ontologies:

						ancestor_i = map__term__ancestors[term_i]
						ancestor_j = map__term__ancestors[term_j]
						
						common_ancestors = ancestor_i.intersection(ancestor_j)
						max_ = 0
						
						for common_ancestor in common_ancestors:
							probability = map__term__probability[common_ancestor]

							if probability > max_:
								max_ = probability

						if max_ > 0.0:
							csv_writer.writerow([term_i, term_j, max_])

		
		return output_file_paths,ontology_domains
	
	def run(self,):
		return self.__create_resnik_similarity__()
		

