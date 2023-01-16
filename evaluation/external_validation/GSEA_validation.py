import gseapy
import os
import pandas as pd
import csv
from random import sample
import time

class GSEAValidation():
	def __init__(self,
		locus_gene_file_path,
		disease_experiment_dir_path,
		gene_set,
		filter_ = [],
		):
		self.locus_gene_file_path = locus_gene_file_path

		self.validation_dir_path = disease_experiment_dir_path + "validation/"
		self.disease_experiment_dir_path = disease_experiment_dir_path
		self.algorithms_file_path = disease_experiment_dir_path + "algorithms/"

		if not os.path.exists(self.validation_dir_path):
			os.makedirs(self.validation_dir_path)

		self.gene_set = gene_set
		if gene_set == "KEGG_2021_Human":
			self.gsea_reactome_validation_dir_path = self.validation_dir_path + "gsea_kegg/"
		else:
			self.gsea_reactome_validation_dir_path = self.validation_dir_path + "gsea_reactome/"
		
		if not os.path.exists(self.gsea_reactome_validation_dir_path):
			os.makedirs(self.gsea_reactome_validation_dir_path)

		self.filter = filter_

		self.map__algorithm__solution = {dir_: self.__load_solution__(self.algorithms_file_path + dir_ + "/solution.txt") for dir_ in os.listdir(self.algorithms_file_path) if dir_[0] != "." and dir_ not in self.filter}
		print(self.map__algorithm__solution)
		self.__load_locus_genes__()

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

	def __load_solution__(self,file_path):

		csv_reader = csv.reader(open(file_path , "r"),delimiter = "\t")
		set_ = set()
		
		for row in csv_reader:
			set_.add(row[0])

		return set_

	def run_gseapy(self):

		dir_path = self.gsea_reactome_validation_dir_path + "standard/"
		
		if not os.path.exists(dir_path):
			os.makedirs(dir_path)
		
		for algorithm, set_ in self.map__algorithm__solution.items():


			result = gseapy.enrichr(gene_list=list(set_),
				gene_sets=[self.gene_set],
				organism='Human', # don't forget to set organism to the one you desired! e.g. Yeast
				description='test_name',
				outdir='test/enrichr_kegg',
				# no_plot=True,
				cutoff= 0.05 # test dataset, use lower value from range(0,1)
		                )
			data_frame = result.results
				
			column_to_drop = ['Gene_set', 'Overlap', 'P-value',
						'Old P-value', 'Old Adjusted P-value', 'Odds Ratio', 'Combined Score',
			]
			
			data_frame = data_frame.drop(columns=column_to_drop)
				
			data_frame.to_csv(dir_path + algorithm + ".tsv", sep = "\t")



	def __compute_random_solution__(self,map__locus__solution):

		solution = set()
		
		for locus, set_ in map__locus__solution.items():

			random_gene = sample(set_,1)[0]

			solution.add(random_gene)

		return solution





	def run_gseapy_on_sample(self, trial = 10): #Reactome_2016

		names = gseapy.get_library_name()

		# show top 20 entries.
		for i in names:
			print(i)

		for algorithm, set_ in self.map__algorithm__solution.items():

			dir_path = self.gsea_reactome_validation_dir_path + algorithm +"/"

			if not os.path.exists(dir_path):
				os.makedirs(dir_path)


			map__locus__solution = {}
			
			for gene in set_:
				
				if gene not in self.map__gene__locus:
					continue
				
				locus = self.map__gene__locus[gene]

				if locus not in map__locus__solution:
					map__locus__solution[locus] = set()

				map__locus__solution[locus].add(gene)


			for i in range(trial):
				
				solution = self.__compute_random_solution__(map__locus__solution)

				result = gseapy.enrichr(gene_list=list(solution),
				
				gene_sets=[self.gene_set],
				organism='Human', # don't forget to set organism to the one you desired! e.g. Yeast
				description='test_name',
				outdir='test/enrichr_kegg',
				# no_plot=True,
				cutoff= 0.05 # test dataset, use lower value from range(0,1)
		                )
				data_frame = result.results
					
				column_to_drop = ['Gene_set', 'Overlap', 'P-value',
							'Old P-value', 'Old Adjusted P-value', 'Odds Ratio', 'Combined Score',
				]
				
				data_frame = data_frame.drop(columns=column_to_drop)
				data_frame.to_csv(dir_path + str(i) + ".tsv", sep = "\t")

				time.sleep(1)


if __name__ == "__main__":


	gsea = GSEAValidation(
		locus_gene_file_path = "../../datasets/GWAS/unbiased_copd_snp_database.txt",
		disease_experiment_dir_path = "../../experiments/algorithm_comparison/",
		filter_ = [],
		gene_set = "Reactome_2016" #KEGG_2021_Human
	)

	gsea.run_gseapy()
	gsea.run_gseapy_on_sample()