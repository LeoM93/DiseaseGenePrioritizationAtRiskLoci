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
		ensembl_db_file_path,
		filter_ = ['dmGWAS','SIGMOID'],
		):
		self.locus_gene_file_path = locus_gene_file_path
		self.ensembl_db = ensembl_db_file_path
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

		self.__load_ensembl_db__()
		self.map__algorithm__solution = {dir_: self.__load_solution__(self.algorithms_file_path + dir_ + "/GCST007692.tsv") for dir_ in os.listdir(self.algorithms_file_path) if dir_[0] != "." and dir_ not in self.filter}
		

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


	def distribution_of_solution(self,):
		self.map__algorithm__locus__score = {}
		
		for algorithm, solution in self.map__algorithm__solution.items():
			self.map__algorithm__locus__score[algorithm] = {}
			for gene in solution:
				if gene in self.map__gene__locus:
					locus =  self.map__gene__locus[gene]

					if locus not in self.map__algorithm__locus__score[algorithm]:
						self.map__algorithm__locus__score[algorithm][locus] = 0
					self.map__algorithm__locus__score[algorithm][locus] += 1

		pd.DataFrame(self.map__algorithm__locus__score).fillna(0.0).to_csv(self.validation_dir_path + "distribution_of_solution.tsv",sep= "\t")




	def __load_solution__(self,file_path):

		csv_reader = csv.reader(open(file_path , "r"),delimiter = "\t")
		set_ = set()
		
		if 'RMM' in file_path:
			for row in csv_reader:
				set_.add(row[0])
		else:
			for row in csv_reader:
				
				if row[0]in self.map__ensembl_id__gene:
					set_.add(self.map__ensembl_id__gene[row[0]])
		
		print(set_)
		return set_

	def run_gseapy(self):

		dir_path = self.gsea_reactome_validation_dir_path + "standard/"
		
		if not os.path.exists(dir_path):
			os.makedirs(dir_path)
		
		for algorithm, set_ in self.map__algorithm__solution.items():


			result = gseapy.enrichr(gene_list=list(set_),
				gene_sets=[self.gene_set],
				organism='Human', # don't forget to set organism to the one you desired! e.g. Yeast
				outdir='test/enrichr_kegg',
				no_plot=True,
				cutoff= 0.05 # test dataset, use lower value from range(0,1)
		                )
			data_frame = result.results
				
			column_to_drop = ['Gene_set', 'Overlap', 'P-value',
						'Old P-value', 'Old Adjusted P-value', 'Odds Ratio', 'Combined Score',
			]
			
			data_frame = data_frame.drop(columns=column_to_drop)
				
			data_frame.to_csv(dir_path + algorithm + ".tsv", sep = "\t")



if __name__ == "__main__":


	gsea = GSEAValidation(
		locus_gene_file_path = "../../datasets/GWAS/unbiased_copd_snp_database.txt",
		disease_experiment_dir_path = "../../experiments/algorithm_comparison/",
		gene_set = "Reactome_2016", #"KEGG_2021_Human"
		ensembl_db_file_path = "../../datasets/curated_db/mart_export.txt",
	)

	gsea.run_gseapy()
