import csv
from scipy.stats import fisher_exact


import numpy as np
import random

def _ecdf(x):
	"""No frills empirical cdf used in fdrcorrection."""
	nobs = len(x)
	return np.arange(1, nobs + 1) / float(nobs)


def fdr_correction(pvals, alpha=0.05, method='indep'):

	"""P-value correction with False Discovery Rate (FDR).
	Correction for multiple comparison using FDR.
	This covers Benjamini/Hochberg for independent or positively correlated and
	Benjamini/Yekutieli for general or negatively correlated tests.
	Parameters
	----------
	pvals : array_like
		set of p-values of the individual tests.
	alpha : float
		error rate
	method : 'indep' | 'negcorr'
		If 'indep' it implements Benjamini/Hochberg for independent or if
		'negcorr' it corresponds to Benjamini/Yekutieli.
	Returns
	-------
	reject : array, bool
		True if a hypothesis is rejected, False if not
	pval_corrected : array
		pvalues adjusted for multiple hypothesis testing to limit FDR
	"""
	pvals = np.asarray(pvals)
	shape_init = pvals.shape
	pvals = pvals.ravel()

	pvals_sortind = np.argsort(pvals)
	pvals_sorted = pvals[pvals_sortind]
	sortrevind = pvals_sortind.argsort()

	if method in ['i', 'indep', 'p', 'poscorr']:
		ecdffactor = _ecdf(pvals_sorted)
	elif method in ['n', 'negcorr']:
		cm = np.sum(1. / np.arange(1, len(pvals_sorted) + 1))
		ecdffactor = _ecdf(pvals_sorted) / cm
	else:
		raise ValueError("Method should be 'indep' and 'negcorr'")

	reject = pvals_sorted < (ecdffactor * alpha)
	if reject.any():
		rejectmax = max(np.nonzero(reject)[0])
	else:
		rejectmax = 0
	reject[:rejectmax] = True

	pvals_corrected_raw = pvals_sorted / ecdffactor
	pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
	pvals_corrected[pvals_corrected > 1.0] = 1.0
	pvals_corrected = pvals_corrected[sortrevind].reshape(shape_init)
	reject = reject[sortrevind].reshape(shape_init)
	return reject, pvals_corrected


class EnrichmentTool():
	def __init__(self, total_universe_flag = False):

		self.ontology_file_path = "/Users/leonardomartini/Documents/network_medicine/Data/db/annotated_db/GO_10_02_2021.txt"
		self.universe_file_path = "/Users/leonardomartini/Documents/network_medicine/Data/db/gwas/COPD/unbiased_copd_snp_database.txt"
		
		self.go_annotation_file_path = "/Users/leonardomartini/Documents/network_medicine/Data/db/annotated_db/ontology/go-basic.txt"
		self.reactome_file_path = "/Users/leonardomartini/Documents/network_medicine/Data/db/annotated_db/GO_KEGG_Reactome_total.txt"
		self.reactom_description_file_path = "/Users/leonardomartini/Documents/network_medicine/Data/db/annotation_description/ReactomePathways.txt"
		self.solution = set()

		self.total_universe_flag = total_universe_flag
		
		self.__build__()


	def __load_universe__(self,):

		self.V = set()
		self.map__locus__genes ={}

		if self.total_universe_flag is True:
			return 

		with open(self.universe_file_path,'r') as fp:
			csv_reader = csv.reader(fp, delimiter = "\t")
			
			for index,row in enumerate(csv_reader):

				if row[1] not in self.map__locus__genes:
					self.map__locus__genes[row[1]] = set()

				self.map__locus__genes[row[1]].add(row[0])

				gene_name = row[0]

				self.V.add(gene_name)

	def __create_random_solution__(self,):

		random_solution = set()

		for locus,genes in self.map__locus__genes.items():
			random_solution.add(random.sample(genes,1)[0])

		return random_solution

	def __load_solution__(self,):
		self.solution = set()

		if type(self.solution_file_path) is set:
			self.solution = self.solution_file_path

		else:
			with open(self.solution_file_path,'r') as fp:
				csv_reader = csv.reader(fp, delimiter = "\t")
				for index,row in enumerate(csv_reader):

					gene_name = row[0]

					self.solution.add(gene_name)


	def __load_reactome_description__(self,):

		with open(self.reactom_description_file_path,"r") as fp:

			csv_reader = csv.reader(fp, delimiter = "\t")
			for index, row in enumerate(csv_reader):
				term_id = row[1]
				description = row[0]

				self.map__term_id__description[term_id] = description



	
	def __load_map__ontology_domain__ontologies__(self, ES = "IEA",current_domain = "biological_process"):
		
		self.map__ensembl__id__ontologies = {}
		self.map__ontology_id__ensembl_ids = {}
		self.map__term_id__description = {}
		
		with open(self.go_annotation_file_path,"r") as fp:
			
			lines = fp.readlines()
			term_id = ""
			description = ""

			for line in lines:

				if line[:3] == "id:":
					term_id = line.split(": ")[-1].replace("\n","")
					
					if "GO:" in term_id:
						self.map__term_id__description[term_id] = ""
					else:
						term_id = ""


				if line[:5] == "name:":
					description = line.split(": ")[-1].replace("\n","")

					if term_id != "":
						self.map__term_id__description[term_id] = description

		
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

				if gene_name in self.V:

					if term_id not in self.map__ontology_id__ensembl_ids:
						self.map__ontology_id__ensembl_ids[term_id] = set()

					self.map__ontology_id__ensembl_ids[term_id].add(gene_name)

					if gene_name not in self.map__ensembl__id__ontologies:
						self.map__ensembl__id__ontologies[gene_name] = set()

					self.map__ensembl__id__ontologies[gene_name].add(term_id)


	
	def __load_reactome__(self,):
		
		self.map__ensembl__id__reactome = {}
		self.map__reactome_id__ensemble = {}

		with open(self.reactome_file_path, "r") as fp:
			csv_reader = csv.reader(fp,delimiter = "\t")
			
			for index, row in enumerate(csv_reader):
				
				if index == 0:
					continue

				gene_name = row[0]
				term_id = row[1]
				database = row[2]

				if database != "Reactome":
					continue

				if len(self.V) == 0:
					
					if gene_name not in self.map__ensembl__id__reactome:
						self.map__ensembl__id__reactome[gene_name] = set()

					if term_id not in self.map__reactome_id__ensemble:
						self.map__reactome_id__ensemble[term_id] = set()

					self.map__ensembl__id__reactome[gene_name].add(term_id)
					self.map__reactome_id__ensemble[term_id].add(gene_name)


				elif gene_name in self.V:

					if gene_name not in self.map__ensembl__id__reactome:
						self.map__ensembl__id__reactome[gene_name] = set()

					if term_id not in self.map__reactome_id__ensemble:
						self.map__reactome_id__ensemble[term_id] = set()

					self.map__ensembl__id__reactome[gene_name].add(term_id)
					self.map__reactome_id__ensemble[term_id].add(gene_name)

		if len(self.V) == 0:
			self.V = set(self.map__ensembl__id__reactome.keys())


	def __reactome_enrichment__(self,print_flag):
		term_id_list = []
		p_vals = []
		sub_solution_list = []

		enrichment = set()


		for term_id in self.reactome_ids:
			term_id_list.append(term_id)

			ensebl_ids_involved_in_term_id = self.map__reactome_id__ensemble[term_id]

			ensebl_ids_involved_in_term_id_in_solution = ensebl_ids_involved_in_term_id.intersection(self.solution)
			ensebl_ids_not_involved_in_term_id_in_solution = self.solution.difference(ensebl_ids_involved_in_term_id)

			universe_minus_solution = self.V.difference(self.solution)
			universe_minus_solution_involved_in_term_id = universe_minus_solution.intersection(ensebl_ids_involved_in_term_id)
			universe_minus_solution_not_involved_in_term_id = universe_minus_solution.difference(universe_minus_solution_involved_in_term_id)

			a = len(ensebl_ids_involved_in_term_id_in_solution)
			b = len(ensebl_ids_not_involved_in_term_id_in_solution)
			c = len(universe_minus_solution_involved_in_term_id)
			d = len(universe_minus_solution_not_involved_in_term_id)
			
			table = [[a,b],[c,d]]

			odd,p_value = fisher_exact(table, alternative='two-sided')

			if term_id == "R-HSA-5099900":
				print(p_value)
			p_vals.append(p_value)
			sub_solution_list.append(ensebl_ids_involved_in_term_id_in_solution)


		result, scores = fdr_correction(p_vals)

		for flag, term_id,score,sub_set in zip(result,term_id_list,scores,sub_solution_list):
			if flag:
				enrichment.add(term_id)

				if print_flag:
					
					if term_id in self.map__term_id__description:
						print(term_id + "\t" + self.map__term_id__description[term_id] +"\t" + str(score) + "\t" + str(sub_set))

					else:
						print(term_id +"\t" + str(score) + "\t" + str(sub_set))
					

		return enrichment


	def __get_solution_go_and_reactome_term_ids__(self,):

		self.term_ids = set()
		self.reactome_ids = set()

		for node in self.solution:
			
			if node in self.map__ensembl__id__ontologies:
				
				for term_id in self.map__ensembl__id__ontologies[node]:
					self.term_ids.add(term_id)


			if node in self.map__ensembl__id__reactome:
				
				for term_id in self.map__ensembl__id__reactome[node]:
					self.reactome_ids.add(term_id)




	def __go_enrichment__(self, print_flag):

		term_id_list = []
		p_vals = []
		sub_solution_list = []
		enrichment = set()

		for term_id in self.term_ids:
			term_id_list.append(term_id)

			ensebl_ids_involved_in_term_id = self.map__ontology_id__ensembl_ids[term_id]
			
			ensebl_ids_involved_in_term_id_in_solution = ensebl_ids_involved_in_term_id.intersection(self.solution)
			ensebl_ids_not_involved_in_term_id_in_solution = self.solution.difference(ensebl_ids_involved_in_term_id)

			universe_minus_solution = self.V.difference(self.solution)
			universe_minus_solution_involved_in_term_id = universe_minus_solution.intersection(ensebl_ids_involved_in_term_id)
			universe_minus_solution_not_involved_in_term_id = universe_minus_solution.difference(universe_minus_solution_involved_in_term_id)

			a = len(ensebl_ids_involved_in_term_id_in_solution)
			b = len(ensebl_ids_not_involved_in_term_id_in_solution)
			c = len(universe_minus_solution_involved_in_term_id)
			d = len(universe_minus_solution_not_involved_in_term_id)
			
			table = [[a,b],[c,d]]
			odd,p_value = fisher_exact(table, alternative='two-sided')

			p_vals.append(p_value)
			sub_solution_list.append(ensebl_ids_involved_in_term_id_in_solution)


		result, scores = fdr_correction(p_vals)
		for flag, term_id,score,sub_set in zip(result,term_id_list,scores,sub_solution_list):
			if flag:

				enrichment.add(term_id)
				if print_flag:
					
					if term_id in self.map__term_id__description:
						print(term_id + "\t" + self.map__term_id__description[term_id] +"\t" + str(score) + "\t" + str(sub_set))

					else:
						print(term_id +"\t" + str(score) + "\t" + str(sub_set))
					

		return enrichment
			

	def __build__(self,):

		self.__load_universe__()	
		self.__load_reactome__()
		self.__load_map__ontology_domain__ontologies__()
		self.__load_reactome_description__()

	def run_on_reactome(self,solution_file_path, print_flag):

		self.solution_file_path = solution_file_path
		self.__load_solution__()
		self.__get_solution_go_and_reactome_term_ids__()

		return self.__reactome_enrichment__(print_flag)

	def run_on_go(self,solution_file_path,print_flag):

		self.solution_file_path = solution_file_path
		self.__load_solution__()
		self.__get_solution_go_and_reactome_term_ids__()
		return self.__go_enrichment__(print_flag)


