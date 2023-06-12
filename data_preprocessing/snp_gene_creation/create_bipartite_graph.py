import csv

class SNPGeneGraph():
	
	def __init__(self,
		rs_id_position_file_path,
		gene_overlap_2mb_windows_file_path,
		output_file_path):
		self.rs_id_position_file_path = rs_id_position_file_path
		self.gene_overlap_2mb_windows_file_path = gene_overlap_2mb_windows_file_path
		self.output_file_path = output_file_path

	def __load_rs_id__position__(self,):
		csv_reader = csv.reader(open(self.rs_id_position_file_path, "r"),delimiter = "\t")
		self.map__rs_id__position = {}
		for index, row in enumerate(csv_reader):
			if index == 0:
				continue
			rs_id = row[2]
			chr_ = row[0]
			position = row[1]
			p_val = float(row[3])

			self.map__rs_id__position[rs_id] = (chr_,int(position),p_val)

	def __load_genes_overlapping_windows__(self,):
		csv_reader = csv.reader(open(self.gene_overlap_2mb_windows_file_path, "r"),delimiter = "\t")
		self.map__gene__positions = {}
		for index, row in enumerate(csv_reader):
			if index == 0:
				continue
			chr_ = row[0].replace("chr","")
			start = int(row[1])
			end = int(row[2])
			gene_name = row[3]

			if gene_name not in self.map__gene__positions:
				self.map__gene__positions[gene_name] = set()

			self.map__gene__positions[gene_name].add((chr_, start,end))


	def run(self, w = 1e06):
		self.__load_rs_id__position__()
		self.__load_genes_overlapping_windows__()
		
		gene_snp_graph = []
		
		self.map__gene__closest_loci = {}
		self.map__snp__closest_genes = {}
		
		map__gene__p_val = {}
		
		for gene, positions in self.map__gene__positions.items():
			self.map__gene__closest_loci[gene] = []
			
			for position in positions:
				chr_,start,end = position

				for rs_id, snp_position in self.map__rs_id__position.items():
					snp_chr, snp_locus, p_val = snp_position

					if snp_chr != chr_:
						continue

					if rs_id not in self.map__snp__closest_genes:
						self.map__snp__closest_genes[rs_id] = []
										
					distance = min(abs(snp_locus - start),abs(snp_locus - end))
					
					if distance < w:
						self.map__gene__closest_loci[gene].append((rs_id, distance))
						self.map__snp__closest_genes[rs_id].append((gene,distance,p_val))
		
		for gene, closest_snps in self.map__gene__closest_loci.items():
			closest_snps.sort(key = lambda x:x[1], reverse = False)
			
			if len(closest_snps) == 0:
				continue

			closest_snp = closest_snps[0]
			gene_snp_graph.append([gene,closest_snp[0]])

		
		map__gene_pool__uniform_pval = {}
		for rs_id, closest_genes in self.map__snp__closest_genes.items():
			closest_genes.sort(key = lambda x:x[1], reverse = False)

			if len(closest_genes) == 0:
				continue
			
			closest_gene = closest_genes[0]
			map__gene__p_val[closest_gene[0]] = p_val

			for gene in closest_genes:
				map__gene_pool__uniform_pval[gene[0]] = 0.05
		
		
		if self.output_file_path != "":
			csv_writer = csv.writer(open(self.output_file_path, "w"),delimiter = "\t")
			csv_writer.writerows(gene_snp_graph)
		
		return map__gene__p_val,map__gene_pool__uniform_pval


