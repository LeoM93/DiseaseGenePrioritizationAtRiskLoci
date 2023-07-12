import csv
import os
class PreprocessAccosications():
    
    def __init__(self,):
        self.GWAS_association_directory = '../experiments/GWAS_associations/association/'
        self.GWAS_associations_files = [self.GWAS_association_directory  + file for file in os.listdir(self.GWAS_association_directory) if file[0] != "."]
        self.GWAS_preprocess_directory = '../experiments/GWAS_associations/GWAS/'
        self.list_of_snps_preprocess_directory = '../experiments/GWAS_associations/list_of_snps/'
        self.input_vegas_2 = '../experiments/GWAS_associations/vegas/input/'

    def __load_and_preprocess_association__(self,file):
        csv_reader = csv.reader(open(file,'r'),delimiter = ",")
        file_name = file.split("/")[-1].split("_")[2].split("-")[0]
        table = [['chrom', 'position','SNP', 'pval']]
        map__chrm__list_of_snps = {}
        for index, row in enumerate(csv_reader):
            if index == 0:
                print(row)
                continue

            try:
                rsId = row[0].split("-")[0]
                if rsId[0] == 'r':
                    chr, position = row[-1].split(":")
                    exponent = float(row[1].split("x")[-1].replace(" ","").split("-")[-1])
                    p_val = 10**-exponent
                    table.append([chr, position, rsId,p_val])
                    if chr not in map__chrm__list_of_snps:
                        map__chrm__list_of_snps[chr] = []
                    map__chrm__list_of_snps[chr].append([rsId])       
            except:
                pass
        csv_writer = csv.writer(open(self.GWAS_preprocess_directory + file_name + ".tsv", "w"), delimiter = "\t")
        csv_writer.writerows(table)
        for k,v in map__chrm__list_of_snps.items():
            csv_writer = csv.writer(open(self.list_of_snps_preprocess_directory + file_name + "_" + str(k) +".tsv", "w"), delimiter = "\t")
            csv_writer.writerows(v) 
    
    def compute_snp_by_chromosome(self,):
        for file in self.GWAS_associations_files:
            self.__load_and_preprocess_association__(file)
    
    def compute_snp_filtered_by_linkage_dis(self,):
        self.list_of_snps_filtered_by_ld = [self.list_of_snps_preprocess_directory  + file for file in os.listdir(self.list_of_snps_preprocess_directory) if file[0] != "." and "_" not in file]

        for file in self.list_of_snps_filtered_by_ld:
            csv_reader = csv.reader(open(file),delimiter = "\t")
            set_of_snp = set()
            table = [['chrom', 'position','SNP', 'pval']]
            vegas_input = []
            for row in csv_reader:
                
                set_of_snp.add(row[0])

            GWAS_reader = csv.reader(open(self.GWAS_preprocess_directory + file.split("/")[-1]),delimiter = "\t")
            for index,row in enumerate(GWAS_reader):
                if index == 0:
                    continue
                if row[2] in set_of_snp:
                    table.append(row)
                    vegas_input.append([row[2], row[3]])
            csv_writer = csv.writer(open(self.GWAS_preprocess_directory + file.split("/")[-1].replace(".tsv","") + "_ld.tsv", "w"), delimiter = "\t")
            csv_writer.writerows(table)
            csv_writer = csv.writer(open(self.input_vegas_2 + file.split("/")[-1].replace(".tsv","") + "_ld.tsv", "w"), delimiter = "\t")
            csv_writer.writerows(vegas_input)

p = PreprocessAccosications()
p.compute_snp_by_chromosome()
p.compute_snp_filtered_by_linkage_dis()