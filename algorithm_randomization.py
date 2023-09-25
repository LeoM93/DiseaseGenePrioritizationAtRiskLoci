import csv
import random
import os
import subprocess

class AlgorithmRandomization():
    def __init__(
        self, 
        rmm_gwas_output_path,
        experiment_dir_path,
        GWAS_dir_path,
        rmm_input_file_path,):


        self.onco_kb = "./datasets/curated_db/onco_kb.tsv"
        self.rmm_gwas_output_path = rmm_gwas_output_path
        self.GWAS_dir_path = GWAS_dir_path
        self.experiment_dir_path = experiment_dir_path
        self.rmm_input_file_path = rmm_input_file_path

        self.random_seed_dir_path =  self.GWAS_dir_path + "random_seeds/"
        self.input_network_dir_path_for_RMM_GWAS = self.GWAS_dir_path + "random_newtork_RMM-GWAS/"
        
        
        if not os.path.exists(self.random_seed_dir_path):
            os.makedirs(self.random_seed_dir_path)
        
        if not os.path.exists(self.experiment_dir_path):
            os.makedirs(self.experiment_dir_path)
        
        if not os.path.exists(self.input_network_dir_path_for_RMM_GWAS):
            os.makedirs(self.input_network_dir_path_for_RMM_GWAS)

        self.rmm_output_dir_path = self.experiment_dir_path + 'RMM-GWAS/'

        if not os.path.exists(self.rmm_output_dir_path):
            os.makedirs(self.rmm_output_dir_path)


        #self.__save_random_seeds__()
        self.__run_rmm_gwas_on_random_seeds__()
        self.__create_distribution_on_onco_kb__()
    
    def __load_solution__(self,file):
        csv_reader = csv.reader(open(file,"r"),delimiter = "\t")
        solution = set()
        for index, row in enumerate(csv_reader):
            solution.add(row[0])
        
        return solution
    def __load_onco_bk__(self,):
        csv_reader = csv.reader(open(self.onco_kb,"r"),delimiter = "\t")
        target_nodes = set()
        for index, row in enumerate(csv_reader):
            
            target_nodes.add(row[0])
        
        return target_nodes
    
    def __load_rmm_input__(self,):
        map__locus__genes = {}
        gene_pool = set()
        csv_reader = csv.reader(open(self.rmm_input_file_path),delimiter = "\t")
        
        for row in csv_reader:
            
            locus = row[1]
            gene = row[0]
            
            gene_pool.add(gene)

            if locus not in map__locus__genes:
                map__locus__genes[locus] = set()
            
            map__locus__genes[locus].add(gene)

        return map__locus__genes, gene_pool

    def __create_random_input__(self,output_file_path, map__locus__genes, gene_pool):
        
        current_gene_pool = gene_pool
        map__locus__random_genes = {}
        
        for k,v in map__locus__genes.items():
            size = len(v)

            sample = set(random.sample(gene_pool, size))
            map__locus__random_genes[k] = sample
            
            current_gene_pool = current_gene_pool.difference(sample)
        
        final_table = []
        
        for k,v in map__locus__random_genes.items():
            for item in v:
                final_table.append([item,k])


        csv_writer = csv.writer(open(output_file_path, "w"),delimiter = "\t")
        csv_writer.writerows(final_table)
    

    def __save_random_seeds__(self, threshold = 1000):
        map__locus__genes, gene_pool = self.__load_rmm_input__()

        for i in range(threshold):
            self.__create_random_input__(self.random_seed_dir_path + str(i) + ".tsv",map__locus__genes,gene_pool )

    def __run__RMM_GWAS__(self,
    file_path, 
    network_file_path = "/Users/leonardomartini/Documents/network_medicine/DiseaseGenePrioritizationAtRiskLoci/datasets/networks/co_regulated_network_1e-5_thresholded.txt"):
        algorithm_name = "RMM-GWAS"
        disease_name = file_path.split("/")[-1]

        preprocessing = 'cd data_preprocessing; python3 main.py'
        command = f'{preprocessing} {network_file_path} {self.input_network_dir_path_for_RMM_GWAS + disease_name} {file_path}'

        print(command)
        subprocess.call(command, shell=True,env = os.environ,stdout=subprocess.PIPE)

        run_rmm_gwas =  f'cd rmm_gwas; python3 run_relations_maximization.py {self.input_network_dir_path_for_RMM_GWAS + disease_name} {self.rmm_output_dir_path + disease_name}'
        print(run_rmm_gwas)
        subprocess.call(run_rmm_gwas, shell=True,env = os.environ,stdout=subprocess.PIPE)


    def __run_rmm_gwas_on_random_seeds__(self,):
        for file in os.listdir(self.random_seed_dir_path):
            if '.tsv' in file:
                print(self.random_seed_dir_path + file)
                self.__run__RMM_GWAS__(self.random_seed_dir_path + file)


    def __create_distribution_on_onco_kb__(self,):
        target_nodes = self.__load_onco_bk__()

        for file in os.listdir(self.rmm_output_dir_path):
            if '.tsv' in file:
                solution = self.__load_solution__(self.rmm_output_dir_path + file)
                precision = len(target_nodes.intersection(solution))/len(solution)
                print(precision)
        
        solution = self.__load_solution__(self.rmm_gwas_output_path )
        precision = len(target_nodes.intersection(solution))/len(solution)
        print(precision)

ar = AlgorithmRandomization(
    rmm_gwas_output_path = "/Users/leonardomartini/Documents/network_medicine/DiseaseGenePrioritizationAtRiskLoci/experiments/algorithm_comparison/algorithms/RMM-GWAS/GCST004988.tsv",
    experiment_dir_path = "/Users/leonardomartini/Documents/network_medicine/DiseaseGenePrioritizationAtRiskLoci/experiments/randomization/",
    GWAS_dir_path = "/Users/leonardomartini/Documents/network_medicine/DiseaseGenePrioritizationAtRiskLoci/experiments/GWAS/", 
    rmm_input_file_path = "/Users/leonardomartini/Documents/network_medicine/DiseaseGenePrioritizationAtRiskLoci/experiments/GWAS/seed_RMM-GWAS/GCST004988.tsv"
)