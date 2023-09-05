import csv
import os
import numpy as np
import pandas as pd
import random
import json

class CosmicOncoKb():

    def __init__(self, 
        
        disease_experiment_dir_path,
        ensembl_db_file_path,
        cosmic_db,
        oncokb_db,
        config_file,
        filter_
        ):
        self.config_file = config_file
        self.map__trait__disease_info = json.load(open(self.config_file))


        self.disease_experiment_dir_path = disease_experiment_dir_path
        self.algorithms_file_path = disease_experiment_dir_path + "algorithms/"
        self.validation_dir_path = disease_experiment_dir_path + "validation/"
        
        self.seed_dir_path =  "../../experiments/GWAS_studies/seed/"
        self.seed_dir_rmm_gwas_path = "../../experiments/GWAS_studies/seed_RMM-GWAS/"
        self.cosmic_db = cosmic_db
        
        if not os.path.exists(self.validation_dir_path):
            os.makedirs(self.validation_dir_path)

        self.filter_ = filter_
        self.ensembl_db = ensembl_db_file_path
        
        self.__load_ensembl_db__()
        self.map__algorithm__solutions = {dir_: self.__load_solutions__(self.algorithms_file_path + dir_ +"/",algorithm_name = dir_) for dir_ in os.listdir(self.algorithms_file_path) if dir_[0] != "." and dir_ not in self.filter_}
        self.__load_cosmic_db__()
        self.__load_disease_seed__()
    
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
    
    
    def __load_solutions__(self,directory_path, algorithm_name):
        
        file_paths = [directory_path + file for file in os.listdir(directory_path) if file[0] != "." ]
        map__disease_name__disease_module = {}
        
        for file_path in file_paths:
            csv_reader = csv.reader(open(file_path , "r"),delimiter = "\t")
            set_ = set()
        
            for row in csv_reader:
                if algorithm_name != 'RMM-GWAS' and algorithm_name != 'RMM-GWAS-OLD' and algorithm_name != 'DSG':
                    if row[0] in self.map__ensembl_id__gene:
                        set_.add(self.map__ensembl_id__gene[row[0]])
                else:
                    set_.add(row[0])
            map__disease_name__disease_module[file_path.split("/")[-1].replace(".tsv","")] = set_

        return map__disease_name__disease_module
    
    def __load_disease_seed__(self):

        file_paths = [file for file in os.listdir(self.seed_dir_path) if file[0] != "."]
        map__disease__seed = {}
        for file in file_paths:
            csv_reader = csv.reader(open(self.seed_dir_path +file, "r"),delimiter = "\t")
            set_ = set()
            for row in csv_reader:
                if row[0] in self.map__ensembl_id__gene:
                    set_.add(self.map__ensembl_id__gene[row[0]])

            map__disease__seed[file.split("/")[-1].replace(".tsv","")] = set_
        return map__disease__seed
    
    def __load_disease_seed_for_RMM_GWAS__(self):

        file_paths = [file for file in os.listdir(self.seed_dir_rmm_gwas_path) if file[0] != "."]
        map__disease__seed = {}
        map__disease__gene__locus = {}
        for file in file_paths:
            map__disease__gene__locus[file.split("/")[-1].replace(".tsv","")] = {}
            csv_reader = csv.reader(open(self.seed_dir_rmm_gwas_path +file, "r"),delimiter = "\t")
            set_ = set()
            for row in csv_reader:
                set_.add(row[0])
                map__disease__gene__locus[file.split("/")[-1].replace(".tsv","")][row[0]] = row[1]

            map__disease__seed[file.split("/")[-1].replace(".tsv","")] = set_
        return map__disease__seed,map__disease__gene__locus
    
    def __load_onco_bk__(self,):
        csv_reader = csv.reader(open(self.cosmic_db,"r"),delimiter = "\t")
        target_nodes = set()
        for index, row in enumerate(csv_reader):
            
            target_nodes.add(row[0])
        
        return target_nodes

    
    def __load_cosmic_db__(self,):

        csv_reader = csv.reader(open(self.cosmic_db,"r"),delimiter = "\t")
        self.map__phenotype__groundtruth = {}
        for index, row in enumerate(csv_reader):
            if index == 0:
                continue
            
            somatic_type = row[9]
            germline_type = row[10]

            for value in self.map__trait__disease_info.values():
                if value not in self.map__phenotype__groundtruth:
                    self.map__phenotype__groundtruth[value] = set()
                
                if value in somatic_type or value in germline_type:
                    self.map__phenotype__groundtruth[value].add(row[0])

    def run(self, trial = 100, validation_dir = "cosmic/", cosmic = True):
        cosmic_validation_dir =self.validation_dir_path + validation_dir
        
        if not os.path.exists(cosmic_validation_dir):
            os.makedirs(cosmic_validation_dir)
        
        map__disease__seed = self.__load_disease_seed__()
        map__disease__seed_RMM_GWAS,map__disease__gene__locus = self.__load_disease_seed_for_RMM_GWAS__()
        precision_table = []
        drug_indication = []
        algorithm_pval = []
        random_distribution = []
        for algorithm, solutions in self.map__algorithm__solutions.items():
            for gwas, solution in solutions.items():
                
                if gwas in self.map__trait__disease_info:
                    if cosmic:
                        target_nodes = self.map__phenotype__groundtruth[self.map__trait__disease_info[gwas]]
                        target_nodes = target_nodes.intersection(map__disease__seed_RMM_GWAS[gwas.split("__")[0]])
                    else:
                        target_nodes = self.__load_onco_bk__()
                        target_nodes = target_nodes.intersection(map__disease__seed_RMM_GWAS[gwas.split("__")[0]])
                    
                    if len(target_nodes) == 0:
                        continue
                    
                    target_nodes_locus = len(set([map__disease__gene__locus[gwas.split("__")[0]][item] for item in target_nodes if item in map__disease__gene__locus[gwas.split("__")[0]]]))

                    phi =  len(solution.intersection(target_nodes))
                    phi_locus = len(set([map__disease__gene__locus[gwas.split("__")[0]][item] for item in solution.intersection(target_nodes) if item in map__disease__gene__locus[gwas.split("__")[0]]]))
                    ld = "Not filtered"    
                    precision_table.append([algorithm,gwas.split("_")[0],phi/len(target_nodes),ld, gwas.split("__")[-1]])
                    

                    counter = 0
                    for i in range(trial):
                        print(gwas,algorithm)
                        
                        if algorithm == 'RMM-GWAS':
                            random_solution = set(random.sample(map__disease__seed_RMM_GWAS[gwas],len(solution)))
                        else:
                            random_solution = set(random.sample(map__disease__seed[gwas],len(solution.intersection(map__disease__seed[gwas]))))

                        
                        phi_random = len(random_solution.intersection(target_nodes))

                        random_distribution.append([algorithm, gwas.split("_")[0],phi_random/len(random_solution),ld,gwas.split("__")[-1]])

                        if phi <= phi_random:
                            counter += 1
                    if counter == 0:
                        counter_str = "p_{val} < 10^{-6}"
                    else:
                        counter_str = "p_{val} ~ " + str(counter/trial) 


                    algorithm_pval.append([algorithm, gwas.split("_")[0],counter_str,ld,gwas.split("__")[-1]])
        
        pd.DataFrame(precision_table, columns = ["Algorithm","GWAS", "Recall", "filter","input type"]).to_csv(cosmic_validation_dir + "metrics.tsv", sep = "\t")
        pd.DataFrame(random_distribution,columns = ["Algorithm","GWAS","Recall", "filter","input type"]).to_csv(cosmic_validation_dir + "metrics_random_distribution.tsv", sep = "\t")
        pd.DataFrame(algorithm_pval,columns = ["Algorithm","GWAS","p_val", "filter","input type"]).to_csv(cosmic_validation_dir + "metrics_p_val.tsv", sep = "\t")






co = CosmicOncoKb(
    disease_experiment_dir_path = "../../experiments/GWAS_study_exp/",
    cosmic_db = "../../datasets/curated_db/cosmic_census.tsv",
    oncokb_db = "../../datasets/curated_db/onco_kb.tsv",
    config_file = "config_files/cosmic.json",
    ensembl_db_file_path = "../../datasets/curated_db/mart_export.txt",
    filter_ = []
)

co.run()
co.run( validation_dir = "onco/", cosmic = False)