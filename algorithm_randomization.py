import csv
import random
import os
import subprocess
import seaborn as sns
sns.set_theme(style="whitegrid")
import matplotlib.pyplot as plt
import matplotlib as mpl

class AlgorithmRandomization():
    def __init__(
        self, 
        rmm_gwas_output_path,
        experiment_dir_path,
        GWAS_dir_path,
        rmm_input_file_path,):


        self.onco_kb = "./datasets/curated_db/onco_kb.tsv"
        self.drug_table_file_path = "./datasets/curated_db/drug_repurposing_hub.tsv"
        self.cosmic_db = "./datasets/curated_db/cosmic_census.tsv"
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
        #self.__run_rmm_gwas_on_random_seeds__()
        self.__create_distribution__()
    
    def __load_cosmic_db__(self,):

        csv_reader = csv.reader(open(self.cosmic_db,"r"),delimiter = "\t")
        ground_truth = set()
        for index, row in enumerate(csv_reader):
            if index == 0:
                continue
            
            somatic_type = row[9]
            germline_type = row[10]
                
            if "breast" in somatic_type or "breast" in germline_type:
                ground_truth.add(row[0])
        return ground_truth
    
    def __load_solution__(self,file):
        csv_reader = csv.reader(open(file,"r"),delimiter = "\t")
        solution = set()
        for index, row in enumerate(csv_reader):
            solution.add(row[0])
        
        return solution
    
    def __load_drug_repurposing_hub__(self,):

        self.map__target__drugs = {}
        self.map__clinical_phase_disease_area_indication__targets = {}

        with open(self.drug_table_file_path, "r") as fp:
            csv_reader = csv.reader(fp,delimiter = "\t")

            for index,row in enumerate(csv_reader):

                if index == 0:
                    continue

                targets = row[3].split("|")
                disease_area = row[4]
                clinical_phase = row[1]
                mou = row[2]
                drug_name = row[0]
                indication = row[5]

                if (clinical_phase,disease_area,indication) not in self.map__clinical_phase_disease_area_indication__targets:
                    self.map__clinical_phase_disease_area_indication__targets[(clinical_phase,disease_area,indication)] = set()


                for target in targets:

                    if target in self.map__target__drugs:
                        self.map__target__drugs[target].append((drug_name, disease_area, mou,clinical_phase,indication))
                    else:
                        self.map__target__drugs[target] = [(drug_name, disease_area, mou,clinical_phase,indication)]

                    self.map__clinical_phase_disease_area_indication__targets[(clinical_phase,disease_area,indication)].add(target)



    def __compute_drug_targets__(self, query_disease_area = None, query_indications = None): 	
        target_nodes = set()
		

        for node in self.map__target__drugs:
            if node == "":
                continue

            target_drugs = self.map__target__drugs[node]

            for drug_name, disease_area, mou,clinical_phase,indication in target_drugs:


                if query_disease_area is not None:
	
                    if query_disease_area in disease_area:


                        if query_indications is not None:
                            query_indications_vector = query_indications.split("|")
                            for query_indication in query_indications_vector:
                                if query_indication in indication:
	                                target_nodes.add(node)
                        else:
                            target_nodes.add(node)

                else:
                    target_nodes.add(node)

	
        return target_nodes
    
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


    def __create_distribution__(self,dir_ = 'internal_validation'):
        target_nodes = self.__load_onco_bk__()
        self.__load_drug_repurposing_hub__()
        self.__load_cosmic_db__()

        drug_targets = self.__compute_drug_targets__("oncology", "breast cancer")
        drug_targets =  self. __load_cosmic_db__()
        onco_gene_distribution = []
        drug_targets_distribution = []

        for file in os.listdir(self.rmm_output_dir_path):
            if '.tsv' in file:
                solution = self.__load_solution__(self.rmm_output_dir_path + file)
                
                precision = len(target_nodes.intersection(solution))/len(solution)
                onco_gene_distribution.append(precision)

                precision = len(drug_targets.intersection(solution))/len(solution)
                drug_targets_distribution.append(precision)
        
        solution = self.__load_solution__(self.rmm_gwas_output_path )
        onco_precision = len(target_nodes.intersection(solution))/len(solution)

        solution = self.__load_solution__(self.rmm_gwas_output_path )
        drug_precision = len(drug_targets.intersection(solution))/len(solution)
        countr = 0
        onco_contr = 0
        for x in drug_targets_distribution:
            if x > drug_precision:
               countr += 1
        p_val = countr/len(drug_targets_distribution)
        f,axes = plt.subplots(2,1,figsize = (5,10))

        for x in onco_gene_distribution:
            if x > onco_precision:
               onco_contr += 1
        p_val_onco = onco_contr/len(onco_gene_distribution)
        f,axes = plt.subplots(2,1,figsize = (5,10))
        
        sns.histplot(onco_gene_distribution,kde=True, ax=axes[0],stat='probability')
        axes[0].set_title('Dataset: OncoKB')
        axes[0].axvline(x=onco_precision, color='red',linestyle='--', label=f'RMM-GWAS recall')
        axes[0].annotate('$p_{val}:$' + str(p_val_onco), xy=(onco_precision, 0.0), xytext=(onco_precision+ 0.001, 0.025 ),
                        arrowprops=dict(arrowstyle="->",color='blue'))

        sns.histplot(drug_targets_distribution,kde=True, bins=10,ax=axes[1],stat='probability')
        axes[1].axvline(x=drug_precision, color='red', linestyle='--', label=f'',)
        
        axes[1].annotate('$p_{val}:$' + str(p_val), xy=(drug_precision, 0.0), xytext=(drug_precision+ 0.001, 0.05 ),
                        arrowprops=dict(arrowstyle="->",color='blue'))
        axes[1].set_title('Dataset: COSMIC')
        axes[0].legend()

        f.savefig('./imgs/' + dir_ + '.pdf',bbox_inches='tight')

        
ar = AlgorithmRandomization(
    rmm_gwas_output_path = "/Users/leonardomartini/Documents/network_medicine/DiseaseGenePrioritizationAtRiskLoci/experiments/algorithm_comparison/algorithms/RMM-GWAS/GCST004988.tsv",
    experiment_dir_path = "/Users/leonardomartini/Documents/network_medicine/DiseaseGenePrioritizationAtRiskLoci/experiments/randomization/",
    GWAS_dir_path = "/Users/leonardomartini/Documents/network_medicine/DiseaseGenePrioritizationAtRiskLoci/experiments/GWAS/", 
    rmm_input_file_path = "/Users/leonardomartini/Documents/network_medicine/DiseaseGenePrioritizationAtRiskLoci/experiments/GWAS/seed_RMM-GWAS/GCST004988.tsv"
)