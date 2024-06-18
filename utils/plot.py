import os
import math
import csv
import pandas as pd
import seaborn as sns
sns.set_theme(style="whitegrid")
import string
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl



class Plot():
    def __init__(self, experiment_dir_path):
        
        self.experiment_dir_path = experiment_dir_path
        self.validation_path = experiment_dir_path + "validation/"
        self.reactome_relationship = "../datasets/curated_db/reactome_pathway_relationship.tsv"
        self.cluster_map = "../imgs/gsea.tsv"

    def __load_pvals__(self, dir_,file_name = "metrics_p_val.tsv"):
        return pd.read_csv(dir_ + file_name, sep = "\t", index_col = 0)
    
    def __load_metrics__(self, dir_, file_name = "metrics.tsv"):
        return pd.read_csv(dir_ + file_name, sep = "\t", index_col = 0)

    
    def __load_random_distribution__(self, dir_,file_name = "metrics_random_distribution.tsv"):
        return pd.read_csv(dir_ + file_name, sep = "\t", index_col = 0)
    
    
    def __load_reactome_graph__(self,):
        
        csv_reader = csv.reader(open(self.reactome_relationship),delimiter = "\t")
        self.G = nx.DiGraph()
        
        for row in csv_reader:
            self.G.add_edge(row[0], row[1])

        self.root_nodes = [node for node in self.G.nodes() if self.G.in_degree(node) == 0]

    def __find_root_path__(self, node):
        if node in self.root_nodes:
            return [node]
        for neighbor in self.G.predecessors(node):
            if neighbor not in self.visited:
                self.visited.add(neighbor)
                path = self.__find_root_path__( neighbor)
                if path:
                    return [node] + path
        return []
    
    def __load_gsea__(self,file_path,threshold_value = 0.05):
        df =  pd.read_csv( file_path, sep = "\t", index_col = 1)
        df = df[df['Adjusted P-value'] <= threshold_value]

        return df
    
    def compute_gsea(self,dir_ = "gsea_reactome/standard" ):
        dfs = []
        self.__load_reactome_graph__()

        for file in [self.validation_path + dir_ + "/" + file  for file in os.listdir(self.validation_path  + dir_ + "/") if file[0] != "."]:
            alg_name = file.split("/")[-1].split(".")[0]
            df = self.__load_gsea__(file)
            
            filtered_df = df.drop(["Unnamed: 0", "Genes"], axis=1)
            filtered_df[alg_name] = df['Adjusted P-value']
            filtered_df = filtered_df.drop(columns=['Adjusted P-value'])
            filtered_df[alg_name] = -np.log(filtered_df[alg_name])

            dfs.append(filtered_df)
        
        concat_df = pd.concat(dfs, axis=1)
        concat_df = concat_df.dropna().append(concat_df[concat_df.isna().any(axis=1)])
        concat_df.fillna(1, inplace=True)
        roots = []
        for term in concat_df.index:
            reactome_id = term.split(" ")[-1]
            try:
                self.visited = set()
                root = self.__find_root_path__(reactome_id)
                roots.append(root[-1])
            except:
                roots.append(np.nan)
        
        concat_df['root'] = roots
        concat_df = concat_df.sort_values(by='root', ascending=True)
        concat_df = concat_df.dropna()
        concat_df.to_csv("../imgs/gsea.tsv",sep="\t")
        
        concat_df = concat_df[concat_df['root'] == "R-HSA-162582"]
        concat_df = concat_df.drop(columns=['root'])
        concat_df.index = concat_df.index.str.replace(r' Homo sapiens R-HSA-\d+', '', regex=True)
        f,axes = plt.subplots(1,1,figsize = (4,8))
        sns_plot = sns.heatmap(concat_df, annot=False,ax = axes)
        colorbar = sns_plot.collections[0].colorbar
        colorbar.set_label(r'$-\log P_{value}$', fontsize=12)
        f.savefig('../imgs/reactome_heatmap.pdf',bbox_inches='tight')
        
    
    def compute_comparison(self, chosen_diseases = []):
        random_distribution_concats  = []
        metrics_concats  = []
       
        for dir_ in [self.validation_path + file + "/" for file in os.listdir(self.validation_path) if file[0] != "."]:
            
            if dir_ == "gsea_reactome" or dir_ == "mouse_phenotypes":
                continue 
            
            data_set = dir_.split("/")[-2]
            
            random_distribution_data_frame = self.__load_random_distribution__(dir_)
            random_distribution_data_frame['Dataset'] = data_set
            metrics_data_frame = self.__load_metrics__(dir_)
            metrics_data_frame['Dataset'] = data_set
            
            metrics_concats.append(metrics_data_frame)
            random_distribution_concats.append(random_distribution_data_frame)

        distribution_df = pd.concat(random_distribution_concats, axis=0)
        metrics_df = pd.concat(metrics_concats,axis=0)

        if (len(chosen_diseases) > 0):
            diseases = chosen_diseases
        else:
            diseases = list(metrics_df["GWAS"].unique())

        data_sets = list(metrics_df["Dataset"].unique())
        f,axes = plt.subplots(3,3,figsize = (20,20))
        
        for i,data_set in enumerate(data_sets):
            
            for j,disease in enumerate(diseases):
                
                try:

                    current_df = distribution_df.loc[(distribution_df['GWAS'] == disease) & (distribution_df['Dataset'] == data_set)]
                    current_metrics = metrics_df.loc[(metrics_df['GWAS'] == disease) & (metrics_df['Dataset'] == data_set)]
                    
                    sns.boxplot(data = current_df,x="Algorithm", y="Recall",showfliers = False,color='white',ax= axes[i][j])
                    sns.scatterplot(data = current_metrics,x = "Algorithm", y = "Recall",s=200, color=".2", hue= "Algorithm", style = "Algorithm",ax= axes[i][j])
                    axes[i][j].set_title('GWAS: ' + disease)
        
        
                
                except:
                    pass

        f.savefig('../imgs/algorithm_comparison.pdf',bbox_inches='tight')
    
    def show_mouse_phenotypes(self, chosen_diseases = []):
        random_distribution_concats  = []
        metrics_concats  = []
        for dir_ in [self.validation_path + file + "/" for file in os.listdir(self.validation_path) if file[0] != "."]:
            print(dir_)
            if "mouse_phenotypes" not in dir_:
                continue 
            
            data_set = dir_.split("/")[-2]
            
            random_distribution_data_frame = self.__load_random_distribution__(dir_)
            random_distribution_data_frame['Dataset'] = data_set
            metrics_data_frame = self.__load_metrics__(dir_)
            metrics_data_frame['Dataset'] = data_set
            
            metrics_concats.append(metrics_data_frame)
            random_distribution_concats.append(random_distribution_data_frame)

        distribution_df = pd.concat(random_distribution_concats, axis=0)
        metrics_df = pd.concat(metrics_concats,axis=0)

        f,axes = plt.subplots(1,1,figsize = (5,5))
        sns.boxplot(data = distribution_df,x="Algorithm", y="Recall",showfliers = False,color='white',ax= axes)
        sns.scatterplot(data = metrics_df,x = "Algorithm", y = "Recall",s=200, color=".2", hue= "Algorithm", style = "Algorithm",ax= axes)
        axes.set_title('Mouse Phenotypes')
        f.savefig('../imgs/mouse_phenotypes.pdf',bbox_inches='tight')

    
p = Plot(
    experiment_dir_path = "../experiments/algorithm_comparison/"
    )


p.compute_gsea()