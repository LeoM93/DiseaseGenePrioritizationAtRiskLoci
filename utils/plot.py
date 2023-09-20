import os
import math
import csv
import pandas as pd
import seaborn as sns
sns.set_theme(style="whitegrid")
import string

import matplotlib.pyplot as plt
import matplotlib as mpl
class Plot():
    def __init__(self, experiment_dir_path):
        self.experiment_dir_path = experiment_dir_path
        self.validation_path = experiment_dir_path + "validation/"
    

    def __load_pvals__(self, file_name = "metrics_p_val.tsv"):
        pass
    
    def __load_metrics__(self, file_name = "metrics.tsv"):
        return pd.read_csv(self.current_validation_path + file_name, sep = "\t", index_col = 0)

    
    def __load_random_distribution__(self, file_name = "metrics_random_distribution.tsv"):
        return pd.read_csv(self.current_validation_path + file_name, sep = "\t", index_col = 0)
    
    def run(self,dir_ = "cosmic"):
        
        self.current_validation_path = self.validation_path + dir_ + "/"
        
        random_distribution_data_frame = self.__load_random_distribution__()
        metrics_df = self.__load_metrics__()
        f,axes = plt.subplots(3,3,figsize = (20,20))
        diseases = list(random_distribution_data_frame["GWAS"].unique())
        
        for i in range(3):
            for j in range(3):
                index = 3*i + j
                try:
                    disease = diseases[index]
                    current_df = random_distribution_data_frame.loc[(random_distribution_data_frame['GWAS'] == disease)]
                    print(current_df)
                    current_metrics = metrics_df.loc[(metrics_df['GWAS'] == disease) ]
                    sns.boxplot(data = current_df,x="Algorithm", y="Recall",showfliers = False,color='white',ax = axes[i][j])
                    sns.scatterplot(data = current_metrics,x = "Algorithm", y = "Recall",s=200, color=".2", hue= "Algorithm", style = "Algorithm",ax= axes[i][j])
                    axes[i][j].set_title(disease)
                except:
                    pass
        f.savefig('../imgs/' + dir_ + '.pdf',bbox_inches='tight')

p = Plot(
    experiment_dir_path = "../experiments/algorithm_comparison/"
    )

p.run()