import pandas as pd
import os
from functools import reduce

def generate_performance_table(validation_path = "../experiments/algorithm_comparison/validation/",):
    
    validation_dir = [validation_path + file + "/" for file in os.listdir(validation_path) if file[0] != "."]
    final_dfs = []
    for dir_ in validation_dir:
        if 'mouse_phenotypes' in dir_ or 'gsea_reactome' in dir_:
            continue
        
        files = [dir_+file for file in os.listdir(dir_) if file[0] != "."]

        db_name = dir_.split("/")[-2]
        dfs = []
        join_columns = ['Algorithm', 'GWAS', 'validation_data_set']


        for file in files:
            if 'metrics_p_val.tsv' not in file and 'metrics.tsv' not in file:
                continue
            df = pd.read_csv(file, index_col = 0,sep="\t")
            df['validation_data_set'] = db_name
            dfs.append(df)
        result = reduce(lambda left, right: pd.merge(left, right, on=join_columns, how='inner'), dfs)
        final_dfs.append(result)
        
    

    result = pd.concat(final_dfs, axis=0, ignore_index=True)
    result.to_csv('./prova.tsv', sep="\t")

generate_performance_table()

