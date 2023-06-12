import argparse
import os
from core.domino import main as domino_main
from core.preprocess_slices import create_slices
from utils.visualize_modules import visualize_modules
import constants as constants

def main_domino():

    parser = argparse.ArgumentParser(description='DOMINO: An active module identification algorithm with reduce rate of false.\n NOTE YOU SHOULD RUN THE SLICES SCRIPT FIRST! (more info, type slicer -h) \n Example input files are available @ https://github.com/Shamir-Lab/DOMINO/tree/master/examples')
    parser.add_argument('-a', '--active_genes_file', dest='active_genes_file', help='Comma delimited list of absolute paths to files, each containing a list of active genes, separated by a new line char (\\n). e.g. /path/to/active_genes_files_1,/path/to/active_genes_files_2.', default="examples/tnfa_active_genes_file.txt")
    parser.add_argument('-n', '--network_file', dest='network_file', help='A path to network file (sif format). e.g. /path/to/network_file.sif', default="examples/huri.sif")
    parser.add_argument('-s', '--slices_file', dest='slices_file', help='A path to slices file (i.e. the output of "slicer" script). e.g., /path/to/slices_file.txt', default="examples/huri_slices.txt")
    parser.add_argument('-o', '--output_folder', dest='output_folder', help='A folder where output files will be written e.g., /path/to/output', default="examples/output")
    parser.add_argument('-f', '--output_file_name', dest='output_file_name', help='A name of the file that will be written e.g., /path/to/output', default="examples/output")    
    parser.add_argument('-c', '--use_cache', dest='use_cache', help='Use auto-generated cache network files (*.pkl) from previous executions with the same network. NOTE: (1) THIS IS NOT THE SLICES FILE! (2) If the content of the file has changed, you should set this option to "false"', default="true")
    parser.add_argument('-p', '--parallelization', dest='parallelization', help='The number of threads allocated to the run (usually single thread is enough)', default="1")
    parser.add_argument('-v', '--visualization', dest='visualization', help='Indicates whether a visualization of the modules ought to be generated', default="true")
    parser.add_argument('-sth', '--slice_threshold', dest='slice_threshold', default="0.3", help='The threshold for considering a slice as relevant')
    parser.add_argument('-mth', '--module_threshold', dest='module_threshold', default="0.05", help='The threshold for considering a putative module as final module')


    args = parser.parse_args()
    cur_ag  = args.active_genes_file
    output_folder = args.output_folder
    
    network_file = args.network_file
    output_file_name = args.output_file_name
    slices_file = args.slices_file
    slice_threshold = float(args.slice_threshold)
    module_threshold = float(args.module_threshold)
    use_cache = args.use_cache=="false"
    parallelization = int(args.parallelization)
    visualization = args.visualization=="true"

    constants.N_OF_THREADS=parallelization
    constants.USE_CACHE=use_cache


    print(network_file)
    print(slices_file)
    print(cur_ag)

    file_name = cur_ag.split("/")[-1].split("__")[-1]
    G_final_modules=domino_main(active_genes_file=cur_ag, network_file=network_file, slices_file=slices_file, slice_threshold=slice_threshold, module_threshold=module_threshold)
        
    out_file= os.path.join(output_folder,output_file_name) 
        
    if len(G_final_modules) !=0:
        open(out_file, 'w+').write("\n".join(['[%s]' % ', '.join(list(m.nodes)) for m in G_final_modules])+"\n")
    else:
        open(out_file, 'w+').write("")

    print(f'{len(G_final_modules)} final modules are reported at {out_file}')
        


if __name__=="__main__":
    main_domino()
