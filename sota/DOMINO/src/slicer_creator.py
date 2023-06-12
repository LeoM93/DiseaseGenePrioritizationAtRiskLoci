import argparse
import os
from core.domino import main as domino_main
from core.preprocess_slices import create_slices
from utils.visualize_modules import visualize_modules
import constants as constants


def main_slicer():

    parser = argparse.ArgumentParser(description='Slicer for DOMINO (step #0): A preprocessing step for the network')
    parser.add_argument('-n', '--network_file', dest='network_file', help='A path to network file (sif format). e.g. /path/to/network_file.sif', default="examples/huri.sif")
    parser.add_argument('-o', '--output_file', dest='output_file', default="examples/huri.sif", help='A path to the output slices file. e.g., /path/to/output/slices_file.txt')


    args = parser.parse_args()
    network_file = args.network_file
    output_file = args.output_file
    create_slices(network_file, output_file)




if __name__=="__main__":
    main_slicer()
