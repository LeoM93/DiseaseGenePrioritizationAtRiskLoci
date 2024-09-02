# A Network-Based Approach to Gene Prioritization at Genome-Wide Association Study Loci

## Citation 

[1] Martini, Leonardo, et al. "Network Based Approach to Gene Prioritization at Genome-Wide Association Study Loci." arXiv preprint arXiv:2210.16292 (2022).

## Author: 

Leonardo Martini, Adriano Fazzone, Michele Gentili, Michael H. Cho, Luca Becchetti and Brian D. Hobbs

Contact:

Leonardo Martini (relma@channing.harvard.edu),  Michele Gentili (remge@channing.harvard.edu)

Department of Computer, Control, and Management Engineering Antonio Ruberti, Sapienza University of Rome, Rome, Italy

Channing Division of Network Medicine, Brigham and Women’s Hospital, Boston, MA,

Harvard Medical School, Boston, MA

CENTAI Institute, Turin, Italy.



## Algorithm Description

In this work, we propose the Relations-Maximization Method, or RMM-GWAS. It leverages co-regulatory networks and biological annotations to identify a functional set of genes likely to be associated with a complex disease’s risk loci. Our method builds on the hypothesis that disease-relevant genes across GWAS loci share biological processes, pathways, and correlated gene expression levels. For this reason, it integrates network topology and gene biological information to derive a weighted co-regulatory network. Then, given a complex disease with multiple associated risk loci, our method prioritizes a set of genes, one per disease locus, so that the overall number of co-regulatory interactions, each weighted by the number of shared biological annotations, is maximized. The predicted set of genes i) covers the set of disease loci and ii) forms a set of genes that, on average, share many co-regulatory interactions.


![alt text](https://github.com/LeoM93/DiseaseGenePrioritizationAtRiskLoci/blob/main/img/algorithm.png?raw=true)


## Python Libraries
Imported libraries and their version:

 - Python version: Python '3.10.4'
 - SciPy version '1.8.1'
 - NumPy version '1.22.4'
 - NetworkX version '3.1'

### Weighting the Gene Co-regulatory network

RMM-GWAS requires the input graph to be weighted. We use Resnik’s method [37] to assign a semantic similarity score to a pair of a and b of GO annotations. This method quantifies the similarity measure of a term pair (a, b) as the information content (IC) of their Most Informative Common Ancestor (MICA),

To compute the weighted graph, run the following commands:

 ```
 cd data_preprocessing
 ```
 
 ```
 python3 main.py <graph_input_file_path> <output_file_path> <seed_RMM-GWAS_file_path>
 ```
The files used in the manuscript can be downloaded at: https://drive.google.com/file/d/12oDaaEs1vso82UXsRe2AWeoGqNccZuLM/view?usp=sharing

## Running Relations-Maximization Method

The directory DiseaseGenePrioritizationAtRiskLoci/example contains the input graph we used in the manuscript

To run RMM-GWAS, pleas follow these commands:

 ```
 cd rmm_gwas python3 run_relations_maximization.py <graph_input_file_name> <output_file_path>
 ```
 


