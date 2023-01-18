# DiseaseGenePrioritizationAtRiskLoci
## Cite the following manuscript
If you desire to run the algorithm, please cite:

Martini, Leonardo, et al. "Network Based Approach to Gene Prioritization at Genome-Wide Association Study Loci." arXiv preprint arXiv:2210.16292 (2022).

## Instruction to recreate pictures and tables in the main manuscript
1. Download Datasets: download datasets.zip from https://drive.google.com/file/d/14CZfizfTPy1I7tZGoCVst2Rp9kp6-PJz/view?usp=sharing and unizip it in the main directory.
2. Run Computational Validation: from the main directory, go to evaluation/computational_validation and run on the terminal odd_vs_even.py and cross_study_validation.py
3. Run External Validation: from the main directory, go to evaluation/external_validation and run on the terminal python3 drug_hub_validation.py and python3 GSEA_validation.py python3 open_target_validation.py
4. Create plots: run on the terminal python3 plot.py

