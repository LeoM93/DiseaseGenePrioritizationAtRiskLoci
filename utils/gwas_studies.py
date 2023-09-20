import requests
import csv
import pandas as pd


class GWASCatalogScraper():
    def __init__(self, gwas_directory):
        self.gwas_dataset = '../datasets/GWAS/gwas_catalog_v1.0.2-studies_r2022-06-15.tsv'
        self.gwas_directory = gwas_directory
        
    
    def __load_filtered_gwas__(self,count = 100):
        csv_reader = csv.reader(open(self.gwas_dataset), delimiter = "\t")
        gwas_studies = set()
        
        for index, row in enumerate(csv_reader):
            if index == 0:
                continue
            accession_number = row[-2]
            accession_count = int(row[-5])
            disease = row[7]
            if 'cancer' in disease and accession_count > count:
                gwas_studies.add(accession_number)
        
        return gwas_studies
    
    def __download_gwas_study__(self,accession_number):
        # Define the URL
        url = 'https://www.ebi.ac.uk/gwas/api/search/downloads?q=accessionId:' + accession_number + '&pvalfilter=&orfilter=&betafilter=&datefilter=&genomicfilter=&genotypingfilter[]=&traitfilter[]=&dateaddedfilter=&facet=association&efo=true'

        try:
            # Send a GET request to the URL
            response = requests.get(url)

            # Check if the request was successful (status code 200)
            if response.status_code == 200:

                text_data = response.text
                file_name = self.gwas_directory + accession_number + '.tsv'
                # Save the text data to a TSV file
                with open(file_name, 'w', encoding='utf-8') as tsv_file:
                    tsv_file.write(text_data)
                
                data_frame = pd.read_csv(file_name, sep ='\t')
                columns_to_keep = ['CHR_ID', 'CHR_POS', 'SNPS', 'P-VALUE', 'UPSTREAM_GENE_ID', 'DOWNSTREAM_GENE_ID', 'SNP_GENE_IDS']
                filtered_df = data_frame[columns_to_keep]
                
                filtered_df.to_csv(file_name, sep='\t', index=False)

            
                
            
            else:
                print(f"Request failed with status code: {response.status_code}")
        except Exception as e:
            print(f"An error occurred: {str(e)}")
    
    def run(self,):
        gwas_studies = self.__load_filtered_gwas__()
        for gwas in gwas_studies:
            self.__download_gwas_study__(gwas)
        

        

c = GWASCatalogScraper(
    gwas_directory = '../experiments/GWAS/studies/'
)
c.run()