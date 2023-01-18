'''
Fetch study information from NCBI GEO and prepare for preprocessing. 

Possible inputs include:
    - GSE accession number
    - List of GSE accession numbers
    - Studies from the last n days
    - search terms (e.g. organism and/or sequencing type)

Accession,Title,Organism,Samples,SRA,Release_Date,All_protocols,seq_types,GSE,GSE_Supplementary,BioProject,PMID,authors,abstract,title,doi,date_published,PMC,journal

'''

import argparse
from Bio import Entrez

def run_search(term: str ='Ribo-Seq[All Fields]', retmax: int=10) -> list:
    '''
    Run a search term on NCBI GEO and return the results.

    Parameters:
        term (str): Search term to run on NCBI GEO
        retmax (int): Maximum number of results to return

    Returns:
        gds_ids (list): List of GDS IDs
    '''
    Entrez.email = "riboseq@gmail.com"
    handle = Entrez.esearch(db="gds", term=term, retmax=retmax)
    gds_ids = Entrez.read(handle)["IdList"]
    
    return gds_ids

def fetch_study_information(gds_ids: list) -> dict:
    '''
    Fetch study information for a list of GDS IDs.

    Parameters:
        gds_ids (list): List of GDS IDs

    Returns:    
        study_information (dict): Dictionary of study information
    '''
    study_information = {"Accession": [], 
                        "Organism": [], 
                        "Title": [], 
                        "Samples": [], 
                        "SRA": [], 
                        "Release_Date": [], 
                        "All_protocols": [], 
                        "seq_types": [], 
                        "GSE": [], 
                        "GSE_Supplementary": [], 
                        "BioProject": [], 
                        "PMID": [], 
                        "authors": [], 
                        "abstract": [], 
                        "title": [], 
                        "doi": [], 
                        "date_published": [], 
                        "PMC": [], 
                        "journal": []}
    for id in gds_ids:
        handle = Entrez.esummary(db="gds", id=id)
        records = Entrez.read(handle)
        for i in records[0]:
            print(i, records[0][i])
        study_information["Accession"].append(records[0]["Accession"])
        study_information["Organism"].append(records[0]["taxon"])
        study_information["Title"].append(records[0]["title"])
        study_information["Samples"].append(int(records[0]["n_samples"]))
        # study_information["SRA"].append(records[0]["TargetFTPLink"])
        study_information["Release_Date"].append(records[0]["PDAT"])
        study_information["All_protocols"].append(records[0]["summary"])
        study_information["seq_types"].append(records[0]["gdsType"])
        study_information["GSE"].append(f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE{records[0]['GSE']}")
        study_information["GSE_Supplementary"].append(f"ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE{records[0]['GSE'][:-3]}nnn/GSE{records[0]['GSE']}/suppl/")
        print(study_information["GSE_Supplementary"])
        break


        


def main(args):
    '''
    Based on inputted arguments from argparse, fetch study information from NCBI GEO and prepare for preprocessing.
    '''
    if args.GSE:
        # Fetch study information for a single GSE accession number
        pass
    elif args.studies_csv:
        # Fetch study information for multiple GSE accession numbers
        pass
    elif args.days:
        # Fetch studies from the last n days
        pass
    elif args.search_terms:
        # Fetch studies using search terms
        pass
    else:
        # No arguments passed
        print('No arguments passed. Please run with -h for help.')
    
    study_information = {"Accession": [], 
                    "Organism": [], 
                    "Title": [], 
                    "Samples": [], 
                    "SRA": [], 
                    "Release_Date": [], 
                    "All_protocols": [], 
                    "seq_types": [], 
                    "GSE": [], 
                    "GSE_Supplementary": [], 
                    "BioProject": [], 
                    "PMID": [], 
                    "authors": [], 
                    "abstract": [], 
                    "title": [], 
                    "doi": [], 
                    "date_published": [], 
                    "PMC": [], 
                    "journal": []}
    gds_list = run_search()
    fetch_study_information(gds_list)


if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Fetch study information from NCBI GEO and prepare for preprocessing.')
    parser.add_argument('--GSE', type=str, help='Add a study to the database from the studies csv.')
    parser.add_argument('--studies_csv', help='Path to the study metadata csv file. Studies table is populated with data from this file.')
    parser.add_argument('--days', type=int, help='Get studies from the last n days')
    parser.add_argument('--search_terms', type=str, help='Search terms for GEO')
    args = parser.parse_args()
    main(args)
