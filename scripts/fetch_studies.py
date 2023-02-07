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
import pandas as pd


def run_search_gds(term: str ='Ribo-Seq[All Fields]', retmax: int=10) -> list:
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


def run_search_pubmed(term: str, retmax: int=1) -> list:
    '''
    Run a search term on NCBI pubmed and return the results.

    Parameters:
        term (str): Search term to run on NCBI pubmed
        retmax (int): Maximum number of results to return

    Returns:
        gds_ids (list): List of GDS IDs
    '''
    Entrez.email = "riboseq@gmail.com"
    handle = Entrez.esearch(db="pubmed", term=term, retmax=retmax)
    pubmed_ids = Entrez.read(handle)["IdList"]

    return pubmed_ids

def fetch_information_pubmed(study_information: dict) -> list:
    '''
    Run a search term on NCBI PubMed and return the results.

    Parameters:
        term (str): Search term to run on NCBI GEO
        retmax (int): Maximum number of results to return

    Returns:
        gds_ids (list): List of GDS IDs
    '''
    Entrez.email = "riboseq@gmail.com"
    for idx, entry in enumerate(study_information["PMID"]):
        if entry != "NA":
            handle = Entrez.esummary(db="pubmed", id=entry, retmax=1)
            records = Entrez.read(handle)
            study_information["authors"][idx] = f'{" ".join(records[0]["AuthorList"])}'
            study_information["title"][idx] = f'{records[0]["Title"]}'
            study_information["doi"][idx] = f'{records[0]["ArticleIds"]["doi"]}'
            study_information["date_published"][idx] = f'{records[0]["PubDate"]}'
            if "pmc" in records[0]["ArticleIds"]:
                study_information["PMC"][idx] = f'{records[0]["ArticleIds"]["pmc"]}'
            else:
                study_information["PMC"][idx] = "NA"
            study_information["journal"][idx] = f'{records[0]["FullJournalName"]}'
        
        else:
            pubmed_ids = run_search_pubmed(study_information['Title'][idx])
            if len(pubmed_ids) > 0:
                print(pubmed_ids[0])
                handle = Entrez.esummary(db="pubmed", id=pubmed_ids[0], retmax=1)
                records = Entrez.read(handle)
                study_information["authors"][idx] = f'{" ".join(records[0]["AuthorList"])}'
                study_information["title"][idx] = f'{records[0]["Title"]}'
                study_information["doi"][idx] = f'{records[0]["ArticleIds"]["doi"]}'
                study_information["date_published"][idx] = f'{records[0]["PubDate"]}'
                if "pmc" in records[0]["ArticleIds"]:
                    study_information["PMC"][idx] = f'{records[0]["ArticleIds"]["pmc"]}'
                else:
                    study_information["PMC"][idx] = "NA"
                study_information["journal"][idx] = f'{records[0]["FullJournalName"]}'
        

    return study_information



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
                        "title": [], 
                        "doi": [], 
                        "date_published": [], 
                        "PMC": [], 
                        "journal": []}

    for id in gds_ids:
        handle = Entrez.esummary(db="gds", id=id)
        records = Entrez.read(handle)
        study_information["Accession"].append(f'{records[0]["Accession"]}')
        study_information["Organism"].append(f'{records[0]["taxon"]}')
        study_information["Title"].append(f'{records[0]["title"]}')
        study_information["Samples"].append(f'{int(records[0]["n_samples"])}')

        if records[0]["ExtRelations"] != []:
            if records[0]["ExtRelations"][0]['RelationType'].casefold() == "sra".casefold():
                study_information["SRA"].append(records[0]["ExtRelations"][0]['TargetObject'])
            else:
                study_information["SRA"].append("NA")
        else:
            study_information["SRA"].append("NA")
        study_information["Release_Date"].append(f'{records[0]["PDAT"]}')
        study_information["All_protocols"].append(f'{records[0]["summary"]}')
        study_information["seq_types"].append(f'{records[0]["gdsType"]}')
        study_information["GSE"].append(f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE{records[0]['GSE']}")
        study_information["GSE_Supplementary"].append(f"{records[0]['FTPLink']}suppl/")
        if records[0]['Projects'] == []:
            study_information["BioProject"].append("NA")
        else:
            study_information["BioProject"].append(f"https://www.ncbi.nlm.nih.gov/bioproject/{records[0]['Projects'][0]}")
        if records[0]["PubMedIds"] == []:
            study_information["PMID"].append("NA")
        else:
            study_information["PMID"].append(f'{int(records[0]["PubMedIds"][0])}')
        
        study_information['authors'].append("NA")
        study_information['title'].append("NA")
        study_information['doi'].append("NA")
        study_information['date_published'].append("NA")
        study_information['PMC'].append("NA")
        study_information['journal'].append("NA")


    return study_information
        


def main(args):
    '''
    Based on inputted arguments from argparse, fetch study information from NCBI GEO and prepare for preprocessing.
    '''
    if args.GSE:
        # Fetch study information for a single GSE accession number
        gds_list = run_search_gds(term=f'{args.GSE}[All Fields]', retmax=1)

    elif args.studies_csv:
        # Fetch study information for multiple GSE accession numbers
        pass
    elif args.days:
        # Fetch studies from the last n days
        pass
    elif args.search_terms:
        gds_list = run_search_gds(term=f'{args.search_terms}[All Fields]', retmax=10)
        # Fetch studies using search terms
        pass
    else:
        gds_list = run_search_gds(retmax=100)

        # No arguments passed
        print('No arguments passed. Please run with -h for help.')

    
    study_information = fetch_study_information(gds_list)

    study_information = fetch_information_pubmed(study_information)


    print(pd.DataFrame(study_information))


if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Fetch study information from NCBI GEO and prepare for preprocessing.')
    parser.add_argument('--GSE', type=str, help='Add a study to the database from the studies csv.')
    parser.add_argument('--studies_csv', help='Path to the study metadata csv file. Studies table is populated with data from this file.')
    parser.add_argument('--days', type=int, help='Get studies from the last n days')
    parser.add_argument('--search_terms', type=str, help='Search terms for GEO')
    args = parser.parse_args()
    main(args)
