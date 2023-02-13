'''
This script is used to retrieve the metadata for a specific sample from Entrez.
 It is used to retrieve the organism, the cell line, the strain, the tissue, the 
 library strategy, the library source, the library selection, the library layout, 
 the platform, the title, the description, the protocol, the tags and the source name.




'''

import argparse
from Bio import Entrez
import pandas as pd

def run_search_biosample(term: str, retmax: int=10) -> list:
    '''
    return sample information for a specific sample
    Parameters:
        term (str): Search term to run on NCBI GEO
        retmax (int): Maximum number of results to return

    Returns:
        gds_ids (list): List of GDS IDs
    '''
    Entrez.email = "riboseq@gmail.com"
    handle = Entrez.esearch(db="bioproject", term=term, retmax=retmax)

    biosample_ids = Entrez.read(handle)["IdList"]
    print(biosample_ids)
    
    for i in biosample_ids:
        handle = Entrez.esummary(db="bioproject", id=i, retmax=1)
        records = Entrez.read(handle)
        print(records)
    return biosample_ids


def main(args):

    if args.GSE:
        run_search_biosample(args.GSE)
    elif args.SRP:
        run_search_biosample(args.SRP)


    return False

if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Fetch sample information')
    parser.add_argument('--GSE', type=str, help='Studies GSE accession number.')
    parser.add_argument('--SRP', type=str, help='Sample SRP accession number.')

    args = parser.parse_args()
    main(args)
