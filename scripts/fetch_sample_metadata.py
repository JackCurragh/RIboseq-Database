'''
This script is used to retrieve the metadata for a specific sample from Entrez.
 It is used to retrieve the organism, the cell line, the strain, the tissue, the 
 library strategy, the library source, the library selection, the library layout, 
 the platform, the title, the description, the protocol, the tags and the source name.


'''

import argparse
from Bio import Entrez
import pandas as pd

def run_search_biosample(term: str, retmax: int=10) -> dict:
    '''
    Return sample information for a specific sample
    Parameters:
        term (str): Search term to run on biosample (GSM ID, SRA ID, etc.)
        retmax (int): Maximum number of results to return

    Returns:
        gds_ids (list): List of GDS IDs
    '''
    Entrez.email = "riboseq@gmail.com"
    handle = Entrez.esearch(db="biosample", term=term, retmax=retmax)

    biosample_ids = Entrez.read(handle)["IdList"]

    for id in biosample_ids:
        handle = Entrez.esummary(db="biosample", id=id, retmax=1)
        records = Entrez.read(handle)

    return dict(records['DocumentSummarySet']['DocumentSummary'][0])


def main(args):

    if args.GSE:
        records = run_search_biosample(args.GSE)
    elif args.SRP:
        records = run_search_biosample(args.SRP)

    if records:
        records = {k: [v] for k, v in records.items()}
        records_df = pd.DataFrame.from_dict(records)
        records_df.to_csv( f'{args.output}', index=True)
    else:
        raise Exception("No records found")

if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Fetch sample information')
    parser.add_argument('--GSE', type=str, help='Studies GSE accession number.')
    parser.add_argument('--SRP', type=str, help='Study SRP accession number.')
    parser.add_argument('--output', type=str, help='Output file name.')

    args = parser.parse_args()
    main(args)
