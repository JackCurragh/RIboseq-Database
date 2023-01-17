'''
This script is used to create the database for the riboseqorg web application via the command line.
 The database is stored in riboseqorg.db but can be changed by passing the name of the database as 
 an argument to the script. The database is managed by the class DatabaseManager. This class contains 
 methods for creating the database, adding data to the database, and querying the database.

'''

from database_manager import DatabaseManager

import pandas as pd
import argparse

def main(args):
    '''
    This function is the main function for the script. It is called when the script is run from the command line.
    '''
    db = DatabaseManager(args.database_path)

    samples_fields = {
        'gse_accession' : '',
        'title' : '',
        'organism' : '',
        'num_samples' : '',
        'sra_accession' : '',
        'release_date' : '',
        'protocols' : '',
        'sequencing_types' : '',
        'gse_url' : '',
        'gse_supplementary' : '',
        'bioproject' : ''
    }

    if args.mode == 'F':
        studies_df = pd.read_csv(args.studies_csv)       
        GSE_row = studies_df[studies_df['Accession'] == args.GSE]
        row_number = GSE_row.index[0]

        samples_fields['gse_accession'] = GSE_row['Accession'][row_number]
        samples_fields['title'] = GSE_row['Title'][row_number]
        samples_fields['organism'] = GSE_row['Organism'][row_number]
        samples_fields['num_samples'] = GSE_row['Samples'][row_number]
        samples_fields['sra_accession'] = GSE_row['SRA'][row_number]
        samples_fields['release_date'] = GSE_row['Release_Date'][row_number]
        samples_fields['protocols'] = GSE_row['All_protocols'][row_number]
        samples_fields['sequencing_types'] = GSE_row['seq_types'][row_number]
        samples_fields['gse_url'] = GSE_row['GSE'][row_number]
        samples_fields['gse_supplementary'] = GSE_row['GSE_Supplementary'][row_number]
        samples_fields['bioproject'] = GSE_row['BioProject'][row_number]

        print(samples_fields['gse_accession'], type(samples_fields['num_samples']))
        
    elif args.mode == 'E':
        samples_fields['gse_accession'] = args.GSE
        samples_fields['title'] = args.title
        samples_fields['organism'] = args.organism
        samples_fields['num_samples'] = args.num_samples
        samples_fields['sra_accession'] = args.sra_accession
        samples_fields['release_date'] = args.release_date
        samples_fields['protocols'] = args.protocols
        samples_fields['sequencing_types'] = args.sequencing_types
        samples_fields['gse_url'] = args.gse_url
        samples_fields['gse_supplementary'] = args.gse_supplementary
        samples_fields['bioproject'] = args.bioproject


    db.add_study(gse_accession=samples_fields['gse_accession'], 
                title=samples_fields['title'], 
                organism=samples_fields['organism'], 
                num_samples=samples_fields['num_samples'], 
                sra_accession=samples_fields['sra_accession'], 
                release_date=samples_fields['release_date'], 
                protocols=samples_fields['protocols'], 
                sequencing_types=samples_fields['sequencing_types'], 
                gse=samples_fields['gse_url'], 
                gse_supplementary=samples_fields['gse_supplementary'], 
                bioproject=samples_fields['bioproject'], 
                samples_csv=args.samples_csv)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='CLI for managing the database for the riboseqorg web application')

    # Add arguments that are relevant to both modes
    parser.add_argument('--database_path', type=str, default='riboseqorg.db', help='The name of the database to create')
    parser.add_argument('--samples_csv', type=str, help='Path to the sample metadata csv file. Samples table is populated with data from this file.')
    parser.add_argument('--GSE', type=str, help='Add a study to the database from the studies csv.')

    # Add a "mode" argument that specifies which mode to use
    parser.add_argument('--mode', choices=['F', 'E'], default='F', help="The mode to use. F for study info from csv files and E explicitly get study info from inputs.")

    # Create argument groups for the two modes
    file_group = parser.add_argument_group('File mode')
    explicit_group = parser.add_argument_group('Explicit mode')

    # Add arguments that are only relevant to the "file" mode
    file_group.add_argument('--studies_csv', help='Path to the study metadata csv file. Studies table is populated with data from this file.')

    explicit_group.add_argument('--title', help="Title of the study as on GEO") 
    explicit_group.add_argument('--organism', help="Organism(s) in the study")
    explicit_group.add_argument('--num_samples', help="Number of samples in the study")
    explicit_group.add_argument('--sra_accession', help="SRA accession number for the study")
    explicit_group.add_argument('--release_date', help="Release date of the study")
    explicit_group.add_argument('--protocols', help="Protocols used in the study")
    explicit_group.add_argument('--sequencing_types', help="Sequencing types used in the study")
    explicit_group.add_argument('--gse_url', help="URL for the study on GEO")
    explicit_group.add_argument('--gse_supplementary', help="URL for the supplementary information on GEO")
    explicit_group.add_argument('--bioproject', help="URL to the BioProject page for the study")

    args = parser.parse_args()
    main(args)