import sys
import pandas as pd


def parse_entry(entry_multiline_string):
    '''
    Parse the 8 lines of GDS information for an individual entry 
    '''
    parsed_information = {}

    entry_lines = entry_multiline_string.split('\n')
    if entry_lines[-1].split('\t')[0] == 'Sample':
        return None

    parsed_information['title'] = ' '.join(entry_lines[0].split(" ")[1:])
    for line in entry_lines[1:]:
        if 'Organism:' in line:
            parsed_information['oragnism'] = line.split('\t')[1]
            
        elif 'Source Name:' in line:
            parsed_information['source_name'] = line.split('\t')[1]

        elif 'Type:' in line:
            parsed_information['dataset_type'] = line.split('\t')[1]

        elif 'Platform:' in line:
            parsed_information['platform'] = line.split(' ')[1]

            if 'Sample' in line:
                parsed_information['number_of_samples'] = line.split(' ')[-2]

            if 'Series:' in line:
                parsed_information['series'] = '; '.join(line.split(' ')[3:-1])

        elif 'Platforms:' in line:
            parsed_information['platform'] = '; '.join(line.split(' ')[1:-3])
            parsed_information['number_of_samples'] = line.split(' ')[-2]


        elif 'Platforms ' in line:
            parsed_information['platform'] = ' '.join(line.split(' ')[0:3])

            if 'Sample' in line:
                parsed_information['number_of_samples'] = line.split(' ')[-2]

        elif 'FTP' in line:
            parsed_information['ftp'] = line.split(' ')[-1]

        elif 'SRA Run Selector:' in line: 
            parsed_information['sra_run_selector'] = line.split(' ')[3]
            parsed_information['bioproject'] = line.split('=')[-1]
        
        elif 'Accession:' in line:
            parsed_information['GSE'] = line.split('\t')[2].split(': ')[-1]
    
    if 'bioproject' not in parsed_information: 
        parsed_information['bioproject'] = '' 
    
    if 'sra_run_selector' not in parsed_information:
        parsed_information['sra_run_selector'] = ''
        
    return parsed_information


def parse_txt_gds_output_to_csv(gds_txt_path, gds_csv_path):
    '''
    Read in the gds text output and output the same information into a csv file
    '''
    with open(gds_txt_path, 'r') as file:
        entry = []
        information_dict = {
            'GSE':[],
            'title':[], 
            'oragnism':[], 
            'dataset_type':[], 
            'platform':[], 
            'number_of_samples':[], 
            'ftp':[], 
            'sra_run_selector':[], 
            'bioproject':[], 
        }
        for line in file.readlines():
            if entry != [] and line == '\n':
                information = parse_entry('\n'.join(entry))
                if information:
                    for element in information:
                        information_dict[element].append(information[element])
                entry = []

            if line != "\n":
                entry.append(line.strip('\n'))
    df = pd.DataFrame(information_dict)
    df.to_csv(gds_csv_path, index=False)


if __name__ == "__main__":
    gds_txt_path = sys.argv[1]
    gds_csv_path = sys.argv[2]

    minimal_df = parse_txt_gds_output_to_csv(gds_txt_path, gds_csv_path)
    # minimal_df.to_csv(f"{gds_csv_path}_minimal_gds_results_table.tsv", sep="\t")



