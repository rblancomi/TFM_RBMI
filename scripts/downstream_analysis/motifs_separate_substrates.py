# -- Import libraries -- # 

import os
import pandas as pd
import numpy as np
from tqdm import tqdm              
import gzip                        
import json
import click
from copy import deepcopy


# -- Auxiliary functions -- #

def load_ESIs(pos_set_dir):
    """
    Reads a dataframe of E3 ligase-substrate interactions (ESIs) that contains a "gene" column
    with the ID of the E3 ligase substrate and generates a ESIs dictionary
    
    Parameters
    ----------
    pos_set_dir: str
           Path to the folder where dataframes of ESIs per motif are stored
    
    Returns
    -------
    ESIs: dict
            Dictionary containing motifs as keys and corresponding sets of substrates as values
    """

    E3s = [E3.split(".")[0] for E3 in os.listdir(pos_set_dir)]

    ESIs = {}

    for E3 in E3s:
        msa_df = pd.read_csv(pos_set_dir+E3+".tsv", sep = "\t")
        ESIs[E3] = set(msa_df.gene)
    
    return ESIs


def motif_sep_substrates(E3, ESIs, scan_dir, substrates_dir, rest_of_seqs_dir, log_file):
    """
    Reads a motif's scan file and divides in two scan files: substrates and no-substrates.
    
    Parameters
    ----------
    E3: str
           E3 ligase motif to process
    ESIs: dict
           Dictionary containing motifs as keys and corresponding sets of substrates as values
    scan_dir: str
           Path to the folder containing the E3 ligase motif scan
    substrates_dir: str
           Path to the folder where the substrates-filtered scan will be stored
    rest_of_seqs_dir: str
           Path to the folder where the scan with the non-substrates will be stored
    log_file: str
           Path to the folder where the log file of the process will be stored

    Returns
    -------
    None
    """

    print(f'E3-ligase AC: {E3}', file = log_file)
    print(f'E3-ligase AC: {E3}')

    # Load E3 scanner file
    with open(scan_dir+E3+".json") as fp:
        E3_scan = json.load(fp)
    
    # Substrates and non-substrates objects
    rest_of_seqs = deepcopy(E3_scan)
    substrates = {}

    print(f'Number of scanned proteins (proteome): {len(E3_scan)}', file = log_file)
    print(f'Number of scanned proteins (proteome): {len(E3_scan)}')
    print(f'Expected number of E3-ligase substrates: {len(ESIs[E3])}', file = log_file)
    print(f'Expected number of E3-ligase substrates: {len(ESIs[E3])}')
    print(f'Expected number of remaining sequences: {len(E3_scan) - len(ESIs[E3])}', 
                                                    file = log_file)
    print(f'Expected number of remaining sequences: {len(E3_scan) - len(ESIs[E3])}')

    # Substrate and non-substrate dictionaries generation
    for substrate in ESIs[E3]:
        for protein in E3_scan:
            if substrate == protein:
                substrates[protein] = E3_scan[protein]
                del rest_of_seqs[protein]
                break

    # Substrate and non-substrate files generation
    print(f'Saved number of E3-ligase substrates {len(substrates)}', file = log_file)
    print(f'Saved number of E3-ligase substrates {len(substrates)}')
    with open(substrates_dir+E3+".json", 'w') as fp:
        json.dump(substrates, fp)
    
    print(f'Saved number of E3-ligase non-substrates: {len(rest_of_seqs)}', file = log_file)
    print(f'Saved number of E3-ligase non-substrates: {len(rest_of_seqs)}')
    with open(rest_of_seqs_dir+E3+".json", 'w') as fp:
        json.dump(rest_of_seqs, fp)
    
    print("----------------------------------------------------\n", file = log_file)
    print("----------------------------------------------------\n")
    
    return None

    

# -- Main function options -- #

@click.command()

@click.option('--scan_dir', 
	 		  '-s',
			  required = True,
			  help = "path to the folder where motifs scans are stored")

@click.option('--pos_set_dir', 
	 		  '-ps',
			  required = True,
			  help = "path to the folder where positive sequences files are stored")

@click.option('--log_dir', 
	 		  '-l',
			  required = True,
			  help = "path to the folder to store logs file")

@click.option('--scan_subs_dir', 
	 		  '-subs',
			  required = True,
			  help = "path to the folder to store substrates scans")

@click.option('--scan_no_subs_dir', 
	 		  '-nosubs',
			  required = True,
			  help = "path to the folder to store rest-of-sequences scans")



# -- Main function  -- #

def motifs_separate_substrates(scan_dir, pos_set_dir, log_dir, scan_subs_dir, 
                            scan_no_subs_dir):
    """
    Takes a PWMs scans and generates substrates and non-substrates scans
    """

    print('\nRetrieving E3-ligases ACs to be analysed\n')
    E3s = [E3.split(".")[0] for E3 in os.listdir(scan_dir)]

    print('Generating E3-substrate interactions (ESIs) dictionary\n')
    ESIs = load_ESIs(pos_set_dir)

    print('STARTING FILES GENERATION\n')
    counter = 1
    with open(log_dir, 'a') as log_file:
        for E3 in E3s:
            print(str(counter)+'/'+str(len(E3s)))
            motif_sep_substrates(E3, ESIs, scan_dir, scan_subs_dir, 
                                        scan_no_subs_dir, log_file)
            counter += 1
    
    return None

if __name__ == "__main__":
    motifs_separate_substrates()





