# execute with an interactive -c 56 (56 cores)

# -- Import libraries -- # 

import os
import sys
import pandas as pd
import gzip                        
import json
import click
from multiprocessing import Pool
from itertools import repeat

## my modules ##
sys.path.append("../Utils/")    # modules folder
from fasta_utils import readFasta_gzip, readFasta_header_gzip
from motif_utils import motif_scan_v1 as motif_scan



# -- Main function options -- #

@click.command()

@click.option('--proteome_dir', 
	 		  '-p',
			  required = True,
			  help = "path to the folder where the proteome file is located")

@click.option('--weight_m_dir', 
			  '-wm',
			  required = True,
			  help = "path to the folder containing position weight matrices motifs")

@click.option('--scan_dir', 
	 		  '-s',
			  required = True,
			  help = "path to the folder to store motifs scans")

@click.option('--cpus', 
	 		  '-c',
			  required = False,
			  default = 1,
			  help = "number of cpus to be used")

@click.option('--proteome_format', 
	 		  '-pf',
			  required = False,
			  default = "fasta_header_gzip",
			  help = "proteome format, which can be json dictionary or fasta compressed")



# -- Main function  -- #


def motifs_scan_proteome(proteome_dir, weight_m_dir, scan_dir, cpus, proteome_format):
    """
    Scans a set of proteins (e.g.: proteome) using a sliding window technnique with a PWM.
    Returns a JSON file of the form {protein: scores}, so that every protein sequence
    has a list of subsequences scores. Proteins shorter than the PWM have an empty list associated.
    """

    print('\nReading the proteome\n')
    if proteome_format == "json":
        with open(proteome_dir) as fp:
            proteome = json.load(fp)
    elif proteome_format == "fasta_gzip":
        proteome = readFasta_gzip(proteome_dir)
    elif proteome_format == "fasta_header_gzip":
        proteome = readFasta_header_gzip(proteome_dir)
    else:
        return print("Incorrect proteome format, only JSON or fasta compressed")

    print('\n## SCAN STARTING ## \n')
    E3_ligases = os.listdir(weight_m_dir)
    counter = 0
    for E3_ligase in E3_ligases:
        counter +=1
        print(f'\t Scanning {E3_ligase}... ({counter}/{len(E3_ligases)})')

        weight_m_dict = {}
        weight_m = pd.read_csv(weight_m_dir+E3_ligase, sep = "\t")
        
        global scanning_proteome

        def scanning_proteome(protein_id):
            protein_seq = proteome[protein_id]
            protein_score = motif_scan(protein_seq, weight_m)
            return [protein_id, protein_score]

        with Pool(cpus) as p:
            weight_m_list = p.map(scanning_proteome, proteome.keys())

        weight_m_dict = {item[0]:item[1] for item in weight_m_list}

        print(f'\t Saving {E3_ligase} scan as JSON file')
        with open(scan_dir+E3_ligase.split(".")[0]+'.json', 'w') as fp:
            json.dump(weight_m_dict, fp)
    
    print('\n## SCAN FINISHED ## \n')


if __name__ == "__main__":
    motifs_scan_proteome()