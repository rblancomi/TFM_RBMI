# -- Import libraries -- # 

import os
import sys
import pandas as pd
import json
import gzip
import click

## my modules ##
sys.path.append("./Utils/")    # modules folder
from fasta_utils import readFasta_gzip, readFasta_header_gzip

# -- Main function options -- #

@click.command()

@click.option('--motifs_metrics_dir', 
	 		  '-metrics',
			  required = True,
			  help = "path to the motifs metrics file")

@click.option('--scan_dir', 
	 		  '-s',
			  required = True,
			  help = "path to the folder where motifs proteome-wide scans are stored")

@click.option('--proteome_dir', 
	 		  '-p',
			  required = True,
			  help = "path to the folder where the proteome file is located")

@click.option('--new_instances_discovered_dir', 
	 		  '-disc',
			  required = True,
			  help = "path to the new instances discovered folder")

@click.option('--proteome_format', 
	 		  '-pf',
			  required = False,
			  default = "fasta_header_gzip",
			  help = "proteome format, which can be json dictionary or fasta compressed")

# -- Main function  -- #

def motifs_new_instances_discovered(motifs_metrics_dir, scan_dir, proteome_dir, new_instances_discovered_dir, proteome_format):

    print("\nLoad motifs IDs\n")
    motifs_IDs = [motif.split(".")[0] for motif in os.listdir(scan_dir)]
    motifs_IDs.remove("APC_DBOX")  # cohort.gz ends up with 9409722 rows after adding discovered degrons information

    print("Load motif metrics file to extract the positivity range per motif\n")
    metrics = pd.read_csv(motifs_metrics_dir, sep = "\t", compression = 'gzip')

    print("Load proteome")
    if proteome_format == "json":
        with open(proteome_dir) as fp:
            proteome = json.load(fp)
    elif proteome_format == "fasta_gzip":
        proteome = readFasta_gzip(proteome_dir)
    elif proteome_format == "fasta_header_gzip":
        proteome = readFasta_header_gzip(proteome_dir)
    else:
        return print("Incorrect proteome format, only JSON or fasta compressed")

    counter = 0

    for motif in motifs_IDs:

        counter += 1
        print(f"Extracting {motif} new discovered instances... ({counter}/{len(motifs_IDs)})\n")

        # Lists to be calculated
        proteins_id = []
        sequences = []
        starts = []
        ends = []

        # Load files
        print(f'\tExtract {motif} positivity range and motif length')
        pos_range = metrics.loc[metrics["motif_id"] == motif, "positivity_range"].values[0] 
        low_pos_range = float(pos_range.split(",")[0][1:])   # tuples are loaded as strings in the df
        motif_length = metrics.loc[metrics["motif_id"] == motif, "motif_length"].values[0]        

        print("\tLoad motif scan (dictionary)")
        with open(scan_dir+motif+".json") as fp:
            scan = json.load(fp)
        
        print("\tExtracting...")
        for protein in scan:
            for idx, score in enumerate(scan[protein]):   # enumerate starts in zero
                if (score >= low_pos_range):
                    
                    # If the header is complete, retrieve ENST, which is in 2nd position
                    if proteome_format == "fasta_header_gzip":
                        proteins_id.append(protein.split("|")[1])
                    else:
                        proteins_id.append(protein)
                    sequences.append(proteome[protein][idx:idx+motif_length])   # checked this is the correct way of indexing
                    starts.append(idx+1)   # adjust to non-pythonic indexing (sequence coordinates)
                    ends.append(idx+motif_length)
        
        print("\tGenerate dataframe with [protein_id, sequence, start, end] for discovered instances\n")
        new_instances_discovered = pd.DataFrame(list(zip(proteins_id,
                                            sequences,
                                            starts,
                                            ends)),
                                            columns = ["protein_id", 
                                            "sequence",
                                            "start",
                                            "end"])

        new_instances_discovered.to_csv(new_instances_discovered_dir+motif+".tsv.gz", sep = "\t", index = False, compression = 'gzip')

if __name__ == "__main__":
    motifs_new_instances_discovered()
