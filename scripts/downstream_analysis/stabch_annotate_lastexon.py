# -- Import libraries -- # 

import os
import pandas as pd
import click
from tqdm import tqdm
tqdm.pandas()

# -- Auxiliary functions -- #

def isinterval(position): 
    """
    Takes a genomic position and determines if it is an 
    interval

    Parameters
    ----------
    position: str
            Posible interval

    Returns
    -------
    True if interval
    False if not interval
    """
    try:
        for p in position:
            if p == "-":
                return True
    except:
        return False

def annot_mut_last_exon(row, last_exons, symb_colname):
    """
    Adds True/False depending on the mutation falling in the
    sequence of the last gene's exon

    Parameters
    ----------
    row: Dataframe row
           Row from the stability change dataframe. 
           Mandatory columns: gene/Hugo_Symbol, Location/bgvep_genomic_loc, Feature
    last_exons: Dataframe
            Contains the genomic coordinates of the last exons of the canonical
            transcripts genes.
            Mandatory columns: ENST, abs_start, abs_end
    
    Returns
    -------
    row with the Mut_in_lastexon column added with a True or False value
    """
    
    # Only process mutations
    if row["Phenotype"] != "WT":
        enst = row["Feature"]

        # Retrieved exons from biomart
        if enst in last_exons["ENST"].values:

            last_exon_abs_start = int(last_exons.loc[last_exons["ENST"] == enst, "abs_start"].item())
            last_exon_abs_end = int(last_exons.loc[last_exons["ENST"] == enst, "abs_end"].item())

            # MS
            if symb_colname == "gene":
                mut_pos = row["Location"]
                # Avoid nan and .
                if isinstance(mut_pos, str):
                    if mut_pos != ".":
                        if isinterval(mut_pos):
                            mut_pos = [int(i) for i in mut_pos.split(":")[1].split("-")]  # Format -> chr#:######-#######
                        else: 
                            mut_pos = int(mut_pos.split(":")[1])  # Format -> chr#:#######

            # RPPA
            elif symb_colname == "Hugo_Symbol":
                # Avoid nan and .
                if isinstance(row["bgvep_genomic_loc"], str):
                    if row["bgvep_genomic_loc"] != ".":
                        mut_pos = int(row["bgvep_genomic_loc"])
            
            # Check if mut location falls in the last exon coordinates
            if isinstance(mut_pos, list): # no SNV
                if (mut_pos[0] >= last_exon_abs_start) and (mut_pos[1] <= last_exon_abs_end):
                    row["Mut_in_lastexon"] = True
                else:
                    row["Mut_in_lastexon"] = False
            else:
                if (mut_pos >= last_exon_abs_start) and (mut_pos <= last_exon_abs_end):
                    row["Mut_in_lastexon"] = True
                else:
                    row["Mut_in_lastexon"] = False
        
        else:
            row["Mut_in_lastexon"] = False
    
    # WT
    else:
        row["Mut_in_lastexon"] = False
    
    return row
        
        
# -- Main function options -- #

@click.command()

@click.option('--stabch_dir', 
	 		  '-stabch',
			  required = True,
			  help = "path to the stability change files directory (ENST annotated + genomic coords in rppa)")

@click.option('--symb_colname', 
	 		  '-symb',
              type = click.Choice(["gene", "Hugo_Symbol"]), 
              default = "gene",
			  required = True,
			  help = "gene symbol column name in the stability change file ('gene' or 'Hugo_Symbol')")

@click.option('--lastexons_dir', 
	 		  '-lexons',
			  required = True,
              default = "/home/rblanco/projects/degrons/data/external/biomart/biomart92_proteome_lastexons_cantranscripts_relcoords_valid.fasta.gz",
			  help = "path to the file containing the last exons with their absolute and relative coordinates") 

@click.option('--stabch_annot_dir', 
	 		  '-annot',
			  required = True,
			  help = "path to the directory where the stability change files with last exon mutations annotated\
                  will be located")

# -- Main function  -- #

def stabch_annotate_lastexon(stabch_dir, symb_colname, lastexons_dir, stabch_annot_dir):
    """
    Adds to the stability change tables a column indicating whether a mutation falls in the last
    protein exon or not (True or False) 
    """

    print("Loading last exons file...")
    last_exons = pd.read_csv(lastexons_dir, sep = "\t", compression = "gzip")

    # CCLE files
    if os.path.isfile(stabch_dir):   
        print("Starting CCLE pipeline (one file: pancancer)\n")

        stabch_df = pd.read_csv(stabch_dir, sep = "\t", compression = "gzip")

        print(f"\tAnnotating mutations happening in the last exon (True/False)")
        stabch_df = stabch_df.progress_apply(lambda row: annot_mut_last_exon(row, last_exons, symb_colname), axis = 1)

        stabch_df.to_csv(stabch_annot_dir, sep = "\t", index = False, compression = "gzip")

    # TCGA and CPTAC files
    elif os.path.isdir(stabch_dir):
        print("Starting TCGA or CPTAC pipeline (several files: cancer types)\n")

        print("\nDefining cancer types...\n")
        ctypes = os.listdir(stabch_dir)       

        c = 0

        print(f"Annotating mutations happening in the last exon (True/False) for {len(ctypes)} cancer types")
        print("##############################################################################################\n")

        for ctype in ctypes:

            c += 1

            print(f"\tPreprocessing of {ctype.split('.')[0]}... ({c}/{len(ctypes)})")
            stabch_df = pd.read_csv(stabch_dir+ctype, sep = "\t", compression = "gzip")

            print(f"\tAnnotating mutations happening in the last exon (True/False)")
            stabch_df = stabch_df.progress_apply(lambda row: annot_mut_last_exon(row, last_exons, symb_colname), axis = 1)

            stabch_df.to_csv(stabch_annot_dir+ctype, sep = "\t", index = False, compression = "gzip")

    print("\n##################")
    print("Analysis finished")


if __name__ == "__main__":
    stabch_annotate_lastexon() 