# -- Import libraries -- # 

import os
import pandas as pd
import numpy as np
import math
import click
from tqdm import tqdm
tqdm.pandas()
pd.set_option('display.max_rows', None)

# -- Auxiliary functions -- #

def annot_symb(row, mapping, rep_symbs):
    """
    Adds the corresponding ENST ID based on the gene symbol in those rows
    where the ENST is missing OR if there is not ENST column in the df

    Parameters
    ----------
    row: Dataframe row
           Row from the mutations dataframe from a specific tumor-type cohort. 
           Mandatory columns: gene, Feature
    mapping: Dataframe
            Contains the correspondance between ENST and gene symbol.
            Mandatory columns: symbol, enst
    
    Returns
    -------
    row with the ENST added if applicable 
    """

    # For df not having Feature column (RPPA)
    if "Feature" not in row.index:

        symb = row["Hugo_Symbol"]  # Tables not having Feature column, have Hugo_Symbol column to store the gene symbol
        
        # Avoid including annotating ENSTs in repetead symbols
        # Also, avoid errors due to symbols not having their counterpart in the mapping file
        if (symb not in rep_symbs) and (symb in mapping.symbol.unique()):
            #print(mapping.loc[mapping["symbol"] == symb, "enst"])
            row["Feature"] = mapping.loc[mapping["symbol"] == symb, "enst"].item()
        
        else:  # we need to add something to the column, a NA, to later remove this row from the df
            row["Feature"] = np.NaN
        
        return row

    
    # For df having Feature column (MS)
    elif "Feature" in row.index:

        # The feature field has a string (100% is not NA)
        if isinstance(row["Feature"], str):
            
            # Make sure the feature field contains a ENST
            if row["Feature"][:4] == "ENST":
                return row
            # Check what the cell contains in case is not a ENST
            else:
                print("The feature is a str but does not start with ENST")
                print(row["Feature"])
                return row
        
        # The feature field contains a missing value
        elif math.isnan(row["Feature"]):
            
            symb = row["gene"]

            # Avoid including annotating ENSTs in repetead symbols
            # Also, avoid errors due to symbols not having their counterpart in the mapping file
            if (symb not in rep_symbs) and (symb in mapping.symbol.unique()):
                row["Feature"] = mapping.loc[mapping["symbol"] == symb, "enst"].item()
            else:  # we need to add something to the column, a NA, to later remove this row from the df
                row["Feature"] = np.NaN
        
            return row
        
        # None of the above happened
        else:
            print("Strange thing 2")
            print(row["Feature"])
            return row

def check_reps(df, colname):
    """
    Counts the number of times a value is repeated in a column
    and counts and returns a list with those values repeated 
    two times or more

    Parameters
    ----------
    df: Dataframe
           Dataframe from where to check the replicates of a specific
           column
    colname: str
            Name of the column from where to check the replicated
    
    Returns
    -------
    reps: list
            Replicated items
    """
    
    counts = df[colname].value_counts(ascending=True)
    reps = counts[counts > 1].index
    
    return reps

# -- Main function options -- #

@click.command()

@click.option('--stabch_dir', 
	 		  '-stabch',
			  required = True,
			  help = "path to the stability change files directory (those in the cluster)")

@click.option('--symb_colname', 
	 		  '-symb',
              type = click.Choice(["gene", "Hugo_Symbol"]), 
			  required = True,
			  help = "gene symbol column name in the stability change file ('gene' or 'Hugo_Symbol')")

@click.option('--mapping_dir', 
	 		  '-map',
			  required = True,
			  help = "path to the file containing the ENST-gene symbol mapping") 

@click.option('--stabch_annot_dir', 
	 		  '-annot',
			  required = True,
			  help = "path to the directory where the stability change files with new ENST annotated\
                  will be located")

@click.option('--log_dir', 
	 		  '-l',
			  required = True,
			  help = "path to the folder to store logs file to record those genes whose ENST is not unique")

# -- Main function  -- #

def stabch_annotate_enst(stabch_dir, symb_colname, mapping_dir, stabch_annot_dir, log_dir):
    """
    Adds to the stability change tables the unique ENST id so that every gene/protein has a 
    unique identifier 
    """

    print("\nLoading mapping file...")
    mapping = pd.read_csv(mapping_dir, sep = "\t")

    print("Detecting repetead SYMBOLS...\n")
    rep_symbs = check_reps(mapping, "symbol")

    # CCLE files
    if os.path.isfile(stabch_dir):   
        print("Starting CCLE pipeline (one file: pancancer)\n")
        
        with open(log_dir, 'w') as log_file:

            stabch_df = pd.read_csv(stabch_dir, sep = "\t", compression = "gzip")

            print(f"\tAnnotating ENSTs")
            stabch_df = stabch_df.progress_apply(lambda row: annot_symb(row, mapping, rep_symbs), axis = 1)

            print(f"\tDiscarding NA genes...")   # those still having NA in Feature column
            reps = stabch_df.loc[stabch_df["Feature"].isna(), symb_colname].value_counts()
            print(f"\t\tDiscarded genes (duplicates or missing in mapping file: {reps}")
            print(reps, file = log_file)
            stabch_df.dropna(subset = ["Feature"], inplace = True)

            stabch_df.to_csv(stabch_annot_dir, sep = "\t", index = False, compression = "gzip")
    
    # TCGA and CPTAC files
    elif os.path.isdir(stabch_dir):
        print("Starting TCGA or CPTAC pipeline (several files: cancer types)\n")

        print("\nDefining cancer types...\n")
        ctypes = os.listdir(stabch_dir)        

        c = 0

        print(f"ENSTs annotation starting for {len(ctypes)} cancer types")
        print("########################################################\n")
        with open(log_dir, 'w') as log_file:

            for ctype in ctypes:

                c += 1

                print(f"\tPreprocessing of {ctype.split('.')[0]}... ({c}/{len(ctypes)})")
                stabch_df = pd.read_csv(stabch_dir+ctype, sep = "\t", compression = "gzip")
                
                # Some ctypes are empty (TCGA for now)
                if stabch_df.empty:
                    print(f'{ctype} removed: empty Dataframe')
                
                else:
                    print(f"\tAnnotating ENSTs")
                    stabch_df = stabch_df.progress_apply(lambda row: annot_symb(row, mapping, rep_symbs), axis = 1)

                    print(f"\tDiscarding NA genes...")   # those still having NA in Feature column
                    reps = stabch_df.loc[stabch_df["Feature"].isna(), symb_colname].value_counts()
                    print(f"\t\tDiscarded genes (duplicates or missing in mapping file): {reps}")
                    print(f"\nCohort: {ctype.split('.')[0]}\n", file = log_file)
                    print(reps, file = log_file)
                    stabch_df.dropna(subset = ["Feature"], inplace = True)

                    stabch_df.to_csv(stabch_annot_dir+ctype, sep = "\t", index = False, compression = "gzip")

    print("\n##################")
    print("Analysis finished")


if __name__ == "__main__":
    stabch_annotate_enst() 