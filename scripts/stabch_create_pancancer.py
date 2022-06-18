# Execute in an interactive -m 64

# -- Import libraries -- # 

import os
import pandas as pd
from tqdm import tqdm
tqdm.pandas()
import click

# -- Auxiliary functions -- # 

def annot_ctype(row, ctype):
    
    row["Cancer_type"] = ctype.split(".")[0]
    
    return row

# -- Main function options -- #

@click.command()

@click.option('--stabch_dir', 
	 		  '-stabch',
			  required = True,
			  help = "path to the stability change files directory (generally,\
                   with discovered instances annotated")

@click.option('--pancancer_dir', 
	 		  '-panc',
			  required = True,
			  help = "path to the mutations pancancer file")

# -- Main function  -- #

def stabch_create_pancancer(stabch_dir, pancancer_dir):

    print(f'\nDefining cancer types...\n')
    ctypes = os.listdir(stabch_dir)

    print(f'Annotating cancer type to each cancer type table:\n')
    ctypes_dfs = {}
    c = 0

    for ctype in ctypes:

        c += 1

        print(f'\t{ctype.split(".")[0]} ({c}/{len(ctypes)})')
        ctype_df = pd.read_csv(stabch_dir+ctype, sep = "\t", compression = "gzip")
        ctype_df = ctype_df.progress_apply(lambda row: annot_ctype(row, ctype), 
        axis = 1)

        ctypes_dfs[ctype] = ctype_df
    
    print("Generating pancancer table...")
    pancancer = pd.concat(ctypes_dfs.values(), ignore_index = True)

    print("Saving pancancer table...")
    pancancer.to_csv(pancancer_dir, sep = "\t", index = False, 
    compression = "gzip")

if __name__ == "__main__":
    stabch_create_pancancer()