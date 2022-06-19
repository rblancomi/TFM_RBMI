# Execute with qmap: qmap submit stabch_annotate_new_instances.qmap

# -- Import libraries -- # 

import os
import pandas as pd
from tqdm import tqdm
tqdm.pandas()
import click

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

def dict_stabch_new_instance(genes, discov_insts):
    """
    Generate a dictionary of the form {gene: [motif, start, end]} from a discovered instances
    dataframe and stability change dataframe from a specific cancer type.

    Parameters
    ----------
    genes: list
           Unique values for all the genes that have been sequenced and included in the stability 
           change dataframe
    discov_insts: dict
            Dictionary of the form {motif: discovered instances dataframe}
    
    
    Returns
    -------
    stabch_new_instances: dict
            Dictionary of the form {gene: [motif, start, end]}
    """

    stabch_new_instances = {}

    for g in genes:

        # Store all the degrons found in the gene
        g_degrons = []

        for E3 in discov_insts:
            if g in discov_insts[E3]["protein_id"].values:  
                
                # Dataframe subset to contain all the degrons falling in the specific gene
                # (Generally, a protein got more than one degron discovered) 
                gene_discov_insts = discov_insts[E3].loc[discov_insts[E3]["protein_id"] == g].copy()
                gene_discov_insts.reset_index(inplace = True, drop = True)

                # If a protein contains more than one degron, the dict value is a list of lists
                for row in gene_discov_insts.itertuples(index = False):
                    g_degrons.append([E3, row.start, row.end])
                
                if g in stabch_new_instances.keys():
                    stabch_new_instances[g] = stabch_new_instances[g] + g_degrons
                else:
                    stabch_new_instances[g] = g_degrons
                    
    return stabch_new_instances      

def annot_stabch_new_instances(row, d):
    """
    Generates the information of degron_status, E3, degron_start, degron_end to add
    to the stability change df (per gene)

    Parameters
    ----------
    row: Dataframe row
           Row from the stability change dataframe from a specific cancer type. 
           Mandatory columns: Protein_position, Consequence, Phenotype, Feature
    d: dict
            Output dictionary of dict_stabch_new_instance
    
    Returns
    -------
    row: Dataframe row
            Row from the stability change dataframe containing 4 additional columns:
            Degron_status, E3, degron_start, degron_end
    """
    
    enst = row["Feature"]

    degron_vs_mut_loc_l = []
    E3_l = []
    degron_start_l = []
    degron_end_l = []
    
    # There is a degron discovered for this protein
    if enst in d.keys():

        row["Degron"] = True

        # Protein's alteration localization with respect to the degron
        # Obtain mutant position (in protein) and transform to int; some mutants do not have a position or the position is an interval
        try:
            mut_pos = int(row["Protein_position"])
        except:
            if row["Protein_position"] == "-":
                mut_pos = row["Protein_position"]
            else:
                if isinterval(row["Protein_position"]):
                    try:
                        mut_pos = [int(i) for i in row["Protein_position"].split("-")]
                    except:
                        mut_pos = row["Protein_position"]  # for intervals with "?" character
                else:
                    mut_pos = row["Protein_position"]
        
        # Localize the mutations position with respect to the degron(s)
        # Account for proteins with more than one degron
        for inst in d[enst]:
            
            E3 = inst[0]
            degron_start = inst[1]
            degron_end = inst[2]
            
            # No mutation, WT form
            if row.Phenotype == "WT":
                degron_vs_mut_loc_l.append("WT")

            # Mutation in UNKNOWN POSITION or unable to recover it
            elif isinstance(mut_pos, str):
                degron_vs_mut_loc_l.append("unknown")

            # Mutation is an interval
            elif isinstance(mut_pos, list):    
                # Mutation INSIDE the degron
                if (degron_start <= mut_pos[0]) and (degron_end >= mut_pos[1]):
                    degron_vs_mut_loc_l.append("inside")
                    
                # Mutation OUTSIDE and BEFORE the degron
                elif degron_start >= mut_pos[1]:
                    degron_vs_mut_loc_l.append("outside_before")
                
                # Mutation OUTSIDE and AFTER the degron
                elif degron_end <= mut_pos[0]:
                    degron_vs_mut_loc_l.append("outside_after")

            
            # Mutation is not an interval
            else:    
                # Mutation INSIDE the degron
                if (degron_start <= mut_pos) and (degron_end >= mut_pos):
                    degron_vs_mut_loc_l.append("inside")
                    
                # Mutation OUTSIDE and BEFORE the degron
                elif degron_start >= mut_pos:
                    degron_vs_mut_loc_l.append("outside_before")
                
                # Mutation OUTSIDE and AFTER the degron
                elif degron_end <= mut_pos:
                    degron_vs_mut_loc_l.append("outside_after")
                
            E3_l.append(E3)
            degron_start_l.append(degron_start)
            degron_end_l.append(degron_end)
        
        row["Loc_mut_degron"] = degron_vs_mut_loc_l
        row["E3"] = E3_l
        row["degron_start"] = degron_start_l
        row["degron_end"] = degron_end_l

        return row
            
    
    # There is no degron discovered for this protein
    else:
        
        row["Degron"] = False
        row["Loc_mut_degron"] = "NA"
        row["E3"] = "NA"
        row["degron_start"] = "NA"
        row["degron_end"] = "NA"

        return row

# -- Main function options -- #

@click.command()

@click.option('--stabch_dir', 
	 		  '-stabch',
			  required = True,
			  help = "path to the stability change table (with ENSTs and last exon info annotated)")

@click.option('--discov_instances_dir', 
	 		  '-disc',
			  required = True,
			  help = "path to the motifs discovered instances directory (overlapping degrons pooled)")

@click.option('--stabch_annot_dir', 
	 		  '-annot',
			  required = True,
			  help = "path to the stability change table + discovered instances annotated")


# -- Main function  -- #

def stabch_annotate_new_instances(stabch_dir, discov_instances_dir, stabch_annot_dir):
    """
    Annotates in the stability change tables the discovered degrons information: 
    ["Degron", "Loc_mut_degron", "E3", "degron_start", "degron_end"]
    """

    # Discovered instances (all motifs)
    print("Defining motifs discovered instances...")
    motifs = os.listdir(discov_instances_dir)
    print("\tGenerating a dictionary of discovered instances dataframes...\n")
    discov_insts = {}
    for motif in motifs:
        discov_insts[motif.split(".")[0]] = pd.read_csv(discov_instances_dir+motif, sep = "\t")
    
    # Process the stability change table
    print(f'Annotating {stabch_dir.split("/")[-1].split(".")[0]} stability change table...\n')

    stabch_df = pd.read_csv(stabch_dir, sep = "\t", compression = "gzip")

    # Filter proteins having a discovered instance
    print("\tGathering sequenced genes in a list...")
    genes_ensts = stabch_df.Feature.unique()

    print("\tGenerating a dictionary with the sequenced genes having a new instance discovered by any motif...")
    stabch_new_instance = dict_stabch_new_instance(genes_ensts, discov_insts)

    print("\tAnnotating new instances in the stability change dataframe...")
    stabch_df = stabch_df.progress_apply(lambda row: annot_stabch_new_instances(row, stabch_new_instance), axis = 1)

    print("\tRemoving those proteins without a discovered instance reported...")
    stabch_df_f1 = stabch_df.loc[stabch_df.Degron == True]

    print("\tExploding table: one row per protein per discovered instance...")
    stabch_df_exploded = stabch_df_f1.explode(["Loc_mut_degron", "E3", "degron_start", "degron_end"], ignore_index = True)

    print("\tSaving...\n")
    stabch_df_exploded.to_csv(stabch_annot_dir, sep = "\t", index = False, compression = 'gzip')


if __name__ == "__main__":
    stabch_annotate_new_instances()