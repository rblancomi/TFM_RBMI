# -- Import libraries -- # 

import os
import sys
import pandas as pd
import numpy as np
from tqdm import tqdm              
import logomaker
import json
import random
import gzip
import click

## my modules ##
sys.path.append("../Utils/")    # modules folder
from fasta_utils import readFasta_gzip
from motif_utils import *
from loads_utils import load_bg_aa


# -- Main function options -- #

@click.command()

@click.option('--aa_bg_dir', 
	 		  '-aa_bg',
			  required = True,
			  help = "path to the folder where the aminoacids background probabilities\
              file is located")

@click.option('--proteome_dir', 
	 		  '-p',
			  required = True,
			  help = "path to the folder where the proteome file is located")

@click.option('--prob_m_dir', 
			  '-pm',
			  required = True,
			  help = "path to the folder containing motifs position probability matrices")

@click.option('--weight_m_dir', 
			  '-wm',
			  required = True,
			  help = "path to the folder containing motifs position weight matrices")

@click.option('--pos_set_dir', 
	 		  '-ps',
			  required = True,
			  help = "path to the folder where positive sequences (real degrons) files are stored")

@click.option('--scan_dir', 
	 		  '-s',
			  required = True,
			  help = "path to the folder where motifs proteome-wide scans are stored")

@click.option('--scan_own_seqs_dir', 
	 		  '-s_pos',
			  required = True,
			  help = "path to the folder where positive sequences (real degrons) scans are stored")

@click.option('--scan_rest_seqs_dir', 
	 		  '-s_nopos',
			  required = True,
			  help = "path to the folder where the non-positive sequences scans are stored")

@click.option('--motifs_metrics_dir', 
	 		  '-metrics',
			  required = True,
			  help = "path to the motifs metrics file")

@click.option('--hits_per_sequence_vs_len_dir', 
	 		  '-hitsvslen',
			  required = True,
			  help = "path to the hits per sequence vs length file")

@click.option('--proteome_format', 
	 		  '-pf',
			  required = False,
			  default = "fasta_gzip",
			  help = "proteome format, which can be json dictionary or fasta compressed")



# -- Main function  -- #

# Execute in an interactive -m 64 (64 Gb memory)

def motif_quality_analysis(aa_bg_dir, proteome_dir, prob_m_dir, weight_m_dir, pos_set_dir,
                            scan_dir, scan_own_seqs_dir, scan_rest_seqs_dir,
                            motifs_metrics_dir, hits_per_sequence_vs_len_dir, proteome_format):

    print("\n## PERFORM GENERAL LOADINGS AND COMPUTATIONS ##\n")
    print("Load aminoacid background probabilities\n")
    aa_probs, aa_names = load_bg_aa(aa_bg_dir)

    print("Load proteome")
    if proteome_format == "json":
        with open(proteome_dir) as fp:
            proteome = json.load(fp)
    elif proteome_format == "fasta_gzip":
        proteome = readFasta_gzip(proteome_dir)
    else:
        return print("Incorrect proteome format, only JSON or fasta compressed")

    print("Calculate proteome proteins length and store this data in a dataframe for future correlation analysis with the number of hits identified in each protein per motif")
    protein_lengths = {}
    for protein in proteome:
        protein_lengths[protein] = len(proteome[protein])
    
    hits_per_sequence_vs_len = pd.DataFrame.from_dict(protein_lengths, orient = 'index', columns=['Protein length'])
    len_df = hits_per_sequence_vs_len.copy() # to be used with motif_positive_hits_vs_len_corr()
    
    print("Load motifs IDs\n")
    motifs_IDs = [motif.split(".")[0] for motif in os.listdir(weight_m_dir)]

    print("Defining structural and activity metrics to be computed (lists)\n")

    # Structural
    motif_lengths = []
    inform_contents = []
    number_seqs_msas = []

    # Activity
    positivity_ranges = [] 
    percent_discovery_own_seqs = []
    n_positive_hits = []
    percent_positive_hits = []
    n_positive_hits_subseq = []
    percent_random_positive_hits = []
    corr_hits_vs_len = []

    print("## ANALYSIS STARTING ##\n")

    counter = 0

    for motif in motifs_IDs:

        counter += 1
        print(f"Computing {motif} metrics... ({counter}/{len(motifs_IDs)})\n")

        # Load files
        
        print("\tLoad probability matrix")
        prob_m = pd.read_csv(prob_m_dir+motif+".tsv", sep = "\t")
        
        print("\tLoad weight matrix")
        weight_m = pd.read_csv(weight_m_dir+motif+".tsv", sep = "\t")

        print("\tLoad MSA matrix: motif sequences (positive sequences)")
        seqs_df = pd.read_csv(pos_set_dir+motif+".tsv", sep = "\t")

        print("\tLoad motif scan (dictionary)")
        with open(scan_dir+motif+".json") as fp:
            scan = json.load(fp)
        
        print("\tLoad motif own sequences scan (dictionary)")
        with open(scan_own_seqs_dir+motif+".json") as fp:
            scan_own_seqs = json.load(fp)

        print("\tLoad motif rest-of-sequences scan (dictionary)\n")
        with open(scan_rest_seqs_dir+motif+".json") as fp:
            scan_rest_seqs = json.load(fp)


        # Compute structural metrics

        print("\tCalculating structural metrics...")

        print("\tMotif length")
        motif_lengths.append(motif_length(weight_m))

        print("\tNumber of sequences which were aligned to generate this motif\n")
        number_seqs_msa = motif_n_sequences(seqs_df)
        number_seqs_msas.append(number_seqs_msa)

        print("\tMotif informational content (length normalized)")
        inform_contents.append(motif_info_content_len_norm(prob_m, aa_probs, number_seqs_msa))

        # Compute activity metrics

        print("\tCalculating activity metrics...")

        print("\tPositivity range")
        positivity_range = motif_positivity_range(weight_m, seqs_df, proteome, seqs_format = "df", correct = True)
        positivity_ranges.append(positivity_range)

        print("\tPercent of sequence which were aligned to generate this motif that are discovered by it")
        percent_discovery_own_seqs.append(motif_percent_discovery_own_seqs(scan_own_seqs, seqs_df, weight_m))

        print("\tPercentage of positive hits in the proteome")
        n, percent = motif_n_and_percent_positive_hits(scan, positivity_range)
        n_positive_hits.append(n)
        percent_positive_hits.append(percent)

        print("\tNumber of positive subsequence hits in the proteome")
        n_positive_hits_subseq.append(motif_n_positive_hits_subseq(scan, positivity_range))
        
        print("\tPercentage of random positive hits in the proteome")
        percent_random_positive_hits.append(motif_percent_random_positive_hits(scan_rest_seqs, seqs_df, positivity_range, 
        min_sampling = 20, random_set = None))

        print("\tCompute hits per sequence to add to correlation dataframes...")
        hits_per_sequence = motif_positive_hits_per_sequence(motif, scan, positivity_range)
        hits_per_sequence_vs_len[motif] =  hits_per_sequence # protein names as indeces are all present because the starting df is the protein's length one, which does not skip any protein

        print("\tCalculate the linear correlation between the number of hits per sequence and its length\n")
        corr_hits_vs_len.append(motif_positive_hits_vs_len_corr(len_df, hits_per_sequence))

    print("Saving motifs metrics as csv...\n")
    motifs_metrics_df = pd.DataFrame(list(zip(motifs_IDs, 
                                            motif_lengths,
                                            inform_contents,
                                            number_seqs_msas,
                                            positivity_ranges,
                                            percent_discovery_own_seqs,
                                            n_positive_hits,
                                            n_positive_hits_subseq,
                                            percent_positive_hits,
                                            percent_random_positive_hits,
                                            corr_hits_vs_len)),
                                            columns = ["motif_id", 
                                            "motif_length",
                                            "motif_info_content",
                                            "n_msa_sequences", 
                                            "positivity_range",
                                            "percent_discovery_msa_seqs",
                                            "number_positive_hits",
                                            "number_positive_hits_subseq",
                                            "percent_positive_hits",
                                            "percent_random_positive_hits",
                                            "correlation_hits_vs_length_sequence"])

    motifs_metrics_df.to_csv(motifs_metrics_dir, sep = "\t", index = False, compression = 'gzip')

    print("Saving hits-per-sequence-vs-length tables as csv...\n")
    hits_per_sequence_vs_len.to_csv(hits_per_sequence_vs_len_dir, sep = "\t", compression = 'gzip')

    print("ANALYSIS FINISHED")

if __name__ == "__main__":
    motif_quality_analysis()






    







