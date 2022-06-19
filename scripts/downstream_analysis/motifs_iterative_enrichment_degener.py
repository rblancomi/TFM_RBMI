# -- Import libraries -- # 

import os
import pandas as pd
import numpy as np
from tqdm import tqdm 
import gzip
import logomaker
from  itertools import chain
import click              # note:click only works with lower case arguments?



# -- Auxiliary functions -- # 

def readFasta_gzip(file):
    """ 
    Reads all sequences of a compressed FASTA file and returns a dictionary with the protein's accession number
    as key and the sequence as value
    
    Parameters
    ----------
    file: str
           Path to the folder where the compressed FASTA file is located
    
    
    Returns
    -------
    d: dict
            Dictionary containing each protein's accession number as key and they sequence as value
    """  
    
    with gzip.open(file, 'rt') as f:         # important to set read mode as 'rt' because we are reading a compressed text file
        
        data = f.readlines()
        d = {} 
        seq_values = []
        
        count_id = 0      # Flag variable
        
        for line in data:
            
            if count_id == 0:  # Flag in cero indicates this line is the first sequence description line in the fasta file
                if line[0] == '>': # Description lines start with “>” character
                    id_key = line.split()[0][1:].split("|")[1] # The first column of the description line is the sequence id
                                                               # and the second element of the description is the AC             
                    count_id = 1  # Change flag to 1 to indicate the following sequences aren't the first
            
            else:
                if line[0] == '>':
                    d[id_key] = ''.join(seq_values) # Adds a new entry to the dictionary before overwritting id_key
                    id_key = line.split()[0][1:].split("|")[1]
                    seq_values = [] # Restarts the variable that stores every sequence
                else:
                    seq_values.append(line.strip('\n'))  # Sequences are divided in several lines in the fasta file
        
        d[id_key] = ''.join(seq_values) # Adds last dictionary entry
    
    print(f'\tNumber of retrieved sequences: {len(d)}')
    return d

# Modified version to account for degeneration
def motif_scan(protein, w_matrix, deg_level = 0):
    """
    Scans a protein's sequence with a motif's weight matrix and returns a list of protein-associated
    matching scores
    
    Parameters
    ----------
    protein: str
                Protein's sequence
    w_matrix: pandas dataframe
                Dataframe containing a motif's weight matrix, whose rows correspond to motif positions 
                and columns to aminoacids. Each cell contains the weight associated to each aminoacid in
                every position of the motif
                   
    Returns
    -------
    substrate scores: list
                List of protein-associated matching scores, of length: len(protein) - len(w_matrix)
    """
    protein_scores = list()             # list to store scores per protein
    #("Analyzing new protein!")
    
    for pos in range(len(protein)-len(w_matrix)+1):                  # scanning the protein avoiding out of index error
        
        window_score = 0
        
        for i in range(len(w_matrix)):                               # score the specific window
            #print(f"Analyzing motif's position: {i}")
            
            if deg_level != 0:
                deg_positions = list(chain.from_iterable((j, len(w_matrix)-1-j) for j in range(deg_level)))
                if i in deg_positions:
                    window_score += w_matrix.iloc[[i]].values.max()
                    pos += 1

            elif protein[pos] == "U" or protein[pos] == "O" or protein[pos] == "X":        # U is selenocysteine; O is pyrrolysine
                                                                                                                # X is unknown 
                window_score += 0
                #print(f"Analyzing substrate's position: {pos}")        
                pos += 1

            else:
                pos_score = w_matrix.loc[i, protein[pos]]              # find the aa puntuation in the weight matrix
                window_score += pos_score
                #print(f"Analyzing substrate's position: {pos}")        
                pos += 1
                
        protein_scores.append(window_score)
        
    if len(protein_scores) == 1:
        for ps in protein_scores:
            score = float(ps)
        return score
    else:
        return protein_scores

def refine_alignment(weight_m, degrons, substrates, proteome, aa, aa_probs, iter, deg_level = 0, enrichment = False):
    """
    """
    # Iteration number
    iter += 1
    print(f'\n### ITERATION {iter} ###\n')
    
    # Positivity score range for the substrates used to generate the motif  
    scores_degrons = []
    
    for idx, row in degrons.iterrows():
        
        scores_degrons.append(motif_scan(row.sequence, weight_m))
    
    positive_range_degrons = (min(scores_degrons), max(scores_degrons))
    if positive_range_degrons[0] < 0:
        #positive_range_degrons = (0, max(scores_degrons))
        positive_range_degrons = (min([s for s in scores_degrons if (s >= 0)]), max(scores_degrons))
    
    print(f'\tPositivity scores range: {positive_range_degrons}\n')

    # Scan each substrate with the PWM to find degrons
    new_degrons = []
    
    for substrate in substrates:
        #print(f'\tSubstrate being scanned: {substrate}')

        if deg_level == 0:
            sequence = proteome[substrate]
            scores = motif_scan(sequence, weight_m)
        else:
            sequence = proteome[substrate]
            scores = motif_scan(sequence, weight_m, deg_level)

        pos_scores = [s for s in scores if (s >= positive_range_degrons[0]) and (s <= positive_range_degrons[1])].sort(reverse = True) # sort to first analyze higher values
        #print(f'\t\tPositive scores: {pos_scores}')
        
        # In case there is a sequence scoring positive (presumed degron), extract its position
        if pos_scores != []:
            idx_pos_scores = [idx for idx, s in enumerate(scores) 
                              if (s >= positive_range_degrons[0]) and (s <= positive_range_degrons[1])]
            for start in idx_pos_scores:
                end = start+len(weight_m)
                degron = proteome[substrate][start:end]
                
                # Check whether this degron was already considered in the motif
                ## This degron was not considered in the motif and the substrate has not contributed to the motif more than 1 time
                if (degrons[(degrons.substrate == substrate) & (degrons.sequence == degron)].empty == True) and\
                    (degrons[degrons.substrate == substrate].shape[0] < 1) and (substrate not in chain(*new_degrons)):
                    #print("\t\tNew degron instance found! Adding to the alignment...\n")
                    new_degrons.append([substrate, degron, start, end, iter])
                    #print(f'{degrons[degrons.substrate == substrate].shape[0]}\n')
                    #print(f'{degrons[degrons.substrate == substrate]}\n')
                    #print("substrate\tsequence\tstart\tend\titeration")
                    #print(f'\t{substrate}\t{degron}\t{start}\t{end}\t{iter}\n')

                # This degron was already considered in the motif
                else:
                    ## Degron's position is unknown
                    if degrons[(degrons.substrate == substrate) & (degrons.sequence == degron) &\
                        (degrons.start == "NA") & (degrons.end == "NA")].empty == False:
                        #print("\tThis degron instance was already considered in the alignment but coordinates were not annotated. Annotating coordinates...\n")
                        degrons.at[(degrons["substrate"] == substrate), 'start'] = start
                        degrons.at[(degrons["substrate"] == substrate), 'end'] = end
    
    # End iterative process (no new degron added in this iteration)
    if new_degrons == []:
        if enrichment or ((deg_level*2) > (len(weight_m)-3)):
            print("## No new degrons were found. Ending iterative process ##\n")
            return weight_m, degrons, iter, deg_level
        else:
            print("## No new degrons were found but no enrichment ocurred. Forzing motif outter degeneration... ##\n")
            deg_level += 1
            return refine_alignment(weight_m, degrons, substrates, proteome, aa, aa_probs, iter, deg_level)
    
    # If new degrons are found, generate new alignment and weight matrix
    else:
        enrichment = True

        print("New degrons were found. Including them in the alignment...\n")
        new_degrons_df = pd.DataFrame(new_degrons, columns = ['substrate', 'sequence', 'start', 'end', 'iteration'])
    
        # Join with previous degrons df
        refined_alignment = pd.concat([degrons, new_degrons_df])
        refined_alignment.reset_index(inplace = True, drop = True)
        #print(f'{refined_alignment}\n')
        
        # Generate new count matrix from old and new degrons
        count_m = logomaker.alignment_to_matrix(list(refined_alignment.sequence), to_type = "counts")
    
        if count_m.shape[1] != 20:
        
            diff_aa = set(aa) - set(count_m.columns)
            for d in diff_aa:
                count_m[d] = 0.0

            count_m.sort_index(axis = 1, inplace = True)
    
        # Generate weight matrix considering bg aa probabilities
        weight_m = logomaker.transform_matrix(count_m, from_type = 'counts', 
                                              to_type = 'weight', background = aa_probs)
        #print(f'New PWM generated: \n{weight_m}\n')
        
        # Perform the scanner again with the refined weight matrix (only until third iteration)
        if iter < 3:
            return refine_alignment(weight_m, refined_alignment, substrates, proteome, aa, aa_probs, iter, deg_level, enrichment)
        
        else:
            print(f'## Max number of iterations ({iter}) met. Ending iterative process ##\n')
            return weight_m, refined_alignment, iter, deg_level



# -- Main function options -- #

# note: to use it with individual E3 ligases, generate temp folders for pre-enrichment alignments and PWMs
# note2: to use it with E3 ligases whose name is not UniProt AC, add AC at the beginning of the name, separated from the rest with an underscore

@click.command()

@click.option('--esis_dir', 
	 		  '-esis',
			  required = True,
			  help = "path to the folder where the human E3-substrate interactions file is located")

@click.option('--proteome_dir', 
	 		  '-p',
			  required = True,
			  help = "path to the folder where the proteome file is located")

@click.option('--aa_bg_dir', 
	 		  '-aa',
			  required = True,
			  help = "path to the folder where the aminoacids background probabilities\
              file is located")

@click.option('--alignments_dir', 
			  '-align',
			  required = True,
			  help = "path to the folder containing degrons alignments")

@click.option('--weight_m_dir', 
			  '-wm',
			  required = True,
			  help = "path to the folder containing motifs weight matrices")

@click.option('--refined_alignments_dir', 
			  '-ralign',
			  required = True,
			  help = "path to the folder in which the refined degrons alignments will be located")

@click.option('--refined_weight_m_dir', 
			  '-rwm',
			  required = True,
			  help = "path to the folder in which the refined motifs weight matrices will be located")

@click.option('--log_file', 
			  '-log',
			  required = True,
			  help = "file to store enrichment progress")



# -- Main function  -- #


def refine_motif(esis_dir, proteome_dir, aa_bg_dir, alignments_dir, weight_m_dir, refined_alignments_dir, refined_weight_m_dir, log_file):

    # Data loading
    print("\nESIs file being loaded...\n")
    ESIs = pd.read_csv(esis_dir, sep = "\t")
    
    print("Proteome being loaded...\n")
    proteome = readFasta_gzip(proteome_dir)

    print("\nLoading aminoacid background probabilities...\n")
    bg_matrix = pd.read_table(aa_bg_dir).sort_values(by = "Aminoacid")
    aa_probs = bg_matrix["Frequency"].to_numpy()
    aa = bg_matrix["Aminoacid"].to_numpy() 

    # E3 ligases for which we have  motifs generated 
    #E3s = list(set([E3.split(".")[0] for E3 in os.listdir(alignments_dir)]))
    E3s = np.unique(np.array([E3.split(".")[0] for E3 in os.listdir(alignments_dir)]))
    counter = 0

    for E3 in E3s:
        counter += 1
        print("===============================================")
        print(f'E3 ligase PWM being enriched: {E3} ({counter}/{len(E3s)})')
        print("===============================================\n")
        
        # Load PWM
        weight_m = pd.read_csv(weight_m_dir+E3+".tsv", sep = "\t")
        
        # Load degrons (1st alignment, not refined)
        degrons = []
        seqs = []
        substrates = []
        with open(alignments_dir+E3+".fasta", 'r') as f:
            for line in f:
                if line[0] != ">":
                    seqs.append(line.strip())
                else:
                    substrates.append(line.strip()[1:].split("_")[0])
        
        for substrate, seq in zip(substrates, seqs):
            degrons.append([substrate, seq, "NA", "NA", 0])
        
        degrons_df = pd.DataFrame(degrons, columns = ['substrate', 'sequence', 'start', 'end', 'iteration'])
        n_original_degrons = len(degrons_df)
        #print(f'Starting degrons which generated first PWD:\n{degrons_df}\n')
            
        # ESIs subset for the specific E3 ligase
        E3_ESIs = set(ESIs[ESIs.E3_AC == E3.split("_")[0]].SUB_AC)
        
        # Iterative enrichment
        #print(refine_alignment(weight_m, degrons_df, E3_ESIs, proteome, aa, aa_probs, 0))
        weight_m, degrons_df, iter, deg_level = refine_alignment(weight_m, degrons_df, E3_ESIs, proteome, aa, aa_probs, 0)
        print("\n\nENRICHMENT ENDED\n")
        #print(f'\nFinal Refined PWM:\n\n{weight_m}')
        if len(degrons_df)-n_original_degrons == 0:
            print("No new sequences were found, the motif was not enriched\n")
        else:
            print("Refined PWM being saved...\n")
            weight_m.to_csv(refined_weight_m_dir+E3+".tsv", sep = "\t", index = False)
            #print(f'\nStarting degrons plus added degrons:\n\n{degrons_df}\n')
            print("Refined alignment being saved...\n")
            degrons_df.to_csv(refined_alignments_dir+E3+".tsv", sep = "\t", index = False)

        with open(log_file, "a") as f:
            print(f'{counter}. E3 ligase {E3}. Refined with {iter} iterations. {len(degrons_df)-n_original_degrons} new sequences added (degeneration level: {deg_level}). Number of UbiNet substrates: {len(E3_ESIs)}\n', 
            file = f)

if __name__ == "__main__":
    refine_motif()

