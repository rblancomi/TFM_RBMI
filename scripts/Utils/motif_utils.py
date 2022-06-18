# -- Import libraries -- # 

import pandas as pd
import logomaker
import matplotlib.pyplot as plt
import random
import sys
import os
from scipy import stats

## my modules ##
sys.path.append("./")    # modules folder
from loads_utils import load_bg_aa

# -- Functions -- #

## -- SCAN (ACTIVITY) -- ##

# For some reason, things do not work if I use the second version of this function
def motif_scan_v1(protein, w_matrix):
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
            
            if protein[pos] == "U" or protein[pos] == "O" or protein[pos] == "X":        # U is selenocysteine; O is pyrrolysine
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
        
    return protein_scores

def motif_scan(sequence, weight_m):
    """
    Scans a protein sequence with a motif's weight matrix and returns a list of sequence-associated
    scores (one per possible subsequence)
    
    Parameters
    ----------
    sequence: str
                Protein sequence
    weight_m: pandas Dataframe
                Dataframe containing a motif's weight matrix, whose rows correspond to motif's positions 
                and columns to aminoacids. Each cell contains the weight associated to each aminoacid in
                each position of the motif
                   
    Returns
    -------
    sequence_scores: list
                List of sequence-associated scores (one per possible subsequence), 
                of length: len(sequence) - len(weight_m)+1
    """

    if len(sequence) == len(weight_m):

        score = 0

        for pos in range(len(sequence)):
            score += weight_m.loc[pos, sequence[pos]] 

        return [score]

    elif len(sequence) > len(weight_m):

        sequence_scores = []

        for pos in range(len(sequence)-(len(weight_m)+1)):                  # iterate through all possible subsequences to be scanned
            subseq_score = 0
            
            for i in range(len(weight_m)):                               # scan the subsequence
                if sequence[pos] == "U" or sequence[pos] == "O" or sequence[pos] == "X":        # U is selenocysteine; O is pyrrolysine
                                                                                                # X is unknown 
                    subseq_score += 0
                    pos += 1
                    
                else:
                    subseq_score += weight_m.loc[i, sequence[pos]]              # find the aa score in the weight matrix
                    pos += 1
                    
            sequence_scores.append(subseq_score)
        
        return sequence_scores

    elif len(sequence) < len(weight_m):
        #print("Error! Sequence length must be equal or greater than motif length")
        return [-1000]  # to avoid entering a NoneType object in another function (negative value is related to the called of this function in positivity_range function)
        

def motif_scan_flex(id, sequence, start, end, weight_m, proteome):
    """
    Allows motif_scan() to perform a flexible scan in case the sequence is
    shorter than the motif. It extends the sequence both to the left and to
    the right as much as possible to meet the motif's length

    Parameters
    ----------
    id: str
                Corresponding protein sequence id (according to proteome keys)
    sequence: str
                Protein sequence
    start: int
                Protein subsequence start in the complete protein sequence
    end: int
                Protein subsequence end in the complete protein sequence
    weight_m: pandas Dataframe
                Dataframe containing a motif's weight matrix, whose rows correspond to motif's positions 
                and columns to aminoacids. Each cell contains the weight associated to each aminoacid in
                each position of the motif
    proteome: dictionary
                Dictionary containing protein IDs as keys and protein sequences as values
    Returns
    -------
    list(max(scores)) OR [] OR motif_scan(sequence, weight_m): list
                List of sequence-associated scores (one per possible subsequence)
    """

    # Sequence shorter than motif
    if len(sequence) < len(weight_m):

        scores = []
        complete_seq = proteome[id]
        diff = len(weight_m) - len(sequence)

        # Complete sequence is equal or longer than motif
        if len(complete_seq) >= len(weight_m):
            # Scan all possible extended subsequences (move to the left)
            for j in range((len(weight_m)-(len(sequence))+1)):  # +1 to consider first scan is performed with the sequence filled to the right, without filling anything to the left

                new_start = start-1-j   # -1 to account for python indexing
                new_end = end+diff-j

                # Finish if the beginning of the complete sequence is reached
                if new_start >= 0:
                    ext_seq = complete_seq[new_start:new_end]
                    scores.append(motif_scan(ext_seq, weight_m))
                else:
                    break
            
            # The subsequence scoring the highest is considered the final score  
            return list(max(scores))

        # Complete sequence is shorter than motif
        else:
            return []

    # Sequence equal or longer than motif
    else:
        return motif_scan(sequence, weight_m)

def motif_positivity_range(weight_m, seqs, proteome, seqs_format = "df", correct = True):
    """
    Calculates the positive range of a motif with the scores of the sequences used to generate such motif.
    It can calculate the scores directly from the subsequences (real scores) or from the whole sequence scan
    (presumed scores) depending on the sequences input (see: Parameters)

    Parameters
    ----------
    weight_m: pandas Dataframe
                Dataframe containing a motif's weight matrix, whose rows correspond to motif's positions 
                and columns to aminoacids. Each cell contains the weight associated to each aminoacid in
                each position of the motif
    seqs: pandas Dataframe or dictionary
                If Dataframe: contains one motif sequence per row. Mandatory to have
                'gene', 'sequence', 'start', 'end' columns. This dataframe is supposed to contain the motif's sequences 
                If dictionary: contains protein IDs as keys and lists of sequence-associated scores (one per possible subsequence)
                as values. This scan is supposed to be the one performed over the motif's sequences
    seqs_format: str (default: "df")
                If "df", calculates from the sequences of a dataframe the score
                If "scan", extracts the highest sequence score (subsequence) 
    correct: boolean (default: True)
                If True, correct the range min to be positive
                If False, do not perform such correction
                   
    Returns
    -------
    (min(max_scores), max(max_scores)): tuple
                Motif's positive range of scores
    """
    max_scores = []

    # Bypass to avoid iterating over the df 'gene' column in case df format is selected
    if seqs_format == "df":
        seqs_ids = seqs["gene"].to_numpy()
    elif seqs_format == "scan":
        seqs_ids = list(seqs.keys())

    for seq in seqs_ids:

        # Calculate the subsequence of interest motif's score directly (real max score) + position
        if seqs_format == "df":

            id = seq
            sequence = seqs.loc[seqs["gene"] == seq, "sequence"].values[0]
            start = seqs.loc[seqs["gene"] == seq, "start"].values[0]
            end = seqs.loc[seqs["gene"] == seq, "end"].values[0]

            max_score = float(str(motif_scan_flex(id, sequence, start, end, weight_m, proteome))[1:-1])

        # Extract the highest motif's score subsequence (presumed max score) + position
        elif seqs_format == "scan":
            max_score = max(seqs[seq])
        
        max_scores.append(max_score)

    # Correct the range to be positive: only keep positive max scores
    if correct:
        max_scores = [s for s in max_scores if (s > 0)]

    return (min(max_scores), max(max_scores))
    
def motif_percent_discovery_own_seqs(seqs_scan, seqs_df, weight_m):
    """
    Returns the percentage of sequences which were aligned to generate a motif which
    can be found by such motif
    
    Parameters
    ----------
    seqs_scan: dictionary
                Contains protein IDs as keys and lists of sequence-associated scores (one per possible subsequence)
                as values. This scan is supposed to be the one performed over the motif's sequences
    seqs_df: pandas Dataframe
                Dataframe containing one motif sequence per row. Mandatory to have
                'gene', 'start' and 'end' columns. This dataframe is supposed to 
                contain the motif's sequences
    weight_m: pandas Dataframe
                Dataframe containing a motif's weight matrix, whose rows correspond to motif's positions 
                and columns to aminoacids. Each cell contains the weight associated to each aminoacid in
                each position of the motif
                   
    Returns
    -------
    matches/len(seqs_df)*100 : float
                Percentage of own sequences the motif can discover
    """

    matches = 0

    for seq in seqs_scan:

        # Real start/end positions
        r_start = (seqs_df.loc[seqs_df["gene"] == seq, "start"].values[0])-1
        r_end = (seqs_df.loc[seqs_df["gene"] == seq, "end"].values[0])
        seq_len = r_end - r_start

        # Presumed start/end positions (note that max value can be repetead)
        p_start = [i for i, score in enumerate(seqs_scan[seq]) if score == max(seqs_scan[seq])]
        p_end = [j+len(weight_m) for j in p_start]

        for p_s, p_e in zip(p_start, p_end):

            # Sequence and motif of same length: a match in the starting positions is enough
            if seq_len == len(weight_m):
                if r_start == p_s:
                    matches += 1
            
            # Sequence longer than motif: presumed positions must be in between real positions
            elif seq_len > len(weight_m):
                if (p_s in range(r_start, r_end)) and (p_e in range(r_start, r_end)):
                        matches += 1 

            # Sequence shorter than motif: presumed positions must contain real positions
            elif seq_len < len(weight_m):
                if (r_start in range(p_s, p_e)) and (r_end in range(p_s, p_e)):
                        matches += 1 

    return matches/len(seqs_df)*100    

def motif_n_and_percent_positive_hits(seqs_scan, positivity_range):
    """
    Returns the percentage of proteome sequences which are scored positive by a motif
    
    Parameters
    ----------
    seqs_scan: dictionary
                Contains protein IDs as keys and lists of sequence-associated scores (one per possible subsequence)
                as values. This scan is supposed to be the one performed over the non-motif's sequences
    positivity_range: tuple
                Motif's positive range of scores              
    Returns
    -------
    n_positive_hits/total_proteins*100 : float
                Percentage of protein sequences falling in the motif's positivity range
    """

    n_positive_hits = 0
    total_proteins = 0  # to have the exact number of proteins with scores 
                        # (as we are avoiding empty lists)

    for seq in seqs_scan:
        if seqs_scan[seq] != []:  # bypass proteins without scores 

            total_proteins += 1

            # Extract the potential hit in the sequence, which corresponds to the max score
            max_score = max(seqs_scan[seq])

            # Determine if it is a hit by assessing whether it falls in the positivity range
            if max_score >= positivity_range[0]:
                n_positive_hits += 1
    
    return (n_positive_hits, n_positive_hits/total_proteins*100)

def motif_n_positive_hits_subseq(seqs_scan, positivity_range):
    """
    Returns the number of proteome subsequences which are scored positive by a motif
    
    Parameters
    ----------
    seqs_scan: dictionary
                Contains protein IDs as keys and lists of sequence-associated scores (one per possible subsequence)
                as values. This scan is supposed to be the one performed over the non-motif's sequences
    positivity_range: tuple
                Motif's positive range of scores              
    Returns
    -------
    n_positive_hits_subseq: int
                Number of protein subsequences falling in the motif's positivity range
    """

    n_positive_hits = 0

    for seq in seqs_scan:
        if seqs_scan[seq] != []:  # bypass proteins without scores 

            # Extract all the sequence hits, which are those higher or equal than the positivity range
            n_positive_hits += len([subseq_score for subseq_score in seqs_scan[seq] if (subseq_score >= positivity_range[0])])
    
    return n_positive_hits
            
def motif_percent_random_positive_hits(seqs_scan, seqs_df, positivity_range, min_sampling = 20, random_set = None):
    """
    Returns the percentage of random proteome sequences which are scored positive by a motif
    
    Parameters
    ----------
    seqs_scan: dictionary
                Contains protein IDs as keys and lists of sequence-associated scores (one per possible subsequence)
                as values. This scan is supposed to be the one performed over the non-motif's sequences
    seqs_df: pandas Dataframe
                Dataframe containing one motif sequence per row
    positivity_range: tuple
                Motif's positive range of scores  
    min_sampling: int (default: 20)
                Minimum number of sample subsequences
    random_set: list (default: None)
                In case of provided, list of subsequences scores            
    Returns
    -------
    n_positive_hits/total_proteins*100 : float
                Percentage of protein sequences falling in the motif's positivity range
    """
    random_seqs_scores = []

    # Gather a random set of subsequences
    if random_set == None:

        if motif_n_sequences(seqs_df) < min_sampling:
            nseqs = min_sampling
        else:
            nseqs = motif_n_sequences(seqs_df)

        sample = random.sample(list(seqs_scan.keys()), nseqs)

        for seq in sample:
            if seqs_scan[seq] != []:
                random_seqs_scores.append(random.choice(seqs_scan[seq]))

    # In case the random set is predefined
    else:

        random_seqs_scores = random_set

    # Percentage of random subsequences scoring in the positivity range
    random_positives = 0
    for score in random_seqs_scores:
        if score >= positivity_range[0]:
            random_positives += 1
                
    return random_positives/len(random_seqs_scores)*100

def motif_positive_hits_per_sequence(motif_id, seqs_scan, positivity_range):
    """
    Returns the percentage of random proteome sequences which are scored positive by a motif
    
    Parameters
    ----------
    motif_id: str
                Motif's ID
    seqs_scan: dictionary
                Contains protein IDs as keys and lists of sequence-associated scores (one per possible subsequence)
                as values. This scan is supposed to be the one performed over the non-motif's sequences
    positivity_range: tuple
                Motif's positive range of scores          
    Returns
    -------
   n_hits_per_protein_df : pandas Dataframe
                Dataframe with sequences ID as index and number of positive hits per sequence in 
                motif_id column
    """
    n_hits_per_protein = {}

    for seq in seqs_scan:

        if seqs_scan[seq] != []:

            n_hits = len([score for score in seqs_scan[seq] 
                            if score >= positivity_range[0]])
            
            n_hits_per_protein[seq] = n_hits
            

    n_hits_per_protein_df = pd.DataFrame.from_dict(n_hits_per_protein, 
                                                    orient = 'index', columns = [motif_id])
        
    return n_hits_per_protein_df

def motif_positive_hits_vs_len_corr(len_df, hits_df):
    """
    Returns the percentage of random proteome sequences which are scored positive by a motif
    
    Parameters
    ----------
    len_df: pandas Dataframe
                Dataframe with sequences ID as index and sequence length in length column
    hits_df: pandas Dataframe
                Dataframe with sequences ID as index and number of positive hits per sequence in 
                motif_id column     
    Returns
    -------
    r_value : float
                Linear correlation coefficient between a sequence length and the number of positive hits
    """
    # Concat both df by index and avoid proteins not having scores
    len_hits_df = pd.concat([len_df, hits_df], axis = 1, join = "inner")

    # Calculate a linear least-squares regression
    cols = len_hits_df.columns
    slope, intercept, r_value, p_value, stderr = stats.linregress(x = len_hits_df.iloc[:, -1].values, y = len_hits_df.iloc[:, 0].values)

    return r_value

## -- LOGOS -- ##

def motif_logo(seqs_dir, aa, logo_dir, font_name = 'Calibri', color_scheme = "skylign_protein", format = "dataframe"):
    """
    Generates a motif's logo using logomaker.Logo() function from a list of sequences stored in a 
    dataframe under a column named 'sequence'. Sequences must have the same length. For a logical 
    motif logo to be created from these sequences, they should have been previously aligned and 
    gaps removed
    
    Parameters
    ----------
    seqs_dir: str
                Path to the folder containing the sequences in a dataframe or in fasta format
    aa: numpy.ndarray
                1 dimensional numpy array containing the one-letter characters of the 20 aminoacids
                alphabetically ordered
    logo_dir: str
                Path to the folder in which the motif's logo will be located
    format: str (default: dataframe)
                Format of the sequences. Can be dataframe or fasta.

                   
    Returns
    -------
    None

    """
    # Gather sequences in a list (logomaker requirement)
    if format == "dataframe":
        seqs_df = pd.read_csv(seqs_dir, sep = "\t")
        seqs = seqs_df['sequence'].tolist()
    elif format == "fasta":
        seqs = []
        with open(seqs_dir, 'r') as f:
            for line in f:
                if line[0] != ">":
                    seqs.append(line.strip())
    
    # Generate a count matrix with the sequences to add those aa not present in any sequence (for next step is necessary)
    count_m = logomaker.alignment_to_matrix(seqs, to_type = "counts")
    ## In the case the number of aminoacids is < 20, a value of 0.0 is added to the count matrix
    if count_m.shape[1] != 20:
        
        diff_aa = set(aa) - set(count_m.columns)
        for d in diff_aa:
            count_m[d] = 0.0

        ## Sort aminoacids alphabetically (don't remember why right now...)    
        count_m.sort_index(axis = 1, inplace = True)
    
    # Generate logo
    logo = logomaker.Logo(count_m, font_name = font_name, color_scheme = color_scheme)
    plt.savefig(logo_dir, dpi = 800, transparent = False)

## -- MOTIF STRUCTURE -- ##

def motif_length(weight_m):
    """
    Calculates a motif's length from its weight matrix length
    
    Parameters
    ----------
    weight_m: pandas dataframe
                Dataframe containing a motif's weight matrix, whose rows correspond to motif's positions 
                and columns to aminoacids. Each cell contains the weight associated to each aminoacid in
                each position of the motif
                   
    Returns
    -------
    len(weight_m): int
                Motif's length

    """

    return len(weight_m)

def motif_info_content_len_norm(prob_m, aa_probs, motif_n_sequences):
    """
    Calculates a motif's informational content normalized by the motif's length.
    Update! Also normalize by the number of sequences in the MSA, as logomaker performs
    a pseudocount correction that biases the results towards higher informational content
    if the MSA is larger
    
    Parameters
    ----------
    prob_m: pandas dataframe
                Dataframe containing a motif's probability matrix, whose rows correspond to motif's positions 
                and columns to aminoacids. Each cell contains the probability associated to each aminoacid in
                each position of the motif
    aa_probs: numpy.ndarray
                1 dimensional numpy array containing the background probabilities of the 20 aminoacids (in 
                alphabetical order according to their one-letter characters)
    motif_n_sequences: int
                Motif's number of sequences (MSA length)
                   
    Returns
    -------
    (info_m.to_numpy().sum()) / (len(prob_m)): int
                Motif's informational content normalized by the motif's length

    """

    info_m = logomaker.transform_matrix(prob_m, from_type = 'probability', 
                                        to_type = 'information', 
                                        background = aa_probs)   # does not accept weight matrix

    return ((info_m.to_numpy().sum()) / (len(prob_m))) / motif_n_sequences

def motif_n_sequences(seqs_df):
    """
    Returns the number of sequences which were aligned to generate a motif
    
    Parameters
    ----------
    seqs_df: pandas Dataframe
                Dataframe containing one motif sequence per row
                   
    Returns
    -------
    len(msa_df): int
                Motif's number of sequences

    """

    return len(seqs_df)

## -- MOTIF MATRICES TRANSFORMATIONS -- ##

def motif_from_w_matrix_to_prob_m(weight_m_dir, aa_bg_dir, prob_m_dir):
    """
    Transforms a motif's position weight matrix (PWM) into a position
    probability matrix (PPM)

    Parameters
    ----------
    weight_m_dir: str
                Path to the folder containing the PWMs
    aa_bg_dir: str
                Path to the folder containing the aminoacids background probabilities file
    prob_m_dir: str
                Path to the folder that will store the PPMs
                   
    Returns
    -------
    None
    """

    # Load motifs
    motifs_IDs = [motif.split(".")[0] for motif in os.listdir(weight_m_dir)]

    # Load aa background probabilities for matrix transformation
    aa_probs, aa_names = load_bg_aa(aa_bg_dir)

    for motif in motifs_IDs:

        # Load motif's PWM
        weight_m = pd.read_csv(weight_m_dir+motif+".tsv", sep = "\t")

        # Transform in PPM
        prob_m = weight_m = logomaker.transform_matrix(weight_m, from_type = 'weight', 
                                          to_type = 'probability', background = aa_probs)
        
        # Save
        prob_m.to_csv(prob_m_dir+motif+".tsv", sep = "\t", header = True, index = False)