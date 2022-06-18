# -- Import libraries -- # 

import pandas as pd
import os
import gzip                       

# -- Functions -- #

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
    
    print(f'Number of retrieved sequences: {len(d)}\n')
    return d

# the version below also avoids * signaling for stop codons when needed
def readFasta_header_gzip(file):
    """ 
    Reads all sequences of a compressed FASTA file and returns a dictionary with the protein's complete header
    as key and the sequence as value. This function's version also implements the avoidance of last-character
    asterisks (*) that indicate STOP codons
    
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
                    id_key = line.strip('\n')[1:]  # avoid ">"          
                    count_id = 1  # Change flag to 1 to indicate the following sequences aren't the first
            
            else:
                if line[0] == '>':
                    if ''.join(seq_values)[-1] == "*":   # avoid asterisk at the end of the seq if present
                        d[id_key] = ''.join(seq_values)[:-1] # Adds a new entry to the dictionary before overwritting id_key
                    else:
                        d[id_key] = ''.join(seq_values)

                    id_key = line.strip('\n')[1:]
                    seq_values = [] # Restarts the variable that stores every sequence
                else:
                    seq_values.append(line.strip('\n'))  # Sequences are divided in several lines in the fasta file
        
        if ''.join(seq_values)[-1] == "*":
            d[id_key] = ''.join(seq_values)[:-1] # Adds a new entry to the dictionary before overwritting id_key
        else:
            d[id_key] = ''.join(seq_values)
    
    print(f'Number of retrieved sequences: {len(d)}\n')
    return d


def from_iterationsDf_to_fasta(df_dir, fasta_dir):

    """ 
    Reads dataframes containing sequence information in the form "substrate, sequence, start, end, iteration" 
    and generates a fasta file with substrate AC and sequence. Functions input is generated from refine_motif.py
    or refine_motif_degener.py functions.
    
    Parameters
    ----------
    df_dir: str
           Path to the folder where the dataframes with sequence information \
               (substrate, sequence, start, end, iteration) is located
    fasta_dir: str
           Path to the folder where sequences fasta files will be located
    
    Returns
    -------
    None
    """  

    E3s = [E3.split(".")[0] for E3 in os.listdir(df_dir)]

    for E3 in E3s:
        
        # Read df containing sequence information in the form:
        #  substrate, sequence, start, end, iteration
        df = pd.read_csv(df_dir+E3+".tsv", sep = "\t")

        # Generate fasta file with substrate and sequence per E3 ligase
        with open(fasta_dir+E3+".fasta", "w") as f:
            for idx, row in df.iterrows():
                f.write(">"+row.substrate+"\n")
                f.write(row.sequence+"\n")


def generate_fasta_from_df(fasta_dir, df, name_col, start_col, end_col, sequence_col):
    """ 
    Reads a dataframe containing at name_row the gene's name or AC, at start_col
    gene's start position, at end_col gene's end position and, at sequence_row, 
    the gene's sequence. Generates a fasta file with gene name or AC + start + end
    and sequence. 

    Note: adding start and end position is to differentiate sequences coming from
    the same substrate, for those programs that need to differenciate between them
    (e.g. clustalo.py)
    
    Parameters
    ----------
    fasta_dir: str
           Path to the folder where the fasta file will be located
    df: pandas dataframe
           Path to the folder where the dataframe with gene name/AC and
           sequence is located
    name_col: str
            Dataframe's column name for gene's name or AC
    start_col: str
            Dataframe's column name for gene's sequence start position
    end_col: str
            Dataframe's column name for gene's sequence end position
    sequence_col
            Dataframe's column name for gene's sequence
    Returns
    -------
    None
    """  

    with open(fasta_dir+".fasta", "w") as f:   # Open fasta file in append mode

        for idx, row in df.iterrows():

            # Gene name+start+end and sequence
            name = row[name_col]+"_"+str(int(row[start_col]))+"_"+str(int((row[end_col])))    
            
            sequence = str(row[sequence_col])  # str transformation because sometimes 
                                               # is interpreted as float

            # Write in fasta format
            f.write(">"+name+"\n")
            f.write(sequence+"\n")   
