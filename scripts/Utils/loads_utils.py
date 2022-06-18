# -- Import libraries -- # 

import pandas as pd
import logomaker
import matplotlib.pyplot as plt
import random
from scipy import stats

# -- Functions -- #

def load_bg_aa(aa_bg_dir):
    """
    Loads aminoacids background probabilities and one-character names separately from a 
    dataframe containing both information
    
    Parameters
    ----------
    aa_bg_dir: str
                Path to the folder where the aminoacids background probabilities file is located                
    Returns
    -------
    aa_probs: numpy array
                1-D array containing aminoacid background probabilities alphabetically
                sorted by one-character aminoacid names
    aa_names: numpy array
                1-D array of one-character aminoacid names alphabetically sorted
    """

    # Load table sorted alphabetically
    bg_aa_matrix = pd.read_table(aa_bg_dir).sort_values(by = "Aminoacid")

    # Separate aminoacids from bg probabilities
    aa_probs = bg_aa_matrix["Frequency"].to_numpy()
    aa_names = bg_aa_matrix["Aminoacid"].to_numpy() 

    return aa_probs, aa_names
