def check_overlap(row, matrix_substrate, which = False):
    """
    Takes the row and a matrix of the same substrate and checks if the sequence in the row has 
    overlap with any of the sequences in the matrix.

    Parameters
    ----------
    row: pandas dataframe row
            Row that contains the protein sequence to check the overlap
    matrix_substrate: pandas dataframe
            Dataframe containing those protein sequences to compare to check for overlaps
    which: boolean (default: False)
            If True, returns a list of T/F it the sequence is included or not in the sequences of 
            matrix_substrate
            If False, returns T/F if there is any overlap
    Returns
    -------
            If which = True: boolean indicating whether there is any overlap
            If which = False: list of booleans indicating whethere there is overlap with all the analyzed
            sequences
    """
    if row['Substrate'] != list(set(matrix_substrate['Substrate']))[0] or len(set(matrix_substrate['Substrate'])) > 1:
        raise ValueError('The matrix is not correct')
        
    seq = set(range(row['Start'], row['End']))
    list_seqs = [set(range(start, end)) for start, end in zip(list(matrix_substrate['Start']),
                                                          list(matrix_substrate['End']))]
            
    # check that the positions are not overlapping
    if not which:
        return (any([(seq.issubset(x) or x.issubset(seq)) for x in list_seqs]))    
    else:
        return ([(seq.issubset(x) or x.issubset(seq)) for x in list_seqs])


def completed_sequence(protein, sequence, start, end, proteome, total_len = 12, half = True):
    """
    Takes as input a sequence (protein, start, end) and the length for the output sequence to have

    Parameters
    ----------
    protein: str
            Protein ID
    sequence: str
            Degron sequence
    start: int
            Starting position of the degron in the complete protein sequence
    end: int
            Ending position of the degron in the complete protein sequence
    proteome: dict
            Dictionary of the form {protein_ID: sequence}. Contains a set of protein sequences (e.g.: proteome)
    total_len: int (default: 12)
            Length to which to extend the degron sequence
    half: boolean or int (default: True)
            If True, adds the same number of amino acids to both degron sequence edges. If the number of amino
            acids to add is odd, adds one more amino acid at the beginning of the sequence.
            If False, it is possible to specify the number of amino acids added at the beginning and at the end
            of the degron sequence.
    
    Returns
    -------
    new_seq, new_start, new_end: str, int, int
            Extended degron sequence and start and end position in the protein sequence

    """
    
    # Extract the complete protein sequence
    complete_seq = proteome[protein]

    # Check whether the extraction is correct
    if complete_seq[start-1:end] != sequence:
        raise ValueError('The sequence or positions may be incorrect.')

    else:
        # Extended the degron sequence equally from the two edges
        if half == True: 
            to_add = total_len-(end-start+1)
            # Add even amino acids
            if to_add%2 == 0:
                new_start = int(start-to_add/2)
                new_end = int(end+to_add/2)
            # Add odd amino acids
            else:
                new_start = int(start-int(to_add/2)-1)
                new_end = int(end+int(to_add/2))

        # Customized extension: the number of amino acids introduced is added at the beginning of the sequence; the rest, after
        else: 
            to_add_before = half
            to_add_after = total_len-(end-start+1)-to_add_before
            if to_add_after < 0:                # check that the number of aa to add is correct
                raise ValueError('Check the number of aa to add.')
            else:
                new_start = start-to_add_before
                new_end = end+to_add_after

    # If we arrived at the beginning of the protein: add all amino acids at the end of the degron sequence        
    if new_start <= 0: 
        return(completed_sequence(protein,sequence, start, end, proteome, total_len, half = 0))
    # If we arrived at the end of the protein: add the remaining amino acids at the beginning of the degron sequence
    elif new_end > len(complete_seq):
        to_add_before = total_len-len(sequence)
        return(completed_sequence(protein, sequence, start, end, proteome, total_len, half = to_add_before))
    
    else:
        new_seq = complete_seq[new_start-1:new_end]
    
    return (new_seq, new_start, new_end)    
