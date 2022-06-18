def check_overlap(row, matrix_substrate, which = False):
    """
    Function that takes the row and a matrix of the same substrate and checks if the sequence in the row has 
    overlap with any of the sequences in the matrix.
    By default (which==False) returns T/F (True if there is an overlap).
    If which==True returns a list of T/F if it is included or not
    """
    if row['Substrate'] != list(set(matrix_substrate['Substrate']))[0] or len(set(matrix_substrate['Substrate'])) > 1:
        raise ValueError('The matrix is not correct')
        
    seq = set(range(row['Start'], row['End']))
    list_seqs = [set(range(start, end)) for start, end in zip(list(matrix_substrate['Start']),
                                                          list(matrix_substrate['End']))]
            
    #check that the positions are not overlapping
    if not which:
        return (any([(seq.issubset(x) or x.issubset(seq)) for x in list_seqs]))    
    else:
        return ([(seq.issubset(x) or x.issubset(seq)) for x in list_seqs])






def completed_sequence(protein,sequence,start,end,total_len=12,half=True,dic_fasta_seqs=False):
    """
    start and end from 1!!!
    Takes as input a sequence (protein,start,end) and the total lenght that is the length of the output sequence.
    
    Returns a sequence of the total_len adding aa before and after the original seq and its start-end positions. 
        -If half==True, it adds hald and half on each side (when the number of aa to add is odd, one more is added before
    than after), 
        -if not you can specify how many aa are added before (and the rest are added after).
    """
    if not dic_fasta_seqs:
        with open('/workspace/users/pgomis/degrons/global_classifier/fasta_sequences.json','rt') as d:
            dic_fasta_seqs=json.load(d)

        
    complete_seq=dic_fasta_seqs[protein]
    if complete_seq[start-1:end]!=sequence:
        raise ValueError('The sequence or positions may be incorrect.')
    else:
        if half==True: #when half is selected, the aa added are distributed equally at both sides
            to_add=total_len-(end-start+1)
            if to_add%2==0:
                new_start=int(start-to_add/2)
                new_end=int(end+to_add/2)
            else:
                new_start=int(start-int(to_add/2)-1)
                new_end=int(end+int(to_add/2))
        else: #when half is set as a number, that number of aa are added before the sequence (and the rest after)
            to_add_before=half
            to_add_after=total_len-(end-start+1)-to_add_before
            if to_add_after<0: #check that the number of aa to add is correct
                raise ValueError('Check the number of aa to add.')
            else:
                new_start=start-to_add_before
                new_end=end+to_add_after
            
    if new_start<=0: #when we have arrived the beginning of the protein
        return(completed_sequence(protein,sequence,start,end,total_len,half=0,dic_fasta_seqs=dic_fasta_seqs))
    elif new_end>len(complete_seq):#when we have arrived the end of the protein
        to_add_before=total_len-len(sequence)
        return(completed_sequence(protein,sequence,start,end,total_len,half=to_add_before,dic_fasta_seqs=dic_fasta_seqs))
    
    else:
        new_seq=complete_seq[new_start-1:new_end]
    
    return(new_seq,new_start,new_end)    
