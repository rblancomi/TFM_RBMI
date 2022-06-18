# -- Import libraries -- # 

import pandas as pd
from scipy import stats

# -- Functions -- #

def process_muttype_df(muttype_df, cols_for_drop, cond, check_normality = True, for_one_degron = False):
    
    if for_one_degron:
        muttype_df = muttype_df.drop_duplicates(subset = cols_for_drop)
        muttype_df["Condition"] = cond
        # may need to add drop duplicates here too
        if check_normality:
            stabch_levels = muttype_df.Stability_Change.values
            try:
                print(f"Normality test {cond}: {stats.normaltest(stabch_levels).pvalue}")
            except ValueError: # ValueError: skewtest is not valid with less than 8 samples
                pass
        
    else:    
        muttype_df = muttype_df.drop_duplicates(subset = cols_for_drop)
        
        if len(muttype_df) >= 9:   # samples <=8 generate recurrence error 
            muttype_df["Condition"] = cond
            
            if check_normality:
                stabch_levels = muttype_df.Stability_Change.values
                try:
                    print(f"Normality test {cond}: {stats.normaltest(stabch_levels).pvalue}")
                except ValueError: # ValueError: skewtest is not valid with less than 8 samples
                    pass
                
    return muttype_df

def prepare_subsets_muttype_dict(stabch, conditions, check_normality = True, Altered_E3_Ligases = False):
    """
    """
    
    # Dictionary to store mut type subsets
    subsets_muttype_dict = {}
    
    # To drop duplicates
    cols_to_remove = ["E3", "degron_start", "degron_end", "Loc_mut_degron"]
    cols = stabch.columns.tolist()
    cols_for_drop = [col for col in cols if col not in cols_to_remove]
    
    # WT
    if "wt" in conditions:
        wt = stabch.loc[(stabch.Phenotype == "WT") & (stabch.Altered_E3_Ligases == Altered_E3_Ligases)]
        subsets_muttype_dict["wt"] = process_muttype_df(wt, cols_for_drop, "wt", check_normality = check_normality)

    # Synonymous variants
    if "syn_muts" in conditions:
        syn_muts = stabch.loc[(stabch.Phenotype == "synonymous_variant") & (stabch.Altered_E3_Ligases == Altered_E3_Ligases)]
        subsets_muttype_dict["syn_muts"] = process_muttype_df(syn_muts, cols_for_drop, "syn_muts", check_normality = check_normality)
    
    if "syn_in_muts" in conditions:
        syn_in_muts = stabch.loc[(stabch.Phenotype == "synonymous_variant") & (stabch.Loc_mut_degron == "inside") &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases)]
        subsets_muttype_dict["syn_in_muts"] = process_muttype_df(syn_in_muts, cols_for_drop, "syn_in_muts", check_normality = check_normality)

    if "syn_out_muts" in conditions:
        syn_out_muts = stabch.loc[(stabch.Phenotype == "synonymous_variant") & (stabch.Loc_mut_degron.str.contains("outside")) &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases)]
        subsets_muttype_dict["syn_out_muts"] = process_muttype_df(syn_out_muts, cols_for_drop, "syn_out_muts", check_normality = check_normality)


    # Missense variants
    if "miss_muts" in conditions:
        miss_muts = stabch.loc[(stabch.Phenotype == "missense_variant") & (stabch.Altered_E3_Ligases == Altered_E3_Ligases)]
        subsets_muttype_dict["miss_muts"] = process_muttype_df(miss_muts, cols_for_drop, "miss_muts", check_normality = check_normality)


    if "miss_in_muts" in conditions:
        miss_in_muts = stabch.loc[(stabch.Phenotype == "missense_variant") & (stabch.Loc_mut_degron == "inside") &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases)]
        subsets_muttype_dict["miss_in_muts"] = process_muttype_df(miss_in_muts, cols_for_drop, "miss_in_muts", check_normality = check_normality)


    if "miss_out_muts" in conditions:
        miss_out_muts = stabch.loc[(stabch.Phenotype == "missense_variant") & (stabch.Loc_mut_degron.str.contains("outside")) &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases)]
        subsets_muttype_dict["miss_out_muts"] = process_muttype_df(miss_out_muts, cols_for_drop, "miss_out_muts", check_normality = check_normality)


    # Inframe indels 
    if "inframe_muts" in conditions:
        inframe_muts = stabch.loc[(stabch.Phenotype.str.contains("inframe")) & (stabch.Altered_E3_Ligases == Altered_E3_Ligases)]
        subsets_muttype_dict["inframe_muts"] = process_muttype_df(inframe_muts, cols_for_drop, "inframe_muts", check_normality = check_normality)

    if "inframe_in_muts" in conditions:
        inframe_in_muts = stabch.loc[(stabch.Phenotype.str.contains("inframe")) & (stabch.Loc_mut_degron == "inside") &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases)]
        subsets_muttype_dict["inframe_in_muts"] = process_muttype_df(inframe_in_muts, cols_for_drop, "inframe_in_muts", check_normality = check_normality)

    if "inframe_out_muts" in conditions:
        inframe_out_muts = stabch.loc[(stabch.Phenotype.str.contains("inframe")) & (stabch.Loc_mut_degron.str.contains("outside")) &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases)]
        subsets_muttype_dict["inframe_out_muts"] = process_muttype_df(inframe_out_muts, cols_for_drop, "inframe_out_muts", check_normality = check_normality)


    # Nonsense variants in last gene's exon (stop gained)
    if "nonsense_muts" in conditions:
        nonsense_muts = stabch.loc[(stabch.Phenotype == "stop_gained") & 
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.Mut_in_lastexon == True)]
        subsets_muttype_dict["nonsense_muts"] = process_muttype_df(nonsense_muts, cols_for_drop, "nonsense_muts", check_normality = check_normality)

    
    if "nonsense_inbf_muts" in conditions:
        nonsense_inbf_muts = stabch.loc[(stabch.Phenotype == "stop_gained") & ((stabch.Loc_mut_degron == "inside") |
        (stabch.Loc_mut_degron == "outside_before")) &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.Mut_in_lastexon == True)]
        subsets_muttype_dict["nonsense_inbf_muts"] = process_muttype_df(nonsense_inbf_muts, cols_for_drop, "nonsense_inbf_muts", check_normality = check_normality)

    
    if "nonsense_aft_muts" in conditions:
        nonsense_aft_muts = stabch.loc[(stabch.Phenotype == "stop_gained") & (stabch.Loc_mut_degron == "outside_after") &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.Mut_in_lastexon == True)]
        subsets_muttype_dict["nonsense_aft_muts"] = process_muttype_df(nonsense_aft_muts, cols_for_drop, "nonsense_aft_muts", check_normality = check_normality)


    # Frameshif variants in last gene's exon 
    if "frameshift_muts" in conditions:
        frameshift_muts = stabch.loc[(stabch.Phenotype == "frameshift_variant") & 
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.Mut_in_lastexon == True)]
        subsets_muttype_dict["frameshift_muts"] = process_muttype_df(frameshift_muts, cols_for_drop, "frameshift_muts", check_normality = check_normality)

        
    if "frameshift_inbf_muts" in conditions:
        frameshift_inbf_muts = stabch.loc[(stabch.Phenotype == "frameshift_variant") & ((stabch.Loc_mut_degron == "inside") |
        (stabch.Loc_mut_degron == "outside_before")) &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.Mut_in_lastexon == True)]
        subsets_muttype_dict["frameshift_inbf_muts"] = process_muttype_df(frameshift_inbf_muts, cols_for_drop, "frameshift_inbf_muts", check_normality = check_normality)


    if "frameshift_aft_muts" in conditions:
        frameshift_aft_muts = stabch.loc[(stabch.Phenotype == "frameshift_variant") & (stabch.Loc_mut_degron == "outside_after") &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.Mut_in_lastexon == True)]
        subsets_muttype_dict["frameshift_aft_muts"] = process_muttype_df(frameshift_aft_muts, cols_for_drop, "frameshift_aft_muts", check_normality = check_normality)


    # Non-truncanting variants: missense and inframe
    if "nontrunc_muts" in conditions:
        nontrunc_muts = stabch.loc[((stabch.Phenotype == "missense_variant") |
        (stabch.Phenotype.str.contains("inframe"))) & (stabch.Altered_E3_Ligases == Altered_E3_Ligases)]
        subsets_muttype_dict["nontrunc_muts"] = process_muttype_df(nontrunc_muts, cols_for_drop, "nontrunc_muts", check_normality = check_normality)

        
    if "nontrunc_in_muts" in conditions:
        nontrunc_in_muts = stabch.loc[((stabch.Phenotype == "missense_variant") |
        (stabch.Phenotype.str.contains("inframe"))) & (stabch.Loc_mut_degron == "inside") &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases)]
        subsets_muttype_dict["nontrunc_in_muts"] = process_muttype_df(nontrunc_in_muts, cols_for_drop, "nontrunc_in_muts", check_normality = check_normality)

    
    if "nontrunc_out_muts" in conditions:
        nontrunc_out_muts = stabch.loc[((stabch.Phenotype == "missense_variant") |
        (stabch.Phenotype.str.contains("inframe"))) & (stabch.Loc_mut_degron.str.contains("outside")) &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases)]
        subsets_muttype_dict["nontrunc_out_muts"] = process_muttype_df(nontrunc_out_muts, cols_for_drop, "nontrunc_out_muts", check_normality = check_normality)


    # Truncating variants (last exon): nonsense and frameshift   
    if "trunc_muts" in conditions:
        trunc_muts = stabch.loc[((stabch.Phenotype == "stop_gained") |
        (stabch.Phenotype == "frameshift_variant")) & (stabch.Altered_E3_Ligases == Altered_E3_Ligases)
        & (stabch.Mut_in_lastexon == True)]
        subsets_muttype_dict["trunc_muts"] = process_muttype_df(trunc_muts, cols_for_drop, "trunc_muts", check_normality = check_normality)

        
    if "trunc_inbf_muts" in conditions:
        trunc_inbf_muts = stabch.loc[((stabch.Phenotype == "stop_gained") |
        (stabch.Phenotype == "frameshift_variant")) & ((stabch.Loc_mut_degron == "inside") |
        (stabch.Loc_mut_degron == "outside_before")) & (stabch.Altered_E3_Ligases == Altered_E3_Ligases)
        & (stabch.Mut_in_lastexon == True)]
        subsets_muttype_dict["trunc_inbf_muts"] = process_muttype_df(trunc_inbf_muts, cols_for_drop, "trunc_inbf_muts", check_normality = check_normality)
    
    if "trunc_aft_muts" in conditions:
        trunc_aft_muts = stabch.loc[((stabch.Phenotype == "stop_gained") |
        (stabch.Phenotype == "frameshift_variant")) & (stabch.Loc_mut_degron == "outside_after") &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases)
        & (stabch.Mut_in_lastexon == True)]
        subsets_muttype_dict["trunc_aft_muts"] = process_muttype_df(trunc_aft_muts, cols_for_drop, "trunc_aft_muts", check_normality = check_normality)
    
    # remove conditions with size =< 9
    subsets_muttype_dict_f = {k:v for (k,v) in subsets_muttype_dict.items() if len(v) >= 9}

    return subsets_muttype_dict_f

def prepare_subsets_muttype_wthE3_dict(stabch, conditions, check_normality = True, Altered_E3_Ligases = False):
    """
    """
    
    # Dictionary to store mut type subsets
    subsets_muttype_dict = {}
    
    # To drop duplicates
    cols_to_remove = ["degron_start", "degron_end"]
    cols = stabch.columns.tolist()
    cols_for_drop = [col for col in cols if col not in cols_to_remove]
    
    # WT
    if "wt" in conditions:
        wt = stabch.loc[(stabch.Phenotype == "WT") & (stabch.Altered_E3_Ligases == Altered_E3_Ligases)]
        subsets_muttype_dict["wt"] = process_muttype_df(wt, cols_for_drop, "wt", check_normality = check_normality)

    # Synonymous variants
    if "syn_muts" in conditions:
        syn_muts = stabch.loc[(stabch.Phenotype == "synonymous_variant") & (stabch.Altered_E3_Ligases == Altered_E3_Ligases)]
        subsets_muttype_dict["syn_muts"] = process_muttype_df(syn_muts, cols_for_drop, "syn_muts", check_normality = check_normality)
    
    if "syn_in_muts" in conditions:
        syn_in_muts = stabch.loc[(stabch.Phenotype == "synonymous_variant") & (stabch.Loc_mut_degron == "inside") &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases)]
        subsets_muttype_dict["syn_in_muts"] = process_muttype_df(syn_in_muts, cols_for_drop, "syn_in_muts", check_normality = check_normality)

    if "syn_out_muts" in conditions:
        syn_out_muts = stabch.loc[(stabch.Phenotype == "synonymous_variant") & (stabch.Loc_mut_degron.str.contains("outside")) &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases)]
        subsets_muttype_dict["syn_out_muts"] = process_muttype_df(syn_out_muts, cols_for_drop, "syn_out_muts", check_normality = check_normality)


    # Missense variants
    if "miss_muts" in conditions:
        miss_muts = stabch.loc[(stabch.Phenotype == "missense_variant") & (stabch.Altered_E3_Ligases == Altered_E3_Ligases)]
        subsets_muttype_dict["miss_muts"] = process_muttype_df(miss_muts, cols_for_drop, "miss_muts", check_normality = check_normality)


    if "miss_in_muts" in conditions:
        miss_in_muts = stabch.loc[(stabch.Phenotype == "missense_variant") & (stabch.Loc_mut_degron == "inside") &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases)]
        subsets_muttype_dict["miss_in_muts"] = process_muttype_df(miss_in_muts, cols_for_drop, "miss_in_muts", check_normality = check_normality)


    if "miss_out_muts" in conditions:
        miss_out_muts = stabch.loc[(stabch.Phenotype == "missense_variant") & (stabch.Loc_mut_degron.str.contains("outside")) &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases)]
        subsets_muttype_dict["miss_out_muts"] = process_muttype_df(miss_out_muts, cols_for_drop, "miss_out_muts", check_normality = check_normality)


    # Inframe indels 
    if "inframe_muts" in conditions:
        inframe_muts = stabch.loc[(stabch.Phenotype.str.contains("inframe")) & (stabch.Altered_E3_Ligases == Altered_E3_Ligases)]
        subsets_muttype_dict["inframe_muts"] = process_muttype_df(inframe_muts, cols_for_drop, "inframe_muts", check_normality = check_normality)

    if "inframe_in_muts" in conditions:
        inframe_in_muts = stabch.loc[(stabch.Phenotype.str.contains("inframe")) & (stabch.Loc_mut_degron == "inside") &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases)]
        subsets_muttype_dict["inframe_in_muts"] = process_muttype_df(inframe_in_muts, cols_for_drop, "inframe_in_muts", check_normality = check_normality)

    if "inframe_out_muts" in conditions:
        inframe_out_muts = stabch.loc[(stabch.Phenotype.str.contains("inframe")) & (stabch.Loc_mut_degron.str.contains("outside")) &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases)]
        subsets_muttype_dict["inframe_out_muts"] = process_muttype_df(inframe_out_muts, cols_for_drop, "inframe_out_muts", check_normality = check_normality)


    # Nonsense variants in last gene's exon (stop gained)
    if "nonsense_muts" in conditions:
        nonsense_muts = stabch.loc[(stabch.Phenotype == "stop_gained") & 
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.Mut_in_lastexon == True)]
        subsets_muttype_dict["nonsense_muts"] = process_muttype_df(nonsense_muts, cols_for_drop, "nonsense_muts", check_normality = check_normality)

    
    if "nonsense_inbf_muts" in conditions:
        nonsense_inbf_muts = stabch.loc[(stabch.Phenotype == "stop_gained") & ((stabch.Loc_mut_degron == "inside") |
        (stabch.Loc_mut_degron == "outside_before")) &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.Mut_in_lastexon == True)]
        subsets_muttype_dict["nonsense_inbf_muts"] = process_muttype_df(nonsense_inbf_muts, cols_for_drop, "nonsense_inbf_muts", check_normality = check_normality)

    
    if "nonsense_aft_muts" in conditions:
        nonsense_aft_muts = stabch.loc[(stabch.Phenotype == "stop_gained") & (stabch.Loc_mut_degron == "outside_after") &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.Mut_in_lastexon == True)]
        subsets_muttype_dict["nonsense_aft_muts"] = process_muttype_df(nonsense_aft_muts, cols_for_drop, "nonsense_aft_muts", check_normality = check_normality)


    # Frameshif variants in last gene's exon 
    if "frameshift_muts" in conditions:
        frameshift_muts = stabch.loc[(stabch.Phenotype == "frameshift_variant") & 
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.Mut_in_lastexon == True)]
        subsets_muttype_dict["frameshift_muts"] = process_muttype_df(frameshift_muts, cols_for_drop, "frameshift_muts", check_normality = check_normality)

        
    if "frameshift_inbf_muts" in conditions:
        frameshift_inbf_muts = stabch.loc[(stabch.Phenotype == "frameshift_variant") & ((stabch.Loc_mut_degron == "inside") |
        (stabch.Loc_mut_degron == "outside_before")) &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.Mut_in_lastexon == True)]
        subsets_muttype_dict["frameshift_inbf_muts"] = process_muttype_df(frameshift_inbf_muts, cols_for_drop, "frameshift_inbf_muts", check_normality = check_normality)


    if "frameshift_aft_muts" in conditions:
        frameshift_aft_muts = stabch.loc[(stabch.Phenotype == "frameshift_variant") & (stabch.Loc_mut_degron == "outside_after") &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.Mut_in_lastexon == True)]
        subsets_muttype_dict["frameshift_aft_muts"] = process_muttype_df(frameshift_aft_muts, cols_for_drop, "frameshift_aft_muts", check_normality = check_normality)


    # Non-truncanting variants: missense and inframe
    if "nontrunc_muts" in conditions:
        nontrunc_muts = stabch.loc[((stabch.Phenotype == "missense_variant") |
        (stabch.Phenotype.str.contains("inframe"))) & (stabch.Altered_E3_Ligases == Altered_E3_Ligases)]
        subsets_muttype_dict["nontrunc_muts"] = process_muttype_df(nontrunc_muts, cols_for_drop, "nontrunc_muts", check_normality = check_normality)

        
    if "nontrunc_in_muts" in conditions:
        nontrunc_in_muts = stabch.loc[((stabch.Phenotype == "missense_variant") |
        (stabch.Phenotype.str.contains("inframe"))) & (stabch.Loc_mut_degron == "inside") &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases)]
        subsets_muttype_dict["nontrunc_in_muts"] = process_muttype_df(nontrunc_in_muts, cols_for_drop, "nontrunc_in_muts", check_normality = check_normality)

    
    if "nontrunc_out_muts" in conditions:
        nontrunc_out_muts = stabch.loc[((stabch.Phenotype == "missense_variant") |
        (stabch.Phenotype.str.contains("inframe"))) & (stabch.Loc_mut_degron.str.contains("outside")) &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases)]
        subsets_muttype_dict["nontrunc_out_muts"] = process_muttype_df(nontrunc_out_muts, cols_for_drop, "nontrunc_out_muts", check_normality = check_normality)


    # Truncating variants (last exon): nonsense and frameshift   
    if "trunc_muts" in conditions:
        trunc_muts = stabch.loc[((stabch.Phenotype == "stop_gained") |
        (stabch.Phenotype == "frameshift_variant")) & (stabch.Altered_E3_Ligases == Altered_E3_Ligases)
        & (stabch.Mut_in_lastexon == True)]
        subsets_muttype_dict["trunc_muts"] = process_muttype_df(trunc_muts, cols_for_drop, "trunc_muts", check_normality = check_normality)

        
    if "trunc_inbf_muts" in conditions:
        trunc_inbf_muts = stabch.loc[((stabch.Phenotype == "stop_gained") |
        (stabch.Phenotype == "frameshift_variant")) & ((stabch.Loc_mut_degron == "inside") |
        (stabch.Loc_mut_degron == "outside_before")) & (stabch.Altered_E3_Ligases == Altered_E3_Ligases)
        & (stabch.Mut_in_lastexon == True)]
        subsets_muttype_dict["trunc_inbf_muts"] = process_muttype_df(trunc_inbf_muts, cols_for_drop, "trunc_inbf_muts", check_normality = check_normality)
    
    if "trunc_aft_muts" in conditions:
        trunc_aft_muts = stabch.loc[((stabch.Phenotype == "stop_gained") |
        (stabch.Phenotype == "frameshift_variant")) & (stabch.Loc_mut_degron == "outside_after") &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases)
        & (stabch.Mut_in_lastexon == True)]
        subsets_muttype_dict["trunc_aft_muts"] = process_muttype_df(trunc_aft_muts, cols_for_drop, "trunc_aft_muts", check_normality = check_normality)
    
    # remove conditions with size =< 9
    subsets_muttype_dict_f = {k:v for (k,v) in subsets_muttype_dict.items() if len(v) >= 9}

    return subsets_muttype_dict_f

def prepare_subsetsE3ligase_muttype_dict(stabch, conditions, E3, check_normality = True, Altered_E3_Ligases = False):
    """
    """
    
    # Dictionary to store mut type subsets of the stabch_levels table
    subsets_muttype_dict = {}
    
    # To drop duplicates
    cols_to_remove = ["E3", "degron_start", "degron_end", "Loc_mut_degron"]
    cols = stabch.columns.tolist()
    cols_for_drop = [col for col in cols if col not in cols_to_remove]
        
    # WT (only targets of the E3 ligase) 
    if "wt" in conditions:
        wt = stabch.loc[(stabch.Phenotype == "WT") & (stabch.E3 == E3) & 
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases)]
        subsets_muttype_dict["wt"] = process_muttype_df(wt, cols_for_drop, "wt", check_normality = check_normality)

    # Synonymous variants
    if "syn_muts" in conditions:
        syn_muts = stabch.loc[(stabch.Phenotype == "synonymous_variant") & (stabch.Altered_E3_Ligases == Altered_E3_Ligases)
        & (stabch.E3 == E3)]
        subsets_muttype_dict["syn_muts"] = process_muttype_df(syn_muts, cols_for_drop, "syn_muts", check_normality = check_normality)
        
    
    if "syn_in_muts" in conditions:
        syn_in_muts = stabch.loc[(stabch.Phenotype == "synonymous_variant") & (stabch.Loc_mut_degron == "inside") &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.E3 == E3)]
        subsets_muttype_dict["syn_in_muts"] = process_muttype_df(syn_in_muts, cols_for_drop, "syn_in_muts", check_normality = check_normality)
        

    if "syn_out_muts" in conditions:
        syn_out_muts = stabch.loc[(stabch.Phenotype == "synonymous_variant") & (stabch.Loc_mut_degron.str.contains("outside")) &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.E3 == E3)]
        subsets_muttype_dict["syn_out_muts"] = process_muttype_df(syn_out_muts, cols_for_drop, "syn_out_muts", check_normality = check_normality)
        

    # Missense variants
    if "miss_muts" in conditions:
        miss_muts = stabch.loc[(stabch.Phenotype == "missense_variant") & (stabch.Altered_E3_Ligases == Altered_E3_Ligases)
        & (stabch.E3 == E3)]
        subsets_muttype_dict["miss_muts"] = process_muttype_df(miss_muts, cols_for_drop, "miss_muts", check_normality = check_normality)
        
            
    if "miss_in_muts" in conditions:
        miss_in_muts = stabch.loc[(stabch.Phenotype == "missense_variant") & (stabch.Loc_mut_degron == "inside") &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.E3 == E3)]
        subsets_muttype_dict["miss_in_muts"] = process_muttype_df(miss_in_muts, cols_for_drop, "miss_in_muts", check_normality = check_normality)
        
            
    if "miss_out_muts" in conditions:
        miss_out_muts = stabch.loc[(stabch.Phenotype == "missense_variant") & (stabch.Loc_mut_degron.str.contains("outside")) &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.E3 == E3)]
        subsets_muttype_dict["miss_out_muts"] = process_muttype_df(miss_out_muts, cols_for_drop, "miss_out_muts", check_normality = check_normality)
        
            
    # Inframe indels 
    if "inframe_muts" in conditions:
        inframe_muts = stabch.loc[(stabch.Phenotype.str.contains("inframe")) & (stabch.Altered_E3_Ligases == Altered_E3_Ligases)
        & (stabch.E3 == E3)]
        subsets_muttype_dict["inframe_muts"] = process_muttype_df(inframe_muts, cols_for_drop, "inframe_muts", check_normality = check_normality)
        
        
    if "inframe_in_muts" in conditions:
        inframe_in_muts = stabch.loc[(stabch.Phenotype.str.contains("inframe")) & (stabch.Loc_mut_degron == "inside") &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.E3 == E3)]
        subsets_muttype_dict["inframe_in_muts"] = process_muttype_df(inframe_in_muts, cols_for_drop, "inframe_in_muts", check_normality = check_normality)
        
            
    if "inframe_out_muts" in conditions:
        inframe_out_muts = stabch.loc[(stabch.Phenotype.str.contains("inframe")) & (stabch.Loc_mut_degron.str.contains("outside")) &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.E3 == E3)]
        subsets_muttype_dict["inframe_out_muts"] = process_muttype_df(inframe_out_muts, cols_for_drop, "inframe_out_muts", check_normality = check_normality)
        
            
    # Nonsense variants in last gene's exon (stop gained)
    if "nonsense_muts" in conditions:
        nonsense_muts = stabch.loc[(stabch.Phenotype == "stop_gained") & 
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.Mut_in_lastexon == True) & (stabch.E3 == E3)]
        subsets_muttype_dict["nonsense_muts"] = process_muttype_df(nonsense_muts, cols_for_drop, "nonsense_muts", check_normality = check_normality)
        
            
    if "nonsense_inbf_muts" in conditions:
        nonsense_inbf_muts = stabch.loc[(stabch.Phenotype == "stop_gained") & ((stabch.Loc_mut_degron == "inside") |
        (stabch.Loc_mut_degron == "outside_before")) &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.Mut_in_lastexon == True) & (stabch.E3 == E3)]
        subsets_muttype_dict["nonsense_inbf_muts"] = process_muttype_df(nonsense_inbf_muts, cols_for_drop, "nonsense_inbf_muts", check_normality = check_normality)
        
    
    if "nonsense_aft_muts" in conditions:
        nonsense_aft_muts = stabch.loc[(stabch.Phenotype == "stop_gained") & (stabch.Loc_mut_degron == "outside_after") &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.Mut_in_lastexon == True) & (stabch.E3 == E3)]
        subsets_muttype_dict["nonsense_aft_muts"] = process_muttype_df(nonsense_aft_muts, cols_for_drop, "nonsense_aft_muts", check_normality = check_normality)
        

    # Frameshif variants in last gene's exon 
    if "frameshift_muts" in conditions:
        frameshift_muts = stabch.loc[(stabch.Phenotype == "frameshift_variant") & 
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.Mut_in_lastexon == True) & (stabch.E3 == E3)]
        subsets_muttype_dict["frameshift_muts"] = process_muttype_df(frameshift_muts, cols_for_drop, "frameshift_muts", check_normality = check_normality)
        
        
    if "frameshift_inbf_muts" in conditions:
        frameshift_inbf_muts = stabch.loc[(stabch.Phenotype == "frameshift_variant") & ((stabch.Loc_mut_degron == "inside") |
        (stabch.Loc_mut_degron == "outside_before")) &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.Mut_in_lastexon == True) & (stabch.E3 == E3)]
        subsets_muttype_dict["frameshift_inbf_muts"] = process_muttype_df(frameshift_inbf_muts, cols_for_drop, "frameshift_inbf_muts", check_normality = check_normality)
        
            
    if "frameshift_aft_muts" in conditions:
        frameshift_aft_muts = stabch.loc[(stabch.Phenotype == "frameshift_variant") & (stabch.Loc_mut_degron == "outside_after") &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.Mut_in_lastexon == True) & (stabch.E3 == E3)]
        subsets_muttype_dict["frameshift_aft_muts"] = process_muttype_df(frameshift_aft_muts, cols_for_drop, "frameshift_aft_muts", check_normality = check_normality)
        
            
    # Non-truncanting variants: missense and inframe
    if "nontrunc_muts" in conditions:
        nontrunc_muts = stabch.loc[((stabch.Phenotype == "missense_variant") |
        (stabch.Phenotype.str.contains("inframe"))) & (stabch.Altered_E3_Ligases == Altered_E3_Ligases) 
        & (stabch.E3 == E3)]
        subsets_muttype_dict["nontrunc_muts"] = process_muttype_df(nontrunc_muts, cols_for_drop, "nontrunc_muts", check_normality = check_normality)
        
            
    if "nontrunc_in_muts" in conditions:
        nontrunc_in_muts = stabch.loc[((stabch.Phenotype == "missense_variant") |
        (stabch.Phenotype.str.contains("inframe"))) & (stabch.Loc_mut_degron == "inside") &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.E3 == E3)]
        subsets_muttype_dict["nontrunc_in_muts"] = process_muttype_df(nontrunc_in_muts, cols_for_drop, "nontrunc_in_muts", check_normality = check_normality)
        
            
    if "nontrunc_out_muts" in conditions:
        nontrunc_out_muts = stabch.loc[((stabch.Phenotype == "missense_variant") |
        (stabch.Phenotype.str.contains("inframe"))) & (stabch.Loc_mut_degron.str.contains("outside")) &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.E3 == E3)]
        subsets_muttype_dict["nontrunc_out_muts"] = process_muttype_df(nontrunc_out_muts, cols_for_drop, "nontrunc_out_muts", check_normality = check_normality)
        
            
    # Truncating variants (last exon): nonsense and frameshift   
    if "trunc_muts" in conditions:
        trunc_muts = stabch.loc[((stabch.Phenotype == "stop_gained") |
        (stabch.Phenotype == "frameshift_variant")) & (stabch.Altered_E3_Ligases == Altered_E3_Ligases)
        & (stabch.E3 == E3) & (stabch.Mut_in_lastexon == True)]
        subsets_muttype_dict["trunc_muts"] = process_muttype_df(trunc_muts, cols_for_drop, "trunc_muts", check_normality = check_normality)
        
            
    if "trunc_inbf_muts" in conditions:
        trunc_inbf_muts = stabch.loc[((stabch.Phenotype == "stop_gained") |
        (stabch.Phenotype == "frameshift_variant")) & ((stabch.Loc_mut_degron == "inside") |
        (stabch.Loc_mut_degron == "outside_before")) & (stabch.Altered_E3_Ligases == Altered_E3_Ligases)
        & (stabch.E3 == E3) & (stabch.Mut_in_lastexon == True)]
        subsets_muttype_dict["trunc_inbf_muts"] = process_muttype_df(trunc_inbf_muts, cols_for_drop, "trunc_inbf_muts", check_normality = check_normality)
        
            
    if "trunc_aft_muts" in conditions:
        trunc_aft_muts = stabch.loc[((stabch.Phenotype == "stop_gained") |
        (stabch.Phenotype == "frameshift_variant")) & (stabch.Loc_mut_degron == "outside_after") &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.E3 == E3) & (stabch.Mut_in_lastexon == True)]
        subsets_muttype_dict["trunc_aft_muts"] = process_muttype_df(trunc_aft_muts, cols_for_drop, "trunc_aft_muts", check_normality = check_normality)
    
    # remove conditions with size =< 9
    subsets_muttype_dict_f = {k:v for (k,v) in subsets_muttype_dict.items() if len(v) >= 9}

    return subsets_muttype_dict_f

def prepare_subsetsDegron_muttype_dict(stabch, conditions, E3, gene, start, end, check_normality = True, Altered_E3_Ligases = False):
    """
    """
    
    # Dictionary to store mut type subsets of the stabch_levels table
    subsets_muttype_dict = {}
    
    # To drop duplicates
    cols_to_remove = ["E3", "degron_start", "degron_end", "Loc_mut_degron"]
    cols = stabch.columns.tolist()
    cols_for_drop = [col for col in cols if col not in cols_to_remove]
        
    # WT
    if "wt" in conditions:
        wt = stabch.loc[(stabch.Phenotype == "WT") & (stabch.Altered_E3_Ligases == Altered_E3_Ligases)
        & (stabch.gene == gene)]
        subsets_muttype_dict["wt"] = process_muttype_df(wt, cols_for_drop, "wt", check_normality = check_normality, for_one_degron = True)

    # Synonymous variants
    if "syn_muts" in conditions:
        syn_muts = stabch.loc[(stabch.Phenotype == "synonymous_variant") & (stabch.Altered_E3_Ligases == Altered_E3_Ligases)
        & (stabch.E3 == E3) & (stabch.gene == gene) & (stabch.degron_start == start) & (stabch.degron_end == end)]
        subsets_muttype_dict["syn_muts"] = process_muttype_df(syn_muts, cols_for_drop, "syn_muts", check_normality = check_normality, for_one_degron = True)
        
    
    if "syn_in_muts" in conditions:
        syn_in_muts = stabch.loc[(stabch.Phenotype == "synonymous_variant") & (stabch.Loc_mut_degron == "inside") &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.E3 == E3) & (stabch.gene == gene) & (stabch.degron_start == start) & (stabch.degron_end == end)]
        subsets_muttype_dict["syn_in_muts"] = process_muttype_df(syn_in_muts, cols_for_drop, "syn_in_muts", check_normality = check_normality, for_one_degron = True)
        

    if "syn_out_muts" in conditions:
        syn_out_muts = stabch.loc[(stabch.Phenotype == "synonymous_variant") & (stabch.Loc_mut_degron.str.contains("outside")) &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.E3 == E3) & (stabch.gene == gene) & (stabch.degron_start == start) & (stabch.degron_end == end)]
        subsets_muttype_dict["syn_out_muts"] = process_muttype_df(syn_out_muts, cols_for_drop, "syn_out_muts", check_normality = check_normality, for_one_degron = True)
        

    # Missense variants
    if "miss_muts" in conditions:
        miss_muts = stabch.loc[(stabch.Phenotype == "missense_variant") & (stabch.Altered_E3_Ligases == Altered_E3_Ligases)
        & (stabch.E3 == E3) & (stabch.gene == gene) & (stabch.degron_start == start) & (stabch.degron_end == end)]
        subsets_muttype_dict["miss_muts"] = process_muttype_df(miss_muts, cols_for_drop, "miss_muts", check_normality = check_normality, for_one_degron = True)
        
            
    if "miss_in_muts" in conditions:
        miss_in_muts = stabch.loc[(stabch.Phenotype == "missense_variant") & (stabch.Loc_mut_degron == "inside") &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.E3 == E3) & (stabch.gene == gene) & (stabch.degron_start == start) & (stabch.degron_end == end)]
        subsets_muttype_dict["miss_in_muts"] = process_muttype_df(miss_in_muts, cols_for_drop, "miss_in_muts", check_normality = check_normality, for_one_degron = True)
        
            
    if "miss_out_muts" in conditions:
        miss_out_muts = stabch.loc[(stabch.Phenotype == "missense_variant") & (stabch.Loc_mut_degron.str.contains("outside")) &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.E3 == E3) & (stabch.gene == gene) & (stabch.degron_start == start) & (stabch.degron_end == end)]
        subsets_muttype_dict["miss_out_muts"] = process_muttype_df(miss_out_muts, cols_for_drop, "miss_out_muts", check_normality = check_normality, for_one_degron = True)
        
            
    # Inframe indels 
    if "inframe_muts" in conditions:
        inframe_muts = stabch.loc[(stabch.Phenotype.str.contains("inframe")) & (stabch.Altered_E3_Ligases == Altered_E3_Ligases)
        & (stabch.E3 == E3) & (stabch.gene == gene) & (stabch.degron_start == start) & (stabch.degron_end == end)]
        subsets_muttype_dict["inframe_muts"] = process_muttype_df(inframe_muts, cols_for_drop, "inframe_muts", check_normality = check_normality, for_one_degron = True)
        
        
    if "inframe_in_muts" in conditions:
        inframe_in_muts = stabch.loc[(stabch.Phenotype.str.contains("inframe")) & (stabch.Loc_mut_degron == "inside") &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.E3 == E3) & (stabch.gene == gene) & (stabch.degron_start == start) & (stabch.degron_end == end)]
        subsets_muttype_dict["inframe_in_muts"] = process_muttype_df(inframe_in_muts, cols_for_drop, "inframe_in_muts", check_normality = check_normality, for_one_degron = True)
        
            
    if "inframe_out_muts" in conditions:
        inframe_out_muts = stabch.loc[(stabch.Phenotype.str.contains("inframe")) & (stabch.Loc_mut_degron.str.contains("outside")) &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.E3 == E3) & (stabch.gene == gene) & (stabch.degron_start == start) & (stabch.degron_end == end)]
        subsets_muttype_dict["inframe_out_muts"] = process_muttype_df(inframe_out_muts, cols_for_drop, "inframe_out_muts", check_normality = check_normality, for_one_degron = True)
        
            
    # Nonsense variants in last gene's exon (stop gained)
    if "nonsense_muts" in conditions:
        nonsense_muts = stabch.loc[(stabch.Phenotype == "stop_gained") & 
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.Mut_in_lastexon == True) & (stabch.E3 == E3)
        & (stabch.gene == gene) & (stabch.degron_start == start) & (stabch.degron_end == end)]
        subsets_muttype_dict["nonsense_muts"] = process_muttype_df(nonsense_muts, cols_for_drop, "nonsense_muts", check_normality = check_normality, for_one_degron = True)
        
            
    if "nonsense_inbf_muts" in conditions:
        nonsense_inbf_muts = stabch.loc[(stabch.Phenotype == "stop_gained") & ((stabch.Loc_mut_degron == "inside") |
        (stabch.Loc_mut_degron == "outside_before")) &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.Mut_in_lastexon == True) & (stabch.E3 == E3)
        & (stabch.gene == gene) & (stabch.degron_start == start) & (stabch.degron_end == end)]
        subsets_muttype_dict["nonsense_inbf_muts"] = process_muttype_df(nonsense_inbf_muts, cols_for_drop, "nonsense_inbf_muts", check_normality = check_normality, for_one_degron = True)
        
    
    if "nonsense_aft_muts" in conditions:
        nonsense_aft_muts = stabch.loc[(stabch.Phenotype == "stop_gained") & (stabch.Loc_mut_degron == "outside_after") &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.Mut_in_lastexon == True) & (stabch.E3 == E3)
        & (stabch.gene == gene) & (stabch.degron_start == start) & (stabch.degron_end == end)]
        subsets_muttype_dict["nonsense_aft_muts"] = process_muttype_df(nonsense_aft_muts, cols_for_drop, "nonsense_aft_muts", check_normality = check_normality, for_one_degron = True)
        

    # Frameshif variants in last gene's exon 
    if "frameshift_muts" in conditions:
        frameshift_muts = stabch.loc[(stabch.Phenotype == "frameshift_variant") & 
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.Mut_in_lastexon == True) & (stabch.E3 == E3)
        & (stabch.gene == gene) & (stabch.degron_start == start) & (stabch.degron_end == end)]
        subsets_muttype_dict["frameshift_muts"] = process_muttype_df(frameshift_muts, cols_for_drop, "frameshift_muts", check_normality = check_normality, for_one_degron = True)
        
        
    if "frameshift_inbf_muts" in conditions:
        frameshift_inbf_muts = stabch.loc[(stabch.Phenotype == "frameshift_variant") & ((stabch.Loc_mut_degron == "inside") |
        (stabch.Loc_mut_degron == "outside_before")) &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.Mut_in_lastexon == True) & (stabch.E3 == E3)
        & (stabch.gene == gene) & (stabch.degron_start == start) & (stabch.degron_end == end)]
        subsets_muttype_dict["frameshift_inbf_muts"] = process_muttype_df(frameshift_inbf_muts, cols_for_drop, "frameshift_inbf_muts", check_normality = check_normality, for_one_degron = True)
        
            
    if "frameshift_aft_muts" in conditions:
        frameshift_aft_muts = stabch.loc[(stabch.Phenotype == "frameshift_variant") & (stabch.Loc_mut_degron == "outside_after") &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.Mut_in_lastexon == True) & (stabch.E3 == E3)
        & (stabch.gene == gene) & (stabch.degron_start == start) & (stabch.degron_end == end)]
        subsets_muttype_dict["frameshift_aft_muts"] = process_muttype_df(frameshift_aft_muts, cols_for_drop, "frameshift_aft_muts", check_normality = check_normality, for_one_degron = True)
        
            
    # Non-truncanting variants: missense and inframe
    if "nontrunc_muts" in conditions:
        nontrunc_muts = stabch.loc[((stabch.Phenotype == "missense_variant") |
        (stabch.Phenotype.str.contains("inframe"))) & (stabch.Altered_E3_Ligases == Altered_E3_Ligases) 
        & (stabch.E3 == E3) & (stabch.gene == gene) & (stabch.degron_start == start) & (stabch.degron_end == end)]
        subsets_muttype_dict["nontrunc_muts"] = process_muttype_df(nontrunc_muts, cols_for_drop, "nontrunc_muts", check_normality = check_normality, for_one_degron = True)
        
            
    if "nontrunc_in_muts" in conditions:
        nontrunc_in_muts = stabch.loc[((stabch.Phenotype == "missense_variant") |
        (stabch.Phenotype.str.contains("inframe"))) & (stabch.Loc_mut_degron == "inside") &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.E3 == E3) & (stabch.gene == gene) & (stabch.degron_start == start) & (stabch.degron_end == end)]
        subsets_muttype_dict["nontrunc_in_muts"] = process_muttype_df(nontrunc_in_muts, cols_for_drop, "nontrunc_in_muts", check_normality = check_normality, for_one_degron = True)
        
            
    if "nontrunc_out_muts" in conditions:
        nontrunc_out_muts = stabch.loc[((stabch.Phenotype == "missense_variant") |
        (stabch.Phenotype.str.contains("inframe"))) & (stabch.Loc_mut_degron.str.contains("outside")) &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.E3 == E3) & (stabch.gene == gene) & (stabch.degron_start == start) & (stabch.degron_end == end)]
        subsets_muttype_dict["nontrunc_out_muts"] = process_muttype_df(nontrunc_out_muts, cols_for_drop, "nontrunc_out_muts", check_normality = check_normality, for_one_degron = True)
        
            
    # Truncating variants (last exon): nonsense and frameshift   
    if "trunc_muts" in conditions:
        trunc_muts = stabch.loc[((stabch.Phenotype == "stop_gained") |
        (stabch.Phenotype == "frameshift_variant")) & (stabch.Altered_E3_Ligases == Altered_E3_Ligases)
        & (stabch.E3 == E3) & (stabch.gene == gene) & (stabch.degron_start == start) & (stabch.degron_end == end)
        & (stabch.Mut_in_lastexon == True)]
        subsets_muttype_dict["trunc_muts"] = process_muttype_df(trunc_muts, cols_for_drop, "trunc_muts", check_normality = check_normality, for_one_degron = True)
        
            
    if "trunc_inbf_muts" in conditions:
        trunc_inbf_muts = stabch.loc[((stabch.Phenotype == "stop_gained") |
        (stabch.Phenotype == "frameshift_variant")) & ((stabch.Loc_mut_degron == "inside") |
        (stabch.Loc_mut_degron == "outside_before")) & (stabch.Altered_E3_Ligases == Altered_E3_Ligases)
        & (stabch.E3 == E3) & (stabch.gene == gene) & (stabch.degron_start == start) & (stabch.degron_end == end)
        & (stabch.Mut_in_lastexon == True)]
        subsets_muttype_dict["trunc_inbf_muts"] = process_muttype_df(trunc_inbf_muts, cols_for_drop, "trunc_inbf_muts", check_normality = check_normality, for_one_degron = True)
        
            
    if "trunc_aft_muts" in conditions:
        trunc_aft_muts = stabch.loc[((stabch.Phenotype == "stop_gained") |
        (stabch.Phenotype == "frameshift_variant")) & (stabch.Loc_mut_degron == "outside_after") &
        (stabch.Altered_E3_Ligases == Altered_E3_Ligases) & (stabch.E3 == E3) & (stabch.gene == gene) & 
        (stabch.degron_start == start) & (stabch.degron_end == end) & (stabch.Mut_in_lastexon == True)]
        subsets_muttype_dict["trunc_aft_muts"] = process_muttype_df(trunc_aft_muts, cols_for_drop, "trunc_aft_muts", check_normality = check_normality, for_one_degron = True)


    # remove empty conditions
    subsets_muttype_dict_f = {k:v for (k,v) in subsets_muttype_dict.items() if not v.empty}

    return subsets_muttype_dict_f

def concat_subsets(subsets_dict, conditions, cols_for_drop, drop_duplicates = True,
                            filter_wts = False):
    """
    To filter WT condition (e.g.: E3 ligase specific dfs)
    """
    
    # Obtain needed mut type subsets
    subsets_dfs = []
    for cond in conditions:
        subsets_dfs.append(subsets_dict[cond])
    
    # Concatenate
    subsets_dfs_concat = pd.concat(subsets_dfs, ignore_index = True)
    
    if drop_duplicates:
        # Drop duplicates always keeping inside condition (if applicable)
        subsets_dfs_concat = subsets_dfs_concat.loc[(subsets_dfs_concat.duplicated(
            subset = cols_for_drop, keep = False) == False) | 
            ((subsets_dfs_concat['Loc_mut_degron'] == 'inside') |
             ((subsets_dfs_concat['Loc_mut_degron'] == 'outside_before') &
             ((subsets_dfs_concat['Phenotype'] == 'stop_gained') | 
              (subsets_dfs_concat["Phenotype"] == "frameshift_variant"))))]
        # Drop duplicates again to correct for several degrons affected by the same inside/outside_before mutation
        subsets_dfs_concat = subsets_dfs_concat.drop_duplicates(subset = cols_for_drop)
    
    if filter_wts:
        # Only keep WT conditions of proteins having a mutation
        mutated_genes = subsets_dfs_concat.loc[subsets_dfs_concat.Phenotype != "WT", "gene"].unique().tolist()
        subsets_dfs_concat = subsets_dfs_concat.loc[subsets_dfs_concat.gene.isin(mutated_genes)]
    
    return subsets_dfs_concat