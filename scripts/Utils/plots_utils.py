# -- Import libraries -- # 

import sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors, colorbar
from matplotlib.lines import Line2D
import matplotlib.patches as patches
import seaborn as sns
import numpy as np
from numpy.polynomial.polynomial import polyfit
import random
from copy import deepcopy
from statannotations.Annotator import Annotator
from copy import deepcopy

## my modules ##
sys.path.append("./")    # modules folder
from motif_utils import motif_scan_flex
from stabch_utils import concat_subsets

# -- Functions -- #

## -- GENERAL PLOT UTILS -- ##

# Palette according to y values
#https://stackoverflow.com/questions/36271302/changing-color-scale-in-seaborn-bar-plot
def colors_from_values(values, palette_name, n_remove_outliers = 0):
    """
    Modifies a plot palette based on an array of values to match palette degradation with 
    them. If necessary, does not consider n outliers.
    
    Parameters
    ----------
    values: numpy array 
           Array of integers or floats to calculate the palette from
    palette_name: str
           Valid matplotlib palette name
    n_remove_outliers: int (default: 0)
           Number of outliers to remove for building the palette degradation scheme. 
           Note: currently only implemented for max 1 outlier.
    
    Returns
    -------
    d: numpy array
           Array with the palette indices to consider in the plot
    """
    # normalize the values to range [0, 1]
    ## consider outliers so that there is a consistent degradation in the palette
    ## currently, this will only work for 1 outlier
    for i in range(n_remove_outliers):
        outlier_idx = np.where(values == np.amax(values))     # max value = outlier
        values_wout_outlier = np.delete(values, outlier_idx)
        np.put(values, outlier_idx, np.amax(values_wout_outlier))
    normalized = (values - min(values)) / (max(values) - min(values))
        
    # convert to indices
    indices = np.round(normalized * (len(values) - 1)).astype(np.int32)
        
    # use the indices to get the colors
    palette = sns.color_palette(palette_name, len(values))
    
    return np.array(palette).take(indices, axis = 0)

def create_all_pairs(conditions):
    """
    Generates all possible conditions pairs for statistical comparison

    Parameters
    ----------
    conditions: list
            List with all possible conditions. No needs to be sorted in any 
            specific manner.
    
    Returns
    -------
    pairs: list of tuples
            List with tuples representing all possible conditions to be compared
    """
    
    pairs = []
    i = 0  
    
    for cond1 in conditions:
        i += 1    # to exclude comparison with itself
        for cond2 in conditions[i:]:
            
            pairs.append((cond1, cond2))
    
    return pairs

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


## -- QUALITY ANALYSIS PLOTS -- ##

def barplot_qa_property(qa_df, property, ylab_property, plot_file, font_y = 15, font_x = 15, n_remove_outliers = 0, 
                        fig_width = 20, fig_height = 6):
    """
    Barplot to display a QA calculated property per motif, using a color palette based on property values.

    Parameters
    ----------
    qa_df: pandas dataframe
            Pandas dataframe of the quality analysis metrics.
    property: str
            Column name of the property to be plotted in qa_df
    ylab_property: str
            Y label for the barplot
    plot_file: str
            Path to the folder where to store the resulting plot
    font_y: int (default: 15)
            Y label font size
    font_x: int (default: 15)
            X label font size
    n_remove_outliers: int (default: 0)
            Number of outliers to remove for building the palette degradation scheme. 
            Note: currently only implemented for max 1 outlier. Passed to colors_from_values
            function
    fig_width: int (default: 20)
            Figure width
    fig_height: int (default: 6)
            Figure height
    
    Returns
    -------
    ax: matplotlib plot
            QA property barplot
            
    """
    
    # Generate property df and subsets
    property_df = qa_df[["motif_id", property]].sort_values(by = ["motif_id"])
    x = property_df["motif_id"].values
    y = property_df[property].values
    y_pal = deepcopy(y)  # to avoid changing the y-array for the barplot when handling color outliers
    
    # Plot settings
    plt.figure(figsize = (fig_width, fig_height))
    sns.set_style("white")
    ## Define palette by y values
    palette = colors_from_values(y_pal, "YlOrRd", n_remove_outliers)

    # Barplot
    ax = sns.barplot(x = x, y = y, palette = palette)
    ax.set_ylabel(ylab_property, fontsize = font_y)
    ax.set_xlabel("Degron motif", fontsize = font_x)
    ax.set_xticklabels(x, rotation = 45, ha = "right")
    ax.tick_params(labelsize = 13)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.yticks(fontsize = 15)

    # Add colorbar legend (customized)
    cmap = plt.cm.get_cmap("YlOrRd")
    sm = plt.cm.ScalarMappable(cmap = cmap, norm = plt.Normalize(qa_df[property].min(), qa_df[property].max()))
    sm.set_array([])
    cbar = plt.colorbar(sm)
    cbar.ax.tick_params(labelsize = 15)

    # Save
    plt.savefig(plot_file, dpi = 800, transparent = True, bbox_inches = "tight")
    
    return ax

# https://matplotlib.org/stable/gallery/subplots_axes_and_figures/broken_axis.html
def barplot_broken_axis_qa_property(qa_df, property, ylab_property, y_lim_high, y_lim_low, 
                                    low_pad_x, high_pad_x, low_pad_y, high_pad_y, plot_file, font_y = 18, 
                                    font_x = 18, property_percentage = None, n_remove_outliers = 0,
                                    fig_width = 20, fig_height = 8):
    """
    Barplot to display a QA calculated property per motif, using a color palette based on property values.
    Includes broken axis to account for very big outliers.

    Parameters
    ----------
    qa_df: pandas dataframe
            Pandas dataframe of the quality analysis metrics.
    property: str
            Column name of the property to be plotted in qa_df
    ylab_property: str
            Y label for the barplot
    y_lim_high: int
            Y-axis limit for the top part of the barplot
    y_lim_low: int
            Y-axis limit for the bottom part of the barplot
    low_pad_x: int
            X-axis padding in the bottom part of the barplot
            for integer or percentages annotation over the bars
    high_pad_x: int
            X-axis padding in the top part of the barplot
            for integer or percentages annotation over the bars
    low_pad_y: int
            Y-axis padding in the bottom part of the barplot
            for integer or percentages annotation over the bars
    high_pad_y: int
            Y-axis padding in the top part of the barplot
            for integer or percentages annotation over the bars
    plot_file: str
            Path to the folder where to store the resulting plot
    font_y: int (default: 18)
            Y label font size
    font_x: int (default: 18)
            X label font size
    property_percentage: str (default: None)
            If None, absolute y values are added over the bars
            If str: str has to be the column name of the property_percentage
            in the qa_df. Percentage y values are then added over the bars.
    n_remove_outliers: int (default: 0)
            Number of outliers to remove for building the palette degradation scheme. 
            Note: currently only implemented for max 1 outlier. Passed to colors_from_values
            function
    fig_width: int (default: 20)
            Figure width
    fig_height: int (default: 8)
            Figure height
    
    Returns
    -------
    ax1, ax2: matplotlib plots
            QA top and bottom property barplot
            
    """
    
    # Generate property df and subsets
    ## To plot percentages over the bars
    if property_percentage == None:  
        property_df = qa_df[["motif_id", property]].sort_values(by = ["motif_id"])
    ## To plot absolute values over the bars
    else:
        property_df = qa_df[["motif_id", property, property_percentage]].sort_values(by = ["motif_id"])
        
    x = property_df["motif_id"].values
    y = property_df[property].values
    y_pal = deepcopy(y)   # to avoid changing the y-array for the barplot when handling color outliers
    
    
    # Plot settings for broken axis
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True)   # generate 2 subplots
    sns.set_style("white")
    palette = colors_from_values(y_pal, "YlOrRd", n_remove_outliers)
    fig.set_size_inches((fig_width, fig_height))
    fig.subplots_adjust(hspace = 0.05)  # adjust space between axes
    
    # Barplot: plot the same data on both axes
    sns.barplot(x = x, y = y, palette = palette, ax = ax1)
    sns.barplot(x = x, y = y, palette = palette, ax = ax2)
    
    ## zoom-in / limit the view to different portions of the data
    ax1.set_ylim(y_lim_high)  # high levels (top part of the plot)
    ax2.set_ylim(y_lim_low)  # low levels (bottom part of the plot)
    
    ## hide the spines between ax and ax2
    ax1.spines.bottom.set_visible(False)
    ax2.spines.top.set_visible(False)
    ax1.spines.top.set_visible(False)
    ax1.tick_params(bottom = False, top = False)
    ax1.tick_params(labeltop = False)  # don't put tick labels at the top
    ax2.xaxis.tick_bottom()
    ax2.spines.right.set_visible(False)
    ax1.spines.right.set_visible(False)
    
    ## generate a big subplot for y label
    ax = fig.add_subplot(111, frameon = False)
    ### hide tick and tick label of the big axis
    ax.tick_params(labelcolor = 'none', which = 'both', top = False, bottom = False, left = False, right = False)
    ax.set_ylabel(ylab_property, fontsize = font_y)
    ax2.set_xlabel("Degron motif", fontsize = font_x)
    ax2.set_xticklabels(x, rotation = 45, ha = "right")
    ax2.tick_params(labelsize = 15) 
    ax1.tick_params(labelsize = 15)
    plt.yticks(fontsize = 15)
    
    ## generate axis break
    d = .5  # proportion of vertical to horizontal extent of the slanted line
    kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
              linestyle="none", color='k', mec='k', mew=1, clip_on=False)
    ax1.plot([0, 1], [0, 0], transform = ax1.transAxes, **kwargs)
    ax2.plot([0, 1], [1, 1], transform = ax2.transAxes, **kwargs)
    
    # Labels over the bars
    ## label with percentages
    if property_percentage != None:
        percents = property_df[property_percentage].values
        # Top plot
        for p, percent in zip(ax1.patches, percents):
            if p.get_height() >= y_lim_low[1]:
                ax1.annotate(str(round(percent, 3)), xy = (p.get_x() + high_pad_x, 
                            p.get_height() + high_pad_y), fontsize = 13)
        # Bottom plot
        for p, percent in zip(ax2.patches, percents):
            if p.get_height() < y_lim_low[1]:
                ax2.annotate(str(round(percent, 3)), xy = (p.get_x() + low_pad_x,
                            p.get_height() + low_pad_y), fontsize = 13)
    
    ## label with absolute values
    else:
        # Top plot
        for p in ax1.patches:
            if p.get_height() >= y_lim_low[1]:
                ax1.annotate(str(round(p.get_height())), xy = (p.get_x() + high_pad_x, 
                            p.get_height() + high_pad_y), fontsize = 13)
        # Bottom plot
        for p in ax2.patches:
            if p.get_height() < y_lim_low[1]:
                ax2.annotate(str(round(p.get_height())), xy = (p.get_x() + low_pad_x,
                            p.get_height() + low_pad_y), fontsize = 13)

    # Add colorbar legend (customized)
    cmap = plt.cm.get_cmap("YlOrRd")
    sm = plt.cm.ScalarMappable(cmap = cmap, norm = plt.Normalize(qa_df[property].min(), qa_df[property].max()))
    sm.set_array([])
    cbaxes = fig.add_axes([0.95, 0.1, 0.01, 0.8])
    cbar = plt.colorbar(sm, cax = cbaxes)
    cbar.ax.tick_params(labelsize = 15)

    # Save
    plt.savefig(plot_file, dpi = 800, transparent = True, bbox_inches = "tight")

    return ax1, ax2

def violinplot_pos_threshold(pwm, pos_set, nopos_scan, E3, proteome, plot_file, n_samples = 100, negatives_to_zero = False, 
                            fig_width = 5, fig_height = 10):
    """
    Violinplot to display the different distribution in a PWM scores over the set of true degrons and a random set of protein
    subsequences. Also, shows the positivity threshold. 

    Parameters
    ----------
    pwm: pandas dataframe
            Dataframe containing a motif's weight matrix, whose rows correspond to motif's positions 
            and columns to aminoacids. Each cell contains the weight associated to each aminoacid in
            each position of the motif
    pos_set: pandas dataframe
            Dataframe containing the set of true degrons for a specific motif. Columns: gene, sequence,
            start, end
    nopos_scan: dict
            Dictionary of the form {protein_ID: [scores]}. Contains the protein-associated scores of a
            motif scan (without true degrons)
    E3: str
            E3 ligase motif's name
    proteome: dict
            Dictionary of the form {protein_ID: sequence}. Contains a set of protein sequences (e.g.: proteome)
    plot_file: str
            Path to the folder where to store the resulting plot
    n_samples: int (default: 100)
            Number of random subsequences to include in the random distribution
    negatives_to_zero: boolean (default: False)
            If True, negative scores are transformed to zero (for visualization purposes)
            If False, negative scores are maintained as negative
    fig_width: int (default: 5)
            Figure width
    fig_height: int (default: 10)
            Figure height
    
    Returns
    -------
    ax: matplotlib plot
            Violinplot of true degrons vs random sequences scores distribution   
    """

    # Calculate scores of positive set
    pos_set_scores = []
    for row in pos_set.itertuples(): 
        pos_set_scores.append(float(str(motif_scan_flex(row.gene, row.sequence, row.start, row.end, 
                                                        pwm, proteome))[1:-1]))

    # Collect random sequences from non-positive set
    nopos_scan_valid_keys = [k for k in nopos_scan if nopos_scan[k] != []]  # some proteins have no scores
    random_proteins = random.sample(nopos_scan_valid_keys, n_samples)
    random_scores = [random.choice(nopos_scan[rand_prot]) for rand_prot in random_proteins]
    
    # Generate df
    pos_set_scores_df = pd.DataFrame(pos_set_scores, columns = ["Score"])
    pos_set_scores_df["Set"] = "positive"
    random_scores_df = pd.DataFrame(random_scores, columns = ["Score"])
    random_scores_df["Set"] = "random"
    scores = pd.concat([pos_set_scores_df, random_scores_df])
    # Correct negatives if necessary
    if negatives_to_zero:
        scores.loc[scores.Score < 0] = 0
    
    # Plot settings
    palette1 = {"positive": "#a1d99b", "random": "#ef6548"}
    palette2 = {"positive": "#006d2c", "random": "#b30000"}
    axis_names = ["True\ndegrons\nN = "+str(len(pos_set_scores_df)),
                  "Random\nsequences\nN = "+str(len(random_scores_df))]
    conditions = ["positive", "random"]
    sns.set_style(style = 'whitegrid')
    plt.figure(figsize = (fig_width, fig_height))

    # Violinplot (+ stripplot)
    ax = sns.violinplot(x = "Set", y = "Score", data = scores,
                        width = 0.7, saturation = 0.7, palette = palette1, fliersize = 0., 
                        linewidth = 1, order = conditions)
    ax = sns.stripplot(x = "Set", y = "Score", data = scores,
                        alpha = 0.7, size = 8, palette = palette2,
                        jitter = True, dodge = True, order = conditions)
    
    ax.set_xticks(np.arange(len(axis_names)), axis_names, fontsize = 12)
    plt.ylabel("PWM score", fontsize = 15)
    plt.xlabel("")
    plt.yticks(fontsize = 15)
    plt.title(E3, size = 18, pad = 10)

    # Plot positivity threshold
    pos_scores = [s for s in pos_set_scores_df.Score.values.tolist() if s >= 0]
    pos_threshold = min(pos_scores)
    plt.axhline(y = pos_threshold, color = 'r', linestyle = '--', alpha = 0.5)

    # Save
    plt.savefig(plot_file, dpi = 800, transparent = True, bbox_inches = "tight")
    
    return ax

def hits_length_scatterplot(qa_df, corr_df, E3, plot_file, font_y = 18, font_x = 18, fig_width = 10, fig_height = 7):
    """
    Scatterplot to display a PWMs relation between the number of degrons discovered in each protein sequence vs the protein's
    length

    Parameters
    ----------
    qa_df: pandas dataframe
            Pandas dataframe of the quality analysis metrics.
    corr_df: pandas dataframe
            Dataframe containing the number of discovered degrons per protein per motif and each protein's
            length
    E3: str
            E3 ligase motif's name
    plot_file: str
            Path to the folder where to store the resulting plot
    font_y: int (default: 18)
            Y label font size
    font_x: int (default: 18)
            X label font size
    fig_width: int (default: 10)
            Figure width
    fig_height: int (default: 7)
            Figure height
    
    Returns
    -------
    ax: matplotlib plot
            Violinplot of true degrons vs random sequences scores distribution   
    """
    
    # Motif subset
    y = corr_df["Protein length"].values
    x = corr_df[E3].values

    # Plot settings
    palette = ["red"]
    plt.figure(figsize = (fig_width, fig_height))
    sns.set_style("whitegrid")

    # Scatterplot with regression line
    sns.regplot(x = x, y = y, scatter_kws = {"color": "red"}, line_kws = {"color": "grey"})
    r_value = qa_df.loc[qa_df.motif_id == E3, "correlation_hits_vs_length_sequence"]   # display correlation index
    plt.title(E3 + "\nR-value = %0.2f" % (r_value) + f"   N = "+str(len(x)), size = 16, pad = 10)
    plt.xlabel("Number of discovered degrons", fontsize = font_y)
    plt.ylabel("Protein length", fontsize = font_x)
    plt.yticks(fontsize = 15)
    plt.xticks(fontsize = 15)
    
    # Save
    plt.savefig(plot_file, dpi = 800, transparent = True, bbox_inches = "tight")

    return None


## -- DATASET DESCRIPTION PLOTS -- ##

def dataset_piechart(stabch_df, samples_col, ctype_col, dataset, palette, plot_file, min_sample = 15, fig_width = 10, fig_height = 7):
    """
    Piechart to display a dataset composition.

    Parameters
    ----------
    stabch_df: pandas dataframe
            Pandas dataframe with the information on the mutations and stability change levels of several proteins. CPTAC or CCLE dataset.
    samples_col: str
            Column name in stabch_df for the sample ID
    ctype_col: str
            Column name in stabch_df for the cancer type
    dataset: str
            Dataset name
    palette: dict
            Dictionary of the form {cancer_type: color} to map each cancer type to a specific color.
    plot_file: str
            Path to the folder where to store the resulting plot
    min_sample: int (default: 15)
            Minimum number of samples in a cancer type to consider such cancer type independendently. If less,
            merged it in "OTHER" category
    fig_width: int (default: 10)
            Figure width
    fig_height: int (default: 7)
            Figure height
    
    Returns
    -------
    None 
    """
    
    # Sample IDs and cancer type subset
    subset = stabch_df[[samples_col, ctype_col]]
    subset = subset.drop_duplicates()
    
    # Group by cancer type and count number of samples
    subset_gpby = subset.groupby(ctype_col, as_index = False).count().sort_values(samples_col, ascending = False)
    
    # In case there are groups < min_sample samples, pull them together
    subset_gpby.loc[subset_gpby[samples_col] < min_sample, ctype_col] = "OTHER"
    subset_gpby = subset_gpby.groupby(subset_gpby[ctype_col]).aggregate({samples_col: 'sum'}).reset_index()
    sizes = subset_gpby[samples_col].tolist()
    labels = subset_gpby[ctype_col].tolist()
    if "cohort" in labels:
        labels[labels.index("cohort")] = "PDAC"  # manually include PDAC
    if "HAEMATOPOIETIC AND LYMPHOID TISSUE" in labels:
        labels[labels.index("HAEMATOPOIETIC AND LYMPHOID TISSUE")] = "HAEMATOPOIETIC AND\nLYMPHOID TISSUE"  # better visualization
    
    # Pieplot 
    def make_autopct(values):
        def my_autopct(pct):
            total = sum(values)
            val = int(round(pct*total/100.0))
            return '{v:d}'.format(v = val)
        return my_autopct 
    
    plt.figure(figsize = (fig_width, fig_height))
    
    plt.pie(sizes, labels = labels, autopct = make_autopct(sizes), shadow = False, colors = palette)
    plt.title(dataset, size = 15)
    total = 'Total = '+str(sum(sizes))
    plt.text(-0.5, -1.4, total, fontsize = 12)

    plt.savefig(plot_file, dpi = 800, transparent = True, bbox_inches = "tight")

    return None

## -- STABILITY CHANGE PLOTS -- ##

def stabch_plot(conditions, pairs, subsets_dict, dataset, stabch, ylim1, ylim2, palette, plot_file,
                degron_features = "all", main_plot = "boxplot",
                do_stats = True, stripplot = True, E3 = None, custom_title = None, annot_n = True,
                pad_stats = 130, fig_width = 10, fig_height = 12):
    """
    Boxplot to compare stability change levels between conditions.

    Parameters
    ----------
    conditions: list
            Conditions to compare.
    pairs: list of tuples
            Conditions to compare with a statistical test
    subsets_dict: dict
            Dictionary of the form {condition: dataframe}. Contains all the samples from a mutations dataframe
            corresponding to each condition already filtered. 
    dataset: str
            Dataset name
    stabch: pandas dataframe
            Dataframe containing mutations and stability change levels of several proteins.
    y_lim1: int
            Y-axis top limit
    y_lim2: int
            Y-axis bottom limit
    palette: dict
            Dictionary of the form {condition: color} to map each condition to a specific color.
    plot_file: str
            Path to the folder where to store the resulting plot
    degron_features: str (default: "all")
            If all: the comparison is for the entire dataset
            If E3_ligase: the comparison is for a specific motif
            If degron_instance: the comparison is for a specific degron
    main_plot: str (default: "boxplot")
            If boxplot: the main plot is a seaborn boxplot
            If violinplot: the main plot is a violinplot
    do_stats: boolean (default: True)
            If True, a Mann Whitney test will be computed between pairs
            If False, no statistical test is computed
    stripplot: boolean (default: True)
            If True, a seaborn stripplot will be also plotted
    E3: str (default: None)
            If None, degron_features must be "all"
            If str, indicated the E3 ligase for which the conditions are being compared
    custom_title: str (default: None)
            If None, degron_features must be "all" or "E3_ligase"
            If str, the plot title 
    annot_n: boolean (default: True)
            If True, sample size is annotated below each condition
    pad_stats: int (default: 130)
            The pad between the statistical annotation and the plot's title. Only considered
            when do_stats = True
    fig_width: int (default: 10)
            Figure width
    fig_height: int (default: 12)
            Figure height
    
    Returns
    -------
    ax: matplotlib plot
            Stability change boxplot comparing selected conditions 
    """
    
    # To drop duplicates
    cols_to_remove = ["E3", "degron_start", "degron_end", "Loc_mut_degron"]
    cols = stabch.columns.tolist()
    cols_for_drop = [col for col in cols if col not in cols_to_remove]
    
    # Axis names dictionaries for all possible conditions
    axis_names = {"wt": "WT",
                  "syn_muts": "Synonymous",
                  "miss_muts": "Missense",
                  "inframe_muts": "Inframe",
                  "nonsense_muts": "Nonsense",
                  "frameshift_muts": "Frameshift",
                  "nontrunc_muts": "Non-truncating",
                  "trunc_muts": "Truncanting",
                  "syn_in_muts": "Synonymous\n inside degron",
                  "syn_out_muts": "Synonymous\n outside degron", 
                  "miss_in_muts": "Missense\n inside degron", 
                  "miss_out_muts": "Missense\n outside degron", 
                  "inframe_in_muts": "Inframe\n inside degron",
                  "inframe_out_muts": "Inframe\n outside degron", 
                  "nonsense_inbf_muts": "Nonsense\n inside/before degron",
                  "nonsense_aft_muts": "Nonsense\n after degron", 
                  "frameshift_inbf_muts": "Frameshift\n inside/before degron",
                  "frameshift_aft_muts":  "Frameshift\n after degron", 
                  "nontrunc_in_muts": "Non-truncating\n inside degron", 
                  "nontrunc_out_muts": "Non-truncating\n outside degron",
                  "trunc_inbf_muts": "Truncating\ninside/before\ndegron", 
                  "trunc_aft_muts": "Truncating\nafter\ndegron"}
    
    # Select conditions and generate subsets
    ## axis and palette
    conditions_axis_names = []
    conditions_palette = {}
    for cond in conditions:
        conditions_axis_names.append(axis_names[cond])
        conditions_palette[cond] = palette[cond]

    ## data subsets: concatenate in unique dataframe
    if degron_features == "all":
        subset = concat_subsets(subsets_dict, conditions, cols_for_drop)
    elif degron_features == "E3_ligase":
        subset = concat_subsets(subsets_dict, conditions, cols_for_drop, filter_wts = True)
    elif degron_features == "degron_instance":
        subset = concat_subsets(subsets_dict, conditions, cols_for_drop,  drop_duplicates = False)
    
    # Edit axis conditions to display N
    if annot_n:
        conditions_axis_names_n = []
        for axis_cond, cond in zip(conditions_axis_names, conditions):
            conditions_axis_names_n.append(axis_cond+"\nN = "+str(len(subset.loc[subset.Condition == cond])))
        conditions_axis_names = conditions_axis_names_n
    
    # Plot settings
    sns.set_style(style = 'whitegrid')
    plt.figure(figsize = (fig_width, fig_height)) 
    plt.ylim([ylim1, ylim2])
    pad = 10
    title_size = 20
    if do_stats:
        pad = pad_stats    
    if degron_features == "all":
        plt.title(dataset+' degrons', pad = pad, size = title_size)
    elif degron_features == "E3_ligase":
        plt.title(dataset+' '+E3+' degrons', pad = pad, size = title_size)
    elif degron_features == "degron_instance":
        plt.title(custom_title, pad = pad, size = title_size)
    
    # Stripplot
    if stripplot:
        ax = sns.stripplot(x = "Condition", y = "Stability_Change", data = subset,
                        alpha = 0.3, size = 4, palette = conditions_palette,
                        jitter = True, dodge = True, order = conditions)
    # Main plot: boxplot or violinplot
    if main_plot == "violinplot":
        ax = sns.violinplot(x = "Condition", y = "Stability_Change", data = subset,
                        width = 0.7, saturation = 0.7, palette = conditions_palette, fliersize = 0., 
                        linewidth = 1, order = conditions)
    elif main_plot == "boxplot":
        ax = sns.boxplot(x = "Condition", y = "Stability_Change", data = subset,
                        width = 0.7, saturation = 0.7, palette = conditions_palette, fliersize = 0., 
                        linewidth = 1, order = conditions)
        
    plt.xticks(np.arange(len(conditions_axis_names)), conditions_axis_names, fontsize = 13)
    plt.ylabel("Stability change", fontsize = 18)
    plt.xlabel("Mutation type", fontsize = 18)
        
    # Perform Mann Whitney test and annotate p values
    if do_stats:
        annotator = Annotator(ax, pairs, data = subset, x = "Condition", y = 'Stability_Change',
                              order = conditions)
        annotator.configure(text_format = "star", loc = 'outside', test = 'Mann-Whitney').apply_and_annotate()
    
    # Save
    plt.savefig(plot_file, dpi = 800, transparent = True, bbox_inches = "tight")

    return ax

def rna_vs_protein_scatterplot(stabch, gene, ctype, dataset, regresline_lim_left, regresline_lim_right, plot_file,
                                fig_width = 10, fig_height = 10):
    """
    Scatterplot to compare RNA vs protein levels of WT forms and altered forms and display the raw residual
    of the latter (i.e.: stability change). Comparison for a specific protein, cancer type and dataset.

    Parameters
    ----------
    stabch: pandas dataframe
            Dataframe containing mutations and stability change levels of several proteins. CPTAC only, specific 
            cancer type (no pancancer). Before degrons annotation (does not correct for duplicates)
    ctype: str
            Cancer type name
    gene: str
            Gene name
    dataset: str
            Dataset name. CPTAC only.
    regresline_lim_left: int
            Left limit of the regression line to be displayed in the plot
    regresline_lim_right: int
            Right limit of the regression line to be displayed in the plot
    plot_file: str
            Path to the folder where to store the resulting plot
    fig_width: int (default: 10)
            Figure width
    fig_height: int (default: 10)
            Figure height
    
    Returns
    -------
    None
    """
    
    # Gene subset
    subset = stabch.loc[(stabch["gene"] == gene)]
    ## wt
    x_wt = subset.loc[subset["Phenotype"] == "WT", "log2fpkm"].values
    y_wt = subset.loc[subset["Phenotype"] == "WT", "norm_protein_expression"].values
    ## non-truncanting (missense and inframe)
    phen_nontrunc = ["missense_variant", "inframe_insertion", "inframe_deletion"]
    x_nontrunc = subset.loc[subset["Phenotype"].isin(phen_nontrunc), "log2fpkm"].values
    y_nontrunc = subset.loc[subset["Phenotype"].isin(phen_nontrunc), "norm_protein_expression"].values
    ## truncanting
    phen_trunc = ["frameshift_variant", "stop_gained"]
    x_trunc = subset.loc[subset["Phenotype"].isin(phen_trunc), "log2fpkm"].values
    y_trunc = subset.loc[subset["Phenotype"].isin(phen_trunc), "norm_protein_expression"].values
    
    # Plot settings
    fig, ax = plt.subplots()
    fig.set_size_inches(fig_width, fig_height)
    sns.set_style("whitegrid")
    
    # Scatterplot
    ax.scatter(x_wt, y_wt, c = "grey")
    ax.scatter(x_nontrunc, y_nontrunc, c = "#fb9a99")
    ax.scatter(x_trunc, y_trunc, c = "#fdbf6f")
    plt.yticks(fontsize = 15)
    plt.xticks(fontsize = 15)
    
    # Regression line for WTs
    b, m = polyfit(x_wt, y_wt, deg = 1)
    ax.plot(np.append(x_wt, [regresline_lim_left, regresline_lim_right]), 
            b+m*np.append(x_wt, [regresline_lim_left, regresline_lim_right]), "-") # add 0 to increase lines length
    ## predicted protein expression according to linear model
    subset["pred_protein_expression"] = subset.apply(
        lambda row: b+m*row["log2fpkm"], axis = 1)
    
    # Plot stability change (predicted vs observed protein levels)
    ## non-truncanting
    nontrunc_df = subset.loc[subset["Phenotype"].isin(phen_nontrunc)]
    list_df = nontrunc_df[['log2fpkm', 'norm_protein_expression', 'pred_protein_expression']].values.tolist()
    x = []
    y = []
    for item in list_df:
        x.append([item[0], item[0]])
        y.append([item[1], item[2]])

    for i in range(0, len(x)):
        ax.plot(x[i], y[i], '--', color = "#fb9a99", linewidth = .7)  
    ## truncanting
    trunc_df = subset.loc[subset["Phenotype"].isin(phen_trunc)]
    list_df = trunc_df[['log2fpkm', 'norm_protein_expression', 'pred_protein_expression']].values.tolist()
    x = []
    y = []
    for item in list_df:
        x.append([item[0], item[0]])
        y.append([item[1], item[2]])

    for i in range(0, len(x)):
        ax.plot(x[i], y[i], '--', color = "#fdbf6f", linewidth = .7)

    ax.set_title(gene+" ("+ctype+" tumors from "+dataset+")\n", fontsize = 20)
    ax.set_xlabel("mRNA (log2fpkm)", fontsize = 18)
    ax.set_ylabel("Protein expression", fontsize = 18)

    # Add customized legend
    legend_cols_dict = [Line2D([0], [0], color = 'w', markerfacecolor = '#fb9a99', marker = 'o',
                        markersize = 8, label = 'Missense or inframe'),
                        Line2D([0], [0] , color = 'w', markerfacecolor = '#fdbf6f', marker = 'o',
                        markersize = 8, label = 'Truncanting')]
    legend1 = ax.legend(handles = legend_cols_dict, loc = "upper left",
                         fontsize = 12, framealpha = 0.3)
    ax.add_artist(legend1)
    
    # Save
    plt.savefig(plot_file, dpi = 800, transparent = True, bbox_inches = "tight")

    return None

def stabch_needle_plot(protein_symb, proteome, stabch, dataset, degron_start, degron_end, plot_file, fig_width = 15, fig_height = 10):
    """
    Needle plot to display a protein's sequence together with the mutation count per position and the stability change
    per position

    Parameters
    ----------
    protein_symb: str
            Gene name
    proteome: dict
            Dictionary of the form {protein_ID: sequence}. Contains a set of protein sequences (e.g.: proteome)
    stabch: pandas dataframe
            Dataframe containing mutations and stability change levels of several proteins. 
    dataset: str
            Dataset name.
    degron_start: int
            Starting position of the degron in the sequence
    degron_end: int
            Ending position of the degron in the sequence
    plot_file: str
            Path to the folder where to store the resulting plot
    fig_width: int (default: 15)
            Figure width
    fig_height: int (default: 10)
            Figure height
    
    Returns
    -------
    None
    """
    
    # Extract protein sequence and length
    enst = stabch.loc[stabch.gene == protein_symb, "Feature"].unique()[0]
    protein_len = len(proteome[enst])
        
    # Filter mutations in protein of interest
    cols_to_keep = ["gene", "Feature", "protein_mutation", "Protein_position", "Phenotype",
                   "Stability_Change", "Altered_E3_Ligases", "Mut_in_lastexon"]
    subset = stabch[cols_to_keep].drop_duplicates()
    phens_trunc = ["stop_gained", "frameshift_variant"]
    phens_nontrunc = ["missense_variant", "inframe_deletion", "inframe_insertion", "synonymous_variant"]
    subset = subset.loc[(subset.gene == protein_symb) &
                       ((subset.Phenotype.isin(phens_nontrunc)) | ((subset.Phenotype.isin(phens_trunc)) &
                       (subset.Mut_in_lastexon == True))) & (subset.Altered_E3_Ligases == False)].drop_duplicates().sort_values(
                        by = ["protein_mutation"])
    
    # Simplify mut type to syn, non-truncating and truncating
    def simplify_conditions(row):
        if row.Phenotype == "synonymous_variant":
            row["Phenotype_reduced"] = "synonymous"
        elif row.Phenotype in ["missense_variant", "inframe_deletion", "inframe_insertion"]:
            row["Phenotype_reduced"] = "non_truncating"
        elif row.Phenotype in ["stop_gained", "frameshift_variant"]:
            row["Phenotype_reduced"] = "truncating"
            
        return row
    
    subset = subset.apply(lambda row: simplify_conditions(row), axis = 1)

    # Convert intervals into point mutation signals
    def break_intervals(row):

        if isinterval(row.Protein_position):
            interval = row.Protein_position.split("-")
            row["Protein_position"] = range(int(interval[0]), int(interval[1])+1)

        return row
    
    subset = subset.apply(lambda row: break_intervals(row), axis = 1)
    subset = subset.explode(["Protein_position"], ignore_index = True)

    # Transform protein position column to dtype int
    subset["Protein_position"] = pd.to_numeric(subset["Protein_position"])

    # Groupby mutation type and position and count freq
    subset_gpby = subset.groupby(["Phenotype_reduced", "Protein_position"]).size().reset_index(
    ).rename(columns = {0:'muts_count'})
    
    # Groupby mutation type and calculate median stability change
    subset_gpby_stabch = subset.groupby("Protein_position")["Stability_Change"].median().reset_index()
    
    # Create mutation freq dicts per mutation type and for stability change 
    y_syn = {}
    y_nontrunc = {}
    y_trunc = {}
    y_stabch = {}
    
    for pos in range(protein_len):
        
        pos = pos+1
        subset_gpby_pos = subset_gpby.loc[subset_gpby.Protein_position == pos]
        subset_gpby_stabch_pos = subset_gpby_stabch.loc[subset_gpby_stabch.Protein_position == pos]
        
        # For the needle plot
        if subset_gpby_pos.empty:
            y_syn[int(pos)] = np.nan
            y_nontrunc[int(pos)] = np.nan
            y_trunc[int(pos)] = np.nan
            
        else:
            for row in subset_gpby_pos.itertuples():
                if row.Phenotype_reduced == "synonymous":
                    y_syn[int(pos)] = row.muts_count
                elif row.Phenotype_reduced == "non_truncating":
                    y_nontrunc[int(pos)] = row.muts_count
                elif row.Phenotype_reduced == "truncating":
                    y_trunc[int(pos)] = row.muts_count
            
            if "synonymous" not in subset_gpby_pos.Phenotype_reduced.to_numpy():
                y_syn[int(pos)] = np.nan
            if "non_truncating" not in subset_gpby_pos.Phenotype_reduced.to_numpy():
                y_nontrunc[int(pos)] = np.nan
            if "truncating" not in subset_gpby_pos.Phenotype_reduced.to_numpy():
                y_trunc[int(pos)] = np.nan
            
        # For the stability change plot
        if subset_gpby_stabch_pos.empty:
            y_stabch[pos] = 0
        else:   
            for row in subset_gpby_stabch_pos.itertuples():
                y_stabch[int(pos)] = row.Stability_Change

            
    # Generate needle plot
    fig, (ax1, ax2, ax3) = plt.subplots(3, sharex = True, gridspec_kw = {'height_ratios': [1, 4, 0.25]})
    fig.set_size_inches(fig_width, fig_height)
    sns.set_style("white")
    
    ax2.stem(np.array(list(y_syn.keys())).astype(int), np.array(list(y_syn.values())),
             linefmt = "red", markerfmt = "ro", label = 'Synonymous', basefmt = "grey")
    ax2.stem(np.array(list(y_nontrunc.keys())).astype(int), np.array(list(y_nontrunc.values())),
             linefmt = "blue", markerfmt = "bo", label = 'Non-truncating', basefmt = "grey")
    ax2.stem(np.array(list(y_trunc.keys())).astype(int), np.array(list(y_trunc.values())),
             linefmt = "green", markerfmt = "go", label = 'Truncating', basefmt = "grey")
    ax2.set_ylabel("Mutation count", fontsize = 13)
    
    legend_cols_dict = [Line2D([0],[0], color='w', markerfacecolor = "red", marker = 'o', 
                               markersize = 8, label = 'Synonymous variant'),
                       Line2D([0],[0], color='w', markerfacecolor = "blue", marker = 'o', 
                               markersize = 8, label = 'Non-truncating variant'),
                       Line2D([0],[0], color='w', markerfacecolor = "green", marker = 'o', 
                               markersize = 8, label = 'Truncating variant')]
                        
    legend = ax2.legend(handles = legend_cols_dict, loc = 'upper left',
                         fontsize = 12, framealpha = 0.3)
    ax2.add_artist(legend)
    
    # Generate stability change plot
    ax1.plot(np.array(list(y_stabch.keys())).astype(int), np.array(list(y_stabch.values())),
             color = "#fd7666")
    ax1.set_ylabel("Stability\nchange", fontsize = 11)
    ax1.set_title(protein_symb+" ("+dataset+")", size = 18)
    
    # Generate sequence plot
    ax3.set_ylim(0, 1)
    rect = patches.Rectangle(xy = (degron_start-1, 0), width = (degron_end-1)-(degron_start-1),
                             height = 10, color = "red", alpha = 0.5, zorder = 2)
    ax3.add_patch(rect)
    ax3.set_xlabel("Protein position", fontsize = 15)
    ax3.set_yticks([])

    # Save
    plt.savefig(plot_file, dpi = 800, transparent = True, bbox_inches = "tight")

    return None

def stabch_E3ligase_boxplot(conditions, subsets_dict, dataset, stabch, ylim1, ylim2, palette, plot_file,
                            main_plot = "boxplot", do_stats = True, stripplot = True, pad_stats = 130, 
                            fig_width = 10, fig_height = 12):
    """
    Multiple boxplot to compare stability change levels between conditions for several motifs.

    Parameters
    ----------
    conditions: list
            Conditions to compare.
    subsets_dict: dict
            Dictionary of the form {condition: dataframe}. Contains all the samples from a mutations dataframe
            corresponding to each condition already filtered. 
    dataset: str
            Dataset name
    stabch: pandas dataframe
            Dataframe containing mutations and stability change levels of several proteins.
    y_lim1: int
            Y-axis top limit
    y_lim2: int
            Y-axis bottom limit
    palette: dict
            Dictionary of the form {condition: color} to map each condition to a specific color.
    plot_file: str
            Path to the folder where to store the resulting plot
    main_plot: str (default: "boxplot")
            If boxplot: the main plot is a seaborn boxplot
            If violinplot: the main plot is a violinplot
    do_stats: boolean (default: True)
            If True, a Mann Whitney test will be computed between pairs
            If False, no statistical test is computed
    stripplot: boolean (default: True)
            If True, a seaborn stripplot will be also plotted
    pad_stats: int (default: 130)
            The pad between the statistical annotation and the plot's title. Only considered
            when do_stats = True
    fig_width: int (default: 10)
            Figure width
    fig_height: int (default: 12)
            Figure height
    
    Returns
    -------
    ax: matplotlib plot
            Stability change boxplot comparing selected conditions for several motifs
    """
    
    # To drop duplicates
    cols_to_remove = ["E3", "degron_start", "degron_end", "Loc_mut_degron"]
    cols = stabch.columns.tolist()
    cols_for_drop = [col for col in cols if col not in cols_to_remove]
    
    # Select conditions and generate subsets
    ## axis and palette
    conditions_palette = {}
    for cond in conditions:
        conditions_palette[cond] = palette[cond]

    ## data subsets: concatenate in unique dataframe   
    subset = concat_subsets(subsets_dict, conditions, cols_for_drop)
    
    # Plot settings
    sns.set_style(style = 'whitegrid')
    plt.figure(figsize = (fig_width, fig_height)) 
    plt.ylim([ylim1, ylim2])
    pad = 10
    title_size = 20
    if do_stats:
        pad = pad_stats
    plt.title(dataset+' degrons', pad = pad, size = title_size)
    
    # Filter: motifs which have before/inside condition
    motifs = subset.E3.unique()
    motifs_f = []
    for motif in motifs:
        motif_conditions = subset.loc[subset.E3 == motif].Condition.unique()
        if (("trunc_inbf_muts" in motif_conditions) and ("trunc_aft_muts" in motif_conditions))\
             or (("nontrunc_in_muts" in motif_conditions) and ("nontrunc_out_muts" in motif_conditions)):
            motifs_f.append(motif)
    subset = subset.loc[subset.E3.isin(motifs_f)]

    # Stripplot
    if stripplot:
        ax = sns.stripplot(x = "E3", y = "Stability_Change", data = subset,
                        alpha = 0.3, size = 4, palette = conditions_palette,
                        jitter = True, dodge = True, hue = "Condition", hue_order = conditions,
                        order = motifs_f)
    # Main plot: boxplot or violinplot
    if main_plot == "violinplot":
        ax = sns.violinplot(x = "E3", y = "Stability_Change", data = subset,
                        width = 0.7, saturation = 0.7, palette = conditions_palette, fliersize = 0., 
                        linewidth = 1, hue = "Condition", hue_order = conditions, order = motifs_f)
    elif main_plot == "boxplot":
        ax = sns.boxplot(x = "E3", y = "Stability_Change", data = subset,
                        width = 0.7, saturation = 0.7, palette = conditions_palette, fliersize = 0., 
                        linewidth = 1, hue = "Condition", hue_order = conditions, order = motifs_f)
        
    plt.xticks(fontsize = 13, rotation = 45, ha = "right")
    plt.ylabel("Stability change", fontsize = 18)
    plt.xlabel("Degron motif", fontsize = 18)

    if "trunc_inbf_muts" in conditions:
        label_in = 'Before/inside degron'
        label_out = 'After degron'
    elif "nontrunc_in_muts" in conditions:
        label_in = 'Inside degron'
        label_out = 'Outside degron'

    # Add customized legend for hue
    legend_cols_dict = [Line2D([0],[0], color ='w', markerfacecolor = "#33a02c", marker = 'o', 
                               markersize = 8, label = 'WT'),
                       Line2D([0],[0], color ='w', markerfacecolor = "#e31a1c", marker = 'o', 
                               markersize = 8, label = label_in),
                       Line2D([0],[0], color ='w', markerfacecolor = "#6a3d9a", marker = 'o', 
                               markersize = 8, label = label_out)]
                        
    legend = ax.legend(handles = legend_cols_dict, loc = 'upper left',
                         fontsize = 12, framealpha = 0.3)
    ax.add_artist(legend)
        
    # Perform Mann Whitney test and annotate p values
    if do_stats:
        
        ## create pairs (with hue)
        pairs = []
        for motif in motifs_f:
            i = 0
            for cond1 in conditions:
                i += 1
                for cond2 in conditions[i:]:
                    pairs.append([(motif, cond1), (motif, cond2)])
                
        ## annotate stats            
        annotator = Annotator(ax, pairs, data = subset, x = "E3", y = 'Stability_Change', hue = "Condition",
        hue_order = conditions, order = motifs_f)
        annotator.configure(text_format = "star", loc = 'outside', test = 'Mann-Whitney').apply_and_annotate()
    
    # Save
    plt.savefig(plot_file, dpi = 800, transparent = True, bbox_inches = "tight")

    return ax