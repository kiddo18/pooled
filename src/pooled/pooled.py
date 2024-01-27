
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties, findfont
import gpplot
from adjustText import adjust_text

# Find the path to the Helvetica font
font_path = findfont(FontProperties(family='Helvetica'))
# Set titles and labels with Helvetica font
font_properties = FontProperties(fname=font_path, size=12)

def lognorm(reads):
    """
    Standardize read counts by calculating reads per million,
    adding a pseudo-count of one, and taking the log2

    :param reads: numpy or pandas array
    :returns: numpy or pandas array
    """
    reads_per_million = (reads/reads.sum())*(10**6)
    lognormed_reads = np.log2(reads_per_million + 1)
    return lognormed_reads

def filter_pDNA(filtered_lognorms, pdna_cols = ['pDNA'], z_low = -3):
	'''
	pdna_cols = ['pDNA'] # Specify the name of the pDNA column in df
	z_low = -3  # minimum z-score threshold

	'''

	# Z-score the pDNA columns
	z_scored_cols = []
	for pdna in pdna_cols:
	    z_col = pdna + '_z' # new column's label
	    # Z-score of sgRNA_i = (lognorm of sgRNA_i - mean of all sgRNA lognorms) / std dev of all sgRNA lognorms
	    filtered_lognorms[z_col] = (filtered_lognorms[pdna] - filtered_lognorms[pdna].mean())/filtered_lognorms[pdna].std()
	    z_scored_cols.append(z_col)
	    
	# Filter by z-score
	filtered_lognorms = filtered_lognorms[filtered_lognorms[z_scored_cols].min(axis = 1) > z_low]

	# Drop z-scored pDNA columns
	filtered_lognorms = filtered_lognorms.drop(z_scored_cols, axis=1)

	# Output
	print('Filtered ' + str(lognorm_df.shape[0] - filtered_lognorms.shape[0]) + ' sgRNAs with low pDNA')

	return filtered_lognorms

def calculate_LFC(lfc_df, filtered_lognorms, ref_map):

	'''
	Input: 
		- lfc_df
		- filtered_lognorms df
		- ref_map dict

	Output: LFC df
	'''

	# Iterate over all contrasts
	for target_col, ref_col in ref_map.items():
	    # Calculate log2 fold change
	    lfc_df[target_col+"_v_"+ref_col] = filtered_lognorms[target_col] - filtered_lognorms[ref_col]

	# Remove reference columns
	for ref_col in set(ref_map.values()):
	    if ref_col not in ref_map.keys(): # not a target condition as well
	        lfc_df = lfc_df.drop(ref_col, axis=1)

	# Remove all columns except for the contrasts "_v_" columns
	for col in lfc_df.columns.to_list():
	    if ("_v_" not in col) and col!= "Sequence":
	        lfc_df = lfc_df.drop(col, axis=1)

	return lfc_df

def pivot(lfc_df, condition_map, rep_map):

	# Transform df from wide to long format
	long_lfcs = (lfc_df.melt(id_vars='Sequence',
	                         var_name='condition', value_name='lfc'))

	# Rename the "condition" column to the target sample name
	long_lfcs['Rep'] = long_lfcs['condition'].copy()
	long_lfcs['Rep'] = long_lfcs['Rep'].map(rep_map)
	long_lfcs['condition'] = long_lfcs['condition'].map(condition_map)

	return long_lfcs

def annotate(long_lfcs, guide_annotations):

	annotated_sgrna_lfcs = (long_lfcs.merge(guide_annotations,
                                       how='inner',
                                       on='Sequence'))
	return annotated_sgrna_lfcs

def plot_histogram(dataframe, col="foo", hue="Rep", title='foo', xlabel='foo', ylabel='foo', figsize=(6,4), bins=100, font_properties=font_properties, title_label_fontsize=14, xlabel_fontsize=12, ylabel_fontsize=12):

	# Set up figure and axis
	fig, ax = plt.subplots(figsize=figsize)

	# Create a histogram plot with adjusted parameters
	sns.histplot(data=dataframe, x=col, kde=True, hue=hue, bins=bins, palette='muted')

	# Title and axis labels
	ax.set_title(title, fontproperties=font_properties, fontsize=title_label_fontsize)
	ax.set_xlabel(xlabel, fontproperties=font_properties, fontsize=xlabel_fontsize)
	ax.set_ylabel(ylabel, fontproperties=font_properties, fontsize=ylabel_fontsize)

	# Remove spines
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)

	# Customize ticks
	ax.tick_params(axis='both', which='both', length=0)
	ax.xaxis.set_tick_params(pad=8)
	ax.yaxis.set_tick_params(pad=8)

	# Remove gridlines
	ax.grid(False)

	# Adjust layout
	plt.tight_layout()

	# Show the plot
	plt.show()

	return None

def calculate_control_sgRNA_stats(annotated_sgrna_lfcs, negative_ctrl_genes):

	# First, filter for rows that are negative control genes
	control_sgrna_lfcs = annotated_sgrna_lfcs[annotated_sgrna_lfcs['Gene'].isin(negative_ctrl_genes)].reset_index(drop=True)

	# Then calculate summary statistics
	control_sgrna_stats = (control_sgrna_lfcs.groupby(['condition','Rep'])
                       .agg(neg_ctl_mean = ('lfc', 'mean'),
                            neg_ctl_sd = ('lfc', 'std'))
                       .reset_index())

	return control_sgrna_lfcs, control_sgrna_stats

def calculate_sgRNA_LFCs(annotated_sgrna_lfcs, control_sgrna_stats, guide_annotations, return_tmp=False):

	# Step 1: Normalize the filtered sgRNA-level LFCs by the mean of LFC of the negative control sgRNAs 
	### (stored in `control_sgrna_stats[neg_ctl_mean]`) within each replicate: `norm_lfc`
	NTC_norm_sgrna_lfcs = annotated_sgrna_lfcs.copy()
	NTC_norm_sgrna_lfcs = (NTC_norm_sgrna_lfcs.merge(control_sgrna_stats,
	                                               how='inner', on=['condition','Rep']))
	NTC_norm_sgrna_lfcs['norm_lfc'] = ((NTC_norm_sgrna_lfcs['lfc'] - NTC_norm_sgrna_lfcs['neg_ctl_mean']))
	NTC_norm_sgrna_lfcs = NTC_norm_sgrna_lfcs.drop(['neg_ctl_sd', 'neg_ctl_mean'], axis=1)

	# Step 2: Rank the filtered, normalized sgRNA-level LFCs within each replicate: `sgRNA_Rank`
	### Add a new column 'Rank' based on 'norm_lfc' within each group (condition, Rep)
	NTC_norm_sgrna_lfcs['sgRNA_Rank'] = NTC_norm_sgrna_lfcs.groupby(['condition','Rep'])['norm_lfc'].rank(ascending=False, method='min').astype(int)


	# Step 3: Combine the replicates:
    #### 3.1) `norm_avg_lfc` = mean across all replicates' `norm_lfc`
    #### 3.2) `sgRNA_sum_Rank` = sum across all replicates' `sgRNA_Rank`
	NTC_norm_avg_sgrna_lfcs = NTC_norm_sgrna_lfcs.copy()

	NTC_norm_avg_sgrna_lfcs = NTC_norm_avg_sgrna_lfcs.groupby(['Sequence','condition']).agg(
	    norm_avg_lfc = ('norm_lfc', 'mean'),
	    sgRNA_sum_Rank = ('sgRNA_Rank', 'sum'),
	).reset_index()

	# Before returning, merge with annotations
	NTC_norm_avg_sgrna_lfcs = (NTC_norm_avg_sgrna_lfcs.merge(guide_annotations,
                                       how='inner',
                                       on='Sequence'))

	# Return dataframe
	if return_tmp:
		return NTC_norm_avg_sgrna_lfcs, NTC_norm_sgrna_lfcs, 
	else:
		return NTC_norm_avg_sgrna_lfcs

def calculate_gene_LFCs(NTC_norm_avg_sgrna_lfcs):

	NTC_norm_gene_lfcs = NTC_norm_avg_sgrna_lfcs.copy()
	NTC_norm_gene_lfcs = (NTC_norm_gene_lfcs.groupby(['Gene', 'condition']).agg(gene_norm_avg_lfc = ('norm_avg_lfc', 'mean'),
                                                                                             gene_Rank = ('sgRNA_sum_Rank', 'median'),
                                                                                            n_sgrnas = ('norm_avg_lfc', 'count')).reset_index())

	return NTC_norm_gene_lfcs

def simulate_p_value(NTC_norm_avg_sgrna_lfcs, n_reps = 2, n_iterations=500, seed = 318):
    
    # Set a seed for reproducibility
    np.random.seed(seed)

    # Prep variables
    null_rank_df0 = NTC_norm_avg_sgrna_lfcs.copy()[['Gene']]
    total_sgRNAs = int(NTC_norm_avg_sgrna_lfcs.groupby('condition').count()['Sequence'])

    # Intialize empty null rank matrix
    null_rank_matrix = np.zeros((total_sgRNAs, n_iterations))

    # Iterate
    for i in range(n_iterations): # For each iteration (column)
        random_ranks = []
        for rep in range(n_reps): # Sum the ranks across replicates
            random_ranks.append(np.random.choice(np.arange(1, total_sgRNAs + 1), total_sgRNAs, replace=False))
        null_rank_matrix[:, i] = np.sum(np.array(random_ranks), axis = 0)

    # Wrangle
    null_rank_df1 = pd.DataFrame(null_rank_matrix)
    null_rank_df = pd.concat([null_rank_df0, null_rank_df1], axis = 1)
    null_rank_df2 = null_rank_df.groupby('Gene').agg('median').reset_index().drop(["Gene"], axis = 1)
    null_rank_vector = null_rank_df2.values.flatten()
    len_null_rank_vector = len(null_rank_vector)


    return null_rank_df2, null_rank_vector, len_null_rank_vector

def calculate_p_value(NTC_norm_gene_lfcs, null_rank_vector, len_null_rank_vector, significant_threshold=0.05):
    

    # Initialize list of p-values for each gene
    p_values = [] 

    for i in range(NTC_norm_gene_lfcs.shape[0]):  # Loop through each gene
        count_extreme = np.sum(null_rank_vector <= NTC_norm_gene_lfcs["gene_Rank"][i])
        p_value = np.divide(count_extreme, len_null_rank_vector)
        #print(NTC_norm_gene_lfcs["Gene"][i], "-", count_extreme, "|", p_value)
    
        # Adjust for 2-tailed test
        if p_value == 0. or p_value == 1.:
            p_value = np.divide(1, len_null_rank_vector) * 2
        elif p_value < 0.5:
            p_value = p_value * 2 
        elif p_value >= 0.5:
            p_value = (1 - p_value) * 2
        
        # Add to list of p-values
        p_values.append(p_value)


    # Merge list of p-values with NTC_norm_gene_lfcs dataframe
    # Calculate the p-value
    NTC_norm_gene_lfcs_pval = NTC_norm_gene_lfcs.copy()
    NTC_norm_gene_lfcs_pval['pval'] = p_values
    
    # Calculate -log10 pval
    NTC_norm_gene_lfcs_pval['nlog10pval'] = -np.log10(NTC_norm_gene_lfcs_pval['pval'])
    
    # Replace Inf if there are any log(0) execution errors
    max_value = np.nanmax(np.where(NTC_norm_gene_lfcs_pval['nlog10pval'] == np.inf, np.nan, NTC_norm_gene_lfcs_pval['nlog10pval']))
    NTC_norm_gene_lfcs_pval['nlog10pval'][NTC_norm_gene_lfcs_pval['nlog10pval'] == np.inf] = max_value


    # Adjust for multiple-hypothesis testing by BH procedure
    ### Use statsmodels function multipletests to adjust FDR
    FDR = multipletests(NTC_norm_gene_lfcs_pval['pval'], alpha=significant_threshold, method='fdr_bh')

    NTC_norm_gene_lfcs_adjpval = NTC_norm_gene_lfcs_pval.copy()

    NTC_norm_gene_lfcs_adjpval["significant"] = FDR[0] # by signficant_threshold = 0.05
    NTC_norm_gene_lfcs_adjpval["FDR"] = FDR[1]
    NTC_norm_gene_lfcs_adjpval["nlog10FDR"] = -np.log10(NTC_norm_gene_lfcs_adjpval['FDR'])

    # Replace Inf if there are any log(0) execution errors
    max_value = np.nanmax(np.where(NTC_norm_gene_lfcs_adjpval['nlog10FDR'] == np.inf, np.nan, NTC_norm_gene_lfcs_adjpval['nlog10FDR']))
    NTC_norm_gene_lfcs_adjpval['nlog10FDR'][NTC_norm_gene_lfcs_adjpval['nlog10FDR'] == np.inf] = max_value

    # Print

    print("Number of significant hits at FDR =", significant_threshold, "is", NTC_norm_gene_lfcs_adjpval["significant"].sum())

    return NTC_norm_gene_lfcs_adjpval

def plot_scatter(df, x, y, labels, significant_genes_to_label, title='foo', xlabel='foo', ylabel='foo', figsize=(6,4), font_properties=font_properties, title_label_fontsize=14, xlabel_fontsize=12, ylabel_fontsize=12):

    # Set up figure and axis
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create a scatter plot with adjusted parameters
    scatter = ax.scatter(x, y, color='gray', s=50, alpha=0.3, edgecolors='none')
    
    # Title and axis labels
    ax.set_title(title, fontproperties=font_properties, fontsize=title_label_fontsize)
    ax.set_xlabel(xlabel, fontproperties=font_properties, fontsize=xlabel_fontsize)
    ax.set_ylabel(ylabel, fontproperties=font_properties, fontsize=ylabel_fontsize)
    
    # Remove spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Customize ticks
    ax.tick_params(axis='both', which='both', length=0)
    ax.xaxis.set_tick_params(pad=8)
    ax.yaxis.set_tick_params(pad=8)
    
    # Customize grid
    #ax.grid(color='lightgray', linestyle='-', linewidth=0.5, alpha=0.5)
    
    # Remove gridlines
    ax.grid(False)
    
    # Set y-axis limits
    #ax.set_ylim(0, y.max() + 0.3 )  # Adjust the limits according to your specific needs
    
    # Label points with their labels
    offset = 0.
    annotations = []
    for i, label in enumerate(labels):
        if label in significant_genes_to_label:
            annotation = ax.text(x[i] + offset, y[i] + offset, label, ha='center', va='bottom', fontproperties=font_properties, fontsize=10, color='black')
            annotations.append(annotation)
            #ax.annotate(label, (x[i], y[i]), textcoords="offset points", xytext=(0,10), ha='center', fontproperties=font_properties, fontsize=6, color='red', bbox=dict(boxstyle='round,pad=0.3', edgecolor='red', facecolor='white'))
            ax.scatter(x[i], y[i], color='red', edgecolor='black')
            
    # Use adjust_text to automatically adjust the position of text labels to avoid overlap
    adjust_text(annotations, arrowprops=dict(arrowstyle='-', color='red'), force_text=(0.5, 0.8), expand=(1.6, 1.6))
    
    # Adjust layout
    plt.tight_layout()
    
    # Show the plot
    plt.show()

    return None

def foo():
	print("You have imported pooled!")
	return None

