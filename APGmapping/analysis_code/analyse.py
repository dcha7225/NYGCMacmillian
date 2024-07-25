import pysam
import pyBigWig
from openpyxl import Workbook
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import scipy.stats as stats
import numpy as np
import ast
import concurrent.futures
from tqdm import tqdm
import copy
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio
import plotly.express as px
from matplotlib.gridspec import GridSpec
from scipy.stats import kruskal
from scipy.stats import gaussian_kde



def chi_square_test(observed_counts):
    """
    Perform a Chi-square test for independence and return the p-value.

    Parameters:
    observed_counts (list of lists): A 2x2 contingency table of observed counts.

    Returns:
    float: p-value from the Chi-square test.
    """
    # Convert the input to a numpy array
    observed = np.array(observed_counts)

    # Perform the Chi-square test

    res = stats.chi2_contingency(observed)
    chi2, p, dof, expected = res
    return p

##### example chi square test #######

'''
observed_counts = [[1961, 67], [107721, 14491]]
p_value = chi_square_test(observed_counts)
print(f"The p-value from the Chi-square test is: {p_value}")

observed_counts2 = [[101517, 12371], [8165, 2187]]
p_value2 = chi_square_test(observed_counts2)
print(f"The p-value from the Chi-square test is: {p_value2}")
'''


def get_histo(values1, values2, filename, title="Value Distribution", step=0.1, label1="Dataset 1", label2="Dataset 2", xlabel="Intervals"):
    """
    Generates overlaid histogram plot with density lines

    Params:
        values1 (float list):  list of values to plot
        values2 (float list):  list of values to plot
        filename (str): name of output file
        title (str): overall plot title 
        step (float): interval steps
        label1 (str): label for first dataset
        label2 (str): label for second dataset
        xlabel (str): label for x axis
    
    Returns:
        a histogram saved as png file
    """



    start = 0
    end = 1

    bins = [x for x in range(int(start * 10), int(end * 10 + 1), int(step * 10))]
    bins = [x / 10 for x in bins]

    plt.clf()
    
    # Plot the first histogram with density=True
    plt.hist(values1, bins=bins, edgecolor='black', alpha=0.5, label=label1, density=True, color = "blue")
    
    # Plot the second histogram with density=True
    plt.hist(values2, bins=bins, edgecolor='black', alpha=0.5, label=label2, density=True, color = "red")
    
    # Calculate and plot KDE for the first dataset
    kde1 = gaussian_kde(values1)
    x1 = np.linspace(start, end, 1000)
    line1, = plt.plot(x1, kde1(x1), color='blue', linestyle='--', label=f'{label1} KDE')
    line1.set_label('_nolegend_')
    # Calculate and plot KDE for the second dataset
    kde2 = gaussian_kde(values2)
    x2 = np.linspace(start, end, 1000)
    line2, = plt.plot(x2, kde2(x2), color='red', linestyle='--', label=f'{label2} KDE')
    line2.set_label('_nolegend_')

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel('Density')
    plt.legend()  
    plt.savefig(filename)


# these following functions generate dictionaries to translate sequence aliases

# Reference sequence -> chromosome (i.e NC_060925.1 -> chr1)
def get_ref_chrom_map(ref_seq_path):
    ref_chrom_map = {}

    with open(ref_seq_path, 'r') as file:
        next(file)

        for line in file:
            result = line.strip().split()
            ref_chrom_map[result[8]] = result[11]
    return ref_chrom_map


# Reference sequence -> geneBank id (i.e NC_060925.1 -> CP068277.2)
def get_ref_genbank_map(ref_seq_path):
    ref_genbank_map = {}

    with open(ref_seq_path, 'r') as file:
        next(file)

        for line in file:
            result = line.strip().split()
            ref_genbank_map[result[8]] = result[6]
    return ref_genbank_map


# Reference sequence -> chromosome (i.e NC_000001.11 -> chr1)
def get_hg38_chrom_names(hg38_alias_path):
    ucsc_chrom_map = {}
    with open(hg38_alias_path, 'r') as file:
        next(file)

        for line in file:
            result = line.strip().split()
            if len(result[-1]) == 1:
                code = result[-2][:-1]
            else:
                code = result[-1]

            ucsc_chrom_map[code] = result[0]
    return ucsc_chrom_map


# accession name -> sequence_name (i.e PDBU01000001.1 -> CAAPA_TwoEndPlaced_1)
def get_contig_name(contig_name_path):
    acc_name_path = {}
    with open(contig_name_path, 'r') as file:
        next(file)

        for line in file:
            result = line.strip().split('\t')
            if result[1] in acc_name_path:
                print("alr found")
            acc_name_path[result[1]] = result[0]

    return acc_name_path

# accession name -> contig length (i.e PDBU01000001.1 -> 2495)
def get_contig_lengths(contig_lengths_path):
    contig_lengths = {}
    with open(contig_lengths_path, 'r') as file:
        for line in file:
            result = line.strip().split(',')
            contig_lengths[result[0]] = int(result[1])
    return contig_lengths


def check_unique(bedFile_path, chrom, a_start, a_end, complete = False):
    """
    checks if query location is T2T unique

    Params:
        bedFile_path (float list):  T2T unique annotation bedfile
        chrom (str): name of query chromosome (i.e "chr1")
        a_start (int): start position of query region
        a_end (int): end position of query region
        complete (bool): flag for complete overlap 
    
    Returns:
        true iff region is unique
    """
    bb = pyBigWig.open(bedFile_path)

    entries = bb.entries(chrom, 0, bb.chroms(chrom), withString=False) #this returns list of all entries  # noqa: E501
    for r_start, r_end in entries:
        if complete:
            if a_start >= r_start and a_end <= r_end:  # completely within unique region  # noqa: E501
                return True
        else:
            if a_start <= r_end and a_end >= r_start:  # partially within unique region  # noqa: E501
                return True

    return False

def get_CPGdata(bedFile_path, chrom, a_start, a_end):
    """
    determines number of CpG islands in query region

    Params:
        bedFile_path (float list):  T2T CpG island annotation bedfile
        chrom (str): name of query chromosome (i.e "chr1")
        a_start (int): start position of query region
        a_end (int): end position of query region
    
    Returns:
        number of CpG islands in region
    """

    bb = pyBigWig.open(bedFile_path)
    entries = bb.entries(chrom, 0, bb.chroms(chrom)) #this returns list of all entries
    island_count = 0
    for r_start, r_end, string in entries:
        if a_start <= r_end and a_end >= r_start:
            name,length,cpgNum,gcNum,perCpg,perGc,obsExp = string.strip().split('\t') # CpG: 36\t390\t36\t258\t18.5\t66.2\t0.85  # noqa: E231, E501
            island_count += 1

    return island_count


def filter(df, coverage_thresh, identity_thresh, cumul=False):
    """
    filters contigs that pass threshold, pass threshold & are in unique regions, pass threshold & contain CpG island
    (useful for creating contigency tables)
    Params:
        df (panda.dataframe): reference specific contig alignment results (output from create_contig_df)
        coverage_thresh (float): coverage threshold
        identity_thresh (float): identity threshold
        cumul (bool): flag for multiple location analysis
    
    Returns:
        three sets according to description above
    """

    passed = set()
    unique_regions = set()
    contains_cpg_island = set()

    for _, row in df.iterrows():

        if cumul:

            if row['cumulative_coverage'] >= coverage_thresh and row['cumulative_identity'] >= identity_thresh:
                passed.add(row['accession'])
                if any(x > 0 for x in row['cpg_count']):
                    contains_cpg_island.add(row['accession'])
                if True in row['unique_region']:
                    unique_regions.add(row['accession'])

        else:
            for i in range(len(row['coverage'])):
                if row['coverage'][i] >= coverage_thresh and row['identity'][i] >= identity_thresh:
                    passed.add(row['accession'])
                    if row['unique_region'][i] == True:
                        unique_regions.add(row['accession'])
                    if row['cpg_count'][i] > 0:
                        contains_cpg_island.add(row['accession'])
                    break

    return passed, unique_regions, contains_cpg_island



def props_all_alignments(df):
    """
    Counts alignments that are in unique regions, contain CpG islands, and both are in unique regions and contain CpG islands

    Params:
        df (panda.dataframe): reference specific contig alignment results (output from create_contig_df)
    
    Returns:
        three counts according to description above
    """

    unique_only = 0
    unique_CpG = 0
    CpG_only = 0
    total = 0
    for _, row in df.iterrows():
        for i in range(len(row['coverage'])):
            total +=1
            unique = False
            CpG = False
            if row['unique_region'][i] == True:
                unique_only += 1
                unique = True
            if row['cpg_count'][i] > 0:
                CpG_only += 1
                CpG = True
            if CpG and unique:
                unique_CpG +=1

    return unique_only, CpG_only, unique_CpG, total



def props_all_contigs(df):
    """
    Counts contigs that are in unique regions, contain CpG islands, and both are in unique regions and contain CpG islands, and # of total contigs
    **(used to create contingency tables)

    Params:
        df (panda.dataframe): reference specific contig alignment results (output from create_contig_df)
    
    Returns:
        4 counts according to description above
    """

    unique_only = 0
    unique_CpG = 0
    CpG_only = 0
    total = 0
    for _, row in df.iterrows():
        total +=1
        unique = False
        CpG = False
        if True in row['unique_region']:
            unique_only += 1
            unique = True
        if any(x > 0 for x in row['cpg_count']):
            CpG_only += 1
            CpG = True
        if CpG and unique:
            unique_CpG +=1

    return unique_only, CpG_only, unique_CpG, total


def cov_id_distribution(df):
    """
    Prepares data for scatterplots.
    Creates lists for coverage and identity for alignments with and w/o CpG islands
    **(used as input for histograms)
    Params:
        df (panda.dataframe): reference specific contig alignment results (output from create_contig_df)
    
    Returns:
        4 lists:
            cpg_identity = list of identity scores for CpG island containing alignments
            cpg_coverage = list of coverage scores for CpG island containing alignments
            noCPG_identity = list of identity scores for alignments w/o CpG islands
            noCpG_coverage = list of coverage scores for alignments w/o CpG islands
    """

    noCpG_coverage = []
    noCPG_identity = []

    cpg_coverage = []
    cpg_identity = []

 
    for _, row in df.iterrows():
        for i in range(len(row['coverage'])):
            if row['cpg_count'][i] > 0:
                cpg_coverage.append(row['coverage'][i])
                cpg_identity.append(row['identity'][i])
            else:
                noCpG_coverage.append(row['coverage'][i])
                noCPG_identity.append(row['identity'][i])

    return  cpg_identity, cpg_coverage, noCpG_coverage, noCPG_identity


def create_contig_df(sam_file_path, ref_seq_path, hg38_alias_path, contig_name_path, contig_lengths_path, CPG_bedFile_path = None, unique_bedFile_path= None, ref=None):
    """
    prepares df containing reference specific contig alignment results.

    Params:
        sam_file_path (str): sam file containing alignment results
        ref_seq_path (str): path for tsv file containing T2T alias names
        hg38_alias_path (str): path for tsv file containing GRCh38 alias names
        contig_name_path (str): path for tsv file containing APG Contig alias names
        contig_lengths_path (str): path for tsv file containing APG Contig lengths
        CPG_bedFile_path (str): path for T2T CpG annotations
        unique_bedFile_path (str): path for T2T Unique region annotations
        ref (str): name of ref sequence (only relevant if ref = "T2T" to allow for annotations)
    Returns:
        a df of format:{    'contig': contig name,
                            'accession': contig accession name,
                            'query_positions': a list of query positions as tuples (query_positions[i] = query position of ith alignment),
                            'contig_length': contig length,
                            'mapped_chromosome': a list of chromosomes (mapped_chromosome[i] = mapped chromosome of ith alignment),
                            'ref_positions': a list of reference positions as tuples (ref_positions[i] = reference position of ith alignment),
                            'unique_region': a list of booleans indicating if alignment is in unique region (unique_region[i] = true iff ith alignment is in unique region),
                            'coverage': a list of coverage scores (coverage[i] = coverage of ith alignment),
                            'identity': a list of identity scores (identity[i] = identity of ith alignment),
                            'cpg_count': a list of CpG island counts (cpg_count[i] = cpg island count of ith alignment),
                            'cumulative_coverage': cumulative coverage across all alignments,
                            'cumulative_identity': cumulative identity accross all alignments
                        }
        note: query position = position relative to query sequence. (i.e the region along the contig that is being mapped)
              reference position = position relative to reference sequence (i.e the region along the reference genome that the contig is mapped to)
            
    """

    ucsc_chrom_map = get_hg38_chrom_names(hg38_alias_path)
    ref_chrom_map = get_ref_chrom_map(ref_seq_path)
    ref_genbank_map = get_ref_genbank_map(ref_seq_path)
    acc_contig_map = get_contig_name(contig_name_path)
    contig_lengths = get_contig_lengths(contig_lengths_path)

    contig_stats = {}

    samfile = pysam.AlignmentFile(sam_file_path, "r")

    for alignment in samfile:
        qname = alignment.query_name
        if qname not in contig_stats:
            contig_stats[qname] = {'contig': acc_contig_map[qname],
                                    'accession': qname,
                                    'query_positions': [],
                                    'contig_length': contig_lengths[qname],
                                    'mapped_chromosome': [],
                                    'ref_positions': [],
                                    'unique_region': [],
                                    'coverage': [],
                                    'identity': [],
                                    'cpg_count': [],
                                    'cumulative_coverage': 0,
                                    'cumulative_identity': 0
                                    }

        if alignment.is_unmapped:
            continue

        #location info
        rname = alignment.reference_name
        chrom = rname

        ref_start = alignment.reference_start
        ref_end = alignment.reference_end

        query_start = alignment.query_alignment_start
        query_end = alignment.query_alignment_end

        if ref == "GRCh38":
            if rname in ucsc_chrom_map:
                chrom = ucsc_chrom_map[rname]
                end = chrom.split('_')[-1]
                if end == "alt" or end == "fix": #skipping alternate loci and patches
                    continue


        #identity per align
        aligned_length =  query_end - query_start
        distance = alignment.get_tag('NM')
        align_idenity = (aligned_length - distance)/aligned_length

        contig_stats[qname]['query_positions'].append((query_start, query_end))
        contig_stats[qname]['ref_positions'].append((ref_start, ref_end))
        contig_stats[qname]['mapped_chromosome'].append(chrom)
        contig_stats[qname]['identity'].append(align_idenity)

    for qname in contig_stats:
        #total cumulative
        weighted_tot_id = 0
        total_weight = 0

        for i in range(len(contig_stats[qname]['identity'])):
            weight = contig_stats[qname]['query_positions'][i][1] - contig_stats[qname]['query_positions'][i][0]
            total_weight += weight
            weighted_tot_id += contig_stats[qname]['identity'][i]*weight

        tot_aligned_length, _ = calculateTotAlignedLength(contig_stats[qname]['query_positions'])

        if tot_aligned_length>0:
            contig_stats[qname]['cumulative_identity'] = weighted_tot_id/total_weight
            contig_stats[qname]['cumulative_coverage'] = tot_aligned_length/contig_stats[qname]['contig_length']


        #chaining together alignments within 300 bp
        groups = groupIntervals(contig_stats[qname]['mapped_chromosome'], contig_stats[qname]['ref_positions'])
        group_chrom=[]
        all_group_ref_positions=[]
        group_identity=[]
        merged_group_query_positions = []

        for group in groups:
            weighted_group_id = 0
            total_weight = 0
            chrom = group['Chromosome']
            raw_group_query_positions = []

            for i in group['OriginalIndices']:
                query_position = contig_stats[qname]['query_positions'][i]
                raw_group_query_positions.append(query_position)
                weight = query_position[1] - query_position[0]
                total_weight += weight
                weighted_group_id += contig_stats[qname]['identity'][i]*weight

            group_aligned_length, merged = calculateTotAlignedLength(raw_group_query_positions)

            group_identity.append(weighted_group_id/total_weight)
            contig_stats[qname]['coverage'].append(group_aligned_length/contig_stats[qname]['contig_length'])
            merged_group_query_positions.append(merged)

            ###T2T annotations
            if ref=="T2T":
                genbank = ref_genbank_map[group['Chromosome']]
                chrom = ref_chrom_map[group['Chromosome']]

                isUnique=False
                cpgCount = 0
                for position in group['group_ref_positions']:
                    start,end = position
                    if not isUnique and check_unique(unique_bedFile_path, chrom, start, end):
                        isUnique = True

                    cpgCount += get_CPGdata(CPG_bedFile_path, genbank, start, end)

                contig_stats[qname]['unique_region'].append(isUnique)
                contig_stats[qname]['cpg_count'].append(cpgCount)

            group_chrom.append(chrom)
            all_group_ref_positions.append(group['group_ref_positions'])

        contig_stats[qname]['query_positions'] = merged_group_query_positions
        contig_stats[qname]['mapped_chromosome'] = group_chrom
        contig_stats[qname]['ref_positions'] = all_group_ref_positions
        contig_stats[qname]['identity'] = group_identity


    samfile.close()
    return pd.DataFrame(contig_stats.values())


### Helper functions for create_contig_df #########
def calculateTotAlignedLength(intervalsOG):
    """
    calculates total aligned length given intervals that may contain overlaps
    """
    # Sort intervals by start position
    intervals = intervalsOG.copy()

    intervals.sort(key=lambda x: x[0])

    # Initialize the list to hold merged intervals
    merged = []

    for interval in intervals:
        if not merged or merged[-1][1] < interval[0]:
            # If merged is empty or there is no overlap, add the interval to merged
            merged.append(interval)
        else:
            # If there is an overlap, merge the intervals
            merged[-1] = (merged[-1][0], max(merged[-1][1], interval[1]))

    total_align_length = 0
    for interval in merged:
        total_align_length += interval[1]-interval[0]

    return total_align_length, merged

def groupIntervals(chromosomes, positions):
    """
    merges overlapping intervals (handles alignment query regions that overlap)
    """
    starts = [pos[0] for pos in positions]
    ends = [pos[1] for pos in positions]

    df = pd.DataFrame({'Chromosome': chromosomes, 'Start': starts, 'End': ends})
    df['OriginalIndex'] = df.index

    # Sort the DataFrame by Chromosome and Start position
    df = df.sort_values(by=['Chromosome', 'Start']).reset_index(drop=True)

    # Function to determine if two intervals are within 300bp
    def within_300bp(row1, row2):
        return row2['Start'] <= row1['End'] + 300

    # Initialize group identifiers
    current_group_id = 0
    group_ids = []

    # Iterate through the sorted dataframe to assign group IDs
    for i in range(len(df)):
        if i == 0:
            group_ids.append(current_group_id)
        else:
            if df['Chromosome'][i] == df['Chromosome'][i - 1] and within_300bp(df.iloc[i - 1], df.iloc[i]):
                group_ids.append(current_group_id)
            else:
                current_group_id += 1
                group_ids.append(current_group_id)

    # Add group IDs to the dataframe
    df['Group'] = group_ids

    # Group the intervals within 500bp and keep track of the original indices
    grouped = df.groupby(['Chromosome', 'Group']).agg({
        'Start': list,
        'End': list,
        'OriginalIndex': list
    }).reset_index()

    # Collect the grouped original indices and group details in a list of dictionaries
    grouped_info = []
    for _, row in grouped.iterrows():
        group_postions = list(zip(row['Start'], row['End']))
        grouped_info.append({
            'Chromosome': row['Chromosome'],
            'group_ref_positions': group_postions,
            'OriginalIndices': row['OriginalIndex']
        })
    return grouped_info


#### CSV functions #######
def create_csv(df,filename="output.csv"):
    df.to_csv(filename, index=False)
    print("dataframe saved to", filename)

def read_csv(csv_path):

    def toLiteral(s):
        try:
            return ast.literal_eval(s)
        except (ValueError, SyntaxError):
            return s
    df = pd.read_csv(csv_path)
    df = df.applymap(toLiteral)

    return df


##### this function is outdated. ##################################
##### Still used for plot.py for generating chromosome annotations ######


def create_alignment_table(sam_file_path, contig_lengths_path, unique_bedFile_path, ref_seq_path, contig_name_path):
    """
    creates alignment data xlsx for T2T reference only. Used for plot.py
    (todo: refactor plot.py to work with alignment csv files)

    Params:
        sam_file_path (str): sam file containing alignment results
        contig_lengths_path (str): path for tsv file containing APG Contig lengths
        unique_bedFile_path (str): path for T2T Unique region annotations
        ref_seq_path (str): path for tsv file containing T2T alias names
        contig_name_path (str): path for tsv file containing APG Contig alias names

    Returns:
       an xlsx file containing alignment results for T2T
    """
    contig_lengths = get_contig_lengths(contig_lengths_path)
    ref_chrom_map = get_ref_chrom_map(ref_seq_path)
    acc_name_path = get_contig_name(contig_name_path)

    wb = Workbook()
    ws = wb.active
    ws.append(['contig_name', 'accession', 'aligned_length', 'mapped_chromosome', 'start_pos', 'end_pos', 'unique_region', 'coverage', 'identity'])

    samfile = pysam.AlignmentFile(sam_file_path, "r")
  
    for alignment in samfile:
        if alignment.is_unmapped:
            continue

        aligned_length = alignment.query_alignment_length
        qname = alignment.query_name
        coverage = aligned_length/contig_lengths[qname]
        distance = alignment.get_tag('NM')
        align_idenity = (aligned_length - distance) / aligned_length

        chrom = ref_chrom_map[alignment.reference_name]
        a_start = alignment.reference_start
        a_end = alignment.reference_end
        unique_region = check_unique(unique_bedFile_path, chrom, a_start, a_end)

        ws.append([acc_name_path[qname], qname, aligned_length, chrom, a_start, a_end, unique_region, coverage, align_idenity])

    samfile.close()


    fname = 'alignments_with_names.xlsx'
    wb.save(fname)

    return fname


def scatterData(csv_path, cumul=False, filename="scatterdata.csv"):
    """
    generates csv file containing data for scatter plots

    Params:
       csv_path (str): path for csv file containing containing reference specific contig alignment results.
                        (csv result of create_contig_df())
        cumul (bool): flag for multiple location analysis
        filename (str): name for output file
    Returns:
        csv file containing scatter data
    """
    df = read_csv(csv_path)

    coverage = np.zeros(124240)
    identity = np.zeros(124240)
    passed = np.zeros(124240)
    if cumul:
        for i, row in df.iterrows():
            coverage[i] = row['cumulative_coverage']
            identity[i] = row['cumulative_identity']
            if row['cumulative_coverage'] >= 0.5 and row['cumulative_identity'] >= 0.8:
                passed[i] = 1
            if row['cumulative_coverage'] >= 0.8 and row['cumulative_identity'] >= 0.9:
                passed[i] = 2

    else:
        for i, row in df.iterrows():
            maxComb = 0
            for j in range(len(row['coverage'])):
                curr = row['coverage'][j] * row['identity'][j]

                if row['coverage'][j] > 1:
                    print(row)
                    break


                if curr > maxComb:
                    maxComb = curr
                    coverage[i] = row['coverage'][j]
                    identity[i] = row['identity'][j]
                if row['coverage'][j] >= 0.5 and row['identity'][j] >= 0.8:
                    passed[i] = 1
                    coverage[i] = row['coverage'][j]
                    identity[i] = row['identity'][j]
                    maxComb = float('inf')

                if row['coverage'][j] >= 0.8 and row['identity'][j] >= 0.9:
                    passed[i] = 2
                    coverage[i] = row['coverage'][j]
                    identity[i] = row['identity'][j]
    data = {
        'coverage': coverage,
        'identity': identity,
        'passed': passed,
    }
    df = pd.DataFrame(data)
    create_csv(df, filename=filename)



def draw_scatterplot(csv_path, title="title", xlabel="X", ylabel="Y", filename='scatterplot.png'):
    """
    generates scatterplots from scatter data

    Params:
        csv_path (str): path for csv file containing scatter data (output of scatterData())
        title (str): overall plot title
        xlabel (str): x axis title
        ylabel (str): y axis title
        filename (str): output file name
    Returns:
        scatterplot with heatmap saved as png
    """

    df = read_csv(csv_path)

    # Create a figure with GridSpec
    fig = plt.figure(figsize=(14, 8))
    gs = GridSpec(1, 2, width_ratios=[3, 1], wspace=0.3)

    ax0 = fig.add_subplot(gs[0])
    ax1 = fig.add_subplot(gs[1])

    # Scatter plot with contour density plot
    scatter = sns.scatterplot(data=df, x='coverage', y='identity', hue='passed', ax=ax0)
    sns.kdeplot(data=df, x='coverage', y='identity', levels=6, cmap='viridis', alpha=0.30, fill=True, ax=ax0, cbar=True)

    ax0.set_xlim(0, 1)
    ax0.set_ylim(0, 1)

    ax0.set_title(title)
    ax0.set_xlabel(xlabel)
    ax0.set_ylabel(ylabel)

    ax0.axhline(y=0.8, color='r', linestyle='--')
    ax0.axvline(x=0.5, color='r', linestyle='--')

    ax0.axhline(y=0.9, color='g', linestyle='--')
    ax0.axvline(x=0.8, color='g', linestyle='--')

    labelMap = {'0.0':'Weak', '1.0':'Reasonably Good', '2.0':'Near-Perfect'}

    handles, labels = scatter.get_legend_handles_labels()
    scatter.legend(handles, [labelMap[l] for l in labels], loc='lower right')

    # Bar plot for percentage of data in each category
    category_counts = df['passed'].value_counts(normalize=True).sort_index()
    categories = category_counts.index.map({0: 'Weak', 1: 'Reasonably Good', 2: 'Near-Perfect'})
    sns.barplot(x=category_counts.values * 100, y=categories, ax=ax1, alpha=0.6)

    ax1.set_xlabel('Percentage of Data (%)')
    ax1.set_xlim(0, 100)

    for i, count in enumerate(category_counts.values * 100):
        ax1.text(count + 1, i, f'{count:.1f}%', va='center')

    plt.savefig(filename)
    plt.clf()

    print("saved to:", filename)


#### outdated bar plotting code, use create_combined_plot_plotly() instead ##########
# this one was bad for many many bars

def plot_grouped_bar_chart(df, x_col, y_col, hue_col, ax, xlabel, ylabel, title='Grouped Bar Chart',palette='hls', bar_width=0.55):
    # Create the bar chart
    barplot = sns.barplot(x=x_col, y=y_col, hue=hue_col, data=df, palette=palette, ax=ax, dodge=True, errorbar=None)


    # Adjust bar width
    for i, bar in enumerate(barplot.patches):
        bar.set_width(bar_width)

    # Add labels to the bars
    for p in barplot.patches:
        barplot.annotate(format(p.get_height(), '.1f'),
                         (p.get_x() + p.get_width() / 2., p.get_height()),
                         ha='center', va='center',
                         xytext=(0, 9))

    # Add title and labels
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    ax.legend(loc='lower left')

def create_combined_plot(df1, df2, x_col, y_col, hue_col, titles, overall_title='Combined Bar Chart',xlabel='Samples', ylabel='Values',filename='combined_bar_chart.png'):

    # Create subplots
    fig, axs = plt.subplots(2, 1, figsize=(55, 30))

    # Plot the first grouped bar chart
    plot_grouped_bar_chart(df1, x_col, y_col, hue_col, axs[0], xlabel, ylabel, title=titles[0])

    # Plot the second grouped bar chart
    plot_grouped_bar_chart(df2, x_col, y_col, hue_col, axs[1], xlabel, ylabel, title=titles[1])

    fig.suptitle(overall_title, fontsize=16)

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])


    plt.savefig(filename)
    plt.clf()



def calcKruskal(df):
    """
    compute Kruskal-Wallis h test for independence:

    Params:
        df (panda.dataframe): a df containing reference comparison data:
                            Columns: 'Ancestry', 'Reference_name', 'Value' 
                            (value represents percentage of contigs passing threshold)
                            **(output of compare_ref())
    Returns:
        H-stat and p value
    """

    African = df[df['Ancestry'] == 'African Ancestry']['Value']
    American = df[df['Ancestry'] == 'American Ancestry']['Value']
    East_asian = df[df['Ancestry'] == 'East Asian Ancestry']['Value']
    European = df[df['Ancestry'] == 'European Ancestry']['Value']
    South_asian = df[df['Ancestry'] == 'South Asian Ancestry']['Value']

    # Perform Kruskal-Wallis H Test
    h_stat, p_val = kruskal(African, American, East_asian, European, South_asian)
    print(f"H-statistic: {h_stat}, P-value: {p_val}")


def create_combined_plot_plotly(reasonable_df, perf_df, overall_title="Title", xlabel="x", ylabel="y", output_file="barPlot.png"):
    """
    Create a Plotly bar chart with four subplots
    
    Parameters:
        reasonable_df (pd.DataFrame): DataFrame with columns 'Ancestry', 'Reference_name', 'Value' for single location analysis.
        perf_df (pd.DataFrame): DataFrame with columns 'Ancestry', 'Reference_name', 'Value' for multiple location analysis.
                                **(value represents percentage of contigs passing threshold)
                                **(output of compare_ref())
        overall_title (str): Overall title of the plot.
        xlabel (str): Label for the x-axis.
        ylabel (str): Label for the y-axis.
        output_file (str): The path and name of the output file to save the chart.
        
    Returns:
        grouped bar chart saved as png
    """
    reasonable_df = reasonable_df.sort_values(by=["Ancestry", 'Value'])
    perf_df = perf_df.sort_values(by=["Ancestry", 'Value'])

    color_sequence = px.colors.qualitative.Plotly  # You can choose other sequences like D3, G10, etc.

    # Define a color map for ancestries
    colors = {
        'African Ancestry': color_sequence[0],
        'American Ancestry': color_sequence[1],
        'East Asian Ancestry': color_sequence[2],
        'European Ancestry': color_sequence[3],
        'South Asian Ancestry': color_sequence[4]
    }
    reasonable_df_filtered = reasonable_df[reasonable_df['Reference_name'] != "GRCh38"]
    perf_df_filtered = perf_df[perf_df['Reference_name'] != "GRCh38"]

    # Calculate average and standard deviation for single data
    reasonable_avg = reasonable_df_filtered.groupby('Ancestry')['Value'].mean().reset_index()
    reasonable_std = reasonable_df_filtered.groupby('Ancestry')['Value'].std().reset_index()

    # Calculate average and standard deviation for cumulative data
    perf_avg = perf_df_filtered.groupby('Ancestry')['Value'].mean().reset_index()
    perf_std = perf_df_filtered.groupby('Ancestry')['Value'].std().reset_index()

    def relabelAvgPlots(s):
        removed = s.split()[:-1]
        res = ' '.join(removed)
        if res == "European":
            res = "European (w/o GRCh38)"
        return res

    reasonable_avg_labels = [relabelAvgPlots(ancestry) for ancestry in reasonable_avg['Ancestry']]
    perf_avg_labels = [relabelAvgPlots(ancestry) for ancestry in perf_avg['Ancestry']]

    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=("Passing Reasonably Good Threshold", "Reasonably Good Threshold Averages", "Passing Near-Perfect Threshold", "Near-Perfect Threshold Averages"),
        vertical_spacing=0.20, horizontal_spacing=0.05,
        specs=[[{"type": "bar", "colspan": 1}, {"type": "bar", "colspan": 1}],
               [{"type": "bar", "colspan": 1}, {"type": "bar", "colspan": 1}]],
        column_widths=[0.8, 0.2]  # Adjust column widths: 70% for main plots, 30% for average plots
    )


    # Track ancestries to avoid duplicate legend entries
    ancestries_tracked = set()

    # Add bar chart for single_data
    for ancestry in reasonable_df['Ancestry'].unique():
        df = reasonable_df[reasonable_df['Ancestry'] == ancestry]
        show_legend = ancestry not in ancestries_tracked
        ancestries_tracked.add(ancestry)
        fig.add_trace(
            go.Bar(
                x=df['Reference_name'],
                y=df['Value'],
                name=f"{ancestry}",
                marker_color=colors.get(ancestry, 'gray'),  # Default to gray if ancestry not in color map
                showlegend=show_legend
            ),
            row=1, col=1
        )

        # Add annotations for single_data
        for i in range(len(df)):
            rounded_value = round(df['Value'].iloc[i], 1)
            fig.add_annotation(
                x=df['Reference_name'].iloc[i],
                y=df['Value'].iloc[i] + max(df['Value']) * 0.02,
                text=str(rounded_value),
                showarrow=False,
                font=dict(size=10),
                row=1, col=1
            )

    # Add bar chart for cumulative data
    for ancestry in perf_df['Ancestry'].unique():
        df = perf_df[perf_df['Ancestry'] == ancestry]
        show_legend = ancestry not in ancestries_tracked
        ancestries_tracked.add(ancestry)
        fig.add_trace(
            go.Bar(
                x=df['Reference_name'],
                y=df['Value'],
                name=f"{ancestry}",
                marker_color=colors.get(ancestry, 'gray'),  # Default to gray if ancestry not in color map
                showlegend=show_legend
            ),
            row=2, col=1
        )
        # Add annotations for cumulative data
        for i in range(len(df)):
            rounded_value = round(df['Value'].iloc[i], 1)
            fig.add_annotation(
                x=df['Reference_name'].iloc[i],
                y=df['Value'].iloc[i] + max(df['Value']) * 0.02,
                text=str(rounded_value),
                showarrow=False,
                font=dict(size=10),
                row=2, col=1
            )


    # Add average bar chart for single_data with error bars
    fig.add_trace(
        go.Bar(
            x=reasonable_avg_labels,
            y=reasonable_avg['Value'],
            showlegend=False,
            error_y=dict(
                type='data',
                array=reasonable_std['Value'],
                visible=True
            ),
            marker_color=[colors.get(ancestry, 'gray') for ancestry in reasonable_avg['Ancestry']]
        ),
        row=1, col=2
    )

    # Add average bar chart for cumulative data with error bars
    fig.add_trace(
        go.Bar(
            x=perf_avg_labels,
            y=perf_avg['Value'],
            showlegend=False,
            error_y=dict(
                type='data',
                array=perf_std['Value'],
                visible=True
            ),
            marker_color=[colors.get(ancestry, 'gray') for ancestry in perf_avg['Ancestry']]
        ),
        row=2, col=2
    )

    fig.update_layout(
        title_text=overall_title,
        title=dict(
            y=0.99,  # Adjust the title position (1 is the top of the plot area, 0 is the bottom)
            font=dict(size=24, color='black')  # Increase title font size and set to black
        ),
        font=dict(color='black'),  # Set all fonts to black
        barmode='group',
        margin=dict(t=80, b=50, l=50, r=50),  # Adjust margins
        height=950,  # Set the height of the entire figure
        width=1300,  # Set the width of the entire figure
        legend=dict(
            orientation="h",  # Horizontal legend
            yanchor="bottom",  # Anchor the legend to the bottom
            y=1.02,  # Position the legend just above the top of the plot
            xanchor="center",  # Center the legend horizontally
            x=0.5,  # Center the legend at the middle of the plot
            font=dict(size=15)  # Adjust the legend font size
        )
    )

    # Update x and y axes labels
    fig.update_xaxes(title_text=xlabel, tickangle=-45, tickfont=dict(size=12), row=1, col=1)  # First subplot x-axis
    fig.update_xaxes(title_text=xlabel, tickangle=-45, tickfont=dict(size=12), row=2, col=1)  # Third subplot x-axis
    fig.update_yaxes(title_text=ylabel, tickfont=dict(size=12), row=1, col=1)  # First subplot y-axis
    fig.update_yaxes(title_text=ylabel, tickfont=dict(size=12), row=2, col=1)  # Third subplot y-axis

    # Update x and y axes labels for average plots
    fig.update_xaxes(title_text="Ancestry", tickangle=-35, tickfont=dict(size=10), row=1, col=2)  # Second subplot x-axis
    fig.update_xaxes(title_text="Ancestry", tickangle=-35, tickfont=dict(size=10), row=2, col=2)  # Fourth subplot x-axis
    fig.update_yaxes(title_text="Average Value", tickfont=dict(size=12), row=1, col=2)  # Second subplot y-axis
    fig.update_yaxes(title_text="Average Value", tickfont=dict(size=12), row=2, col=2)  # Fourth subplot y-axis

    # Save the plot as a PNG file
    pio.write_image(fig, output_file)





def compare_ref(file_path):
    """
    generate a df containing thresholding results for different reference assemblies in paralel
    
    Parameters:
       file_path (str): path of a tsv file containing columns: "Ancestry", "Path", "Name" 
                        where path = path to alignment csv file for a reference assembly with name = "Name" of ancestry = "Ancestry"
                        (i.e data contained in "csvFilePaths.txt")
                        ** update file paths prior to using
    Returns:
        4 dfs containing reference comparison data:
                            Columns: 'Ancestry', 'Reference_name', 'Value' 
                            (value represents percentage of contigs passing threshold)
                            **(used for group bar plots)
        0.5_0.8CumulPassed.csv = reasonable threshold multiple location analysis
        0.5_0.8SinglePassed.csv = reasonable threshold single location analysis
        0.8_0.9CumulPassed.csv = near-perfect threshold multiple location analysis
        0.8_0.9SinglePassed.csv = near-perfect threshold single location analysis
    """


    ## helper for parallelization
    def process_path(region, path, cov, ident, path_to_name):
        single_entry = {"Ancestry": region, "Reference_name": path_to_name[path]}
        cumul_entry = {"Ancestry": region, "Reference_name": path_to_name[path]}

        df = read_csv(path)
        singlePassed,_,_ = filter(df, cov, ident, cumul=False)
        cumulPassed,_,_ = filter(df, cov, ident, cumul=True)

        single_entry["Value"] = len(singlePassed)/124240 * 100
        cumul_entry["Value"] = len(cumulPassed)/124240 * 100
        return single_entry, cumul_entry

    # Read the file into a DataFrame
    df_paths = pd.read_csv(file_path, sep="\t", header=None, names=["Ancestry", "Path", "Name"])

    # Group by Ancestry and create dictionaries
    group_paths = df_paths.groupby("Ancestry")["Path"].apply(list).to_dict()
    path_to_name = df_paths.set_index("Path")["Name"].to_dict()

    thresholdDict = {(0.5, 0.8): "Reasonably Good", (0.8, 0.9): "Near-Perfect"}
    thresholds = [(0.5, 0.8), (0.8, 0.9)]

    for cov, ident in thresholds:
        single_data = {
            "Ancestry": [],
            "Reference_name": [],
            "Value": []
        }
        cumul_data = {
            "Ancestry": [],
            "Reference_name": [],
            "Value": []
        }

        with concurrent.futures.ThreadPoolExecutor() as executor:
            futures = []
            for region in group_paths:
                for path in group_paths[region]:
                    futures.append(executor.submit(process_path, region, path, cov, ident, path_to_name))

            for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures)):
                single_entry, cumul_entry = future.result()
                for key in single_entry:
                    single_data[key].append(single_entry[key])
                    cumul_data[key].append(cumul_entry[key])

        single_df = pd.DataFrame(single_data).sort_values(by=["Ancestry", 'Value'])
        cumul_df = pd.DataFrame(cumul_data).sort_values(by=["Ancestry", 'Value'])

        create_csv(single_df,filename=str(cov)+'_'+str(ident)+'SinglePassed.csv')
        create_csv(cumul_df,filename=str(cov)+'_'+str(ident)+'CumulPassed.csv')



def generateComparisonBars():
    """
    wrapper function to calculate grouped bar plots based on compare_ref results
    
    returns:
        two grouped bar plots:
            CumulPassedchart.png = multiple location analysis plots
            SinglePassedchart.png = single location analysis plots

    """


    thresholds = [(0.5, 0.8), (0.8, 0.9)]
    methods = ["SinglePassed", "CumulPassed"]
    methodsDict = {"SinglePassed": "Single Location Analysis", "CumulPassed":"Multiple Location Analysis"}

    for method in methods:
        dfs = []
        for cov, ident in thresholds:
            dfs.append(read_csv("/grid/chambwe/home/dcha/APGmapping/figures/barcharts/"+str(cov)+'_'+str(ident)+method+".csv"))

        reasonable_df, perf_df = dfs
        calcKruskal(reasonable_df)
        calcKruskal(perf_df)
        create_combined_plot_plotly(reasonable_df, perf_df, overall_title='Percentage of Contigs Passing Thresholds ({})'.format(methodsDict[method]), xlabel='Reference Assemblies', ylabel='Percentage of Passing Contigs', output_file=method + 'chart.png')



def generateRefCSV(sam_file_path, output_filename):
    """
    wrapper function to create csv of create_contig_df() results
    
    returns:
        csv file containing containing reference specific contig alignment results.
    """
    df = create_contig_df(sam_file_path, ref_seq_path, hg38_alias_path, contig_name_path, contig_lengths_path, CPG_bedFile_path_masked, unique_bedFile_path)
    create_csv(df, filename=output_filename)



# update these file paths prior to usage

contigs_path = "/grid/chambwe/macmillan/african_pangenome/apg_contigs/PDBU01.1.fsa_nt"
contig_lengths_path = "/grid/chambwe/macmillan/african_pangenome/apg_contigs/lengths.csv"
unique_bedFile_path = "/grid/chambwe/macmillan/T2T/uniqueT2T/hgUnique.hg38.bb"
ref_seq_path = "/grid/chambwe/macmillan/T2T/sequence_report.tsv"
contig_name_path = "/grid/chambwe/macmillan/african_pangenome/apg_contigs/PDBU01_contigs.tsv"
hg38_alias_path = "/grid/chambwe/macmillan/GRCh38/hg38.chromAlias.txt"

CPG_bedFile_path_masked = "/grid/chambwe/macmillan/T2T/cpgIslandExt.bb"
CPG_bedFile_path_unmasked = "/grid/chambwe/macmillan/T2T/cpgIslandExtUnmasked.bb"





#################################### example usage #####################################



################### generate CSV Files in parallel ################
'''
sam_file_paths = []
with open('/grid/chambwe/home/dcha/APGmapping/APG_HPRC/HPRCSamPaths.txt', 'r') as file:
    sam_file_paths = [line.strip() for line in file if line.strip()]

# Create contig_data with corresponding CSV file names
contig_data = [(sam_path, f"{sam_path.split('/')[-1].split('.')[0]}.csv") for sam_path in sam_file_paths]

print(contig_data)

# Use ProcessPoolExecutor for multiAprocessing
with concurrent.futures.ProcessPoolExecutor() as executor:
    futures = [executor.submit(generateRefCSV, sam_path, output) for sam_path, output in contig_data]

    # Use tqdm to show the progress bar
    for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures)):
        try:
            future.result()  # This will raise an exception if the callable raised
        except Exception as e:
            print(f"An error occurred: {e}")
'''


######## generate scatter plot csv and plots ################
'''
scatterData('/grid/chambwe/home/dcha/APGmapping/alignmentCSVs/GRCh38.csv', cumul=False, filename = "GRCh38scatterData.csv")
scatterData('/grid/chambwe/home/dcha/APGmapping/alignmentCSVs/T2T.csv', cumul=False, filename = "T2TscatterData.csv")

draw_scatterplot('T2TscatterData.csv', title="APG Contigs Alignment Scores to T2T-CHM13(Single Location)", xlabel="Coverage", ylabel="Identity", filename='Heat_T2T_Contig_Scatter.png')
draw_scatterplot('GRCh38scatterData.csv', title="APG Contigs Alignment Scores to GRCh38(Single Location)", xlabel="Coverage", ylabel="Identity", filename='Heat_GRCh38_Contig_Scatter.png')

'''

####### generate histogram ################################
'''
df = read_csv("/grid/chambwe/home/dcha/APGmapping/alignmentCSVs/T2T.csv")
cpg_identity, cpg_coverage, noCpG_coverage, noCPG_identity = cov_id_distribution(df)

get_histo(cpg_identity, noCPG_identity, "identity_histo.png", title="Distribution of Alignment Identities", step=0.1, label1="Alignments Containing CpG Islands", label2="Alignments w/o CpG Islands", xlabel="Identity Intervals")
get_histo(cpg_coverage, noCpG_coverage, "coverage_histo.png", title="Distribution of Alignment Coverages", step=0.1, label1="Alignments Containing CpG Islands", label2="Alignments w/o CpG Islands", xlabel="Coverage Intervals")


'''
