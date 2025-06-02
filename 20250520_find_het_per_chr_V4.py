#### Script to calculate heterozgyosity based on variants in vcf file ####
## Date: 20250520 (May 20th, 2025)
## Author: Amanda Gardiner
## Version 4 (V3 is 20250325_find_het_whole_genome_V3.py)
## NOTES: This is one part of V3 -- split the original script into two different scripts
##        Did to make code usable within Snakemake snakefile
##        Partner script is 20250325_find_het_per_chr_V3.py
##        UPDATING THIS TO WORK FOR ONE CHROMOSOME AT A TIME instead of looping over all chromosomes
## GOAL: Take input of txt file containing chromosomes and variant positions and use it to calculate heterozygosity per chromosome and over whole genome
####

###############################################################################

## Import necessary libraries
import numpy as np
import pandas as pd
import re as re
import sys
from sys import stdin
import io
import os

## Load in variables from shell
dat = sys.argv[1] # File containing variant information
chr_file = sys.argv[2] # File containing chromosome names and their lengths
clade = sys.argv[3] # Clade name
species = sys.argv[4] #Species name/reference genome name -- directory for storing results
current_window_length = int(sys.argv[5]) # Window length defined for calculating het
current_window_interval = int(sys.argv[6]) # How often windows start, defines their level of overlap
chrom = sys.argv[7] # Number of autosomal chromosomes in the genome
roh_data = sys.argv[8] # File containing ROH start, end, and length for whole genome



#### FUNCTIONS ####
def filter_var_by_chrom(file_of_variants, chromosome):
    with open(file_of_variants, 'r') as file:
        dat = [line.split() for line in file if len(line.split()) > 1 and line.split()[1] == chromosome]
    return pd.DataFrame(dat)

def filter_chr_file_by_chrom(chromosome_file, chromosome):
    with open(chromosome_file, 'r') as file:
        dat = [word for line in file if line.split()[0] == chromosome for word in line.split()]
    return dat

def filter_roh_by_chrom(roh_data, chromosome):
    dat=[]
    with open(roh_data, 'r') as file:
        for line in file:
            splitline = line.split(',')
            if splitline[0] == chromosome:
                splitline[3] = splitline[3].replace('\n', '')
                dat.append(splitline)
    return pd.DataFrame(dat)

def calc_het(variant_file, roh_file, chromosome_length_file, window_length, window_interval, chromosome):
    single_chr_df = filter_var_by_chrom(variant_file, chromosome)
    chr_info = filter_chr_file_by_chrom(chromosome_length_file, chromosome)
    chr_length = int(chr_info[1])
    if chr_length < window_length:
        window_length_alt = 5000
        window_interval_alt = window_length_alt/2
        window_starts = np.arange(0, (chr_length-window_length_alt+1), window_interval_alt)
        window_ends = window_starts + window_length_alt
        window_sizes = np.repeat(window_length_alt, len(window_starts))
    else:
        window_starts = np.arange(0, (chr_length-window_length+1), window_interval)
        window_ends = window_starts + window_length 
        window_sizes = np.repeat(window_length, len(window_starts))
    ## End if statement
    het = np.zeros(len(window_starts)) ## Create heterozygosity vector
    het_wo_roh = np.zeros(len(window_starts)) ## heterozygosity vector excluding ROH
    single_chr_results = pd.DataFrame({
        'Start': window_starts, 
        'End': window_ends, 
        'Het': het, 
        'Het_excl_ROH': het_wo_roh, 
        'Window_Size': window_sizes, 
        'Window_Size_excl_ROH': window_sizes
    }) ## Create dataframe to store results for the given chromosome
    ## Count the number of variants in each window
    starts = single_chr_results.iloc[:, 0].astype(int)
    ends = single_chr_results.iloc[:, 1].astype(int)
    variant_positions = single_chr_df.iloc[:, 2].astype(int)
    single_chr_results.loc[:,'Het'] = [(variant_positions.between(start, end - 1)).sum() for start, end in zip(starts, ends)] ## Sum all variants in each window, not worrying about ROH
    ## Count number of variants in each window, excluding those found in ROH
    ## Get ROH positions
    roh_dat = filter_roh_by_chrom(roh_file, chromosome)
    ## If there are ROH present on chromosome -- Calculate het
    if roh_dat.empty == False:
        ROH_starts = (roh_dat.iloc[:, 1].values).astype(int) 
        ROH_ends = (roh_dat.iloc[:, 2].values).astype(int)    
        variant_positions = (single_chr_df.iloc[:, 2].values).astype(int)  
        for k in range(len(single_chr_results)):  
            start, end = single_chr_results.iloc[k, [0, 1]].astype(int) ## Start and end of a given window
            var_in_win = variant_positions[(variant_positions >= start) & (variant_positions < end)] ## Find all variants within a given window
            ## Find the first ROH that starts after 'start'
            # l = int(np.searchsorted(ROH_starts, start, side='right') - 1)
            l = np.where(ROH_starts > start)[0]
            if l.size > 0:
                l = int(l[0])
            else:
                l = -1
            if l < 0 or end < ROH_starts[l]:  
                ## No ROH overlap, count all variants in the range
                sum_variants = len(var_in_win)
                single_chr_results.loc[k,'Het_excl_ROH'] = sum_variants
            elif start >= ROH_starts[l] and end <= ROH_ends[l]:  
                ## Window fully inside ROH
                sum_variants = 0  
                single_chr_results.loc[k,'Het_excl_ROH'] = sum_variants
            elif start < ROH_starts[l] and end > ROH_starts[l] and end <= ROH_ends[l]:  
                ## Window overlaps start of ROH
                var_in_win_excl_ROH = var_in_win[(var_in_win < ROH_starts[l])]
                sum_variants = len(var_in_win_excl_ROH)
                single_chr_results.loc[k,'Het_excl_ROH'] = sum_variants  
                single_chr_results.loc[k,'Window_Size_excl_ROH'] = ROH_starts[l] - start
            elif start >= ROH_starts[l] and start <= ROH_ends[l] and end > ROH_ends[l]:  
                ## Window overlaps end of ROH
                var_in_win_excl_ROH = var_in_win[(var_in_win > ROH_ends[l])]
                sum_variants = len(var_in_win_excl_ROH)
                single_chr_results.loc[k,'Het_excl_ROH'] = sum_variants
                single_chr_results.loc[k,'Window_Size_excl_ROH'] = end - ROH_ends[l]
            elif start < ROH_starts[l] and end > ROH_ends[l]:  
                ## Window Completely contains ROH
                variants_pre_ROH = var_in_win[(var_in_win >= start) & (var_in_win < ROH_starts[l])]
                sum_variants_pre_ROH = len(variants_pre_ROH)
                variants_post_ROH = var_in_win[(var_in_win > ROH_ends[l]) & (var_in_win < end)]
                sum_variants_post_ROH = len(variants_post_ROH)
                sum_variants = sum_variants_pre_ROH + sum_variants_post_ROH
                single_chr_results.loc[k,'Het_excl_ROH'] = sum_variants
                single_chr_results.loc[k,'Window_Size_excl_ROH'] = (ROH_starts[l] - start) + (end - ROH_ends[l])
        ## Finish for loop
    elif roh_dat.empty == True: ## Calculate het if there are NO ROH on chr
        single_chr_results.loc[:,'Het_excl_ROH'] = single_chr_results.loc[:,'Het'] ## Het_excl_ROH is same as het given that there are no ROH on chromosome
    ## Calculate heterozygosity per kb
    single_chr_results['Het_Per_Kb'] = (single_chr_results['Het']/single_chr_results['Window_Size'])*1000
    single_chr_results['Het_Per_Kb_excl_ROH'] = (single_chr_results['Het_excl_ROH']/single_chr_results['Window_Size_excl_ROH'])*1000
    ## Add chromosome column
    chrom_list=[chromosome]*single_chr_results.shape[0]
    single_chr_results.insert(loc = 0, column = 'chrom', value=chrom_list)
    ## Save results for chromosome
    chr_file_name = os.path.join(clade, species, chromosome + '_het.txt')
    single_chr_results.to_csv(chr_file_name, index=False, header=True)
## End function

run_het_calculations = calc_het(dat, roh_data, chr_file, current_window_length, current_window_interval, chrom)
