#### Script to calculate FROH across the whole genome (Inbreeding Coefficient from ROH Durbin Calculations ####
## Date: 20250617 (June 17th, 2025)
## Author: Amanda Gardiner
## Version 3 (Version 2 is 20250603_FROH_Calc_Whole_Genome_V2.R)
## NOTES: Based on rds/users/ag2427/hpc-work/ROH_analyses/20241211_FROH_ROH_Durbin_Calc.R
## NOTES: Doing this to calculate it based on chromosome length from only aligned bases 

####

#### ---- Import necessary libraries ---- ####
import numpy as np
import pandas as pd
import re as re
import sys
from sys import stdin
import io
import os

#### ----- Load in variables from the shell ---- ####
Aln_file = sys.argv[1]
num_aut_chr = int(sys.argv[2])
chromosome_list = sys.argv[3]
roh_data_file = sys.argv[4]
output_file_name = sys.argv[5]

#### ---- FUNCTIONS ---- ####
## Define function to read the alignment file and store it in an array
def read_aln_file(Alignment_file, chromosome):
    with open(Alignment_file, 'r') as file:
        dat = [line.split() for line in file if line.split()[1] == chromosome]
    return dat

## Define function to parse through whole chromosome and map each individual base and whether it is aligned or not
def aln_map(Alignment_list, chrom_length):
    base_sums = np.zeros(chrom_length, dtype=int) ## Create an array the length of the chromosome
    for start, end in ((int(line[2]), int(line[3])) for line in Alignment_list):
        base_sums[start:end] += 1 ## Calculate depth of alingments at each individual base
    base_gross_cov = np.where(base_sums > 0, 1, 0) ## Convert depth to binary of whether it is covered or not
    aln_sums = np.cumsum(base_gross_cov) ## Get alingment sum at each given base
    return np.vstack((np.arange(1, chrom_length + 1), base_sums, aln_sums.astype(int)))
## Finish function

def read_chromosome_dat(chromosome_list_file, number_autosomal):
    Chrom = pd.read_csv(chromosome_list_file, header=None, sep='\t')
    Chrom.columns = ['Chrom_name', 'Chrom_end']
    Chrom_autosomal = Chrom.iloc[0:num_aut_chr, 0:2]
    return Chrom_autosomal
## Finish function
chrom_dat = read_chromosome_dat(chromosome_list, num_aut_chr)

def get_roh_sum(roh_file_name):
    try:
        dat = pd.read_csv(roh_file_name, header=None)
        dat.columns = ['chrom', 'start', 'end', 'length']
        dat['length'] = pd.to_numeric(dat['length'], errors='coerce').astype('Int64')
        roh_sum_all_lengths = dat['length'].sum()
        roh_sum_med_long = dat[dat['length'] >= 500000]['length'].sum()  ## Sum of medium and long ROH lengths
        roh_sum_long = dat[dat['length'] >= 1000000]['length'].sum() ## Sum of long ROH lengths
    except pd.errors.EmptyDataError: ## If there are no ROH found
        print('Note: ROH csv was empty. Skipping.')
        dat = [0]
        roh_sum_all_lengths = 0
        roh_sum_med_long = 0  
        roh_sum_long = 0 
    all_roh_sums = [roh_sum_all_lengths, roh_sum_med_long, roh_sum_long]
    return all_roh_sums
## Finish function
Lroh = get_roh_sum(roh_data_file)


def calc_FROH_aut(sum_roh_lengths, list_all_chromosomes, Aln_file, file_name):
    all_aln_length = []
    for row in range(len(list_all_chromosomes)):
        chrom = list_all_chromosomes.iloc[row,0]
        chrom_length = list_all_chromosomes.iloc[row,1]
        aln_list = read_aln_file(Aln_file, chrom)
        base_aln_map = aln_map(aln_list, chrom_length)
        aln_length = base_aln_map[2,-1]
        all_aln_length.append(aln_length)
    sum_aln_length = sum(all_aln_length)
    ## Divide each Lroh by sum
    Froh_all = sum_roh_lengths[0]/sum_aln_length
    Froh_all_percent = Froh_all*100
    Froh_med_long = sum_roh_lengths[1]/sum_aln_length
    Froh_med_long_percent = Froh_med_long*100
    Froh_long = sum_roh_lengths[2]/sum_aln_length
    Froh_long_percent = Froh_long*100
    ## Multiply each by 100 to get the percentage
    ## Save in a file
    with open(file_name, 'wt') as file:
        file.write(f"The Froh value for all ROH (short, med, long) is: \n")
        file.write(f"{Froh_all} \n")
        file.write(f"The Froh in percent for all ROH (short, med, long) \n")
        file.write(f"{Froh_all_percent} % \n \n")
        file.write(f"The Froh value for medium and long ROH (med, long) is: \n")
        file.write(f"{Froh_med_long} \n")
        file.write(f"The Froh in percent for medium and long ROH (med, long) is: \n")
        file.write(f"{Froh_med_long_percent} % \n \n")
        file.write(f"The Froh value for long ROH(long) is: \n")
        file.write(f"{Froh_long} \n")
        file.write(f"The Froh in percent for long ROH (long) \n")
        file.write(f"{Froh_long_percent} % \n")
        file.close()
## Finish Function

calc_FROH_aut(Lroh, chrom_dat, Aln_file, output_file_name)
