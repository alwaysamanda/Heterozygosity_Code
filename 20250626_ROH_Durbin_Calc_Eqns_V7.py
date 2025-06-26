#### WRITING DURBIN'S EQUATIONS TO CALCULATE ROH FOR VGP GENOMES ####
## Date: 20250626 (June 26th, 2025)
## Author: Amanda Gardiner
## Version 7 (V5 is 20250320_ROH_Durbin_Calc_Eqns_V6.py)
## NOTES: 
##          Start with X0 = 0, a0 = b0 = c0 = S0 = 0
##          ai = start, bi = end, ci = difference
##          Move along chromosome between SNPs
##          Si = Si-1 + (Xi - Xi-1) - P
##          Si-1 = score before going to this region
##          Xi - Xi-1 = distance between this SNP and previous SNP
##          P = set penalty (try 1,000,000 and then 100,000)
##          If Si ≤ 0, then set ai = Xi, bi = Xi, ci = 0
##          If Si > 0, then if Si > ci, make ci = Si, and bi = Xi
##          Repeat along whole chromosome to get regions of interest

## NOTES: Doing this to modify input from ALN file to account for swapped position
##        of Ref Chr names in Aln files using Richard Durbin's FastGA alingments of VGP data freeze
####

#### ---- Import necessary libraries ---- ####
import numpy as np
import pandas as pd
import re as re
import sys
from sys import stdin
import io
import os

#### ----  Load in variables from the shell script ---- ####
chrom = sys.argv[1] ## Chromosome ROH is being calculated for
chrom_length = int(sys.argv[2]) ## Chromosome length in bases

ref_name = sys.argv[3] ## Reference fasta sequence name 
clade = sys.argv[4] ## Clade that species belongs to

Aln_file = sys.argv[5]

Var_file = sys.argv[6]

spec_name = sys.argv[7]

#### ----  FUNCTIONS ---- ####
## Define function to read the alignment file and store it in an array
def read_aln_file(Alignment_file, chromosome):
    with open(Alignment_file, 'r') as file:
        dat = [line.split() for line in file if line.split()[4] == chromosome]
    return dat

aln_list = read_aln_file(Aln_file, chrom)


## Define function to parse through whole chromosome and map each individual base and whether it is aligned or not
def aln_map(Alignment_list, chrom_length):
    base_sums = np.zeros(chrom_length, dtype=int) ## Create an array the length of the chromosome
    for start, end in ((int(line[2]), int(line[3])) for line in Alignment_list):
        base_sums[start:end] += 1 ## Calculate depth of alingments at each individual base
    base_gross_cov = np.where(base_sums > 0, 1, 0) ## Convert depth to binary of whether it is covered or not
    aln_sums = np.cumsum(base_gross_cov) ## Get alingment sum at each given base
    return np.vstack((np.arange(1, chrom_length + 1), base_sums, aln_sums.astype(int)))
## Finish function

base_aln_map = aln_map(aln_list, chrom_length)

## Define function to read the variant file and store it in an array
def read_var_file(variant_file, chromosome):
    with open(variant_file, 'r') as file:
        var_lines = [line.split() for line in file if line.split()[8] == chromosome]
    var_pos = [int(line[2]) for line in var_lines]
    return var_lines, var_pos

var_list, var_pos = read_var_file(Var_file, chrom)


## Define function to calculate ROH
def calculate_ROH(var_pos, base_aln_map, chrom):
    score = np.zeros(len(var_pos), dtype=int) ## Create array to store the score at each base position
    aln_sums = np.zeros(len(var_pos), dtype=int) ## Create array to store the total number of aligned bases at each position
    df = pd.DataFrame({'variant': var_pos, 
                       'var_aln_sum': aln_sums, 
                       'score': score})
    indices = np.searchsorted(base_aln_map[0, :], df['variant']) - 1 ## Find the variants in the base_aln_map and subtract one to get them to 0-indexing
    df['var_aln_sum'] = base_aln_map[2, indices] ## var_aln_sum == sum of all alingned bases at the variant's point
    ## Run equation to calculate score based on aligned chr length (var_aln_sum)
    Penalty = 1e+05 ## Define penalty -- currently 100,000
    df['score'] = np.maximum(0, (np.roll(df['score'].astype(int), 1) + (df['var_aln_sum'].astype(int) - np.roll(df['var_aln_sum'], 1)) - Penalty))
    ## Extract the start and end points for each ROH
    df['ROH_true'] = np.maximum(df['score'], 0) > 0
    mask = df['ROH_true']
    groups = (mask != mask.shift()).cumsum()
    true_groups = groups[mask].unique()
    result = pd.Series(true_groups).apply(
        lambda group: (
            df[groups == group]['variant'].iloc[0] 
            if df[groups == group].index[0] == 0 
            else df.loc[df[groups == group].index[0] - 1, 'variant'], 
            df[groups == group]['variant'].iloc[-1]
        )
    )
    ## Use start and end points to create a data frame containing ROH start, end, and length
    ROH_start = [roh[0] for roh in result]
    ROH_end = [roh[1] for roh in result]
    df_roh = pd.DataFrame({
        'start': ROH_start, 
        'end': ROH_end
    })
    df_roh['length'] = df_roh['end'] - df_roh['start']
    chrom_list = [chrom]*df_roh.shape[0]
    df_roh.insert(loc = 0, column = 'chrom', value=chrom_list)
    return df_roh

ROH_report = calculate_ROH(var_pos, base_aln_map, chrom)
# print(ROH_report)

if ROH_report.empty == True:
    print('No ROH detected in', chrom)

output_name = os.path.join(clade, spec_name, chrom + "_ROH_Results.txt")
ROH_report.to_csv(output_name, index=False, header=True)
