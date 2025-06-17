#### Script to calculate FROH (Inbreeding Coefficient from ROH Durbin Calculations ####
## Date: 20250617 (June 17th, 2025)
## Author: Amanda Gardiner
## Version 2 (Version 1 is 20250210_FROH_per_chr_calc.R)
## GOAL: Calculate The Inbreeding coefficient per chromosome for an individual,
##       and report results both in orignial format and normalized for the size of the chromosome so that they can be compared
## NOTES: Initial script pulled from 20250106_FROH_Calc.R
## NOTES: Each time I am running this script, I am running it for a specific chromosome -- as such I do not need to loop through all of them
## NOTES: Altering this from the previous script to get the array of aligned bases and calculate l_aut and chr length based on aligned bases only


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
Roh_file = sys.argv[2]
clade = sys.argv[3]
spec_name = sys.argv[4]
chrom = sys.argv[5]
chrom_length = int(sys.argv[6])

output_file_name = os.path.join(clade, spec_name, chrom+"_FROH_results.txt")


#### ---- FUNCTIONS ---- ####
## Define function to read the alignment file and store it in an array
def read_aln_file(Alignment_file, chromosome):
    with open(Alignment_file, 'r') as file:
        dat = [line.split() for line in file if line.split()[1] == chromosome]
    return dat

aln_list = read_aln_file(Aln_file, chrom)
# print(aln_list)

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
# print(base_aln_map)

## Define function to read the ROH data
def read_ROH_dat(roh_file_name):
    with open(roh_file_name, 'r') as file:
        dat = [line.strip().split(',') for line in file]
        raw_data = pd.DataFrame(dat[1:], columns=dat[0])
        raw_data['length'] = pd.to_numeric(raw_data['length'], errors='coerce').astype('Int64')
    return raw_data
## Finish function
roh_dat = read_ROH_dat(Roh_file)


## Define function to calculate FROH
def calc_FROH(alingment_map, Roh_data, chromosome, file_name):
    chr_roh_dat = Roh_data[Roh_data["chrom"] == chromosome]
    roh_sum_all = chr_roh_dat['length'].sum() ## Sum of all ROH of all lengths
    roh_sum_med_long = chr_roh_dat[chr_roh_dat['length'] >= 500000]['length'].sum() ## Sum of medium and long ROH lengths
    roh_sum_long = chr_roh_dat[chr_roh_dat['length'] >= 1000000]['length'].sum() ## Sum of long ROH lengths
    aln_length = alingment_map[2,-1] ## Get length of chromosome with only aligned bases
    Froh_all = roh_sum_all/aln_length
    Froh_all_percent = Froh_all*100
    Froh_med_long = roh_sum_med_long/aln_length
    Froh_med_lon_percent = Froh_med_long*100
    Froh_long = roh_sum_long/aln_length
    Froh_long_percent = Froh_long*100
    with open(file_name, 'wt') as file:
        file.write(f"Froh_aut value for all autosomal ROH (short, med, long) is: \n")
        file.write(f"{Froh_all} \n")
        file.write(f"Froh_aut in percent for all autosomal ROH (short, med, long) is: \n")
        file.write(f"{Froh_all_percent} \n \n")
        file.write(f"Froh_aut value for medium and long autosomal ROH is: \n")
        file.write(f"{Froh_med_long} \n")
        file.write(f"Froh_aut in percent for medium and long autosomal ROH is: \n")
        file.write(f"{Froh_med_lon_percent} \n \n")
        file.write(f"Froh_aut value for long autosomal ROH is: \n")
        file.write(f"{Froh_long} \n")
        file.write(f"Froh_aut in percent for long autosomal ROH is: \n")
        file.write(f"{Froh_long_percent} \n")
        file.close()
## End function

calc_FROH(base_aln_map, roh_dat, chrom, output_file_name)
