#### Script to calculate heterozgyosity based on variants in vcf file ####
## Date: 20250319 (March 19th, 2025)
## Author: Amanda Gardiner
## Version 3 (V2 is 20250319_find_het_calc_V2.R)
## NOTES: This is one part of V3 -- split the original script into two different scripts
##        Did to make code usable within Snakemake snakefile
##        Partner script is 20250325_find_het_per_chr_V3.py
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
chr_file = sys.argv[1] # File containing chromosome names and their lengths
clade = sys.argv[2] # Clade name
species = sys.argv[3] #Species name/reference genome name -- directory for storing results
date = sys.argv[4] # Date of running file
num_aut_chr = int(sys.argv[5]) # Number of autosomal chromosomes in the genome

#### FUNCTIONS ####
def filter_var_by_chrom(file_of_variants, chromosome):
    with open(file_of_variants, 'r') as file:
        dat = [line.split() for line in file if len(line.split()) > 1 and line.split()[1] == chromosome]
    return pd.DataFrame(dat)

def filter_chr_file_by_chrom(chromosome_file, chromosome):
    with open(chromosome_file, 'r') as file:
        dat = [word for line in file if line.split()[0] == chromosome for word in line.split()]
    return dat

def read_chromosome_file(chromosome_file):
    with open(chromosome_file, 'r') as file:
        dat = [line.split()[:] for line in file]
    return pd.DataFrame(dat, columns=["chr_name", "chr_length"])

def filter_roh_by_chrom(roh_data, chromosome):
    with open(roh_data, 'r') as file:
        dat = [line.split() for line in file if line.split()[0] == chromosome]
    return pd.DataFrame(dat)

def calc_het(chromosome_length_file, clade, species, date, num_aut_chr):
    ## Read in chromosomes and their lengths
    chr_df = read_chromosome_file(chromosome_length_file)

    ## Create dataframe to track mean results for each chromosome
    all_chromosomes = chr_df['chr_name']
    chr_mean_het = np.zeros(len(all_chromosomes))
    chr_mean_het_excl_ROH = np.zeros(len(all_chromosomes))
    mean_het_df = pd.DataFrame({
        'chr': all_chromosomes, 
        'mean_het': chr_mean_het, 
        'mean_het_excl_ROH': chr_mean_het_excl_ROH, 
    })

    all_het = []
    all_auto_het = []

    all_het_excl_ROH = []
    all_auto_het_excl_ROH = []

    ## Find autosomal chromosomes
    aut_chr = chr_df.iloc[0:num_aut_chr, 0]
    # print(aut_chr)

    for chr_name in chr_df.iloc[:,0]: ## Read in chromosome files
        file_path = os.path.join(clade, species, date + '_' + chr_name + '_het.txt')
        chrom_file = pd.read_csv(file_path, sep=',', header=0)

        ## Calculate mean het per kb for the chromosome and store it in mean_het_df
        mean_het_df.loc[mean_het_df['chr']==chr_name,'mean_het'] = np.mean(chrom_file['Het_Per_Kb'])
        mean_het_df.loc[mean_het_df['chr']==chr_name,'mean_het_excl_ROH'] = np.mean(chrom_file['Het_Per_Kb_excl_ROH'])

        ## Calculate mean het per kb window for the whole genome and just autosomal genome
        all_het.extend(chrom_file['Het_Per_Kb'].tolist())
        all_het_excl_ROH.extend(chrom_file['Het_Per_Kb_excl_ROH'].tolist())
    ## Finish for loop

    ## Save file containing mean heterozygosty per kb window per chromosome
    mean_het_file_name = os.path.join(clade, species, date + '_' + species + '_per_chr_mean_heterozygosity.txt')
    mean_het_df.to_csv(mean_het_file_name, header=True, index=False)
    # print(mean_het_df)

    total_mean = np.mean(all_het)
    total_mean_excl_ROH = np.mean(all_het_excl_ROH)

    for chr_name in aut_chr:
        file_path = os.path.join(clade, species, date + '_' + chr_name + '_het.txt')
        chrom_file = pd.read_csv(file_path, sep=',')
        all_auto_het.extend(chrom_file['Het_Per_Kb'].tolist())
        all_auto_het_excl_ROH.extend(chrom_file['Het_Per_Kb_excl_ROH'].tolist())
    ## Finish for loop
    total_auto_mean = np.mean(all_auto_het)
    total_auto_mean_excl_ROH = np.mean(all_auto_het_excl_ROH)

    ## Save results
    results = [
        'The mean heterozygosity per 1kb over the whole genome: \n', total_mean, '\n', 
        'The mean heterozygosity per 1kb over the whole autosomal genome: \n', total_auto_mean, '\n', 
        'The mean heterozygosity per 1kb over the whole genome excluding ROH: \n', total_mean_excl_ROH, '\n',
        'The mean heterozygosity per 1kb over the whole autosomal genome excluding ROH: \n', total_auto_mean_excl_ROH 
    ]
    total_mean_file_name = os.path.join(clade, species, date  + '_' + species + '_whole_genome_mean_heterozygosity.txt')

    result_file = open(total_mean_file_name, 'w')
    result_file.writelines(str(element) for element in results)
    result_file.close()
## End function

run_het_calculations = calc_het(chr_file, clade, species, date, num_aut_chr)
