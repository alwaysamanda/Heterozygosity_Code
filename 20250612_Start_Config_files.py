#### Script to create config files for all VGP species ####
## Date: 20250612 (June 12th, 2025)
## Author: Amanda Gardiner
## Version 1
## GOAL: Create template of config files for all VGP species with alternate and secondary haplotypes
####

###############################################################################

#### ---- Import necessary libraries ---- ####
import numpy as np
import pandas as pd
import re as re
import sys
from sys import stdin
import io
import os

#### ---- Load in data ---- ####
dat_file = sys.argv[1]

def READ_FILE(file_name):
    return pd.read_csv(file_name, sep='\t')

data = READ_FILE(dat_file)


#### ---- Create config files ---- ####
def GET_SPEC_NAME(string):
    sep_name=string.split("_")
    name1=sep_name[0][0:3]
    name2=sep_name[1][0:3].capitalize()
    fullname = name1+name2
    return(fullname)

def CREATE_CONFIG(data):
    folder = "config_files"
    for row in range(len(data)):
        clade=data.iloc[row,1]
        species_name=data.iloc[row,2]
        shortname=GET_SPEC_NAME(species_name)
        ref_seq=data.iloc[row,3]
        sec_seq=data.iloc[row,4]
        output_file_name = os.path.join(folder, clade, shortname + "_config.yml")
        if os.path.exists(output_file_name):
            print(f"File already exists: {output_file_name} — Skipping.")
            continue
        with open(output_file_name, "wt") as file:
            file.write(f"    CLADE : {clade} \n")
            file.write(f"    SPEC_NAME : {shortname} \n")
            file.write(f"    REF_NAME : {ref_seq} \n")
            file.write(f"    ALT_NAME : {sec_seq} \n")
            file.write(f"    CHROM_START_CHR : \n")
            file.write(f"    WINDOW_INTERVAL : 500000 \n")
            file.write(f"    WINDOW_LENGTH : 1000000 \n")
            file.write(f"    NUM_AUT_CHROMOSOMES : \n")
            file.write(f"    NUM_ALL_CHR: \n")
            file.write(f"    MUTATION_RATE: 1.25e-8 \n")
            file.write(f"    GENERATION_TIME: \n")
            file.write(f"    ALL_CHROMOSOMES : \n")
            file.write(f"        - \n")
            file.write(f"    AUTOSOMAL_CHROMOSOMES: \n")
            file.write(f"        - \n")
            file.write(f"    BOOTSTRAPPING_VALUES: \n")
            for i in range(30):
                file.write(f"        - {i} \n")
            file.close()

CREATE_CONFIG(data)
