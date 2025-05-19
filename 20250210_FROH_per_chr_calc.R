#### Script to calculate FROH (Inbreeding Coefficient from ROH Durbin Calculations ####
## Date: 20250210 (February 10th, 2025)
## Author: Amanda Gardiner
## Version 1
## GOAL: Calculate The Inbreeding coefficient per chromosome for an individual,
##       and report results both in orignial format and normalized for the size of the chromosome so that they can be compared
## NOTES: Initial script pulled from 20250106_FROH_Calc.R
## NOTES: Each time I am running this script, I am running it for a specific chromosome -- as such I do not need to loop through all of them


####

## Load in necessary packages
library(dplyr)
args <- commandArgs()

## Load in required variables
## Load in csv file containing ROH
ROH_dat <- read.csv(args[6], header=FALSE)
colnames(ROH_dat) <- c('Chr', 'Start', 'End', 'Length')
# print(head(ROH_dat))

Today_date <- args[8]
Ref_name <- args[9] ## Reference genome name
clade <- args[10]
num_aut_chr <- as.numeric(args[11]) ## Number of autosomal chromosomes
chrom_length <- as.numeric(args[12]) ## Length of chromosome from first to last measured variant (will be slightly shorter than total lenth of chromosome)
chr <- args[13]

## Calculate FROH for whole chromosome, all ROH 
chr_ROH <- ROH_dat %>% filter(ROH_dat$Chr == chr) ## Select ROH for that chromosome
L_ROH <- sum(chr_ROH$Length)

FROH_all <- L_ROH/chrom_length
# Report as a percentage
FROH_all_percent <- FROH_all*100

print("Froh_aut value for all autosomal ROH (short, med, long) is:")
print(FROH_all)
print(" ")
print("Froh_aut in percent for all autosomal ROH (short, med, long)")
print(paste((FROH_all_percent), "%"))
print("")

## Calculate FROH for medium and long ROH only
chr_ROH_med_long <- chr_ROH %>% filter(chr_ROH$Length >= 500000)
L_ROH_med_lon <- sum(chr_ROH_med_long$Length)
FROH_med_long <- L_ROH_med_lon/chrom_length
# Report as a percentage
FROH_med_long_percent <- FROH_med_long*100

print("Froh_aut value for medium and long autosomal ROH is:")
print(FROH_med_long)
print(" ")
print("Froh_aut in percent for medium and long autosomal ROH")
print(paste((FROH_med_long_percent), "%"))
print("")

## Calculate FROH for long ROH only
chr_ROH_long <- chr_ROH %>% filter(chr_ROH$Length >= 1000000)
L_ROH_long <- sum(chr_ROH_long$Length)
FROH_long <- L_ROH_long/chrom_length
# Report as a percentage
FROH_long_percent <- FROH_long*100

print("Froh_aut value for long autosomal ROH is:")
print(FROH_long)
print(" ")
print("Froh_aut in percent for long autosomal ROH")
print(paste((FROH_long_percent), "%"))


## Calculate LROH -- the length of all the ROH on that chromosome
## Calculate Laut -- the length of the chromosome from first to last variant
## Calculate FROH for 1Mb, 500kb, and 100kb --> LROH/Laut
## Normalize FROH for chromosome length
## Report -- save in a .txt file









