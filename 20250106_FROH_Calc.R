#### Script to calculate FROH (Inbreeding Coefficient from ROH Durbin Calculations ####
## Date: 20250106 (January 6th, 2025)
## Author: Amanda Gardiner
## Version 1
## NOTES: Based on rds/users/ag2427/hpc-work/ROH_analyses/20241211_FROH_ROH_Durbin_Calc.R

####
library(dplyr)
args <- commandArgs()

num_aut_chr <- args[11]

Species <- args[8]
clade <- args[9]
# vcf_name <- paste(clade, "/", Ref_name, "/", Ref_name, '_aligned.mm2.vcf', sep='')
# vcf <- read.table(vcf_name, header=TRUE)

# Load in file with chromosomes and their lengths
Chrom <- read.table(args[7], header=FALSE)
Chrom <- as.data.frame(Chrom)
colnames(Chrom) <- c('Chrom_name', 'Chrom_end')
Chrom_autosomal <- Chrom[1:num_aut_chr, 1:2]
Chrom_sum <- sum(Chrom_autosomal$Chrom_end)
# num_chr <- nrow(Chrom)

## Steps
## Load in csv file containing ROH
ROH_dat <- read.csv(args[6], header=FALSE)
colnames(ROH_dat) <- c('Chr', 'Start', 'End', 'Length')
L_ROH <- sum(ROH_dat$Length)

## Calculate traditional FROH (Froh = LROH/Laut)
FROH_trad_all <- L_ROH/Chrom_sum

print("The traditional Froh value for all ROH (short, med, long) is:")
print(FROH_trad_all)
print(" ")
print("Traditional Froh in percent for all ROH (short, med, long)")
print(paste((FROH_trad_all*100), "%"))
print(" ")

## Calculate FROHaut (Froh, aut = LROH/Laut, where Laut = length of genome from first to last variant detected on each chromosome summed over all chr)
Laut <- as.numeric(args[10])
FROH_aut_all <- L_ROH/Laut

print("Froh_aut value for all ROH (short, med, long) is:")
print(FROH_aut_all)
print(" ")
print("Froh_aut in percent for all ROH (short, med, long)")
print(paste((FROH_aut_all*100), "%"))

## Solely autosomal FROH from here on out
## FROHaut autosomal
Laut_autosomal <- as.numeric(args[12])
FROH_aut_autosomal_all <- L_ROH/Laut_autosomal

print("Froh_aut value for all autosomal ROH (short, med, long) is:")
print(FROH_aut_autosomal_all)
print(" ")
print("Froh_aut in percent for all autosomal ROH (short, med, long)")
print(paste((FROH_aut_autosomal_all*100), "%"))

## Calculate FROHaut for medium length ROH
ROH_dat_med <- ROH_dat %>% filter(ROH_dat$Length >= 500000)
L_ROH_med <- sum(ROH_dat_med$Length)

FROH_aut_med <- L_ROH_med/Laut_autosomal

print("Froh_aut value for medium and long autsomal ROH is:")
print(FROH_aut_med)
print(" ")
print("Froh_aut in percent for medium and long autosomal ROH is")
print(paste((FROH_aut_med*100), "%"))

## Calculate FROHaut for only long ROH
ROH_dat_long <- ROH_dat %>% filter(ROH_dat$Length >= 1000000)
L_ROH_long <- sum(ROH_dat_long$Length)

FROH_aut_long <- L_ROH_long/Laut_autosomal

print("Froh_aut value for long autosomal ROH is:")
print(FROH_aut_long)
print(" ")
print("Froh_aut in percent for long autosomal ROH is")
print(paste((FROH_aut_long*100), "%"))
