#### Script to plot ROH on chromosomes from ROH Durbin Calculations ####
## Date: 20250106 (January 6th, 2025)
## Author: Amanda Gardiner
## Version 1
## NOTES: Based on rds/users/ag2427/hpc-work/ROH_analyses/20241211_Plot_ROH_Durbin_Calc.R

####

## Load in necessary packages
library(dplyr)
library(ggplot2)

## Load in data
args <- commandArgs()
dat <- read.csv(args[6], header=FALSE) ## Pull from the ROH_results.csv file
colnames(dat) <- c('Chr', 'Start', 'End', 'Length')

date <- args[8]
ref_name <- args[9]
clade <- args[10]
num_auto_chr <- as.numeric(args[11])
spec_name <- args[12]
num_all_chr <- as.numeric(args[13])

## Add in column to data to categorize based on size of ROH
ROH_type <- c()

for (i in 1:nrow(dat)) { 
        if (dat[i,4] >= 1000000) { 
                ROH_type[i] <- 'Long' 
        } else if (dat[i,4] < 1000000 & dat[i,4] >= 500000) {
                ROH_type[i] <- 'Medium'
        } else {
                ROH_type[i] <- 'Short'
        }
}
dat$ROH_Type <- ROH_type

# Load in file with chromosomes and their lengths
Chrom <- read.table(args[7], header=FALSE)
Chrom <- as.data.frame(Chrom)
colnames(Chrom) <- c('Chrom_name', 'Chrom_end')
Chrom_autosomal <- Chrom[1:num_auto_chr, 1:2]
autosomal_names <- seq(1, num_auto_chr, 1)
chrom_numbers <- c(autosomal_names, rep("Sex", each=(num_all_chr - num_auto_chr)))

## Plot results
file_name <- paste(clade, "/", spec_name, "/", date, "_", spec_name, "_ROH_Map.pdf", sep = "")
# pdf(file = file_name)
# svg(file = "out.svg")

plot_ROH <- function(){
        ## Plot base chromosomes
        genome_base = ggplot() +
        geom_bar(data = Chrom, 
                aes(x = chrom_numbers, y = Chrom_end) , 
                stat='identity', 
                fill='grey80', 
                colour='grey80', 
                width=.2) +
        labs(y = "Base pair position", 
                x = "Chromosome") +
        theme_minimal() +
        theme(panel.border = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1),
                axis.ticks.x = element_blank())

        ## Plot ROH
        seg_example = genome_base + 
                geom_segment(data = dat, 
                        aes(x=Chr, xend=Chr, y=Start, yend=End, color = ROH_Type), 
                        linewidth = 1.5
                ) + 
                scale_color_manual(values=c('#000066', '#0033FF', '#33CCFF')) + 
                guides(color = guide_legend(title = "ROH Length"))

        pdf(file = file_name, onefile=FALSE)
        print(seg_example)
        dev.off()
}

plot_ROH()

# pdf(file = file_name, width = 8, height = 6)
# print(seg_example)
# # plot <- plot_func()

# dev.off()

