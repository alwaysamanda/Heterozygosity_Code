#### Script to plot ROH on chromosomes from ROH Durbin Calculations ####
## Date: 20250620 (June 20th, 2025)
## Author: Amanda Gardiner
## Version 2 (Version 1 is 20250106_Plot_ROH.R)
## GOAL: Plot ROH on each chromosome in the genome of a species
## NOTES: Based on rds/users/ag2427/hpc-work/ROH_analyses/20241211_Plot_ROH_Durbin_Calc.R
## NOTES: Modifying script to include cases where no ROH are detected
####

#### ---- Load in necessary packages ---- ####
library(dplyr)
library(ggplot2)


#### ---- Load in data and variables ---- ####
args <- commandArgs()
Roh_file <- args[6]
chrom_length_file <- args[7]
clade <- args[8]
num_auto_chr <- as.numeric(args[9])
spec_name <- args[10]
num_all_chr <- as.numeric(args[11])


#### ---- Read ROH data and add column to categorize based on ROH size ---- ####
read_roh_data <- function(roh_data_file){
    if (file.size(roh_data_file) == 0){
        dat <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), 
            c('Chr', 'Start', 'End', 'Length', 'ROH_Type'))
    }
    else{
        dat <- read.csv(roh_data_file, header=FALSE) ## Pull from the ROH_results.csv file
        colnames(dat) <- c('Chr', 'Start', 'End', 'Length')
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
    }
    return(dat)
}

#### ---- Read in chromosomes and their lengths ---- ####
read_chrom_dat <- function(chrom_length_file, num_all_chr, num_auto_chr){
    Chrom <- read.table(chrom_length_file, header=FALSE)
    Chrom <- as.data.frame(Chrom[1:num_all_chr,])
    colnames(Chrom) <- c('Chrom_name', 'Chrom_end')
    return(Chrom)
}

#### ---- Plot results ---- ####
plot_ROH <- function(clade, spec_name, Roh_file, chrom_length_file, num_all_chr, num_auto_chr){
    dat <- read_roh_data(Roh_file)
    Chrom <- read_chrom_dat(chrom_length_file, num_all_chr, num_auto_chr)
    Chrom_autosomal <- Chrom[1:num_auto_chr, 1:2]
    autosomal_names <- seq(1, num_auto_chr, 1)
    chrom_numbers <- c(autosomal_names, rep("Sex", each=(num_all_chr - num_auto_chr)))

    chrom_numbers <- as.character(chrom_numbers)
    Chrom$Chrom_numbers <- chrom_numbers
    dat$Chrom_Label <- Chrom$Chrom_numbers[match(dat$Chr, Chrom$Chrom_name)]

    chrom_map <- setNames(seq_along(chrom_numbers), chrom_numbers)
    Chrom$Chrom_pos <- chrom_map[as.character(Chrom$Chrom_numbers)]
    dat$Chrom_pos <- chrom_map[as.character(dat$Chrom_Label)]

    file_name <- paste(clade, "/", spec_name, "/", spec_name, "_ROH_Map.pdf", sep = "")
     ## Plot base chromosomes
     genome_base = ggplot() +
     geom_bar(data = Chrom, 
        aes(x = Chrom_pos, y = Chrom_end) , 
        stat='identity', 
        fill='grey80', 
        colour='grey80', 
        width=.2) +
    labs(y = "Base pair position", 
        x = "Chromosome") +
    theme_minimal() +
    theme(panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            axis.text.x = element_text(angle = 90, hjust = 1), 
            axis.ticks.x = element_blank()
    ) + 
    scale_x_continuous(breaks = chrom_map, labels = names(chrom_map))

    ## Plot ROH
    seg_example = genome_base + 
            geom_segment(data = dat, 
                    aes(x=Chrom_pos, xend=Chrom_pos, y=Start, yend=End, color = ROH_Type), 
                    linewidth = 1.5
            ) + 
            scale_color_manual(values=c('#000066', '#0033FF', '#33CCFF')) + 
            guides(color = guide_legend(title = "ROH Length"))

    pdf(file = file_name, onefile=FALSE)
    print(seg_example)
    dev.off()
}

plot_ROH(clade, spec_name, Roh_file, chrom_length_file, num_all_chr, num_auto_chr)
