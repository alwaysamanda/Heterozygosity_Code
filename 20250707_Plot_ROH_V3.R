#### Script to plot ROH on chromosomes from ROH Durbin Calculations ####
## Date: 20250620 (June 20th, 2025)
## Author: Amanda Gardiner
## Version 3 (Version 1 is 20250620_Plot_ROH.R)
## GOAL: Plot ROH on each chromosome in the genome of a species
## NOTES: Based on rds/users/ag2427/hpc-work/ROH_analyses/20241211_Plot_ROH_Durbin_Calc.R
## NOTES: Modifying script to include cases where no ROH are detected
## NOTES: Modifying script to plot alingments of each chromosome next to them on the map
####

#### ---- Load in necessary packages ---- ####
library(dplyr)
library(ggplot2)

#### ---- Load in data and variables ---- ####
args <- commandArgs()
Roh_file <- args[6]
chrom_length_file <- args[7]
aln_file <- args[8]
clade <- args[9]
num_auto_chr <- as.numeric(args[10])
spec_name <- args[11]
num_all_chr <- as.numeric(args[12])
output_file <- args[13]

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

#### ---- Read in alignment file ---- ####
read_aln_dat <- function(aln_file_name){
    aln_dat <- read.table(aln_file_name, header=FALSE)
    dat <- as.data.frame(aln_dat[, 2:4])
    colnames(dat) <- c("Chrom", "Start", "End")
    return(dat)
}

#### ---- Plot results ---- ####
plot_ROH <- function(clade, spec_name, Roh_file, chrom_length_file, aln_file_name, num_all_chr, num_auto_chr, file_name){
    dat <- read_roh_data(Roh_file)
    Chrom <- read_chrom_dat(chrom_length_file, num_all_chr, num_auto_chr)
    aln_dat <- read_aln_dat(aln_file_name)

    Chrom_autosomal <- Chrom[1:num_auto_chr, 1:2]
    autosomal_names <- seq(1, num_auto_chr, 1)
    chrom_numbers <- c(autosomal_names, rep("Sex", each=(num_all_chr - num_auto_chr)))

    chrom_numbers <- as.character(chrom_numbers)
    Chrom$Chrom_numbers <- chrom_numbers
    dat$Chrom_Label <- Chrom$Chrom_numbers[match(dat$Chr, Chrom$Chrom_name)]

    aln_dat$Chrom_Label <- Chrom$Chrom_numbers[match(aln_dat$Chrom, Chrom$Chrom_name)]

    chrom_map <- setNames(seq_along(chrom_numbers), chrom_numbers)
    Chrom$Chrom_pos <- chrom_map[as.character(Chrom$Chrom_numbers)]
    dat$Chrom_pos <- chrom_map[as.character(dat$Chrom_Label)]
    aln_dat$Chrom_pos <- chrom_map[as.character(aln_dat$Chrom_Label)]

     ## Plot base chromosomes
     genome_base = ggplot() +
     geom_bar(data = Chrom, 
        aes(x = Chrom_pos, y = Chrom_end) , 
        stat='identity', 
        fill='grey80', 
        colour='grey80', 
        width=.4) +
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


    ## Plot roh and alingment coverage
    aln_plot = genome_base + 
    geom_segment(data = aln_dat,
        aes(x = Chrom_pos+0.1, xend = Chrom_pos+0.1, 
        y = Start, yend = End), 
        color = "#DB5461", 
        linewidth = 1
    )

    roh_seg_plot = aln_plot + 
            geom_segment(data = dat, 
                    aes(x=Chrom_pos-0.1, xend=Chrom_pos-0.1, y=Start, yend=End, color = ROH_Type), 
                    linewidth = 1.5
            ) + 
            scale_color_manual(values=c(Long = '#000066', 
            Medium = '#0033FF', 
            Short = '#33CCFF')) + 
            guides(color = guide_legend(title = "ROH Length"))

    png(file = file_name, width = 2200, height = 1600, res = 300)
    print(roh_seg_plot)
    dev.off()
}

plot_ROH(clade, spec_name, Roh_file, chrom_length_file, aln_file, num_all_chr, num_auto_chr, output_file)


