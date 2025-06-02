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

ref_name <- args[8]
clade <- args[9]
num_auto_chr <- as.numeric(args[10])
spec_name <- args[11]
num_all_chr <- as.numeric(args[12])

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
Chrom <- as.data.frame(Chrom[1:num_all_chr,])
colnames(Chrom) <- c('Chrom_name', 'Chrom_end')
Chrom_autosomal <- Chrom[1:num_auto_chr, 1:2]
autosomal_names <- seq(1, num_auto_chr, 1)
chrom_numbers <- c(autosomal_names, rep("Sex", each=(num_all_chr - num_auto_chr)))

chrom_numbers <- as.character(chrom_numbers)
Chrom$Chrom_numbers <- chrom_numbers
dat$Chrom_Label <- Chrom$Chrom_numbers[match(dat$Chr, Chrom$Chrom_name)]

chrom_map <- setNames(seq_along(chrom_numbers), chrom_numbers)
Chrom$Chrom_pos <- chrom_map[as.character(Chrom$Chrom_numbers)]
dat$Chrom_pos <- chrom_map[as.character(dat$Chrom_Label)]

## Plot results
file_name <- paste(clade, "/", spec_name, "/", spec_name, "_ROH_Map.pdf", sep = "")
# pdf(file = file_name)
# svg(file = "out.svg")


plot_ROH <- function(){
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

plot_ROH()

# pdf(file = file_name, width = 8, height = 6)
# print(seg_example)
# # plot <- plot_func()

# dev.off()

