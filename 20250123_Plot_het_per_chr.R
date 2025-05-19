#### Script to plot heterozygosity ####
## Date: 20250123 (January 23rd, 2025)
## Author: Amanda Gardiner
## Version 1
## GOAL: To plot heterozygosity calculated by 20250114_find_het_per_chr.sh 
##       I want to plot heterozygosity for each individual chromosome, save it, and then combine them for a whole genome graph
## NOTES: Want plots to follow those seen in Stanhope et al. 2023 Fig. 2
##        In Fig. 2A: Chromosomes on x axis across whole genome, and heterozygosity on y-axis
##        Het for each window in a chromosome will be plotted
##        In Fig. 2B: Histogram of windows across whole genome -- counting the number of windows with each value of heterozygosity

####

## Load in necessary packages
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(scales)
library(patchwork)
library(svglite)

## Load in data
args <- commandArgs()
dat <- read.table(args[6], header=FALSE) ## dat will be a table of all het for windows across each chromosome, compiled into a single file
colnames(dat) <- c('Chr', 'Start', 'End', 'Het', 'Het_excl_ROH', 'Window_Size', 'Window_Size_excl_ROH', 'Het_Per_KB', 'Het_Per_KB_excl_ROH')
dat$Midpoint <- as.numeric(dat$Start) + 500000

date <- args[7]
ref_name <- args[8]
clade <- args[9]
num_chr <- args[10]
spec_name <- args[11]

## Get all chromosomes
all_chr <- unique(dat$Chr)

plot_list = list()

## Create loop to plot heterozygosity across each of the autosomal chromosomes
for (i in 1:num_chr) {
        single_chr <- all_chr[i] ## Select Chromosome 
        single_chr_dat <- dat %>% filter(dat$Chr == single_chr) ## Filter data for only heterozygosity on that chromosome

        single_chr_plot_file_name <- paste(clade, "/", spec_name, "/", date, "_", spec_name, "_" , single_chr, "_Het_Map.svg", sep = "") 
        if (i == 1) {
                single_chr_map = ggplot() + 
                        geom_area(data = single_chr_dat, 
                        aes(x = Midpoint, y = Het_Per_KB_excl_ROH), 
                        fill="#003C66"
                        ) + 
                        labs(y = "Heterozygosity per Kb") +
                        ylim(0,6) + 
                        theme_minimal() +
                        theme(panel.border = element_blank(),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),  
                        axis.ticks.x = element_blank(), 
                        axis.text.x = element_text(angle = 90),
                        axis.title.x = element_blank(), 
                        plot.margin = margin(t = 0.5, r = 0.01, l = 0.25, b = 0.5) 
                        ) + 
                        scale_x_continuous(labels = label_scientific(digits=2))
        }
        else {
                if ((i %% 2) == 0) {
                        fill_color = "#707070"
                }
                if ((i %% 2) != 0) {
                        fill_color = "#003C66"
                }
                single_chr_map = ggplot() + 
                        geom_area(data = single_chr_dat, 
                        aes(x = Midpoint, y = Het_Per_KB_excl_ROH), 
                        fill=fill_color
                        ) + 
                        ylim(0,6) + 
                        theme_minimal() +
                        theme(panel.border = element_blank(),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),  
                        axis.ticks.x = element_blank(), 
                        axis.text.y = element_blank(),
                        axis.text.x = element_text(angle = 90), 
                        axis.title = element_blank(), 
                        plot.margin = margin(t = 0.5, r = 0.01, l = 0.01, b = 0.5) 
                        ) + 
                        scale_x_continuous(labels = label_scientific(digits=2))
        }
        
        ggsave(single_chr_plot_file_name, single_chr_map, device = 'svg')

        # plot_list[[i]] <- single_chr_map

}

#### NOTE: COMMENTED OUT ON 20250325 TO MAKE CODE COMPATIBLE WITH SNAKEMAKE ####
####       COPIED INTO SISTER SCRIPT 20250325_Plot_het_whole_genome.R to get whole genome plot ####
# ## Plot all the chromosomes together using a single command
# num_plots <- length(plot_list)
# # print(num_plots)

# all_plots <- wrap_plots(plot_list, nrow=1)

# all_plot_file_name <- paste(clade, "/", spec_name, "/", date, "_", ref_name, "_Whole_Genome" , "_Het_Map.svg", sep = "")

# ggsave(all_plot_file_name, all_plots, width = num_plots, units = 'in', device = 'svg', limitsize = FALSE)



########## TEST BELOW ##########
## Plot all the chromosomes together using a single ggplot command
# test <- ggplot(dat, aes(Midpoint, Het_Per_KB, color = Chr)) +
#                 geom_area() + 
#                 facet_wrap(~Chr) + 
#                 labs(y = "Heterozygosity per Kb", 
#                 x = "Chromosome") +
#                 ylim(0,6) + 
#                 theme_minimal() +
#                 theme(panel.border = element_blank(),
#                 panel.grid.major = element_blank(),
#                 panel.grid.minor = element_blank(),  
#                 axis.ticks.x = element_blank(), 
#                 legend.position="none", 
#                 axis.text.x = element_blank(), 
#                 plot.title = element_blank(), 
#                 panel.spacing.x = unit(0.5, "lines"), 
#                 strip.text.x = element_blank()
#                 ) + 
#                 scale_x_continuous(labels = label_scientific(digits=2))

# ggsave(paste(clade, "/", spec_name, "/", date, "_", ref_name, '_test_group.pdf', sep=""), test, device = 'pdf', width = 21)


# ## Create plot to show heterozygosity on a specific chromosome
# spec_chrom_dat <- dat %>% filter(dat$Chr == 'NC_083444.1')
# spec_chrom_dat_2 <- dat %>% filter(dat$Chr == 'NC_083401.1')

# spec_chrom_map_file_name <- paste(clade, "/", spec_name, "/", date, "_", ref_name, "_Specific_Chr_Genome_Map.pdf", sep = "")
# pdf(file = spec_chrom_map_file_name)

# spec_chrom_map = ggplot() + 
#         geom_area(data = spec_chrom_dat, 
#                 aes(x = Midpoint, y = Het_Per_KB), 
#                 fill="#003C66"
#         ) + 
#         labs(y = "Heterozygosity per Kb", 
#         x = "Chromosome") +
#         ylim(0,6) + 
#         theme_minimal() +
#         theme(panel.border = element_blank(),
#             panel.grid.major = element_blank(),
#             panel.grid.minor = element_blank(),  
#             axis.ticks.x = element_blank(), 
#             legend.position = "none") + 
#         scale_x_continuous(labels = label_scientific(digits=2))

# # spec_chrom_map

# spec_chrom_map_2 = ggplot() + 
#         geom_area(data = spec_chrom_dat_2, 
#                 aes(x = Midpoint, y = Het_Per_KB), 
#                 fill="#3B4454"
#         ) + 
#         labs(y = "Heterozygosity per Kb", 
#         x = "Chromosome") +
#         ylim(0,6) + 
#         theme_minimal() +
#         theme(panel.border = element_blank(),
#             panel.grid.major = element_blank(),
#             panel.grid.minor = element_blank(), 
#             axis.text.y = element_blank(), 
#             axis.title.y = element_blank(), 
#             axis.ticks.x = element_blank(), 
#             legend.position = "none"
#             ) + 
#         scale_x_continuous(labels = label_scientific(digits=2))

# grid.arrange(spec_chrom_map, spec_chrom_map_2,ncol=2)

# dev.off()





# ## Create plot showing heterozygosity along genome by window and chromosome
# het_map_file_name <- paste(clade, "/", spec_name, "/", date, "_", ref_name, "_Het_Genome_Map.pdf", sep = "")
# pdf(file = het_map_file_name)

# het_map = ggplot() +
#     geom_bar(data = dat, 
#             aes(x = Chr, y = Het_Per_KB, fill=factor(Start)), 
#             position='dodge', 
#             stat='identity', 
#             width=1) + 
#     labs(y = "Heterozygosity per Kb", 
#             x = "Chromosome") +
#     theme_minimal() +
#     theme_light() +
#     theme(panel.border = element_blank(),
#             panel.grid.major = element_blank(),
#             panel.grid.minor = element_blank()) + 
#     theme(axis.text.x = element_text(angle = 90, hjust = 1),
#         axis.ticks.x = element_blank(), 
#         legend.position = "none")

# het_map

# dev.off()
