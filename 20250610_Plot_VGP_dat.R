#### Script to take tsv file containing all data and plot it for introductory figure ####
## Date: 20250610 (June 10th, 2025)
## Author: Amanda Gardiner
## Version 1
## GOAL: Create stacked bar charts to plot data for all species
####

###############################################################################

#### ---- Import necessary packages ---- ####
library(dplyr)
library(ggplot2)

#### ---- Load in data ---- ####
prim_dat_file <- "VGP_data/primary.taxon.metaData.tsv"
sec_dat_file <- "VGP_data/secondary.taxon.metaData.tsv"
all_dat_file <- "VGP_data/VGPPhase1-freeze-1.0.tsv"

## Read files
prim_dat <- read.csv(prim_dat_file, header=TRUE, sep='\t')
sec_dat <- read.csv(sec_dat_file, header=TRUE, sep='\t')
all_dat <- read.csv(all_dat_file, header=TRUE, sep='\t')

#### ---- Prep data for finding hap1/2 to align in analysis ---- ####
## Create an empty dataframe to store results
hap_df <- data.frame(
    taxId = c(),
    vgpClade = c(), 
    Scientific_name = c(), 
    Prim_accession = c(),
    Sec_accession = c(),
    NCBI_taxon_string = c(),
    stringsAsFactors=FALSE
)

## Use loop to find all species with a single seconary haplotype
for (row in 1:nrow(prim_dat)) {
    spec_name <- prim_dat[row,5] ## Get scientific name for species
    sec_index <- grep(spec_name, sec_dat$Scientific_name) ## Find secondary haplotype(s) for that species
    if (length(sec_index)== 0) {
        cat(paste0(spec_name, " has no secondary \n"))
    }
    else if (length(sec_index) > 1) {
        cat(paste0(spec_name, " has multiple secondary haplotypes \n"))
    }
    else {
        sec_row <- sec_index[1]
        sec_access <- sec_dat[sec_row, 1]
        newrow <- c(prim_dat[row,3], prim_dat[row,4], spec_name, prim_dat[row,1], sec_access, prim_dat[row,7])
        hap_df <- rbind(hap_df, newrow)
    }
}
colnames(hap_df) <- c("taxId", "vgpClade", "Scientific_name", "Prim_accession", "Sec_accession", "NCBI_taxon_string")

## Save to table so that I can access these to create more config files later
write.table(hap_df, "20250611_Prim_Sec_table.tsv", sep="\t")

#### ---- Prep data for plotting ---- ####
prim_df <- all_dat[all_dat$Accession...for.main.haplotype != "",] ## Get plotting data for primary assemblies
prim_sum_df <- prim_df %>% 
    group_by(Superorder, Lineage) %>% 
    summarise(sum = n(), .groups="drop")
colnames(prim_sum_df) = c("superorder", "lineage", "sum")

lineage_order <- prim_sum_df %>%
  group_by(lineage) %>%
  summarise(total = sum(sum), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(lineage)

prim_sum_df$lineage <- factor(prim_sum_df$lineage, levels = lineage_order)

superorder_order <- prim_sum_df %>%
  group_by(superorder) %>%
  summarise(total = sum(sum), .groups = "drop") %>%
  arrange(total) %>%
  pull(superorder)

prim_sum_df$superorder <- factor(prim_sum_df$superorder, levels = superorder_order)

#### ---- Plot Primary data only ---- ####
primary_output_file <- "test.svg"
svg(file = primary_output_file, height=60, width=100)

ggplot(prim_sum_df, aes(fill = superorder, y = sum, x = lineage)) +
  geom_bar(position = "stack", stat = "identity") +
  geom_text(
    aes(label = superorder),
    position = position_stack(vjust = 0.5),
    size = 20,  # Adjust this for legibility in your SVG units
    color = "black"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 80),
    axis.title = element_blank()
  ) +
  scale_fill_manual(values = c(
    "Actinopterygii (Ray-finned Fishes)" = "#d0d1e6", "Chondrichthyes (Cartilaginous Fishes)" = "#a6bddb", "Agnatha (Jawless Fishes)" = "#74a9cf", "Sarcopterygii (Lobe-finned fishes)" = "#3690c0",
    "Laurasiatheria" = "#e5f5e0", "Marsupials" = "#a1d99b", "Monotremes" = "#238b45", "Supraprimates (Euarchontoglires)" = "#c7e9c0", "Afrotheria" = "#74c476", "Xenarthra" = "#41ab5d", 
    "Squamata" = "#fdd0a2", "Testudines (Cryptodira)" = "#fdae6b", "Testudines (Pleurodira)" = "#fd8d3c", "Crocodylia" = "#f16913", 
    "Core landbirds" = "#fee0d2", "Core waterbirds" = "#fcbba1", "Columbea" = "#ef3b2c", "GalloAnserformes" = "#fc9272", "Otidimorphae" = "#a50f15", "Caprimulgimorphae" = "#fb6a4a", "Palaeognathae" = "#cb181d",
    "Anura" = "#bcbddc", "Caudata" = "#9e9ac8", "Gymnophiona" = "#807dba", 
    "Tunicata" = "#d9d9d9", "Echinodermata" = "#737373", "Cephalochordata" = "#bdbdbd", "Hemichordata" = "#969696"
  ))

dev.off()

