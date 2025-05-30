#### SNAKEFILE FOR AUTOMATED PIPLELINE OF RUNNING PIPELINE ####

"""
login to icelake node using login-icelake.hpc.cam.ac.uk
start a screen terminal (screen -RD)
make sure to use snakemake version 7.8.5

check how many jobs are required AT EACH STAGE, this command also shows the shell commands it would use:
- snakemake -n --printshellcmds -p

to run locally use:
    snakemake --configfile config.yml
to run on slurm use:
    snakemake --executor slurm --jobs 30 --workflow-profile profiles --latency-wait 60 --keep-going --configfile config.yml
For changing input parameters, change in the config.yml file

Generate simplified plot to show rule flow with:
    snakemake --executor slurm --workflow-profile profiles --jobs 10 --configfile config.yml --rulegraph | dot -Tpdf > dag.pdf
Show full complicated path for each chromosome with:
    snakemake --executor slurm --workflow-profile profiles --jobs 10 --configfile config.yml --dag | dot -Tpdf > dag_complex.pdf

NOTES:
    Based on Snakefile made by Bettina Fisher (https://github.com/bef22/vcf_snakemake/blob/main/Snakefile)

AUTHORS and UPDATES:
Amanda Gardiner, 20250313

"""

###############################################################################

user = "ag2427"

### set global parameters
CLADE=config['CLADE']
SPEC_NAME=config['SPEC_NAME']
REF_NAME=config['REF_NAME']
ALT_NAME=config['ALT_NAME']
TODAY_DATE=config['TODAY_DATE']
CHROM_START_CHR=config['CHROM_START_CHR']
CHROM_LIST_FILE=config['CHROM_LIST_FILE']
WINDOW_INTERVAL=config['WINDOW_INTERVAL']
WINDOW_LENGTH=config['WINDOW_LENGTH']
NUM_AUT_CHROMOSOMES=config['NUM_AUT_CHROMOSOMES']
NUM_ALL_CHR=config['NUM_ALL_CHR']
CHROMS=config['ALL_CHROMOSOMES']
AUTO_CHROMS=config['AUTOSOMAL_CHROMOSOMES']
BOOTSTRAPS=config['BOOTSTRAPPING_VALUES']
MU=config['MUTATION_RATE']
GEN_TIME=config['GENERATION_TIME']

### paths to dependent scripts
ROH_CALC_PY = "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/20250320_ROH_Durbin_Calc_Eqns_V6.py"
ROH_PLOT_R = "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/20250106_Plot_ROH.R"
FROH_CALC_PER_CHR_R = "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/20250210_FROH_per_chr_calc.R"
FROH_CALC_R = "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/20250106_FROH_Calc.R"
# CALC_HET_PER_CHR_PY = "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/20250325_find_het_per_chr_V3.py"
CALC_HET_PER_CHR_PY = "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/20250520_find_het_per_chr_V4.py"
CALC_HET_WHOLE_GENOME_PY = "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/20250325_find_het_whole_genome_V3.py"
PLOT_HET_PER_CHR_R = "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/20250123_Plot_het_per_chr.R"
PLOT_WHOLE_HET_R = "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/20250325_Plot_het_whole_genome.R"
PLOT_MSMC_R = "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/20250508_Plot_MSMC.R"
PLOT_MSMC_BOOSTRAP_R = "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/20250529_Plot_MSMC_Bootstrap.R"
ROH_MASKER="/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/20250528_ROH_Masker.py"
GENERATE_MULTIHETSEP="/rds/project/rds-8b3VcZwY7rY/users/ag2427/hpc-work/msmc-tools/generate_multihetsep.py" ##change path, the script can be found in "MSMC-tools"
MASK_FILE_GENERATOR="/rds/project/rds-8b3VcZwY7rY/users/ag2427/hpc-work/20241011_jaskaran_script_msmc_fileprep/mask_file.py" ##In this directory
BOOTSTRAPPING_GENERATOR="/rds/project/rds-8b3VcZwY7rY/users/ag2427/hpc-work/20241011_jaskaran_script_msmc_fileprep/block_bootstrap_mhs.py" ## From "Cobraa-Trevor Cousins" also in this github directory
BAM_CALLER="/rds/project/rds-8b3VcZwY7rY/users/ag2427/hpc-work/msmc-tools/bamCaller.py" ##change path, the script can be found in "MSMC-tools"
INDEL_MASKER="/rds/project/rds-8b3VcZwY7rY/users/ag2427/hpc-work/20241011_jaskaran_script_msmc_fileprep/indel_mask.py" ##In this directory

###############################################################################
### the following shouldn't be changed to ensure working software versions on icelake are used
bcfPath = f'/home/{user}/rds/rds-durbin-group-8b3VcZwY7rY/software_RHEL8/bin/bcftools'
vcfPath = f'/home/{user}/rds/rds-durbin-group-8b3VcZwY7rY/software_RHEL8/bin/vcftools'
tabixPath = f'/home/{user}/rds/rds-durbin-group-8b3VcZwY7rY/software_RHEL8/bin/tabix'

###############################################################################

import pdb
import subprocess
import os
import re

wildcard_constraints:
    CLADE=CLADE,
    SPEC_NAME=SPEC_NAME, 
    REF_NAME=REF_NAME,
    ALT_NAME=ALT_NAME, 
    TODAY_DATE=TODAY_DATE, 
    CHROM_START_CHR=CHROM_START_CHR,
    WINDOW_INTERVAL=WINDOW_INTERVAL, 
    WINDOW_LENGTH=WINDOW_LENGTH, 
    NUM_AUT_CHROMOSOMES=NUM_AUT_CHROMOSOMES, 
    NUM_ALL_CHR=NUM_ALL_CHR, 
    CHROMS=CHROMS, 
    AUTO_CHROMS=AUTO_CHROMS,
    BOOTSTRAPS=BOOTSTRAPS

rule all: 
    input:
        expand("{CLADE}/{SPEC_NAME}/{SPEC_NAME}_FASTGA.chain.chr.fltr.srt.paf", CLADE=CLADE, SPEC_NAME=SPEC_NAME), 
        expand("{SPEC_NAME}_ALN.chain.pdf", SPEC_NAME=SPEC_NAME), 
        expand("{CLADE}/chrom_lists/{REF_NAME}_chroms.txt", CLADE=CLADE, REF_NAME=REF_NAME), 
        expand("{CLADE}/{SPEC_NAME}/temp/{SPEC_NAME}_Var_Only.txt", CLADE=CLADE, SPEC_NAME=SPEC_NAME), 
        expand("{CLADE}/{SPEC_NAME}/temp/{SPEC_NAME}_Aln_Only.txt", CLADE=CLADE, SPEC_NAME=SPEC_NAME), 
        expand("{CLADE}/{SPEC_NAME}/temp/{TODAY_DATE}_{CHR}_Aln_Only.txt", CLADE=CLADE, SPEC_NAME=SPEC_NAME, TODAY_DATE=TODAY_DATE, CHR=CHROMS), 
        expand("{CLADE}/{SPEC_NAME}/temp/{TODAY_DATE}_{CHR}_Var_Only.txt", CLADE=CLADE, SPEC_NAME=SPEC_NAME, TODAY_DATE=TODAY_DATE, CHR=CHROMS), 
        expand("{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{CHR}_ROH_Results.txt", CLADE=CLADE, SPEC_NAME=SPEC_NAME, TODAY_DATE=TODAY_DATE, CHR=CHROMS), 
        expand("{CLADE}/{SPEC_NAME}/{SPEC_NAME}_ROH_Results.csv", CLADE=CLADE, SPEC_NAME=SPEC_NAME), 
        expand("{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{SPEC_NAME}_ROH_Map.pdf", CLADE=CLADE, SPEC_NAME=SPEC_NAME, TODAY_DATE=TODAY_DATE), 
        expand("{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{CHR}_{REF_NAME}_FROH.txt", CHR=CHROMS, CLADE=CLADE, SPEC_NAME=SPEC_NAME, TODAY_DATE=TODAY_DATE, REF_NAME=REF_NAME), 
        expand("{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{SPEC_NAME}_FROH.txt", CLADE=CLADE, SPEC_NAME=SPEC_NAME, TODAY_DATE=TODAY_DATE), 
        expand("{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{CHROM}_het.txt", CHROM=AUTO_CHROMS, CLADE=CLADE, SPEC_NAME=SPEC_NAME, TODAY_DATE=TODAY_DATE, REF_NAME=REF_NAME), 
        expand("{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{SPEC_NAME}_per_chr_mean_heterozygosity.txt", CLADE=CLADE, SPEC_NAME=SPEC_NAME, TODAY_DATE=TODAY_DATE), 
        expand("{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{SPEC_NAME}_whole_genome_mean_heterozygosity.txt", CLADE=CLADE, SPEC_NAME=SPEC_NAME, TODAY_DATE=TODAY_DATE), 
        expand("{CLADE}/{SPEC_NAME}/{SPEC_NAME}_Het_Compiled.tsv", CLADE=CLADE, SPEC_NAME=SPEC_NAME), 
        expand("{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{SPEC_NAME}_{CHROM}_Het_Map.svg", CHROM=AUTO_CHROMS, CLADE=CLADE, SPEC_NAME=SPEC_NAME, TODAY_DATE=TODAY_DATE, REF_NAME=REF_NAME),
        expand("{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{SPEC_NAME}_Het_Whole_Genome_Map.svg", CLADE=CLADE, SPEC_NAME=SPEC_NAME, TODAY_DATE=TODAY_DATE), 
        expand("{CLADE}/{SPEC_NAME}/MSMC/paf_files/{SPEC_NAME}_{CHROM}.paf", CLADE=CLADE, SPEC_NAME=SPEC_NAME, CHROM=AUTO_CHROMS), 
        expand("{CLADE}/{SPEC_NAME}/MSMC/vcf_files/{SPEC_NAME}_{CHROM}.vcf.gz", CLADE=CLADE, SPEC_NAME=SPEC_NAME, CHROM=AUTO_CHROMS), 
        # expand("{CLADE}/{SPEC_NAME}/MSMC/Output_Mask/{SPEC_NAME}_{CHROM}.bed.gz", CLADE=CLADE, SPEC_NAME=SPEC_NAME, CHROM=AUTO_CHROMS), 
        # expand("{CLADE}/{SPEC_NAME}/MSMC/Output_VCF/{SPEC_NAME}_{CHROM}_only_snps.vcf.gz", CLADE=CLADE, SPEC_NAME=SPEC_NAME, CHROM=AUTO_CHROMS), 
        # expand("{CLADE}/{SPEC_NAME}/MSMC/ROH_Negative_Mask/{SPEC_NAME}_{CHROM}_ROH_Mask.bed.gz", CLADE=CLADE, SPEC_NAME=SPEC_NAME, CHROM=AUTO_CHROMS),
        expand("{CLADE}/{SPEC_NAME}/MSMC/Output_NEW_multihetsep/{CHROM}_multihet.txt", CLADE=CLADE, SPEC_NAME=SPEC_NAME, CHROM=AUTO_CHROMS), 
        # expand("{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{SPEC_NAME}_MSMC2_test.svg", CLADE=CLADE, SPEC_NAME=SPEC_NAME, TODAY_DATE=TODAY_DATE), 




        # expand("{CLADE}/{SPEC_NAME}/{SPEC_NAME}.sam", CLADE=CLADE, SPEC_NAME=SPEC_NAME), 
        # expand("{CLADE}/{SPEC_NAME}/{SPEC_NAME}.bam", CLADE=CLADE, SPEC_NAME=SPEC_NAME), 
        # expand("{CLADE}/{SPEC_NAME}/{SPEC_NAME}.bam.bai", CLADE=CLADE, SPEC_NAME=SPEC_NAME), 
        # expand("{CLADE}/{SPEC_NAME}/{SPEC_NAME}_{CHROM}.bam", CLADE=CLADE, SPEC_NAME=SPEC_NAME, CHROM=AUTO_CHROMS), 
        # expand("{CLADE}/{SPEC_NAME}/{SPEC_NAME}_{CHROM}.bam.bai", CLADE=CLADE, SPEC_NAME=SPEC_NAME, CHROM=AUTO_CHROMS), 
        # expand("{CLADE}/{SPEC_NAME}/temp/{REF_NAME}.fa", CLADE=CLADE, SPEC_NAME=SPEC_NAME, REF_NAME=REF_NAME), 
        # expand("{CLADE}/{SPEC_NAME}/temp/{CHROM}.fasta", CLADE=CLADE, SPEC_NAME=SPEC_NAME, CHROM=AUTO_CHROMS), 
        # expand("{CLADE}/{SPEC_NAME}/temp/{CHROM}_restructured.vcf.gz", CLADE=CLADE, SPEC_NAME=SPEC_NAME, CHROM=AUTO_CHROMS), 
        # expand("{CLADE}/{SPEC_NAME}/temp/{CHROM}_pileup.vcf.gz", CLADE=CLADE, SPEC_NAME=SPEC_NAME, CHROM=AUTO_CHROMS), 
        # expand("{CLADE}/{SPEC_NAME}/{CHROM}.recode.vcf", CLADE=CLADE, SPEC_NAME=SPEC_NAME, CHROM=AUTO_CHROMS), 
        # expand("{CLADE}/{SPEC_NAME}/{CHROM}.recode.vcf.gz", CLADE=CLADE, SPEC_NAME=SPEC_NAME, CHROM=AUTO_CHROMS), 
        # expand("{CLADE}/{SPEC_NAME}/temp/{CHROM}_only_snps.vcf.bgz", CLADE=CLADE, SPEC_NAME=SPEC_NAME, CHROM=AUTO_CHROMS), 
        # expand("{CLADE}/{SPEC_NAME}/temp/{CHROM}_only_snps.vcf.bgz.tbi", CLADE=CLADE, SPEC_NAME=SPEC_NAME, CHROM=AUTO_CHROMS), 
        # expand("{CLADE}/{SPEC_NAME}/MSMC/Output_indel_BED_file/{CHROM}_indel_bed.txt", CLADE=CLADE, SPEC_NAME=SPEC_NAME, CHROM=AUTO_CHROMS), 
        # expand("{CLADE}/{SPEC_NAME}/MSMC/Output_bam_caller_BED/{CHROM}_unmodified_bed.txt.gz", CLADE=CLADE, SPEC_NAME=SPEC_NAME, CHROM=AUTO_CHROMS),
        # expand("{CLADE}/{SPEC_NAME}/MSMC/Output_bam_caller_BED/{CHROM}_modified_bed.txt.gz", CLADE=CLADE, SPEC_NAME=SPEC_NAME, CHROM=AUTO_CHROMS),  
        # expand("{CLADE}/{SPEC_NAME}/MSMC/Output_multihet_filtered_95/{CHROM}_multihet.txt", CHROM=AUTO_CHROMS, CLADE=CLADE, SPEC_NAME=SPEC_NAME), 
        # expand("{CLADE}/{SPEC_NAME}/MSMC/Output_Multihet_Bootstrapped/{BOOT}/{CHROM}_multihet.txt", CHROM=AUTO_CHROMS, BOOT=BOOTSTRAPS, CLADE=CLADE, SPEC_NAME=SPEC_NAME),
        # expand("{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{SPEC_NAME}.msmc2.final.txt", CLADE=CLADE, SPEC_NAME=SPEC_NAME, TODAY_DATE=TODAY_DATE), 
        # expand("{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{SPEC_NAME}_Bootstrapping_{BOOT}.msmc2", CLADE=CLADE, SPEC_NAME=SPEC_NAME, TODAY_DATE=TODAY_DATE, BOOT=BOOTSTRAPS), 
        # expand("{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{SPEC_NAME}_MSMC2_test.svg", CLADE=CLADE, SPEC_NAME=SPEC_NAME, TODAY_DATE=TODAY_DATE), 
        # # # f"{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{SPEC_NAME}_BAM_Coverage.txt",  
        # f"{CLADE}/{SPEC_NAME}/{SPEC_NAME}_FASTGA_collinear_plot.pdf"




#### RULES FOR FASTGA ALINGMENT AND PAF FILE GENERATION ####
rule FAtoGDB_REF:
    input: 
        f"/rds/project/rds-p67MZilb2eQ/projects/VGP/241117.UCSC-hubs-VGP-alignment/alignment/reference/{CLADE}/{REF_NAME}.fa.gz"
    output:
        f"{CLADE}/{SPEC_NAME}/temp/{SPEC_NAME}_{REF_NAME}.1gdb",
    shell:
        """
        FAtoGDB -v {input} {output}
        """

rule FAtoGDB_ALT:
    input: 
        f"/rds/project/rds-p67MZilb2eQ/projects/VGP/241117.UCSC-hubs-VGP-alignment/alignment/alternate/{CLADE}/{ALT_NAME}.fa.gz"
    output:
        f"{CLADE}/{SPEC_NAME}/temp/{SPEC_NAME}_{ALT_NAME}.1gdb",
    shell:
        """
        FAtoGDB -v {input} {output}
        """

rule GIXmake_REF:
    input: 
        f"{CLADE}/{SPEC_NAME}/temp/{SPEC_NAME}_{REF_NAME}.1gdb"
    output:
        f"{CLADE}/{SPEC_NAME}/temp/{SPEC_NAME}_{REF_NAME}.gix"
    shell:
        """
        GIXmake -v -P. -T8 {input}
        """

rule GIXmake_ALT:
    input: 
        f"{CLADE}/{SPEC_NAME}/temp/{SPEC_NAME}_{ALT_NAME}.1gdb"
    output:
        f"{CLADE}/{SPEC_NAME}/temp/{SPEC_NAME}_{ALT_NAME}.gix"
    shell:
        """
        GIXmake -v -P. -T8 {input}
        """

rule FASTGA:
    input:
        REF=expand("{CLADE}/{SPEC_NAME}/temp/{SPEC_NAME}_{REF_NAME}.gix", CLADE=CLADE, SPEC_NAME=SPEC_NAME, REF_NAME=REF_NAME),
        ALT=expand("{CLADE}/{SPEC_NAME}/temp/{SPEC_NAME}_{ALT_NAME}.gix", CLADE=CLADE, SPEC_NAME=SPEC_NAME, ALT_NAME=ALT_NAME)
    output:
        f"{CLADE}/{SPEC_NAME}/{SPEC_NAME}_ALN.1aln"
    resources:
        runtime="12h"
    shell:
        """
        FastGA -v -P. -T8 -1:{output} {input.ALT} {input.REF}
        """
    
rule ALNCHAIN:
    input:
        f"{CLADE}/{SPEC_NAME}/{SPEC_NAME}_ALN.1aln"
    output:
        f"{CLADE}/{SPEC_NAME}/{SPEC_NAME}_ALN.chain.1aln"
    shell:
        """
        ALNchain {input}
        """

rule ALNtoPAF:
    input:
        f"{CLADE}/{SPEC_NAME}/{SPEC_NAME}_ALN.chain.1aln"
    output:
        f"{CLADE}/{SPEC_NAME}/{SPEC_NAME}_FASTGA.chain.paf"
    shell:
        """
        ALNtoPAF -s -T8 {input} > {output}
        """

rule FILTER_PAF_VARIANCE:
    input:
        f"{CLADE}/{SPEC_NAME}/{SPEC_NAME}_FASTGA.chain.paf"
    output:
        f"{CLADE}/{SPEC_NAME}/{SPEC_NAME}_FASTGA.chain.fltr.paf"
    shell:
        """
        awk -v OFS="\t" '{{if(substr($13,6)<=0.1) print}}' {input} > {output}
        """

rule SORT_PAF:
    input:
        f"{CLADE}/{SPEC_NAME}/{SPEC_NAME}_FASTGA.chain.fltr.paf"
    output:
        f"{CLADE}/{SPEC_NAME}/{SPEC_NAME}_FASTGA.chain.fltr.srt.paf"
    shell:
        """
        sort -k6,6V -k8,8n {input} > {output}
        """


rule FILTER_PAF_CHR_ONLY:
    input:
        PAF=f"{CLADE}/{SPEC_NAME}/{SPEC_NAME}_FASTGA.chain.fltr.srt.paf", 
        ALL_CHROMS=f"{CLADE}/chrom_lists/{REF_NAME}_chroms.txt"
    output:
        f"{CLADE}/{SPEC_NAME}/{SPEC_NAME}_FASTGA.chain.chr.fltr.srt.paf"
    shell:
        """
        awk 'BEGIN {{ while (getline < "{input.ALL_CHROMS}") list[$0] }} $6 in list' {input.PAF} > {output}
        """

rule ALNPLOT:
    input:
        expand("{CLADE}/{SPEC_NAME}/{SPEC_NAME}_ALN.chain.1aln", CLADE=CLADE, SPEC_NAME=SPEC_NAME)
    output:
        "{SPEC_NAME}_ALN.chain.pdf"
    shell:
        """
        ALNplot -p -H500 {input}
        """
####

#### RULES FOR GETTING NECESSARY FILES AND INFO  ####
rule GET_CHROM_LISTS:
    input:
    "../250430.VGP-Phase1/alignment/reference/{CLADE}/{REF_NAME}.fa.gz"
    output:
    "{CLADE}/chrom_lists/{REF_NAME}_chroms.txt"
    params:
        CHROM_START_CHR=CHROM_START_CHR
    shell:
        """
        zcat < {input} | grep '>{params.CHROM_START_CHR}' > {output}
        """

rule CHROM_LENGTH_CALC:
    input:
        FASTA=f"../250430.VGP-Phase1/alignment/reference/{CLADE}/{REF_NAME}.fa.gz", 
        ALL_CHROMS=f"{CLADE}/chrom_lists/{REF_NAME}_chroms.txt",
    output:
        f"{CLADE}/{SPEC_NAME}/Reference_{SPEC_NAME}_Chroms_Lengths.txt"
    params:
        CHROM_START_CHR=CHROM_START_CHR
    shell:
        """
        zcat {input.FASTA} \
        | awk '
            $0 ~ ">" {{if (NR > 1) {{print name "\t" len;}}
                name = substr($0, 2);
                len = 0;
                next;
            }}
            {{
                len += length($0);
            }}
            END {{
                print name "\t" len;
            }}
        ' \
        | grep -F -f {input.ALL_CHROMS} > {output}
        """


rule GET_WHOLE_VAR:
    input:
        f"{CLADE}/{SPEC_NAME}/{SPEC_NAME}_FASTGA.chain.chr.fltr.srt.paf"
    output:
        f"{CLADE}/{SPEC_NAME}/temp/{SPEC_NAME}_Var_Only.txt"
    shell:
        """
        k8 paftools.js call {input} | awk '$1 == "V" && $2 !~ /^JA/' > {output}
        """

rule GET_WHOLE_ALN:
    input:
        f"{CLADE}/{SPEC_NAME}/{SPEC_NAME}_FASTGA.chain.chr.fltr.srt.paf"
    output:
        f"{CLADE}/{SPEC_NAME}/temp/{SPEC_NAME}_Aln_Only.txt"
    shell:
        """
        k8 paftools.js call {input} | awk '$1 == "R" && $2 !~ /^JA/' > {output}
        """

rule GET_ALN_ONLY_PER_CHROM:
    input:
        f"{CLADE}/{SPEC_NAME}/temp/{SPEC_NAME}_Aln_Only.txt"
    output:
        chr_aln_only="{CLADE}/{SPEC_NAME}/temp/{TODAY_DATE}_{CHR}_Aln_Only.txt"
    shell:
        """
        awk '$1 == "R" && $2 == "{wildcards.CHR}"' {input} > {output.chr_aln_only}
        """

rule GET_VAR_ONLY_PER_CHROM:
    input:
        f"{CLADE}/{SPEC_NAME}/temp/{SPEC_NAME}_Var_Only.txt"
    output:
        chr_var_only="{CLADE}/{SPEC_NAME}/temp/{TODAY_DATE}_{CHR}_Var_Only.txt"
    shell:
        """
        awk '$1 == "V" && $2 == "{wildcards.CHR}"' {input} > {output.chr_var_only}
        """


####

#### RULES FOR CALCULATING ROH ####
rule ROH_CALC: # Calculations to find ROH
    input:
        chr_aln_only="{CLADE}/{SPEC_NAME}/temp/{TODAY_DATE}_{CHR}_Aln_Only.txt", 
        chr_var_only="{CLADE}/{SPEC_NAME}/temp/{TODAY_DATE}_{CHR}_Var_Only.txt", 
        CHROM_LENGTH_FILE=f"{CLADE}/{SPEC_NAME}/Reference_{SPEC_NAME}_Chroms_Lengths.txt"
    output:
        ROH_outfiles="{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{CHR}_ROH_Results.txt"
    params: 
        TODAY_DATE=TODAY_DATE, 
        REF_NAME=REF_NAME, 
        CLADE=CLADE, 
        SPEC_NAME=SPEC_NAME
    shell:
        """
        CHROM_LENGTH="$(awk -v var='{wildcards.CHR}' '$1 == var {{print $2; exit}}' {input.CHROM_LENGTH_FILE})"
        python {ROH_CALC_PY} {wildcards.CHR} $CHROM_LENGTH {params.TODAY_DATE} {params.REF_NAME} {params.CLADE} {input.chr_aln_only} {input.chr_var_only} {params.SPEC_NAME}
        """

rule COMPILE_ROH:
    input:
        expand("{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{CHR}_ROH_Results.txt", CLADE=CLADE, SPEC_NAME=SPEC_NAME, TODAY_DATE=TODAY_DATE, CHR=CHROMS)
    output:
        f"{CLADE}/{SPEC_NAME}/{SPEC_NAME}_ROH_Results.csv"
    params:
        CLADE=CLADE, 
        SPEC_NAME=SPEC_NAME, 
        TODAY_DATE=TODAY_DATE
    shell:
        """
        awk -v var='chrom' -F',' '{{if ($1!=var) {{ print $0 }}}}' {input} >> {output}
        """

####

#### RULES FOR PLOTTING ROH ####
rule Plot_ROH:
    input:
        ROH="{CLADE}/{SPEC_NAME}/{SPEC_NAME}_ROH_Results.csv",
        CHROM_LENGTH_FILE="{CLADE}/{SPEC_NAME}/Reference_{SPEC_NAME}_Chroms_Lengths.txt"
    output:
        "{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{SPEC_NAME}_ROH_Map.pdf"
    params:
        TODAY_DATE=TODAY_DATE,
        REF_NAME=REF_NAME,
        CLADE=CLADE,
        NUM_AUT_CHROMOSOMES=NUM_AUT_CHROMOSOMES, 
        SPEC_NAME=SPEC_NAME, 
        NUM_ALL_CHR=NUM_ALL_CHR
    shell:
        """
        Rscript {ROH_PLOT_R} {input.ROH} {input.CHROM_LENGTH_FILE} {params.TODAY_DATE} {params.REF_NAME} {params.CLADE} {params.NUM_AUT_CHROMOSOMES} {params.SPEC_NAME} {params.NUM_ALL_CHR} > {output}
        """
####

#### RULES FOR CALCULATING FROH ####
rule FROH_PER_AUT_CHR:
    input:
        ROH="{CLADE}/{SPEC_NAME}/{SPEC_NAME}_ROH_Results.csv", 
        CHROM_LENGTH_FILE="{CLADE}/{SPEC_NAME}/Reference_{SPEC_NAME}_Chroms_Lengths.txt", 
        ALL_CHROMS="{CLADE}/chrom_lists/{REF_NAME}_chroms.txt",
        VAR_FILE="{CLADE}/{SPEC_NAME}/temp/{TODAY_DATE}_{CHR}_Var_Only.txt"
    output:
        outfiles="{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{CHR}_{REF_NAME}_FROH.txt"
    params:
        TODAY_DATE=TODAY_DATE,
        REF_NAME=REF_NAME,
        CLADE=CLADE,
        NUM_AUT_CHROMOSOMES=NUM_AUT_CHROMOSOMES, 
        SPEC_NAME=SPEC_NAME
    shell:
        """
        start="$(head -n  1 {input.VAR_FILE} | awk ' {{print $3}} ')"
        end="$(tail -n  1 {input.VAR_FILE} | awk ' {{print $4}} ')"
        chrom_length="$((end - start))"

        Rscript {FROH_CALC_PER_CHR_R} {input.ROH} {input.CHROM_LENGTH_FILE} {params.TODAY_DATE} {params.REF_NAME} {params.CLADE} {params.NUM_AUT_CHROMOSOMES} "$chrom_length" {wildcards.CHR} > {output.outfiles}
        """

rule WHOLE_FROH:
    input:
        ROH="{CLADE}/{SPEC_NAME}/{SPEC_NAME}_ROH_Results.csv", 
        CHROM_LENGTH_FILE="{CLADE}/{SPEC_NAME}/Reference_{SPEC_NAME}_Chroms_Lengths.txt", 
        ALL_CHROMS=expand("{CLADE}/chrom_lists/{REF_NAME}_chroms.txt", CLADE=CLADE, REF_NAME=REF_NAME), 
        VAR_FILE="{CLADE}/{SPEC_NAME}/temp/{SPEC_NAME}_Var_Only.txt"
    output:
        "{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{SPEC_NAME}_FROH.txt"
    params:
        TODAY_DATE=TODAY_DATE,
        REF_NAME=REF_NAME,
        CLADE=CLADE,
        NUM_AUT_CHROMOSOMES=NUM_AUT_CHROMOSOMES, 
        SPEC_NAME=SPEC_NAME
    shell:
        """
        set -e 
        
        Laut_autosomal=0
        Laut=0

        while read -r line; do
            start="$(head -n 1 {params.CLADE}/{params.SPEC_NAME}/temp/{params.TODAY_DATE}_"$line"_Var_Only.txt | awk ' {{print $3}} ')"
            end="$(tail -n 1 {params.CLADE}/{params.SPEC_NAME}/temp/{params.TODAY_DATE}_"$line"_Var_Only.txt | awk ' {{print $4}} ')"
            chrom_length="$((end - start))"
            Laut_autosomal=$(("$Laut_autosomal" + "$chrom_length"))
        done < <(head -n {params.NUM_AUT_CHROMOSOMES} "{input.ALL_CHROMS}") 

        while read -r line; do
            start="$(head -n 1 {params.CLADE}/{params.SPEC_NAME}/temp/{params.TODAY_DATE}_"$line"_Var_Only.txt | awk ' {{print $3}} ')"
            end="$(tail -n 1 {params.CLADE}/{params.SPEC_NAME}/temp/{params.TODAY_DATE}_"$line"_Var_Only.txt | awk ' {{print $4}} ')"
            chrom_length="$((end - start))"
            Laut=$(("$Laut" + "$chrom_length"))
        done < "{input.ALL_CHROMS}"

        Rscript {FROH_CALC_R} {input.ROH} {input.CHROM_LENGTH_FILE} {params.TODAY_DATE} {params.REF_NAME} {params.CLADE} "$Laut" {params.NUM_AUT_CHROMOSOMES} "$Laut_autosomal" > {output}
        """

####


#### RULES FOR CALCULATING HETEROZYGOSITY ####
rule CALC_HET_PER_CHR:
    input:
        VAR_FILE="{CLADE}/{SPEC_NAME}/temp/{SPEC_NAME}_Var_Only.txt", 
        CHROM_LENGTH_FILE="{CLADE}/{SPEC_NAME}/Reference_{SPEC_NAME}_Chroms_Lengths.txt", 
        ROH="{CLADE}/{SPEC_NAME}/{SPEC_NAME}_ROH_Results.csv"
    output:
        "{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{CHROM}_het.txt"
    params:
        TODAY_DATE=TODAY_DATE,
        REF_NAME=REF_NAME, 
        CLADE=CLADE, 
        NUM_AUT_CHROMOSOMES=NUM_AUT_CHROMOSOMES, 
        SPEC_NAME=SPEC_NAME, 
        WINDOW_INTERVAL=WINDOW_INTERVAL,  
        WINDOW_LENGTH=WINDOW_LENGTH
    shell:
        """
        python {CALC_HET_PER_CHR_PY} {input.VAR_FILE} {input.CHROM_LENGTH_FILE} {params.CLADE} {params.SPEC_NAME} {params.TODAY_DATE} {params.WINDOW_LENGTH} {params.WINDOW_INTERVAL} {wildcards.CHROM} {input.ROH}
        """

rule CALC_HET_WHOLE_GENOME:
    input:
        CHROM_LENGTH_FILE="{CLADE}/{SPEC_NAME}/Reference_{SPEC_NAME}_Chroms_Lengths.txt",
        PER_CHR_FILES=expand("{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{CHROM}_het.txt", CLADE=CLADE, SPEC_NAME=SPEC_NAME, TODAY_DATE=TODAY_DATE, CHROM=AUTO_CHROMS)
    output:
        WHOLE_PER_CHR="{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{SPEC_NAME}_per_chr_mean_heterozygosity.txt",
        WHOLE_GENOME="{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{SPEC_NAME}_whole_genome_mean_heterozygosity.txt", 
    params:
        TODAY_DATE=TODAY_DATE,
        CLADE=CLADE, 
        SPEC_NAME=SPEC_NAME,
        NUM_AUT_CHROMOSOMES=NUM_AUT_CHROMOSOMES
    shell:
        """
        python {CALC_HET_WHOLE_GENOME_PY} {input.CHROM_LENGTH_FILE} {params.CLADE} {params.SPEC_NAME} {params.TODAY_DATE} {params.NUM_AUT_CHROMOSOMES}
        """

rule COMPILE_HET:
    input:
        expand("{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{CHROM}_het.txt", CHROM=AUTO_CHROMS, CLADE=CLADE, SPEC_NAME=SPEC_NAME, TODAY_DATE=TODAY_DATE, REF_NAME=REF_NAME)
    output:
        "{CLADE}/{SPEC_NAME}/{SPEC_NAME}_Het_Compiled.tsv"
    shell:
        """
        awk -v var='chrom' -F',' '{{if ($1!=var) {{ print $1, $2, $3, $4, $5, $6, $7, $8, $9 }}}}' {input} >> {output}
        """
####


#### RULES FOR PLOTTING HETEROZYGOSITY ####

rule PLOT_HET_PER_CHR:
    input:
        "{CLADE}/{SPEC_NAME}/{SPEC_NAME}_Het_Compiled.tsv"
    output:
        "{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{SPEC_NAME}_{CHR}_Het_Map.svg"
    params:
        TODAY_DATE=TODAY_DATE,
        REF_NAME=REF_NAME, 
        CLADE=CLADE, 
        SPEC_NAME=SPEC_NAME, 
        NUM_ALL_CHR=NUM_ALL_CHR
    shell:
        """
        Rscript {PLOT_HET_PER_CHR_R} {input} {params.TODAY_DATE} {params.REF_NAME} {params.CLADE} {params.NUM_ALL_CHR} {params.SPEC_NAME}
        """

rule PLOT_WHOLE_HET:
    input:
        "{CLADE}/{SPEC_NAME}/{SPEC_NAME}_Het_Compiled.tsv"
    output:
        "{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{SPEC_NAME}_Het_Whole_Genome_Map.svg"
    params:
        TODAY_DATE=TODAY_DATE,
        REF_NAME=REF_NAME, 
        CLADE=CLADE, 
        NUM_AUT_CHROMOSOMES=NUM_AUT_CHROMOSOMES, 
        SPEC_NAME=SPEC_NAME, 
    shell:
        """
        Rscript {PLOT_WHOLE_HET_R} {input} {params.TODAY_DATE} {params.REF_NAME} {params.CLADE} {params.NUM_AUT_CHROMOSOMES} {params.SPEC_NAME}
        """
####


#### RULES FOR CREATING MSMC INPUT FILES ####
# rule SEP_PAF_BY_CHR:
#     input:
#         "{CLADE}/{SPEC_NAME}/{SPEC_NAME}_FASTGA.chain.chr.fltr.srt.paf"
#     output:
#         "{CLADE}/{SPEC_NAME}/MSMC/paf_files/{SPEC_NAME}_{CHROM}.paf"
#     shell:
#         """
#         mkdir -p "{CLADE}/{SPEC_NAME}/MSMC/paf_files"
#         awk -v var={wildcards.CHROM} '{{if($6 == var) print}}' {input} > {output}
#         """

# rule GEN_CHROM_VCF:
#     input:
#         paf="{CLADE}/{SPEC_NAME}/MSMC/paf_files/{SPEC_NAME}_{CHROM}.paf",
#         refseq=expand("/rds/project/rds-p67MZilb2eQ/projects/VGP/250430.VGP-Phase1/alignment/reference/{CLADE}/{REF_NAME}.fa.gz", CLADE=CLADE, REF_NAME=REF_NAME)
#     output:
#         "{CLADE}/{SPEC_NAME}/MSMC/vcf_files/{SPEC_NAME}_{CHROM}.vcf.gz"
#     shell:
#         """
#         mkdir -p "{wildcards.CLADE}/{wildcards.SPEC_NAME}/MSMC/vcf_files"
#         k8 paftools.js call -s {wildcards.SPEC_NAME} -f {input.refseq} {input.paf} > {wildcards.CLADE}/{wildcards.SPEC_NAME}/MSMC/vcf_files/{wildcards.SPEC_NAME}_{wildcards.CHROM}.vcf
#         gzip {wildcards.CLADE}/{wildcards.SPEC_NAME}/MSMC/vcf_files/{wildcards.SPEC_NAME}_{wildcards.CHROM}.vcf
#         """

# rule GEN_ROH_NEGATIVE_MASK:
#     input:
#         lambda wildcards: f"{wildcards.CLADE}/{wildcards.SPEC_NAME}/{TODAY_DATE}_{wildcards.CHROM}_ROH_Results.txt"
#     output:
#         "{CLADE}/{SPEC_NAME}/MSMC/ROH_Negative_Mask/{SPEC_NAME}_{CHROM}_ROH_Mask.bed.gz"
#     shell:
#         """
#         mkdir -p "{wildcards.CLADE}/{wildcards.SPEC_NAME}/MSMC/ROH_Negative_Mask"
#         python {ROH_MASKER} {wildcards.CLADE} {wildcards.SPEC_NAME} {wildcards.CHROM} {input}
#         gzip {wildcards.CLADE}/{wildcards.SPEC_NAME}/MSMC/ROH_Negative_Mask/{wildcards.SPEC_NAME}_{wildcards.CHROM}_ROH_Mask.bed
#         """

rule MAKE_MASK:
    input:
        "{CLADE}/{SPEC_NAME}/temp/{TODAY_DATE}_{CHR}_Aln_Only.txt"
    output:

    shell:
        """
        """

# rule MAIN_MULTIHETSEP:
#     input:
#         vcf="{CLADE}/{SPEC_NAME}/MSMC/vcf_files/{SPEC_NAME}_{CHROM}.vcf.gz", 
#         ROH_negative_mask="{CLADE}/{SPEC_NAME}/MSMC/ROH_Negative_Mask/{SPEC_NAME}_{CHROM}_ROH_Mask.bed.gz"
#     output:
#         "{CLADE}/{SPEC_NAME}/MSMC/Output_NEW_multihetsep/{CHROM}_multihet.txt"
#     shell:
#         """
#         mkdir -p "{CLADE}/{SPEC_NAME}/MSMC/Output_NEW_multihetsep"
#         python {GENERATE_MULTIHETSEP} {input.vcf} --negative_mask={input.ROH_negative_mask} > {output}
#         """















# rule GEN_INPUT_MASK_VCF:
#     input:
#         "{CLADE}/{SPEC_NAME}/MSMC/vcf_files/{SPEC_NAME}_{CHROM}.vcf.gz"
#     output:
#         mask="{CLADE}/{SPEC_NAME}/MSMC/Output_Mask/{SPEC_NAME}_{CHROM}.bed.gz", 
#         vcf="{CLADE}/{SPEC_NAME}/MSMC/Output_VCF/{SPEC_NAME}_{CHROM}_only_snps.vcf.gz"
#     shell:
#         """
#         mkdir -p "{wildcards.CLADE}/{wildcards.SPEC_NAME}/MSMC/Output_Mask"
#         mkdir -p "{wildcards.CLADE}/{wildcards.SPEC_NAME}/MSMC/Output_VCF"
#         bcftools view -V indels -m2 -M2 {input} | {BAM_CALLER} 1 {output.mask} | bgzip -c > {output.vcf}
#         """
### bcftools call --ploidy 2 -c -V indels {input} | {BAM_CALLER} 1 {output.mask} | bgzip -c > {output.vcf}

# rule MAIN_MULTIHETSEP:
#     input:
#         vcf="{CLADE}/{SPEC_NAME}/MSMC/Output_VCF/{SPEC_NAME}_{CHROM}_only_snps.vcf.gz", 
#         ROH_negative_mask="{CLADE}/{SPEC_NAME}/MSMC/ROH_Negative_Mask/{SPEC_NAME}_{CHROM}_ROH_Mask.bed.gz", 
#         mask="{CLADE}/{SPEC_NAME}/MSMC/Output_Mask/{SPEC_NAME}_{CHROM}.bed.gz"
#     output:
#         "{CLADE}/{SPEC_NAME}/MSMC/Output_NEW_multihetsep/{CHROM}_multihet.txt"
#     shell:
#         """
#         mkdir -p "{CLADE}/{SPEC_NAME}/MSMC/Output_NEW_multihetsep"
#         python {GENERATE_MULTIHETSEP} {input.vcf} --negative_mask={input.ROH_negative_mask} --mask={input.mask} > {output}
#         """

# rule BOOTSTRAPPING_MULTIHET_FILE:
#     input:
#         "{CLADE}/{SPEC_NAME}/MSMC/Output_multihet_filtered_95/{CHROM}_multihet.txt"
#     output:
#         "{CLADE}/{SPEC_NAME}/MSMC/Output_Multihet_Bootstrapped/{BOOT}/{CHROM}_multihet.txt"
#     shell:
#         """
#         mkdir -p {CLADE}/{SPEC_NAME}/MSMC/Output_Multihet_Bootstrapped/{wildcards.BOOT}
#         python {BOOTSTRAPPING_GENERATOR} --inmhs {input} --windowsize 5e+05 --outmhs {wildcards.CLADE}/{wildcards.SPEC_NAME}/MSMC/Output_Multihet_Bootstrapped/{wildcards.BOOT}/{wildcards.CHROM}_multihet.txt
#         """

#### 

#### RULES FOR RUNNING MSMC ####
rule RUN_PRIMARY_MSMC:
    input:
        expand("{CLADE}/{SPEC_NAME}/MSMC/Output_multihet_filtered_95/{CHROM}_multihet.txt", CHROM=AUTO_CHROMS, CLADE=CLADE, SPEC_NAME=SPEC_NAME)
    output:
        loop="{CLADE}/{SPEC_NAME}/MSMC/Primary_Results/{TODAY_DATE}_{SPEC_NAME}.msmc2.loop.txt", 
        log="{CLADE}/{SPEC_NAME}/MSMC/Primary_Results/{TODAY_DATE}_{SPEC_NAME}.msmc2.log", 
        final="{CLADE}/{SPEC_NAME}/MSMC/Primary_Results/{TODAY_DATE}_{SPEC_NAME}.msmc2.final.txt"
    shell:
        """
        mkdir -p "{CLADE}/{SPEC_NAME}/MSMC/Primary_Results"
        build/release/msmc2 -t 12 --fixedRecombination -o {CLADE}/{SPEC_NAME}/{TODAY_DATE}_{SPEC_NAME}.msmc2 {input}
        """


# rule PLOT_MSMC:
#     input:
#         "{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{SPEC_NAME}.msmc2.final.txt"
#     output:
#         "{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{SPEC_NAME}_MSMC2_test.svg"
#     params:
#         CLADE=CLADE, 
#         TODAY_DATE=TODAY_DATE, 
#         SPEC_NAME=SPEC_NAME, 
#         MU=MU, 
#         GEN_TIME=GEN_TIME
#     shell:
#         """
#         Rscript {PLOT_MSMC_R} {CLADE} {TODAY_DATE} {SPEC_NAME} {MU} {GEN_TIME} {output} {input}
#         """

# rule PLOT_MSMC_BOOTSTRAP:
#     input:
#         "{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{SPEC_NAME}.msmc2.final.txt"
#     output:
#         "{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{SPEC_NAME}_MSMC2_Bootstrap.svg"
#     params:
#         CLADE=CLADE, 
#         TODAY_DATE=TODAY_DATE, 
#         SPEC_NAME=SPEC_NAME, 
#         MU=MU, 
#         GEN_TIME=GEN_TIME
#     shell:
#         """
#         Rscript {PLOT_MSMC_BOOTSTRAP_R} {CLADE} {TODAY_DATE} {SPEC_NAME} {MU} {GEN_TIME} {output} {input}
#         """




















# rule MAKE_SAM:
#     input:
#         REF=expand("/rds/project/rds-p67MZilb2eQ/projects/VGP/250430.VGP-Phase1/alignment/reference/{CLADE}/{REF_NAME}.fa.gz", CLADE=CLADE, REF_NAME=REF_NAME), 
#         ALT=expand("/rds/project/rds-p67MZilb2eQ/projects/VGP/250430.VGP-Phase1/alignment/alternate/{CLADE}/{ALT_NAME}.fa.gz", CLADE=CLADE, ALT_NAME=ALT_NAME)
#     output:
#         "{CLADE}/{SPEC_NAME}/{SPEC_NAME}.sam"
#     resources:
#         runtime="12h", 
#         mem_mb=200000
#     params:
#         BOOTSTRAPS=BOOTSTRAPS
#     shell:
#         """
#         mkdir -p {CLADE}/{SPEC_NAME}/MSMC
#         mkdir -p {CLADE}/{SPEC_NAME}/MSMC/Output_BED_file
#         mkdir -p {CLADE}/{SPEC_NAME}/MSMC/Output_indel_BED_file
#         mkdir -p {CLADE}/{SPEC_NAME}/MSMC/Output_bam_caller_BED_file
#         mkdir -p {CLADE}/{SPEC_NAME}/MSMC/Output_multihet_filtered_95
#         mkdir -p {CLADE}/{SPEC_NAME}/MSMC/Output_Multihet_Bootstrapped/{params.BOOTSTRAPS}
#         minimap2 -L -t 16 -ax asm5 {input.REF} {input.ALT} > {output}
#         """

# rule MAKE_BAM:
#     input:
#         "{CLADE}/{SPEC_NAME}/{SPEC_NAME}.sam"
#     output:
#         "{CLADE}/{SPEC_NAME}/{SPEC_NAME}.bam"
#     shell:
#         """
#         samtools sort {input} -o {output}
#         """

# rule INDEX_WHOLE_BAM:
#     input:
#         "{CLADE}/{SPEC_NAME}/{SPEC_NAME}.bam"
#     output:
#         "{CLADE}/{SPEC_NAME}/{SPEC_NAME}.bam.bai"
#     shell:
#         """
#         samtools index {input}
#         """

# rule SEPARATE_BAM_BY_CHR:
#     input:
#         "{CLADE}/{SPEC_NAME}/{SPEC_NAME}.bam"
#     output:
#         "{CLADE}/{SPEC_NAME}/{SPEC_NAME}_{CHROM}.bam"
#     shell:
#         """
#         samtools view -b {input} {wildcards.CHROM} > {output}
#         """

# rule UNZIP_REFERENCE:
#     input:
#         "/rds/project/rds-p67MZilb2eQ/projects/VGP/250430.VGP-Phase1/alignment/reference/{CLADE}/{REF_NAME}.fa.gz"
#     output:
#         "{CLADE}/{SPEC_NAME}/temp/{REF_NAME}.fa"
#     shell:
#         """
#         gunzip -c {input} > {output}
#         """

# rule EXTRACT_CHR_FASTA:
#     input:
#         expand("{CLADE}/{SPEC_NAME}/temp/{REF_NAME}.fa", CLADE=CLADE, SPEC_NAME=SPEC_NAME, REF_NAME=REF_NAME)
#     output:
#         "{CLADE}/{SPEC_NAME}/temp/{CHROM}.fasta"
#     params:
#         CLADE=CLADE, 
#         SPEC_NAME=SPEC_NAME, 
#         REF_NAME=REF_NAME
#     shell:
#         """
#         samtools faidx {input} {wildcards.CHROM} > {output}
#         samtools faidx {output}
#         """

# rule INDEX_BAM_FILE:
#     input:
#         "{CLADE}/{SPEC_NAME}/{SPEC_NAME}_{CHROM}.bam"
#     output:
#         "{CLADE}/{SPEC_NAME}/{SPEC_NAME}_{CHROM}.bam.bai"
#     shell:
#         """
#         samtools index {input}
#         """

# rule PILEUP_VARIANTS:
#     input:
#         FASTA="{CLADE}/{SPEC_NAME}/temp/{CHROM}.fasta",
#         BAM="{CLADE}/{SPEC_NAME}/{SPEC_NAME}_{CHROM}.bam"
#     output:
#         "{CLADE}/{SPEC_NAME}/temp/{CHROM}_pileup.vcf.gz"
#     shell:
#         """
#         bcftools mpileup -f {input.FASTA} {input.BAM} --config pacbio-ccs -Oz -o {output}
#         """

# rule INDEL_MASKER:
#     input:
#         "{CLADE}/{SPEC_NAME}/temp/{CHROM}_pileup.vcf.gz"
#     output:
#         "{CLADE}/{SPEC_NAME}/MSMC/Output_indel_BED_file/{CHROM}_indel_bed.txt"
#     shell:
#         """
#         python {INDEL_MASKER} --input_gen_vcf {input} --Chromosome_name {wildcards.CHROM} --output_BED {output}
#         """

# rule CALL_PLOIDY:
#     input:
#         "{CLADE}/{SPEC_NAME}/temp/{CHROM}_pileup.vcf.gz"
#     output:
#         output_bam_caller="{CLADE}/{SPEC_NAME}/MSMC/Output_bam_caller_BED/{CHROM}_unmodified_bed.txt.gz", 
#         only_snps="{CLADE}/{SPEC_NAME}/temp/{CHROM}_only_snps.vcf.bgz"
#     shell:
#         """
#         bcftools call --ploidy 2 -c -V indels {input} | {BAM_CALLER} 1 {output.output_bam_caller} | bgzip -c > {output.only_snps}
#         """

# rule BAM_CALLER_BED:
#     input:
#         "{CLADE}/{SPEC_NAME}/MSMC/Output_bam_caller_BED/{CHROM}_unmodified_bed.txt.gz"
#     output:
#         "{CLADE}/{SPEC_NAME}/MSMC/Output_bam_caller_BED/{CHROM}_modified_bed.txt.gz"
#     shell:
#         """
#         zcat {input}  | \
#         awk -v chrom="{wildcards.CHROM}" 'BEGIN{{OFS="\t"}} {{ $1=chrom; print $0 }}' | \
#         gzip -c > {output}
#         """

# rule TABIX:
#     input:
#         "{CLADE}/{SPEC_NAME}/temp/{CHROM}_only_snps.vcf.bgz"
#     output:
#         "{CLADE}/{SPEC_NAME}/temp/{CHROM}_only_snps.vcf.bgz.tbi"
#     shell:
#         """
#         tabix -p vcf {input}
#         """

# rule RESTRUCTURE:
#     input:
#         "{CLADE}/{SPEC_NAME}/temp/{CHROM}_only_snps.vcf.bgz"
#     output:
#         "{CLADE}/{SPEC_NAME}/temp/{CHROM}_restructured.vcf.gz"
#     shell:
#         """
#         zcat {input} | sed  's/1\/1/1\/0/g' | bgzip -c > {output}
#         """

# rule REMOVE_INDELS:
#     input:
#         "{CLADE}/{SPEC_NAME}/temp/{CHROM}_restructured.vcf.gz"
#     output:
#         "{CLADE}/{SPEC_NAME}/{CHROM}.recode.vcf"
#     params:
#         CLADE=CLADE,
#         SPEC_NAME=SPEC_NAME
#     shell:
#         """
#         vcftools --gzvcf {input} \
#          --remove-indels \
#          --recode \
#          --out {params.CLADE}/{params.SPEC_NAME}/{wildcards.CHROM}
#         """

# rule BGZIP_VCF:
#     input:
#         "{CLADE}/{SPEC_NAME}/{CHROM}.recode.vcf"
#     output:
#         "{CLADE}/{SPEC_NAME}/{CHROM}.recode.vcf.gz"
#     shell:
#         """
#         bgzip {input}
#         """

# rule GENERATE_MULTIHETSEP:
#     input:
#         "{CLADE}/{SPEC_NAME}/{CHROM}.recode.vcf.gz"
#     output:
#         "{CLADE}/{SPEC_NAME}/MSMC/Output_Multihet/{CHROM}_multihet.txt"
#     shell:
#         """
#         python {GENERATE_MULTIHETSEP} {input} > {output}
#         """

# rule GENERATE_MASK:
#     input:
#         "{CLADE}/{SPEC_NAME}/MSMC/Output_Multihet/{CHROM}_multihet.txt"
#     output:
#         "{CLADE}/{SPEC_NAME}/MSMC/Output_Bed/{CHROM}/{CHROM}_95.txt"
#     params:
#         CLADE=CLADE, 
#         SPEC_NAME=SPEC_NAME
#     shell:
#         """
#         python {MASK_FILE_GENERATOR} --input_file {input} --Chromosome_name {wildcards.CHROM} --output_BED {params.CLADE}/{params.SPEC_NAME}/MSMC/Output_Bed/{wildcards.CHROM}
#         """

# rule MAIN_MULTIHETSEP:
#     input:
#         vcf="{CLADE}/{SPEC_NAME}/{CHROM}.recode.vcf.gz", 
#         negative_mask="{CLADE}/{SPEC_NAME}/MSMC/Output_Bed/{CHROM}/{CHROM}_95.txt", 
#         mask="{CLADE}/{SPEC_NAME}/MSMC/Output_bam_caller_BED/{CHROM}_modified_bed.txt.gz", 
#         negative_mask_indel="{CLADE}/{SPEC_NAME}/MSMC/Output_indel_BED_file/{CHROM}_indel_bed.txt"
#     output:
#         "{CLADE}/{SPEC_NAME}/MSMC/Output_multihet_filtered_95/{CHROM}_multihet.txt"
#     shell:
#         """
#         python {GENERATE_MULTIHETSEP} {input.vcf} --negative_mask={input.negative_mask} --mask={input.mask} --negative_mask={input.negative_mask_indel} > {output}
#         """

# rule BOOTSTRAPPING_MULTIHET_FILE:
#     input:
#         "{CLADE}/{SPEC_NAME}/MSMC/Output_multihet_filtered_95/{CHROM}_multihet.txt"
#     output:
#         "{CLADE}/{SPEC_NAME}/MSMC/Output_Multihet_Bootstrapped/{BOOT}/{CHROM}_multihet.txt"
#     shell:
#         """
#         mkdir -p {CLADE}/{SPEC_NAME}/MSMC/Output_Multihet_Bootstrapped/{wildcards.BOOT}
#         python {BOOTSTRAPPING_GENERATOR} --inmhs {input} --windowsize 5e+05 --outmhs {wildcards.CLADE}/{wildcards.SPEC_NAME}/MSMC/Output_Multihet_Bootstrapped/{wildcards.BOOT}/{wildcards.CHROM}_multihet.txt
#         """
# ###

# #### RULES FOR RUNNING MSMC ANALYSIS ####
# rule RUN_PRIMARY_MSMC:
#     input:
#         expand("{CLADE}/{SPEC_NAME}/MSMC/Output_multihet_filtered_95/{CHROM}_multihet.txt", CHROM=AUTO_CHROMS, CLADE=CLADE, SPEC_NAME=SPEC_NAME)
#     output:
#         loop="{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{SPEC_NAME}.msmc2.loop.txt", 
#         log="{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{SPEC_NAME}.msmc2.log", 
#         final="{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{SPEC_NAME}.msmc2.final.txt"
#     shell:
#         """
#         build/release/msmc2 -t 12 --fixedRecombination -o {CLADE}/{SPEC_NAME}/{TODAY_DATE}_{SPEC_NAME}.msmc2 {input}
#         """

# rule RUN_BOOTSTRAPPING_MSMC:
#     input:
#         expand("{CLADE}/{SPEC_NAME}/MSMC/Output_Multihet_Bootstrapped/{BOOT}/{CHROM}_multihet.txt", CLADE=CLADE, SPEC_NAME=SPEC_NAME, BOOT=BOOTSTRAPS, CHROM=AUTO_CHROMS)
#     output:
#         "{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{SPEC_NAME}_Bootstrapping_{BOOT}.msmc2"
#     resources:
#         runtime="6h", 
#         mem_mb=200000
#     shell:
#         """
#         build/release/msmc2 -t 12 --fixedRecombination -o {output} {input}
#         """

# rule PLOT_MSMC:
#     input:
#         "{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{SPEC_NAME}.msmc2.final.txt"
#     output:
#         "{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{SPEC_NAME}_MSMC2_Bootstrap.svg"
#     params:
#         CLADE=CLADE, 
#         TODAY_DATE=TODAY_DATE, 
#         SPEC_NAME=SPEC_NAME, 
#         MU=MU, 
#         GEN_TIME=GEN_TIME
#     shell:
#         """
#         Rscript {PLOT_MSMC_R} {CLADE} {TODAY_DATE} {SPEC_NAME} {MU} {GEN_TIME} {output} {input}
#         """


# # # # #### RULES FOR CREATING BAM FILES AND CALCULATING COVERAGE ####
# # rule MAKE_SAM:
# #     input:
# #         REF=expand("/rds/project/rds-p67MZilb2eQ/projects/VGP/241117.UCSC-hubs-VGP-alignment/alignment/reference/{CLADE}/{REF_NAME}.fa.gz", CLADE=CLADE, REF_NAME=REF_NAME), 
# #         ALT=expand("/rds/project/rds-p67MZilb2eQ/projects/VGP/241117.UCSC-hubs-VGP-alignment/alignment/alternate/{CLADE}/{ALT_NAME}.fa.gz", CLADE=CLADE, ALT_NAME=ALT_NAME)
# #     output:
# #         "{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{SPEC_NAME}.sam"
# #     shell:
# #         """
# #         minimap2 -t 16 -ax asm5 {input.REF} {input.ALT} > {output}
# #         """

# # rule MAKE_BAM:
# #     input:
# #         "{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{SPEC_NAME}.sam"
# #     output:
# #         "{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{SPEC_NAME}.bam"
# #     shell:
# #         """
# #         samtools sort {input} -o {output}
# #         """

# # rule BAM_COVERAGE:
# #     input:
# #         "{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{SPEC_NAME}.bam"
# #     output:
# #         "{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{SPEC_NAME}_BAM_Coverage.txt"
# #     shell:
# #         """
# #         samtools coverage -D {input} > {output}
# #         """

# # ####


# # #### RULES FOR GENERATING PAF DOTPLOT ####
# # rule PAF_DOTPLOT:
# #     input:
# #         f"{CLADE}/{SPEC_NAME}/{SPEC_NAME}_FASTGA.fltr.srt.paf",
# #     output:
# #         f"{CLADE}/{SPEC_NAME}/{SPEC_NAME}_FASTGA_collinear_plot.pdf"
# #     shell:
# #         """
# #         ALNplot -p {input} > {output}
# #         """

# # ####

# # # # #### REMOVE TEMPORARY FILES ####
# # # # rule RM_TEMP_ROH_FILES:
# # # #     input:
# # # #         f"{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{chrom}_filtered.paf",
# # # #         f"{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{chrom}_Aln_Var.txt", 
# # # #         f"{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{chrom}_Aln_Only.txt", 
# # # #         f"{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{chrom}_Var_Only.txt", 
# # # #         f"{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{SPEC_NAME}.sam", 
# # # #         f"{CLADE}/{SPEC_NAME}/{SPEC_NAME}_Aln_Var.txt", 
# # # #         f"{CLADE}/{SPEC_NAME}/{SPEC_NAME}_Var_Only.txt"
# # # #     shell:
# # # #         """
# # # #         rm -f -r {input}
# # # #         """
# # # # ####
