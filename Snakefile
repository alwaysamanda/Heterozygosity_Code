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
    snakemake --executor slurm --jobs 30 --workflow-profile profiles --latency-wait 60 --keep-going --configfile config_files/{SPEC_NAME}_config.yml
Add in the command --touch for when I want to make sure it doesn't regenerate previously made files
For changing input parameters, change in the config.yml file

Generate simplified plot to show rule flow with:
    snakemake --executor slurm --workflow-profile profiles --jobs 10 --configfile config.yml --rulegraph | dot -Tpdf > dag.pdf
Show full complicated path for each chromosome with:
    snakemake --executor slurm --workflow-profile profiles --jobs 10 --configfile config.yml --dag | dot -Tpdf > dag_complex.pdf

NOTES:
    Based on Snakefile made by Bettina Fisher (https://github.com/bef22/vcf_snakemake/blob/main/Snakefile)

AUTHORS and UPDATES:
Amanda Gardiner, 20250313

UPDATES: 20250530 -- Redoing pipeline to streamline it

"""

###############################################################################

user = "ag2427"

### set global parameters
CLADE=config['CLADE']
SPEC_NAME=config['SPEC_NAME']
REF_NAME=config['REF_NAME']
ALT_NAME=config['ALT_NAME']
CHROM_START_CHR=config['CHROM_START_CHR']
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
        expand("{CLADE}/{SPEC_NAME}/temp/{SPEC_NAME}_{REF_NAME}.1gdb", CLADE=CLADE, SPEC_NAME=SPEC_NAME, REF_NAME=REF_NAME),
        expand("{CLADE}/{SPEC_NAME}/temp/{SPEC_NAME}_{ALT_NAME}.1gdb", CLADE=CLADE, SPEC_NAME=SPEC_NAME, ALT_NAME=ALT_NAME),
        expand("{CLADE}/{SPEC_NAME}/temp/{SPEC_NAME}_{REF_NAME}.gix", CLADE=CLADE, SPEC_NAME=SPEC_NAME, REF_NAME=REF_NAME), 
        expand("{CLADE}/{SPEC_NAME}/temp/{SPEC_NAME}_{ALT_NAME}.gix", CLADE=CLADE, SPEC_NAME=SPEC_NAME, ALT_NAME=ALT_NAME),
        expand("{CLADE}/{SPEC_NAME}/{SPEC_NAME}_ALN.1aln", CLADE=CLADE, SPEC_NAME=SPEC_NAME), 
        f"{CLADE}/{SPEC_NAME}/{SPEC_NAME}_FASTGA.chain.paf", 
        f"{CLADE}/chrom_lists/{SPEC_NAME}_chroms.txt", 
        f"{CLADE}/{SPEC_NAME}/{SPEC_NAME}_FASTGA.chain.chr.fltr.srt.paf", 
        f"{CLADE}/{SPEC_NAME}/{SPEC_NAME}_ALN.chain.pdf", 
        ### End alignment section
        ### Start chromosome information section
        f"{CLADE}/chrom_lists/{SPEC_NAME}_chroms.txt", 
        f"{CLADE}/{SPEC_NAME}/{SPEC_NAME}_Chroms_Lengths.txt", 
        f"{CLADE}/{SPEC_NAME}/temp/{SPEC_NAME}_Var_Only.txt", 
        f"{CLADE}/{SPEC_NAME}/temp/{SPEC_NAME}_Aln_Only.txt", 
        expand("{CLADE}/{SPEC_NAME}/temp/{CHR}_Aln_Only.txt", CLADE=CLADE, SPEC_NAME=SPEC_NAME, CHR=CHROMS), 
        expand("{CLADE}/{SPEC_NAME}/temp/{CHR}_Var_Only.txt", CLADE=CLADE, SPEC_NAME=SPEC_NAME, CHR=CHROMS), 
        ### End chromosome information section
        ### Start ROH section
        expand("{CLADE}/{SPEC_NAME}/{CHR}_ROH_Results.txt", CLADE=CLADE, SPEC_NAME=SPEC_NAME, CHR=CHROMS), 
        f"{CLADE}/{SPEC_NAME}/{SPEC_NAME}_ROH_Results.csv", 
        f"{CLADE}/{SPEC_NAME}/{SPEC_NAME}_ROH_Map.pdf", 
        expand("{CLADE}/{SPEC_NAME}/{CHR}_{SPEC_NAME}_FROH.txt", CLADE=CLADE, SPEC_NAME=SPEC_NAME, CHR=CHROMS), 
        f"{CLADE}/{SPEC_NAME}/{SPEC_NAME}_FROH.txt", 
        expand("{CLADE}/{SPEC_NAME}/{CHROM}_het.txt", CLADE=CLADE, SPEC_NAME=SPEC_NAME, CHROM=AUTO_CHROMS), 
        f"{CLADE}/{SPEC_NAME}/{SPEC_NAME}_per_chr_mean_heterozygosity.txt", 
        f"{CLADE}/{SPEC_NAME}/{SPEC_NAME}_whole_genome_mean_heterozygosity.txt", 
        f"{CLADE}/{SPEC_NAME}/{SPEC_NAME}_Het_Compiled.tsv", 
        expand("{CLADE}/{SPEC_NAME}/{SPEC_NAME}_{CHROM}_Het_Map.png", CLADE=CLADE, SPEC_NAME=SPEC_NAME, CHROM=AUTO_CHROMS), 
        f"{CLADE}/{SPEC_NAME}/{SPEC_NAME}_Het_Whole_Genome_Map.png"

####---------- RULES FOR FASTGA ALINGMENT AND PAF FILE GENERATION ----------####
rule FAtoGDB_REF:
    input: 
        "/rds/project/rds-p67MZilb2eQ/projects/VGP/250430.VGP-Phase1/alignment/reference/{CLADE}/{REF_NAME}.fa.gz"
    output:
        "{CLADE}/{SPEC_NAME}/temp/{SPEC_NAME}_{REF_NAME}.1gdb",
    shell:
        """
        mkdir -p {wildcards.CLADE}
        mkdir -p {wildcards.CLADE}/{wildcards.SPEC_NAME}
        mkdir -p {wildcards.CLADE}/{wildcards.SPEC_NAME}/temp
        FAtoGDB -v {input} {output}
        """

rule FAtoGDB_ALT:
    input: 
        "/rds/project/rds-p67MZilb2eQ/projects/VGP/250430.VGP-Phase1/alignment/alternate/{CLADE}/{ALT_NAME}.fa.gz"
    output:
        "{CLADE}/{SPEC_NAME}/temp/{SPEC_NAME}_{ALT_NAME}.1gdb",
    shell:
        """
        FAtoGDB -v {input} {output}
        """

rule GIXmake_REF:
    input: 
        "{CLADE}/{SPEC_NAME}/temp/{SPEC_NAME}_{REF_NAME}.1gdb"
    output:
        "{CLADE}/{SPEC_NAME}/temp/{SPEC_NAME}_{REF_NAME}.gix"
    shell:
        """
        GIXmake -v -P. -T8 {input}
        """

rule GIXmake_ALT:
    input: 
        "{CLADE}/{SPEC_NAME}/temp/{SPEC_NAME}_{ALT_NAME}.1gdb"
    output:
        "{CLADE}/{SPEC_NAME}/temp/{SPEC_NAME}_{ALT_NAME}.gix"
    shell:
        """
        GIXmake -v -P. -T8 {input}
        """

rule FASTGA:
    input:
        REF=expand("{CLADE}/{SPEC_NAME}/temp/{SPEC_NAME}_{REF_NAME}.gix", CLADE=CLADE, SPEC_NAME=SPEC_NAME, REF_NAME=REF_NAME),
        ALT=expand("{CLADE}/{SPEC_NAME}/temp/{SPEC_NAME}_{ALT_NAME}.gix", CLADE=CLADE, SPEC_NAME=SPEC_NAME, ALT_NAME=ALT_NAME)
    output:
        "{CLADE}/{SPEC_NAME}/{SPEC_NAME}_ALN.1aln"
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

rule GET_CHROM_LISTS:
    input:
        expand("../250430.VGP-Phase1/alignment/reference/{CLADE}/{REF_NAME}.fa.gz", CLADE=CLADE, REF_NAME=REF_NAME)
    output:
        "{CLADE}/chrom_lists/{SPEC_NAME}_chroms.txt"
    params:
        CHROM_START_CHR=CHROM_START_CHR
    shell:
        """
        mkdir -p {wildcards.CLADE}/chrom_lists
        zcat < {input} | grep '>{params.CHROM_START_CHR}' | sed 's/^>//' > {output}
        """

rule FILTER_PAF_CHR_ONLY:
    input:
        PAF="{CLADE}/{SPEC_NAME}/{SPEC_NAME}_FASTGA.chain.fltr.srt.paf", 
        ALL_CHROMS="{CLADE}/chrom_lists/{SPEC_NAME}_chroms.txt"
    output:
        "{CLADE}/{SPEC_NAME}/{SPEC_NAME}_FASTGA.chain.chr.fltr.srt.paf"
    shell:
        """
        awk 'BEGIN {{ while (getline < "{input.ALL_CHROMS}") list[$0] }} $6 in list' {input.PAF} > {output}
        """

rule ALNPLOT:
    input:
        "{CLADE}/{SPEC_NAME}/{SPEC_NAME}_ALN.chain.1aln"
    output:
        "{CLADE}/{SPEC_NAME}/{SPEC_NAME}_ALN.chain.pdf"
    shell:
        """
        ALNplot -p -H500 {input}
        mv {wildcards.SPEC_NAME}_ALN.chain.pdf {wildcards.CLADE}/{wildcards.SPEC_NAME}/{wildcards.SPEC_NAME}_ALN.chain.pdf
        """
####---------- END ----------####

####---------- RULES TO GET CHROMOSOME INFORMATION ----------####
rule CHROM_LENGTH_CALC:
    input:
        FASTA=expand("../250430.VGP-Phase1/alignment/reference/{CLADE}/{REF_NAME}.fa.gz", CLADE=CLADE, REF_NAME=REF_NAME), 
        ALL_CHROMS="{CLADE}/chrom_lists/{SPEC_NAME}_chroms.txt",
    output:
        "{CLADE}/{SPEC_NAME}/{SPEC_NAME}_Chroms_Lengths.txt"
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
        chr_aln_only="{CLADE}/{SPEC_NAME}/temp/{CHR}_Aln_Only.txt"
    shell:
        """
        awk '$1 == "R" && $2 == "{wildcards.CHR}"' {input} > {output.chr_aln_only}
        """

rule GET_VAR_ONLY_PER_CHROM:
    input:
        f"{CLADE}/{SPEC_NAME}/temp/{SPEC_NAME}_Var_Only.txt"
    output:
        chr_var_only="{CLADE}/{SPEC_NAME}/temp/{CHR}_Var_Only.txt"
    shell:
        """
        awk '$1 == "V" && $2 == "{wildcards.CHR}"' {input} > {output.chr_var_only}
        """
####---------- END ----------####

####---------- RULES FOR CALCULATING ROH ----------####
rule ROH_CALC:
    input:
        chr_aln_only="{CLADE}/{SPEC_NAME}/temp/{CHR}_Aln_Only.txt", 
        chr_var_only="{CLADE}/{SPEC_NAME}/temp/{CHR}_Var_Only.txt", 
        CHROM_LENGTH_FILE=f"{CLADE}/{SPEC_NAME}/{SPEC_NAME}_Chroms_Lengths.txt"
    output:
        ROH_outfiles="{CLADE}/{SPEC_NAME}/{CHR}_ROH_Results.txt"
    params:  
        REF_NAME=REF_NAME, 
        CLADE=CLADE, 
        SPEC_NAME=SPEC_NAME
    shell:
        """
        CHROM_LENGTH="$(awk -v var='{wildcards.CHR}' '$1 == var {{print $2; exit}}' {input.CHROM_LENGTH_FILE})"
        python {ROH_CALC_PY} {wildcards.CHR} $CHROM_LENGTH {params.REF_NAME} {params.CLADE} {input.chr_aln_only} {input.chr_var_only} {params.SPEC_NAME}
        """

rule COMPILE_ROH:
    input:
        expand("{CLADE}/{SPEC_NAME}/{CHR}_ROH_Results.txt", CLADE=CLADE, SPEC_NAME=SPEC_NAME, CHR=CHROMS)
    output:
        "{CLADE}/{SPEC_NAME}/{SPEC_NAME}_ROH_Results.csv"
    params:
        CLADE=CLADE, 
        SPEC_NAME=SPEC_NAME, 
    shell:
        """
        awk -v var='chrom' -F',' '{{if ($1!=var) {{ print $0 }}}}' {input} >> {output}
        """
####---------- END ----------####

####---------- RULES FOR PLOTTING ROH ----------####
rule Plot_ROH:
    input:
        ROH="{CLADE}/{SPEC_NAME}/{SPEC_NAME}_ROH_Results.csv",
        CHROM_LENGTH_FILE="{CLADE}/{SPEC_NAME}/{SPEC_NAME}_Chroms_Lengths.txt"
    output:
        "{CLADE}/{SPEC_NAME}/{SPEC_NAME}_ROH_Map.pdf"
    params:
        REF_NAME=REF_NAME,
        CLADE=CLADE,
        NUM_AUT_CHROMOSOMES=NUM_AUT_CHROMOSOMES, 
        SPEC_NAME=SPEC_NAME, 
        NUM_ALL_CHR=NUM_ALL_CHR
    shell:
        """
        Rscript {ROH_PLOT_R} {input.ROH} {input.CHROM_LENGTH_FILE} {params.REF_NAME} {params.CLADE} {params.NUM_AUT_CHROMOSOMES} {params.SPEC_NAME} {params.NUM_ALL_CHR} > {output}
        """
####---------- END ----------####

####---------- RULES FOR CALCULATING FROH ----------####
rule FROH_PER_AUT_CHR:
    input:
        ROH="{CLADE}/{SPEC_NAME}/{SPEC_NAME}_ROH_Results.csv", 
        CHROM_LENGTH_FILE="{CLADE}/{SPEC_NAME}/{SPEC_NAME}_Chroms_Lengths.txt", 
        ALL_CHROMS="{CLADE}/chrom_lists/{SPEC_NAME}_chroms.txt",
        VAR_FILE="{CLADE}/{SPEC_NAME}/temp/{CHR}_Var_Only.txt"
    output:
        outfiles="{CLADE}/{SPEC_NAME}/{CHR}_{SPEC_NAME}_FROH.txt"
    params:
        REF_NAME=REF_NAME,
        CLADE=CLADE,
        NUM_AUT_CHROMOSOMES=NUM_AUT_CHROMOSOMES, 
        SPEC_NAME=SPEC_NAME
    shell:
        """
        start="$(head -n  1 {input.VAR_FILE} | awk ' {{print $3}} ')"
        end="$(tail -n  1 {input.VAR_FILE} | awk ' {{print $4}} ')"
        chrom_length="$((end - start))"

        Rscript {FROH_CALC_PER_CHR_R} {input.ROH} {input.CHROM_LENGTH_FILE} {params.REF_NAME} {params.CLADE} {params.NUM_AUT_CHROMOSOMES} "$chrom_length" {wildcards.CHR} > {output.outfiles}
        """

rule WHOLE_FROH:
    input:
        ROH="{CLADE}/{SPEC_NAME}/{SPEC_NAME}_ROH_Results.csv", 
        CHROM_LENGTH_FILE="{CLADE}/{SPEC_NAME}/{SPEC_NAME}_Chroms_Lengths.txt", 
        ALL_CHROMS=expand("{CLADE}/chrom_lists/{SPEC_NAME}_chroms.txt", CLADE=CLADE, SPEC_NAME=SPEC_NAME), 
        VAR_FILE="{CLADE}/{SPEC_NAME}/temp/{SPEC_NAME}_Var_Only.txt"
    output:
        "{CLADE}/{SPEC_NAME}/{SPEC_NAME}_FROH.txt"
    params:
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
            start="$(head -n 1 {params.CLADE}/{params.SPEC_NAME}/temp/"$line"_Var_Only.txt | awk ' {{print $3}} ')"
            end="$(tail -n 1 {params.CLADE}/{params.SPEC_NAME}/temp/"$line"_Var_Only.txt | awk ' {{print $4}} ')"
            chrom_length="$((end - start))"
            Laut_autosomal=$(("$Laut_autosomal" + "$chrom_length"))
        done < <(head -n {params.NUM_AUT_CHROMOSOMES} "{input.ALL_CHROMS}") 

        while read -r line; do
            start="$(head -n 1 {params.CLADE}/{params.SPEC_NAME}/temp/"$line"_Var_Only.txt | awk ' {{print $3}} ')"
            end="$(tail -n 1 {params.CLADE}/{params.SPEC_NAME}/temp/"$line"_Var_Only.txt | awk ' {{print $4}} ')"
            chrom_length="$((end - start))"
            Laut=$(("$Laut" + "$chrom_length"))
        done < "{input.ALL_CHROMS}"

        Rscript {FROH_CALC_R} {input.ROH} {input.CHROM_LENGTH_FILE} {params.REF_NAME} {params.CLADE} "$Laut" {params.NUM_AUT_CHROMOSOMES} "$Laut_autosomal" > {output}
        """
####---------- END ----------####

####---------- RULES FOR CALCULATING HETEROZYGOSITY ----------####
rule CALC_HET_PER_CHR:
    input:
        VAR_FILE="{CLADE}/{SPEC_NAME}/temp/{SPEC_NAME}_Var_Only.txt", 
        CHROM_LENGTH_FILE="{CLADE}/{SPEC_NAME}/{SPEC_NAME}_Chroms_Lengths.txt", 
        ROH="{CLADE}/{SPEC_NAME}/{SPEC_NAME}_ROH_Results.csv"
    output:
        "{CLADE}/{SPEC_NAME}/{CHROM}_het.txt"
    params:
        REF_NAME=REF_NAME, 
        CLADE=CLADE, 
        NUM_AUT_CHROMOSOMES=NUM_AUT_CHROMOSOMES, 
        SPEC_NAME=SPEC_NAME, 
        WINDOW_INTERVAL=WINDOW_INTERVAL,  
        WINDOW_LENGTH=WINDOW_LENGTH
    shell:
        """
        python {CALC_HET_PER_CHR_PY} {input.VAR_FILE} {input.CHROM_LENGTH_FILE} {params.CLADE} {params.SPEC_NAME} {params.WINDOW_LENGTH} {params.WINDOW_INTERVAL} {wildcards.CHROM} {input.ROH}
        """

rule CALC_HET_WHOLE_GENOME:
    input:
        CHROM_LENGTH_FILE="{CLADE}/{SPEC_NAME}/{SPEC_NAME}_Chroms_Lengths.txt",
        PER_CHR_FILES=expand("{CLADE}/{SPEC_NAME}/{CHROM}_het.txt", CLADE=CLADE, SPEC_NAME=SPEC_NAME, CHROM=AUTO_CHROMS)
    output:
        WHOLE_PER_CHR="{CLADE}/{SPEC_NAME}/{SPEC_NAME}_per_chr_mean_heterozygosity.txt",
        WHOLE_GENOME="{CLADE}/{SPEC_NAME}/{SPEC_NAME}_whole_genome_mean_heterozygosity.txt", 
    params:
        CLADE=CLADE, 
        SPEC_NAME=SPEC_NAME,
        NUM_AUT_CHROMOSOMES=NUM_AUT_CHROMOSOMES
    shell:
        """
        python {CALC_HET_WHOLE_GENOME_PY} {input.CHROM_LENGTH_FILE} {params.CLADE} {params.SPEC_NAME} {params.NUM_AUT_CHROMOSOMES}
        """

rule COMPILE_HET:
    input:
        expand("{CLADE}/{SPEC_NAME}/{CHROM}_het.txt", CHROM=AUTO_CHROMS, CLADE=CLADE, SPEC_NAME=SPEC_NAME)
    output:
        "{CLADE}/{SPEC_NAME}/{SPEC_NAME}_Het_Compiled.tsv"
    shell:
        """
        awk -v var='chrom' -F',' '{{if ($1!=var) {{ print $1, $2, $3, $4, $5, $6, $7, $8, $9 }}}}' {input} >> {output}
        """
####---------- END ----------####

####---------- RULES FOR PLOTTING HETEROZYGOSITY ----------####
rule PLOT_HET_PER_CHR:
    input:
        "{CLADE}/{SPEC_NAME}/{SPEC_NAME}_Het_Compiled.tsv"
    output:
        "{CLADE}/{SPEC_NAME}/{SPEC_NAME}_{CHROM}_Het_Map.png"
    params:
        REF_NAME=REF_NAME, 
        CLADE=CLADE, 
        SPEC_NAME=SPEC_NAME, 
        NUM_ALL_CHR=NUM_ALL_CHR
    shell:
        """
        Rscript {PLOT_HET_PER_CHR_R} {input} {params.REF_NAME} {params.CLADE} {params.NUM_ALL_CHR} {params.SPEC_NAME}
        """

rule PLOT_WHOLE_HET:
    input:
        "{CLADE}/{SPEC_NAME}/{SPEC_NAME}_Het_Compiled.tsv"
    output:
        "{CLADE}/{SPEC_NAME}/{SPEC_NAME}_Het_Whole_Genome_Map.png"
    params:
        REF_NAME=REF_NAME, 
        CLADE=CLADE, 
        NUM_AUT_CHROMOSOMES=NUM_AUT_CHROMOSOMES, 
        SPEC_NAME=SPEC_NAME, 
    shell:
        """
        Rscript {PLOT_WHOLE_HET_R} {input} {params.REF_NAME} {params.CLADE} {params.NUM_AUT_CHROMOSOMES} {params.SPEC_NAME}
        """
####---------- END ----------####

####---------- RULES FOR MSMC DATA PREP ----------####
####---------- END ----------####

####---------- RULES FOR RUNNING AND PLOTTING MSMC ----------####
####---------- END ----------####
