#### SNAKEFILE FOR AUTOMATED PIPLELINE OF RUNNING ANALYSES ####

"""
login to icelake node using login-icelake.hpc.cam.ac.uk
start a screen terminal (screen -RD)
make sure to use snakemake version 7.8.5

check how many jobs are required AT EACH STAGE, this command also shows the shell commands it would use:
- snakemake -n --printshellcmds -p

to run locally use:
    snakemake --configfile config_files/sharks/{SPEC_NAME}_config.yml
    
    --touch
    --keep-going
to run on slurm use:
    snakemake --executor slurm --jobs 30 --workflow-profile profiles --configfile config_files/fishes/XXXXXX_config.yml
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
UPDATES: 20250618 -- Redoing pipeline to use fastga alingments which Richard created, and put output files in the VGP/groups/CLADE/SPECIES directories
UPDATES: 20250626 -- Updating pipeline for FASTGA alignments which Richard created -- having to use new files because of the alternate format of the paftools output
"""

###############################################################################

user = "ag2427"

pass_aln_paf = "yes"
pass_var = "yes"
pass_roh = "yes"
pass_het = "yes"
pass_primary_msmc = "no"

### set global parameters
CLADE=config['CLADE']
SPEC_NAME=config['SPEC_NAME']
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
FULL_SPEC_NAME=config['FULL_SPEC_NAME']
ALN_NAME=config['ALN_NAME']


### paths to dependent scripts
ROH_CALC_PY = "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/20250626_ROH_Durbin_Calc_Eqns_V7.py"
ROH_PLOT_R = "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/20250707_Plot_ROH_V3.R" ## 20250620_Plot_ROH.R
FROH_CALC_PER_CHR = "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/20250627_FROH_per_chr_calc_V3.py"
FROH_CALC_WHOLE = "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/20250617_FROH_Calc_Whole_Genome_V3.py"
CALC_HET_PER_CHR_PY = "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/20250827_find_het_per_chr_V5.py"
CALC_HET_WHOLE_GENOME_PY = "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/20250728_find_het_whole_genome_V4.py"
PLOT_HET_PER_CHR_R = "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/20250123_Plot_het_per_chr.R"
PLOT_WHOLE_HET_R = "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/20250325_Plot_het_whole_genome.R"
PLOT_MSMC_R = "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/20250508_Plot_MSMC.R"
PLOT_MSMC_BOOSTRAP_R = "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/20250529_Plot_MSMC_Bootstrap.R"
ROH_MASKER = "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/20250528_ROH_Masker.py"
ALN_MASK_MAKER_PY = "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/20250529_Mask_Maker.py"
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
    CHROM_START_CHR=CHROM_START_CHR,
    WINDOW_INTERVAL=WINDOW_INTERVAL, 
    WINDOW_LENGTH=WINDOW_LENGTH, 
    NUM_AUT_CHROMOSOMES=NUM_AUT_CHROMOSOMES, 
    NUM_ALL_CHR=NUM_ALL_CHR, 
    CHROMS=CHROMS, 
    AUTO_CHROMS=AUTO_CHROMS,
    BOOTSTRAPS=BOOTSTRAPS, 
    FULL_SPEC_NAME=FULL_SPEC_NAME, 
    ALN_NAME=ALN_NAME


#### ---- Create list of output files for each stage after passing checkpoints ---- ####
output_files = list()
## Initial output from FastGA aligned file
output_files.append(f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{ALN_NAME}.chain.chr.fltr.srt.paf")
output_files.append(f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{SPEC_NAME}_chroms.txt")
output_files.append(f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{ALN_NAME}.chain.pdf")

## Outputs for Aln and Var Files
if pass_aln_paf == 'yes':
    output_files.append(f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{SPEC_NAME}_Chroms_Lengths.txt")
    output_files.append(f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/temp/{SPEC_NAME}_Var_Only.txt")
    output_files.append(f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/temp/{SPEC_NAME}_Aln_Only.txt")
    output_files.append(expand("../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/temp/{CHR}_Var_Only.txt", CLADE=CLADE, FULL_SPEC_NAME=FULL_SPEC_NAME, CHR=CHROMS))
    output_files.append(expand("../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/temp/{CHR}_Aln_Only.txt", CLADE=CLADE, FULL_SPEC_NAME=FULL_SPEC_NAME, CHR=CHROMS))

if pass_var == 'yes':
    output_files.append(expand("../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/ROH/{CHR}_ROH_Results.txt", CLADE=CLADE, FULL_SPEC_NAME=FULL_SPEC_NAME, CHR=CHROMS))
    output_files.append(f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/ROH/{SPEC_NAME}_ROH_Results.csv")
    output_files.append(f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{SPEC_NAME}_ROH_Map.png")
    output_files.append(expand("../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/ROH/{CHROM}_FROH_results.txt", CLADE=CLADE, FULL_SPEC_NAME=FULL_SPEC_NAME, CHROM=AUTO_CHROMS))
    output_files.append(f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{SPEC_NAME}_FROH.txt")

if pass_roh == 'yes':
    output_files.append(expand("../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/Het/{CHROM}_het.txt", CLADE=CLADE, FULL_SPEC_NAME=FULL_SPEC_NAME, CHROM=AUTO_CHROMS))
    output_files.append(f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{SPEC_NAME}_per_chr_mean_heterozygosity.txt")
    output_files.append(f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{SPEC_NAME}_whole_genome_mean_heterozygosity.txt")
    output_files.append(f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/Het/{SPEC_NAME}_Het_Compiled.tsv")
    output_files.append(f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{SPEC_NAME}_Het_Whole_Genome_Map.png")

if pass_het == 'yes':
    output_files.append(expand("../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/paf_files/{SPEC_NAME}_{CHROM}.paf", CLADE=CLADE, FULL_SPEC_NAME=FULL_SPEC_NAME, SPEC_NAME=SPEC_NAME, CHROM=AUTO_CHROMS))
    output_files.append(expand("../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/vcf_files/{SPEC_NAME}_{CHROM}.vcf.gz", CLADE=CLADE, FULL_SPEC_NAME=FULL_SPEC_NAME, SPEC_NAME=SPEC_NAME, CHROM=AUTO_CHROMS))
    output_files.append(expand("../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/filtered_vcf_files/{SPEC_NAME}_{CHROM}_filtered.vcf.gz", CLADE=CLADE, FULL_SPEC_NAME=FULL_SPEC_NAME, SPEC_NAME=SPEC_NAME, CHROM=AUTO_CHROMS))
    output_files.append(expand("../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/ROH_Negative_Mask/{SPEC_NAME}_{CHROM}_ROH_Mask.bed.gz", CLADE=CLADE, FULL_SPEC_NAME=FULL_SPEC_NAME, SPEC_NAME=SPEC_NAME, CHROM=AUTO_CHROMS))
    output_files.append(expand("../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/Aln_Mask/{CHROM}_Aln_Mask.bed", CLADE=CLADE, FULL_SPEC_NAME=FULL_SPEC_NAME, CHROM=AUTO_CHROMS))
    output_files.append(expand("../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/Output_primary_multihetsep/{SPEC_NAME}_{CHROM}_multihet_new.txt", SPEC_NAME=SPEC_NAME, CLADE=CLADE, FULL_SPEC_NAME=FULL_SPEC_NAME, CHROM=AUTO_CHROMS))
    # output_files.append(f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/Primary_Results/{SPEC_NAME}.msmc2.loop.txt")
    # output_files.append(f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/Primary_Results/{SPEC_NAME}.msmc2.log")
    # output_files.append(f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/Primary_Results/{SPEC_NAME}.msmc2.final.txt")
    # output_files.append(f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/{SPEC_NAME}_MSMC2.png")

if pass_primary_msmc == 'yes':
    output_files.append(expand("../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/Output_Multihet_Bootstrapped/{BOOT}/{SPEC_NAME}_{CHROM}_multihet.txt", SPEC_NAME=SPEC_NAME, CLADE=CLADE, FULL_SPEC_NAME=FULL_SPEC_NAME, BOOT=BOOTSTRAPS, CHROM=AUTO_CHROMS))
    output_files.append(expand("../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/Bootstrap_results/{SPEC_NAME}_Bootstrapping_{BOOT}.msmc2.loop.txt", CLADE=CLADE, FULL_SPEC_NAME=FULL_SPEC_NAME, BOOT=BOOTSTRAPS, CHROM=AUTO_CHROMS, SPEC_NAME=SPEC_NAME))
    output_files.append(expand("../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/Bootstrap_results/{SPEC_NAME}_Bootstrapping_{BOOT}.msmc2.log", CLADE=CLADE, FULL_SPEC_NAME=FULL_SPEC_NAME, BOOT=BOOTSTRAPS, CHROM=AUTO_CHROMS, SPEC_NAME=SPEC_NAME))
    output_files.append(expand("../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/Bootstrap_results/{SPEC_NAME}_Bootstrapping_{BOOT}.msmc2.final.txt", CLADE=CLADE, FULL_SPEC_NAME=FULL_SPEC_NAME, BOOT=BOOTSTRAPS, CHROM=AUTO_CHROMS, SPEC_NAME=SPEC_NAME))
    output_files.append(f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/{SPEC_NAME}_MSMC2_Bootstrapped.png")

#### ---- End ---- ####

rule all: 
    input: 
        output_files

#### ---------- RULES FOR FASTGA ALINGMENT AND PAF FILE GENERATION ---------- ####
rule ALNCHAIN:
    input:
        f"../groups/{CLADE}/{FULL_SPEC_NAME}/{ALN_NAME}-12.1aln"
    output:
        f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{ALN_NAME}.chain.1aln"
    shell:
        """
        mkdir -p ../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity
        ALNchain -o../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{ALN_NAME}.chain.1aln {input}
        """

rule ALNtoPAF:
    input:
        f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{ALN_NAME}.chain.1aln"
    output:
        f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{ALN_NAME}.chain.paf"
    shell:
        """
        ALNtoPAF -s -w -T8 {input} > {output}
        """

rule FILTER_PAF_VARIANCE:
    input:
        f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{ALN_NAME}.chain.paf"
    output:
        f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{ALN_NAME}.chain.fltr.paf"
    shell:
        """
        awk -v OFS="\t" '{{if(substr($13,6)<=0.1) print}}' {input} > {output}
        """

rule SORT_PAF:
    input:
        f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{ALN_NAME}.chain.fltr.paf"
    output:
        f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{ALN_NAME}.chain.fltr.srt.paf"
    shell:
        """
        sort -k6,6V -k8,8n {input} > {output}
        """

rule GET_CHROM_LISTS:
    input:
        f"../groups/{CLADE}/{FULL_SPEC_NAME}/{ALN_NAME}.fa.gz"
    output:
        f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{SPEC_NAME}_chroms.txt"
    params:
        CHROM_START_CHR=CHROM_START_CHR,
        NUM_ALL_CHR=NUM_ALL_CHR
    shell:
        """
        zcat < {input} | grep '>{params.CHROM_START_CHR}' | head -n {params.NUM_ALL_CHR} | sed 's/^>//' > {output}
        """

rule FILTER_PAF_CHR_ONLY:
    input:
        PAF=f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{ALN_NAME}.chain.fltr.srt.paf", 
        ALL_CHROMS=f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{SPEC_NAME}_chroms.txt"
    output:
        f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{ALN_NAME}.chain.chr.fltr.srt.paf"
    shell:
        """
        awk 'BEGIN {{ while (getline < "{input.ALL_CHROMS}") list[$0] }} $6 in list' {input.PAF} > {output}
        """

rule ALNPLOT:
    input:
        f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{ALN_NAME}.chain.1aln"
    output:
        f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{ALN_NAME}.chain.pdf"
    params:
        ALN_NAME=ALN_NAME, 
        CLADE=CLADE, 
        FULL_SPEC_NAME=FULL_SPEC_NAME
    shell:
        """
        ALNplot -p -H500 {input}
        mv {params.ALN_NAME}.chain.pdf ../groups/{params.CLADE}/{params.FULL_SPEC_NAME}/heterozygosity/{params.ALN_NAME}.chain.pdf
        """
#### ---- END ---- ####

#### ---------- RULES TO GET CHROMOSOME INFORMATION ---------- ####
rule CHROM_LENGTH_CALC:
    input:
        FASTA=expand("../groups/{CLADE}/{FULL_SPEC_NAME}/{ALN_NAME}.fa.gz", CLADE=CLADE, FULL_SPEC_NAME=FULL_SPEC_NAME, ALN_NAME=ALN_NAME), 
        ALL_CHROMS=expand("../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{SPEC_NAME}_chroms.txt", CLADE=CLADE, FULL_SPEC_NAME=FULL_SPEC_NAME, SPEC_NAME=SPEC_NAME),
    output:
        f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{SPEC_NAME}_Chroms_Lengths.txt"
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
        f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{ALN_NAME}.chain.chr.fltr.srt.paf"
    output:
        f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/temp/{SPEC_NAME}_Var_Only.txt"
    shell:
        """
        k8 paftools.js call {input} | awk '$1 == "V"' > {output}
        """

rule GET_WHOLE_ALN:
    input:
        f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{ALN_NAME}.chain.chr.fltr.srt.paf"
    output:
        f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/temp/{SPEC_NAME}_Aln_Only.txt"
    shell:
        """
        k8 paftools.js call {input} | awk '$1 == "R"' > {output}
        """

rule GET_VAR_ONLY_PER_CHROM:
    input:
        var_file=f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/temp/{SPEC_NAME}_Var_Only.txt", 
        paf_file=f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{ALN_NAME}.chain.chr.fltr.srt.paf"
    output:
        "../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/temp/{CHR}_Var_Only.txt"
    shell:
        """
        awk '$1 == "V" && $2 == "{wildcards.CHR}"' {input.var_file} > {output}
        if [ "$(stat -c %s {output})" -eq 0 ]; then
            k8 paftools.js call -L10000 {input.paf_file} | awk '$1 == "V" && $2 == "{wildcards.CHR}"' > {output}
        fi
        """

rule GET_ALN_ONLY_PER_CHROM:
    input:
        aln_file=f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/temp/{SPEC_NAME}_Aln_Only.txt", 
        paf_file=f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{ALN_NAME}.chain.chr.fltr.srt.paf"
    output:
        "../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/temp/{CHR}_Aln_Only.txt"
    shell:
        """
        awk '$1 == "R" && $2 == "{wildcards.CHR}"' {input} > {output}
        if [ "$(stat -c %s {output})" -eq 0 ]; then
            k8 paftools.js call -L10000 {input.paf_file} | awk '$1 == "R" && $2 == "{wildcards.CHR}"' > {output}
        fi
        """

#### ---- END ---- ####

####---------- RULES FOR CALCULATING ROH ----------####
rule ROH_CALC:
    input:
        chr_aln_only="../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/temp/{CHR}_Aln_Only.txt", 
        chr_var_only="../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/temp/{CHR}_Var_Only.txt", 
        CHROM_LENGTH_FILE=f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{SPEC_NAME}_Chroms_Lengths.txt"
    output:
        ROH_outfiles="../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/ROH/{CHR}_ROH_Results.txt"
    params:  
        CLADE=CLADE, 
        SPEC_NAME=SPEC_NAME
    resources:
        mem_mb=100000
    shell:
        """
        mkdir -p ../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/ROH/
        python {ROH_CALC_PY} {wildcards.CHR} {input.CHROM_LENGTH_FILE} {params.CLADE} {input.chr_aln_only} {input.chr_var_only} {params.SPEC_NAME} {output}
        """

rule COMPILE_ROH:
    input:
        expand("../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/ROH/{CHR}_ROH_Results.txt", CLADE=CLADE, FULL_SPEC_NAME=FULL_SPEC_NAME, CHR=CHROMS)
    output:
        f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/ROH/{SPEC_NAME}_ROH_Results.csv"
    shell:
        """
        awk -v var='chrom' -F',' '{{if ($1!=var) {{ print $0 }}}}' {input} >> {output}
        """
#### ---- END ---- ####

####---------- RULES FOR PLOTTING ROH ----------####
rule Plot_ROH:
    input:
        ROH="../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/ROH/{SPEC_NAME}_ROH_Results.csv",
        CHROM_LENGTH_FILE="../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{SPEC_NAME}_Chroms_Lengths.txt", 
        ALN_FILE="../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/temp/{SPEC_NAME}_Aln_Only.txt",
    output:
        "../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{SPEC_NAME}_ROH_Map.png"
    params:
        CLADE=CLADE,
        NUM_AUT_CHROMOSOMES=NUM_AUT_CHROMOSOMES, 
        SPEC_NAME=SPEC_NAME, 
        NUM_ALL_CHR=NUM_ALL_CHR
    shell:
        """
        Rscript {ROH_PLOT_R} {input.ROH} {input.CHROM_LENGTH_FILE} {input.ALN_FILE} {params.CLADE} {params.NUM_AUT_CHROMOSOMES} {params.SPEC_NAME} {params.NUM_ALL_CHR} {output}
        """
#### ---- END ---- ####

####---------- RULES FOR CALCULATING FROH ----------####
rule FROH_PER_AUT_CHR:
    input:
        ALN_FILE="../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/temp/{CHROM}_Aln_Only.txt",
        ROH=expand("../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/ROH/{SPEC_NAME}_ROH_Results.csv", CLADE=CLADE, FULL_SPEC_NAME=FULL_SPEC_NAME, SPEC_NAME=SPEC_NAME),
        CHROM_LENGTH_FILE=expand("../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{SPEC_NAME}_Chroms_Lengths.txt", CLADE=CLADE, FULL_SPEC_NAME=FULL_SPEC_NAME, SPEC_NAME=SPEC_NAME)
    output:
        "../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/ROH/{CHROM}_FROH_results.txt"
    params:
        CLADE=CLADE, 
        SPEC_NAME=SPEC_NAME, 
    resources:
        mem_mb=100000
    shell:
        """
        python {FROH_CALC_PER_CHR} {input.ALN_FILE} {input.ROH} {params.CLADE} {params.SPEC_NAME} {wildcards.CHROM} {input.CHROM_LENGTH_FILE} {output}
        """

rule WHOLE_FROH:
    input:
        ALN_FILE="../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/temp/{SPEC_NAME}_Aln_Only.txt",
        ROH="../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/ROH/{SPEC_NAME}_ROH_Results.csv", 
        CHROM_LENGTH_FILE="../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{SPEC_NAME}_Chroms_Lengths.txt", 
    output:
        "../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{SPEC_NAME}_FROH.txt"
    params:
        NUM_AUT_CHROMOSOMES=NUM_AUT_CHROMOSOMES, 
    resources:
       mem_mb=100000
    shell:
        """
        python {FROH_CALC_WHOLE} {input.ALN_FILE} {params.NUM_AUT_CHROMOSOMES} {input.CHROM_LENGTH_FILE} {input.ROH} {output}
        """
####---------- END ----------####

####---------- RULES FOR CALCULATING HETEROZYGOSITY ----------####
rule CALC_HET_PER_CHR:
    input:
        VAR_FILE="../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/temp/{CHROM}_Var_Only.txt",
        CHROM_LENGTH_FILE=expand("../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{SPEC_NAME}_Chroms_Lengths.txt", CLADE=CLADE, FULL_SPEC_NAME=FULL_SPEC_NAME, SPEC_NAME=SPEC_NAME), 
        ROH=expand("../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/ROH/{SPEC_NAME}_ROH_Results.csv", CLADE=CLADE, FULL_SPEC_NAME=FULL_SPEC_NAME, SPEC_NAME=SPEC_NAME)
    output:
        "../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/Het/{CHROM}_het.txt"
    params:
        CLADE=CLADE, 
        NUM_AUT_CHROMOSOMES=NUM_AUT_CHROMOSOMES, 
        SPEC_NAME=SPEC_NAME, 
        FULL_SPEC_NAME=FULL_SPEC_NAME,
        WINDOW_INTERVAL=WINDOW_INTERVAL,  
        WINDOW_LENGTH=WINDOW_LENGTH
    shell:
        """
        mkdir -p ../groups/{wildcards.CLADE}/{wildcards.FULL_SPEC_NAME}/heterozygosity/Het
        python {CALC_HET_PER_CHR_PY} {input.CHROM_LENGTH_FILE} {params.CLADE} {params.FULL_SPEC_NAME} {params.WINDOW_LENGTH} {params.WINDOW_INTERVAL} {input.ROH} {wildcards.CHROM} {output}
        """

rule CALC_HET_WHOLE_GENOME:
    input:
        CHROM_LENGTH_FILE="../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{SPEC_NAME}_Chroms_Lengths.txt",
        PER_CHR_FILES=expand("../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/Het/{CHROM}_het.txt", CLADE=CLADE, FULL_SPEC_NAME=FULL_SPEC_NAME, CHROM=AUTO_CHROMS)
    output:
        WHOLE_PER_CHR="../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{SPEC_NAME}_per_chr_mean_heterozygosity.txt",
        WHOLE_GENOME="../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{SPEC_NAME}_whole_genome_mean_heterozygosity.txt", 
    params:
        CLADE=CLADE, 
        SPEC_NAME=SPEC_NAME,
        NUM_AUT_CHROMOSOMES=NUM_AUT_CHROMOSOMES, 
        FULL_SPEC_NAME=FULL_SPEC_NAME
    shell:
        """
        python {CALC_HET_WHOLE_GENOME_PY} {input.CHROM_LENGTH_FILE} {params.CLADE} {params.SPEC_NAME} {params.FULL_SPEC_NAME} {params.NUM_AUT_CHROMOSOMES} {output.WHOLE_PER_CHR} {output.WHOLE_GENOME}
        """

rule COMPILE_HET:
    input:
        expand("../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/Het/{CHROM}_het.txt", CHROM=AUTO_CHROMS, CLADE=CLADE, FULL_SPEC_NAME=FULL_SPEC_NAME)
    output:
        "../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/Het/{SPEC_NAME}_Het_Compiled.tsv"
    shell:
        """
        awk -v var='chrom' -F',' '{{if ($1!=var) {{ print $1, $2, $3, $4, $5, $6, $7, $8, $9 }}}}' {input} >> {output}
        """
####---------- END ----------####

####---------- RULES FOR PLOTTING HETEROZYGOSITY ----------####
rule PLOT_WHOLE_HET:
    input:
        "../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/Het/{SPEC_NAME}_Het_Compiled.tsv"
    output:
        "../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{SPEC_NAME}_Het_Whole_Genome_Map.png"
    params:
        CLADE=CLADE, 
        NUM_AUT_CHROMOSOMES=NUM_AUT_CHROMOSOMES, 
        SPEC_NAME=SPEC_NAME, 
    shell:
        """
        Rscript {PLOT_WHOLE_HET_R} {input} {params.CLADE} {params.NUM_AUT_CHROMOSOMES} {params.SPEC_NAME} {output}
        """
####---------- END ----------####

####---------- RULES FOR MSMC DATA PREP ----------####
rule SEP_PAF_BY_CHR:
    input:
        expand("../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{ALN_NAME}.chain.chr.fltr.srt.paf", CLADE=CLADE, FULL_SPEC_NAME=FULL_SPEC_NAME, ALN_NAME=ALN_NAME)
    output:
        "../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/paf_files/{SPEC_NAME}_{CHROM}.paf"
    shell:
        """
        mkdir -p ../groups/{wildcards.CLADE}/{wildcards.FULL_SPEC_NAME}/heterozygosity/MSMC/paf_files
        awk -v var={wildcards.CHROM} '{{if($6 == var) print}}' {input} > {output}
        """

rule GEN_CHROM_VCF:
    input:
        paf="../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/paf_files/{SPEC_NAME}_{CHROM}.paf",
        refseq=expand("../groups/{CLADE}/{FULL_SPEC_NAME}/{ALN_NAME}.fa.gz", CLADE=CLADE, FULL_SPEC_NAME=FULL_SPEC_NAME, ALN_NAME=ALN_NAME), 
    output:
        "../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/vcf_files/{SPEC_NAME}_{CHROM}.vcf.gz"
    resources:
        mem_mb=50000
    shell:
        """
        mkdir -p "../groups/{wildcards.CLADE}/{wildcards.FULL_SPEC_NAME}/heterozygosity/MSMC/vcf_files"
        k8 paftools.js call -s {wildcards.SPEC_NAME} -f {input.refseq} {input.paf} > ../groups/{wildcards.CLADE}/{wildcards.FULL_SPEC_NAME}/heterozygosity/MSMC/vcf_files/{wildcards.SPEC_NAME}_{wildcards.CHROM}.vcf
        gzip ../groups/{wildcards.CLADE}/{wildcards.FULL_SPEC_NAME}/heterozygosity/MSMC/vcf_files/{wildcards.SPEC_NAME}_{wildcards.CHROM}.vcf
        """

rule FILTER_VCF:
    input:
        "../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/vcf_files/{SPEC_NAME}_{CHROM}.vcf.gz"
    output:
        "../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/filtered_vcf_files/{SPEC_NAME}_{CHROM}_filtered.vcf.gz"
    shell:
        """
        mkdir -p ../groups/{wildcards.CLADE}/{wildcards.FULL_SPEC_NAME}/heterozygosity/MSMC/filtered_vcf_files
        bcftools view -V indels {input} -Oz -o ../groups/{wildcards.CLADE}/{wildcards.FULL_SPEC_NAME}/heterozygosity/MSMC/vcf_files/{wildcards.CHROM}_no_indels.vcf.gz
        zcat ../groups/{wildcards.CLADE}/{wildcards.FULL_SPEC_NAME}/heterozygosity/MSMC/vcf_files/{wildcards.CHROM}_no_indels.vcf.gz | sed  's/1\/1/0\/1/g' | bgzip -c > {output}
        """

rule GEN_ROH_NEGATIVE_MASK:
    input:
        "../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/ROH/{CHROM}_ROH_Results.txt"
    output:
        "../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/ROH_Negative_Mask/{SPEC_NAME}_{CHROM}_ROH_Mask.bed.gz"
    shell:
        """
        mkdir -p ../groups/{wildcards.CLADE}/{wildcards.FULL_SPEC_NAME}/heterozygosity/MSMC/ROH_Negative_Mask
        python {ROH_MASKER} {wildcards.CLADE} {wildcards.SPEC_NAME} {wildcards.CHROM} {input} ../groups/{wildcards.CLADE}/{wildcards.FULL_SPEC_NAME}/heterozygosity/MSMC/ROH_Negative_Mask/{wildcards.SPEC_NAME}_{wildcards.CHROM}_ROH_Mask.bed
        gzip ../groups/{wildcards.CLADE}/{wildcards.FULL_SPEC_NAME}/heterozygosity/MSMC/ROH_Negative_Mask/{wildcards.SPEC_NAME}_{wildcards.CHROM}_ROH_Mask.bed
        """

rule MAKE_MASK:
    input:
        ALN_FILE="../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/temp/{CHROM}_Aln_Only.txt", 
        CHROM_LENGTHS=expand("../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/{SPEC_NAME}_Chroms_Lengths.txt", CLADE=CLADE, FULL_SPEC_NAME=FULL_SPEC_NAME, SPEC_NAME=SPEC_NAME)
    output:
        "../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/Aln_Mask/{CHROM}_Aln_Mask.bed"
    resources:
        mem_mb=100000
    shell:
        """
        mkdir -p ../groups/{wildcards.CLADE}/{wildcards.FULL_SPEC_NAME}/heterozygosity/MSMC/Aln_Mask
        python {ALN_MASK_MAKER_PY} {wildcards.CHROM} {input.CHROM_LENGTHS} {wildcards.CLADE} {input.ALN_FILE} {output}
        """

#expand( , CLADE=CLADE, FULL_SPEC_NAME=FULL_SPEC_NAME, SPEC_NAME=SPEC_NAME, CHROM=AUTO_CHROMS), 
rule MAIN_MULTIHETSEP:
    input:
        VCF="../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/filtered_vcf_files/{SPEC_NAME}_{CHROM}_filtered.vcf.gz", 
        ROH_negative_mask="../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/ROH_Negative_Mask/{SPEC_NAME}_{CHROM}_ROH_Mask.bed.gz", 
        MASK="../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/Aln_Mask/{CHROM}_Aln_Mask.bed"
    output:
        "../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/Output_primary_multihetsep/{SPEC_NAME}_{CHROM}_multihet_new.txt"
    shell:
        """
        mkdir -p ../groups/{wildcards.CLADE}/{wildcards.FULL_SPEC_NAME}/heterozygosity/MSMC/Output_primary_multihetsep
        if [ "$(stat -c %s {input.ROH_negative_mask})" -eq 51 ] || [ "$(stat -c %s {input.ROH_negative_mask})" -eq 52 ]; then
            python {GENERATE_MULTIHETSEP} {input.VCF} --mask={input.MASK} > {output}
        else
            python {GENERATE_MULTIHETSEP} {input.VCF} --mask={input.MASK} --negative_mask={input.ROH_negative_mask} > {output}
        fi
        """
####---------- END ----------####

####---------- RULES FOR RUNNING AND PLOTTING MSMC ----------####
## Define function to get only multihetsep files with 1k+ variants
def get_large_multihetsep(wildcards):
    files = expand(
            "../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/Output_primary_multihetsep/{SPEC_NAME}_{CHROM}_multihet_new.txt",
            FULL_SPEC_NAME=wildcards.FULL_SPEC_NAME, 
            CLADE=wildcards.CLADE,
            SPEC_NAME=wildcards.SPEC_NAME,
            CHROM=AUTO_CHROMS
        )
    large_files = []
    for f in files:
        num_lines = sum(1 for _ in open(f))
        if num_lines >= 1000:
            large_files.append(f)
    return large_files

rule RUN_PRIMARY_MSMC:
    input:
        get_large_multihetsep
    output:
        loop="../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/Primary_Results/{SPEC_NAME}.msmc2.loop.txt", 
        log="../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/Primary_Results/{SPEC_NAME}.msmc2.log", 
        final="../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/Primary_Results/{SPEC_NAME}.msmc2.final.txt"
    resources:
        runtime="1.25h"
    shell:
        """
        mkdir -p ../groups/{wildcards.CLADE}/{wildcards.FULL_SPEC_NAME}/heterozygosity/MSMC/Primary_Results
        msmc2/build/release/msmc2 -t 12 --fixedRecombination -o ../groups/{wildcards.CLADE}/{wildcards.FULL_SPEC_NAME}/heterozygosity/MSMC/Primary_Results/{wildcards.SPEC_NAME}.msmc2 {input}
        """

rule PLOT_MSMC:
    input:
        "../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/Primary_Results/{SPEC_NAME}.msmc2.final.txt"
    output:
        "../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/{SPEC_NAME}_MSMC2.png"
    params:
        CLADE=CLADE, 
        SPEC_NAME=SPEC_NAME, 
        MU=MU, 
        GEN_TIME=GEN_TIME
    shell:
        """
        Rscript {PLOT_MSMC_R} {params.CLADE} {params.SPEC_NAME} {params.MU} {params.GEN_TIME} {output} {input}
        """
####---------- END ----------####

#### ---------- RULES FOR BOOTSTRAPPING MSMC ---------- ####
rule BOOTSTRAPPING_MULTIHET_FILE:
    input:
        "../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/Output_primary_multihetsep/{SPEC_NAME}_{CHROM}_multihet_new.txt"
    output:
        "../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/Output_Multihet_Bootstrapped/{BOOT}/{SPEC_NAME}_{CHROM}_multihet.txt"
    shell:
        """
        mkdir -p ../groups/{wildcards.CLADE}/{wildcards.FULL_SPEC_NAME}/heterozygosity/MSMC/Output_Multihet_Bootstrapped/{wildcards.BOOT}
        python {BOOTSTRAPPING_GENERATOR} --inmhs {input} --windowsize 1e+05 --outmhs ../groups/{wildcards.CLADE}/{wildcards.FULL_SPEC_NAME}/heterozygosity/MSMC/Output_Multihet_Bootstrapped/{wildcards.BOOT}/{wildcards.SPEC_NAME}_{wildcards.CHROM}_multihet.txt
        """
#### ---------- END ---------- ####

####---------- RULES FOR RUNNING AND PLOTTING BOOTSTRAPPED MSMC ----------####
rule RUN_BOOTSTRAPPING_MSMC:
    input:
        lambda wildcards: expand(
            "../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/Output_Multihet_Bootstrapped/{BOOT}/{SPEC_NAME}_{CHROM}_multihet.txt",
            FULL_SPEC_NAME=wildcards.FULL_SPEC_NAME, 
            CLADE=wildcards.CLADE,
            SPEC_NAME=wildcards.SPEC_NAME,
            BOOT=wildcards.BOOT,
            CHROM=AUTO_CHROMS
        )
    output:
        loop="../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/Bootstrap_results/{SPEC_NAME}_Bootstrapping_{BOOT}.msmc2.loop.txt", 
        log="../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/Bootstrap_results/{SPEC_NAME}_Bootstrapping_{BOOT}.msmc2.log",
        final="../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/Bootstrap_results/{SPEC_NAME}_Bootstrapping_{BOOT}.msmc2.final.txt"
    resources:
        runtime="1h", 
        mem_mb=100000
    shell:
        """
        mkdir -p ../groups/{wildcards.CLADE}/{wildcards.FULL_SPEC_NAME}/heterozygosity/MSMC/Bootstrap_results
        msmc2/build/release/msmc2 -t 12 --fixedRecombination -o ../groups/{wildcards.CLADE}/{wildcards.FULL_SPEC_NAME}/heterozygosity/MSMC/Bootstrap_results/{wildcards.SPEC_NAME}_Bootstrapping_{wildcards.BOOT}.msmc2 {input}
        """

rule PLOT_MSMC_BOOSTRAP:
    input:
        f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/Primary_Results/{SPEC_NAME}.msmc2.final.txt"
    output:
        f"../groups/{CLADE}/{FULL_SPEC_NAME}/heterozygosity/MSMC/{SPEC_NAME}_MSMC2_Bootstrapped.png"
    params:
        CLADE=CLADE, 
        SPEC_NAME=SPEC_NAME, 
        MU=MU, 
        GEN_TIME=GEN_TIME
    shell:
        """
        Rscript {PLOT_MSMC_BOOSTRAP_R} {params.CLADE} {params.SPEC_NAME} {params.MU} {params.GEN_TIME} {output} {input}
        """
####---------- END ----------####
