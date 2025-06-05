README

####
20241231 (December 31st, 2024)

Amanda Gardiner made this README file within this directory for doing ROH and heterozygosity analyses on VGP reference and alternate genomes
Made directories for storing results of analysis for each taxon, with directories based on groupings in VGP/clades directory

Created 20241231_ROH_Calc.sh script, which will calculate ROH for each chromosome for each organism


#### UPDATE ####
20250101 (January 1st, 2025)

Created 20241231_ROH_Durbin_Calc_Eqns.R script to be called by 20241231_ROH_Calc.sh to calculate ROH
Finished writing 20241231_ROH_Calc.sh

In VGP/20241117.UCSC-hubs-BGP-alignment/reference/invertebrate:
    zcat < GCA_035083965.1.fa.gz
    zcat < GCA_035083965.1.fa.gz | grep '>CM'
    zcat < GCA_035083965.1.fa.gz | grep '>CM' | wc --> 19

Created other/chroms_lists/GCA_035083965.1_chroms.txt file
Pasted chromosome names for GCA_035083965.1.fa.gz into this file and removed '>'
    sbatch 20241231_ROH_Calc.sh 
Submitted batch job 2751863

Checked initial results for first runs in the array, and it was not working, so I canceled the rest of the jobs
    scancel 2751863
Error:
    Error in read.table(args[6]) : no lines available in input
    Execution halted

Realized (very silly) that I need to align the fasta files to generate a vcf file for each species before I can run ROH analysis
Created 20250101_mm2_alignment.sh file
This will create a species directory within /heterozygosity for each species and align the reference and alternate fasta files to create a vcf file
Ran with:
    ../241117.UCSC-hubs-VGP-alignment/alignment/reference/invertebrate/GCA_035083965.1.fa.gz
    ../241117.UCSC-hubs-VGP-alignment/alignment/alternate/echinoderm/GCA_035149815.1.fa.gz 
Submitted batch job 2751921

Canceled job partway through --> path file for alternate genome was incorrect
Fixed:
    ../241117.UCSC-hubs-VGP-alignment/alignment/alternate/invertebrate/GCA_035149815.1.fa.gz
Also had to fix paftools syntax given that command was not recognized
Added paftools.js file
Confirmed commands can be called
    k8 paftools.js
Added that command to script 
Submitted batch job 2751974


#### UPDATE ####
20250102 (January 2nd, 2025)

Job 2751974 completed
Looking at the error files, it appears that the code worked!
Output files went into heterozygosity directory, moved them into the other/GCA_035083965.1 directory

Updated 20250101_mm2_alignment.sh script to have variable for clade directory (e.g. other, sharks) and have the alignment output files put into the clade/GC... directory
Couldn't find matching alternate genomes for otherChordates or other invertebrates, so moved onto sharks directory
Ran for Hemiscyllium ocellatum (Epaulette shark):
    ../241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/GCF_020745735.1.fa.gz
    ../241117.UCSC-hubs-VGP-alignment/alignment/alternate/sharks/GCA_020745765.1.fa.gz
Submitted batch job 2760190

Updated 20241231_ROH_Calc.sh 
Update includes refining directory variables and using vcf file
Submitted test to make sure everything was working
Submitted batch job 2763138


#### UPDATE ####
20250103 (January 3rd)

Job 2760190 to do minimap alingment for GCF_020745735.1 -- Hemiscyllium ocellatum (Epaulette shark) completed
It worked!

Updated 20250101_mm2_alignment.sh script for next species:
Ran for Mobula birostris (Giant manta):
    REFERENCE=../241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/GCA_030028105.1.fa.gz ## Finish filling in depending on genome and clade 
    ALTERNATE=../241117.UCSC-hubs-VGP-alignment/alignment/alternate/sharks/GCA_030035685.1.fa.gz ## Finish filling in depending on genome and clade
Submitted batch job 2814094

Job 2763138 to test ROH calculations for GCA_035083965.1 worked!
Ran for Branchiostoma lanceolatum (amphioxus):
    REF_NAME=GCA_035083965.1
    REFERENCE=../241117.UCSC-hubs-VGP-alignment/alignment/reference/invertebrate/${REF_NAME}.fa.gz ## Finish filling in depending on genome and clade
    ALTERNATE=../241117.UCSC-hubs-VGP-alignment/alignment/alternate/invertebrate/GCA_035149815.1.fa.gz ## Finish filling in depending on genome and clade
Uncommented out all commands in 20241231_ROH_Durbin_Calc_Eqns.R
Submitted batch job 2814115


#### UPDATE ####
20250106 (January 6th, 2025)

Job 2814115 to calculate ROH for GCA_035083965.1 (amphioxus) appears to have worked!

Created 20250106_Plot_ROH_FROH.sh to feed ROH.txt files into R scripts to plot ROH and calculate 20240106_Plot_ROH_FROH
Created 20250106_Plot_ROH.R to create pdf map of ROH on chromosomes for each individual
Created 20250106_FROH_Calc.R to calculate FROH (inbreeding coefficient for each individual)

Job 2814094 to do mm2 alignment for GCA_030028105.1 Mobula birostris (Giant manta) failed
Due to being out of memory
Upped memory to 30GB
Updated dates for out/err files to 20250106
Submitted batch job 3028223

Prepped 20250106_Plot_ROH_FROH.sh for GCA_035083965.1
Tested this command to get all chromosomes and their lengths for each reference species and save it as a txt file
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/GCA_030028105.1.fa.gz | awk '$0 ~ ">CM" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }'


#### UPDATE ####
20250108 (January 8th, 2025)

Finished writing 20250106_FROH_Calc.R and 20250106_Plot_ROH.R
Submitted 20240106_Plot_ROH_FROH.sh for GCA_035083965.1 (Branchiostoma lanceolatum (amphioxus))
Submitted batch job 3184218

Job 3028223 to do mm2 alignment for GCA_030028105.1 Mobula birostris (Giant manta) failed
Again it failed due to being out of memory
Upped memory to 120GB
Deleted GCA_030028105.1 directory within shark directory to remove corrupted alignment files
Submitted batch job 3185029

Set up 20241231_ROH_Calc.sh script to find ROH for GCF_020745735.1 Hemiscyllium ocellatum (Epaulette shark)
Ran commands to find chromosome names for Hem oce in 241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/
    zcat < GCF_020745735.1.fa.gz | grep '>'
    zcat < GCF_020745735.1.fa.gz | grep '>NC'
    zcat < GCF_020745735.1.fa.gz | grep '>NC' >../../../../heterozygosity/sharks/chrom_lists/GCF_020745735.1_chroms.txt
    zcat < GCF_020745735.1.fa.gz | grep '>NC' | wc ## This led to an output of 54 chromosomes
Manually went in and removed the ">" before "NC" in the GCF_020745735.1_chroms.txt file
Updated dates for output files in this script to 20250108
Submitted batch job 3188436

Job 3184218 to plot ROH and calculate FROH failed
    /var/spool/slurm/slurmd/job3184218/slurm_script: line 32: ../241117.UCSC-hubs-VGP-alignment/alignment/reference/other/GCA_035083965.1.fa.gz: No such file or directory
    awk: fatal: cannot open file `20250103_other/chrom_lists/GCA_035083965.1_chroms.txt_100kb_Results_ROH_Durbin_Calc.txt' for reading (No such file or directory)
    Error in args[9] : object of type 'closure' is not subsettable
    Execution halted
Problem in first error message was that the clade directory was 'invertebrate' instead of 'other' --> fixed in 20250106_Plot_ROH_FROH.sh script
For awk error in second line, fixed by changing for statement line to:
    for i in ${cat ${ALL_CHROMS}}
Instead of 
     for i in ${ALL_CHROMS[@]}
Changed TODAY_DATE variable and dates for out/err files
Submitted batch job 3203285


#### UPDATE ####
20250110 (January 10th, 2025)

Job 3203285 to plot ROH and calculate FROH for GCA_035083965.1 Branchiostoma lanceolatum (amphioxus) failed
Errors:
    /var/spool/slurm/slurmd/job3203285/slurm_script: line 39: ${cat ${ALL_CHROMS}}: bad substitution
    Error in args[9] : object of type 'closure' is not subsettable
Fixed the first error by changing loop syntax from if statement to:
    if [ -f ${SPECIES}/${TODAY_DATE}_${REF_NAME}_ROH_Results.tsv ] 
    then
        echo 'ROH summation file exists'
    else 
        while read -r line; do 
            awk -v var=${line} 'NR!=1 && !/[a-zA-Z]/ { print var, $2, $3, $4 }' ${ROH_CALC_DATE}_${line}_100kb_Results_ROH_Durbin_Calc.txt >> ${SPECIES}/${TODAY_DATE}_${REF_NAME}_ROH_Results.tsv            
        done < "${ALL_CHROMS}"
    fi
Realized second error was from TODAY_DATE and REF_NAME not being included as input variables for the 20250106_Plot_ROH.R script
Added them as inputs in the .sh script
Updated TODAY_DATE and date for .out/.err files
Submitted batch job 3308959

Job 3188436 to find ROH for GCF_020745735.1 Hemiscyllium ocellatum (Epaulette shark) completed successfully
However looking at results I realized that the ROH didn't separate by chromosome --> this is because I am doing it with the whole vcf file rather than subsetting it by chromosome
Went back to look at ROH txt files for GCA_035083965.1 Branchiostoma lanceolatum (amphioxus) and realized the same problem is present
Removed the .fa and ROH txt files for both species

Job 3308959 failed
Went back to 20241231_ROH_Calc.sh and altered script so that the vcf file for the whole genome is subset by chromosome and ROH are calculated just for a single chromosome per txt file
    ## Subset VCF file to get one for just that chromosome 
    vcftools --vcf ${VCF} --chr ${CHROM} --out ${SPECIES}/temp/${TODAY_DATE}_${CHROM}.vcf
Added command to remove CHROM vcf file at end of script to save storage space
    ## Remove temp files to save storage space
    rm -r -f ${SPECIES}/temp/${TODAY_DATE}_${CHROM}.vcf
Removed command to subset fasta script by chromosome since it wasn't necessary
Prepped script for GCA_035083965.1 Branchiostoma lanceolatum (amphioxus)
Submitted batch job 3311282

Job completed but didn't work
chromosome vcf files were not generated - as such txt files were empty
I think the reason the files were not generated is because I didn't include --recode in the vcftools command
Updated:
    vcftools --vcf ${VCF} --chr ${CHROM} --recode --out ${SPECIES}/temp/${TODAY_DATE}_${CHROM}.vcf
For GCA_035083965.1 Branchiostoma lanceolatum (amphioxus)
Submitted batch job 3312849

Job 3312849 completed
chromosome vcf files were generated but with the file names ending in .vcf.recode.vcf
Manually fixed using mv command
Fixed script so that file names owuld end in recode.vcf, and the Rscript and rm commands would look for files ending in .recode.vcf
Added if statement to 20241231_ROH_Calc.sh to avoid wasting time rerunning the vcftools command
    if [ -f ${SPECIES}/temp/${TODAY_DATE}_${CHROM}.recode.vcf ] 
        then
            echo 'Chromosome vcf file exists'
        else 
            vcftools --vcf ${VCF} --chr ${CHROM} --recode --out ${SPECIES}/temp/${TODAY_DATE}_${CHROM}
        done
    fi
    # vcftools --vcf ${VCF} --chr ${CHROM} --recode --out ${SPECIES}/temp/${TODAY_DATE}_${CHROM}
Submitted batch job 3313718

Job 3313718 failed due to syntax error 
    /var/spool/slurm/slurmd/job3314052/slurm_script: line 50: syntax error near unexpected token `done'
    /var/spool/slurm/slurmd/job3314052/slurm_script: line 50: `    done'
Updated if statement to fix this error
    if [ -f ${SPECIES}/temp/${TODAY_DATE}_${CHROM}.recode.vcf ] 
        then
            echo 'Chromosome vcf file exists'
        else 
            vcftools --vcf ${VCF} --chr ${CHROM} --recode --out ${SPECIES}/temp/${TODAY_DATE}_${CHROM}
    fi
Removed empty txt files that were generated
Submitted batch job 3314564

It worked!
Prepped 20241231_ROH_Calc.sh for GCF_020745735.1 Hemiscyllium ocellatum (Epaulette shark)
    #SBATCH --array=1-54 (based on number of chromosomes)
    CLADE=sharks

    ## Reference genome name
    REF_NAME=GCF_020745735.1
Submitted batch job 3315207

Job 3315207 completed!
It worked!
Prepped 20241231_ROH_Calc.sh for GCA_030028105.1 Mobula birostris (Giant manta)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/GCA_030028105.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/GCA_030028105.1.fa.gz | grep '>CM' >sharks/chrom_lists/GCA_030028105.1_chroms.txt
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/GCA_030028105.1.fa.gz | grep '>CM' | wc ## 33 chromosomes
Submitted batch job 3320060

Resubmitting 20240106_Plot_ROH_FROH.sh for GCA_035083965.1 Branchiostoma lanceolatum (amphioxus) now that ROH txt files have been fixed
Submitted batch job 3320175

Prepping 202250101_mm2_alignment.sh script for GCF_030144855.1 Hypanus sabinus (Atlantic stingray)
Submitted batch job 3320367


#### UPDATE ####
20250111 (January 11th, 2025)

Job 3320175 for 20240106_Plot_ROH_FROH.sh for GCA_035083965.1 Branchiostoma lanceolatum (amphioxus) failed
Errors:
    awk: fatal: cannot open file `20250103_CM068865.1_100kb_Results_ROH_Durbin_Calc.txt' for reading (No such file or directory) ## This error for every single chromosome
    Error in args[9] : object of type 'closure' is not subsettable
First, removed generated files from this run since they were empty and could interfere with next run
    rm -f -r 20250110_GCA_035083965.1_ROH_Results.tsv
    rm -f -r 20250110_GCA_035083965.1_ROH_Map.pdf
    rm -f -r 20250110_GCA_035083965.1_FROH.txt
First error (awk: fatal) due to date for ROH txt files being wrong (had in code as 20250103 whereas it is actually 20250110)
Fixed in 20250106_Plot_ROH_FROH.sh script 
    ROH_CALC_DATE=20250110
Second error is coming from args[9], potentially in both R scripts
Realized errors was coming from 20250106_FROH_Calc.R script --> had not defined the args variable at start of script but further down
Fixed by moving to top of script:
    args <- commandArgs()
Updated TODAY_DATE in 20240106_Plot_ROH_FROH.sh
Submitted batch job 3415685

Job 3415685 failed --> realized it was due to me not specifying the directory that these files are found in
Fixed in if statement by adding ${SPECIES}/:
    awk -v var=${line} 'NR!=1 && !/[a-zA-Z]/ { print var, $2, $3, $4 }' ${SPECIES}/${ROH_CALC_DATE}_${line}_100kb_Results_ROH_Durbin_Calc.txt >> ${SPECIES}/${TODAY_DATE}_${REF_NAME}_ROH_Results.tsv     
Removed empty generated files from this failed run:
    rm -f -r 20250111_GCA_035083965.1_ROH_Results.tsv
    rm -r -f 20250111_GCA_035083965.1_FROH.txt
    rm -f -r 20250111_GCA_035083965.1_ROH_Map.pdf
Submitted batch job 3415819

Job 3415819 failed
Errors:
    awk: fatal: cannot open file `other/GCA_035083965.1/20250110_CM068865.1_100kb_Results_ROH_Durbin_Calc.txt' for reading (No such file or directory)
    Error in `[.data.frame`(vcf, , .I[c(1, .N)], by = vcf[1, ]) : 
    unused argument (by = vcf[1, ])
    Calls: [
    Execution halted
    ...
    Error in `geom_segment()`:
    ! Problem while computing aesthetics.
    ℹ Error occurred in the 2nd layer.
    Caused by error:
    ! object 'Chr' not found
First error specifically for this single chromosome -- CM068865.1 -- no txt file is present for this chromosome
I manually verified there are files for the other chromosomes
There is an empty txt file called 20250110__100kb_Results_ROH_Durbin_Calc.txt -- possible this is the one for this chromosome, but not sure why the chromosome name was not attached
Will manually re-run the formation of the txt file for this one chromosome to see what happens
    vcftools --vcf other/GCA_035083965.1/GCA_035083965.1_aligned.mm2.vcf --chr CM068865.1 --recode --out other/GCA_035083965.1/temp/20250110_CM068865.1
Successfully generated VCF file
    Rscript 20241231_ROH_Durbin_Calc_Eqns.R "other/GCA_035083965.1/temp/20250110_CM068865.1.recode.vcf" > other/GCA_035083965.1/20250110_CM068865.1_100kb_Results_ROH_Durbin_Calc.txt
Successfully generated txt file
Removed generated .tsv, .pdf, and .txt files from last run to avoid any corruption of files
    rm -f -r 20250111_GCA_035083965.1_ROH_Results.tsv
    rm -f -r 20250111_GCA_035083965.1_ROH_Map.pdf
    rm -f -r 20250111_GCA_035083965.1_FROH.txt
Resubmitted 20250106_Plot_ROH_FROH.sh for GCA_035083965.1 Branchiostoma lanceolatum (amphioxus)
Submitted batch job 3417164

Job 3417164 failed
Error:
    Error in `[.data.frame`(vcf, , .I[c(1, .N)], by = vcf[1, ]) : 
        unused argument (by = vcf[1, ])
    Calls: [
    Execution halted
    ...
    Error in `geom_segment()`:
    ! Problem while computing aesthetics.
    ℹ Error occurred in the 2nd layer.
    Caused by error:
    ! object 'Chr' not found
The first error is coming from 20250106_FROH_Calc.R script
    indx <- vcf[, .I[c(1, .N)], by = vcf[1,]]
Based on online recommendations - fixed with modification:
    indx <- vcf[, .I[c(1, .N)], by = vcf[1,]]$V1 ; vcf[indx]
The second error is coming from headers not being present in the .tsv file containing the ROH
Added lines into 20250106_Plot_ROH.R and 20250106_FROH_Calc.R scripts respectively to give column names
    colnames(dat) <- c('Chr', 'Start', 'End', 'Length')
    colnames(ROH_1mb) <- c('Chr', 'Start', 'End', 'Length')
Submitted batch job 3418274

Job 3418274 says it comleted successfull but same error popped up in .err file
    Error in `[.data.frame`(vcf, , .I[c(1, .N)], by = vcf[1, ]) : 
    unused argument (by = vcf[1, ])
    Calls: [
    Execution halted
Something went wrong with job, as .pdf file was only 26bytes and FROH txt file was empty
Removed .tsv, FROH txt, and pdf file
    rm -f -r 20250111_GCA_035083965.1_ROH_Results.tsv
    rm -f -r 20250111_GCA_035083965.1_ROH_Map.pdf
    rm -f -r 20250111_GCA_035083965.1_FROH.txt
Problem is with working with vcf file in R to get locations of first and last variant on each chromosome
I think I am going to try to get vcf file handled with information for FROH in shell script instead of r


#### UPDATE ####
20250112 (January 12th, 2025)

Finished editing 20250106_Plot_ROH_FROH.sh
Defined Laut variable with handling of vcf file in shell script
    while read -r line; do
        vcftools --vcf ${VCF} --chr ${line} --recode --out ${TEMP}/${TODAY_DATE}_${line}
        start="$(bcftools query -f '%POS' ${TEMP}/${TODAY_DATE}_${line}.recode.vcf  | head -n 1)"
        echo ${start}
        end="$(bcftools query -f '%POS' ${TEMP}/${TODAY_DATE}_${line}.recode.vcf   | tail -n 1)"
        echo ${end}
        chrom_length=$((${end}-${start}))
        echo ${chrom_length}
        Laut=$(($Laut + ${chrom_length}))
    done < "${ALL_CHROMS}"
Resubmitted for GCA_035083965.1 Branchiostoma lanceolatum (amphioxus)
Updated TODAY_DATE and date for .out/.err files
Submitted batch job 3446639

Job 3446639 completed but failed
Error:
    Error in L_ROH/Laut : non-numeric argument to binary operator
Changed input for Laut in 20250106_FROH_Calc.R script to:
    Laut <- is.numeric(args[11])
Reran
Submitted batch job 3446647

Job 3446647 completed but FROH_aut was inf -- not correct
Changed input for Laut in 20250106_FROH_Calc.R script to:
    Laut <- As.numeric(args[11])


#### UPDATE ####
20250113 (January 13th, 2025)

Resubmitted 20250106_Plot_ROH_FROH.sh for GCA_035083965.1 Branchiostoma lanceolatum (amphioxus)
Submitted batch job 3453300

Job 3453300 completed, with FROH successfully calculated
However PDF for maps for ROH is corrupted and cannot be opened
Updated TODAY_DATE and date in .out/.err files for 20250106_Plot_ROH_FROH.sh
Commented out dev.off() command and resubmitted
Submitted batch job 3453786

Removed old temp/...recode.vcf and .log files to save storage space
    rm 20250110_*
    rm: remove regular file '20250110_CM068865.1.log'? y
    rm: remove regular file '20250110_CM068865.1.recode.vcf'? y

    rm -f -r 20250111_CM0*
    rm -f -r 20250112_CM0*
Only temp/...recode.vcf files from 20250113 remain
Altered script to save computation time/memory so that I don't have to generated more temp/...recode.vcf files
    # while read -r line; do
    #     vcftools --vcf ${VCF} --chr ${line} --recode --out ${TEMP}/${TODAY_DATE}_${line}
    # done < "${ALL_CHROMS}"

    while read -r line; do
        start="$(bcftools query -f '%POS' ${TEMP}/${TODAY_DATE}_${line}.recode.vcf  | head -n 1)"
        end="$(bcftools query -f '%POS' ${TEMP}/${TODAY_DATE}_${line}.recode.vcf   | tail -n 1)"
        chrom_length=$((${end}-${start}))
        Laut=$(($Laut + ${chrom_length}))
    done < "${ALL_CHROMS}"
Commented out rm command to remove .recode.vcf files -- will keep files I generated today until I can get this working and then remove them to save on storage space

Job 3453786 completed successfully, but again the pdf file was not correctly created (0 byte size)
Manually entered R workspace and copied/pasted script from 20250106_Plot_ROH.R to run using files generated in last job
Successfully generated 20250113_GCA_035083965.1_ROH_Map.pdf

Added dev.off() command to 20250106_Plot_ROH.R script right before pdf and ggplot commands -- possibly could be source of issue when running code before

Job 3320060 to find ROH for GCA_030028105.1 Mobula birostris (Giant manta) completed but failed
These are the errors that kept appearing in all of the .err files for each chromosome:
    After filtering, kept 0 out of a possible 3714560 Sites
    No data left for analysis!
    ...
    In file(file, "rt") :
    cannot open file 'sharks/GCA_030028105.1/temp/20250110_.recode.vcf': No such file or directory
When checking in /sharks/GCA_030028105.1 directory --> the .txt files for ROH were generated but the names are strings contained in quotation marks and have a file size of 0
    e.g. '20250110_>CM057557.1_100kb_Results_ROH_Durbin_Calc.txt'
Removed these files with rm -f -r command
    rm -f -r 20250110_*
I checked sharks/chrom_lists/GCA_030028105.1_chroms.txt --> I think issue is the carat symbols before chromosome names so that the names are not exact matches when running vcftools
Manually removed ">" symbols from sharks/chrom_lists/GCA_030028105.1_chroms.txt
Will have to find a way to automate this in the future
Submitted batch job 3459127 for GCA_030028105.1 Mobula birostris (Giant manta)


Job 3320367 of 202250101_mm2_alignment.sh script for GCF_030144855.1 Hypanus sabinus (Atlantic stingray) still hasn't run
Canceled job
    scancel 3320367
Changed from cclake to icelake-himem node
Resubmitted using sbatch
Submitted batch job 3459408

Prepping chrom_lists for other species to do analysis on --> going through altnerate chromosomes in each group and comparing them to the reference genomes available
Prep for GCA_035084275.1 Hydrolagus colliei (Spotted ratfish)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/GCA_035084275.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/GCA_035084275.1.fa.gz | grep '>CM'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/GCA_035084275.1.fa.gz | grep '>CM' | wc ## --> 40 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/GCA_035084275.1.fa.gz | grep '>CM' | tr -d '>' > sharks/chrom_lists/GCA_035084275.1_chroms.txt
tr -d '>' command removes the carat symbol from the chromosome names to avoid issues during 20241231_ROH_Calc.sh

Prep for GCA_035084215.1 Heptranchias perlo (sharpnose sevengill shark)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/GCA_035084215.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/GCA_035084215.1.fa.gz | grep '>CM' | wc ## --> 47 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/GCA_035084215.1.fa.gz | grep '>CM' | tr -d '>' > sharks/chrom_lists/GCA_035084215.1_chroms.txt

Prep for GCA_036365495.1 Heterodontus francisci (horn shark)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/GCA_036365495.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/GCA_036365495.1.fa.gz | grep '>CM' | wc ## --> 51 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/GCA_036365495.1.fa.gz | grep '>CM' | tr -d '>' > sharks/chrom_lists/GCA_036365495.1_chroms.txt

Prep for GCA_036971175.1 Narcine bancroftii (Spotted torpedo ray/Caribbean electric ray)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/GCA_036971175.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/GCA_036971175.1.fa.gz | grep '>CM' | wc ## --> 14 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/GCA_036971175.1.fa.gz | grep '>CM' | tr -d '>' > sharks/chrom_lists/GCA_036971175.1_chroms.txt

Prep for GCF_025201925.1 Gopherus flavomarginatus (Mexican gopher tortoise) (alternate = GCA_025201965.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCF_025201925.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCF_025201925.1.fa.gz | grep '>NC' | wc ## --> 25 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCF_025201925.1.fa.gz | grep '>NC' | tr -d '>' > reptiles/chrom_lists/GCF_025201925.1_chroms.txt

NOTE: Could not find genome in reference directory for Gopherus evgoodei (Goodes thornscrub tortoise)
        Found in alternate directory: GCF_007399415.2
        According to NCBI the other haplotype is GCA_007399395.1 --> cannot find in reptile reference directory
        Ran this command and found nothing:
            find -name GCA_007399395.1*

Prep for GCF_030867095.1 Alligator mississippiensis (American alligator) (alternate = GCA_030867065.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCF_030867095.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCF_030867095.1.fa.gz | grep '>NC' | wc ## --> 17 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCF_030867095.1.fa.gz | grep '>NC' | tr -d '>' > reptiles/chrom_lists/GCF_030867095.1_chroms.txt

Prep for GCF_027887155.1 Malaclemys terrapin pileata (Diamondback terrapin) (alternate = GCA_027887205.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCF_027887155.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCF_027887155.1.fa.gz | grep '>NC' | wc ## --> 25 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCF_027887155.1.fa.gz | grep '>NC' | tr -d '>' > reptiles/chrom_lists/GCF_027887155.1_chroms.txt

Prep for GCA_028017835.1 Emys orbicularis (European pond turtle) (alternate = GCA_028017845.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCA_028017835.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCA_028017835.1.fa.gz | grep '>CM' | wc ## --> 25 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCA_028017835.1.fa.gz | grep '>CM' | tr -d '>' > reptiles/chrom_lists/GCA_028017835.1_chroms.txt

Prep for GCF_029931775.1.fa.gz Euleptes europaea (tarantolino) (alternate = GCA_029931755.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCF_029931775.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCF_029931775.1.fa.gz | grep '>NC' | wc ## --> 21 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCF_029931775.1.fa.gz | grep '>NC' | tr -d '>' > reptiles/chrom_lists/GCF_029931775.1_chroms.txt

Prep for GCA_030020295.1.fa.gz Gavialis gangeticus (Gharial) (altnerate = GCA_030020385.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCA_030020295.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCA_030020295.1.fa.gz | grep '>CM' | wc ## --> 16 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCA_030020295.1.fa.gz | grep '>CM' | tr -d '>' > reptiles/chrom_lists/GCA_030020295.1_chroms.txt

Prep for GCF_030035675.1.fa.gz Rhineura floridana (Florida worm lizard) (alternate = GCA_030035715.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCF_030035675.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCF_030035675.1.fa.gz | grep '>NC' | wc ## --> 22 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCF_030035675.1.fa.gz | grep '>NC' | tr -d '>' > reptiles/chrom_lists/GCF_030035675.1_chroms.txt

Prep for GCA_030412105.1.fa.gz Cyclura pinguis (alternate = GCA_030412085.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCA_030412105.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCA_030412105.1.fa.gz | grep '>CM' | wc ## --> 19 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCA_030412105.1.fa.gz | grep '>CM' | tr -d '>' > reptiles/chrom_lists/GCA_030412105.1_chroms.txt

Prep for GCA_030867105.1.fa.gz Liasis olivaceus (alternate = GCA_030867085.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCA_030867105.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCA_030867105.1.fa.gz | grep '>CM' | wc ## --> 17 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCA_030867105.1.fa.gz | grep '>CM' | tr -d '>' > reptiles/chrom_lists/GCA_030867105.1_chroms.txt

Prep for GCA_031021105.1.fa.gz Erythrolamprus reginae (royal ground snake) (alternates = GCA_031001685.1.fa.gz -- hap3, GCA_031021085.1.fa.gz -- hap2)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCA_031021105.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCA_031021105.1.fa.gz | grep '>CM' | wc ## --> 14 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCA_031021105.1.fa.gz | grep '>CM' | tr -d '>' > reptiles/chrom_lists/GCA_031021105.1_chroms.txt

Prep for GCA_033349115.1.fa.gz Macrochelys suwanniensis (Suwannee snapping turtle) (alternate = GCA_033296515.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCA_033349115.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCA_033349115.1.fa.gz | grep '>CM' | wc ## --> 26 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCA_033349115.1.fa.gz | grep '>CM' | tr -d '>' > reptiles/chrom_lists/GCA_033349115.1_chroms.txt

Prep for GCA_035046505.1.fa.gz Tiliqua scincoides (alternate = GCA_035046495.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCA_035046505.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCA_035046505.1.fa.gz | grep '>CM' | wc ## --> 16 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCA_035046505.1.fa.gz | grep '>CM' | tr -d '>' > reptiles/chrom_lists/GCA_035046505.1_chroms.txt

Prep for GCF_035149785.1.fa.gz Candoia aspera (alternate = GCA_035125265.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCF_035149785.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCF_035149785.1.fa.gz | grep '>NC' | wc ## --> 18 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCF_035149785.1.fa.gz | grep '>NC' | tr -d '>' > reptiles/chrom_lists/GCF_035149785.1_chroms.txt

Prep for GCA_037176765.1.fa.gz Anolis sagrei (Brown anole) (alternate = GCA_037176775.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCA_037176765.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCA_037176765.1.fa.gz | grep '>CM' | wc ## --> 15 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCA_037176765.1.fa.gz | grep '>CM' | tr -d '>' > reptiles/chrom_lists/GCA_037176765.1_chroms.txt

Finished getting chrom_lists for shark/invertebrate/reptile directories

Job 3459127 for GCA_030028105.1 Mobula birostris (Giant manta) completed but failed
Same error as before:
    After filtering, kept 0 out of a possible 3714560 Sites
    No data left for analysis!
    Run Time = 8.00 seconds
    Error in file(file, "rt") : cannot open the connection
    Calls: read.table -> file
    In addition: Warning message:
    In file(file, "rt") :
    cannot open file 'sharks/GCA_030028105.1/temp/20250113_.recode.vcf': No such file or directory
    Execution halted
When checking results, ROH .txt files were generated for each chromosome, with lists of ROH included in them
One empty ROH .txt file was generated
    20250113__100kb_Results_ROH_Durbin_Calc.txt
Given that ROH files were still generated, I'm going to manually go through process for one chromosome and cross check it to see if the results worked despite what the error messages saying
    Chr = CM057556.1
    vcftools --vcf sharks/GCA_030028105.1/GCA_030028105.1_aligned.mm2.vcf --chr CM057556.1 --recode --out sharks/GCA_030028105.1/temp/20250113_CM057556.1

        VCFtools - 0.1.17
    (C) Adam Auton and Anthony Marcketta 2009

    Parameters as interpreted:
            --vcf sharks/GCA_030028105.1/GCA_030028105.1_aligned.mm2.vcf
            --chr CM057556.1
            --out sharks/GCA_030028105.1/temp/20250113_CM057556.1
            --recode

    After filtering, kept 1 out of 1 Individuals
    Outputting VCF file...
    After filtering, kept 243493 out of a possible 3714560 Sites
    Run Time = 3.00 seconds

    Rscript 20241231_ROH_Durbin_Calc_Eqns.R "sharks/GCA_030028105.1/temp/20250113_CM057556.1.recode.vcf" > sharks/GCA_030028105.1/TEST_CM057556.1_100kb_Results_ROH_Durbin_Calc.txt
I compared sharks/GCA_030028105.1/TEST_CM057556.1_100kb_Results_ROH_Durbin_Calc.txt with sharks/GCA_030028105.1/20250113_CM057556.1_100kb_Results_ROH_Durbin_Calc.txt and they were identical
Given this, despite the errors, I will take the ROH files generated from this and move them forward in the pipeline

Prepped 20250106_Plot_ROH_FROH.sh for GCA_030028105.1 Mobula birostris (Giant manta)
Submitted batch job 3469514

Chrom_list prep for amphibians from 241117.UCSC-hubs-VGP-alignment:
Prep for GCA_027789765.1.fa.gz Dendropsophus ebraccatus (hourglass treefrog) (alternate = GCA_027789725.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/amphibians/GCA_027789765.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/amphibians/GCA_027789765.1.fa.gz | grep '>CM' | wc ## --> 15 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/amphibians/GCA_027789765.1.fa.gz | grep '>CM' | tr -d '>' > amphibians/chrom_lists/GCA_027789765.1_chroms.txt

Prep for GCA_027917425.1.fa.gz Gastrophryne carolinensis (eastern narrow-mouthed toad) (alternate = GCA_027917415.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/amphibians/GCA_027917425.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/amphibians/GCA_027917425.1.fa.gz | grep '>CM' | wc ## --> 11 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/amphibians/GCA_027917425.1.fa.gz | grep '>CM' | tr -d '>' > amphibians/chrom_lists/GCA_027917425.1_chroms.txt

Prep for GCF_028390025.1.fa.gz Pseudophryne corroboree (corroboree frog) (alternate = GCA_028390055.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/amphibians/GCF_028390025.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/amphibians/GCF_028390025.1.fa.gz | grep '>NC' | wc ## --> 12 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/amphibians/GCF_028390025.1.fa.gz | grep '>NC' | tr -d '>' > amphibians/chrom_lists/GCF_028390025.1_chroms.txt

Prep for GCF_029499605.1.fa.gz Hyla sarda (Sardinian treefrog) (alternate = GCA_029493135.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/amphibians/GCF_029499605.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/amphibians/GCF_029499605.1.fa.gz | grep '>NC' | wc --> 13 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/amphibians/GCF_029499605.1.fa.gz | grep '>NC' | tr -d '>' > amphibians/chrom_lists/GCF_029499605.1_chroms.txt

(alternate = GCA_002915635.3.fa.gz) --> NOTE: COULD NOT FIND REFERENCE GENOME FOR THIS ALTERNATE IN MATCHING VGP DIRECTORY --> Ambystoma mexicanum (axolotl)

Prep for GCA_031893055.1.fa.gz Leptodactylus fuscus (rufous frog) (alternate = GCA_031893025.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/amphibians/GCA_031893055.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/amphibians/GCA_031893055.1.fa.gz | grep '>CM' | wc ## --> 12 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/amphibians/GCA_031893055.1.fa.gz | grep '>CM' | tr -d '>' > amphibians/chrom_lists/GCA_031893055.1_chroms.txt

Prep for GCA_035609145.1.fa.gz Eleutherodactylus coqui (Puerto Rican coqui) (alternate = GCA_035609135.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/amphibians/GCA_035609145.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/amphibians/GCA_035609145.1.fa.gz | grep '>CM' | wc ## --> 14 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/amphibians/GCA_035609145.1.fa.gz | grep '>CM' | tr -d '>' > amphibians/chrom_lists/GCA_035609145.1_chroms.txt

Prep for GCA_038501925.1.fa.gz Xenopus petersii (Peter's clawed frog) (alternate = GCA_038501915.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/amphibians/GCA_038501925.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/amphibians/GCA_038501925.1.fa.gz | grep '>CM' | wc ## --> 19 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/amphibians/GCA_038501925.1.fa.gz | grep '>CM' | tr -d '>' > amphibians/chrom_lists/GCA_038501925.1_chroms.txt

Prep for GCA_038048845.1.fa.gz Mixophyes fleayi (Fleay's barred frog) (alternate = GCA_038048865.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/amphibians/GCA_038048845.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/amphibians/GCA_038048845.1.fa.gz | grep '>CM' | wc ## --> 12 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/amphibians/GCA_038048845.1.fa.gz | grep '>CM' | tr -d '>' > amphibians/chrom_lists/GCA_038048845.1_chroms.txt

Prep for GCA_037306005.1.fa.gz Rhinophrynus dorsalis (Mexican burrowing toad) (alternate = GCA_037306015.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/amphibians/GCA_037306005.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/amphibians/GCA_037306005.1.fa.gz | grep '>CM' | wc ## --> 12 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/amphibians/GCA_037306005.1.fa.gz | grep '>CM' | tr -d '>' > amphibians/chrom_lists/GCA_037306005.1_chroms.txt

Finished amphibians chrom lists


#### UPDATE ####
20250114 (January 14th, 2024)

Job 3459408 for mm2 alignment of GCF_030144855.1 Hypanus sabinus (Atlantic stingray) completed
It appears to have completed successfully!

Prep chrom_list for GCF_030144855.1.fa.gz Hypanus sabinus (Atlantic stingray) (alternate = GCA_030144785.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/GCF_030144855.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/GCF_030144855.1.fa.gz | grep '>NC' | wc ## --> 35 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/GCF_030144855.1.fa.gz | grep '>NC' | tr -d '>' > sharks/chrom_lists/GCF_030144855.1_chroms.txt

Prepared 20250101_mm2_alignment.sh script for GCA_035084275.1.fa.gz Hydrolagus colliei (spotted ratfish)
    REF_NAME=GCA_035084275.1
    REFERENCE=../241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/GCA_035084275.1.fa.gz ## Finish filling in depending on genome and clade
    ALTERNATE=../241117.UCSC-hubs-VGP-alignment/alignment/alternate/sharks/GCA_035084065.1.fa.gz ## Finish filling in depending on genome and clade
Updated date for .out/.err files to 20250114
Submitted batch job 3566567

Prepared 20241231_ROH_Calc.sh script for GCF_030144855.1 Hypanus sabinus (Atlantic stingray)
    less sharks/chrom_lists/GCF_030144855.1_chroms.txt | wc ## --> 35 chromosomes
Updated TODAY_DATE and date for .out/.err files
Updated array size to 35 for number of chromosomes
Changed REF_NAME
    REF_NAME=GCF_030144855.1
Submitted batch job 3567616


Job 3469514 for 20250106_Plot_ROH_FROH.sh for GCA_030028105.1 Mobula birostris (Giant manta) completed
It appears to have completed successfully with the exception of the pdf file again being generated with a size of 0 bytes
Ran Rscript in command line to see if result would change
    Rscript 20250106_Plot_ROH.R "sharks/GCA_030028105.1/20250113_GCA_030028105.1_ROH_Results.tsv" "sharks/GCA_030028105.1/GCA_030028105.1_Chroms_Lengths.txt" "20250113" "GCA_030028105.1" > sharks/GCA_030028105.1/20250113_GCA_030028105.1_ROH_Map.pdf
Result came out the same --> output file generated but with file size of 0 bytes
Commented out def.off() commands in 20250106_Plot_ROH.R script and reran
    Rscript 20250106_Plot_ROH.R "sharks/GCA_030028105.1/20250113_GCA_030028105.1_ROH_Results.tsv" "sharks/GCA_030028105.1/GCA_030028105.1_Chroms_Lengths.txt" "20250113" "GCA_030028105.1" > sharks/GCA_030028105.1/20250113_GCA_030028105.1_ROH_Map.pdf
Result came out the same --> output file generated but with file size of 0 bytes
Re-run:
  xvfb-run -a Rscript 20250106_Plot_ROH.R "sharks/GCA_030028105.1/20250113_GCA_030028105.1_ROH_Results.tsv" "sharks/GCA_030028105.1/GCA_030028105.1_Chroms_Lengths.txt" "20250113" "GCA_030028105.1" > sharks/GCA_030028105.1/20250113_GCA_030028105.1_ROH_Map.pdf
This time file came out with a size of 210bytes but could not be opened 
Removed file
    rm -f -r 20250113_GCA_030028105.1_ROH_Map.pdf
Will enter R from command line and run code from 20250106_Plot_ROH.R script 
    dat <- read.table("sharks/GCA_030028105.1/20250113_GCA_030028105.1_ROH_Results.tsv", header=FALSE)
    Chrom <- read.table("sharks/GCA_030028105.1/GCA_030028105.1_Chroms_Lengths.txt", header=FALSE)
    pdf(file = paste("sharks/GCA_030028105.1/", date, "_", ref_name, "_ROH_Map.pdf", sep = "")) 
Successfully generated pdf file 20250113_GCA_030028105.1_ROH_Map.pdf -- I think problem was due to not having dev.off() at end of script to close the pdf() device
Removed earlier dev.off command before pdf -- possibly was prematurely stopping R code there

Job 3566567 of 20250101_mm2_alignment.sh script for GCA_035084275.1.fa.gz Hydrolagus colliei (spotted ratfish) completed
It appears the job completed successfully!

Prepared 20250101_mm2_alignment.sh script for GCA_035084215.1 Heptranchias perlo (sharpnose sevengill shark)
    REF_NAME=GCA_035084215.1 
    REFERENCE=../241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/GCA_035084215.1.fa.gz 
    ALTERNATE=../241117.UCSC-hubs-VGP-alignment/alignment/alternate/sharks/GCA_035084135.1.fa.gz
Submitted batch job 3569083

Prepping chrom_lists for fishes -- finding matches for alternate genomes in references given that there are fewer alternates
Prep for  GCA_944039275.1.fa.gz Danio rerio (zebrafish) (alternate = GCA_903684865.1.fa.gz, GCA_903684855.2.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/fish/GCA_944039275.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/fish/GCA_944039275.1.fa.gz | grep '>OX' | wc ## --> 26 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/fish/GCA_944039275.1.fa.gz | grep '>OX' | tr -d '>' > fishes/chrom_lists/GCA_944039275.1_chroms.txt

Prep for GCA_028022725.1.fa.gz Lycodopsis pacificus (blackbelly eelpout) (alternate = GCA_028021495.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/fish/GCA_028022725.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/fish/GCA_028022725.1.fa.gz | grep '>CM' | wc ## --> 24 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/fish/GCA_028022725.1.fa.gz | grep '>CM' | tr -d '>' > fishes/chrom_lists/GCA_028022725.1_chroms.txt

Prep for GCF_902713425.1.fa.gz Acipenser ruthenus (sterlet) (alternate = GCA_902713435.2.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/fish/GCF_902713425.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/fish/GCF_902713425.1.fa.gz | grep '>NC' | wc ## --> 61 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/fish/GCF_902713425.1.fa.gz | grep '>NC' | tr -d '>' > fishes/chrom_lists/GCF_902713425.1_chroms.txt

Prep for GCF_029633865.1.fa.gz Lampris incognitus (smalleye Pacific opah) (alternate = GCA_029633845.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/fish/GCF_029633865.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/fish/GCF_029633865.1.fa.gz | grep '>NC' | wc ## --> 21 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/fish/GCF_029633865.1.fa.gz | grep '>NC' | tr -d '>' > fishes/chrom_lists/GCF_029633865.1_chroms.txt

Prep for GCF_009769545.1.fa.gz Cyclopterus lumpus (lumpfish) (alternate = GCA_963457625.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/fish/GCF_009769545.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/fish/GCF_009769545.1.fa.gz | grep '>NC' | wc ## --> 25 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/fish/GCF_009769545.1.fa.gz | grep '>NC' | tr -d '>' > fishes/chrom_lists/GCF_009769545.1_chroms.txt

Prep for GCA_029633875.1.fa.gz Hoplias malabaricus (trahira) (alternate = GCA_029633855.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/fish/GCA_029633875.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/fish/GCA_029633875.1.fa.gz | grep '>CM' | wc ## --> 19 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/fish/GCA_029633875.1.fa.gz | grep '>CM' | tr -d '>' > fishes/chrom_lists/GCA_029633875.1_chroms.txt

Prep for GCF_030014385.1.fa.gz Trichomycterus rosablanca (alternate = GCA_030015355.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/fish/GCF_030014385.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/fish/GCF_030014385.1.fa.gz | grep '>NC' | wc ## --> 27 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/fish/GCF_030014385.1.fa.gz | grep '>NC' | tr -d '>' > fishes/chrom_lists/GCF_030014385.1_chroms.txt

Prep for GCA_030463535.1.fa.gz Salminus brasiliensis (dorado) (alternate = GCA_030448965.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/fish/GCA_030463535.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/fish/GCA_030463535.1.fa.gz | grep '>CM' | wc ## --> 26 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/fish/GCA_030463535.1.fa.gz | grep '>CM' | tr -d '>' > fishes/chrom_lists/GCA_030463535.1_chroms.txt

Prep for GCF_901709675.1.fa.gz Syngnathus acus (greater pipefish) (alternate = GCA_948146105.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/fish/GCF_901709675.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/fish/GCF_901709675.1.fa.gz | grep '>NC' | wc ## --> 22 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/fish/GCF_901709675.1.fa.gz | grep '>NC' | tr -d '>' > fishes/chrom_lists/GCF_901709675.1_chroms.txt

Prep for GCA_036373705.1.fa.gz Amia calva (bowfin) (alternate = GCA_036365475.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/fish/GCA_036373705.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/fish/GCA_036373705.1.fa.gz | grep '>CM' | wc ## --> 24 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/fish/GCA_036373705.1.fa.gz | grep '>CM' | tr -d '>' > fishes/chrom_lists/GCA_036373705.1_chroms.txt

Prep for GCA_037039145.1.fa.gz Fundulus diaphanus (banded killifish) (alternate = GCA_037038625.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/fish/GCA_037039145.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/fish/GCA_037039145.1.fa.gz | grep '>CM' | wc ## --> 25 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/fish/GCA_037039145.1.fa.gz | grep '>CM' | tr -d '>' > fishes/chrom_lists/GCA_037039145.1_chroms.txt

Prep for GCA_038024135.1.fa.gz Cyprinella venusta (blacktail shiner) (alternate = GCA_038021265.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/fish/GCA_038024135.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/fish/GCA_038024135.1.fa.gz | grep '>CM' | wc
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/fish/GCA_038024135.1.fa.gz | grep '>CM' | tr -d '>' > fishes/chrom_lists/GCA_038024135.1_chroms.txt

Job 3569083 for 20250101_mm2_alignment.sh script for GCA_035084215.1 Heptranchias perlo (sharpnose sevengill shark) completed
It worked!
Prepared 20250101_mm2_alignment.sh script for GCA_036365495 Heterodontus francisci (horn shark)
    REF_NAME=GCA_036365495.1
    REFERENCE=../241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/GCA_036365495.1.fa.gz ## Finish filling in depending on genome and clade
    ALTERNATE=../241117.UCSC-hubs-VGP-alignment/alignment/alternate/sharks/GCA_036365525.1.fa.gz ## Finish filling in depending on genome and clade
Submitted batch job 3580721


#### UPDATE ####
20250115 (January 15th, 2025)

Job 3580721 for 20250101_mm2_alignment.sh script for GCA_036365495 Heterodontus francisci (horn shark) completed
It worked!
Prepared 20250101_mm2_alignment.sh script for GCA_036971175.1 Narcine bancroftii (Spotted torpedo ray/Caribbean electric ray)
    REF_NAME=GCA_036971175.1
    REFERENCE=../241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/GCA_036971175.1.fa.gz ## Finish filling in depending on genome and clade
    ALTERNATE=../241117.UCSC-hubs-VGP-alignment/alignment/alternate/sharks/GCA_036971175.1.fa.gz ## Finish filling in depending on genome and clade
Submitted batch job 3631327


Job 3631327 for 20250101_mm2_alignment.sh script for GCA_036971175.1 Narcine bancroftii (Spotted torpedo ray/Caribbean electric ray) completed
It worked!
I have now completed mm2 alignments for other/shark diploid genomes for VGP -- I will now move onto reptiles
Prepared 20250101_mm2_alignment.sh script for GCF_025201925.1 Gopherus flavomarginatus (Mexican gopher tortoise) (alternate = GCA_025201965.1.fa.gz)
    REF_NAME=GCF_025201925.1
    REFERENCE=../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCF_025201925.1.fa.gz ## Finish filling in depending on genome and clade
    ALTERNATE=../241117.UCSC-hubs-VGP-alignment/alignment/alternate/reptiles/GCA_025201965.1.fa.gz ## Finish filling in depending on genome and clade
Submitted batch job 3646868

Prepared 20250106_Plot_ROH_FROH.sh script for GCF_020745735.1 Hemiscyllium ocellatum (Epaulette shark)
    REF_NAME=GCF_020745735.1 
    TODAY_DATE=20250115
    ROH_CALC_DATE=20250110
Submitted batch job 3647558

Job 3646868 completed for GCF_025201925.1 Gopherus flavomarginatus (Mexican gopher tortoise) (alternate = GCA_025201965.1.fa.gz)
It worked!
Prepared 20250101_mm2_alignment.sh script for GCF_030867095.1 Alligator mississippiensis (American alligator) (alternate = GCA_030867065.1.fa.gz)
    REF_NAME=GCF_030867095.1 ## Change for reference genome of each species
    REFERENCE=../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCF_030867095.1.fa.gz ## Finish filling in depending on genome and clade
    ALTERNATE=../241117.UCSC-hubs-VGP-alignment/alignment/alternate/reptiles/GCA_030867065.1.fa.gz ## Finish filling in depending on genome and clade
Submitted batch job 3648801

Created 20250115_find_het_calc.R script 
Will use this script to calculate heterozygosity for genomes
Starting with simulated chromosome with controlled heterozygosity to make sure calculations are correct


#### UPDATE ####
20250116 (January 16th, 2024)

Job 3648801 completed to do mm2 alingment for GCF_030867095.1 Alligator mississippiensis (American alligator) (alternate = GCA_030867065.1.fa.gz)
It worked!
Prepared 20250101_mm2_alignment.sh script for GCF_027887155.1 Malaclemys terrapin pileata (Diamondback terrapin) (alternate = GCA_027887205.1.fa.gz)
    REF_NAME=GCF_027887155.1
    REFERENCE=../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCF_027887155.1.fa.gz ## Finish filling in depending on genome and clade
    ALTERNATE=../241117.UCSC-hubs-VGP-alignment/alignment/alternate/reptiles/GCA_027887205.1.fa.gz ## Finish filling in depending on genome and clade
Submitted batch job 3685140

Job 3685140 completed of mm2 alignment for GCF_027887155.1 Malaclemys terrapin pileata (Diamondback terrapin) (alternate = GCA_027887205.1.fa.gz)
It worked!
Prepared 20250101_mm2_alignment.sh script for GCA_028017835.1 Emys orbicularis (European pond turtle) (alternate = GCA_028017845.1.fa.gz)
    REF_NAME=GCA_028017835.1
    REFERENCE=../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCA_028017835.1.fa.gz ## Finish filling in depending on genome and clade
    ALTERNATE=../241117.UCSC-hubs-VGP-alignment/alignment/alternate/reptiles/GCA_028017845.1.fa.gz ## Finish filling in depending on genome and clade
Submitted batch job 3685529

Job 3647558 completed of 20250106_Plot_ROH_FROH.sh script for GCF_020745735.1 Hemiscyllium ocellatum (Epaulette shark)
Technically completed but had some errors in running:
    awk: fatal: cannot open file `sharks/GCF_020745735.1/20250110_NC_083454.1_100kb_Results_ROH_Durbin_Calc.txt' for reading (No such file or directory)

    /var/spool/slurm/slurmd/job3647558/slurm_script: line 63: -: syntax error: operand expected (error token is "-")
    Error in names(x) <- value : 
    'names' attribute [2] must be the same length as the vector [1]
    Calls: colnames<-
Removed generated files from this run to avoid errors
    rm -f -r GCF_020745735.1_Chroms_Lengths.txt
    rm -f -r 20250115_GCF_020745735.1_ROH_Results.tsv
    rm -f -r 20250115_GCF_020745735.1_FROH.txt
    rm -f -r 20250115_GCF_020745735.1_ROH_Map.pdf
When checking errors -- for first error not finding a txt file for chromosome NC_083454.1 -- there is not text file generated for this chromosome
There is one empty txt file w/o a chromosome associated with it --> 20250110__100kb_Results_ROH_Durbin_Calc.txt
Manually re-running for missing chromosome:
    vcftools --vcf sharks/GCF_020745735.1/GCF_020745735.1_aligned.mm2.vcf --chr NC_083454.1 --recode --out sharks/GCF_020745735.1/temp/20250110_NC_083454.1

    VCFtools - 0.1.17
    (C) Adam Auton and Anthony Marcketta 2009

    Parameters as interpreted:
            --vcf sharks/GCF_020745735.1/GCF_020745735.1_aligned.mm2.vcf
            --chr NC_083454.1
            --out sharks/GCF_020745735.1/temp/20250110_NC_083454.1
            --recode

    After filtering, kept 1 out of 1 Individuals
    Outputting VCF file...
    After filtering, kept 3558 out of a possible 5127864 Sites
    Run Time = 3.00 seconds

    Rscript 20241231_ROH_Durbin_Calc_Eqns.R "sharks/GCF_020745735.1/temp/20250110_NC_083454.1.recode.vcf" > sharks/GCF_020745735.1/20250110_NC_083454.1_100kb_Results_ROH_Durbin_Calc.txt

    ls -ltr | grep .1_100kb_Results_ROH_Durbin_Calc.txt | wc --> confirmed 54 files are present, one for each chromosome
Removed the empty chromosome ROH file:
    rm -f -r 20250110__100kb_Results_ROH_Durbin_Calc.txt
Addressed second error by changing syntax in line 63 to:
    chrom_length=$((end - start))
Resubmitted 20250106_Plot_ROH_FROH.sh for GCF_020745735.1 Hemiscyllium ocellatum (Epaulette shark) 
Submitted batch job 3685749

Prepared 20241231_ROH_Calc.sh script for GCA_035084275.1 Hydrolagus colliei (Pacific ratfish)
    TODAY_DATE=20250116
    REF_NAME=GCA_035084275.1
Submitted batch job 3685800

Job 3685529 completed for mm2 alignment of GCA_028017835.1 Emys orbicularis (European pond turtle) (alternate = GCA_028017845.1.fa.gz)
It worked!
Prepared 20250101_mm2_alignment.sh script for GCF_029931775.1.fa.gz Euleptes europaea (tarantolino) (alternate = GCA_029931755.1.fa.gz)
    REF_NAME=GCF_029931775.1
    REFERENCE=../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCF_029931775.1.fa.gz ## Finish filling in depending on genome and clade
    ALTERNATE=../241117.UCSC-hubs-VGP-alignment/alignment/alternate/reptiles/GCA_029931755.1.fa.gz ## Finish filling in depending on genome and clade
Submitted batch job 3686881

Job 3686881 completed for mm2 alignment of GCF_029931775.1.fa.gz Euleptes europaea (tarantolino) (alternate = GCA_029931755.1.fa.gz)
It worked!
Prepared 20250101_mm2_alignment.sh for GCA_030020295.1.fa.gz Gavialis gangeticus (Gharial) (altnerate = GCA_030020385.1.fa.gz)
    REF_NAME=GCA_030020295.1 ## Change for reference genome of each species
    REFERENCE=../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCA_030020295.1.fa.gz ## Finish filling in depending on genome and clade
    ALTERNATE=../241117.UCSC-hubs-VGP-alignment/alignment/alternate/reptiles/GCA_030020385.1.fa.gz ## Finish filling in depending on genome and clade
Submitted batch job 3687513

Job 3685800 for ROH calculation of GCA_035084275.1 Hydrolagus colliei (Pacific ratfish) completed
Realized that I didn't change the array number to match the number of chromosomes this species hasn
Removed .err/.out files in Out_Err_Files directory to save space
    rm -f -r 20250116_ROH_Calc_GCA_035084275.1*
Removed all output ROH txt files in sharks/GCA_035084275.1/ directory
    rm -f -r *_100kb_Results_ROH_Durbin_Calc.txt
Changed sbatch array to 1-40 in 20241231_ROH_Calc.sh
Resubmitted
Submitted batch job 3687765

Job 3685749 completed for plotting ROH/calculating FROH for GCF_020745735.1 Hemiscyllium ocellatum (Epaulette shark) 
Job says it completed successfully but the .err files contains these errors:
    Error in names(x) <- value : 
    'names' attribute [2] must be the same length as the vector [1]
    Calls: colnames<-
    Execution halted
Checked sharks/GCF_020745735.1 directory and files did not generate Successfully
GCF_020745735.1_Chroms_Lengths.txt has size of 11bytes and only contains length of 1 chromosome
FROH .txt file and ROH_Map.pdf both have file size of 0
Removing generated files so that they do not mess up future runs of code
    rm -f -r GCF_020745735.1_Chroms_Lengths.txt
    rm -f -r 20250116_GCF_020745735.1_ROH_Results.tsv
    rm -f -r 20250116_GCF_020745735.1_FROH.txt
    rm -f -r 20250116_GCF_020745735.1_ROH_Map.pdf
Looking at zcat command on line 37 to generate the Chroms_Lengths .txt file, I think it was due to the fact that it was looking for chromosome starting with 'CM' instead of 'NC'
Manually running command in command line to test theory
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/GCF_020745735.1.fa.gz | awk '$0 ~ ">NC" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' >sharks/GCF_020745735.1/GCF_020745735.1_Chroms_Lengths.txt
Command ran successfully with file properly generated
Resubmitted 20250106_Plot_ROH_FROH.sh
Submitted batch job 3688103

Job 3687513 completed for mm2 alignment of GCA_030020295.1.fa.gz Gavialis gangeticus
It worked!
Prepared 20250101_mm2_alignment.sh for GCF_030035675.1.fa.gz Rhineura floridana (Florida worm lizard) (alternate = GCA_030035715.1.fa.gz)
    REF_NAME=GCF_030035675.1
    REFERENCE=../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCF_030035675.1.fa.gz ## Finish filling in depending on genome and clade
    ALTERNATE=../241117.UCSC-hubs-VGP-alignment/alignment/alternate/reptiles/GCA_030035715.1.fa.gz ## Finish filling in depending on genome and clade
Submitted batch job 3688353

Job 3688103 completed for plotting ROH/calculating FROH for GCF_020745735.1 Hemiscyllium ocellatum (Epaulette shark) 
Successfully generated Chroms_Lengths.txt file, .tsv file, and FROH.txt file
Once again the ROH_Map.pdf file was not successfully generated, it was made but only has a size of 26 bytes and does not open (says invalid/corrupted)
Looking at 20250106_Plot_ROH.R, it's possible it's due to pdf() name in R script not matching output in shell script
Fixed this
Will rerun Rscript separately in command line
    Rscript 20250106_Plot_ROH.R "sharks/GCF_020745735.1/20250116_GCF_020745735.1_ROH_Results.tsv" "sharks/GCF_020745735.1/GCF_020745735.1_Chroms_Lengths.txt" 20250116 "GCF_020745735.1" "sharks" > sharks/GCF_020745735.1/20250116_GCF_020745735.1_ROH_Map.pdf
ROH_Map.pdf file was successfully generated!

Job 3687765 completed for creating ROH .txt files for each chromosome of GCA_035084275.1 Hydrolagus colliei (Pacific ratfish)
It was successful with the exception of not creating a txt file for the last chromosome CM068781.1
Again, also a txt file without a chromosome in the name and a size of 0bytes: 20250116__100kb_Results_ROH_Durbin_Calc.txt
Attmped to fix this by editing line 39 in 20241231_ROH_Calc.sh, where chromosome looked at in each array is specified
    CHROM=$(awk "NR==$SLURM_ARRAY_TASK_ID" ${ALL_CHROMS}) ## removed the -1 after $SLURM_ARRAY_TASK_ID
Removed ROH .txt files and resubmitted script to see if this workspace
Submitted batch job 3692500


#### UPDATE ####
20250117 (January 17th, 2025)

Job 3688353 completed for mm2 alignment of GCF_030035675.1.fa.gz Rhineura floridana (Florida worm lizard) (alternate = GCA_030035715.1.fa.gz)
It worked!
NOTE: While going through the sharks results to find the next sp. to run through 20250106_Plot_ROH_FROH.sh script, 
      I found out that the mm2 files for GCA_036971175.1 Narcine bancroftii (Spotted torpedo ray/Caribbean electric ray) are incorrect
      Will have to re-run mm2 alignment
Prepped 20250101_mm2_alignment.sh for GCA_036971175.1 Narcine bancroftii (Spotted torpedo ray/Caribbean electric ray)
    REF_NAME=GCA_036971175.1 
    REFERENCE=../241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/GCA_036971175.1.fa.gz ## Finish filling in depending on genome and clade
    ALTERNATE=../241117.UCSC-hubs-VGP-alignment/alignment/alternate/sharks/GCA_036971445.1.fa.gz ## Finish filling in depending on genome and clade
Realized problem from before was that I had a syntax error in the alternate genome file name
Removed old mm2 files to avoid any issues with new file generation
    rm -f -r GCA_036971175.1.mm2.paf
    rm -f -r GCA_036971175.1.mm2.srt.paf
    rm -f -r GCA_036971175.1_aligned.mm2.vcf
Submitted batch job 3725725

Job 3692500 completed for ROH calculation for GCA_035084275.1 Hydrolagus colliei (Pacific ratfish)
It worked! 
ROH .txt files were properly generated
NOTE: Will want to rerun 20241231_ROH_Calc.sh on GCF_030144855.1 Hypanus sabinus (Atlantic stingray)
      Two of the ROH .txt files are empty and another doesn't have a chromosome associated with it 
Prepared 20241231_ROH_Calc.sh for ROH calculation for GCF_030144855.1 Hypanus sabinus (Atlantic stingray)
    #SBATCH --array=1-35 ##Number of chromosomes your species have --> changed for 35 chromosomes for this sp.
    TODAY_DATE=20250117
    REF_NAME=GCF_030144855.1 ## Change for reference genome of each species
Removed old ROH .txt files to avoid them corrupting re-running the ROH calculation
    rm -f -r *.1_100kb_Results_ROH_Durbin_Calc.txt
    rm -f -r *_100kb_Results_ROH_Durbin_Calc.txt
Submitted batch job 3726045

Job 3688103 completed for plotting ROH/calculating FROH for GCF_020745735.1 Hemiscyllium ocellatum (Epaulette shark) 
Job was completed successfully!
ROH_Map.pdf, _ROH_Results.tsv, and FROH.txt file were all generated
Prepared 20250106_Plot_ROH_FROH.sh for GCA_035084275.1 Hydrolagus colliei (Pacific ratfish)
    TODAY_DATE=20250117
    ROH_CALC_DATE=20250116
    REF_NAME=GCA_035084275.1 
Changed awk in zcat command on line 37 for chromosomes starting with 'CM'
Submitted batch job 3725926

Job 3725926 completed for plotting ROH/calculating FROH for GCA_035084275.1 Hydrolagus colliei (Pacific ratfish)
It worked!
ROH_Map.pdf, _ROH_Results.tsv, and FROH.txt file were all generated

Job 3726045 completed for plotting ROH/calculating FROH for GCF_030144855.1 Hypanus sabinus (Atlantic stingray)
On some, but not all, of the .err files these errors appeared:
    No data left for analysis!
    Error in read.table(args[6]) : no lines available in input
    Execution halted
Looking at the ROH .txt files -- it appears that three of them came out empty:
    20250117_NC_082739.1_100kb_Results_ROH_Durbin_Calc.txt
    20250117_NC_082740.1_100kb_Results_ROH_Durbin_Calc.txt
    20250117_NC_082738.1_100kb_Results_ROH_Durbin_Calc.txt
When cross referencing these with the .err files which show the errors
    20250117_ROH_Calc_GCF_030144855.1_3726814.err --> NC_082739.1
    20250117_ROH_Calc_GCF_030144855.1_3726045.err --> NC_082740.1
    20250117_ROH_Calc_GCF_030144855.1_3726729.err --> NC_082738.1
All three are the ones with the empty txt file results
Will attempt to manually re-run the script for these three chromosomes to see if the errors still arise
    vcftools --vcf sharks/GCF_030144855.1/GCF_030144855.1_aligned.mm2.vcf --chr NC_082739.1 --recode --out sharks/GCF_030144855.1/temp/20250117_NC_082739.1

    After filtering, kept 0 out of a possible 5243846 Sites
    No data left for analysis!
After getting this result, I looked through the entire vcf file for something from this chromosome:
    less sharks/GCF_030144855.1/GCF_030144855.1_aligned.mm2.vcf | grep NC_082739.1

    ##contig=<ID=NC_082739.1,length=45031449>
It appears that there are no variants present on this chromosome
    less sharks/GCF_030144855.1/GCF_030144855.1_aligned.mm2.vcf | grep NC_082740.1
    ##contig=<ID=NC_082740.1,length=1576620>

    less sharks/GCF_030144855.1/GCF_030144855.1_aligned.mm2.vcf | grep NC_082738.1
    ##contig=<ID=NC_082738.1,length=69719156>
The same with the other two chromosomes --> there are no variants present.
Given this, I think this species is ready to move onto to plotting ROH/calculating FROH

Prepared 20250106_Plot_ROH_FROH.sh for GCF_030144855.1 Hypanus sabinus (Atlantic stingray)
    TODAY_DATE=20250117
    ROH_CALC_DATE=20250117
    REF_NAME=GCF_030144855.1 
Changed awk in zcat command on line 37 for chromosomes starting with 'NC'
Submitted batch job 3731568

Preparing chrom_lists for the primates in the VGP/20241117.UCSC-hubs-BGP-alignment/reference and alternate directories

Job 3725725 completed for mm2 alignment of GCA_036971175.1 Narcine bancroftii (Spotted torpedo ray/Caribbean electric ray)
It worked!
Prepared 20250101_mm2_alignment.sh for GCA_030412105.1.fa.gz Cyclura pinguis (alternate = GCA_030412085.1.fa.gz)
    REF_NAME=GCA_030412105.1 
    CLADE=reptiles
    REFERENCE=../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCA_030412105.1.fa.gz ## Finish filling in depending on genome and clade
    ALTERNATE=../241117.UCSC-hubs-VGP-alignment/alignment/alternate/reptiles/GCA_030412085.1.fa.gz ## Finish filling in depending on genome and clade
Submitted batch job 3731771

Prepared 20241231_ROH_Calc.sh for GCA_035084215.1 Heptranchias perlo (sharpnose sevengill shark)
    #SBATCH --array=1-46 ##Number of chromosomes your species have
    TODAY_DATE=20250117
    CLADE=sharks
    REF_NAME=GCA_035084215.1 
Submitted batch job 3731971

Job 3731568 completed for plotting ROH/calculating FROH for GCF_030144855.1 Hypanus sabinus (Atlantic stingray)
It worked!

Job 3731971 completed for calculating the number of ROH for GCA_035084215.1 Heptranchias perlo (sharpnose sevengill shark)
It technically worked, but 15 chromosomes have empty ROH .txt files
Prepared 20250106_Plot_ROH_FROH.sh script for GCA_035084215.1 Heptranchias perlo (sharpnose sevengill shark)
    CLADE=sharks
    REF_NAME=GCA_035084215.1
    REF_NAME=GCA_035084215.1
Changed awk in zcat command on line 37 for chromosomes starting with 'CM'
Submitted batch job 3733419

Job 3731771 completed for mm2 alignment of GCA_030412105.1.fa.gz Cyclura pinguis (alternate = GCA_030412085.1.fa.gz)
It worked!
Preparing 20250101_mm2_alignment.sh for GCA_030867105.1.fa.gz Liasis olivaceus (alternate = GCA_030867085.1.fa.gz)
    REF_NAME=GCA_030867105.1 
    CLADE=reptiles
    REFERENCE=../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCA_030867105.1.fa.gz ## Finish filling in depending on genome and clade
    ALTERNATE=../241117.UCSC-hubs-VGP-alignment/alignment/alternate/reptiles/GCA_030867085.1.fa.gz ## Finish filling in depending on genome and clade
Submitted batch job 3733126

Prepared 20241231_ROH_Calc.sh for GCA_036365495.1 Heterodontus francisci (horn shark)
    #SBATCH --array=1-51 ##Number of chromosomes your species have
    TODAY_DATE=20250117 
    CLADE=sharks
    REF_NAME=GCA_036365495.1
Submitted batch job 3733394

Job 3733126 completed for mm2 alignment of GCA_030867105.1.fa.gz Liasis olivaceus (alternate = GCA_030867085.1.fa.gz)
It worked!
Prepared 20250101_mm2_alignment.sh for GCA_033349115.1.fa.gz Macrochelys suwanniensis (Suwannee snapping turtle) (alternate = GCA_033296515.1.fa.gz)
    REF_NAME=GCA_033349115.1
    CLADE=reptiles
    REFERENCE=../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCA_033349115.1.fa.gz ## Finish filling in depending on genome and clade
    ALTERNATE=../241117.UCSC-hubs-VGP-alignment/alignment/alternate/reptiles/GCA_033296515.1.fa.gz ## Finish filling in depending on genome and clade
Submitted batch job 3733533


#### UPDATE ####
20250118 (January 18th, 2025)

Job 3733533 completed for mm2 alignment of GCA_033349115.1.fa.gz Macrochelys suwanniensis (Suwannee snapping turtle) (alternate = GCA_033296515.1.fa.gz)
It worked!
Prepared 20250101_mm2_alignment.sh for GCA_035046505.1.fa.gz Tiliqua scincoides (alternate = GCA_035046495.1.fa.gz)
    REF_NAME=GCA_035046505.1 
    CLADE=reptiles
    REFERENCE=../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCA_035046505.1.fa.gz ## Finish filling in depending on genome and clade
    ALTERNATE=../241117.UCSC-hubs-VGP-alignment/alignment/alternate/reptiles/GCA_035046495.1.fa.gz ## Finish filling in depending on genome and clade
Submitted batch job 3797899

Job 3733394 completed for ROH calculation of GCA_036365495.1 Heterodontus francisci (horn shark)
It worked!

Prepared 20250106_Plot_ROH_FROH.sh for GCA_036365495.1 Heterodontus francisci (horn shark)
    TODAY_DATE=20250118
    ROH_CALC_DATE=20250117
    CLADE=sharks
    REF_NAME=GCA_036365495.1
awk in zcat command on line 37 already set for chromosomes starting with 'CM'
Submitted batch job 3797955

Prepared 20241231_ROH_Calc.sh for GCA_036971175.1 Narcine bancroftii (Spotted torpedo ray/Caribbean electric ray)
    #SBATCH --array=1-14 
    TODAY_DATE=20250118
    CLADE=sharks
    REF_NAME=GCA_036971175.1
Submitted batch job 3798000

Job 3797955 completed for plotting ROH/calculating FROH for GCA_036365495.1 Heterodontus francisci (horn shark)
It worked!
ROH_Map.pdf, _ROH_Results.tsv, and FROH.txt file were all generated

Job 3797899 completed for mm2 alignment of GCA_035046505.1.fa.gz Tiliqua scincoides (alternate = GCA_035046495.1.fa.gz)
It worked!

Job 3798000 completed for ROH caluclation for GCA_036971175.1 Narcine bancroftii (Spotted torpedo ray/Caribbean electric ray)
It worked!

Prepared 20250106_Plot_ROH_FROH.sh for GCA_036971175.1 Narcine bancroftii (Spotted torpedo ray/Caribbean electric ray)
    TODAY_DATE=20250118
    ROH_CALC_DATE=20250117
    CLADE=sharks
    REF_NAME=GCA_036971175.1
awk in zcat command on line 37 already set for chromosomes starting with 'CM'
Submitted batch job 3798474

Prepared 20250101_mm2_alignment.sh for GCF_035149785.1.fa.gz Candoia aspera (alternate = GCA_035125265.1.fa.gz)
    REF_NAME=GCF_035149785.1
    CLADE=reptiles
    REFERENCE=../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCF_035149785.1.fa.gz ## Finish filling in depending on genome and clade
    ALTERNATE=../241117.UCSC-hubs-VGP-alignment/alignment/alternate/reptiles/GCA_035125265.1.fa.gz ## Finish filling in depending on genome and clade
Submitted batch job 3798519

Prepared 20241231_ROH_Calc.sh for GCF_025201925.1 Gopherus flavomarginatus (Mexican gopher tortoise) (alternate = GCA_025201965.1.fa.gz)
    #SBATCH --array=1-25 
    TODAY_DATE=20250118
    CLADE=reptiles
    REF_NAME=GCF_025201925.1
Submitted batch job 3798607

Job 3798474 completed for plotting ROH/calculating FROH for GCA_036971175.1 Narcine bancroftii (Spotted torpedo ray/Caribbean electric ray)
FROH .txt file, ROH .tsv file, and ROH .pdf map all have a file size of 0 bytes
Realized error was an incorrect ROH_CALC_DATE (20250117 when it should be 20250118)
Fixed error and removed empty files
    rm -f -r 20250118_GCA_036971175.1_ROH_Results.tsv
    rm -f -r 20250118_GCA_036971175.1_FROH.txt
    rm -f -r 20250118_GCA_036971175.1_ROH_Map.pdf
Resubmitted
Submitted batch job 3798755

Job 3798519 completed for mm2 alignment of GCF_035149785.1.fa.gz Candoia aspera (alternate = GCA_035125265.1.fa.gz)
It worked!
Prepared 20250101_mm2_alignment.sh for GCA_037176765.1.fa.gz Anolis sagrei (Brown anole) (alternate = GCA_037176775.1.fa.gz)
    REF_NAME=GCA_037176765.1
    CLADE=reptiles
    REFERENCE=../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCA_037176765.1.fa.gz ## Finish filling in depending on genome and clade
    ALTERNATE=../241117.UCSC-hubs-VGP-alignment/alignment/alternate/reptiles/GCA_037176775.1.fa.gz ## Finish filling in depending on genome and clade
Submitted batch job 3798988

Job 3798755 completed for plotting ROH/calculating FROH for GCA_036971175.1 Narcine bancroftii (Spotted torpedo ray/Caribbean electric ray)
It worked!
With this completed, all sharks will have had ROH mapped and FROH calculated 

Prepared 20250106_Plot_ROH_FROH.sh for GCF_025201925.1 Gopherus flavomarginatus (Mexican gopher tortoise) (alternate = GCA_025201965.1.fa.gz)
    TODAY_DATE=20250118
    ROH_CALC_DATE=20250118
    CLADE=reptiles
    REF_NAME=GCF_025201925.1
Changed awk in zcat command on line 37 to chromosomes starting with 'NC', also changed reference file directory to go to reptiles
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/${REF_NAME}.fa.gz | awk '$0 ~ ">NC" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' >${SPECIES}/${REF_NAME}_Chroms_Lengths.txt
Submitted batch job 3799202

Prepared 20241231_ROH_Calc.sh for GCF_030867095.1 Alligator mississippiensis (American alligator) (alternate = GCA_030867065.1.fa.gz)
    #SBATCH --array=1-17
    TODAY_DATE=20250118
    CLADE=reptiles
    REF_NAME=GCF_025201925.1
Submitted batch job 3799248

Job 3799202 completed for plotting ROH/calculating FROH for  GCF_025201925.1 Gopherus flavomarginatus (Mexican gopher tortoise) (alternate = GCA_025201965.1.fa.gz)
Job appears to be successful! Although some ROH .txt files appear to have been modified after the output files were formed -- odd


#### UPDATE ####
20250119 (January 19th, 2025)

Job 3798988 finished for mm2 alignment of GCA_037176765.1.fa.gz Anolis sagrei (Brown anole) (alternate = GCA_037176775.1.fa.gz)
Did not complete but was canceled due to the time limit
Upped time limit to 24 hours
Updated dates for .out/.err files to 20250119
Removed .paf file created to avoid any errors due to job being interrupted
    rm -f -r GCA_037176765.1.mm2.paf
Resubmitted job
Submitted batch job 3855997

Job 3799248 completed for ROH calculation of GCF_030867095.1 Alligator mississippiensis (American alligator) (alternate = GCA_030867065.1.fa.gz)
Realized I incorrectly ran job with reference genome name from Gopherus flavomarginatus
Remove all generated files from both GCF_030867095.1 Alligator mississippiensis (American alligator) and GCF_025201925.1 Gopherus flavomarginatus (Mexican gopher tortoise) directories
In reptiles/GCF_025201925.1/
    rm -f -r *_100kb_Results_ROH_Durbin_Calc.txt
    rm -f -r 20250118_GCF_025201925.1_ROH_Results.tsv
    rm -f -r 20250118_GCF_025201925.1_FROH.txt
    rm -f -r 20250118_GCF_025201925.1_ROH_Map.pdf
    rm -f -r GCF_025201925.1_Chroms_Lengths.txt
Double checked that names were correct for mm2 alignment of GCF_030867095.1 Alligator mississippiensis 
Prepared 20241231_ROH_Calc.sh for GCF_030867095.1 Alligator mississippiensis (American alligator) (alternate = GCA_030867065.1.fa.gz)
    #SBATCH --array=1-17
    TODAY_DATE=20250119
    CLADE=reptiles
    REF_NAME=GCF_030867095.1 
Resubmitted job
Submitted batch job 3856427

Job 3856427 completed for ROH calculation of GCF_030867095.1 Alligator mississippiensis (American alligator) (alternate = GCA_030867065.1.fa.gz)
It worked!
Prepared 20250106_Plot_ROH_FROH.sh for GCF_030867095.1 Alligator mississippiensis (American alligator) (alternate = GCA_030867065.1.fa.gz)
    TODAY_DATE=20250119
    ROH_CALC_DATE=20250119
    CLADE=reptiles
    REF_NAME=GCF_030867095.1
awk in zcat command on line 37 already set for chromosomes starting with 'NC'
Submitted batch job 3856590

Prepared 20241231_ROH_Calc.sh for GCF_027887155.1 Malaclemys terrapin pileata (Diamondback terrapin) (alternate = GCA_027887205.1.fa.gz)
    #SBATCH --array=1-25 
    TODAY_DATE=20250119
    CLADE=reptiles
    REF_NAME=GCF_027887155.1
Submitted batch job 3856719

Job 3856590 completed for plotting ROH/calculating FROH for GCF_030867095.1 Alligator mississippiensis (American alligator) (alternate = GCA_030867065.1.fa.gz)
It worked!

Job 3856719 completed for calculating ROH for GCF_027887155.1 Malaclemys terrapin pileata (Diamondback terrapin) (alternate = GCA_027887205.1.fa.gz)
It worked!
Prepared 20250106_Plot_ROH_FROH.sh for GCF_027887155.1 Malaclemys terrapin pileata (Diamondback terrapin) (alternate = GCA_027887205.1.fa.gz)
    TODAY_DATE=20250119
    ROH_CALC_DATE=20250119
    CLADE=reptiles
    REF_NAME=GCF_027887155.1
awk in zcat command on line 37 already set for chromosomes starting with 'NC'
Submitted batch job 3857075

Prepared 20241231_ROH_Calc.sh for GCA_028017835.1 Emys orbicularis (European pond turtle) (alternate = GCA_028017845.1.fa.gz)
    #SBATCH --array=1-25 
    TODAY_DATE=20250119
    CLADE=reptiles
    REF_NAME=GCA_028017835.1 
Submitted batch job 3857200

Job 3857075 completed for plotting ROH/calculating FROH for GCF_027887155.1 Malaclemys terrapin pileata (Diamondback terrapin) (alternate = GCA_027887205.1.fa.gz)
It worked!

Job 3857200 completed for calculating ROH for GCA_028017835.1 Emys orbicularis (European pond turtle) (alternate = GCA_028017845.1.fa.gz)
It worked!
Prepared 20250106_Plot_ROH_FROH.sh for GCA_028017835.1 Emys orbicularis (European pond turtle) (alternate = GCA_028017845.1.fa.gz)
    TODAY_DATE=20250119
    ROH_CALC_DATE=20250119
    CLADE=reptiles
    REF_NAME=GCA_028017835.1
Changed awk in zcat command on line 37 for chromosomes starting with 'CM'
Submitted batch job 3857448

Prepared 20241231_ROH_Calc.sh for GCF_029931775.1.fa.gz Euleptes europaea (tarantolino) (alternate = GCA_029931755.1.fa.gz)
    #SBATCH --array=1-21
    TODAY_DATE=20250119
    CLADE=reptiles
    REF_NAME=GCF_029931775.1 
Submitted batch job 3857725

Job 3857448 completed for plotting ROH/calculating FROH for GCA_028017835.1 Emys orbicularis (European pond turtle) (alternate = GCA_028017845.1.fa.gz)
It worked!

Job 3857725 completed for calculating ROH for GCF_029931775.1.fa.gz Euleptes europaea (tarantolino) (alternate = GCA_029931755.1.fa.gz)
It worked!
Prepared 20250106_Plot_ROH_FROH.sh for GCF_029931775.1.fa.gz Euleptes europaea (tarantolino) (alternate = GCA_029931755.1.fa.gz)
    TODAY_DATE=20250119
    ROH_CALC_DATE=20250119
    CLADE=reptiles
    REF_NAME=GCF_029931775.1
Changed awk in zcat command on line 37 for chromosomes starting with 'NC'
Submitted batch job 3867377

Prepared 20241231_ROH_Calc.sh for GCA_030020295.1.fa.gz Gavialis gangeticus (Gharial) (altnerate = GCA_030020385.1.fa.gz)
    #SBATCH --array=1-16
    TODAY_DATE=20250119
    CLADE=reptiles
    REF_NAME=GCA_030020295.1 
Submitted batch job 3858344

Job 3858344 completed for ROH calculation for GCA_030020295.1.fa.gz Gavialis gangeticus (Gharial) (altnerate = GCA_030020385.1.fa.gz)
It worked!

Job 3867377 completed for plotting ROH/calculating FROH for GCF_029931775.1.fa.gz Euleptes europaea (tarantolino) (alternate = GCA_029931755.1.fa.gz)
It worked!

Prepared 20250106_Plot_ROH_FROH.sh for GCA_030020295.1.fa.gz Gavialis gangeticus (Gharial) (altnerate = GCA_030020385.1.fa.gz)
    TODAY_DATE=20250119
    ROH_CALC_DATE=20250119
    CLADE=reptiles
    REF_NAME=GCA_030020295.1
Changed awk in zcat command on line 37 for chromosomes starting with 'CM'
Submitted batch job 3867510

Prepared 20241231_ROH_Calc.sh for GCF_030035675.1.fa.gz Rhineura floridana (Florida worm lizard) (alternate = GCA_030035715.1.fa.gz)
    #SBATCH --array=1-22
    TODAY_DATE=20250119
    CLADE=reptiles
    REF_NAME=GCF_030035675.1 ## Change for reference genome of each species
Submitted batch job 3867559


#### UPDATE ####
20250120 (January 20th, 2024)

Job 3867559 completed for ROH calculation of GCF_030035675.1.fa.gz Rhineura floridana (Florida worm lizard) (alternate = GCA_030035715.1.fa.gz)
It worked!

Job 3867510 completed for plotting ROH/calculating FROH for GCA_030020295.1.fa.gz Gavialis gangeticus (Gharial) (altnerate = GCA_030020385.1.fa.gz)
It worked!

Job 3855997 completed for mm2 alignment of GCA_037176765.1.fa.gz Anolis sagrei (Brown anole) (alternate = GCA_037176775.1.fa.gz)
It worked!

Prepared 20250106_Plot_ROH_FROH.sh for GCF_030035675.1.fa.gz Rhineura floridana (Florida worm lizard) (alternate = GCA_030035715.1.fa.gz)
    TODAY_DATE=20250120
    ROH_CALC_DATE=20250119
    CLADE=reptiles
    REF_NAME=GCF_030035675.1 
Changed awk in zcat command on line 37 for chromosomes starting with 'NC'
Submitted batch job 3894997

Prepared 20241231_ROH_Calc.sh for GCA_037176765.1.fa.gz Anolis sagrei (Brown anole) (alternate = GCA_037176775.1.fa.gz)
    #SBATCH --array=1-15
    TODAY_DATE=20250120
    CLADE=reptiles
    REF_NAME=GCA_037176765.1 
Submitted batch job 3895680

Prepared 20250101_mm2_alignment.sh for GCA_027789765.1.fa.gz Dendropsophus ebraccatus (hourglass treefrog) (alternate = GCA_027789725.1.fa.gz)
    REF_NAME=GCA_027789765.1 
    CLADE=amphibians
    REFERENCE=../241117.UCSC-hubs-VGP-alignment/alignment/reference/amphibians/GCA_027789765.1.fa.gz ## Finish filling in depending on genome and clade
    ALTERNATE=../241117.UCSC-hubs-VGP-alignment/alignment/alternate/amphibians/GCA_027789725.1.fa.gz ## Finish filling in depending on genome and clade
Submitted batch job 3898002

Job 3895680 completed for ROH calculation of GCA_037176765.1.fa.gz Anolis sagrei (Brown anole) (alternate = GCA_037176775.1.fa.gz)
It worked!

Job 3894997 completed for plotting ROH/calculating FROH for GCF_030035675.1.fa.gz Rhineura floridana (Florida worm lizard) (alternate = GCA_030035715.1.fa.gz)
It worked!

Prepared 20250106_Plot_ROH_FROH.sh for GCA_037176765.1.fa.gz Anolis sagrei (Brown anole) (alternate = GCA_037176775.1.fa.gz)
    TODAY_DATE=20250120
    ROH_CALC_DATE=20250120
    CLADE=reptiles
    REF_NAME=GCA_037176765.1
Changed awk in zcat command on line 37 for chromosomes starting with 'CM'
Submitted batch job 3898304

Prepared 20241231_ROH_Calc.sh for GCA_030412105.1.fa.gz Cyclura pinguis (alternate = GCA_030412085.1.fa.gz)
    #SBATCH --array=1-19 
    TODAY_DATE=20250120
    CLADE=reptiles
    REF_NAME=GCA_030412105.1
Submitted batch job 3898314

Finished writing 20250114_find_het_per_chr.sh and 20250115_find_het_calc.R scripts
These scripts will be used to calculate heterozygosity per 1Mb window along each chromosome, along with mean per chromosome and mean over whole genome

Prepared 20250114_find_het_per_chr.sh for GCF_020745735.1 Hemiscyllium ocellatum (Epaulette shark) 
    TODAY_DATE=20250120
    CLADE=sharks
    REF_NAME=GCF_020745735.1 
Submitted batch job 3911276

Job 3898002 completed for mm2 alingment for GCA_027789765.1.fa.gz Dendropsophus ebraccatus (hourglass treefrog) (alternate = GCA_027789725.1.fa.gz)
It worked!

Job 3898314 completed for ROH calculation of GCA_030412105.1.fa.gz Cyclura pinguis (alternate = GCA_030412085.1.fa.gz)
It worked!

Job 3898304 completed for plotting ROH/calculating FROH for GCA_037176765.1.fa.gz Anolis sagrei (Brown anole) (alternate = GCA_037176775.1.fa.gz)
It worked!

Prepared 20250106_Plot_ROH_FROH.sh for GCA_030412105.1.fa.gz Cyclura pinguis (alternate = GCA_030412085.1.fa.gz)
    TODAY_DATE=20250120
    ROH_CALC_DATE=20250120
    CLADE=reptiles
    REF_NAME=GCA_030412105.1 
awk in zcat command on line 37 for chromosomes starting with 'CM'
Submitted batch job 3910937

Prepared 20241231_ROH_Calc.sh for GCA_030867105.1.fa.gz Liasis olivaceus (alternate = GCA_030867085.1.fa.gz)
    #SBATCH --array=1-17 
    TODAY_DATE=20250120
    CLADE=reptiles
    REF_NAME=GCA_030867105.1
Submitted batch job 3911155

Preapred 20250101_mm2_alignment.sh for GCA_027917425.1.fa.gz Gastrophryne carolinensis (eastern narrow-mouthed toad) (alternate = GCA_027917415.1.fa.gz)
    REF_NAME=GCA_027917425.1
    CLADE=amphibians
    REFERENCE=../241117.UCSC-hubs-VGP-alignment/alignment/reference/amphibians/GCA_027917425.1.fa.gz ## Finish filling in depending on genome and clade
    ALTERNATE=../241117.UCSC-hubs-VGP-alignment/alignment/alternate/amphibians/GCA_027917415.1.fa.gz ## Finish filling in depending on genome and clade
Submitted batch job 3911273

Job 3910937 completed for plotting ROH/calculating FROH for GCA_030412105.1.fa.gz Cyclura pinguis (alternate = GCA_030412085.1.fa.gz)
It worked!

Job 3911276 finished for calculating heterozygosity of GCF_020745735.1 Hemiscyllium ocellatum (Epaulette shark) 
Did not work -- stopped with this error:
    Error in `filter()`:
    ℹ In argument: `data_frame$Chr == chr`.
    Caused by error:
    ! `..1` must be of size 5127864 or 1, not size 0.
Realized column names for the variant position and chromosome length files were not defined
Fixed by adding:
    colnames(dat) <- c('Chr', 'Pos')
    colnames(chr_file) <- c('Chr', 'Length')
Resubmitted
Submitted batch job 3912022

Job 3912022 finished for calculating heterozygosity of GCF_020745735.1 Hemiscyllium ocellatum (Epaulette shark)
Did not worked -- stopped with this error:
    Error in chr_length - window_length : 
    non-numeric argument to binary operator
    Calls: calc_het -> seq -> seq.default
Fixed by making sure to define window_length, window_interval, and chr_length as integers
    current_window_length <- as.integer(args[10])
    current_window_interval <- as.integer(args[11])

    chr_length <- as.integer(chromosome_length_file[i,2])
Resubmitted
Submitted batch job 3912751

Job 3911155 completed for calculating ROH of GCA_030867105.1.fa.gz Liasis olivaceus (alternate = GCA_030867085.1.fa.gz)
It worked!

Prepared 20250106_Plot_ROH_FROH.sh for GCA_030867105.1.fa.gz Liasis olivaceus (alternate = GCA_030867085.1.fa.gz)
    TODAY_DATE=20250120
    ROH_CALC_DATE=20250120
    CLADE=reptiles
    REF_NAME=GCA_030867105.1 
awk in zcat command on line 37 for chromosomes starting with 'CM'
Submitted batch job 3914206

Prepared 20241231_ROH_Calc.sh for GCA_033349115.1.fa.gz Macrochelys suwanniensis (Suwannee snapping turtle) (alternate = GCA_033296515.1.fa.gz)
    #SBATCH --array=1-26 ##Number of chromosomes your species have
    TODAY_DATE=20250120
    CLADE=reptiles
    REF_NAME=GCA_033349115.1 
Submitted batch job 3914455


#### UPDATE ####
20250121 (January 21st, 2025)

Job 3911273 completed for mm2 alignment of GCA_027917425.1.fa.gz Gastrophryne carolinensis (eastern narrow-mouthed toad) (alternate = GCA_027917415.1.fa.gz)
Job had an error, was not successful:
    /var/spool/slurm/slurmd/job3911273/slurm_script: line 35: 2368403 Killed                  minimap2 -t24 -cx asm5 -L --cs ${REFERENCE} ${ALTERNATE} > ${SPECIES}/${REF_NAME}.mm2.paf
    slurmstepd: error: Detected 1 oom_kill event in StepId=3911273.batch. Some of the step tasks have been OOM Killed.
Resulting .paf and .srt.paf files had file size of 0, .vcf had file sie of only 45581 (too small)
Removed incorrect files
    rm -f -r GCA_027917425.1.mm2.paf
    rm -f -r GCA_027917425.1.mm2.srt.paf
    rm -f -r GCA_027917425.1_aligned.mm2.vcf
Double check 20250101_mm2_alignment.sh script -- REF name and alternate name were correct
Everything in the script seemed correct, so I will resubmit
Just in case it was a memory error, I upped the memory for the task to 200GB
Updated dates for .out/.err files
Resubmitted
Submitted batch job 3998270

Job 3914206 completed for plotting ROH/calculating FROH of GCA_030867105.1.fa.gz Liasis olivaceus (alternate = GCA_030867085.1.fa.gz)
It worked!

Job 3914455 completed for ROH calculations of GCA_033349115.1.fa.gz Macrochelys suwanniensis (Suwannee snapping turtle) (alternate = GCA_033296515.1.fa.gz)
It worked!

Prepared 20250106_Plot_ROH_FROH.sh for GCA_033349115.1.fa.gz Macrochelys suwanniensis (Suwannee snapping turtle) (alternate = GCA_033296515.1.fa.gz)
    TODAY_DATE=20250121
    ROH_CALC_DATE=20250120
    CLADE=reptiles
    REF_NAME=GCA_033349115.1
awk in zcat command on line 37 for chromosomes starting with 'CM'
Submitted batch job 3998847

Prepared 20241231_ROH_Calc.sh for GCA_035046505.1.fa.gz Tiliqua scincoides (alternate = GCA_035046495.1.fa.gz)
    #SBATCH --array=1-16 
    TODAY_DATE=20250121
    CLADE=reptiles
    REF_NAME=GCA_035046505.1
Submitted batch job 3998908

Job 3912751 finished for calculating heterozygosity of GCF_020745735.1 Hemiscyllium ocellatum (Epaulette shark)
It did not work, this error was generated:
    Error in seq.default(0, (chr_length - window_length), by = window_interval) : 
    wrong sign in 'by' argument
    Calls: calc_het -> seq -> seq.default
Updated date in .out/.err files and TODAY_DATE
    TODAY_DATE=20250121
Error is in line 42 of 20250115_find_het_calc.R:
    window_starts <- seq(0, chr_length-window_length, by = window_interval)
When I double check was results were generated, heterozygosity .txt files were generated for the first 50 chromosomes, but are missing for CHR 52,X, and Y 
No output files were created with the means for each chromosome or over the whole genome either
The problem is that NC_083452.1 is less than 1Mb in length, so chr_length - window_length results in a negative value
Based on methods used in Stanhope et al. 2023, when they had scaffolds smaller than 1Mb, they used 50kbp sliding windows instead
Expanded on window_starts command to include an if statement to change window size based on chromosome size
            if(chr_length < window_length){
            window_length_alt <- 50000
            window_interval_alt <- window_length_alt/2
            window_starts <- seq(0, (chr_length-window_length_alt), by = window_interval_alt)
            window_ends <- window_starts + window_length_alt
        }
        else{
            window_starts <- seq(0, (chr_length-window_length), by = window_interval)
            window_ends <- window_starts + window_length 
        }
Removed previously generated heterozygosity .txt files
    rm -r -f *.1_heterozygosity.txt
Resubmitted 20250114_find_het_per_chr.sh
Submitted batch job 4000901

Job 3998847 completed for plotting ROH/calculating FROH of  GCA_033349115.1.fa.gz Macrochelys suwanniensis (Suwannee snapping turtle) (alternate = GCA_033296515.1.fa.gz)
It worked!

Job 3998908 completed for calculating ROH of GCA_035046505.1.fa.gz Tiliqua scincoides (alternate = GCA_035046495.1.fa.gz)
It worked!

Preapred 20250106_Plot_ROH_FROH.sh for GCA_035046505.1.fa.gz Tiliqua scincoides (alternate = GCA_035046495.1.fa.gz)
    TODAY_DATE=20250121
    ROH_CALC_DATE=20250121
    CLADE=reptiles
    REF_NAME=GCA_035046505.1
awk in zcat command on line 37 for chromosomes starting with 'CM'
Submitted batch job 4005316

Prepared 20241231_ROH_Calc.sh for GCF_035149785.1.fa.gz Candoia aspera (alternate = GCA_035125265.1.fa.gz)
    #SBATCH --array=1-18
    TODAY_DATE=20250121
    CLADE=reptiles
    REF_NAME=GCF_035149785.1 
Submitted batch job 4005447

Job 4005316 completed for plotting ROH/calculating FROH for GCA_035046505.1.fa.gz Tiliqua scincoides (alternate = GCA_035046495.1.fa.gz)
It worked!

Job 4005447 completed for ROH calculation for GCF_035149785.1.fa.gz Candoia aspera (alternate = GCA_035125265.1.fa.gz)
It worked!

Prepared 20250106_Plot_ROH_FROH.sh for GCF_035149785.1.fa.gz Candoia aspera (alternate = GCA_035125265.1.fa.gz)
    TODAY_DATE=20250121
    ROH_CALC_DATE=20250121
    CLADE=reptiles
    REF_NAME=GCF_035149785.1
Changed awk in zcat command on line 37 for chromosomes starting with 'NC'
Submitted batch job 4010659

Job 4000901 completed for calculating heterozygosity of GCF_020745735.1 Hemiscyllium ocellatum (Epaulette shark)
It worked!
While results were successful, I realized that I needed to change the script so that it reports consistent numbers (e.g. all heterozygosity in #variants/1Mb), and I need to report window sizes for chromosomes
Modified code in 20250115_find_het_calc.R starting at line 42 to include column in results data frame with window sizes specified
        if(chr_length < window_length){
            window_length_alt <- 50000
            window_interval_alt <- window_length_alt/2
            window_starts <- seq(0, (chr_length-window_length_alt), by = window_interval_alt)
            window_ends <- window_starts + window_length_alt
            window_sizes <- rep(window_length_alt, length(window_starts))
        }
        else{
            window_starts <- seq(0, (chr_length-window_length), by = window_interval)
            window_ends <- window_starts + window_length 
            window_sizes <- rep(window_length, length(window_starts))
        }
        het <- rep(0, length(window_starts))
        single_chr_results <- data.frame(window_starts, window_ends, het, window_sizes)
        colnames(single_chr_results) <- c('Start', 'End', 'Het', 'Window_Size')
Added column to dataframe after heterozygosity is calculated to calculate het per 1kb
        single_chr_results$Het_Per_KB <- (single_chr_results$Het/single_chr_results$Window_Size)*1000
Also modified such that 20250121_per_chr_mean_heterozygosity.txt and 20250121_whole_genome_mean_heterozygosity.txt will both have heterozygosity per 1kb 
I am making results in per 1kb to compare to the results of Stanhope et al 2023

Given that sex chromosomes can sometimes have odd results, I modified the .sh script to include the variable NUM_AUTOSOMAL_CHROMOSOMES
This variable will go into the R script and be used to separately calculate the total heterozygosity per kb over the whole autosomal genome, to be reported alongside the het over the whole genome (autosomal + sex)

To run this, while saving results of last run just in case this doesn't work, I created a directory to store the results of this run in
    mkdir test_het_files
    mv *.1_heterozygosity.txt test_het_files/
    mv 20250121_whole_genome_mean_heterozygosity.txt test_het_files/
    mv 20250121_per_chr_mean_heterozygosity.txt test_het_files/

Prepared 20250114_find_het_per_chr.sh for GCF_020745735.1 Hemiscyllium ocellatum (Epaulette shark)  
    NUM_AUTOSOMAL_CHROMOSOMES=52
Submitted batch job 4014276

Job 4010659 completed for plotting ROH/calculating FROH for GCF_035149785.1.fa.gz Candoia aspera (alternate = GCA_035125265.1.fa.gz)
It worked!

Preapred 20241231_ROH_Calc.sh for GCF_025201925.1 Gopherus flavomarginatus (Mexican gopher tortoise) (alternate = GCA_025201965.1.fa.gz)
I needed to recalculate this after the mix-up with the Alligator mississippiensis genome
    #SBATCH --array=1-25 
    TODAY_DATE=20250121
    CLADE=reptiles
    REF_NAME=GCF_025201925.1
Submitted batch job 4015117

Job 4015117 completed for ROH calculation of GCF_025201925.1 Gopherus flavomarginatus (Mexican gopher tortoise) (alternate = GCA_025201965.1.fa.gz)
It worked!

Prepared 20250106_Plot_ROH_FROH.sh for GCF_025201925.1 Gopherus flavomarginatus (Mexican gopher tortoise) (alternate = GCA_025201965.1.fa.gz)
    TODAY_DATE=20250121
    ROH_CALC_DATE=20250121
    CLADE=reptiles
    REF_NAME=GCF_025201925.1
awk in zcat command on line 37 for chromosomes starting with 'NC'
Submitted batch job 4020788


#### UPDATE ####
20250122 (January 22nd, 2025)

Job 4014276 finished for finding heterozygosity of GCF_020745735.1 Hemiscyllium ocellatum (Epaulette shark)
Job did not work due to error:
    Error in calc_het(data_frame = dat, chromosome_length_file = chr_file,  : 
    argument "number_autosomal_chromosomes" is missing, with no default
    Execution halted
Realized I had not put arugment in for the number of autosomal chromosomes when I actually ran the function
Fixed:
    run_function <- calc_het(data_frame=dat, chromosome_length_file=chr_file, window_length=current_window_length, window_interval=current_window_interval, number_autosomal_chromosomes=num_autosomal_chr)
Updated today's date and date for .out/.err files on 20250114_find_het_per_chr.sh
    TODAY_DATE=20250122
Resubmitted
Submitted batch job 4101181

Job 4020788 completed to plot ROH/calculate FROH of  GCF_025201925.1 Gopherus flavomarginatus (Mexican gopher tortoise) (alternate = GCA_025201965.1.fa.gz)
It worked!

Job 3998270 finished for mm2 alignment of GCA_027917425.1.fa.gz Gastrophryne carolinensis (eastern narrow-mouthed toad) (alternate = GCA_027917415.1.fa.gz)
It did not completed, it timed out due to time limit of 12 hours
Given how long this is taking, I might revisit it later after aligning some other genomes -- gives me time to think on what is causing this to take so long
Removed .paf file in case it was incomplete/corrupted
    rm -f -r GCA_027917425.1.mm2.paf

Prepared 20241231_ROH_Calc.sh for GCA_027789765.1.fa.gz Dendropsophus ebraccatus (hourglass treefrog) (alternate = GCA_027789725.1.fa.gz)
    #SBATCH --array=1-15
    TODAY_DATE=20250122
    CLADE=amphibians
    REF_NAME=GCA_027789765.1
Submitted batch job 4102140

Prepared 20250101_mm2_alignment.sh for GCA_031893055.1.fa.gz Leptodactylus fuscus (rufous frog) (alternate = GCA_031893025.1.fa.gz)
    REF_NAME=GCA_031893055.1
    CLADE=amphibians
    REFERENCE=../241117.UCSC-hubs-VGP-alignment/alignment/reference/amphibians/GCA_031893055.1.fa.gz ## Finish filling in depending on genome and clade
    ALTERNATE=../241117.UCSC-hubs-VGP-alignment/alignment/alternate/amphibians/GCA_031893025.1.fa.gz ## Finish filling in depending on genome and clade
Submitted batch job 4102524

Job 4101181 completed for heterozygosity calculation of  GCF_020745735.1 Hemiscyllium ocellatum (Epaulette shark)
It worked!

Prepared 20240114_find_het_per_chr.sh for GCA_030028105.1 Mobula birostris (Giant manta)
    TODAY_DATE=20250122
    CLADE=sharks
    REF_NAME=GCA_030028105.1 
    NUM_AUTOSOMAL_CHROMOSOMES=32 
    WINDOW_LENGTH=1000000
    WINDOW_INTERVAL=500000
Submitted batch job 4108012

Preparing chrom_lists for birds directory for analysis
Prep for GCA_036013475.2.fa.gz Columba livia (rock pigeon) (alternate = GCA_036010775.2.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_036013475.2.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_036013475.2.fa.gz | grep '>CM' | wc ## --> 42 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_036013475.2.fa.gz | grep '>CM' | tr -d '>' > birds/chrom_lists/GCA_036013475.2_chroms.txt

Prep for GCA_037962945.1.fa.gz Sarcoramphus papa (king vulture) (alternate = GCA_037950955.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_037962945.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_037962945.1.fa.gz | grep '>CM' | wc ## --> 42 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_037962945.1.fa.gz | grep '>CM' | tr -d '>' > birds/chrom_lists/GCA_037962945.1_chroms.txt

Prep for GCA_036873955.1.fa.gz Mergus octosetaceus (Brazilian merganser) (alternate = GCA_036850655.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_036873955.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_036873955.1.fa.gz | grep '>CM' | wc ## --> 37 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_036873955.1.fa.gz | grep '>CM' | tr -d '>' > birds/chrom_lists/GCA_036873955.1_chroms.txt

Prep for GCF_036370855.1.fa.gz Dromaius novaehollandiae (emu) (alternate = GCA_036417515.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCF_036370855.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCF_036370855.1.fa.gz | grep '>NC' | wc ## --> 41 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCF_036370855.1.fa.gz | grep '>NC' | tr -d '>' > birds/chrom_lists/GCF_036370855.1_chroms.txt

Prep for GCF_028858725.1.fa.gz Colius striatus (speckled mousebird) (alternate = GCA_028858625.2.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCF_028858725.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCF_028858725.1.fa.gz | grep '>NC' | wc ## --> 33 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCF_028858725.1.fa.gz | grep '>NC' | tr -d '>' > birds/chrom_lists/GCF_028858725.1_chroms.txt

Prep for GCA_034619465.1.fa.gz Leptosomus discolor (cuckoo roller) (alternate = GCA_034619455.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_034619465.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_034619465.1.fa.gz | grep '>CM' | wc ## --> 29 chromosomes based on hap1
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_034619465.1.fa.gz | grep '>CM' | tr -d '>' > birds/chrom_lists/GCA_034619465.1_chroms.txt

Prep for GCF_963924245.1.fa.gz Chroicocephalus ridibundus (black-headed gull) (alternate = GCA_030820635.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCF_963924245.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCF_963924245.1.fa.gz | grep '>NC' | wc ## --> 33 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCF_963924245.1.fa.gz | grep '>NC' | tr -d '>' > birds/chrom_lists/GCF_963924245.1_chroms.txt

Prep for GCA_036013445.1.fa.gz Caloenas nicobarica (Nicobar pigeon) (alternate = GCA_036010745.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_036013445.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_036013445.1.fa.gz | grep '>CM' | wc ## --> 40 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_036013445.1.fa.gz | grep '>CM' | tr -d '>' > birds/chrom_lists/GCA_036013445.1_chroms.txt

Prep for GCA_036169615.1.fa.gz Heliangelus exortis (alternate = GCA_036172335.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_036169615.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_036169615.1.fa.gz | grep '>CM' | wc ## --> 33 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_036169615.1.fa.gz | grep '>CM' | tr -d '>' > birds/chrom_lists/GCA_036169615.1_chroms.txt

Prep for GCA_036417665.1.fa.gz Passer domesticus (house sparrow) (alternate = GCA_036417895.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_036417665.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_036417665.1.fa.gz | grep '>CM' | wc ## --> 39 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_036417665.1.fa.gz | grep '>CM' | tr -d '>' > birds/chrom_lists/GCA_036417665.1_chroms.txt

Prep for GCF_036250125.1.fa.gz Pseudopipra pipra (alternate = GCA_036250135.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCF_036250125.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCF_036250125.1.fa.gz | grep '>NC' | wc ## --> 33 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCF_036250125.1.fa.gz | grep '>NC' | tr -d '>' > birds/chrom_lists/GCF_036250125.1_chroms.txt

Prep for GCA_036417845.1.fa.gz Apteryx mantelli (North Island brown kiwi) (alternate = GCA_036417975.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_036417845.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_036417845.1.fa.gz | grep '>CM' | wc ## --> 44 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_036417845.1.fa.gz | grep '>CM' | tr -d '>' > birds/chrom_lists/GCA_036417845.1_chroms.txt

Prep for GCA_036417535.1.fa.gz Chlamydotis macqueenii (Macqueen's bustard) (alternate = GCA_036418225.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_036417535.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_036417535.1.fa.gz | grep '>CM' | wc ## --> 39 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_036417535.1.fa.gz | grep '>CM' | tr -d '>' > birds/chrom_lists/GCA_036417535.1_chroms.txt


#### UPDATE #### 
20250123 (January 23rd, 2025)

Job 4108012 completed for finding heterozygosity of GCA_030028105.1 Mobula birostris (Giant manta)
It worked!

Prepared 20250114_find_het_per_chr.sh for GCA_035084275.1 Hydrolagus colliei (Pacific ratfish)
    TODAY_DATE=20250123
    CLADE=sharks
    REF_NAME=GCA_035084275.1
    NUM_AUTOSOMAL_CHROMOSOMES=40 
    WINDOW_LENGTH=1000000 ## How long we want each window that we are using to measure heterozygosity
    WINDOW_INTERVAL=500000 
Submitted batch job 4175337

Job 4102140 completed for ROH calculation of GCA_027789765.1.fa.gz Dendropsophus ebraccatus (hourglass treefrog) (alternate = GCA_027789725.1.fa.gz)
It worked!

Prepared 20250106_Plot_ROH_FROH.sh for GCA_027789765.1.fa.gz Dendropsophus ebraccatus (hourglass treefrog) (alternate = GCA_027789725.1.fa.gz)
    TODAY_DATE=20250123
    ROH_CALC_DATE=20250123
    CLADE=amphibians
    REF_NAME=GCA_027789765.1
Changed awk in zcat command on line 37 for chromosomes starting with 'CM'
Submitted batch job 4175716

Contuining to prepare chrom_lists for birds directory for analysis
NOTE: Alternate bird assembly GCA_024206055.2.fa.gz for Gallus Gallus is haploid, has no matching haplotype easily found
NOTE: Alternate bird assembly GCA_008822115.3.fa.gz for Taeniopygia guttata (zebra finch) -- depsite having a linked haplotype, it was not found in reference directory
NOTE: Alternate bird assembly GCA_027557775.1.fa.gz for Gallus Gallus has a linked haplotype, but it cannot be found within the reference directory
NOTE: Alternate bird assembly GCA_027408225.1.fa.gz for Gallus Gallus has a linked haplotype, but it cannot be found within the reference directory

Prep for GCA_031877795.1.fa.gz Strix aluco (Tawny owl) (alternate = GCA_031877785.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_031877795.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_031877795.1.fa.gz | grep '>CM' | wc ## --> 42 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_031877795.1.fa.gz | grep '>CM' | tr -d '>' > birds/chrom_lists/GCA_031877795.1_chroms.txt   

Prep for GCF_012275295.1.fa.gz Melopsittacus undulatus (budgerigar) (alternate = GCA_012275275.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCF_012275295.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCF_012275295.1.fa.gz | grep '>NC' | wc ## --> 32 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCF_012275295.1.fa.gz | grep '>NC' | tr -d '>' > birds/chrom_lists/GCF_012275295.1_chroms.txt   

Prep for GCA_020800305.1.fa.gz Porphyrio hochstetteri (South Island takahe) (alternate = GCA_020801775.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_020800305.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_020800305.1.fa.gz | grep '>CM' | wc ## --> 35 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_020800305.1.fa.gz | grep '>CM' | tr -d '>' > birds/chrom_lists/GCA_020800305.1_chroms.txt   

Prep for GCF_003957565.2.fa.gz Taeniopygia guttata (zebra finch) (alternate = GCF_008822105.2.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCF_003957565.2.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCF_003957565.2.fa.gz | grep '>NC' | wc ## --> 42 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCF_003957565.2.fa.gz | grep '>NC' | tr -d '>' > birds/chrom_lists/GCF_003957565.2_chroms.txt   

Prep for GCF_016700215.2.fa.gz Gallus gallus (chicken) (alternate = GCF_016699485.2.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCF_016700215.2.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCF_016700215.2.fa.gz | grep '>NC' | wc ## --> 41 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCF_016700215.2.fa.gz | grep '>NC' | tr -d '>' > birds/chrom_lists/GCF_016700215.2_chroms.txt   

Prep for GCA_031468815.1.fa.gz Morus bassanus (northern gannet) (alternate = GCA_031468805.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_031468815.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_031468815.1.fa.gz | grep '>CM' | wc ## --> 35 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_031468815.1.fa.gz | grep '>CM' | tr -d '>' > birds/chrom_lists/GCA_031468815.1_chroms.txt   

Prep for GCF_030936135.1.fa.gz Gavia stellata (red-throated loon) (alternate = GCA_030936125.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCF_030936135.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCF_030936135.1.fa.gz | grep '>NC' | wc ## --> 44 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCF_030936135.1.fa.gz | grep '>NC' | tr -d '>' > birds/chrom_lists/GCF_030936135.1_chroms.txt   

Prep for GCA_030867145.1.fa.gz Opisthocomus hoazin (hoatzin) (alternate = GCA_030867165.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_030867145.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_030867145.1.fa.gz | grep '>CM' | wc ## --> 41 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_030867145.1.fa.gz | grep '>CM' | tr -d '>' > birds/chrom_lists/GCA_030867145.1_chroms.txt   

Prep for GCF_030490865.1.fa.gz Poecile atricapillus (Black-capped chickadee) (alternate = GCA_030490855.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCF_030490865.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCF_030490865.1.fa.gz | grep '>NC' | wc ## --> 42 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCF_030490865.1.fa.gz | grep '>NC' | tr -d '>' > birds/chrom_lists/GCF_030490865.1_chroms.txt   

Prep for GCF_017639655.2.fa.gz Falco naumanni (lesser kestrel) (alternate = GCA_017639645.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCF_017639655.2.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCF_017639655.2.fa.gz | grep '>NC' | wc ## --> 28 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCF_017639655.2.fa.gz | grep '>NC' | tr -d '>' > birds/chrom_lists/GCF_017639655.2_chroms.txt

Prep for GCA_028858755.1.fa.gz Ara ararauna (blue-and-yellow macaw) (alternate = GCA_028858555.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_028858755.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_028858755.1.fa.gz | grep '>CM' | wc ## --> 35 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCA_028858755.1.fa.gz | grep '>CM' | tr -d '>' > birds/chrom_lists/GCA_028858755.1_chroms.txt

Prep for GCF_028858705.1.fa.gz Grus americana (Whooping crane) (alternate = GCA_028858595.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCF_028858705.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCF_028858705.1.fa.gz | grep '>NC' | wc ## --> 41 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCF_028858705.1.fa.gz | grep '>NC' | tr -d '>' > birds/chrom_lists/GCF_028858705.1_chroms.txt

Prep for GCF_028500815.1.fa.gz Rissa tridactyla (black-legged kittiwake) (alternate = GCA_028501385.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCF_028500815.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCF_028500815.1.fa.gz | grep '>NC' | wc ## --> 32 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/birds/GCF_028500815.1.fa.gz | grep '>NC' | tr -d '>' > birds/chrom_lists/GCF_028500815.1_chroms.txt

Finished preparing chrom_lists for bird directory

Job 4175337 completed for finding heterozygosity of GCA_035084275.1 Hydrolagus colliei (Pacific ratfish)
It worked!

Prepared 20250114_find_het_per_chr.sh for GCF_030144855.1 Hypanus sabinus (Atlantic stingray)
    TODAY_DATE=20250123
    CLADE=sharks
    REF_NAME=GCF_030144855.1 
    NUM_AUTOSOMAL_CHROMOSOMES=32
    WINDOW_LENGTH=1000000 ## How long we want each window that we are using to measure heterozygosity
    WINDOW_INTERVAL=500000
Submitted batch job 4179456

Created scripts 20250123_Plot_het.sh and 20240123_Plot_het_per_chr.R
These will create a single file containing heterozygosity of all windows across whole genome, and plot it similar to Stanhope et al. 2023 Fig. 2

Job 4179456 completed for finding heterozygosity of GCF_030144855.1 Hypanus sabinus (Atlantic stingray)
It worked!

Prepared 20250114_find_het_per_chr.sh for GCA_035084215.1 Heptranchias perlo (sharpnose sevengill shark)
    TODAY_DATE=20250123
    CLADE=sharks
    REF_NAME=GCA_035084215.1
    NUM_AUTOSOMAL_CHROMOSOMES=45 
    WINDOW_LENGTH=1000000 ## How long we want each window that we are using to measure heterozygosity
    WINDOW_INTERVAL=500000
Submitted batch job 4191586

Job 4175716 finished with errors for plotting ROH/calculating FROH of GCA_027789765.1.fa.gz Dendropsophus ebraccatus (hourglass treefrog) (alternate = GCA_027789725.1.fa.gz)
Errors were all in this same vein:
    awk: fatal: cannot open file `amphibians/GCA_027789765.1/20250123_CM051056.1_100kb_Results_ROH_Durbin_Calc.txt' for reading (No such file or directory)
Error was in the ROH_CALC_DATE -- had it listed as 20250123 when it should've been 20250122
Fixed
    ROH_CALC_DATE=20250122
Resubmitted
Submitted batch job 4191658

Job 4102524 to do mm2 alignment of GCA_031893055.1.fa.gz Leptodactylus fuscus (rufous frog) (alternate = GCA_031893025.1.fa.gz) still hasn't run
Canceled job
Updated date for .out/.err files
Lowered memory usage to 120GB and resubmitted
    #SBATCH --mem=120GB
Submitted batch job 4194558

Job 4191658 finished with errors for plotting ROH/calculating FROH of GCA_027789765.1.fa.gz Dendropsophus ebraccatus (hourglass treefrog) (alternate = GCA_027789725.1.fa.gz)
Error occuring when generating Chroms_Lengths file -- this file is empty
Realized it was because the directory for the reference genoome hadn't been updated from reptiles to amphibians
Fixed by using clade variable instead of typing it out
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/${CLADE}/${REF_NAME}.fa.gz.....
Removed generated files since they are incorrect
    rm -r -f GCA_027789765.1_Chroms_Lengths.txt
    rm -r -f 20250123_GCA_027789765.1_ROH_Results.tsv
    rm -f -r 20250123_GCA_027789765.1_FROH.txt
    rm -f -r 20250123_GCA_027789765.1_ROH_Map.pdf
Resubmitted
Submitted batch job 4194797

Job 4194797 finished for plotting ROH/calculating FROH of GCA_027789765.1.fa.gz Dendropsophus ebraccatus (hourglass treefrog) (alternate = GCA_027789725.1.fa.gz)
It worked!

Job 4191586 finished for finding heterozygosity of GCA_035084215.1 Heptranchias perlo (sharpnose sevengill shark)
It worked!

Prepared 20250114_find_het_per_chr.sh for GCA_036365495.1 Heterodontus francisci (horn shark)
    TODAY_DATE=20250123
    CLADE=sharks
    REF_NAME=GCA_036365495.1
    NUM_AUTOSOMAL_CHROMOSOMES=50
    WINDOW_LENGTH=1000000
    WINDOW_INTERVAL=500000
Submitted batch job 4203062


#### UPDATE ####
20250124 (January 24th, 2025)

Job 4203062 completed for calculating heterozygosity of GCA_036365495.1 Heterodontus francisci (horn shark)
It worked!

Prepared 20250114_find_het_per_chr.sh for GCA_036971175.1 Narcine bancroftii (Spotted torpedo ray/Caribbean electric ray)
    TODAY_DATE=20250124
    CLADE=sharks
    REF_NAME=GCA_036971175.1
    NUM_AUTOSOMAL_CHROMOSOMES=14
    WINDOW_LENGTH=1000000
    WINDOW_INTERVAL=500000
Submitted batch job 4224768

Job 4224768 completed for calculating heterozygosity of GCA_036971175.1 Narcine bancroftii (Spotted torpedo ray/Caribbean electric ray)
It worked!

Prepared 20250114_find_het_per_chr.sh for GCF_030867095.1 Alligator mississippiensis (American alligator) (alternate = GCA_030867065.1.fa.gz)
    TODAY_DATE=20250124
    CLADE=reptiles
    REF_NAME=GCF_030867095.1
    NUM_AUTOSOMAL_CHROMOSOMES=16
    WINDOW_LENGTH=1000000
    WINDOW_INTERVAL=500000 
Submitted batch job 4241621

Job 4194558 submitted yesterday for mm2 alingment of GCA_031893055.1.fa.gz Leptodactylus fuscus (rufous frog) (alternate = GCA_031893025.1.fa.gz) still hasn't run
Canceled the job
    scancel 4194558
Made some adjustments to the script -- including changing the node to cclake from icelake-himem, lowering memory requirement to 100GB, and lowering time to 10hrs from 12hrs
    #SBATCH -J mm2_alignment
    #SBATCH -p cclake
    #SBATCH --mem=100GB 
    #SBATCH --time=10:00:00
Submitted batch job 4242061

Preparing chrom_lists for mammals directory for analysis
Prep for GCF_003369695.1.fa.gz Bos indicus x Bos taurus (hybrid cattle) (alternate = GCA_003369685.2.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_003369695.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_003369695.1.fa.gz | grep '>NC' | wc ## --> 30 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_003369695.1.fa.gz | grep '>NC' | tr -d '>' > mammals/chrom_lists/GCF_003369695.1_chroms.txt

Prep for GCF_011762595.1.fa.g Tursiops truncatus (common bottlenose dolphin) (alternate = GCA_011762515.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_011762595.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_011762595.1.fa.gz | grep '>NC' | wc ## --> 23 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_011762595.1.fa.gz | grep '>NC' |  tr -d '>' > mammals/chrom_lists/GCF_011762595.1_chroms.txt

Prep for GCF_922984935.1.fa.gz Meles meles (Eurasian badger) (alternate = GCA_922990625.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_922984935.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_922984935.1.fa.gz | grep '>NC' | wc ## --> 24 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_922984935.1.fa.gz | grep '>NC' | tr -d '>' > mammals/chrom_lists/GCF_922984935.1_chroms.txt

Prep for GCF_026419965.1.fa.gz Kogia breviceps (pygmy sperm whale) (alternate = GCA_026419985.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_026419965.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_026419965.1.fa.gz | grep '>NC' | wc ## --> 22 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_026419965.1.fa.gz | grep '>NC' | tr -d '>' > mammals/chrom_lists/GCF_026419965.1_chroms.txt

Prep for GCF_030028045.1.fa.gz Hippopotamus amphibius kiboko (alternate = GCA_030028035.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_030028045.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_030028045.1.fa.gz | grep '>NC' | wc ## --> 19 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_030028045.1.fa.gz | grep '>NC' | tr -d '>' > mammals/chrom_lists/GCF_030028045.1_chroms.txt

Prep for GCF_028564815.1.fa.gz Eubalaena glacialis (North Atlantic right whale) (alternate = GCA_028571275.1.fa.gz) 
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_028564815.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_028564815.1.fa.gz | grep '>NC' | wc  ## --> 23 chromosomes     
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_028564815.1.fa.gz | grep '>NC' | tr -d '>' > mammals/chrom_lists/GCF_028564815.1_chroms.txt

Prep for GCF_028023285.1.fa.gz Balaenoptera ricei (Rice's whale) (alternate = GCA_028017805.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_028023285.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_028023285.1.fa.gz | grep '>NC' | wc ## --> 22 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_028023285.1.fa.gz | grep '>NC' | tr -d '>' > mammals/chrom_lists/GCF_028023285.1_chroms.txt

Prep for GCF_030020395.1.fa.gz Manis pentadactyla (Chinese pangolin) (alternate = GCA_030020945.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_030020395.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_030020395.1.fa.gz | grep '>NC' | wc ## --> 21 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_030020395.1.fa.gz | grep '>NC' | tr -d '>' > mammals/chrom_lists/GCF_030020395.1_chroms.txt

Prep for GCF_030014295.1.fa.gz Loxodonta africana (African savanna elephant) (alternate = GCA_030020305.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_030014295.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_030014295.1.fa.gz | grep '>NC' | wc ## --> 29 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_030014295.1.fa.gz | grep '>NC' | tr -d '>' > mammals/chrom_lists/GCF_030014295.1_chroms.txt

Prep for GCA_030035585.1.fa.gz Dugong dugon (dugong) (alternate = GCA_030020955.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCA_030035585.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCA_030035585.1.fa.gz | grep '>CM' | wc
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCA_030035585.1.fa.gz | grep '>CM' | tr -d '>' > mammals/chrom_lists/GCA_030035585.1_chroms.txt

Prep for GCF_030435755.1.fa.gz Ochotona princeps (American pika) (alternate = GCA_030435715.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_030435755.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_030435755.1.fa.gz | grep '>NC' | wc ## --> 35 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_030435755.1.fa.gz | grep '>NC' | tr -d '>' > mammals/chrom_lists/GCF_030435755.1_chroms.txt

Job 4242061 completed for mm2 alignment of GCA_031893055.1.fa.gz Leptodactylus fuscus (rufous frog) (alternate = GCA_031893025.1.fa.gz)
It worked!

Prepared 20250101_mm2_alignment.sh for GCF_029499605.1.fa.gz Hyla sarda (Sardinian treefrog) (alternate = GCA_029493135.1.fa.gz)
    REF_NAME=GCF_029499605.1
    CLADE=amphibians
    REFERENCE=../241117.UCSC-hubs-VGP-alignment/alignment/reference/amphibians/GCF_029499605.1.fa.gz ## Finish filling in depending on genome and clade
    ALTERNATE=../241117.UCSC-hubs-VGP-alignment/alignment/alternate/amphibians/GCA_029493135.1.fa.gz ## Finish filling in depending on genome and clade
Submitted batch job 4251432

Prepared 20241231_ROH_Calc.sh for GCA_031893055.1.fa.gz Leptodactylus fuscus (rufous frog) (alternate = GCA_031893025.1.fa.gz)
    #SBATCH --array=1-11
    TODAY_DATE=20250124
    CLADE=amphibians
    REF_NAME=GCA_031893055.1 
Submitted batch job 4251504

Job 4241621 completed for calculating heterozygosity of GCF_030867095.1 Alligator mississippiensis (American alligator) (alternate = GCA_030867065.1.fa.gz)
It worked!

Prepared 20250114_find_het_per_chr.sh for GCF_027887155.1 Malaclemys terrapin pileata (Diamondback terrapin) (alternate = GCA_027887205.1.fa.gz)
    TODAY_DATE=20250124
    CLADE=reptiles
    REF_NAME=GCF_027887155.1
    NUM_AUTOSOMAL_CHROMOSOMES=25
    WINDOW_LENGTH=1000000
    WINDOW_INTERVAL=500000 
Submitted batch job 4252138

Continued preparing chrom_lists for non-primate mammals
NOTE: Alternate haplotype GCA_024803745.1.fa.gz Thomomys bottae (Botta's pocket gopher) has no matching reference haplotype in the reference/mammals/ directory

Prep for GCA_031878705.1.fa.gz Tapirus indicus (Asiatic tapir) (alternate = GCA_031878655.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCA_031878705.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCA_031878705.1.fa.gz | grep '>CM' | wc ## --> 27 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCA_031878705.1.fa.gz | grep '>CM' | tr -d '>' > mammals/chrom_lists/GCA_031878705.1_chroms.txt

Prep for GCF_020826845.1.fa.gz Diceros bicornis minor (alternate = GCA_020826835.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_020826845.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_020826845.1.fa.gz | grep '>NC' | wc ## --> 43 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_020826845.1.fa.gz | grep '>NC' | tr -d '>' > mammals/chrom_lists/GCF_020826845.1_chroms.txt

Prep for GCF_020740685.1.fa.gz Jaculus jaculus (lesser Egyptian jerboa) (alternate = GCA_020740715.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_020740685.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_020740685.1.fa.gz | grep '>NC' | wc ## --> 26 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_020740685.1.fa.gz | grep '>NC' | tr -d '>' > mammals/chrom_lists/GCF_020740685.1_chroms.txt

Prep for GCA_031878675.1.fa.gz Thomomys bottae (Botta's pocket gopher) (alternate = GCA_031878665.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCA_031878675.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCA_031878675.1.fa.gz | grep '>CM' | wc ## --> 41 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCA_031878675.1.fa.gz | grep '>CM' | tr -d '>' > mammals/chrom_lists/GCA_031878675.1_chroms.txt

Job 4251504 completed for calculating ROH of GCA_031893055.1.fa.gz Leptodactylus fuscus (rufous frog) (alternate = GCA_031893025.1.fa.gz)
It worked!

Prepared 20250106_Plot_ROH_FROH.sh for GCA_031893055.1.fa.gz Leptodactylus fuscus (rufous frog) (alternate = GCA_031893025.1.fa.gz)
    TODAY_DATE=20250124
    ROH_CALC_DATE=20250124
    CLADE=amphibians
    REF_NAME=GCA_031893055.1
awk in zcat command on line 37 for chromosomes starting with 'CM'
Submitted batch job 4262742

Job 4252138 completed for calculating heterozygosity of GCF_027887155.1 Malaclemys terrapin pileata (Diamondback terrapin) (alternate = GCA_027887205.1.fa.gz)
It worked!

Prepared 20250114_find_het_per_chr.sh for GCA_028017835.1 Emys orbicularis (European pond turtle) (alternate = GCA_028017845.1.fa.gz)
    TODAY_DATE=20250124
    CLADE=reptiles
    REF_NAME=GCA_028017835.1
    NUM_AUTOSOMAL_CHROMOSOMES=25
    WINDOW_LENGTH=1000000
    WINDOW_INTERVAL=500000
Submitted batch job 4263173


#### UPDATE ####
20250126 (January 26th, 2025)

Job 4262742 completed for plotting ROH/calculating FROH of GCA_031893055.1.fa.gz Leptodactylus fuscus (rufous frog) (alternate = GCA_031893025.1.fa.gz)
It worked!

Job 4263173 finished for finding heterozygosity of GCA_028017835.1 Emys orbicularis (European pond turtle) (alternate = GCA_028017845.1.fa.gz)
It did not work, error:
    Error in names(x) <- value : 
        'names' attribute [2] must be the same length as the vector [1]
    Calls: colnames<-
When checking the GCA_028017835.1 directory, I also found that the plotting ROH hadn't worked for this species
GCA_028017835.1_Chroms_Lengths.txt was only 11bytes in size
Both FROH .txt file and ROH Map .pdf file were 0bytes in size

Preparing 20250106_Plot_ROH_FROH.sh for GCA_028017835.1 Emys orbicularis (European pond turtle) (alternate = GCA_028017845.1.fa.gz)
    TODAY_DATE=20250126
    ROH_CALC_DATE=20250119
    CLADE=reptiles
    REF_NAME=GCA_028017835.1 
Removed the files generated from the last run of the code to avoid any issues
    rm -f -r GCA_028017835.1_Chroms_Lengths.txt
    rm -f -r 20250119_GCA_028017835.1_FROH.txt
    rm -f -r 20250119_GCA_028017835.1_ROH_Map.pdf
    rm -f -r 20250119_GCA_028017835.1_ROH_Results.tsv
Resubmitted
Submitted batch job 4364366

Job 4251432 stopped of mm2 alignment of GCF_029499605.1.fa.gz Hyla sarda (Sardinian treefrog) (alternate = GCA_029493135.1.fa.gz)
Did not completed due to an out of memory error:
    slurmstepd: error: Detected 1 oom_kill event in StepId=4251432.batch. Some of the step tasks have been OOM Killed.
Upped memory to 120GB
    #SBATCH --mem=120GB
Removed previously generated files
    rm -f -r GCF_029499605.1.mm2.paf
    rm -f -r GCF_029499605.1.mm2.srt.paf
    rm -f -r GCF_029499605.1_aligned.mm2.vcf
Updated date for .err/.out files to today
Resubmitted
Submitted batch job 4368184

Job 4364366 completed for plotting ROH/calculating FROH of GCA_028017835.1 Emys orbicularis (European pond turtle) (alternate = GCA_028017845.1.fa.gz)
It worked!
Checked generated files, they do not look corrupted

Preparing 20250114_find_het_per_chr.sh for GCA_028017835.1 Emys orbicularis (European pond turtle) (alternate = GCA_028017845.1.fa.gz)
    TODAY_DATE=20250126
    CLADE=reptiles
    REF_NAME=GCA_028017835.1 
    NUM_AUTOSOMAL_CHROMOSOMES=25 
    WINDOW_LENGTH=1000000 
    WINDOW_INTERVAL=500000
Submitted batch job 4368677


#### UPDATE ####
20250127 (January 27th, 2025)

Job 4251432 failed for mm2 alingment of GCF_029499605.1.fa.gz Hyla sarda (Sardinian treefrog) (alternate = GCA_029493135.1.fa.gz)
    slurmstepd: error: Detected 1 oom_kill event in StepId=4368184.batch. Some of the step tasks have been OOM Killed.
Upped memory to 180GB
    #SBATCH --mem=180GB
Removed files generated from previous run
    rm -f -r GCF_029499605.1.mm2.paf
    rm -f -r GCF_029499605.1.mm2.srt.paf
    rm -f -r GCF_029499605.1_aligned.mm2.vcf
Updated date for .out/.err files
Resubmitted
Submitted batch job 4400110

Job 4368677 completed for finding heterozygosity of GCA_028017835.1 Emys orbicularis (European pond turtle) (alternate = GCA_028017845.1.fa.gz)
It worked!

Prepared 20250114_find_het_per_chr.sh for GCF_029931775.1.fa.gz Euleptes europaea (tarantolino) (alternate = GCA_029931755.1.fa.gz)
    TODAY_DATE=20250127
    CLADE=reptiles
    REF_NAME=GCF_029931775.1
    NUM_AUTOSOMAL_CHROMOSOMES=21 
    WINDOW_LENGTH=1000000
    WINDOW_INTERVAL=500000
Submitted batch job 4404896

Continued preparing chrom_lists for non-primate mammals
Prep for GCF_036321535.1.fa.gz Camelus dromedarius (Arabian camel) (alternate = GCA_036321565.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_036321535.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_036321535.1.fa.gz | grep '>NC' | wc ## --> 38 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_036321535.1.fa.gz | grep '>NC' | tr -d '>' > mammals/chrom_lists/GCF_036321535.1_chroms.txt

Prep for GCF_022682495.1.fa.gz Desmodus rotundus (common vampire bat) (alternate = GCA_022682395.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_022682495.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_022682495.1.fa.gz | grep '>NC' | wc ## --> 15 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_022682495.1.fa.gz | grep '>NC' | tr -d '>' > mammals/chrom_lists/GCF_022682495.1_chroms.txt

Prep for GCA_036426135.1.fa.gz Equus caballus (horse) (alternate = GCA_036418255.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCA_036426135.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCA_036426135.1.fa.gz | grep '>CM' | wc ## --> 34 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCA_036426135.1.fa.gz | grep '>CM' | tr -d '>' > mammals/chrom_lists/GCA_036426135.1_chroms.txt

Prep for GCA_036417435.1.fa.gz Inia geoffrensis (boutu) (alternate = GCA_036417475.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCA_036417435.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCA_036417435.1.fa.gz | grep '>CM' | wc ## --> 24 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCA_036417435.1.fa.gz | grep '>CM' | tr -d '>' > mammals/chrom_lists/GCA_036417435.1_chroms.txt

Prep for GCA_037038545.1.fa.gz Rhynchonycteris naso (proboscis bat) (alternate = GCA_037038555.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCA_037038545.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCA_037038545.1.fa.gz | grep '>CM' | wc ## --> 13 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCA_037038545.1.fa.gz | grep '>CM' | tr -d '>' > mammals/chrom_lists/GCA_037038545.1_chroms.txt

Prep for GCF_030445035.1.fa.gz Dasypus novemcinctus (nine-banded armadillo) (alternate = GCA_030445055.2.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_030445035.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_030445035.1.fa.gz | grep '>NC' | wc ## --> 32 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCF_030445035.1.fa.gz | grep '>NC' | tr -d '>' > mammals/chrom_lists/GCF_030445035.1_chroms.txt

Prep for GCA_037038515.1.fa.gz Microtus pennsylvanicus (meadow vole) (alternate = GCA_037039175.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCA_037038515.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCA_037038515.1.fa.gz | grep '>CM' | wc ## --> 25 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCA_037038515.1.fa.gz | grep '>CM' | tr -d '>' > mammals/chrom_lists/GCA_037038515.1_chroms.txt

Prep for GCA_037157525.1.fa.gz Molossus alvarezi (Alvarez's mastiff bat) (alternate = GCA_037176705.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCA_037157525.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCA_037157525.1.fa.gz | grep '>CM' | wc ## --> 25 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCA_037157525.1.fa.gz | grep '>CM' | tr -d '>' > mammals/chrom_lists/GCA_037157525.1_chroms.txt

Prep for GCA_038363145.1.fa.gz Artibeus intermedius (alternate = GCA_038363225.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCA_038363145.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCA_038363145.1.fa.gz | grep '>CM' | wc ## --> 17 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/mammals/GCA_038363145.1.fa.gz | grep '>CM' | tr -d '>' > mammals/chrom_lists/GCA_038363145.1_chroms.txt
## Finished chroms list for 241117.UCSC-hubs-VGP-alignment/alignment/alternate/mammals directory

Preparing chrom_lists for primates
Prep for GCA_016695395.2.fa.gz Homo sapiens (human) (alternate = GCA_016700455.2.fa.gz)
NOTE: BOTH OF THESE GENOMES WERE FOUND IN THE ALTERNATE DIRECTORY
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/alternate/primates/GCA_016695395.2.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/alternate/primates/GCA_016695395.2.fa.gz | grep '>CM' | wc ## --> 24 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/alternate/primates/GCA_016695395.2.fa.gz | grep '>CM' | tr -d '>' > primates/chrom_lists/GCA_016695395.2_chroms.txt

Prep for GCF_011100555.1.fa.gz Callithrix jacchus (white-tufted-ear marmoset) (alternate = GCA_011078405.1.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/primates/GCF_011100555.1.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/primates/GCF_011100555.1.fa.gz | grep '>NC' | wc ## --> 24 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/primates/GCF_011100555.1.fa.gz | grep '>NC' | tr -d '>' > primates/chrom_lists/GCF_011100555.1_chroms.txt

Prep for GCF_028878055.2.fa.gz Symphalangus syndactylus (siamang) (alternate = GCA_028878085.2.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/primates/GCF_028878055.2.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/primates/GCF_028878055.2.fa.gz | grep '>NC' | wc ## --> 27 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/primates/GCF_028878055.2.fa.gz | grep '>NC' | tr -d '>' > primates/chrom_lists/GCF_028878055.2_chroms.txt

Prep for GCF_028885655.2.fa.gz Pongo abelii (Sumatran orangutan) (alternate = GCA_028885685.2.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/primates/GCF_028885655.2.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/primates/GCF_028885655.2.fa.gz | grep '>NC' | wc ## --> 26 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/primates/GCF_028885655.2.fa.gz | grep '>NC' | tr -d '>' > primates/chrom_lists/GCF_028885655.2_chroms.txt

Prep for GCF_028885625.2.fa.gz Pongo pygmaeus (Bornean orangutan) (alternate = GCA_028885525.2.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/primates/GCF_028885625.2.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/primates/GCF_028885625.2.fa.gz | grep '>NC' | wc ## --> 26 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/primates/GCF_028885625.2.fa.gz | grep '>NC' | tr -d '>' > primates/chrom_lists/GCF_028885625.2_chroms.txt

Prep for GCF_028858775.2.fa.gz Pan troglodytes (chimpanzee) (alternate = GCA_028858805.2.fa.gz)
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/primates/GCF_028858775.2.fa.gz | grep '>'
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/primates/GCF_028858775.2.fa.gz | grep '>NC' | wc ## --> 26 chromosomes
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/primates/GCF_028858775.2.fa.gz | grep '>NC' | tr -d '>' > primates/chrom_lists/GCF_028858775.2_chroms.txt

Moved primates directory to within mammals to match Durbin's wants for organization
    mv primates mammals

Job 4404896 completed for calculating heterozygosity of GCF_029931775.1.fa.gz Euleptes europaea (tarantolino) (alternate = GCA_029931755.1.fa.gz)
It worked!

Prepared 20250114_find_het_per_chr.sh for GCA_030020295.1.fa.gz Gavialis gangeticus (Gharial) (altnerate = GCA_030020385.1.fa.gz)
    TODAY_DATE=20250127
    CLADE=reptiles
    REF_NAME=GCA_030020295.1 
    NUM_AUTOSOMAL_CHROMOSOMES=16
    WINDOW_LENGTH=1000000
    WINDOW_INTERVAL=500000 
Submitted batch job 4411973


#### UPDATE ####
20250128 (January 28th, 2025)

Job 4400110 failed for mm2 alingment of GCF_029499605.1.fa.gz Hyla sarda (Sardinian treefrog) (alternate = GCA_029493135.1.fa.gz)
It failed due to timing out after 10 hours
It only generated the .paf file
Removed the incomplete .paf file
    rm -f -r GCF_029499605.1.mm2.paf

Prepared 20250101_mm2_alignment.sh for GCA_038048845.1.fa.gz Mixophyes fleayi (Fleay's barred frog) (alternate = GCA_038048865.1.fa.gz)
    REF_NAME=GCA_038048845.1
    CLADE=amphibians
    REFERENCE=../241117.UCSC-hubs-VGP-alignment/alignment/reference/amphibians/GCA_038048845.1.fa.gz
    ALTERNATE=../241117.UCSC-hubs-VGP-alignment/alignment/alternate/amphibians/GCA_038048865.1.fa.gz
Submitted batch job 4438986

Job 4411973 completed for calculating heterozygosity of GCA_030020295.1.fa.gz Gavialis gangeticus (Gharial) (altnerate = GCA_030020385.1.fa.gz)
It worked!

Prepared 20250114_find_het_per_chr.sh for GCF_030035675.1.fa.gz Rhineura floridana (Florida worm lizard) (alternate = GCA_030035715.1.fa.gz)
    TODAY_DATE=20250128
    CLADE=reptiles
    REF_NAME=GCF_030035675.1
    NUM_AUTOSOMAL_CHROMOSOMES=22 
    WINDOW_LENGTH=1000000 
    WINDOW_INTERVAL=500000
Submitted batch job 4439512

Job 4439512 completed for calculating heterozygosity of GCF_030035675.1.fa.gz Rhineura floridana (Florida worm lizard) (alternate = GCA_030035715.1.fa.gz)
It worked!

Prepared 20250114_find_het_per_chr.sh for GCA_030412105.1.fa.gz Cyclura pinguis (alternate = GCA_030412085.1.fa.gz)
    TODAY_DATE=20250128
    CLADE=reptiles
    REF_NAME=GCA_030412105.1 
    NUM_AUTOSOMAL_CHROMOSOMES=17
    WINDOW_LENGTH=1000000 
    WINDOW_INTERVAL=500000 
Submitted batch job 4443068

Job 4438986 completed for mm2 alignment of GCA_038048845.1.fa.gz Mixophyes fleayi (Fleay's barred frog) (alternate = GCA_038048865.1.fa.gz)
It worked!

Prepared 20241231_ROH_Calc.sh for GCA_038048845.1.fa.gz Mixophyes fleayi (Fleay's barred frog) (alternate = GCA_038048865.1.fa.gz)
    #SBATCH --array=1-12
    TODAY_DATE=20250128
    CLADE=amphibians
    REF_NAME=GCA_038048845.1
Submitted batch job 4449142

Job 4449142 completed for ROH calculation of GCA_038048845.1.fa.gz Mixophyes fleayi (Fleay's barred frog) (alternate = GCA_038048865.1.fa.gz)
It worked!

Prepared 20250106_Plot_ROH_FROH.sh for GCA_038048845.1.fa.gz Mixophyes fleayi (Fleay's barred frog) (alternate = GCA_038048865.1.fa.gz)
    TODAY_DATE=20250128
    ROH_CALC_DATE=20250128
    CLADE=amphibians
    REF_NAME=GCA_038048845.1 
awk in zcat command on line 37 for chromosomes starting with 'CM'
Submitted batch job 4449494


#### UPDATE ####
20250129 (January 29th, 2025)

Job 4443068 completed for calculating heterozygosity of GCA_030412105.1.fa.gz Cyclura pinguis (alternate = GCA_030412085.1.fa.gz)
It worked!

Prepare 20250114_find_het_per_chr.sh for GCA_030867105.1.fa.gz Liasis olivaceus (alternate = GCA_030867085.1.fa.gz)
    TODAY_DATE=20250129
    CLADE=reptiles
    REF_NAME=GCA_030867105.1
    NUM_AUTOSOMAL_CHROMOSOMES=17
    WINDOW_LENGTH=1000000
    WINDOW_INTERVAL=500000
Submitted batch job 4473737
Canceled job due to it taking so long in the queue 
Switched to cclake node
    #SBATCH -p cclake
Resubmitted
Submitted batch job 4495634

Job 4449494 completed for plotting ROH/calculating FROH of GCA_038048845.1.fa.gz Mixophyes fleayi (Fleay's barred frog) (alternate = GCA_038048865.1.fa.gz)
It worked!

Created 20250129_FASTGA_alignment.sh script to do FASTGA alignments of reference and alternate haplotypes in place of minimap2
    Pulled base FastGA command from hpc-work/slurm_submit.peta4-icelake-fastga-2 script
Prepared script to run test with GCF_020745735.1 Hemiscyllium ocellatum (Epaulette shark)
Submitted batch job 4481376

Job 4481376 for fastGA alingment failed
    FAtoGDB: Cannot open ../241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/.GCF_020745735.1.bps for 'w+'
    FastGA: Call to FAtoGDB failed
Realized error is because I didn't specify file name in the reference and alternate variables
Fixed
    REFERENCE=../241117.UCSC-hubs-VGP-alignment/alignment/reference/${CLADE}/${REF_NAME}.fa.gz 
    ALTERNATE=../241117.UCSC-hubs-VGP-alignment/alignment/alternate/${CLADE}/${ALT_NAME}.fa.gz
Resubmitted
Submitted batch job 4482246

Job 4482246 failed for fastGA alignment
Same error repeated:
    FAtoGDB: Cannot open ../241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/.GCF_020745735.1.bps for 'w+'
    FastGA: Call to FAtoGDB failed
Created HemOce.1aln with a size of 0 bytes
Removed generated file and directory to avoid issues in further runs
    rm -f -r HemOce.1aln
    rm -f -r HemOce
Slightly altered syntax of FastGA command
    FastGA -P ${TEMP} -1: ${SPECIES} -v ${REFERENCE} ${ALTERNATE}
Resubmitted
Submitted batch job 4482473

Job 4482473 failed for fastGA alingment
Same error as before:
    FAtoGDB: Cannot open ../241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/.GCF_020745735.1.bps for 'w+'
    FastGA: Call to FAtoGDB failed

Tried pulling out individual FAtoGDB command and running it in the command line 
    FAtoGDB -v ../241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/GCF_020745735.1.fa.gz 20250129_HemOce_GCF_020745735.1.1gdb
It worked!
Removed generated file from heterozygosity directory
Going to alter 20250129_FASTGA_alignment.sh script to run sub-functions instead of FastGA command
Finsished and submitted
Submitted batch job 4488958

Job 4488958 completed but failed
Error:
    ONEcode file open error sharks/HemOce/20250129_HemOce_ALN.1aln: file is empty
    ALNtoPAF: Failed to open .1aln file sharks/HemOce/20250129_HemOce_ALN.1aln
Removed empty .1aln file
    rm -f -r HemOce.1aln
Error came from two places -- the gix files were not made and I had a syntax error in the FastGA command
Fixed -- changed directory for files to be placed in GIXmake command
    GIXmake -v -P. ${TEMP}/${TODAY_DATE}_${SPEC_NAME}_${REF_NAME}.1gdb ${REFERENCE} ${TEMP}/${TODAY_DATE}_${SPEC_NAME}_${REF_NAME}.gix
    GIXmake -v -P. ${TEMP}/${TODAY_DATE}_${SPEC_NAME}_${ALT_NAME}.1gdb ${ALTERNATE} ${TEMP}/${TODAY_DATE}_${SPEC_NAME}_${ALT_NAME}.gix
    FastGA -v -P. -1:${SPECIES}/${TODAY_DATE}_${SPEC_NAME}_ALN ${REFERENCE} ${TEMP}/${TODAY_DATE}_${SPEC_NAME}_${REF_NAME}.gix ${ALTERNATE} ${TEMP}/${TODAY_DATE}_${SPEC_NAME}_${ALT_NAME}.gix
Resubmitted
Submitted batch job 4495346

Job 4495634 completed for calculating heterozygosity of  GCA_030867105.1.fa.gz Liasis olivaceus (alternate = GCA_030867085.1.fa.gz)
It worked!

Prepared 20250114_find_het_per_chr.sh for GCA_033349115.1.fa.gz Macrochelys suwanniensis (Suwannee snapping turtle) (alternate = GCA_033296515.1.fa.gz)
    TODAY_DATE=20250129
    CLADE=reptiles
    REF_NAME=GCA_033349115.1 
    NUM_AUTOSOMAL_CHROMOSOMES=26
    WINDOW_LENGTH=1000000 
    WINDOW_INTERVAL=500000
Submitted batch job 4500134

Job 4500134 completed for calculating heterozygosity of GCA_033349115.1.fa.gz Macrochelys suwanniensis (Suwannee snapping turtle) (alternate = GCA_033296515.1.fa.gz)
It worked!

Prepared 20250114_find_het_per_chr.sh for GCF_025201925.1 Gopherus flavomarginatus (Mexican gopher tortoise) (alternate = GCA_025201965.1.fa.gz)
    TODAY_DATE=20250129
    CLADE=reptiles
    REF_NAME=GCF_025201925.1
    NUM_AUTOSOMAL_CHROMOSOMES=25
    WINDOW_LENGTH=1000000
    WINDOW_INTERVAL=500000
Submitted batch job 4502171


#### UPDATE ####
20250131 (January 25th, 2025)

Job 4495346 for FASTGA alignment finished but failed
Same error as before
    ONEcode file open error sharks/HemOce/20250129_HemOce_ALN.1aln: file is empty
    ALNtoPAF: Failed to open .1aln file sharks/HemOce/20250129_HemOce_ALN.1aln
I think the GIXmake commands are not running
Will run for reference GIX manually in command line to see if command workspace
    GIXmake -v sharks/HemOce/temp/20250129_HemOce_GCF_020745735.1.1gdb sharks/HemOce/temp/20250131_HemOce_GCF_020745735.1.gix
Received this error 
    GIXmake: Cannot create a .gix with a different location and root name than its .gdb
Went into the sharks/HemOce/temp directory and changed command
    GIXmake -v 20250129_HemOce_GCF_020745735.1.1gdb
Command ran, but then failed due to this error:
    GIXmake: IO write to file /tmp/.post.2241996.17.ktb failed
Tried to specify temp directory
    GIXmake -v -P. 20250129_HemOce_GCF_020745735.1.1gdb
New error arose:
    GIXmake: Cannot open part file ./.20250129_HemOce_GCF_020745735.1.ktab.8 for writing
Tried making it with a new temp directory from heterozygosity directory
    mkdir tmp
    GIXmake -v sharks/HemOce/temp/20250129_HemOce_GCF_020745735.1.1gdb
Failed with old error
    GIXmake: IO write to file /tmp/.post.2504976.22.ktb failed

    GIXmake -v ../241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/GCF_020745735.1.fa.gz 20250131_HemOce_reference
    GIXmake: IO write to file /tmp/.post.2765658.19.ktb failed
    GIXmake: IO write to file /tmp/.post.2765658.21.ktb failed


#### UPDATE ####
20250203 (February 3rd, 2025)

Removed .1gdb files to see if the problem was something related to them
    rm -f -r 20250129_HemOce_GCF_020745735.1.1gdb
    rm -f -r 20250129_HemOce_GCA_020745765.1.1gdb
Re-made the .1gdb files
    FAtoGDB -v ../241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/GCF_020745735.1.fa.gz sharks/HemOce/temp/20250203_HemOce_GCF_020745735.1.1gdb
.1gdb and associated .bps file were created for reference haplotype
    FAtoGDB -v ../241117.UCSC-hubs-VGP-alignment/alignment/alternate/sharks/GCA_020745765.1.fa.gz sharks/HemOce/temp/20250203_HemOce_GCA_020745765.1.1gdb
.1gdb and associated .bps file were created for alternate haplotype

Ran GIXmake command to see if it works creating GIX from .1gdb
    GIXmake -v -P./tmp/ sharks/HemOce/temp/20250203_HemOce_GCF_020745735.1.1gdb
Came out with error:
    GIXmake: Cannot open part file sharks/HemOce/temp/.20250203_HemOce_GCF_020745735.1.ktab.8 for writing
Will run it with temporary directory the same as the one for .1gdb
    GIXmake -v -P./sharks/HemOce/temp/ sharks/HemOce/temp/20250203_HemOce_GCF_020745735.1.1gdb
Got error for saying it was out of memory
Submitting 20250129_FASTGA_alignment.sh script with updated TODAY_DATE and only commands for GIXmake running
Submitted batch job 4628464

Job 4628464 didn't work - due to syntax error in GIXmake commands
Fixed and resubmitted
Submitted batch job 4628853

Job 4628853 failed due to saying that gix and .1gdb files would be in different locations
Changed syntax of temporary directory and resubmitted
Submitted batch job 4630084

Ran command within sharks/HemOce/temp direcotry
GIXmake -v -P. 20250203_HemOce_GCF_020745735.1.1gdb
    GIXmake: Cannot open part file ./.20250203_HemOce_GCF_020745735.1.ktab.8 for writing


IN MEETING WITH RICHARD
We updated the version of FASTGA I have within the src directory within my base directory

Richard ran GIXmake command on GCF sequence and it worked
I ran it for the GCA sequence -- it failed due to same error as before
    (base) [ag2427@login-p-2 temp]$ GIXmake -v -P. 20250203_HemOce_GCA_020745765.1

    Creating genome index (GIX) 20250203_HemOce_GCA_020745765.1.gix in the current directory

    Partitioning K-mers via pos-lists into 8 parts
    Starting sort & index output of each part
        Outputing part 8
        Concat'ing parts
    GIXmake: Cannot open part file ./.20250203_HemOce_GCA_020745765.1.post.7 for writing

Wondered if error was due to syntax error in not specifying .1gdb -- corrected and re-ran
    (base) [ag2427@login-p-2 temp]$ GIXmake -v -P. 20250203_HemOce_GCA_020745765.1.1gdb 

    Creating genome index (GIX) 20250203_HemOce_GCA_020745765.1.gix in the current directory

    Partitioning K-mers via pos-lists into 8 parts
    Starting sort & index output of each part
        Outputing part 8
        Concat'ing parts
    GIXmake: Cannot open part file ./.20250203_HemOce_GCA_020745765.1.post.7 for writing

Checked FastGA source
    which FastGA
        ~/bin/FastGA
It's the same as it was before 
Need to change the path for the command


#### UPDATE ####
20250204 (February 4th, 2025)

Used an alias to change which GIXmake executable I am using
    alias GIXmake=~/src/FASTGA/GIXmake

    (base) [ag2427@login-q-1 FASTGA]$ which GIXmake
        alias GIXmake='/home/ag2427/src/FASTGA/GIXmake'
        ~/src/FASTGA/GIXmake
It appears to have worked -- re-run command
    GIXmake -v -P. 20250203_HemOce_GCA_020745765.1.1gdb
Job didn't work:
    Creating genome index (GIX) 20250203_HemOce_GCA_020745765.1.gix in the current directory

    Partitioning K-mers via pos-lists into 8 parts
    Starting sort & index output of each part
    Outputing part 8
    Concat'ing parts
    GIXmake: Cannot open part file ./.20250203_HemOce_GCA_020745765.1.post.6 for writing

Tried with a different syntax for alias command
    alias GIXmake="~/src/FASTGA/GIXmake"

    which GIXmake
        alias GIXmake='~/src/FASTGA/GIXmake'
            ~/src/FASTGA/GIXmake

(base) [ag2427@login-p-3 temp]$ GIXmake -v -P. 20250203_HemOce_GCA_020745765.1.1gdb

  Creating genome index (GIX) 20250203_HemOce_GCA_020745765.1.gix in the current directory

  Partitioning K-mers via pos-lists into 8 parts
  Starting sort & index output of each part
    Outputing part 8
    Concat'ing parts
GIXmake: Cannot open part file ./.20250203_HemOce_GCA_020745765.1.post.6 for writing


Modified 20250106_FROH_Calc.R, 20250106_Plot_ROH_FROH.sh, and 20250106_Plot_ROH.R scripts 
Changed 20250106_FROH_Calc.R to calculate inbreeding coefficient for just long and long/medium ROH and report all three values
Also calculating FROH for autosomal chromosomes in addition to all chromosomes
Added number of autosomal chromosomes variables in 20250106_Plot_ROH_FROH.sh -- use this to only calculate results for autosomal chromosomes
Changed 20250106_Plot_ROH.R to only plot ROH on autosomal chromosomes

Prepared 20250106_Plot_ROH_FROH.sh script for GCF_020745735.1 Hemiscyllium ocellatum (Epaulette shark)
    TODAY_DATE=20250204
    ROH_CALC_DATE=20250110
    CLADE=sharks
    REF_NAME=GCF_020745735.1
    NUM_AUT_CHROMOSOMES=52
awk command changed for chromosome names starting with NC
Submitted batch job 4703322


#### UPDATE ####
20250205 (February 5th, 2025)

Job 4703322 for plotting ROH/calculating FROH failed
    /var/spool/slurm/slurmd/job4703322/slurm_script: line 102: syntax error: unexpected end of file
I believe this error came because I had a syntax error and forgot 'fi' at the end of my first if statement in the file
Fixed and resubmitted
    TODAY_DATE=20250205
    ROH_CALC_DATE=20250110
    CLADE=sharks
    REF_NAME=GCF_020745735.1 
    NUM_AUT_CHROMOSOMES=52
Submitted batch job 4758172

Started working on 20250115_find_het_calc.R script to find a way to calculate heterozygosity without ROH included


#### UPDATE ####
20250206 (February 6th, 2025)

Finished updated 20250115_find_het_calc.R script to calculate heterozygosity excluding ROH 
Prepared 20250114_find_het_per_chr.sh script for GCF_020745735.1 Hemiscyllium ocellatum (Epaulette shark)
    TODAY_DATE=20250206
    ROH_DATE=20250116
    CLADE=sharks
    REF_NAME=GCF_020745735.1
    NUM_AUTOSOMAL_CHROMOSOMES=52 
    WINDOW_LENGTH=1000000
    WINDOW_INTERVAL=500000
Submitted batch job 4804450

Neither job 4804450 submitted today or job 4758172 submitted yesterday have run
Both were submitted on the cclake
Canceled both jobs and resubmitted on icelake-himem node

Submitted 20250114_find_het_per_chr.sh
Submitted batch job 4808121

Submitted 20250106_Plot_ROH_FROH.sh
Submitted batch job 4808131


#### UPDATE ####
20250207 (February 7th, 2025)

Job 4808121 and job 4808131 still haven't run
Canceled jobs and moved both to cclake-himem

Submitted 20250114_find_het_per_chr.sh 
Submitted batch job 4824384

Submitted 20250106_Plot_ROH_FROH.sh
Submitted batch job 4824382


#### UPDATE ####
20250208 (February 8th, 2025)

Job 4824382 completed for plotting ROH/calculating FROH of GCF_020745735.1 -- Hemiscyllium ocellatum (Epaulette shark)
Updated FROH calculation code worked!

Job 4824384 failed due to an error
Error:
    Error in calc_het(data_frame = dat, ROH_df = ROH_data, chromosome_length_file = chr_file,  : 
    object 'Het_excl_ROH' not found
    Execution halted
This error can from a syntax error when defining Het_excl_ROH
Fixed and resubmitted
Submitted batch job 4873787


#### UPDATE ####
20250210 (February 10th, 2025)

Job 4873787 completed for calculating heterozygosity of GCF_020745735.1 -- Hemiscyllium ocellatum (Epaulette shark)
It worked!

Prepared 20250114_find_het_per_chr.sh for GCA_035084275.1 Hydrolagus colliei (Spotted ratfish)
    TODAY_DATE=20250210
    ROH_DATE=20250116
    CLADE=sharks
    REF_NAME=GCA_035084275.1 
    NUM_AUTOSOMAL_CHROMOSOMES=40
    WINDOW_LENGTH=1000000
    WINDOW_INTERVAL=500000
Submitted batch job 4939060

Created 20250210_FROH_per_chr_calc.R 
This script will be eventually added to 20250106_Plot_ROH_FROH.sh
This script will be used to calculate FROH (the inbreeding coefficient) for each chromosome
Testing script with GCF_020745735.1 -- Hemiscyllium ocellatum (Epaulette shark)


Job 4939060 failed for finding heterozygosity of GCA_035084275.1 Hydrolagus colliei (Spotted ratfish)
    Error in file(file, "rt") : cannot open the connection
    Calls: read.table -> file
    In file(file, "rt") :
    cannot open file 'sharks/GCA_035084275.1/20250116_GCA_035084275.1_ROH_Results.tsv': No such file or directory
    Execution halted    
Error was due to incorrect date in ROH_DATE -- corrected and resubmitted
    ROH_DATE=20250117
Submitted batch job 4946048

Job 4946048 completed for finding heterozygosity of GCA_035084275.1 Hydrolagus colliei (Spotted ratfish)
It worked!

Prepared 20250114_find_het_per_chr.sh for GCA_035084215.1 Heptranchias perlo (sharpnose sevengill shark)
    TODAY_DATE=20250210
    ROH_DATE=20250117
    CLADE=sharks
    REF_NAME=GCA_035084215.1
    NUM_AUTOSOMAL_CHROMOSOMES=45
    WINDOW_LENGTH=1000000
    WINDOW_INTERVAL=500000 
Submitted batch job 4960501


#### UPDATE ####
20250211 (February 11th, 2025)

Job 4960501 completed for finding heterozygosity of GCA_035084215.1 Heptranchias perlo (sharpnose sevengill shark)
It worked!

Prepared 20250114_find_het_per_chr.sh for GCA_036365495.1 Heterodontus francisci (horn shark)
    TODAY_DATE=20250211
    ROH_DATE=20250118
    CLADE=sharks
    REF_NAME=GCA_036365495.1
    NUM_AUTOSOMAL_CHROMOSOMES=50
    WINDOW_LENGTH=1000000
    WINDOW_INTERVAL=500000
Submitted batch job 5002608

Finished writing up 20250210_FROH_per_chr_calc.R script

Prepared 20250101_mm2_alignment.sh for GCA_031021105.1.fa.gz Erythrolamprus reginae (royal ground snake) (alternates = GCA_031001685.1.fa.gz -- hap3, GCA_031021085.1.fa.gz -- hap2)
Will be doing analysis for this species twice, once with hap1/hap2 and another time with hap1/hap3 -- will be comparing results
Prepared for hap1/2 analysis
    REF_NAME=GCA_031021105.1 
    CLADE=reptiles
    REFERENCE=../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCA_031021105.1.fa.gz ## Finish filling in depending on genome and clade
    ALTERNATE=../241117.UCSC-hubs-VGP-alignment/alignment/alternate/reptiles/GCA_031021085.1.fa.gz ## Finish filling in depending on genome and clade
Output files have name ${REF_NAME}_hap12...
Changed node to cclake-himem
Submitted batch job 5008782

Created script 20250211_make_bam_plot_coverage.sh to create a bam file alignment for a species and then generate the coverage .txt file for it
Prepared 20250211_make_bam_plot_coverage.sh for GCF_020745735.1 -- Hemiscyllium ocellatum (Epaulette shark)
    TODAY_DATE=20250211
    CLADE=sharks
    SPECIES_NAME=HemOce
    REF_NAME=GCF_020745735.1
    ALT_NAME=GCA_020745765.1
Submitted batch job 5006700


#### UPDATE ####
20250212 (February 12th, 2025)

Job 5008782 completed for mm2 alignment of hap1/2 of GCA_031021105.1.fa.gz Erythrolamprus reginae (royal ground snake) (alternates = GCA_031001685.1.fa.gz -- hap3, GCA_031021085.1.fa.gz -- hap2)
It worked!

Prepared 20250101_mm2_alignment.sh script for hap1/3 of GCA_031021105.1.fa.gz Erythrolamprus reginae (royal ground snake) (alternates = GCA_031001685.1.fa.gz -- hap3, GCA_031021085.1.fa.gz -- hap2)
    REF_NAME=GCA_031021105.1 ## Change for reference genome of each species
    CLADE=reptiles
    REFERENCE=../241117.UCSC-hubs-VGP-alignment/alignment/reference/reptiles/GCA_031021105.1.fa.gz ## Finish filling in depending on genome and clade
    ALTERNATE=../241117.UCSC-hubs-VGP-alignment/alignment/alternate/reptiles/GCA_031001685.1.fa.gz ## Finish filling in depending on genome and clade
Output files have name ${REF_NAME}_hap13...
Submitted batch job 5036742

Prepared 20241231_ROH_Calc.sh for hap1/2 of GCA_031021105.1.fa.gz Erythrolamprus reginae (royal ground snake) (alternates = GCA_031001685.1.fa.gz -- hap3, GCA_031021085.1.fa.gz -- hap2)
    #SBATCH --array=1-14
    TODAY_DATE=20250212
    REF_NAME=GCA_031021105.1 
Updated vcf file names and output file names to have ..._hap12...
Submitted batch job 5036799

Job 5006700 completed for creating BAM plot and generating alignment coverage of GCF_020745735.1 -- Hemiscyllium ocellatum (Epaulette shark)
It worked!
Prepared 20250211_make_bam_plot_coverage.sh for GCA_035084275.1 Hydrolagus colliei (Spotted ratfish)
    TODAY_DATE=20250212
    CLADE=sharks
    SPECIES_NAME=HydCol
    REF_NAME=GCA_035084275.1 ## Change for reference genome of each species
    ALT_NAME=GCA_035084065.1
Submitted batch job 5037603

Job 5002608 completed for heterozygosity calculation of GCA_036365495.1 Heterodontus francisci (horn shark)
It worked!
Prepared 20250114_find_het_per_chr.sh for GCA_030028105.1 Mobula birostris (Giant manta)
    TODAY_DATE=20250212
    ROH_DATE=20250113
    CLADE=sharks
    REF_NAME=GCA_030028105.1
    NUM_AUTOSOMAL_CHROMOSOMES=32
    WINDOW_LENGTH=1000000 ## How long we want each window that we are using to measure heterozygosity
    WINDOW_INTERVAL=500000
Submitted batch job 5037750

Job 5037603 completed for creating BAM and calculating alignment coverage of GCA_035084275.1 Hydrolagus colliei (Spotted ratfish)
It worked!

Talked with Pío and he recommended using dotplots to check for places where there is not alignment between reference and altnerate in comparison to using the samtools coverage
Installed the pafr R package to create dotplots from the generated PAF files for each species
Testing with GCF_020745735.1 -- Hemiscyllium ocellatum (Epaulette shark)
It worked!
Generated sharks/HemOce/20250212_HemOce_DotPlot.pdf

Created 20250212_pafr_dotplot.sh to be able to submit scripts to create dotplots for each species through 20250212_pafr_dotplot.r
Prepared script for GCA_027789765.1.fa.gz Dendropsophus ebraccatus (hourglass treefrog) (alternate = GCA_027789725.1.fa.gz)
    TODAY_DATE=20250212
    REF_NAME=GCA_027789765.1
    CLADE=amphibians
Submitted batch job 5041419

Job 5041419 completed for making the dotplot of GCA_027789765.1.fa.gz Dendropsophus ebraccatus (hourglass treefrog) (alternate = GCA_027789725.1.fa.gz)
It worked!

Job 5036742 completed for mm2 alignment of hap1/3 of GCA_031021105.1.fa.gz Erythrolamprus reginae (royal ground snake) (alternates = GCA_031001685.1.fa.gz -- hap3, GCA_031021085.1.fa.gz -- hap2)
It worked!

Job 5036799 completed for calculating ROH of hap1/2 of GCA_031021105.1.fa.gz Erythrolamprus reginae (royal ground snake) (alternates = GCA_031001685.1.fa.gz -- hap3, GCA_031021085.1.fa.gz -- hap2)
It worked!

Prepared 20241231_ROH_Calc.sh for hap1/3 of GCA_031021105.1.fa.gz Erythrolamprus reginae (royal ground snake) (alternates = GCA_031001685.1.fa.gz -- hap3, GCA_031021085.1.fa.gz -- hap2)
    TODAY_DATE=20250212
    CLADE=reptiles
    REF_NAME=GCA_031021105.1
Updated vcf file names and output file names to have ..._hap13...
Submitted batch job 5049615


#### UPDATE ####
20250213 (February 13th, 2025)

Job 5037750 completed for calculating heterozygosity of GCA_030028105.1 Mobula birostris (Giant manta)
It worked!
Prepared 20250114_find_het_per_chr.sh for GCA_036971175.1 Narcine bancroftii (Spotted torpedo ray/Caribbean electric ray)
    TODAY_DATE=20250213
    ROH_DATE=20250118
    CLADE=sharks
    REF_NAME=GCA_036971175.1 
    NUM_AUTOSOMAL_CHROMOSOMES=14 
    WINDOW_LENGTH=1000000 
    WINDOW_INTERVAL=500000 
Submitted batch job 5070911

Job 5049615 completed for calculating ROH of hap1/3 of GCA_031021105.1.fa.gz Erythrolamprus reginae (royal ground snake) (alternates = GCA_031001685.1.fa.gz -- hap3, GCA_031021085.1.fa.gz -- hap2)
It worked!

Prepared 20250106_Plot_ROH_FROH.sh for hap1/2 of GCA_031021105.1.fa.gz Erythrolamprus reginae (royal ground snake) (alternates = GCA_031001685.1.fa.gz -- hap3, GCA_031021085.1.fa.gz -- hap2)
    TODAY_DATE=20250213
    ROH_CALC_DATE=20250112
    CLADE=sharks
    REF_NAME=GCA_031021105.1
    NUM_AUT_CHROMOSOMES=14
    VCF=${SPECIES}/${REF_NAME}_hap12_aligned.mm2.vcf
    ${ROH_CALC_DATE}_hap12_${line}_100kb_Results_ROH_Durbin_Calc.txt
    ${TODAY_DATE}_hap12_${REF_NAME}_ROH_Results.tsv 
    ${TODAY_DATE}_hap12_${line}.recode.vcf 
    ${TODAY_DATE}_hap12_autosomal_${line}
Submitted batch job 5071233


I want to run the ROH and heterozygosity analyses with human and mouse genomes to see if the results make SPECIES_NAME
Human --> GCF_009914755.1.fa.gz Homo sapiens (Alternate provided by VGP = GCA_016700455.2.fa.gz)
Mouse --> GCA_949316315.1.fa.gz Mus musculus (house mouse) (alternate = ) --> could use reference T2T generated by Sanger GCA_964188535.1

Prepared 20250101_mm2_alignment.sh for GCF_009914755.1.fa.gz Homo sapiens (Alternate provided by VGP = GCA_016700455.2.fa.gz)
    REF_NAME=GCF_009914755.1
    CLADE=mammals/primates
    REFERENCE=../241117.UCSC-hubs-VGP-alignment/alignment/reference/primates/GCA_031021105.1.fa.gz ## Finish filling in depending on genome and clade
    ALTERNATE=../241117.UCSC-hubs-VGP-alignment/alignment/alternate/primates/GCA_016700455.1.fa.gz ## Finish filling in depending on genome and clade
Submitted batch job 5075987


Job 5071233 failed for plotting ROH/calculating FROH of hap1/2 of GCA_031021105.1.fa.gz Erythrolamprus reginae (royal ground snake) (alternates = GCA_031001685.1.fa.gz -- hap3, GCA_031021085.1.fa.gz -- hap2)
Realized it was the clade was incorrectly put as 'sharks' -- corrected
    CLADE=reptiles
Resubmitted
Submitted batch job 5076449

Job 5076449 completed for plotting ROH/calculating FROH of hap1/2 of GCA_031021105.1.fa.gz Erythrolamprus reginae (royal ground snake) (alternates = GCA_031001685.1.fa.gz -- hap3, GCA_031021085.1.fa.gz -- hap2)
It finished but job did not work due to typo in the ROH_CALC_DATE, incorrectly stating it as 20250112, corrected error
    ROH_CALC_DATE=20250212
Resubmitted
Submitted batch job 5079897

Job 5070911 completed for calculating heterozygosity of GCA_036971175.1 Narcine bancroftii (Spotted torpedo ray/Caribbean electric ray)
It worked!


#### UPDATE ####
20250217 (February 17th, 2025)

Job 5079897 completed for plotting ROH/calculating FROH of hap1/2 of GCA_031021105.1.fa.gz Erythrolamprus reginae (royal ground snake) (alternates = GCA_031001685.1.fa.gz -- hap3, GCA_031021085.1.fa.gz -- hap2)
Error:
    Error in read.table(args[6], header = FALSE) : 
    no lines available in input
When checking output files, the problem seems to be with getting the Chroms_Lengths.txt and ROH_Results.tsv files, which have file sizes of 11 and 0 bytes respectively
Problem for generating Chroms_Lengths file is due to command looking for chromosomes starting with NC instead of CM -- Fixed
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/${CLADE}/${REF_NAME}.fa.gz | awk '$0 ~ ">CM" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' >${SPECIES}/Reference_${REF_NAME}_Chroms_Lengths.txt
Removed empty ROH_results and Chroms_Lengths files 
    rm -f -r 20250213_hap12_GCA_031021105.1_ROH_Results.tsv
    rm -f -r Reference_GCA_031021105.1_Chroms_Lengths.txt
Removed empty output files
    rm -f -r 20250213_hap12_CM*
    rm -f -r 20250213_hap12_GCA_031021105.1_ROH_Map.pdf
    rm -f -r 20250213_hap12_GCA_031021105.1_FROH.txt
Updated TODAY_DATE 
    TODAY_DATE=20250217
Resubmitted
Submitted batch job 5197519

Job 5075987 failed for mm2 alignment of GCF_009914755.1.fa.gz Homo sapiens (Alternate provided by VGP = GCA_016700455.2.fa.gz)
Error:
    failed to open file '../241117.UCSC-hubs-VGP-alignment/alignment/reference/primates/GCA_031021105.1.fa.gz': No such file or directory
Error was due to me failing to update the reference name in the reference path -- fixed
    REFERENCE=../241117.UCSC-hubs-VGP-alignment/alignment/reference/primates/GCF_009914755.1.fa.gz ## Finish filling in depending on genome and clade
I confirmed the altnerate was correct
Resubmitted
Submitted batch job 5197735

Job 5197519 completed but failed for plotting ROH/calculating FROH of hap1/2 of GCA_031021105.1.fa.gz Erythrolamprus reginae (royal ground snake) (alternates = GCA_031001685.1.fa.gz -- hap3, GCA_031021085.1.fa.gz -- hap2)
Error:
    Error in L_ROH/chrom_length : non-numeric argument to binary operator
    Execution halted
Checked 20250210_FROH_per_chr_calc.sh script -- appears to be working fine
Removed all temporary files in the reptiles/GCA_031021105.1/temp directory in case they were corrupted and somehow affecting the results
Resubmitted
Submitted batch job 5203661

Working on HemOce GIXmake command
Removed all aliases from bashrc, removed old FASTGA commands from bin, and copied FASTGA commands from src/FASTGA to bin using cp command
     which GIXmake
        ~/bin/GIXmake
Tested with more than 2 threads to see if error is fixed
Tested in sharks/HemOce/temp directory where the .1gdb file is located
    GIXmake -v -P. -T6 20250203_HemOce_GCA_020745765.1.1gdb
Failed with error:
    GIXmake: Cannot open part file ./.20250203_HemOce_GCA_020745765.1.post.5 for writing
Tried:
    GIXmake -v -P. -T4 20250203_HemOce_GCA_020745765.1.1gdb
Failed with error:
    GIXmake: Cannot open part file ./.20250203_HemOce_GCA_020745765.1.post.3 for writing
After adding the error commands I re-ran it 
    GIXmake -v -P. -T8 20250203_HemOce_GCA_020745765.1.1gdb
Failed with error:
    GIXmake: Cannot open part file ./.20250203_HemOce_GCA_020745765.1.post.7 for writing. Error: Too many open files


Loop in GIXmake.c before I made any modifications
  for (p = 0; p < NTHREADS; p++)
    { carm[p].tout = open(Catenate(TPATH,"/.",TROOT,
                          Numbered_Suffix(".ktab.",p+1,"")),
                          O_WRONLY|O_CREAT|O_TRUNC,S_IRWXU);
      if (carm[p].tout < 0)
        { fprintf(stderr,"%s: Cannot open part file %s/.%s.ktab.%d for writing\n",
                         Prog_Name,TPATH,TROOT,p+1);
          goto remove_parts;
        }
      carm[p].pout = open(Catenate(TPATH,"/.",TROOT,
                          Numbered_Suffix(".post.",p+1,"")),
                          O_WRONLY|O_CREAT|O_TRUNC,S_IRWXU);
      if (carm[p].pout < 0)
        { fprintf(stderr,"%s: Cannot open part file %s/.%s.post.%d for writing\n",
                         Prog_Name,TPATH,TROOT,p+1);
          goto remove_parts;
        }
    }

FOUND SOLUTION: HAVE TO RUN ALL FASTGA COMMANDS IN TERMINAL
FastGA -v -P. -T8 20250203_HemOce_GCF_020745735.1.gix 20250203_HemOce_GCA_020745765.1.gix -1:20250203_HemOce_FASTGA
ALNtoPAF -x -T8 20250203_HemOce_FASTGA.1aln > 20250217_HemOce_FASTGA.paf

Job 5197735 failed for mm2 alignment of GCF_009914755.1.fa.gz Homo sapiens (Alternate provided by VGP = GCA_016700455.2.fa.gz)
Error:
    ERROR: failed to open file '../241117.UCSC-hubs-VGP-alignment/alignment/alternate/primates/GCA_016700455.1.fa.gz': No such file or directory
Error is that accession number ends in .2, NOT .1 -- fixed and resubmitted
Submitted batch job 5214436

Prepared 20250211_make_bam_plot_coverage.sh for GCA_035084215.1 Heptranchias perlo (sharpnose sevengill shark)
    TODAY_DATE=20250217
    CLADE=sharks
    SPECIES_NAME=HepPer
    REF_NAME=GCA_035084215.1 ## Change for reference genome of each species
    ALT_NAME=GCA_035084135.1
Submitted batch job 5215045

Job 5203661 for plotting ROH/calculating FROH of hap1/2 of GCA_031021105.1.fa.gz Erythrolamprus reginae (royal ground snake) (alternates = GCA_031001685.1.fa.gz -- hap3, GCA_031021085.1.fa.gz -- hap2)
Same error as before:
    Error in L_ROH/chrom_length : non-numeric argument to binary operator
Altered chrom_length variable in 20250210_FROH_per_chr_calc.sh 
    chrom_length <- as.numeric(args[12])
Resubmitted
Submitted batch job 5215703

Job 5214436 completed for mm2 alignment of GCF_009914755.1.fa.gz Homo sapiens (Alternate provided by VGP = GCA_016700455.2.fa.gz)
It worked! 

Prepared 20241231_ROH_Calc.sh for GCF_009914755.1.fa.gz Homo sapiens (Alternate provided by VGP = GCA_016700455.2.fa.gz)
    #SBATCH --array=1-24
    TODAY_DATE=20250217
    CLADE=mammals/primates
    REF_NAME=GCF_009914755.1
Submitted batch job 5219322

Job 5215703 finished incorrectly for plotting ROH/calculating FROH of hap1/2 of GCA_031021105.1.fa.gz Erythrolamprus reginae (royal ground snake) (alternates = GCA_031001685.1.fa.gz -- hap3, GCA_031021085.1.fa.gz -- hap2)
Error:
    In file(file, "rt") :
    cannot open file 'reptiles/GCA_031021105.1/GCA_031021105.1_aligned.mm2.vcf': No such file or directory
This error was in 20250106_FROH_Calc.R -- didn't account for hap12 inclusion in file name -- fixed
    vcf_name <- paste(clade, "/", Ref_name, "/", Ref_name, '_hap12_aligned.mm2.vcf', sep='')
Re-ran using bash command
It worked!

Prepared 20250106_Plot_ROH_FROH.sh for hap1/3 of GCA_031021105.1.fa.gz Erythrolamprus reginae (royal ground snake) (alternates = GCA_031001685.1.fa.gz -- hap3, GCA_031021085.1.fa.gz -- hap2)
Changed all ...hap12... to ...hap13....
Submitted batch job 5225013

Job 5215045 completed for creating BAM and plotting coverage of for GCA_035084215.1 Heptranchias perlo (sharpnose sevengill shark)
It worked!

Preparing 20250211_make_bam_plot_coverage.sh for GCA_036365495.1 Heterodontus francisci (horn shark)
    TODAY_DATE=20250217
    CLADE=sharks
    SPECIES_NAME=HetFra
    REF_NAME=GCA_036365495.1 ## Change for reference genome of each species
    ALT_NAME=GCA_036365525.1
Submitted batch job 5226270

Job 5225013 completed for plotting ROH/calculating FROH of hap1/3 of GCA_031021105.1.fa.gz Erythrolamprus reginae (royal ground snake) (alternates = GCA_031001685.1.fa.gz -- hap3, GCA_031021085.1.fa.gz -- hap2)
It worked!

Prepared 20250106_Plot_ROH_FROHs.sh for GCA_035084275.1 Hydrolagus colliei (Spotted ratfish)
    TODAY_DATE=20250217
    ROH_CALC_DATE=20250117
    CLADE=sharks
    REF_NAME=GCA_035084275.1 
    NUM_AUT_CHROMOSOMES=40
zcat line in line 43 prepped for chromosomes starting with CM
Submitted batch job 5226772


#### UPDATE ####
20250218 (February 18th, 2025)

Job 5226772 failed for plotting ROH/calculating FROH of GCA_035084275.1 Hydrolagus colliei (Spotted ratfish)
Could not read ROH text files -- incorrect ROH_CALC_DATE
Fixed and resubmitted
    ROH_CALC_DATE=20250116
Submitted batch job 5244288

Job 5226270 failed for creating BAM/calculating coverage of GCA_036365495.1 Heterodontus francisci (horn shark)
Failed due to an out of memory error
Updated memory to 100GB and updated the date to 20250218
Removed output txt file
     rm -f -r 20250217_HetFra.txt
Resubmitted
Submitted batch job 5244388

Job 5244288 completed for plotting ROH/calculating FROH of GCA_035084275.1 Hydrolagus colliei (Spotted ratfish)
It worked!

Prepared 20250106_Plot_ROH_FROH.sh for GCA_035084215.1 Heptranchias perlo (sharpnose sevengill shark)
    TODAY_DATE=20250218
    ROH_CALC_DATE=20250117
    CLADE=sharks
    REF_NAME=GCA_035084215.1 
    NUM_AUT_CHROMOSOMES=45
zcat/awk on line 43 prepared for chromosomes starting with CM...
Submitted batch job 5249193

Job 5249193 completed for plotting ROH/calculating FROH of GCA_035084215.1 Heptranchias perlo (sharpnose sevengill shark)
It worked!

Prepared 20250106_Plot_ROH_FROH.sh for GCA_036365495.1 Heterodontus francisci (horn shark)
    TODAY_DATE=20250218
    ROH_CALC_DATE=20250117
    CLADE=sharks
    REF_NAME=GCA_036365495.1 
    NUM_AUT_CHROMOSOMES=50
zcat/awk command on line 43 prepared for chromosomes starting with CM...
Submitted batch job 5249640

Job 5249640 completed for plotting ROH/calculating FROH of GCA_036365495.1 Heterodontus francisci (horn shark)
It worked!

Prepared 20250106_Plot_ROH_FROH.sh for GCA_030028105.1 Mobula birostris (Giant manta)
    TODAY_DATE=20250218
    ROH_CALC_DATE=20250113
    CLADE=sharks
    REF_NAME=GCA_030028105.1 
    NUM_AUT_CHROMOSOMES=32
zcat/awk command on line 43 prepared for chromosomes starting with CM...
Submitted batch job 5250378

Job 5250378 completed for plotting ROH/calculating FROH of GCA_030028105.1 Mobula birostris (Giant manta)
It worked!

Prepared 20250106_Plot_ROH_FROH.sh for GCA_036971175.1 Narcine bancroftii (Spotted torpedo ray/Caribbean electric ray)
    TODAY_DATE=20250218
    ROH_CALC_DATE=20250118
    CLADE=sharks
    REF_NAME=GCA_036971175.1 
    NUM_AUT_CHROMOSOMES=14
zcat/awk command on line 43 prepared for chromosomes starting with CM...
Submitted batch job 5251577

Job 5251577 completed for plotting ROH/calculating FROH of GCA_036971175.1 Narcine bancroftii (Spotted torpedo ray/Caribbean electric ray)
It worked!

Prepped 20250106_Plot_ROH_FROH.sh for GCF_030144855.1 Hypanus sabinus (Atlantic stingray)
    TODAY_DATE=20250218
    ROH_CALC_DATE=20250117
    CLADE=sharks
    REF_NAME=GCF_030144855.1 
    NUM_AUT_CHROMOSOMES=32
zcat/awk command on line 43 prepped for chromosomes starting with NC...
Submitted batch job 5254078

Job 5244388 failed for creating BAM/calculating coverage of GCA_036365495.1 Heterodontus francisci (horn shark)
Error:
    [E::parse_cigar] CIGAR length too long at position 1 (268907262H)
    [W::sam_read1_sam] Parse error at line 1513
    samtools sort: truncated file. Aborting
    [E::hts_open_format] Failed to open file "sharks/GCA_036365495.1/20250218_HetFra.bam" : No such file or directory
    samtools coverage: Could not open "sharks/GCA_036365495.1/20250218_HetFra.bam": No such file or directory

Job 5254078 completed for plotting ROH/calculating FROH of GCF_030144855.1 Hypanus sabinus (Atlantic stingray)
It worked!


#### UPDATE ####
20250219 (February 19th, 2025)

Modified zcat/awk command in line 43 of 20250106_Plot_ROH_FROH.sh to make sure that the last chromosome is the correct length and doesn't have unplaced scaffolds added onto it
    zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/GCA_035084275.1.fa.gz | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' | grep 'CM' > sharks/GCA_035084275.1/test.txt

Prepared 20250106_Plot_ROH_FROH.sh for GCA_035084275.1 Hydrolagus colliei (Spotted ratfish)
    TODAY_DATE=20250219
    ROH_CALC_DATE=20250116
    CLADE=sharks
    REF_NAME=GCA_035084275.1
    NUM_AUT_CHROMOSOMES=40
    CHROM_START_CHR=CM
Ran with bash to make sure it worked and I could troubleshoot as needed
It works!

Created 20250219_ROH_Calc_V2.sh and 20250219_ROH_Durbin_Calc_Eqns_V2.R
These will calculate ROH when only looking at areas of genome that are successfully aligned


#### UPDATE ####
20250220 (February 20th, 2025)

Finished writing 20250219_ROH_Calc_V2.sh and 20250219_ROH_Durbin_Calc_Eqns_V2.R
Prepared 20250219_ROH_Calc_V2.sh for GCA_035084275.1 Hydrolagus colliei (Spotted ratfish)
    #SBATCH --array=1-40
    TODAY_DATE=20250220
    CLADE=sharks
    REF_NAME=GCA_035084275.1
Submitted batch job 5373292

Canceled job 5373292 after it ran through 2 chromosomes and saw that it was not working
    tat error: No such file or directory
    Error: Can't determine file type of sharks/GCA_035084275.1/GCA_035084275.1_aligned.mm2.vc

    NAs introduced by coercion 
    Error in file(con, "r") : invalid 'description' argument
    Calls: read_paf -> readLines -> file
    Execution halted
Due to syntax error in ilne 40 of 20250219_ROH_Calc_V2.sh file -- fixed
    VCF=${SPECIES}/${REF_NAME}_aligned.mm2.vcf
Removed the empty generated txt files from within sharks/GCA_035084275.1
    rm -f -r 20250220_*
Resubmitted
Submitted batch job 5373930

Started FASTGA alignment for GCA_035084275.1 Hydrolagus colliei (Spotted ratfish)
Will attempt to do it through a shell script (submitted through terminal) instead of manually
Prepared 20250129_FASTGA_alignment.sh for GCA_035084275.1 Hydrolagus colliei (Spotted ratfish)
    TODAY_DATE=20250220
    REF_NAME=GCA_035084275.1
    ALT_NAME=GCA_035084065.1
    CLADE=sharks
    SPEC_NAME=HydCol
Submitted batch job 5374983

Canceled job 5373930 after it ran through 6 chromosomes and I saw that it was not working
Error:
    NAs introduced by coercion 
    Error in file(con, "r") : invalid 'description' argument
    Calls: read_paf -> readLines -> file
    Execution halted
Realized issue (A rather silly one at that) was that I had forgotten to define a PAF variable in the shell script
Fixed
    PAF=${SPECIES}/${REF_NAME}_aligned.mm2.srt.paf
Resubmitted
Submitted batch job 5375191

Job 5374983 failed for FASTGA alignment of GCA_035084275.1 Hydrolagus colliei (Spotted ratfish)
Error:
    GIXmake: Cannot create a .gix with a different location and root name than its .gdb
    GIXmake: Cannot create a .gix with a different location and root name than its .gdb
GDB files were successfully made
Altering first GIXmake command to not specify an output file name to see if that resolves the error
    GIXmake -v -P. -T8 ${TEMP}/${TODAY_DATE}_${SPEC_NAME}_${REF_NAME}.1gdb 
Resubmitted
Submitted batch job 5375630

Job 5375630 completed
Changing the GIXmake command worked!
Changed the other GIXmake file
    GIXmake -v -P. -T8 ${TEMP}/${TODAY_DATE}_${SPEC_NAME}_${ALT_NAME}.1gdb
Resubmitted
Submitted batch job 5375952

Canceled job 5375191 due to it not working on the first several chromosomes
Error:
    Warning message:
    NAs introduced by coercion 
    Warning messages:
    1: In .make_numeric(as.data.frame(res), c(2, 3, 4, 7:12)) :
    NAs introduced by coercion
    2: In .make_numeric(as.data.frame(res), c(2, 3, 4, 7:12)) :
    NAs introduced by coercion
    3: In .make_numeric(as.data.frame(res), c(2, 3, 4, 7:12)) :
    NAs introduced by coercion
    4: In .make_numeric(as.data.frame(res), c(2, 3, 4, 7:12)) :
    NAs introduced by coercion
    5: In .make_numeric(as.data.frame(res), c(2, 3, 4, 7:12)) :
    NAs introduced by coercion
    6: In .make_numeric(as.data.frame(res), c(2, 3, 4, 7:12)) :
    NAs introduced by coercion
    7: In .make_numeric(as.data.frame(res), c(2, 3, 4, 7:12)) :
    NAs introduced by coercion
    Error in file(file, "rt") : invalid 'description' argument
    Calls: read.table -> file
    Execution halted


#### UPDATE ####
20250221 (February 21st, 2025)

Job 5375952 completed for FASTGA alignment of  GCA_035084275.1 Hydrolagus colliei (Spotted ratfish)
The GIX files were successfully made, but FASTGA did not run
To confirm if this is a syntax error or error with running FastGA in a shell script, I will run FastGA command in terminal
Ran in sharks/HydCol
    FastGA -v -P. -1:20250220_HydCol_ALN temp/20250220_HydCol_GCA_035084275.1.gix temp/20250220_HydCol_GCA_035084065.1.gix
This worked!
Upon double checking syntax, I updated syntax of FastGA command in shell script
    FastGA -v -P. -T8 -1:${SPECIES}/${TODAY_DATE}_${SPEC_NAME}_ALN ${TEMP}/${TODAY_DATE}_${SPEC_NAME}_${REF_NAME}.gix ${TEMP}/${TODAY_DATE}_${SPEC_NAME}_${ALT_NAME}.gix
Manually ran ALNtoPAF in sharks/HydCol through terminal
    ALNtoPAF -x 20250220_HydCol_ALN.1aln
Canceled running this as it was outputting results into command line -- updated command to specify output file
    ALNtoPAF -x 20250220_HydCol_ALN.1aln > 20250220_HydCol_FASTGA.paf
It appears to have worked successfully!
Updated ALNtoPAF script in 20250129_FASTGA_alignment.sh
    ALNtoPAF -x ${SPECIES}/${TODAY_DATE}_${SPEC_NAME}_ALN.1aln > ${TODAY_DATE}_${SPEC_NAME}_FASTGA.paf

Preparing 20250129_FASTGA_alignment for GCA_035084215.1 Heptranchias perlo (sharpnose sevengill shark)
    TODAY_DATE=20250221
    REF_NAME=GCF_035084215.1
    ALT_NAME=GCA_035084135.1
    CLADE=sharks
    SPEC_NAME=HepPer
Submitted batch job 5427575

Alterend some syntax within 20250219_ROH_Durbin_Calc_Eqns_V2.R to see if that fixed errors
Running 20250219_ROH_Calc_V2.sh for GCA_035084275.1 Hydrolagus colliei (Spotted ratfish)
    #SBATCH --array=1-40
    TODAY_DATE=20250221
    CLADE=sharks
    REF_NAME=GCA_035084275.1 
Submitted batch job 5428108

Job 5427575 finished for FASTGA alingment of GCA_035084215.1 Heptranchias perlo (sharpnose sevengill shark)
Job did NOT complete
Generated .1gdb and .gix files for the alternate genome, but not the reference
Realized syntax error was the REF_NAME being erronenously labeled with GCF instead of GCA
Fixed
    REF_NAME=GCA_035084215.1
Removed empty .1aln file from sharks/HepPer directory
    rm -f -r 20250221_HepPer_ALN.1aln
Resubmitted
Submitted batch job 5428307

Canceled job 5428108 due to same error coming up
Error:
Warning message:
    NAs introduced by coercion 
    Warning messages:
    1: In .make_numeric(as.data.frame(res), c(2, 3, 4, 7:12)) :
    NAs introduced by coercion
    2: In .make_numeric(as.data.frame(res), c(2, 3, 4, 7:12)) :
    NAs introduced by coercion
    3: In .make_numeric(as.data.frame(res), c(2, 3, 4, 7:12)) :
    NAs introduced by coercion
    4: In .make_numeric(as.data.frame(res), c(2, 3, 4, 7:12)) :
    NAs introduced by coercion
    5: In .make_numeric(as.data.frame(res), c(2, 3, 4, 7:12)) :
    NAs introduced by coercion
    6: In .make_numeric(as.data.frame(res), c(2, 3, 4, 7:12)) :
    NAs introduced by coercion
    7: In .make_numeric(as.data.frame(res), c(2, 3, 4, 7:12)) :
    NAs introduced by coercion
    Error in file(file, "rt") : invalid 'description' argument
    Calls: read.table -> file
    Execution halted
I found an error (possible causing this?) -- Chrom_length variable in the shell script never defined
Fixed
    CHROM_START_CHR=CM

    ## Generating file with names of all chromosomes and their lengths (if not already made)
    if [ -f ${SPECIES}/Reference_${REF_NAME}_Chroms_Lengths.txt ]
        then
            echo 'Chrom_length file exists'
        else
            zcat < ../241117.UCSC-hubs-VGP-alignment/alignment/reference/${CLADE}/${REF_NAME}.fa.gz | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' | grep ${CHROM_START_CHR} > ${SPECIES}/Reference_${REF_NAME}_Chroms_Lengths.txt
    fi

    CHROM_LENGTH="awk -v var=${CHROM} '{ if ($1 == var) { print $2 }}' ${SPECIES}/Reference_${REF_NAME}_Chroms_Lengths.txt"
Resubmitted with fixed script
Submitted batch job 5432508


Job 5428307 finished for FASTGA alingment of GCA_035084215.1 Heptranchias perlo (sharpnose sevengill shark)
It worked!

Preparing 20250129_FASTGA_alignment.sh for GCA_036365495.1 Heterodontus francisci (horn shark)
    TODAY_DATE=20250221
    REF_NAME=GCA_036365495.1
    ALT_NAME=GCA_036365525.1
    CLADE=sharks
    SPEC_NAME=HetFra
Submitted batch job 5432048


#### UPDATE ####
20250224 (February 24th, 2025)

Job 5432048 for FASTGA alingment of GCA_036365495.1 Heterodontus francisci (horn shark) completed
It worked!

Job 5432508 failed for ROH calculation of GCA_035084275.1 Hydrolagus colliei (Spotted ratfish)
ROH_Results.txt files were generated but with file size of 0
Error:
    Warning message:
    NAs introduced by coercion 
    Warning message:
    NAs introduced by coercion 
    Error in file(con, "r") : cannot open the connection
    Calls: read_paf -> readLines -> file
    In addition: Warning message:
    In file(con, "r") : cannot open file 'if': No such file or directory
    Execution halted
Removed the empty files
    rm -f -r 20250221_*
Changed by having quotation marks around ${PAF} in the shell script input for the R script like I had for the vcf file name
Updated TODAY_DATE
    TODAY_DATE=20250224
Resubmitted
Submitted batch job 5573160

Prepared 20250129_FASTGA_alignment.sh for GCA_030028105.1 Mobula birostris (Giant manta)
    TODAY_DATE=20250224
    REF_NAME=GCA_030028105.1
    ALT_NAME=GCA_030035685.1
    CLADE=sharks
    SPEC_NAME=MobBir
Submitted batch job 5573452

Canceled job 5573160 due to same error as before
    Warning message:
    NAs introduced by coercion 
    Warning message:
    NAs introduced by coercion 
    Error in file(con, "r") : cannot open the connection
    Calls: read_paf -> readLines -> file
    In addition: Warning message:
    In file(con, "r") : cannot open file 'if': No such file or director
Problem is with reading the PAF file --> args[11] is the string 'if' instead of paf file name
Printed args
Actual problem is with defining the CHROM_LENGTH variable
     [1] "/usr/lib64/R/bin/exec/R"                                            
    [2] "--no-echo"                                                          
    [3] "--no-restore"                                                       
    [4] "--file=20250219_ROH_Durbin_Calc_Eqns_V2.R"                          
    [5] "--args"                                                             
    [6] "CM068742.1"                                                         
    [7] "awk"                                                                
    [8] "-v"                                                                 
    [9] "var=CM068742.1"                                                     
    [10] "'{"                                                                 
    [11] "if"                                                                 
    [12] "("                                                                  
    [13] "=="                                                                 
    [14] "var)"                                                               
    [15] "{"                                                                  
    [16] "print"                                                              
    [17] "}}'"                                                                
    [18] "sharks/GCA_035084275.1/Reference_GCA_035084275.1_Chroms_Lengths.txt"
    [19] "20250224"                                                           
    [20] "GCA_035084275.1"                                                    
    [21] "sharks"                                                             
    [22] "sharks/GCA_035084275.1/GCA_035084275.1_aligned.mm2.srt.paf"         
    [23] "sharks/GCA_035084275.1/temp/20250224_CM068742.1.recode.vcf"      
Fixed by changing syntax of variable
    CHROM_LENGTH=$(awk -v var=${CHROM} '{ if ($1 == var) { print $2 }}' ${SPECIES}/Reference_${REF_NAME}_Chroms_Lengths.txt)
     [1] "/usr/lib64/R/bin/exec/R"                                   
    [2] "--no-echo"                                                 
    [3] "--no-restore"                                              
    [4] "--file=20250219_ROH_Durbin_Calc_Eqns_V2.R"                 
    [5] "--args"                                                    
    [6] "CM068742.1"                                                
    [7] "156096171"                                                 
    [8] "20250224"  
    [9] "GCA_035084275.1"                                           
    [10] "sharks"                                                    
    [11] "sharks/GCA_035084275.1/GCA_035084275.1_aligned.mm2.srt.paf"
    [12] "sharks/GCA_035084275.1/temp/20250224_CM068742.1.recode.vcf"
Resubmitted
Submitted batch job 5574418

Canceled Job 5574418 due to error 
Error:
    Error in file(con, "r") : cannot open the connection
    Calls: read_paf -> readLines -> file
    In addition: Warning message:
    In file(con, "r") :
    cannot open file 'sharks/GCA_035084275.1/GCA_035084275.1_aligned.mm2.srt.paf': No such file or directory
    Execution halted
Fixed name in shell script
    PAF=${SPECIES}/${REF_NAME}.mm2.srt.paf
Removed empty files generated
    rm -f -r 20250224_*
Resubmitted
Submitted batch job 5574763

Job 5573452 completed for FASTGA alingment of GCA_030028105.1 Mobula birostris (Giant manta)
It worked!

Prepared 20250129_FASTGA_alignment.sh for GCA_036971175.1 Narcine bancroftii (Spotted torpedo ray/Caribbean electric ray)
    TODAY_DATE=20250224
    REF_NAME=GCA_036971175.1
    ALT_NAME=GCA_036971445.1
    CLADE=sharks
    SPEC_NAME=NarBan
Submitted batch job 5576812


#### UPDATE ####
20250225 (February 25th, 2025)

Job 5574763 ended for ROH calculation of GCA_035084275.1 Hydrolagus colliei (Spotted ratfish)
Job did not finish due to error:
    Error:
    Loading required package: ggplot2
    /var/spool/slurm/slurmd/job5578403/slurm_script: line 70: 708519 Killed                  Rscript 20250219_ROH_Durbin_Calc_Eqns_V2.R ${CHROM} ${CHROM_LENGTH} ${TODAY_DATE} ${REF_NAME} ${CLADE} "${PAF}" "${TEMP}/${TODAY_DATE}_${CHROM}.recode.vcf" > ${SPECIES}/${TO
    DAY_DATE}_${CHROM}_ROH_Results.txt
    slurmstepd: error: Detected 1 oom_kill event in StepId=5578403.batch. Some of the step tasks have been OOM Killed.
Updated TODAY_DATE
    TODAY_DATE=20250225
Updated memory allotment
    #SBATCH --mem=120GB
Resubmitted
Submitted batch job 5631290

Job 5576812 completed for FASTGA alingment of GCA_036971175.1 Narcine bancroftii (Spotted torpedo ray/Caribbean electric ray)
It worked!
Prepared 20250129_FASTGA_alignment.sh for GCF_030144855.1 Hypanus sabinus (Atlantic stingray)
    TODAY_DATE=20250225
    REF_NAME=GCF_030144855.1
    ALT_NAME=GCA_030144785.1
    CLADE=sharks
    SPEC_NAME=HypSab
Submitted batch job 5631592


#### UPDATE ####
20250226 (February 26th, 2025)

Job 5631592 completed for FASTGA alignment of GCF_030144855.1 Hypanus sabinus (Atlantic stingray)
It worked!
Prepared 20250129_FASTGA_alignment.sh for GCA_027789765.1.fa.gz Dendropsophus ebraccatus (hourglass treefrog) (alternate = GCA_027789725.1.fa.gz)
    TODAY_DATE=20250226
    REF_NAME=GCA_027789765.1
    ALT_NAME=GCA_027789725.1
    CLADE=amphibians
    SPEC_NAME=DenEbr
Submitted batch job 5699517

Job 5631290 finished for ROH calculation of GCA_035084275.1 Hydrolagus colliei (Spotted ratfish)
It did not work -- only one ROH results txt file was generated and it had a file size of zero
All chromosome runs had this error:
Error:
    Error in unlist(apply(cg_df, MARGIN = 1, count_cg_type)) : 
    long vectors not supported yet: ../../src/include/Rinlinedfuns.h:537
    Calls: create_read_df -> unlist
Updated the count_cg_type function to remove the for loop to avoid this error
Removed the empty ROH_results.txt file that was made
    rm -f -r 20250225_CM068742.1_ROH_Results.txt
Resubmitted
Submitted batch job 5703727


#### UPDATE ####
20250227

Job 5699517 completed for FASTGA alignment of GCA_027789765.1.fa.gz Dendropsophus ebraccatus (hourglass treefrog) (alternate = GCA_027789725.1.fa.gz)
It worked!
Prepared 20250129_FASTGA_alignment.sh for GCA_031893055.1.fa.gz Leptodactylus fuscus (rufous frog) (alternate = GCA_031893025.1.fa.gz)
    TODAY_DATE=20250227
    REF_NAME=GCA_031893055.1
    ALT_NAME=GCA_031893025.1
    CLADE=amphibians
    SPEC_NAME=LepFus
Submitted batch job 5749589

Job 5703727 failed for ROH calculation of GCA_035084275.1 Hydrolagus colliei (Spotted ratfish)
Error:
    Error in (function (..., row.names = NULL, check.rows = FALSE, check.names = TRUE,  : 
    arguments imply differing number of rows: 1524855, 1532852, 1532872, 1537339
    Calls: create_read_df ... as.data.frame -> as.data.frame.list -> do.call -> <Anonymous>
    Execution halted

    slurmstepd: error: Detected 1 oom_kill event in StepId=5737668.batch. Some of the step tasks have been OOM Killed.
In the error for differing number of rows: these numbers chnaged depending on the chromosome

Prepared 20250212_pafr_dotplot.sh for GCA_035084275.1 Hydrolagus colliei (Spotted ratfish)
    TODAY_DATE=20250227
    REF_NAME=GCA_035084275.1
    CLADE=sharks
Submitted batch job 5752247

Prepared 20250211_make_bam_plot_coverage.sh for GCA_030028105.1 Mobula birostris (Giant manta)
    TODAY_DATE=20250227
    CLADE=sharks
    SPECIES_NAME=MobBir
    REF_NAME=GCA_030028105.1 ## Change for reference genome of each species
    ALT_NAME=GCA_030035685.1
Submitted batch job 5752590

Job 5752247 completed for creating dotplot of GCA_035084275.1 Hydrolagus colliei (Spotted ratfish)
It worked!

Prepared 20250212_pafr_dotplot.sh for GCA_035084215.1 Heptranchias perlo (sharpnose sevengill shark)
    TODAY_DATE=20250227
    REF_NAME=GCA_035084215.1
    CLADE=sharks
Submitted batch job 5753301

Job 5752590 completed for BAM plot coverage of GCA_030028105.1 Mobula birostris (Giant manta)
It worked!
Prepared 20250211_make_bam_plot_coverage.sh for GCA_036971175.1 Narcine bancroftii (Spotted torpedo ray/Caribbean electric ray)
    TODAY_DATE=20250227
    CLADE=sharks
    SPECIES_NAME=NarBan
    REF_NAME=GCA_036971175.1 ## Change for reference genome of each species
    ALT_NAME=GCA_036971445.1
Submitted batch job 5765892

Job 5753301 completed for paf dotplot of GCA_035084215.1 Heptranchias perlo (sharpnose sevengill shark)
It worked!
Prepared 20250212_pafr_dotplot.sh for GCA_036365495.1 Heterodontus francisci (horn shark)
    TODAY_DATE=20250227
    REF_NAME=GCA_036365495.1
    CLADE=sharks
Submitted batch job 5765953

Job 5765953 completed for paf dotplot of GCA_036365495.1 Heterodontus francisci (horn shark)
It worked!
Prepared 20250212_pafr_dotplot.sh for GCA_030028105.1 Mobula birostris (Giant manta)
    TODAY_DATE=20250227
    REF_NAME=GCA_030028105.1
    CLADE=sharks
Submitted batch job 5768824


#### UPDATE ####
20250228 (February 28th, 2025)

Job 5749589 completed for FASTGA alignment of GCA_031893055.1.fa.gz Leptodactylus fuscus (rufous frog) (alternate = GCA_031893025.1.fa.gz)
It worked!
Prepared 20250129_FASTGA_alignment.sh for GCA_038048845.1.fa.gz Mixophyes fleayi (Fleay's barred frog) (alternate = GCA_038048865.1.fa.gz)
    TODAY_DATE=20250228
    REF_NAME=GCA_038048845.1
    ALT_NAME=GCA_038048865.1
    CLADE=amphibians
    SPEC_NAME=MixFle
Submitted batch job 5824537

Job 5768824 completed for generating a dotplot of GCA_030028105.1 Mobula birostris (Giant manta)
It worked!
Prepared 20250212_pafr_dotplot.sh for GCA_036971175.1 Narcine bancroftii (Spotted torpedo ray/Caribbean electric ray)
    TODAY_DATE=20250228
    REF_NAME=GCA_036971175.1
    CLADE=sharks
Submitted batch job 5824668

Job 5765892 failed for looking at BAM plot coverage of GCA_036971175.1 Narcine bancroftii (Spotted torpedo ray/Caribbean electric ray)
It didn't work -- output file was generated with a file size of zero
Error:
    [E::parse_cigar] CIGAR length too long at position 1 (345549512S)
    [W::sam_read1_sam] Parse error at line 448
    samtools sort: truncated file. Aborting
    [E::hts_open_format] Failed to open file "sharks/GCA_036971175.1/20250227_NarBan.bam" : No such file or directory
    samtools coverage: Could not open "sharks/GCA_036971175.1/20250227_NarBan.bam": No such file or directory

Prepared 20250211_make_bam_plot_coverage.sh for GCF_030144855.1 Hypanus sabinus (Atlantic stingray)
    TODAY_DATE=20250228
    CLADE=sharks
    SPECIES_NAME=HypSab
    REF_NAME=GCF_030144855.1 ## Change for reference genome of each species
    ALT_NAME=GCA_030144785.1
Submitted batch job 5825146

Altered 20250219_ROH_Calc_V2.sh to submit a test version of 20250227_ROH_Durbin_Calc_Eqns_V3.R to see if it works
If it does I will implement it in the array
Submitted batch job 5825265 

Job 5824668 completed for dotplot of GCA_036971175.1 Narcine bancroftii (Spotted torpedo ray/Caribbean electric ray)
It worked!
Prepared 20250212_pafr_dotplot.sh for GCF_030144855.1 Hypanus sabinus (Atlantic stingray)
    TODAY_DATE=20250228
    REF_NAME=GCF_030144855.1
    CLADE=sharks
Submitted batch job 5825675

Job 5824537 completed for FASTGA alignment of GCA_038048845.1.fa.gz Mixophyes fleayi (Fleay's barred frog) (alternate = GCA_038048865.1.fa.gz)
It worked!

Job 5825675 completed for dotplot of GCF_030144855.1 Hypanus sabinus (Atlantic stingray)
It worked!
Prepared 20250212_pafr_dotplot.sh for GCA_027789765.1.fa.gz Dendropsophus ebraccatus (hourglass treefrog) (alternate = GCA_027789725.1.fa.gz)
    TODAY_DATE=20250228
    REF_NAME=GCA_027789765.1
    CLADE=amphibians
Submitted batch job 5839547

Job 5825265 failed for ROH calculation
Due to out of memory error
Upped memory to 200GB and resubmitted
Submitted batch job 5839970

Job 5825146 completed for BAM plot coverage of GCF_030144855.1 Hypanus sabinus (Atlantic stingray)
It worked!
Prepared 20250211_make_bam_plot_coverage.sh for GCA_027789765.1.fa.gz Dendropsophus ebraccatus (hourglass treefrog) (alternate = GCA_027789725.1.fa.gz)
    TODAY_DATE=20250228
    CLADE=amphibians
    SPECIES_NAME=DenEbr
    REF_NAME=GCA_027789765.1 ## Change for reference genome of each species
    ALT_NAME=GCA_027789725.1
Submitted batch job 5840158

Prepared 20250129_FASTGA_alignment.sh for GCA_027917425.1.fa.gz Gastrophryne carolinensis (eastern narrow-mouthed toad) (alternate = GCA_027917415.1.fa.gz)
    TODAY_DATE=20250228
    REF_NAME=GCA_027917425.1
    ALT_NAME=GCA_027917415.1
    CLADE=amphibians
    SPEC_NAME=GasCar
Submitted batch job 5840308

Job 5839547 completed for dotplot of GCA_027789765.1.fa.gz Dendropsophus ebraccatus (hourglass treefrog) (alternate = GCA_027789725.1.fa.gz
It worked!
Prepared 20250212_pafr_dotplot.sh for GCA_031893055.1.fa.gz Leptodactylus fuscus (rufous frog) (alternate = GCA_031893025.1.fa.gz)
    TODAY_DATE=20250228
    REF_NAME=GCA_031893055.1
    CLADE=amphibians
Submitted batch job 5841075

Job 5841075 completed for dotplot of GCA_031893055.1.fa.gz Leptodactylus fuscus (rufous frog) (alternate = GCA_031893025.1.fa.gz)
It worked!
Prepared 20250212_pafr_dotplot.sh for GCA_038048845.1.fa.gz Mixophyes fleayi (Fleay's barred frog) (alternate = GCA_038048865.1.fa.gz)
    TODAY_DATE=20250228
    REF_NAME=GCA_038048845.1
    CLADE=amphibians
Submitted batch job 5844044

Job 5844044 completed for dotplot of GCA_038048845.1.fa.gz Mixophyes fleayi (Fleay's barred frog) (alternate = GCA_038048865.1.fa.gz)
It worked!

#### UPDATE ####
20250303 (March 3rd, 2025)

Job 5839970 failed for ROH calculation
Failed due to out of memory error
    slurmstepd: error: Detected 1 oom_kill event in StepId=5839970.batch. Some of the step tasks have been OOM Killed.

Job 5840158 completed for BAM file and coverage of GCA_027789765.1.fa.gz Dendropsophus ebraccatus (hourglass treefrog) (alternate = GCA_027789725.1.fa.gz)
It worked!

Prepared 20250211_make_bam_plot_coverage.sh for GCA_031893055.1.fa.gz Leptodactylus fuscus (rufous frog) (alternate = GCA_031893025.1.fa.gz)
    TODAY_DATE=20250303
    CLADE=amphibians
    SPECIES_NAME=LepFus
    REF_NAME=GCA_031893055.1 ## Change for reference genome of each species
    ALT_NAME=GCA_031893025.1
Submitted batch job 6012303

Job 5840308 completed for FASTGA alingment of GCA_027917425.1.fa.gz Gastrophryne carolinensis (eastern narrow-mouthed toad) (alternate = GCA_027917415.1.fa.gz)
It worked!

Prepared 20250129_FASTGA_alignment.sh for GCF_028390025.1.fa.gz Pseudophryne corroboree (corroboree frog) (alternate = GCA_028390055.1.fa.gz)
    TODAY_DATE=20250303
    REF_NAME=GCF_028390025.1
    ALT_NAME=GCA_028390055.1
    CLADE=amphibians
    SPEC_NAME=PseCor
Submitted batch job 6013109

Prepared 20250212_pafr_dotplot.sh for GCF_030867095.1 Alligator mississippiensis (American alligator) (alternate = GCA_030867065.1.fa.gz)
    TODAY_DATE=20250303
    REF_NAME=GCF_030867095.1
    CLADE=reptiles
Submitted batch job 6013481

Created 20250303_PAF_to_VCF.sh to convert FASTGA generated PAF files to VCF files for ROH calculations
Prepared for GCF_020745735.1 -- Hemiscyllium ocellatum (Epaulette shark)
    TODAY_DATE=20250303
    PAF_DATE=20250217
    CLADE=sharks
    SPECIES=HemOce
    REF_NAME=GCF_020745735.1 
Submitted batch job 6021129

Job 6021129 failed due to the inability to open file
Error:
    paftools.js:413: Error: k8_open: failed to open file
        var file = args[getopt.ind] == '-'? new File() : new File(args[getopt.ind]);
I believe error was due to syntax error in PAF variable, with file type written in uppercase as .PAF
Fixed
    PAF=${SPECIES_DIR}/${PAF_DATE}_${SPECIES}_FASTGA.paf
Removed empty output file
    rm -f -r 20250303_HemOce_FASTGA.vcf
Resubmitted
Submitted batch job 6021591

Job 6021591 finished for FASTGA alignment VCF file of GCF_020745735.1 -- Hemiscyllium ocellatum (Epaulette shark)
Job did not work! It did not find any variants
The problem is due to paftools.js pulling variants from cs string, which FASTGA does not generate when going from ALNtoPAF
Tested to see if this changed based on cigar string syntax used
    ALNtoPAF -m temp/20250203_HemOce_FASTGA.1aln > 20250303_TEST_HemOce_FASTGA.paf
This PAF file also does not have cs strings

Prepared 20241231_ROH_Calc.sh to run for FASTGA generated alignment of GCF_020745735.1 -- Hemiscyllium ocellatum (Epaulette shark)
Added in code to create VCF file from FASTGA generated PAF
    TODAY_DATE=20250303
    VCF_DATE=20250303
    CLADE=sharks
    REF_NAME=GCF_020745735.1 
    SPECIES=HemOce


#### UPDATE ####
20250304 (March 4th, 2025)

Job 6013109 completed for FASTGA alingment of GCF_028390025.1.fa.gz Pseudophryne corroboree (corroboree frog) (alternate = GCA_028390055.1.fa.gz)
It worked!

Prepared 20250129_FASTGA_alignment.sh for GCF_029499605.1.fa.gz Hyla sarda (Sardinian treefrog) (alternate = GCA_029493135.1.fa.gz)
    TODAY_DATE=20250304
    REF_NAME=GCF_029499605.1
    ALT_NAME=GCA_029493135.1
    CLADE=amphibians
    SPEC_NAME=HylSar
Submitted batch job 6068739

Job 6068739 completed for FastGA alignment of GCF_029499605.1.fa.gz Hyla sarda (Sardinian treefrog) (alternate = GCA_029493135.1.fa.gz)
It worked!

Job 6012303 failed for BAM file and coverage of GCA_031893055.1.fa.gz Leptodactylus fuscus (rufous frog) (alternate = GCA_031893025.1.fa.gz)
Error:
    [E::parse_cigar] CIGAR length too long at position 1 (343113863H)
    [W::sam_read1_sam] Parse error at line 1271
    samtools sort: truncated file. Aborting
    [E::hts_open_format] Failed to open file "amphibians/GCA_031893055.1/20250303_LepFus.bam" : No such file or directory
    samtools coverage: Could not open "amphibians/GCA_031893055.1/20250303_LepFus.bam": No such file or directory

Prepared 20250211_make_bam_plot_coverage.sh for GCA_038048845.1.fa.gz Mixophyes fleayi (Fleay's barred frog) (alternate = GCA_038048865.1.fa.gz)
    TODAY_DATE=20250304
    CLADE=amphibians
    SPECIES_NAME=MixFle
    REF_NAME=GCA_038048845.1 ## Change for reference genome of each species
    ALT_NAME=GCA_038048865.1
Submitted batch job 6068887

Job 6013481 completed for making a dotplot of GCF_030867095.1 Alligator mississippiensis (American alligator) (alternate = GCA_030867065.1.fa.gz)
It worked!

Prepared 20250212_pafr_dotplot.sh for GCA_030020295.1.fa.gz Gavialis gangeticus (Gharial) (altnerate = GCA_030020385.1.fa.gz)
    TODAY_DATE=20250304
    REF_NAME=GCA_030020295.1
    CLADE=reptiles
Submitted batch job 6069084

Job 6069084 completed for the dotplot of GCA_030020295.1.fa.gz Gavialis gangeticus (Gharial) (altnerate = GCA_030020385.1.fa.gz)
It worked!

Prepared 20250212_pafr_dotplot.sh for GCA_033349115.1.fa.gz Macrochelys suwanniensis (Suwannee snapping turtle) (alternate = GCA_033296515.1.fa.gz)
    TODAY_DATE=20250304
    REF_NAME=GCA_033349115.1
    CLADE=reptiles
Submitted batch job 6075810


#### UPDATE ####
20250307 (March 7th, 2025)

Submitted 20250305_ROH_Calc_V4.sh 
Script calls on 20250304_ROH_Durbin_Calc_Eqns_V4.py
Running as a test on a single chromosome -- CM068742.1
    TODAY_DATE=20250307
    CLADE=sharks ## Clade directory
    REF_NAME=GCA_035084275.1 ## Reference genome name
    CHROM_START_CHR=CM
Submitted batch job 6266957

Canceled 6266957 due to the amount of time it was taking -- Need to make it more efficient


#### UPDATE ####
20250309 (March 9th, 2025)

Commented out ROH calculation function in 20250304_ROH_Durbin_Calc_Eqns_V4.py
Bash running 20250304_ROH_Durbin_Calc_Eqns_V4.py
Evalulating how long the base sums function takes
Started running at 2:00PM 
Job only takes ~2-4 minutes, not very long -- time seems to be coming from ROH calculation script
length_vcf = 287968
testing calculate_ROH function time to run by setting for loop to only run for 100 variants
Submitted at 2:26PM
Finished at 2:29PM -- finished  with error
    ValueError: Expected 1D or 2D array, got 0D array instead
Tested and it appeared to work, but just found no ROH in that area
Changed to look at chromosome CM071011.1 in GCA_036365495.1 Heterodontus francisci (horn shark)
This chromosome is significantly shorter and appears to be rife with ROH from previous calculations
Submitted batch job 6352500

Job 6352500 failed
    Failed out with error:
        base_sums[1,k] = base_sums[1,k]+1
    IndexError: index 9100128 is out of bounds for axis 1 with size 9100128
Changed 
    num_bases = range(0,chrom_length, 1) to 
    num_bases = range(1,chrom_length+1, 1)
Submitted batch job 6352567

Job 6352567 failed due to the same error
      File "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/20250304_ROH_Durbin_Calc_Eqns_V4.py", line 108, in cg_to_alingment
        base_sums[1,k] = base_sums[1,k]+1
    IndexError: index 9100128 is out of bounds for axis 1 with size 9100128


#### UPDATE ####
20250310 (March 10th, 2025)

Troubleshooting 20250305_ROH_Calc_V4.sh
Confirmed everything uptil for k in range... loop is working as it should
Problem is coming somewhere from counting base_pos in the if aln_type.../for k in range... statement
Position of an alingment is going beyond end of reference chromosome in two specific cases: cigar 303 and cg_array 283 and 284
Hypothetically have it working -- have had to add an clause in line 116 which I'm not happy about 
                if position > chrom_length or next_position > chrom_length:
                print(i, k)
Started running at 12:52-12:53 with ROH calculation function added to see how long it takes
Returned with an empty array -- I'm thinking this is because of the low alingment amount that it didn't find any ROH


#### UPDATE ####
20250313 (March 13th, 2025)

Got 20250311_ROH_Calc_V5.sh to work, with 20250311_ROH_Durbin_Calc_Eqns_V5.py
Worked for single chromosome on test with GCA_036365495.1 Heterodontus francisci (horn shark)
    TODAY_DATE=20250313
    CLADE=sharks ## Clade directory
    REF_NAME=GCA_036365495.1 ## Reference genome name
    CHROM_START_CHR=CM
    CHROM=CM071011.1

Now will test for whole genome of GCA_036365495.1 Heterodontus francisci (horn shark)
    #SBATCH --array=1-51
    TODAY_DATE=20250313
    CLADE=sharks ## Clade directory
    REF_NAME=GCA_036365495.1 ## Reference genome name
    CHROM_START_CHR=CM
Submitted batch job 6561561

Job 6561561 canceled due to error on first run:
    FileNotFoundError: [Errno 2] No such file or directory: '/20250313/_/CM071011.1/_ROH_Results.txt'
Changed output_name line syntax to fix:
    output_name = os.path.join(clade, ref_name, today_dat + "_" + chrom + "_ROH_Results.txt")
Resubmitted
Submitted batch job 6561862

Job 6561862 canceled due to error on first array run
    from readpaf import parse_paf
    ModuleNotFoundError: No module named 'readpaf'
Removed command due to the fact that it is not needed in this version of ROH calculation script
Resubmitted
Submitted batch job 6562086

Starting to work on snakemake for workflow automation
    module load ceuadmin/snakemake/7.19.1
Working on snakemake tutorial (https://snakemake.readthedocs.io/en/stable/tutorial/setup.html)
First downloaded:
    conda create -c conda-forge -c bioconda -n snakemake snakemake
Then updated Conda
    conda update -n base -c defaults conda
Created SnakeFile
Will test if snakefile works by creating FASTGA alingment for GCA_035609145.1.fa.gz Eleutherodactylus coqui (Puerto Rican coqui) (alternate = GCA_035609135.1.fa.gz)
    rule fastga_alignment:
    output:
        "amphibians/EleCoq/temp/EleCoq_GCA_035609145.1.1gdb"
    shell:
        """
        FAtoGDB -v ../241117.UCSC-hubs-VGP-alignment/alignment/reference/amphibians/GCA_035609145.1.fa.gz \
        amphibians/EleCoq/temp/EleCoq_GCA_035609145.1.1gdb
        """
Got this test to work -- Snakemake functions
    

#### UPDATE ####
20250314 (March 17th, 2025)

Job 6562086 ran out of time -- did not successfully generate ROH files for the whole genome
Updated to run for a single chromosome, to make sure that pulling files and information in shell script works
    #SBATCH --array=1
    CHROM=CM071010.1

    conda activate msmc_env
Submitted batch job 6738625

Job 6738625 completed successfully
No ROH detected on this chromosome

Given that this script appears to be working, I think the issue was computation time/memory with a large genome
Prepare 20250311_ROH_Calc_V5.sh for GCA_035084275.1 Hydrolagus colliei (Spotted ratfish)
This genome is smaller than that for HetFra (Horn Shark) both in chromosome size and number of chromosome
    #SBATCH --array=1-40 
    TODAY_DATE=20250317
    CLADE=sharks ## Clade directory
    REF_NAME=GCA_035084275.1 ## Reference genome name
    CHROM_START_CHR=CM
Upped memory to 150GB and runtime to 16 hours
Submitted batch job 6739060


Working on snakemake file
Set up rules to run FASTGA alingment on GCA_035609145.1.fa.gz Eleutherodactylus coqui (Puerto Rican coqui) (alternate = GCA_035609135.1.fa.gz)
It worked!


#### UPDATE ####
20250318 (March 18th, 2025)

Job 6739060 completed for ROH calculation of GCA_035084275.1 Hydrolagus colliei (Spotted ratfish)
It worked!

Prepared 20250106_Plot_ROH_FROH.sh for GCA_035084275.1 Hydrolagus colliei (Spotted ratfish)
    TODAY_DATE=20250318
    ROH_CALC_DATE=20250317
    CLADE=sharks
    REF_NAME=GCA_035084275.1 ## Change for reference genome of each species
    NUM_AUT_CHROMOSOMES=40
    CHROM_START_CHR=CM
Submitted batch job 6794682

Job completed, but did not generate any output due to syntax error when reading in ROH files
Fixed syntax error and resubmitted
Submitted batch job 6795070

Job finished, but ran incorrectly
.tsv file containing ROH results only contains the names of chromosomes and no ROH
Error when run in bash:
    'names' attribute [4] must be the same length as the vector [1]


#### UPDATE ####
20250319 (March 19th, 2025)

Working on adding rules for calculating FROH to the snakefile
Testing with chromosome CM068742.1 from GCA_035084275.1 Hydrolagus colliei (Spotted ratfish)
Working with minimap2 generated paf file -- GCA_035084275.1.mm2.srt.paf in sharks/GCA_035084275.1 directory
    awk -v var=CM068742.1 '{if ($6 == var && $12 >= 60) {print $0}}' GCA_035084275.1.mm2.srt.paf > 20250319_CM068742.1_filtered.paf
    k8 paftools.js call sharks/GCA_035084275.1/20250319_CM068742.1_filtered.paf > sharks/GCA_035084275.1/20250319_CM068742.1_Aln_Var.txt
    awk -v var='V' '{if ($1 == var && $5 == 1) {print $0}}' 20250319_CM068742.1_Aln_Var.txt > 20250319_CM068742.1_Var_Only.txt

    head -n 1 20250319_CM068742.1_Var_Only.txt | awk ' {print $3} ' ## This command prints the start of the first variant
    tail -n 1 20250319_CM068742.1_Var_Only.txt | awk ' {print $4} ' ## This command prints the end of the last variant

    awk -v var='R' '{if ($1 == var) {print $0}}' 20250319_CM068742.1_Aln_Var.txt > 20250319_CM068742_Aln_Only.txt


#### UPDATE ####
20250320 (March 20th, 2025)

Created 20250320_ROH_Durbin_Calc_Eqns_V6.py -- based on Version 5 but hopefully running more efficiently for large genomes
Running with 20250311_ROH_Calc_V5.sh shell script
Testing with GCA_035084275.1 Hydrolagus colliei (Spotted ratfish)
Submitted batch job 7009199

Potentially will need to optimize 20250115_find_het_calc.R script
Currently rerunning with no changes to see how it works
Submitted 20250114_find_het_per_chr.sh for GCA_036971175.1 Narcine bancroftii (Spotted torpedo ray/Caribbean electric ray)
Submitted batch job 7009693

Canceled job 7009693 since it was not running due to syntax errors in the shell script


#### UPDATE ####
20250321 (March 21st, 2025)

Fixed syntax errors in 20250311_ROH_Calc_V5.sh and Resubmitted
Submitted batch job 7062643

Canceled job since it hadn't started running after several hours
Switched to icelake-himem node and resubmitted
Submitted batch job 7077974

Job 7077974 finished with mixed results
17 ROH results files were generated, and no error files
Possible that mixed results were that no ROH were found on some chromosomes
I need to have something to indicate that instead of a legit error
Added if statement to 20250320_ROH_Durbin_Calc_Eqns_V6.py to print in out file if no ROH are found
    if ROH_report.empty == True:
        print('No ROH detected in', chrom)
Also added command in shell script to remove temporary files
Resubmitted
Submitted batch job 7083343

Job 7083343 completed without errors
It appears to have worked!
Came back with only 17 ROH files again

Finished writing 20250319_find_het_calc_V2.py and shell script to execute: 20250319_find_het_calc_V2.sh
This should run faster than version 1
Testing with GCA_035084275.1 Hydrolagus colliei (Spotted ratfish)
    TODAY_DATE=20250321
    ROH_DATE=20250118
    CLADE=sharks
    REF_NAME=GCA_035084275.1 ## Change for reference genome of each species
    NUM_AUTOSOMAL_CHROMOSOMES=40 ## Change for each species -- will have to manually check on NCBI
    WINDOW_LENGTH=1000000 ## How long we want each window that we are using to measure heterozygosity
    WINDOW_INTERVAL=500000
Submitted with bash
Need to fix error getting thrown in line 118:
    within_roh = roh_in_range & (starts >= roh_starts[roh_indices]) & (ends <= roh_ends[roh_indices]) ## FIX ERROR GETTING THROWN HERE

    ValueError: Can only compare identically-labeled Series objects


#### UPDATE ####
20250324 (March 24th, 2025)

Troubleshooting 20250319_find_het_calc_V2.py
Error finding variants on other chromosomes -- I have confirmed that they are present
Hypothetically have finished troubleshooting script and it is now working

Added heterozygosity calculation rule to snakefile


#### UPDATE ####
20250325 (March 25th, 2025)

Due to rules of snakemake, had to split heterozygosity calculations per chromosome and over whole genome into two separate files
Created 20250325_find_het_per_chr_V3.py and 20250325_find_het_whole_genome_V3.py


#### UPDATE ####
20250326 (March 26th, 2025)

I ran the snakefile, starting with just running the commands to generate the .paf file
    snakemake ALNtoPAF

    CLADE="amphibians"
    SPEC_NAME="EleCoq_TEST"
    REF_NAME="GCA_035609145.1"
    ALT_NAME="GCA_035609135.1"
    TODAY_DATE="20250325"
    CHROM_START_CHR="CM"
    CHROM_LIST_FILE="/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/amphibians/chrom_lists/GCA_035609145.1_chroms.txt"
    WINDOW_INTERVAL=500000
    WINDOW_LENGTH=1000000
    NUM_AUT_CHROMOSOMES=13
It worked!

Created config.yml yaml file to move configuration settings there instead of in snakefile
Have to add all chromosome in there manually

Testing snakefile with 6 chromosomes

Command to print DAG
    snakemake --configfile config.yml --dag | dot -Tpdf > dag.pdf

Having error where no variants are getting tracked in the files generated by k8 paftools.js
I'm going to try running that same command on the mm2 and FASTGA alingments of GCA_035084275.1 Hydrolagus colliei (Spotted ratfish) to see if I can find the error
Running it on the mm2 generated .paf file (GCA_035084275.1.mm2.srt.paf) I find variants
Running it on the FASTGA .paf file (HydCol/20250220_HydCol_FASTGA.paf) generates the same issue

Will continue testing snakefile working using GCA_035084275.1 Hydrolagus colliei (Spotted ratfish)
Using first five chromosomes, and will use minimap2 alingment
    CLADE : "sharks"
    SPEC_NAME : "HydCol_Test"
    REF_NAME : "GCA_035084275.1"
    ALT_NAME : "GCA_035084065.1"
    TODAY_DATE : "20250326"
    CHROM_START_CHR : "CM"
    CHROM_LIST_FILE : "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/sharks/chrom_lists/GCA_035084275.1_chroms.txt"
    WINDOW_INTERVAL : 500000
    WINDOW_LENGTH : 1000000
    NUM_AUT_CHROMOSOMES : 40
    Chromosomes : 
    - CM068742.1
    - CM068743.1
    - CM068744.1
    - CM068745.1
    - CM068746.1

Manually changed input of FILTER_PAF rule on line 219 to work with mm2 alingment for testing
    input:
        f"{CLADE}/{SPEC_NAME}/GCA_035084275.1.mm2.srt.paf"

#### UPDATE ####
20250331 (March 31st, 2025)

Did some troubleshooting of WHOLE_FROH rule and altered code in 20250106_Plot_ROH.R to see if it will save a pdf properly
Ran snakemake file

Realized that PAF_DOTPLOT rule was not working because FASTGA generated .paf files do not have tp strings for primary vs secondary alignments
Changed for mm2 paf file and retried
Dotplot successfully made!

Got ROH plot to successfully save!
Deleted HydCol_Test directory and remade it, copying over mm2.srt.paf file for HydCol from GCA_035084275.1
Will now test for all chromosomes
Updated config.yml file with all chromosomes
    CLADE : "sharks"
    SPEC_NAME : "HydCol_Test"
    REF_NAME : "GCA_035084275.1"
    ALT_NAME : "GCA_035084065.1"
    TODAY_DATE : "20250331"
    CHROM_START_CHR : "CM"
    CHROM_LIST_FILE : "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/sharks/chrom_lists/GCA_035084275.1_chroms.txt"
    WINDOW_INTERVAL : 500000
    WINDOW_LENGTH : 1000000
    NUM_AUT_CHROMOSOMES : 40
    Chromosomes : 
    - CM068742.1
    - CM068743.1
    - CM068744.1
    - CM068745.1
    - CM068746.1
    - CM068747.1
    - CM068748.1
    - CM068749.1
    - CM068750.1
    - CM068751.1
    - CM068752.1
    - CM068753.1
    - CM068754.1
    - CM068755.1
    - CM068756.1
    - CM068757.1
    - CM068758.1
    - CM068759.1
    - CM068760.1
    - CM068761.1
    - CM068762.1
    - CM068763.1
    - CM068764.1
    - CM068765.1
    - CM068766.1
    - CM068767.1
    - CM068768.1
    - CM068769.1
    - CM068770.1
    - CM068771.1
    - CM068772.1
    - CM068773.1
    - CM068774.1
    - CM068775.1
    - CM068776.1
    - CM068777.1
    - CM068778.1
    - CM068779.1
    - CM068780.1
    - CM068781.1 
Ran with snakemake --configfile config.yml

It worked!!!


#### UPDATE ####
20250401 (April 1st, 2025)

When looking at FROH results for HydCol_Test, I noticed the results were odd for Froh_aut for all ROH
    "Froh_aut value for all ROH (short, med, long) is:"
    399.4996
    "Froh_aut in percent for all ROH (short, med, long)"
    "39949.9589344003 %"

Need to look into snakefile, and find cause for these odd results and fix it
Additionally realized value of Froh_aut and Froh_aut for autosomal ROH are different
In this case, no sex chromosomes were listed so they should be identical for HydCol, but here they are not
    echo "$Laut_autosomal"
    echo "$Laut"
Returns 937287621 and 47485, respectively
Laut is much smaller than it should be, Laut_autosomal looks ~correct
Fixed by calculating them the same way -- having Laut calculated by going through each chromosome and finding start and end and summing them in a while loop

Created 20250401_Reference_Clade_Plot.py --> will use this script to plot reference genome clade and order affiliation
The script needs further work when I have an updated taxa list, especially to make the graph better visually, but baseline is functional


#### UPDATE ####
20250402 (April 2nd, 2025)

Created 20250402_cs_string_test.py
Will use this script to test my ability to make cs strings from cigar strings and store them in paf files
Created sample paf file with only four entries to use in testing this
    sed 4q HydCol_Test_FASTGA.paf > cs_creation_test.paf
Got it working -- created test output output_cs_creation_test.paf
Next tested with the HydCol_Test_FASTGA.paf file
    It worked!
Testing to see if it generates variants when run through paftools.js
    k8 paftools.js call sharks/HydCol_Test/output_cs_FASTGA.paf > sharks/HydCol_Test/FASTGA_Test_Aln_Var.txt

I noticed in 20250331_HydCol_Test_whole_genome_mean_heterozygosity.txt that the heterozygosity values for with and without ROH were identical -- need to look into this and see what is going on
Checked the other output file from CALC_HET_WHOLE_GENOME rule and results are not coming out correctly
Problem with ...._per_chr_mean_heterozygosity.txt file is that the dataframe is appended rather than filling in the empty cells

Fixed error in generating ...._per_chr_mean_heterozygosity.txt file -- Updated python code
    mean_het_df.loc[mean_het_df['chr']==chr_name,'mean_het'] = np.mean(chrom_file['Het_Per_Kb'])
    mean_het_df.loc[mean_het_df['chr']==chr_name,'mean_het_excl_ROH'] = np.mean(chrom_file['Het_Per_Kb_excl_ROH'])

I think error for het and het_excl_ROH coming out the same is being sourced from calculating heterozygosity per chromosome
Looking at 20250325_find_het_per_chr_V3


#### UPDATE ####
20250403 (April 3rd, 2025)

Created 20250403_cs_string_generator.py
This is the cleaned up version of the script to make cs strings for paf files
This will be used in snakefile
Based on final set of testing from 20250402_cs_string_test.py

Added PAF_CS rule to snakefile
This rule will take FASTGA generated .paf file and create cs strings for it
Updated rule_all to need the output of PAF_CS   
    f"{CLADE}/{SPEC_NAME}/{SPEC_NAME}_FASTGA_Updated.paf"

Updated config.yml file to run test and see if cs string generation works
    SPEC_NAME : "HydCol_Test2"
    TODAY_DATE : "20250403"
Ran with snakemake --configfile config.yml

Job stopped due to error with CALC_HET_WHOLE_GENOME
CALC_HET_WHOLE_GENOME tries to run before CALC_HET_PER_CHR -- cannot work without necessary files
I think I will add output of CALC_HET_PER_CHR as input for CALC_HET_WHOLE_GENOME, and just not call it in the shell
    PER_CHR_FILES=expand("{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{CHR}_het.txt", CLADE=CLADE, SPEC_NAME=SPEC_NAME, TODAY_DATE=TODAY_DATE, CHR=CHROMS)
Resubmitted

Altered the 20250212_pafr_dotplot.R script so that if the paf file doesn't have the tp tag, it will run without filtering for primary alignments
Did this because currently FASTGA generated PAF files do not have the tp tag


#### UPDATE ####
20250404 (April 4th, 2025)

In some but not all of the ROH calculations, I'm getting an index out of bounds error
Happening in line 100 of 20250320_ROH_Durbin_Calc_Eqns_V6.py 
    df['var_aln_sum'] = base_aln_map[2, indices]
Hypothetically fixed by subtracting one from the indices in the line above
    indices = np.searchsorted(base_aln_map[0, :], df['variant']) - 1

Analysis seems to be taking quite a long time -- it's the generation of the whole genome Aln_Var.txt file which is taking a long time
The reason it's taking so long is because it's also mapping things to the scaffolds
Need to alter script so it is only mapping alingments and variants to chromosomes as reference sequences

Adding rule to snakefile to filter FASTGA file before we add CS strings to only include chromosomes as reference for alingments
Added rule FILTER_PAF_CHR_ONLY

Removed HydCol_Test2 directory and resubmitted SnakeMake, first only to generate updated PAF file to see if it worked
It seems to have worked!

Testing snakemake on next set of code -- to get Var and Aln files
Step to get whole genome Aln_Var file is still taking a long time -- the file still has scaffolds as reference sequences in addition to chromosomes
Trying an alteration of the script to run faster 
    k8 paftools.js call {input} | awk '$1 == "V" && $2 !~ /^JAY/' > {output}

Decided to alter snakefile to avoid making unnecessary files

Noticed when re-running that 10/40 Chromosome specific Aln_Only files have a file size of zero
This indicates that there were no alingments found on these chromosomes, yet variants are found
I also think something else odd is going on because some of the variant files for a chromosome are larger than that for the whole genome

Tested GET_ALN_ONLY_PER_CHROM using mm2 generated paf file 
    input = f"{CLADE}/GCA_035084275.1/GCA_035084275.1.mm2.srt.paf"
It went incredibly quickly, and all chromosomes had alignments

I found the problem with the heterozygosity code --> the ROH data is not getting read properly
Fixed this in the filter_roh_by_chrom function in 20250325_find_het_per_chr_V3.py
However results now show Het_excl_ROH is the same for every single window along a chromosome
Fixed this by altering for loop starting at line 122 in 20250325_find_het_per_chr_V3.py


#### UPDATE ####
20250407 (April 7th, 2025)

Will run test of snakemake file on GCA_036971175.1 Narcine bancroftii (Spotted torpedo ray/Caribbean electric ray)
See if I get the same issues with the alingment and variant txt file generation as I did with HydCol_Test
Prepared config.yml for NarBan_Test
    CLADE : "sharks"
    SPEC_NAME : "NarBan_Test"
    REF_NAME : "GCA_036971175.1"
    ALT_NAME : "GCA_036971445.1"
    TODAY_DATE : "20250407"
    CHROM_START_CHR : "CM"
    CHROM_LIST_FILE : "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/sharks/chrom_lists/GCA_036971175.1_chroms.txt"
    WINDOW_INTERVAL : 500000
    WINDOW_LENGTH : 1000000
    NUM_AUT_CHROMOSOMES : 14
    Chromosomes : 
    - CM072956.1
    - CM072957.1
    - CM072958.1
    - CM072959.1
    - CM072960.1
    - CM072961.1
    - CM072962.1
    - CM072963.1
    - CM072964.1
    - CM072965.1
    - CM072966.1
    - CM072967.1
    - CM072968.1
    - CM072969.1

Started with only the rules to generate the paf file with cs strings added.
Stopped at the PAF_CS rule due to a memory error

Recompiled FASTGA after Chenxi pushed updates to package 
Update has ALNtoPAF now have option to create cs strings
Will test this with sharks/HydCol_Test/HydCol_Test_ALN.1aln
    ALNtoPAF -s -T8 HydCol_Test_ALN.1aln > 20250407_TEST_CS.paf
it worked! It generated the shortform CS strings but with no cigar strings
Will now test if paftools.js can work on this
    k8 ../../paftools.js call 20250407_TEST_CS.paf > 20250407_TEST_CM068782.1_ALN.txt

If this doesn't work -- try generating paf file again with cigar strings as well as cs strings -- both in short form


#### UPDATE ####
20250408 (April 8th, 2025)

paftools call on the paf file with the cs strings did not complete, but it appears to have worked
Included finding variants for scaffolds (JAYJKW.....) -- Will need to see if I can stop that from happening to make it run faster and reduce file size
    awk -v var="CM068743.1" '{if ($2 == var && $1 == "R") {print $0}}' 20250407_TEST_CM068782.1_ALN.txt
File appears to look normal and as expected
    awk '{if ($1 !~ /^CM0/) {print $0}}' 20250407_TEST_CS.paf
The paf file contains alignments to scaffolds which I don't think I need
I will filter them out and then run the paftools command to see if this works
    awk '{{ if ($6 !~ /^JAY/) {{print $0}}}}' 20250407_TEST_CS.paf > 20250407_TEST_CS_filtered.paf
It worked!
Testing command from GET_ALN_ONLY_PER_CHROM
    k8 ../../paftools.js call 20250407_TEST_CS_filtered.paf | awk '$1 == "R" && $2 == "CM068742.1"' > 20250407_TEST_CM068782_ALN.txt
Started running at 9:42AM
Canceled job at 9:49AM
Wondering if it would run faster to generate whole txt file with all aln and var and then filter it instead of piping for each
    k8 ../../paftools.js call 20250407_TEST_CS_filtered.paf > TEST_Total_Aln_Var.txt
Submitting at 9:50AM
Canceled at 10:04AM -- Wondering if it will run faster if supported
    sort -k6,6V -k8,8n 20250407_TEST_CS_filtered.paf > 20250407_TEST_CS_filtered.srt.paf
It finished
    k8 ../../paftools.js call 20250407_TEST_CS_filtered.srt.paf > TEST_Total_Aln_Var.txt
started running at 10:09AM
Canceled at 10:33AM
    k8 ../../paftools.js call 20250407_TEST_CS_filtered.srt.paf | awk '$1 == "R"' > TEST_Total_Aln.txt
Submitted at 10:34AM
It ran until 10:42, when the remote connection timed out to cluster
Got to chromosome CM068754.1

Also testing submitting it through a shell script to see if it runs faster
Created paftools.sh
Submitted batch job 7983878
Job timed out after 2 hours and did not complete

Attempting to run with vcf file generation to see if it will run faster
    ALNtoPAF -s -T8 -m sharks/HydCol_Test/HydCol_Test_ALN.1aln > TEST.paf
Ran at 12:24PM
Finished at 12:25PM
    k8 paftools.js call -s HAP2 -f ../241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/GCA_035084275.1.fa.gz TEST.paf > TEST.vcf
Started running at 12:36PM
Canceled at 12:49 when I saw that it hadn't created any new entries in the .txt file after 12:37PM

Updated Snakemake file rules for getting var and aln
First have GET_WHOLE_VAR and GET_WHOLE_ALN -- these rules generate txt files with all variants and alignments for whole genome

Then will filter these for specific chromosomes in GET_ALN_ONLY_PER_CHROM and GET_VAR_ONLY_PER_CHROM
Canceled job after it had been running for 5 hours and not finished
The alingment file sharks/HydCol_Test2/temp/HydCol_Test2_Aln_Only.txt finished a few minutes after job started
Both were stuck running on CM068756.1
Specifically variants kept getting added to this chromosome

I'm going to remove this chromosome from the PAF file and rerun jobs to see if something funky is going on with that specific chromosome
    awk '$6 != "CM068756.1"' HydCol_Test2_FASTGA_Filtered.srt.paf > HydCol_Test2_wo_CM068756.1_TEST.srt.paf
Now will run paftools on this
    k8 ../../paftools.js call HydCol_Test2_wo_CM068756.1_TEST.srt.paf > TEST_wo_CM068756.1_aln_var.txt
Started running at 6:48PM
Finished at 7:30 -- it appears to have worked?


#### UPDATE ####
20250409 (April 9th, 2025)

See if I can filter TEST_wo_CM068756.1_aln_var.txt for just Aln or Var
    awk '$1 == "R"' TEST_wo_CM068756.1_aln_var.txt > TEST_wo_CM068756.1_aln_only.txt
Took only a few seconds to run
    awk '$1 == "V"' TEST_wo_CM068756.1_aln_var.txt > TEST_wo_CM068756.1_var_only.txt
Took only a few seconds to run
Most efficient is generating whole file and then splitting it after the fact -- not running pipeline

Need to figure out what is going on with CM068756.1
Filtered paf file to only contain this chr
    awk '$6 == "CM068756.1"' HydCol_Test2_FASTGA.srt.paf > CM068756.1_only_PAF.paf
5383 lines in .paf file
Ran this file to generate aln_var.txt
    k8 ../../paftools.js call CM068756.1_only_PAF.paf > CM068756.1_only_aln_var.txt
Started at 11:21AM
Line count for aln_var.txt file for all but CM068756.1 --> 7578851
Line count for aln_var.txt file for just CM068756.1 --> at least 518433 (when job canceled)
Canceled job at 12:16PM, it was still running
    awk '$1 == "R"' CM068756.1_only_aln_var.txt | wc -l 
110 alingments
    awk '$1 == "R" && $2 == "CM068755.1"' TEST_wo_CM068756.1_aln_only.txt | wc -l
The chromosome immediately before (and of similar length) has 101 alingments -- not significantly different
    awk '$1 == "V"' CM068756.1_only_aln_var.txt | wc -l 
518324 variants in chr
    awk '$1 == "V" && $2 == "CM068755.1"' TEST_wo_CM068756.1_var_only.txt | wc -l
36438 variants in chromosome immediately before
Why are there significantly more variants found in CM068756.1? Something doesn't seem right

Went into GCA_035084275.1 directory to see how these numbers compare to mm2 generated aln and var files
    awk '$1 == "V" && $2 == "CM068756.1"' 20250321_Var_Only.txt | wc -l
38791 variants
    awk '$1 == "V" && $2 == "CM068755.1"' 20250321_Var_Only.txt | wc -l
35636 variants
    awk '$1 == "R" && $2 == "CM068756.1"' 20250321_Aln_Var.txt | wc -l
81 alingments
    awk '$1 == "R" && $2 == "CM068755.1"' 20250321_Aln_Var.txt | wc -l
70 alingments


Testing with another species to see if this same issue with a singular chromosome comes up
Set config.yml for NarBan (as in line 4282) and submitted just to generate the .paf files

Create aln only txt file for the abberant chromosome 
    awk '$1 == "R"' CM068756.1_only_aln_var.txt > CM068756.1_aln.txt
Confirm it has 110 lines -- same number as above
Last alignment is this one:
    R       CM068756.1      17527895        17528460
See if I can find this region in the .paf file for this chromosome 
    awk '{print $7}' CM068756.1_only_PAF.paf | uniq
This returns a single value:
    17630497
Is this the same in other chromosomes?
    awk '{print $7}' HydCol_Test2_wo_CM068756.1_TEST.srt.paf | uniq | wc -l 
    39


#### UPDATE ####
20250410 (April 10th, 2024)

NOTE: Have to find a way to update scaffold names in the FILTER_PAF_CHR_ONLY rule
Start of scaffold names was different between HydCol and NarBan

Successfully generated NarBan_Test_FASTGA_Filtered.srt.paf
Will test if same chromosome issue occurs
    k8 ../../paftools.js call NarBan_Test_FASTGA_Filtered.srt.paf > TEST_aln_var.txt
Started running at 9:46AM
Canceled job at 9:56AM -- it was still on the first chromosome and seemingly not moving forward, just generating variants in the same region
Remove first chromosome
    awk '$6 != "CM072956.1"' NarBan_Test_FASTGA_Filtered.srt.paf > TEST_wo_CM072956.1.paf
Rerun paftools
    k8 ../../paftools.js call TEST_wo_CM072956.1.paf > TEST_wo_CM072956_aln_var.txt
Canceled job after it ran for ~40min and was still on the first chromosome
    awk '$1 == "V"' TEST_wo_CM072956_aln_var.txt | wc -l
    188705

    awk '$1 == "R"' TEST_wo_CM072956_aln_var.txt | wc -l
    3

While it's running, will generate same file from mm2 file sharks/GCA_036971175.1/GCA_036971175.1.mm2.srt.paf
    k8 ../../paftools.js call GCA_036971175.1.mm2.srt.paf > TEST_mm2_aln_var.txt
finished in a few seconds
    awk '$1 == "V" && $2 == "CM072957.1"' TEST_mm2_aln_var.txt | wc -l
    318399
    awk '$1 == "R" && $2 == "CM072957.1"' TEST_mm2_aln_var.txt | wc -l
    826
Possible that the chromosomes are just incredibly long and I didn't give it enough time to finish

However I went back and compared the number of variants for CM072956.1 from the intial aln_var file to the one made by mm2
    awk '$1 == "V"' TEST_aln_var.txt | wc -l
902219
    awk '$1 == "V" && $2 == "CM072956.1"' TEST_mm2_aln_var.txt | wc -l
410169
Again an insanely inflated number of variants
    awk '$1 == "R"' TEST_aln_var.txt > TEST_aln_only_CM068756.1.txt
Last alingment in file:
    R       CM072956.1      13129647        13158556
Will I find same results as in HydCol below?
    awk '$6 == "CM072956.1" && $8 <= 13129647 && $9 >= 13158556' NarBan_Test_FASTGA.srt.paf | wc -l
1 alignment found
    CM072942.1      494489450       14668111        14910361        +       CM072956.1      495097020       12916635        13158556     215738  248483  255 
Maybe something is going odd with the cs strings?
    awk '{if ($6 == "CM072956.1" && $8 <= 13129647 && $9 >= 13158556) {print $15}}' NarBan_Test_FASTGA.srt.paf > aln_cs_string.txt
Looking through this string in python, filtering it for unique characters, it doesn't appear that any unusual things pop up
Just normal expected (+, -, :, *)



For finding problem alingment in HydCol:
Last recorded alignment in CM068756.1_aln.txt:
    R       CM068756.1      17527895        17528460
Looking in paf file 
    awk '$8 == 17527895' CM068756.1_only_PAF.paf
Nothing Returned
    awk '$8 == 17528460' CM068756.1_only_PAF.paf
Returned a single alignment
    JAYJKX010000405.1       420528  10754   112810  -       CM068756.1      17630497        17528460        17630493        82183104049 

    awk '$9 == 17528460' CM068756.1_only_PAF.paf
Nothing Returned
    awk '$8 <= 17527895 && $9 >= 17528460' CM068756.1_only_PAF.paf | wc -l
3 alignments are returned
    CM068809.1      9187931 1453649 1461102 +       CM068756.1      17630497        17520983        17528472        7358    7523    255
    CM068813.1      3656718 2112593 2120057 -       CM068756.1      17630497        17520983        17528472        7328    7546    255 
    CM068796.1      18076791        17945120        18015602        +       CM068756.1      17630497        17521698        17592184        58755   71730   255
I am wondering if issue comes from first two --> They are listed as having the same start and ending positions, but are opposite relative strands to each other
    tail -20 CM068756.1_only_aln_var.txt
Doesn't appear that the issue is from a specific type of variant (e.g. indel vs point mutation) 
    awk '{print $9}' CM068756.1_only_aln_var.txt | uniq | wc -l 
233 supposedly unique things alinging to chromosome
    awk '{print $9}' CM068756.1_only_aln_var.txt | uniq | awk '$1 == "CM068796.1"' | wc -l
There were repeats of a single chromosome, so I counted it
111 repeats of this in the 'unique' filter

Looked for variants starting at same position
    awk '$1 == "V" && $3 == 17587213' CM068756.1_only_aln_var.txt
Got these results:
    V       CM068756.1      17587213        17587214        562     255     c       t       CM068784.1      109974261       109974262       +
    V       CM068756.1      17587213        17587214        562     255     t       g       JAYJKX010000169.1       28173   28174   -
    V       CM068756.1      17587213        17587214        562     255     g       t       JAYJKX010000405.1       230877  230878  +
    V       CM068756.1      17587213        17587214        562     255     t       g       JAYJKX010000405.1       232354  232355  +
    V       CM068756.1      17587213        17587214        562     255     g       t       JAYJKX010000405.1       233674  233675  +
5 variants at same position on same chromosome, 3 of which are coming from the same scaffold
Tried at another point:
    awk '$1 == "V" && $3 == 17587193' CM068756.1_only_aln_var.txt
Got these results:
    V       CM068756.1      17587193        17587194        562     255     a       t       JAYJKX010000405.1       226766  226767  +
    V       CM068756.1      17587193        17587194        562     255     a       g       JAYJKX010000562.1       69918   69919   +
    V       CM068756.1      17587193        17587194        562     255     t       a       JAYJKX010000562.1       67919   67920   +
    V       CM068756.1      17587193        17587194        562     255     t       c       JAYJKX010000405.1       232914  232915  +
    V       CM068756.1      17587193        17587193        562     255     -       c       JAYJKX010000405.1       233654  233655  +
    V       CM068756.1      17587193        17587194        562     255     g       a       JAYJKX010000405.1       233655  233656  +
5 of the 6 variants start and end at the same position -- 3 of which are from one scaffold and the other 2 from another scaffold
What if I were to remove scaffolds from the paf file before calling paftools?
    awk '{if ($1 !~ /^JA/) {print $0}}' HydCol_Test2_FASTGA_Filtered.srt.paf > HydCol_Test2_wo_Scaffolds_TEST.srt.paf
Re-run paftools on this
    k8 ../../paftools.js call HydCol_Test2_wo_Scaffolds_TEST.srt.paf > TEST_wo_scaffolds_aln_var.txt
started at 1:03PM
Canceled job at 1:14PM -- Same problem as before
Removed HydCol_Test2_wo_Scaffolds_TEST.srt.paf & TEST_wo_scaffolds_aln_var.txt to save on file space

Looking at the last alignment listed versus last few variants listed:
    R       CM068756.1      17527895        17528460

    V       CM068756.1      17587193        17587193        562     255     -       c       JAYJKX010000405.1       233654  233655  +
    V       CM068756.1      17587193        17587194        562     255     g       a       JAYJKX010000405.1       233655  233656  +
    V       CM068756.1      17587197        17587198        562     255     g       -       JAYJKX010000405.1       233659  233659  +
    V       CM068756.1      17587200        17587201        562     255     a       g       JAYJKX010000405.1       233661  233662  +
    V       CM068756.1      17587204        17587205        562     255     c       a       JAYJKX010000405.1       233665  233666  +
    V       CM068756.1      17587207        17587208        562     255     g       t       JAYJKX010000405.1       233668  233669  +
    V       CM068756.1      17587213        17587214        562     255     g       t       JAYJKX010000405.1       233674  233675  +
The location of the variants is after the alingment

For NarBan, the last alingment listed in TEST_aln_only_CM068756.1.txt is:
    R       CM072956.1      13129647        13158556
The last few variants were:
    V       CM072956.1      13055212        13055213        106     255     t       c       CM072942.1      14785447        14785448        +
    V       CM072956.1      13055214        13055215        106     255     t       a       CM072942.1      14785449        14785450        +
    V       CM072956.1      13055216        13055217        106     255     t       g       CM072942.1      14785451        14785452        +
    V       CM072956.1      13055235        13055236        106     255     t       a       CM072942.1      14785470        14785471        +
These do not go beyond the end of the alignment

Preparing config.yml for GCA_035084215.1 Heptranchias perlo (sharpnose sevengill shark)
    CLADE : "sharks"
    SPEC_NAME : "HepPer_Test"
    REF_NAME : "GCA_035084215.1"
    ALT_NAME : "GCA_035084135.1"
    TODAY_DATE : "20250410"
    CHROM_START_CHR : "CM"
    CHROM_LIST_FILE : "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/sharks/chrom_lists/GCA_035084215.1_chroms.txt"
    WINDOW_INTERVAL : 500000
    WINDOW_LENGTH : 1000000
    NUM_AUT_CHROMOSOMES : 45
    Chromosomes : 
        - CM068637.1
        - CM068638.1
        - CM068639.1
        - CM068640.1
        - CM068641.1
        - CM068642.1
        - CM068643.1
        - CM068644.1
        - CM068645.1
        - CM068646.1
        - CM068647.1
        - CM068648.1
        - CM068649.1
        - CM068650.1
        - CM068651.1
        - CM068652.1
        - CM068653.1
        - CM068654.1
        - CM068655.1
        - CM068656.1
        - CM068657.1
        - CM068658.1
        - CM068659.1
        - CM068660.1
        - CM068661.1
        - CM068662.1
        - CM068663.1
        - CM068664.1
        - CM068665.1
        - CM068666.1
        - CM068667.1
        - CM068668.1
        - CM068669.1
        - CM068670.1
        - CM068671.1
        - CM068672.1
        - CM068673.1
        - CM068674.1
        - CM068675.1
        - CM068676.1
        - CM068677.1
        - CM068678.1
        - CM068679.1
        - CM068680.1
        - CM068681.1
        - CM068682.1
Running snakemake just to get FASTGA generated .1aln and .paf files for HepPer
Finished
Have updated FILTER_PAF_CHR_ONLY shell rule to:
    awk '{{ if ($6 !~ /^JA/) {{print $0}}}}' {input} > {output}
NOTE: Found that it included an extra chromosome from the alternate haplotype and not the primary
      Need to update FILTER_PAF_CHR_ONLY so that it only looks at alingments to reference

Run paftools and see what happens
    k8 ../../paftools.js call HepPer_Test_FASTGA_Filtered.srt.paf > TEST_aln_var.txt
Started running at 12:41PM
Canceled at 1:58PM after getting stuck on the first chromosome
    awk '$3 == 158632220' TEST_aln_var.txt
    V       CM068637.1      158632220       158632221       115     255     c       a       JAYJKZ010000001.1       39520909        39520910        +
    V       CM068637.1      158632220       158632220       115     255     -       t       JAYJKZ010000001.1       39531671        39531672        +
    V       CM068637.1      158632220       158632221       115     255     g       a       JAYJKZ010000001.1       39532293        39532294        +
    V       CM068637.1      158632220       158632221       115     255     t       a       JAYJKZ010000001.1       39533535        39533536        +
    V       CM068637.1      158632220       158632221       115     255     t       a       JAYJKZ010000001.1       39542846        39542847        +
    V       CM068637.1      158632220       158632221       115     255     a       t       JAYJKZ010000001.1       39434408        39434409        +
    V       CM068637.1      158632220       158632221       115     255     c       t       JAYJKZ010000001.1       39437504        39437505        +
Multiple variants from same alignment getting placed at the exact same spotted

Got help from Chenxi
He helped me solve the problem, which is two-fold:
    -- many repeats of alignments close to each other and just slightly shifted off from each other
    -- Very high variance in these repeats, leading to an exceedingly high number of variants for regions of the chromosome
Turns out code was working correctly the whole time, there were just some odd alignments
Solution is to filter out alignments with exceedingly high variance:
Command that Chenxi ran
    "k8 ../../paftools.js call <(awk -v OFS="\t" '{if(substr($13,6)<=0.1) print}' HydCol_Test2_FASTGA.srt.paf) >HydCol_Test2_FASTGA.srt_var.txt" >script.sh

Testing this solution with HydCol_Test2

k8 ../../paftools.js call <(awk -v OFS="\t" '{if(substr($13,6)<=0.1) print}' HydCol_Test2_FASTGA.srt.paf) > 20250410_HydCol_All_Chr_aln_var.txt
Started at 6:44PM
It finished at 6:49PM!!! It worked!!! (Hallelujah)


Added rule to Snakefile FILTER_PAF_VARIANCE to filter out alignments with high variance
Updated config.yml file for HydCol, and going into HydCol directory
    SPEC_NAME : "HydCol"
    TODAY_DATE : "20250331" ## NOTE -- DATE IS OLD...I forgot to update it before submitting. Will change retroactively
Running snakemake on this, first with only rules to generate .paf file

Paf file succesfully generated! Next uncommented commands to generate the var and aln files for the whole genome
submitted

It worked!!! And it gave a similar number of variants for the problem chromosome as from mm2
39572 variants for Chr CM068756.1

Generated everything but paf alingment and BAM coverage


#### UPDATE ####
20250411 (April 11th, 2025)

Updated config.yml for HepPer
     SPEC_NAME : "HepPer"
     TODAY_DATE : "20250411"
Otherwise the config is same as above

Successfully ran until GET_WHOLE_VAR -- Same problem as before was encountered, only making it to the second chromosome before issues arose
Canceled job after GET_WHOLE_VAR was running for ~30minutes without finishing
Confirmed that it was the same problem as before:
3039076 variants for CM068637.1 from FASTGA compared to 545765 from mm2

To fix this I will need a stricter cutoff for dv, less than 0.1 which is what Chenxi and I set yesterday
will test with CM068637.1.paf
    k8 ../../paftools.js call CM068637.1.paf | awk -v OFS="\t" '{{if(substr($13,6)<=0.8) print}}' > CM068637.1_var.txt
Submitted at 11:09AM
Canceled job at 11:29AM after running for 20 minutes
Had 767024 variants -- still encountering the problem

Realized I was off by 1 sig fig in the filtering (should've been 0.08 not 0.8)
Fixed and resubmitted
    k8 ../../paftools.js call CM068637.1.paf | awk -v OFS="\t" '{{if(substr($13,6)<=0.08) print}}' > CM068637.1_var.txt
Submitted at 11:32AM
Canceld at 11:50AM with hte same problem

We need more stringent criteria
    rm -f -r CM068637.1_var.txt
    k8 ../../paftools.js call CM068637.1.paf | awk -v OFS="\t" '{{if(substr($13,6)<=0.07) print}}' > CM068637.1_var.txt
Started running at 11:53AM

In the meeantime, updated the PAF_DOTPLOT rule to use ALNplot instead of pafr

Canceled job at 12:12PM
    (base) [ag2427@login-p-1 HepPer]$ awk -v OFS="\t" '{{if(substr($13,6)<=0.06) print}}' CM068637.1.paf | wc -l
    31407
    (base) [ag2427@login-p-1 HepPer]$ awk -v OFS="\t" '{{if(substr($13,6)<=0.1) print}}' CM068637.1.paf | wc -l
    63824
Reducing it to 0.06 would half the number of alingments used in comparison to the initial criteria
Will try this, and if it doesn't work then possibly reach out to Chenxi and Richard
    rm -f -r CM068637.1_var.txt
    k8 ../../paftools.js call CM068637.1.paf | awk -v OFS="\t" '{{if(substr($13,6)<=0.06) print}}' > CM068637.1_var.txt
Started at ~12:15PM
Canceled at 12:31PM after the SAME PROBLEM occurred
    awk -v OFS="\t" '{{if(substr($13,6)<=0.06 && $1 == "JAYJKZ010000002.1") print}}' CM068637.1.paf > CM068637.1_0.06_fltr.paf
    ALNplot -p CM068637.1_0.06_fltr.paf
This one scaffold in the filter above is responsible for the majority of the alingment to the chromosome
alingment shows that there is one problem region with all the repeats which are messy, the rest is clean

Need to see if this is a species specific problem
Updated config.yml for NarBan criteria
    CLADE : "sharks"
    SPEC_NAME : "NarBan"
    REF_NAME : "GCA_036971175.1"
    ALT_NAME : "GCA_036971445.1"
    TODAY_DATE : "20250411"
    CHROM_START_CHR : "CM"
    CHROM_LIST_FILE : "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/sharks/chrom_lists/GCA_036971175.1_chroms.txt"
    WINDOW_INTERVAL : 500000
    WINDOW_LENGTH : 1000000
    NUM_AUT_CHROMOSOMES : 14
    Chromosomes : 
    - CM072956.1
    - CM072957.1
    - CM072958.1
    - CM072959.1
    - CM072960.1
    - CM072961.1
    - CM072962.1
    - CM072963.1
    - CM072964.1
    - CM072965.1
    - CM072966.1
    - CM072967.1
    - CM072968.1
    - CM072969.1
Submitted snakemake at 12:44PM

While this is running and I am troubleshooting this, I can start generating the paf files for the other species
Config for HemOce
    CLADE : "sharks"
    SPEC_NAME : "HemOce"
    REF_NAME : "GCF_020745735.1"
    ALT_NAME : "GCA_020745765.1"
    TODAY_DATE : "20250411"
    CHROM_START_CHR : "NC"
    CHROM_LIST_FILE : "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/sharks/chrom_lists/GCF_020745735.1_chroms.txt"
    WINDOW_INTERVAL : 500000
    WINDOW_LENGTH : 1000000
    NUM_AUT_CHROMOSOMES : 52
    Chromosomes : 
    - NC_083401.1
    - NC_083402.1
    - NC_083403.1
    - NC_083404.1
    - NC_083405.1
    - NC_083406.1
    - NC_083407.1
    - NC_083408.1
    - NC_083409.1
    - NC_083410.1
    - NC_083411.1
    - NC_083412.1
    - NC_083413.1
    - NC_083414.1
    - NC_083415.1
    - NC_083416.1
    - NC_083417.1
    - NC_083418.1
    - NC_083419.1
    - NC_083420.1
    - NC_083421.1
    - NC_083422.1
    - NC_083423.1
    - NC_083424.1
    - NC_083425.1
    - NC_083426.1
    - NC_083427.1
    - NC_083428.1
    - NC_083429.1
    - NC_083430.1
    - NC_083431.1
    - NC_083432.1
    - NC_083433.1
    - NC_083434.1
    - NC_083435.1
    - NC_083436.1
    - NC_083437.1
    - NC_083438.1
    - NC_083439.1
    - NC_083440.1
    - NC_083441.1
    - NC_083442.1
    - NC_083443.1
    - NC_083444.1
    - NC_083445.1
    - NC_083446.1
    - NC_083447.1
    - NC_083448.1
    - NC_083449.1
    - NC_083450.1
    - NC_083451.1
    - NC_083452.1
    - NC_083453.1
    - NC_083454.1

Config for HetFra
    CLADE : "sharks"
    SPEC_NAME : "HetFra"
    REF_NAME : "GCA_036365495.1"
    ALT_NAME : "GCA_036365525.1"
    TODAY_DATE : "20250411"
    CHROM_START_CHR : "CM"
    CHROM_LIST_FILE : "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/sharks/chrom_lists/GCA_036365495.1_chroms.txt"
    WINDOW_INTERVAL : 500000
    WINDOW_LENGTH : 1000000
    NUM_AUT_CHROMOSOMES : 50
    Chromosomes : 
    - CM070962.1
    - CM070963.1
    - CM070964.1
    - CM070965.1
    - CM070966.1
    - CM070967.1
    - CM070968.1
    - CM070969.1
    - CM070970.1
    - CM070971.1
    - CM070972.1
    - CM070973.1
    - CM070974.1
    - CM070975.1
    - CM070976.1
    - CM070977.1
    - CM070978.1
    - CM070979.1
    - CM070980.1
    - CM070981.1
    - CM070982.1
    - CM070983.1
    - CM070984.1
    - CM070985.1
    - CM070986.1
    - CM070987.1
    - CM070988.1
    - CM070989.1
    - CM070990.1
    - CM070991.1
    - CM070992.1
    - CM070993.1
    - CM070994.1
    - CM070995.1
    - CM070996.1
    - CM070997.1
    - CM070998.1
    - CM070999.1
    - CM071000.1
    - CM071001.1
    - CM071002.1
    - CM071003.1
    - CM071004.1
    - CM071005.1
    - CM071006.1
    - CM071007.1
    - CM071008.1
    - CM071009.1
    - CM071010.1
    - CM071011.1
    - CM071012.1

Config for MobBir
    CLADE : "sharks"
    SPEC_NAME : "MobBir"
    REF_NAME : "GCA_030028105.1"
    ALT_NAME : "GCA_030035685.1"
    TODAY_DATE : "20250411"
    CHROM_START_CHR : "CM"
    CHROM_LIST_FILE : "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/sharks/chrom_lists/GCA_030028105.1_chroms.txt"
    WINDOW_INTERVAL : 500000
    WINDOW_LENGTH : 1000000
    NUM_AUT_CHROMOSOMES : 32
    Chromosomes : 
    - CM057556.1
    - CM057557.1
    - CM057558.1
    - CM057559.1
    - CM057560.1
    - CM057561.1
    - CM057562.1
    - CM057563.1
    - CM057564.1
    - CM057565.1
    - CM057566.1
    - CM057567.1
    - CM057568.1
    - CM057569.1
    - CM057570.1
    - CM057571.1
    - CM057572.1
    - CM057573.1
    - CM057574.1
    - CM057575.1
    - CM057576.1
    - CM057577.1
    - CM057578.1
    - CM057579.1
    - CM057580.1
    - CM057581.1
    - CM057582.1
    - CM057583.1
    - CM057584.1
    - CM057585.1
    - CM057586.1
    - CM057587.1
    - CM057588.1

Config for HypSab
    CLADE : "sharks"
    SPEC_NAME : "HypSab"
    REF_NAME : "GCF_030144855.1"
    ALT_NAME : "GCA_030144785.1"
    TODAY_DATE : "20250411"
    CHROM_START_CHR : "NC"
    CHROM_LIST_FILE : "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/sharks/chrom_lists/GCF_030144855.1_chroms.txt"
    WINDOW_INTERVAL : 500000
    WINDOW_LENGTH : 1000000
    NUM_AUT_CHROMOSOMES : 32
    Chromosomes : 
    - NC_082706.1
    - NC_082707.1
    - NC_082708.1
    - NC_082709.1
    - NC_082710.1
    - NC_082711.1
    - NC_082712.1
    - NC_082713.1
    - NC_082714.1
    - NC_082715.1
    - NC_082716.1
    - NC_082717.1
    - NC_082718.1
    - NC_082719.1
    - NC_082720.1
    - NC_082721.1
    - NC_082722.1
    - NC_082723.1
    - NC_082724.1
    - NC_082725.1
    - NC_082726.1
    - NC_082727.1
    - NC_082728.1
    - NC_082729.1
    - NC_082730.1
    - NC_082731.1
    - NC_082732.1
    - NC_082733.1
    - NC_082734.1
    - NC_082735.1
    - NC_082736.1
    - NC_082737.1
    - NC_082738.1
    - NC_082739.1
    - NC_082740.1


Flipping the order of the filtering and sorting rules in snakemake to make it more efficient
Successfully generated paf file for NarBan
Running snakemake to generate Var_only and Aln_only txt files
Started running at 5:43PM
Canceled at 6:27PM when it was still on the first chromosome 
Checked region where variants were being generated -- dv value was 0.08-0.09 range
Will make criteria more strict -- .07 cutoff
    awk -v OFS="\t" '{{if(substr($13,6)<=0.07 && $6 == "CM072956.1") print}}' NarBan_FASTGA.fltr.srt.paf > CM072956.1.paf
Now will run this to see if the same error is occuring:
    k8 ../../paftools.js call CM072956.1.paf > test.txt
Canceled this one, because it appeared to be running in a circle -- will wait for results of NarBan snakefile with updated criteria to see

Realized I was running the filtering at the incorrect step this morning on HepPer -- will attempt to re-run it
    awk -v OFS="\t" '{{if(substr($13,6)<=0.07 && $6 == "CM068637.1") print}}' HepPer_FASTGA.fltr.srt.paf > CM068637.1.paf
    k8 ../../paftools.js call CM068637.1.paf > test.txt
Started running at 7:29PM
Finished at 7:32PM
It worked?!

Given that this worked, I updated the shell in the FILTER_PAF_VARIANCE rule:
    awk -v OFS="\t" '{{if(substr($13,6)<=0.0.07) print}}' {input} > {output}
Removed old files
    (snakemake) [ag2427@login-p-3 HepPer]$ rm -f -r HepPer_FASTGA.srt.paf
    (snakemake) [ag2427@login-p-3 HepPer]$ rm -f -r HepPer_FASTGA.fltr.srt.paf
    (snakemake) [ag2427@login-p-3 HepPer]$ rm -f -r HepPer_FASTGA.fltr.srt.pdf
Resubmitted on snakemake to get updated files with stricter filtering and new variants


#### UPDATE ####
20250412 (April 12th, 2025)

NarBan stopped at SORT_PAF due to running out of memory
Realized there was a syntax error in FILTER_PAF_VARIANCE and fixed
    awk -v OFS="\t" '{{if(substr($13,6)<=0.07) print}}' {input} > {output}
Successfully generated paf files -- variants still had same problem as before -- got stuck on the first chromosome


#### UPDATE ####
20250414 (April 14th, 2025)

Testing NarBan first chromosome to see if lowering it to 0.06 will fix the problem
    awk -v OFS="\t" '{{if(substr($13,6)<=0.05 && $6 == "CM072956.1") print}}' NarBan_FASTGA.fltr.srt.paf > CM072956.1.paf
    k8 ../../paftools.js call CM072956.1.paf 
Running into the same problem -- even with more stringent criteria. 
I'm thinking I need to filter for removing duplicates as well or by another filter than just dv
Approximately where it starts:
    V       CM072956.1      69270141        69270142        2       255     c       a       CM072942.1      70907869        70907870        +
    V       CM072956.1      69270176        69270177        2       255     a       t       CM072942.1      70907904        70907905        +
    V       CM072956.1      69270183        69270184        2       255     c       a       CM072942.1      70907911        70907912        +
    V       CM072956.1      69270194        69270195        2       255     g       c       CM072942.1      70907922        70907923        +
    V       CM072956.1      69270218        69270219        3       255     a       g       CM072942.1      70907946        70907947        +
Need to find where this is in the paf file
    awk '$8 <= 69270141 && $9 >= 69270141 && $1 == "CM072942.1"' CM072956.1.paf | less -SN
Two lines were shown:
      1 CM072942.1      494489450       70540089        71137278        +       CM072956.1      495097020       68902425        69499504        587098  597662  255     dv:f:0.0169     df:i:10564      cs:Z::4-a:2*ga:6+g:1+gg:10*ga:1+ggg:16+gaa:2*cg:8+a:7+g:6*ga:>
      2 CM072942.1      494489450       70907398        71137276        +       CM072956.1      495097020       69269964        69499795        220126  230216  255     dv:f:0.0424     df:i:10090  
Possibly see if removing these two lines fixes the problem
    awk '$13 == "dv:f:0.0169" && $14 == "df:i:10564"' CM072956.1.paf | less -SN
This indeed only returns the first line --> remove these and re-run
    awk '$13 != "dv:f:0.0169" && $14 != "df:i:10564"' CM072956.1.paf > test.paf
    awk '$13 == "dv:f:0.0169" && $14 == "df:i:10564"' test.paf | less -SN 
Nothing returned in second awk --> removal worked
    k8 ../../paftools.js call test.paf
This did not fix the problem, but it did change the query depth from 130 to 129
    V       CM072956.1      69392316        69392317        129     255     c       t       CM072942.1      71028617        71028618        +
    V       CM072956.1      69392319        69392320        129     255     c       a       CM072942.1      71028620        71028621        +
    V       CM072956.1      69392327        69392328        129     255     c       t       CM072942.1      71028628        71028629        +
    V       CM072956.1      69392357        69392358        129     255     a       c       CM072942.1      71028658        71028659        +
    V       CM072956.1      69392383        69392384        129     255     a       t       CM072942.1      71028684        71028685        +
Will try looking into these variants
    awk '$8 <= 69392383 && $9 >= 69392383 && $1 == "CM072942.1"' test.paf | less -SN
This returns 129 alingments, the same number as the number of query depth for the variants
Would removing the repeat alignments fix the problem?
    awk '$8 <= 69392383 && $9 >= 69392383 && $1 == "CM072942.1"' CM072956.1.paf | less -SN
130 alingments are returned here
    awk '$8 <= 69392383 && $9 >= 69392383 && $1 == "CM072942.1"' CM072956.1.paf > alns_130_depth.paf
Will remove all alingments except the first one, which is the longest, and see if this fixes the problem
    awk '$8 <= 69392383 && $9 >= 69392383 && $1 == "CM072942.1"' test.paf > alns_to_remove.paf
    awk '{if (f==1) { r[$0] } else if (! ($0 in r)) { print $0 } } ' f=1 alns_to_remove.paf f=2 CM072956.1.paf
    awk '$8 <= 69392383 && $9 >= 69392383 && $1 == "CM072942.1"' CM072956.1_alns_removed.paf | less -SN
The last awk only shows the first line from alns_130_depth.paf -- it worked
    k8 ../../paftools.js call CM072956.1_alns_removed.paf

It worked! It made it past the region where it was spinning endlessly
I'm going to try this same fix on HepPer and see if it works there as well
Problem chromosome in HepPer is CM068637.1
    awk -v OFS="\t" '{{if(substr($13,6)<=0.07 && $6 == "CM068637.1") print}}' HepPer_FASTGA.paf > CM068637.1.paf
    sort -k6,6V -k8,8n CM068637.1.paf > CM068637.1.srt.paf
    k8 ../../paftools.js call CM068637.1.srt.paf
Canceled when it got stuck around base 54905907
    awk '$8 <= 54905907 && $9 >= 54905907' CM068637.1.srt.paf | less -SN
Again, matching query number of where it got stuck, there are 154 alingments at this point
One has a significantly higher number of matching bases in comparison to the others, retain it and see if it fixes the problem
    awk '$8 <= 54905907 && $9 >= 54905907' CM068637.1.srt.paf > alns_test.paf
    tail -n +2 alns_test.paf > alns_to_remove.paf
    awk '{if (f==1) { r[$0] } else if (! ($0 in r)) { print $0 } } ' f=1 alns_to_remove.paf f=2 CM068637.1.srt.paf
    awk '{if (f==1) { r[$0] } else if (! ($0 in r)) { print $0 } } ' f=1 alns_to_remove.paf f=2 CM068637.1.srt.paf > CM068637.1_aln_removed.paf
    k8 ../../paftools.js CM068637.1_aln_removed.paf > test_aln_var_results.txt
This again appears to have worked!

How could I make this quantitative?
Maybe do it for all query depths >= 100
Steps:
    -- Load in sorted and filtered PAF file
    -- Use paftools.js call function
    -- When query depth >= 100
        -- stop function and save the last variant line of the file
        -- Find all alns in the paf file which match that chromosome and which have start < var and end > var
        -- Find aln out of that list with the highest number of base matches
        -- Keep that aln, delete the others out of the paf file
    -- Repeat paftools.js call function but this time on the updated paf file
    -- Finish when paftools.js call function finishes without query depth > 100

Created 20250414_filter_repeat_alns.sh to try and do this

Heard from Richard regarding issue -- apparently ALNchain should fix issue
Testing
    ALNchain 20250221_HepPer_ALN.1aln
created output 20250221_HepPer_ALN.chain.1aln
    ALNtoPAF -s -T8 20250221_HepPer_ALN.chain.1aln > chain_test.paf
It worked!
Now will test if same number of repeats is found
    awk '$1 == "CM068637.1" && $3 <= 54905907 && $4 >= 54905907' chain_test.paf | less -SN
Only 3 are found for this position
Test to see if paftools gets stuck
    k8 ../../paftools.js call chain_test.paf
The query and targets swapped places in the paf file, and as a result also in the output file, and was causing issues
Will re-run fastga with the reference first and alt second, and then do ALNchain and ALNtoPAF to see if this fixes the issue
    FastGA -v -P. -T8 -1:20240414_test 20250221_HepPer_GCA_035084215.1.gix 20250221_HepPer_GCA_035084135.1.gix
    ALNchain 20240414_test.1aln
    ALNtoPAF -s -T8 20240414_test.chain.1aln > 20240414_test.paf


#### UPDATE #### 
20250415 (April 15th 2025)

First going to test running paftools on this to see if running ALNchain fixes the issue
    k8 ../../../paftools.js call 20240414_test.paf > test_out_var_aln.txt
It appears to work, but it is also running alingments on scaffolds -- need to remove them first
    awk '{ if ($1 ~ /^JA/) print}' 20240414_test.paf | wc -l
    164208
    awk '{ if ($1 ~ /^CM/) print}' 20240414_test.paf | wc -l
    100448
Will filter out alingments to scaffolds in query
    awk '{{ if ($1 !~ /^JA/) {{print $0}}}}' 20240414_test.paf > 20240414_test.fltr.paf
Re-run paftools
    k8 ../../../paftools.js call 20240414_test.fltr.paf > test_out_var_aln.txt
Same issue as before -- realized this is because the paf folder still was improperly formatted for paftools
    k8 ../../../paftools.js call -f ../../../../241117.UCSC-hubs-VGP-alignment/alignment/reference/sharks/GCA_035084215.1.fa.gz 20240414_test.fltr.paf


#### UPDATE ####
20250416 (April 16th, 2025)

First, I realized I ran FASTGA on outdated files from February -- will clean out temp directory in HepPer and then re-run to see if I can fix the formatting issue
    rm -f -r 20250221*
    rm -f -r 20240414*
    FastGA -v -P. -T8 -1:20250416_FASTGA_Ref_first HepPer_GCA_035084215.1.gix HepPer_GCA_035084135.1.gix
    FastGA -v -P. -T8 -1:20250416_FASTGA_Alt_first HepPer_GCA_035084135.1.gix HepPer_GCA_035084215.1.gix
Running FastGA with both to see if that makes a difference in order, and which once I should do
    ALNtoPAF -s -T8 20250416_FASTGA_Ref_first.1aln > 20250416_FASTGA_Ref_first.paf
    ALNtoPAF -s -T8 20250416_FASTGA_Alt_first.1aln > 20250416_FASTGA_Alt_first.paf
Started running these at 12:40PM
Finished at 12:50 and 12:52
Proper formatting found for 20250416_FASTGA_Alt_first.paf
Will now try running ALNchain on 20250416_FASTGA_Alt_first.1aln and then converting to paf and seeing what happens
    ALNchain -v 20250416_FASTGA_Alt_first.1aln
    Segmentation fault
This generates part of a file, but with the error: Segmentation fault
    ALNchain -v 20250416_FASTGA_Ref_first.1aln
This runs normally
Will attempt to convert 20250416_FASTGA_Alt_first.chain.1aln to paf and see if formatting is fixed
    ALNtoPAF -s -T8 20250416_FASTGA_Alt_first.chain.1aln > 20250416_FASTGA_Alt_first.chain.paf
Cannot run - I get this error
    FATAL ERROR: ONE file error: can't seek to start of footer
Do not get this error when running with 20250416_FASTGA_Ref_first.chain.1aln

k8 ../../../paftools.js call test.paf > test_wrong_orientation.txt
Confimed this didn't work

Emailed Chenxi -- He confirmed that Alt then Ref are listed for FASTGA
Need to see if the issue is occurring with other species
Went to sharks/NarBan
    ALNchain NarBan_ALN.1aln
    ALNchain: retained 98382 alignments in 27795 chains
This worked -- Generating NarBan_ALN.chain.1aln
    ALNtoPAF -s -T8 NarBan_ALN.chain.1aln > 20250416_NarBan_ALN.chain.paf
    k8 ../../paftools.js call 20250416_NarBan_ALN.chain.paf
This appears to work fine?

Set config.yml for NarBan
Running with ALNchain rule to see if I can get variant and alignment files
It worked! Was able to generate Aln and Var files without any issues for all chromosomes
Now running Snakemake for other results -- 
    Running into memory issues for my ROH calculation script -- will work on that later
    Want to test if the variant generation works for other elasmobranchs

Set up config.yml for HemOce
Running snakemake only to get paf file, and variant and aln files
Submitted
It worked! Note, it also generated alingments and variants for unplaced scaffolds
I don't think this should be an issue for downstream analyses, as those will only look for relevant chromosomes based on list in config.yml

Set up config.yml for HetFra
Running snakemake only to get paf file, and variant and aln files
Submitted
Job stopped due to memory consumption on GIXmake for HetFra_GCA_036365495.1.gix


#### UPDATE ####
20250417 (April 17th, 2025)

Downloaded plugin to run snakemake on slurm
    pip install snakemake-executor-plugin-slurm

Added new section to config file describing slurm resources:
    default-resources:
        slurm_account: "DURBIN-SL2-CPU"
        slurm_partition: "cclake"
        mem_mb: 100000
        runtime: "2h"

Will attempt to run snakemake on slurm for HetFra to see if this can fix memory consumption problem
    snakemake --executor slurm --configfile config.yml
    Error: maximum number of parallel jobs/used nodes has to be specified for remote execution (use --jobs N with N being a number >= 1)
Try again
    snakemake --executor slurm --workflow-profile profiles --jobs 10 --configfile config.yml
profiles/config.yaml == file that describes slurm running parameters (e.g. mem, time per job)
config.yml == file that describes config parameters for the species being analysed 


#### UPDATE ####
20250418 (April 18th, 2025)

Check on update for snakemake -- no jobs currently running.
Generated GIX files, but not .1aln file
I will up the time for running FASTGA rule to 12hrs and resubmit


#### UPDATE ####
20250422 (April 22nd, 2025)

Recompiled snakemake to get fix for ALNchain which Chenxi implemented
    git clone https://github.com/thegenemyers/FASTGA.git
    cd src/FASTGA
    make && make install

Checked on results of running snakemake on HetFra
.1aln file was generated, but not any of the var or aln files
Resubmitted

Will also re-run ALNchain on HepPer to see if ALNchain fix worked for issue
    ALNchain HepPer_ALN.1aln
It finished, and worked!!!!
Generated HepPer_ALN.chain.1aln
Once snakemake finishes running on HetFra, will have it run on HepPer

Canceled job running on HetFra after it had been running for 30 minutes
The aln and var files hadn't updated in 23 minutes
It got stuck on the X sex chromosome CM071012.1

Updated config.yml for HepPer, and will run it while trying to find a fix for HetFra sex chromosome
Generated var and aln files, both for whole genome and per-chr
It worked!
Running for ROH and Het calculations


#### UPDATE ####
20250423 (April 23rd, 2025)

Snakemake run for HepPer stopped yesterday as a result of a bad wifi connection
Re-ran

Generated simplified DAG of snakemake rules
snakemake --executor slurm --workflow-profile profiles --jobs 10 --configfile config.yml --rulegraph | dot -Tpdf > dag.pdf


#### UPDATE ####
20250424 (April 24th, 2025)

Some of the previous job ran -- generated ROH plot, heterozygosity values, and FROH for all chromosomes
Issue is on the WHOLE_FROH rule 

Re-ran locally and same issue resulted 
Issue in line 389:
    head: cannot open 'sharks/HepPer/temp/20250422__Var_Only.txt' for reading: No such file or directory

Checked sharks/chrom_lists/GCA_035084215.1_chroms.txt
Problem was that there was an empty line at the end of the file after the chromosomes, which was being read in shell.
Removed extra line and re-ran
This fixed the issue

New issue arose with CALC_HET_PER_CHR rule -- Job failing out on slurm for 6 chromosomes
Potential issue -- the Chromosome length file contains one extra chromosome than is offically listed
    CM068729.1      19484
Need to fix this 
Original command
        zcat < {input} | awk '$0 ~ ">" {{if (NR > 1) {{print c;}} c=0;printf substr($0,2,100) "\\t"; }} $0 !~ ">" {{c+=length($0);}} END {{ print c; }}' | grep "{params.CHROM_START_CHR}" > {output}
Altered in CHROM_LENGTH_CALC Rule:
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
Removed file and re-ran
    rm -f -r Reference_HepPer_Chroms_Lengths.txt
It worked to successfully generate a reference file

Rest of job is working with the exception of an error when plotting heterozygosity
    MissingOutputException in rule PLOT_HET_PER_CHR in file /rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/Snakefile, line 481:
    Job 293  completed successfully, but some output files are missing. Missing files after 5 seconds. This might be due to filesystem latency. If that is the case, consider to increase the wait time with --latency-wait:
Problem is not generating sharks/HepPer/20250422_HepPer_CM068682.1_Het_Map.svg
This problem persists even after increasing latency wait to 60 seconds

Figured out issue -- only plots for all autosomal chromosomes but snakemake is expecting svgs for all chromosomes, including sex chromosomes
Altered 20250123_Plot_het_per_chr.R to plot for all chromosomes
Added NUM_ALL_CHR to config.yml
Got it to work!

Will now set up for HetFra
    NUM_ALL_CHR: 51
Submitted to cluster

Downloaded VGP freeze (will likely have to re-download) to work on introductory figure
Created script 20250424_Taxon_Figure.R for this purpose


#### UPDATE ####
20250425 (April 25th, 2024)

Submission of Snakefile with HetFra ran into issues generating the var and aln files at the sex chromosome
Thinking I should alter the rules so that it only looks for autosomal chromosomes to recorded
Brought back FILTER_PAF_CHR_ONLY rule and altered it so that I only have the lines in the paf file for the relevant chromosomes
    awk 'BEGIN { while (getline < {input.ALL_CHROMS}) list[$0] } $6 in list' {input.PAF} > {output}
Tested this and it worked, including on the sex chromosome
Job failed on rule COMPILE_ROH, due to this error:
    MissingOutputException in rule COMPILE_ROH in file /rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/Snakefile, line 323:
    Job 167  completed successfully, but some output files are missing. Missing files after 5 seconds. This might be due to filesystem latency. If that is the case, consider to increase the wait time with --latency-wait:    
When I look for ROH_Results.txt files, it appears that some are missing, I only find 42
Re-ran, and got to 48 out of 51 needed.
Commented out COMPILE_ROH rule and everything below -- wondering if it is running before ROH are finished getting made
Re-ran ROH_CALC rule, increased --latency-wait to 60 seconds
It appears to have worked -- The issues were due to slurm issues
Finished running, and it worked!

Set up config.yml for MobBir
    TODAY_DATE : "20250425"
    NUM_ALL_CHR: 33
Changed slurm partition in config.yaml to icelake given node failures on cclake
Submitted

Will try creating an additional config file to see if I can run multiple species at the same time
Created NarBan_config.yml
Running by changing config file input
    snakemake --executor slurm --workflow-profile profiles --jobs 10 --configfile NarBan_config.yml
So far it appears to be working! Submitted jobs to SLURM

Job succesfully finished for NarBan

Set up config file for HypSab
Created HypSab_config.yml
Ran by changing config file input
    snakemake --executor slurm --workflow-profile profiles --jobs 10 --configfile HypSab_config.yml


#### UPDATE ####
20250428 (April 28th, 2025)

.1aln files generated for both HypSab and MobBir
Going to run Snakefile to generate other output files from these
    snakemake --executor slurm --workflow-profile profiles --jobs 10 --configfile HypSab_config.yml
    snakemake --executor slurm --workflow-profile profiles --jobs 10 --configfile config.yml

Run Snakefile for HemOce:
    snakemake --executor slurm --workflow-profile profiles --jobs 10 --configfile HemOce_config.yml


#### UPDATE ####
20250430 (April 30th, 2025)

Run Snakefile for HemOce:
    snakemake --executor slurm --workflow-profile profiles --jobs 10 --configfile HemOce_config.yml
Jobs didn't run due to HPC being down for maintenance


#### UPDATE ####
20250501 (May 1st, 2025)

Need to compare ROH locations and low heterozygosity locations in HepPer to see if they match
Comparing the ROH and Whole Genome Het Map, the chromosomes which are almost entirely ROH do match the ones with low heterozygosity
Need to more deeply review heterozygosity code for potential causes of problems in Het w/o ROH being lower than total Het --
    Find out what is writing results in columns 5 and 8 of the HepPer_Het_Compiled.tsv file


#### UPDATE ####
20250502 (May 2nd, 2025)

When calculating heterozygosity per chromosome, the Het_excl_ROH is the same for all windows, which doesn't make sense
I think I found the issue -- in 20250325_find_het_per_chr_V3.py, for loop in lines 122 t0 149, the heterozygosity was populating the whole Het_excl_ROH column instead of a single cell due to a syntax error
Fixed error in an interactive sesssion
Removed all heterozygosity files from HepPer and created HepPer_config.yml
    CLADE : "sharks"
    SPEC_NAME : "HepPer"
    REF_NAME : "GCA_035084215.1"
    ALT_NAME : "GCA_035084135.1"
    TODAY_DATE : "20250422"
    CHROM_START_CHR : "CM"
    CHROM_LIST_FILE : "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/sharks/chrom_lists/GCA_035084215.1_chroms.txt"
    WINDOW_INTERVAL : 500000
    WINDOW_LENGTH : 1000000
    NUM_AUT_CHROMOSOMES : 45
    NUM_ALL_CHR: 46
    Chromosomes : 
        - CM068637.1
        - CM068638.1
        - CM068639.1
        - CM068640.1
        - CM068641.1
        - CM068642.1
        - CM068643.1
        - CM068644.1
        - CM068645.1
        - CM068646.1
        - CM068647.1
        - CM068648.1
        - CM068649.1
        - CM068650.1
        - CM068651.1
        - CM068652.1
        - CM068653.1
        - CM068654.1
        - CM068655.1
        - CM068656.1
        - CM068657.1
        - CM068658.1
        - CM068659.1
        - CM068660.1
        - CM068661.1
        - CM068662.1
        - CM068663.1
        - CM068664.1
        - CM068665.1
        - CM068666.1
        - CM068667.1
        - CM068668.1
        - CM068669.1
        - CM068670.1
        - CM068671.1
        - CM068672.1
        - CM068673.1
        - CM068674.1
        - CM068675.1
        - CM068676.1
        - CM068677.1
        - CM068678.1
        - CM068679.1
        - CM068680.1
        - CM068681.1
        - CM068682.1
Resubmitted to run locally
Realized I will have to do this on teh cluster once it is back up and running -- new config file means that files are being remade

Removed all heterozygosity files for NarBan and resubmitted since that already had a separate config file.
It ran succesfully, but Het_excl_ROH was significantly lower than it was before
There is another error
I realized the other error is that I am not accounting for the shorter length of the genome where ROH are present
E.g. for  Het_Per_Kb_excl_ROH I am diving the number of variants by the window size, not accounting for less of the window being counted due to ROH being present


#### UPDATE ####
20250505 (May 5th, 2025)

Problems leading to odd heterozygosity results"
    -- Not accounting for less of the window being counted due to ROH being present
    -- Varaints being counted within a run of homozygosity --> need to dig into the script to see what is going on

Going to run an interactive session with NarBan data to work through ROH calculation code to see if I can find something going screwy
Possibly just need to increase the penalty score
Run within snakemake environment in NarBan directory
chrom = "CM072956.1"
chrom_length = 495097020
ref_name = "GCA_036971175.1"
clade = "sharks"
Aln_file = open("temp/NarBan_Aln_Only.txt", 'r')
Var_file = open("temp/NarBan_Var_Only.txt", 'r')
spec_name = "NarBan"

I think I found issue #1
The base_aln_map from the aln_map function is not working properly,
The aln_sums are higher than the entire length of the chromosome
Updated the aln_map function to fix this.
    def aln_map(Alignment_list, chrom_length):
    base_sums = np.zeros(chrom_length, dtype=int) ## Create an array the length of the chromosome
    for start, end in ((int(line[2]), int(line[3])) for line in Alignment_list):
        base_sums[start:end] += 1 ## Calculate depth of alingments at each individual base
    base_gross_cov = np.where(base_sums > 0, 1, 0) ## Convert depth to binary of whether it is covered or not
    aln_sums = np.cumsum(base_gross_cov) ## Get alingment sum at each given base
    return np.vstack((np.arange(1, chrom_length + 1), base_sums, aln_sums.astype(int)))
Will now resubmit codes for HepPer to see what differences this makes
Removed all files except those for FASTGA and generating the paf files within HepPer
Resubmitted
    snakemake --executor slurm --workflow-profile profiles --jobs 10 --configfile HepPer_config.yml

As for the issue in 20250325_find_het_per_chr_V3.py, need to account for 'shorter' chromosome when ROH are removed
Will add commands to calc_het function to sum all variants and divide by length of chromosome
Do for both including and excluding ROH to see difference between them and to what I am currently doing
I will need to account for overlapping windows


#### UPDATE ####
20250506 (May 6th, 2025)

Rerunning of HepPer yesterday didn't finish, so I resubmitted it
Finished generating ROH and FROH results, but getting errors with heterozygosity
First Het file generated for CM068637.1 shows no variants when excl ROH after base 96500000, which is nearly half of the chromosome
But checking ROH file shows a singular ROH from 96329383 to 96451101

Something is definitely going wrong with the heterozygosity calculations
Running interactive session to see what is going wrong
Will use HepPer
    dat = "temp/HepPer_Var_Only.txt"
    chr_file = "Reference_HepPer_Chroms_Lengths.txt"
    clade = "sharks"
    species = "HepPer"
    current_window_length = 1000000
    current_window_interval = 500000
    num_aut_chr = 45
    roh_data = "HepPer_ROH_Results.csv"

mean_het_df = pd.DataFrame({
    'chr': all_chromosomes, 
    'mean_het': chr_mean_het, 
    'mean_het_excl_ROH': chr_mean_het_excl_ROH, 
})
    chr = 'CM068637.1'

window_starts = np.arange(0, (chr_length-current_window_length+1), current_window_interval)
window_ends = window_starts + current_window_length
window_sizes = np.repeat(current_window_length, len(window_starts))
single_chr_results = pd.DataFrame({
    'Start': window_starts, 
    'End': window_ends, 
    'Het': het, 
    'Het_excl_ROH': het_wo_roh, 
    'Window_Size': window_sizes, 
    'Window_Size_excl_ROH': window_sizes
})
Verified that sum of variants in each window is done correctly (line 110)
Possible that the issue is still coming back to ROH calculation
E.g. in HepPer chr CM068637.1, there is a single ROH detected from 
96329383 to 96451101, with length 121,718
supposedly 940 variants are found in this -- this seems quite high for something we're calling an ROH
Yet when I filter for variants in this range
    test = variant_positions[variant_positions > 96300000]
    test = test[test < 96500000]
There are only 89 variants found

I think I realized one of the issues -- not filtering for variants just in the window before doing counts
Fixed this by filtering for variants just in the window and updating Het_excl_ROH calculation
Updated single_chr_results df to include window size excluding ROH, which will be used to calculate Het_excl_ROH per kb
Removed heterozygosity result files from HepPer and resubmitted
Despite changes, the same issues persisted
Issue is with line 129:
    l = int(np.searchsorted(ROH_starts, start, side='right') - 1)

Appears that data freeze is finally done and the data is downloaded!
In VGP/250430.VGP-Phase1


#### UPDATE ####
20250507 (May 7th, 2025)

Getting closer to fixing heterozygosity code, but not 100% 
Het_excl_ROH is still lower than Het for everything on a handful of chromosomes in HepPer
    E.g. CM068641.1,0.16436131386861316,0.022619524390398278
I think it is still an issue of indexing with ROH, 
E.g. first het window in the chromosome:
    CM068641.1,0,1000000,4781.0,0.0,1000000,252294,4.781,0.0
ROH in this region:
    CM068641.1,252294,1584956,1332662

chr = "CM068641.1" 
It appears the issue with this chromosome is the heterozygous loci being found within the ROH
All 4781 variants are found in bases 304413 - 334159, which is within the ROH
Also 4781 variants wtihin a 29746 base range, which is much too high to be considered part of an ROH
The next variant isn't until base 1753267, which is over 1.4 million bases later
Will check ROH calculation script for this chromosome
chrom = "CM068641.1"
chrom_length = 137863001
Var_file = open("temp/HepPer_Var_Only.txt", 'r')
Aln_file = open("temp/HepPer_Aln_Only.txt", 'r')
spec_name = "HepPer"


../250430.VGP-Phase1/vgp.alignment.set.metaData.txt


#### UPDATE ####
20250408 (May 8th, 2025)

Modified Snakefile to add rules for doing MSMC analysis on species
Created 20250508_Plot_MSMC.R to use in Snakefile to plot MSMC results for each individual species
Rules for MSMC analysis are not completed -- Need to continue working on them.

For now I'm going to switch to troubleshooting ROH for chromosome on HepPer described above in lines 5572 to 5582
Running python interactive session
Coverage using aln_map function looks like it makes sense for bases bases 304413 - 334159
Wondering if the issue is due to the variants not being properly filtered for each chromosome
I updated the function for reading variants and re-ran it with the new function:
    def read_var_file(variant_file, chrom):
        var_lines = [line.split() for line in variant_file.readlines() if line.split()[1] == chrom]
        var_pos = [int(line[2]) for line in var_lines]
        return var_lines, var_pos
Results were different for chromosome CM068641.1, with the first ROH starting at 334158 instead of 252294
Given this, I will keep new function and re-run to see differences across genome
Removed het and ROH files
    rm -f -r *_Het_Map.svg
    rm -f -r 20250422_HepPer_Het_Whole_Genome_Map.svg
    rm -f -r HepPer_Het_Compiled.tsv
    rm -f -r *mean_heterozygosity.txt
    rm -f -r *_het.txt
    rm -f -r *_FROH.txt
    rm -f -r *ROH_Map.pdf
    rm -f -r *ROH_Results*
Resubmitted analysis
    snakemake --executor slurm --workflow-profile profiles --jobs 10 --configfile HepPer_config.yml


#### UPDATE ####
20250509 (May 9th, 2025)

ROH results were generated, but for chr CM068641.1 they were the same as they were before, not what I was getting in the interactive session
What is going on?
Removed ROH results
    rm -f -r *_ROH_Results.txt
Resubmitted -- was it possible the ROH code hadn't updated yet?

chrom = "CM068641.1"
chrom_length = 137863001 ## awk -v var='CM068641.1' '$1 == var {{print $2; exit}}' Reference_HepPer_Chroms_Lengths.txt
Var_file = "temp/20250422_CM068641.1_Var_Only.txt"
Aln_file = "temp/20250422_CM068641.1_Aln_Only.txt"
spec_name = "HepPer"
    df = pd.DataFrame({'variant': var_pos, 
        'var_aln_sum': aln_sums, 
        'score': score})
    result = pd.Series(true_groups).apply(
        lambda group: (
            df[groups == group]['var_aln_sum'].iloc[0] 
            if df[groups == group].index[0] == 0 
            else df.loc[df[groups == group].index[0] - 1, 'var_aln_sum'], 
            df[groups == group]['var_aln_sum'].iloc[-1]
        )
    )
    df_roh = pd.DataFrame({
        'start': ROH_start, 
        'end': ROH_end
    })
I think I found the problem -- It's tracking variants properly but I think it's because I'm using the alingment sums up to that point rather than the actual variant positions
result = pd.Series(true_groups).apply(
    lambda group: (
        df[groups == group]['variant'].iloc[0] 
        if df[groups == group].index[0] == 0 
        else df.loc[df[groups == group].index[0] - 1, 'variant'], 
        df[groups == group]['variant'].iloc[-1]
    )
)
Switched from 'var_aln_sum' to 'variant' in the dataframe to get true position of variants when tracking results
Resubmitted
It worked!
Some slight changes, but overall plot looks the same, with only a few ROH shifted
Now running the rest of the snakefile for HepPer, including heterozygosity, to plot results

Job finished, and it worked with the exception of the heterozygosity plots not properly being generated
I think the problem is with adding a column to the heterozygosity files for window length when excl ROH
Updated 20250123_Plot_het_per_chr.R and 20250325_Plot_het_whole_genome.R to add the extra column and also plot Het_Per_Kb_excl_ROH instead of just Het_Per_KB
I also updated the COMPILE_HET rule to make sure it also prints the new window size excl ROH and het per kb excl ROH


#### UPDATE ####
20250510 (May 10th, 2025)

Updated rules in Snakefile to create all rules for prepping files for MSMC and running MSMC


#### UPDATE ####
20250512 (May 12th, 2025)

Troubleshooting MSMC prep and run rules in Snakefile
Added outputs of last two rules in MSMC file prep, for creating main multihetsep file and the bootstraps, to rule_all
Updated HepPer_config.yml to include BOOTSTRAPPING_VALUES
Dry run appears to be successful
Submitted

Job to create sam file did not finish as it ran out of runtime
Specified in resources for that rule to have a runtime of 6hrs
Resubmitted

Modified 20250325_Plot_het_whole_genome.R to get the chromosomes listed instead of their lengths on the x-axis
Modified 20250106_Plot_ROH.R to clean up the legend name and change from chromosome names into sequential numbers

Submitted HypSab_config to redo jobs for plotting and run jobs for MSMC analysis


#### UPDATE ####
20250513 (May 13th, 2025)

Snakemake lists that the HepPer.sam file is incomplete and needs to be regenerated
Also only 13 fasta and 14 .fai files got made
Need to figure out what is making slurm fail on these
Similar situation with HypSab

I upped the runtime for creating the .SAM files to 12 hours and resubmitted
JOBIDs for SAM file jobs: 9329411 and 9329494

Updated 20250424_Taxon_Figure.R code to generate stacked barchart for the number of reference genomes


#### UPDATE ####
20250514 (May 14th, 2025)

It appears that the .SAM file for HypSab was succesfully made, it just needed higher running time
Have issues where some jobs are consistently failing, with error SLURM status is: 'FAILED'
Found the error in PLOT_ROH rule -- syntax error in update I had made, and subsequently fixed
Error for  EXTRACT_CHR_FASTA:
    [E::fai_load3_core] Failed to open FASTA file sharks/HemOce/temp/GCF_020745735.1.fa
Consistently this way with some chromosomes in HemOce, HepPer, and HypSab
Potentially due to not defining parameters for this rule
Added params and resubmitted
Slurm is still failing out on some iterations of this rule with the same error
Updated EXTRACT_CHR_FASTA into two rules, UNZIP_REFERENCE and EXTRACT_CHR_FASTA, and used temp command to remove .fa file when done to see if this fixes the error

CALC_HET_PER_CHR rule is failing for a few chromosome in HypSab due to following error:
 File "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/20250325_find_het_per_chr_V3.py", line 109, in calc_het
     41     variant_positions = single_chr_df.iloc[:, 2].astype(int)
     42                         ~~~~~~~~~~~~~~~~~~^^^^^^
     IndexError: single positional indexer is out-of-bounds

PILEUP VARIANTS rule is failing for HepPer with following error:
[E::hts_open_format] Failed to open file "sharks/HepPer/temp/CM068638.1.vcf.gz" 


#### UPDATE ####
20250515 (May 15th, 2025)

Fixed syntax problems resulting in error for HepPer and resubmitted

For HypSab, error is due to 2 chromosomes having no alingments or variants present
NC_082738.1
NC_082740.1
Both of these chromosomes are sex chromosomes -- X1 and Y
    awk -v var="NC_082738.1" '{if ($6 == var) {print $0}}' HypSab_FASTGA.chain.chr.fltr.srt.paf
Awk shows 47 alignments in paf file for NC_082738.1 and 2 alingments for NC_082740.1
Removed the var and aln files from HypSab/temp and re-ran on snakemake to only regenerate those files in case something happened
Re-made the ALn and Var files and the same errors exist. It's not due to filtering $2 for scaffolds, as most alingments of NC_082738.1 are not from unplaced scaffolds
    awk -v var="NC_082738.1" '{if ($6 == var) {print $0}}' HypSab_FASTGA.chain.chr.fltr.srt.paf > test.paf
    k8 ../../paftools.js call test.paf
    0 reference bases covered by exactly one contig
    0 substitutions; ts/tv = NaN
    0 1bp deletions
    0 1bp insertions
    0 2bp deletions
    0 2bp insertions
    0 [3,50) deletions
    0 [3,50) insertions
    0 [50,1000) deletions
    0 [50,1000) insertions
    0 >=1000 deletions
    0 >=1000 insertions

Temporarily fixing issue with HypSab by only working with autosomal chromosomes
Change line 76 in 20250325_find_het_per_chr_V3.py from 
    for i in range(len(all_chromosomes)):
to
    for i in range(num_aut_chr):
so that it only calculates for autosomal chromosomes
Re-submitted

Rule INDEL_MASKER for HepPer failed due to error:
    ModuleNotFoundError: No module named 'matplotlib'
Install maplotlib
    pip install matplotlib
Re-ran

Rule BAM_CALLER_BED failed for HypSab and HepPer due to syntax error
Fixed syntax error.

Rule REMOVE_INDELS failed for HypSab and HepPer due to syntax error
Fixed syntax error

awk -v var="Chondrichthyes" '{if ($3 == var) {print $0}}' VGPPhase1-freeze-1.0.tsv | cut -f16 | grep -c GC
awk -v var="Chondrichthyes" '{if ($3 == var) {print $0}}' VGPPhase1-freeze-1.0.tsv | cut -f22 | tr ',' '\n' | tr -d ' ' | grep GC | wc -l
    GCA_035084065.1
    GCA_036365495.1
    GCA_020745765.1
    GCA_030684295.1
    GCA_964194165.1
    GCA_964214025.1
    GCA_035084135.1
    GCA_044704965.1
    GCA_036971175.1
    GCA_030144785.1
    GCA_030035685.1


#### UPDATE ####
20250516 (May 16th, 2025)

Working through submitting snakemake for HepPer, HypSab, and HemOce
Fixing syntax errors in Snakefile

Got one error for GENERATE_MASK rule:
    ModuleNotFoundError: No module named 'msprime'
Installed msprime in snakemake environment
    pip install msprime

Modified the 20250325_Plot_het_whole_genome.R script to plot chromosome names instead of length along chromosome

Ran into error when running for HepPer:
    FileNotFoundError: [Errno 2] No such file or directory: 'sharks/HepPer/Output_Bed/CM068653.1_95.txt'
Error in GENERATE_MASK rule
Updated shell in GENERATE_MASK:
    python {MASK_FILE_GENERATOR} --input_file {input} --Chromosome_name {wildcards.CHR} --output_BED {params.CLADE}/{params.SPEC_NAME}/MSMC/Output_Bed/{wildcards.CHR}

Getting an error in CALL_PLOIDY with HypSab
    log /rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/.snakemake/slurm_logs/rule_CALL_PLOIDY/sharks_HypSab_NC_082740.1/9436857.log
This worked fine with HepPer, and error is saying output file is missing that was generated with HepPer
I will try increasing latency-wait to 60s and resubmitting
    snakemake --executor slurm --jobs 25 --workflow-profile profiles --configfile HypSab_config.yml --latency-wait 60
This did not work, the same error persisted
    36 Waiting at most 60 seconds for missing files:
    37 sharks/HypSab/MSMC/Output_bam_caller_BED/NC_082740.1_unmodified_bed.txt.gz (missing locally)

Getting an error in BOOTSTRAPPING_MULTIHET_FILE for HepPer:
    .snakemake/log/2025-05-16T115347.359293.snakemake.log
Was able to fix this error by changing syntax in rule_all input
Resubmitted HepPer
Failed with error:
    0 Traceback (most recent call last):
     41   File "/rds/project/rds-8b3VcZwY7rY/users/ag2427/hpc-work/20241011_jaskaran_script_msmc_fileprep/block_bootstrap_mhs.py", line 44, in <module>
     42     data = pd.read_csv(inmhs, header = None,sep='\t') # load data
     43            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     44   File "/home/ag2427/.conda/envs/snakemake/lib/python3.11/site-packages/pandas/io/parsers/readers.py", line 1026, in read_csv
     45     return _read(filepath_or_buffer, kwds)
     46            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     47   File "/home/ag2427/.conda/envs/snakemake/lib/python3.11/site-packages/pandas/io/parsers/readers.py", line 620, in _read
     48     parser = TextFileReader(filepath_or_buffer, **kwds)
     49              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     50   File "/home/ag2427/.conda/envs/snakemake/lib/python3.11/site-packages/pandas/io/parsers/readers.py", line 1620, in __init__
     51     self._engine = self._make_engine(f, self.engine)
     52                    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     53   File "/home/ag2427/.conda/envs/snakemake/lib/python3.11/site-packages/pandas/io/parsers/readers.py", line 1898, in _make_engine
     54     return mapping[engine](f, **self.options)
     55            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     56   File "/home/ag2427/.conda/envs/snakemake/lib/python3.11/site-packages/pandas/io/parsers/c_parser_wrapper.py", line 93, in __init__
     57     self._reader = parsers.TextReader(src, **kwds)
     58                    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     59   File "parsers.pyx", line 581, in pandas._libs.parsers.TextReader.__cinit__
     60 pandas.errors.EmptyDataError: No columns to parse from file

Getting an error in CALC_HET_PER_CHR for HemOce:
    Waiting at most 5 seconds for missing files:
     37 sharks/HemOce/20250428_NC_083454.1_het.txt (missing locally)
NC_083454.1
This is the y-chromosome


#### UPDATE ####
20250519 (May 19th, 2025)

Removed MSMC directory and all files within for HepPer, given that some of them were file size 0 and I think this was causing issues
Resubmitted

Set up Git extension for VScode and looked at how to push scripts to Github repo
    echo "# Heterozygosity_Code" >> README.md
    git init
    git add README.md
    git commit -m "first commit"
    git branch -M main
    git remote add origin https://github.com/alwaysamanda/Heterozygosity_Code.git
    git push -u origin main

Resubmitted code for HemOce -- Same step had an error as before with the y-chromosome
Problem is that the code is running only for autosomal chromosomes, but snakemake is expecting output files for all chromosomes including sex chromosomes
Modified HemOce_config.yml to include all chromosomes and a separate list of just the autosomal chromosomes

Updated all MSMC data prep rules to only look at autosomal chromosomes
Reran for HepPer

Write rule in Snakefile to create ALNplot for visual alignment of species analysed 

Created HydCol_config.yml and submitted to get additional files made

Had issues with ALNplot for HypSab
    MissingOutputException in rule ALNPLOT in file /rds/project/rds-p67MZilb2eQ/projects/VGP/heterozyg>
    39 Job 0 completed successfully, but some output files are missing. Missing files after 60 seconds.

Updated NarBan_config.yml and submitted

Created MobBir_config.yml and submitted

HydCol log: .snakemake/log/2025-05-19T184328.515320.snakemake.log
HepPer log: .snakemake/log/2025-05-19T185625.397439.snakemake.log
HemOce log: .snakemake/log/2025-05-19T184420.839445.snakemake.log
HypSab log: .snakemake/log/2025-05-19T185353.752621.snakemake.log
NarBan log: .snakemake/log/2025-05-19T212706.097623.snakemake.log
MobBir log: .snakemake/log/2025-05-19T212933.156237.snakemake.log


#### UPDATE ####
20250520 (May 20th, 2025)

Resubmitted HydCol and HepPer

ALNplot rule is not working 
Consistently getting this error:
    38 MissingOutputException in rule ALNPLOT in file /rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/Snakefile, line 237:
     39 Job 0 completed successfully, but some output files are missing. Missing files after 60 seconds. This might be due to filesystem latency. If that is the case, consider to increase>
     40 sharks/HydCol/HydCol_ALN.chain.pdf (missing locally, parent dir contents: CM068742.1.recode.vcf, C...
     )

Submitted ALNplot manually
    ALNplot -p -H500 HydCol_ALN.chain.1aln
It worked

Error in HydCol BOOTSTRAPPING_MULTIHET_FILE
    44   File "/rds/project/rds-8b3VcZwY7rY/users/ag2427/hpc-work/20241011_jaskaran_script_msmc_fileprep/block_bootstrap_mhs.py", line 55, in <module>
     45     block_indices_to_use = np.random.randint(1,num_windows,num_windows)
     46                            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     47   File "numpy/random/mtrand.pyx", line 782, in numpy.random.mtrand.RandomState.randint
     48   File "numpy/random/_bounded_integers.pyx", line 1334, in numpy.random._bounded_integers._rand_int64
     49 ValueError: low >= high
     50 RuleException:
     51 CalledProcessError in file /rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/Snakefile, line 737:

Error in HepPer BOOTSTRAPPING_MULTIHET_FILE
    num_windows=0
     43 Traceback (most recent call last):
     44   File "/rds/project/rds-8b3VcZwY7rY/users/ag2427/hpc-work/20241011_jaskaran_script_msmc_fileprep/block_bootstrap_mhs.py", line 87>
     45     if np.min([newmhs_np[i+1,0]-newmhs_np[i,0] for i in range(0,len(newmhs_np[:,0])-1)])<1:
     46        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     47   File "/home/ag2427/.conda/envs/snakemake/lib/python3.11/site-packages/numpy/core/fromnumeric.py", line 2953, in min
     48     return _wrapreduction(a, np.minimum, 'min', axis, None, out,
     49            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     50   File "/home/ag2427/.conda/envs/snakemake/lib/python3.11/site-packages/numpy/core/fromnumeric.py", line 88, in _wrapreduction
     51     return ufunc.reduce(obj, axis, dtype, out, **passkwargs)
     52            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     53 ValueError: zero-size array to reduction operation minimum which has no identity
This error is because the window size is bigger than the seq length

Updated BOOTSTRAPPING_MULTIHET_FILE to make the window size smaller (5e+05)
Will resubmit for HydCol and HepPer

Still getting an error on some bootstrapping for HepPer:
    line 44 of block_bootstrap_mhs.py
    data = pd.read_csv(inmhs, header = None,sep='\t') # load data
    pandas.errors.EmptyDataError: No columns to parse from file

Set up new rule in snakemake to get chrom lists:

Created SteTig_config.yml, CetMax_config.yml, MusAst_config.yml, PriJap_config.yml
However, despite these having alternate, they do not have reference haplotypes yet

Created HetFra_config.yml since that had not been created yet and ran to get ROH and heterozygosity plot

Getting errors in plotting Het for HetFra and NarBan
    Error in names(x) <- value : 
     55   'names' attribute [9] must be the same length as the vector [8]
     56 Calls: colnames<-
I think this error is because the het for these species was generated before the new code I put in to update Het
Will remove all files after generating .1aln and resubmitted

Errors for calculating heterozygosity for HemOce and MobBir
    38 MissingOutputException in rule CALC_HET_PER_CHR in file /rds/project/rds-p67MZilb2eQ/projects/VGP/>
     39 Job 0 completed successfully, but some output files are missing.

Updated MSMC rules for them to run after MSMC file prep is done

Error for HypSab in CALC_HET_PER_CHR:
     37 sharks/HypSab/20250425_NC_082740.1_het.txt (missing locally)
     38 MissingOutputException in rule CALC_HET_PER_CHR in file /rds/project/rds-p67MZilb2eQ/projects/VGP/heteroz>
     39 Job 0 completed successfully, but some output files are missing. Missing files after 60 seconds. This mig>
     40 sharks/HypSab/20250425_NC_082740.1_het.txt 

I updated the code to calculate heterozygosity to 20250520_find_het_per_chr_V4.py
Instead of running a loop through all chromosomes it will run separately for each autosomal chromosome
Hopefully this will run faster and reduce error

Error with BOOTSTRAPPING_MULTIHET_FILE for HepPer:
     40 Traceback (most recent call last):
     41   File "/rds/project/rds-8b3VcZwY7rY/users/ag2427/hpc-work/20241011_jaskaran_script_msmc_fileprep/block_bootstrap_mhs.py", line 44, in <module>
     42     data = pd.read_csv(inmhs, header = None,sep='\t') # load data
     43            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     44   File "/home/ag2427/.conda/envs/snakemake/lib/python3.11/site-packages/pandas/io/parsers/readers.py", line 1026, in read_csv
     45     return _read(filepath_or_buffer, kwds)
     46            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     47   File "/home/ag2427/.conda/envs/snakemake/lib/python3.11/site-packages/pandas/io/parsers/readers.py", line 620, in _read
     48     parser = TextFileReader(filepath_or_buffer, **kwds)
     49              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     50   File "/home/ag2427/.conda/envs/snakemake/lib/python3.11/site-packages/pandas/io/parsers/readers.py", line 1620, in __init__
     51     self._engine = self._make_engine(f, self.engine)
     52                    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     53   File "/home/ag2427/.conda/envs/snakemake/lib/python3.11/site-packages/pandas/io/parsers/readers.py", line 1898, in _make_engine
     54     return mapping[engine](f, **self.options)
     55            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     56   File "/home/ag2427/.conda/envs/snakemake/lib/python3.11/site-packages/pandas/io/parsers/c_parser_wrapper.py", line 93, in __init__
     57     self._reader = parsers.TextReader(src, **kwds)
     58                    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     59   File "parsers.pyx", line 581, in pandas._libs.parsers.TextReader.__cinit__
     60 pandas.errors.EmptyDataError: No columns to parse from file
It's because some of the input files are empty and as such cannot be read properly

#### UPDATE ####
20250521 (May 21st, 2025)

Error in MAKE_BAM for NarBan:
      36 [E::parse_cigar] CIGAR length too long at position 1 (345549512S)
     37 [W::sam_read1_sam] Parse error at line 448
     38 samtools sort: truncated file. Aborting
How do I get around this?


Error in CALC_HET_PER_CHR for HypSab for Chr CHR=NC_082738.1:
     36 Traceback (most recent call last):
     37   File "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/20250520_find_het_per_chr_V4.py", line 168, in <module>
     38     run_het_calculations = calc_het(dat, roh_data, chr_file, current_window_length, current_window_interval, chrom, date)
     39                            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     40   File "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/20250520_find_het_per_chr_V4.py", line 95, in calc_het
     41     variant_positions = single_chr_df.iloc[:, 2].astype(int)
     42                         ~~~~~~~~~~~~~~~~~~^^^^^^
     43   File "/home/ag2427/.conda/envs/snakemake/lib/python3.11/site-packages/pandas/core/indexing.py", line 1184, in __getitem__
     44     return self._getitem_tuple(key)
     45            ^^^^^^^^^^^^^^^^^^^^^^^^
     46   File "/home/ag2427/.conda/envs/snakemake/lib/python3.11/site-packages/pandas/core/indexing.py", line 1690, in _getitem_tuple
     47     tup = self._validate_tuple_indexer(tup)
     48           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     49   File "/home/ag2427/.conda/envs/snakemake/lib/python3.11/site-packages/pandas/core/indexing.py", line 966, in _validate_tuple_indexer
     50     self._validate_key(k, i)
     51   File "/home/ag2427/.conda/envs/snakemake/lib/python3.11/site-packages/pandas/core/indexing.py", line 1592, in _validate_key
     52     self._validate_integer(key, axis)
     53   File "/home/ag2427/.conda/envs/snakemake/lib/python3.11/site-packages/pandas/core/indexing.py", line 1685, in _validate_integer
     54     raise IndexError("single positional indexer is out-of-bounds")
     55 IndexError: single positional indexer is out-of-bounds
This error is because no variants were found for this chromosome
This is also a sex chromosome -- should not be running for this

Commented out all rules in relation to MSMC file prep and ran snakemake on all 7 elamsobranch species for testing
(HydCol, HepPer, HemOce, HetFra, HypSab, NarBan, MobBir) and got all FROH, ROH, and HET graphs and calculations for them

Started rules on MSMC prep for all 7 species
Got an error in the MAKE_BAM rule for NarBan:
    36 [E::parse_cigar] CIGAR length too long at position 1 (345549512S)
     37 [W::sam_read1_sam] Parse error at line 448
     38 samtools sort: truncated file. Aborting
This same error appears when I'm trying to use samtools view instead of samtools sort

Some rules failed for all species, so commented out rules after INDEX_BAM_FILE and re-ran for HepPer, HypSab, and HemOce

Job for HemOce finished with 100% completion
Added bake in all the rest of the rules for MSMC file formatting through bootstrapping the data_frame
Ran for HemOce, HepPer, and HypSab


#### UPDATE ####
20250522 (May 22nd, 2025)

Had to resubmit jobs for HemOce, HepPer, and HypSab since my computer stopped them running to do a software update overnight

Got advice from Chenxi about how to fix issues with NarBan cigar string
He recommended re-doing the SAM file generation for NarBan and adding the -L to the minimap2 command to address this
    minimap2 -t -L 16 -ax asm5 {input.REF} {input.ALT} > {output}
Removed NarBan.sam file
Switched over config.yaml to icelake-himem

CALL_PLOIDY jobs for NarBan are failing:
    36 Waiting at most 60 seconds for missing files:
     37 sharks/NarBan/MSMC/Output_bam_caller_BED/CM072967.1_unmodified_bed.txt.gz (missing locally)
     38 MissingOutputException in rule CALL_PLOIDY in file /rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/Snakefile, line 645:
     39 Job 0 completed successfully, but some output files are missing.
When reviewing rule and comparing it to Jaskaran's original code -- discovered that final command had bgzip instead of gzip -- fixed and resubmitted
Started manually reviewing the rule, and found out that the issue is in the .vcf.gz files for NarBan
They are empty so nothing gets called from bcftools call
Problem is traced back to teh chromosome specific .bam files, which are empty
.fasta files appear correct
Deleted all files back to the .bam file for the whole chromosome and resubmitted
Canceled jobs for NarBan as the chromosome specific .bam files are still showing up empty

Running files for HypSab locally since none of the jobs are intensive computing wise

I'm going to test whether it would be faster to run MSMC data prep in a shell script rather than breaking it up into snakemake rules
I will create a shell script to run for HepPer and submit it separately
Created 20250522_msmc_prep_HepPer.sh
Submitted batch job 9760677


#### UPDATE ####
20250523 (May 23rd, 2025)

Job 9760677 still hadn't run, so I canceled it and lowered the memory and run time and resubmitted to cclake-himem
    Submitted batch job 9798729
Job ran and finished with exit code 0
Checked error files and job didn't work because modules weren't loading
There was also an odd error at the top of the error file:
    [W::bam_hdr_read] EOF marker is absent. The input is probably truncated

Resubmitted HypSab with snakemake to cluster
Finished, but said that it wasn't done
Realized it was because the "{CLADE}/{SPEC_NAME}/{CHROM}.recode.vcf" file was getting deleted
when it was getting compressed into a .gz file and it was running in a loop
Commented out in rule all:
    expand("{CLADE}/{SPEC_NAME}/{CHROM}.recode.vcf", CLADE=CLADE, SPEC_NAME=SPEC_NAME, CHROM=AUTO_CHROMS), 
Submitted HypSab
Submitted HemOce
Submitted MobBir

Boostrapping jobs failed for HepPer because:
    pandas.errors.EmptyDataError: No columns to parse from file
The multihet files for several chromosomes are empty
    ./CM068673.1_multihet.txt
    ./CM068644.1_multihet.txt
    ./CM068679.1_multihet.txt
    ./CM068662.1_multihet.txt
    ./CM068681.1_multihet.txt
    ./CM068639.1_multihet.txt
    ./CM068659.1_multihet.txt
    ./CM068667.1_multihet.txt
    ./CM068663.1_multihet.txt
    ./CM068670.1_multihet.txt
    ./CM068666.1_multihet.txt
    ./CM068660.1_multihet.txt
    ./CM068641.1_multihet.txt
All of these chromosomes are chromosomes with long ROH spanning almost the whole chromosome
But not all ROH chromosomes have empty files -- e.g.
    CM068648.1 (chr12)
    CM068671.1 (chr35)
How do I work around this problem?

Realized that all the primary bootstrapping files have been generated for HypSab, just waiting for bootstrapping
Updating rules RUN_PRIMARY_MSMC and PLOT_MSMC to make sure they run
Adding generation time and mutation rate to the HypSab_config.yml file as parameters for plotting MSMC
    MUTATION_RATE: 1.25e-8
    GENERATION_TIME: 18.8
Transferred build (for msmc2 software) from hpc-work to heterozygosity folder using mv command
Ran for HypSab, commenting out boostrapping just to get primary results
PLOT_MSMC for HypSab failed with error:
     47 Error in dat$left_time_boundary/mu : 
     48   non-numeric argument to binary operator
Fixed by loading in mu and gen_time with as.numeric()
It worked!

Will now submit for HemOce and MobBir 
For now will keep bootstrapping commented out and focus on getting primary plotted for the report


#### UPDATE ####
20250525 (May 25th, 2025)

Re-ran PLOT_MSMC for HypSab, and HydCol with new y-axis limit
Submitted for MobBir and HemOce to get PLOT_MSMC run

Submitted MSMC Prep snakemake rules for HetFra
MAKE_BAM failed for HetFra for the same reason that it did for NarBan -- Cigar length too long
Removed SAM file and resubmitted with MAKE_SAM

The Call plidy jobs for HetFra are failing due to missing file error:
    sharks/HetFra/MSMC/Output_bam_caller_BED/CM070964.1_unmodified_bed.txt.gz (missing locally)
     38 MissingOutputException in rule CALL_PLOIDY in file 


#### UPDATE ####
20250526 (May 26th, 2025)

Running commands in CALL_PLOIDY individually for HetFra to see if I can get it working:
    bcftools call --ploidy 2 -c -V indels sharks/HetFra/temp/CM070962.1_pileup.vcf.gz
This came back empty
It appears the chromosome specific bam files are empty
the whole genome bam file is fine
    samtools view -H "sharks/HetFra/HetFra.bam" | grep "@SQ" | cut -f 2 | cut -d ':' -f 2
The HetFra bam file doesn't have any @SQ lines in the header -- no chromosomes are defined
I think the problem is in the minimap2 command to make the SAM file using the -L command
Potentially the issue is that "-t -L 16" broke up the threads specification instead of "-L -t 16"
Ran MAKE_SAM rule with:
    minimap2 -L -t 16 -ax asm5 {input.REF} {input.ALT} > {output}
Failed with non-zero exit status 137 -- I believe this is out of memory
Specified memory in MAKE_SAM resources
    mem_mb=100000
Resubmitted

It appears to be working, and generating @SQ lines
Also submitted NarBan and HepPer to get the sam files made


#### UPDATE ####
20250527 (May 27th, 2025)

It appears that the SAM files for NarBan, HepPer, and HetFra have been made successfully
Will continue them in Snakemake pipeline
Submitted HetFra
Removed previously made files for HepPer and NarBan
    .bam, .bam.bai files, vcf files, fasta and fa files, etc
    All files which had previously been made during MSMC file prep

Failed MAKE_BAM rule for HetFra with error:
    36 [E::parse_cigar] CIGAR length too long at position 1 (268907262H)
     37 [W::sam_read1_sam] Parse error at line 1513
     38 samtools sort: truncated file. Aborting
     39 RuleException:

Failed MAKE_BAM for NarBan with error:
    36 [E::parse_cigar] CIGAR length too long at position 1 (345549512S)
     37 [W::sam_read1_sam] Parse error at line 448
     38 samtools sort: truncated file. Aborting

This error was the exact same as before for NarBan, despite having added the -L command to minimap2

Error for separating BAM file by CHR for HepPer
     36 [E::idx_find_and_load] Could not retrieve index file for 'sharks/HepPer/HepPer.bam'
     37 samtools view: Random alignment retrieval only works for indexed SAM.gz, BAM or CRAM files.
    
Attempting to run boostrapping for HydCol MSMC2
Missing bootstrapping multihet files for just chr40 (CM068781.1)
Generating bootstrap files and then running boostrapping analysis

Got error for bootstrapping chr40:
    Traceback (most recent call last):
  File "/rds/project/rds-8b3VcZwY7rY/users/ag2427/hpc-work/20241011_jaskaran_script_msmc_fileprep/block_bootstrap_mhs.py", line 55, in <module>
    block_indices_to_use = np.random.randint(1,num_windows,num_windows)
                           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "numpy/random/mtrand.pyx", line 782, in numpy.random.mtrand.RandomState.randint
  File "numpy/random/_bounded_integers.pyx", line 1334, in numpy.random._bounded_integers._rand_int64
    ValueError: low >= high
I think the error is due to how short the chromosome is: 812,877 bases
Changed window size from 5e+05 to 2.5e+05
Resubmitted
This appears to have worked!

Failed out on running bootstrapping MSMC however due to an out of memory error
Added resource parameters for RUN_BOOTSTRAPPING_MSMC and resubmitted
Submitting to run Boostrapping MSMC for HypSab

I believe I found the link for the reads for HepPer:
https://www.ncbi.nlm.nih.gov/sra/SRX24176190[accn]
Need to download them and start mapping reads to reference

Richard stated that I need to create a script to create the MSMC2 input file and mask etc starting from the paf file I get from FASTGA
What does the end file need to look like?
Steps:
    Generate a mask and VCF file
        mask specifies what regions should be called (3 columns, chr, start, and end)
    Run generate_multihetsep.py, which merges VCF and mask files together, and also performs simple trio-phasing.
    Should output multihepsep file


#### UPDATE ####
20250528 (May 28th, 2025)

Finished generating Bootstrapping multihet files for HydCol
Submitted to run Bootstrapping MSMC for HydCol

Started writing new rules for doing my own pipeline to prep MSMC input from FASTGA generated files
Created 20250528_ROH_Masker.py to create mask files which exclude all ROH regions
rule SEP_PAF_BY_CHR:
    input:
        "{CLADE}/{SPEC_NAME}/{SPEC_NAME}_FASTGA.chain.chr.fltr.srt.paf"
    output:
        "{CLADE}/{SPEC_NAME}/MSMC/paf_files/{SPEC_NAME}_{CHROM}.paf"
    shell:
        """
        mkdir -p "{CLADE}/{SPEC_NAME}/MSMC/paf_files"
        awk -v var={wildcards.CHROM} '{{if($6 == var) print}}' {input} > {output}
        """

rule GEN_CHROM_VCF:
    input:
        paf="{CLADE}/{SPEC_NAME}/MSMC/paf_files/{SPEC_NAME}_{CHROM}.paf",
        refseq=expand("/rds/project/rds-p67MZilb2eQ/projects/VGP/250430.VGP-Phase1/alignment/reference/{CLADE}/{REF_NAME}.fa.gz", CLADE=CLADE, REF_NAME=REF_NAME)
    output:
        "{CLADE}/{SPEC_NAME}/MSMC/vcf_files/{SPEC_NAME}_{CHROM}.vcf.gz"
    shell:
        """
        mkdir -p "{CLADE}/{SPEC_NAME}/MSMC/vcf_files"
        k8 paftools.js call -s {SPEC_NAME} -f {input.refseq} {input.paf} > {CLADE}/{SPEC_NAME}/MSMC/vcf_files/{SPEC_NAME}_{CHROM}.vcf
        gzip {CLADE}/{SPEC_NAME}/MSMC/vcf_files/{SPEC_NAME}_{CHROM}.vcf
        """

rule GEN_ROH_NEGATIVE_MASK:
    input:
        "{CLADE}/{SPEC_NAME}/{TODAY_DATE}_{CHROM}_ROH_Results.txt"
    output:
        "{CLADE}/{SPEC_NAME}/MSMC/ROH_Negative_Mask/{SPEC_NAME}_{CHROM}_ROH_Mask.bed.gz"
    shell:
        """
        mkdir -p "{CLADE}/{SPEC_NAME}/MSMC/ROH_Negative_Mask"
        python {ROH_MASKER}
        gzip {CLADE}/{SPEC_NAME}/MSMC/ROH_Negative_Mask/{SPEC_NAME}_{CHROM}_ROH_Mask.bed
        """

rule GEN_INPUT_MASK_VCF:
    input:
        "{CLADE}/{SPEC_NAME}/MSMC/vcf_files/{SPEC_NAME}_{CHROM}.vcf.gz"
    output:
        mask="{CLADE}/{SPEC_NAME}/MSMC/Output_Mask/{SPEC_NAME}_{CHROM}.bed.gz", 
        vcf="{CLADE}/{SPEC_NAME}/MSMC/vcf_files/{SPEC_NAME}_{CHROM}_only_snps.vcf.gz"
    shell:
        """
        mkdir -p "{CLADE}/{SPEC_NAME}/MSMC/Output_Mask"
        bcftools call --ploidy 2 -c -V indels {input} | {BAM_CALLER} 1 {output.mask} | bgzip -c > {output.vcf}
        """

rule MAIN_MULTIHETSEP:
    input:
        vcf="{CLADE}/{SPEC_NAME}/MSMC/vcf_files/{SPEC_NAME}_{CHROM}_only_snps.vcf.gz", 
        ROH_negative_mask="{CLADE}/{SPEC_NAME}/MSMC/ROH_Negative_Mask/{SPEC_NAME}_{CHROM}_ROH_Mask.bed.gz", 
        mask="{CLADE}/{SPEC_NAME}/MSMC/Output_Mask/{SPEC_NAME}_{CHROM}.bed.gz"
    output:
        "{CLADE}/{SPEC_NAME}/MSMC/Output_NEW_multihetsep/{CHROM}_multihet.txt"
    shell:
        """
        mkdir -p "{CLADE}/{SPEC_NAME}/MSMC/Output_NEW_multihetsep"
        python {GENERATE_MULTIHETSEP} {input.vcf} --negative_mask={input.ROH_negative_mask} --mask={input.mask} > {output}
        """

Added rules into Snakefile, and commented out old rules
Commented out all rules but SEP_PAF_BY_CHR and submitted for HydCol as a test run 
This worked succesfully!

Now will try GEN_CHROM_VCF rule for HydCol
Failed with error:
    36 ERROR: fail to open file 'paftools'.
     37 ERROR: failed to read file 'paftools'
Syntax error saying paftools instead of paftools.js
Fixed and resubmitted
It worked!

Submitted for GEN_INPUT_MASK_VCF 
Failed with error: 
    Wrong number of PL fields? nals=2 npl=-1


#### UPDATE ####
20250529 (May 29th, 2025)

The error I was getting for GEN_INPUT_MASK_VCF is because the VCF wasn't created from a BAM file, and therefore doesn't contain GL or PL columns (genotype likelihoods)
I think I will use bcftools view instead of call, and filter to not use variants
Original:
    bcftools call --ploidy 2 -c -V indels {input} | {BAM_CALLER} 1 {output.mask} | bgzip -c > {output.vcf}
New:
    bcftools view -V indels -m2 -M2 -v snps {input} | {BAM_CALLER} 1 {output.mask} | bgzip -c > {output.vcf}
Failed with this error:
    Error: only supply one of --include-types, --exclude-types options
     37 RuleException:
Fixed to only exclude indels
    bcftools view -V indels -m2 -M2 {input} | {BAM_CALLER} 1 {output.mask} | bgzip -c > {output.vcf}
Resubmitted
It appears to have worked!

Modified 20250528_ROH_Masker.py to account for if no ROH are detected
Submitted for GEN_ROH_NEGATIVE_MASK
Failed with error:
    gzip: sharks/HydCol/MSMC/ROH_Negative_Mask/HydCol_CM068756.1_ROH_Mask.bed: No such file or directory
Error was due to syntax error in output file name in 20250528_ROH_Masker.py
Fixed and resubmitted
Appears to have worked!

Will now attempt to run primary multihepsep file generation
Says it finished without error but all the files are empty
Problem is tracing back to the VCF files generated in GEN_INPUT_MASK_VCF,
They have headers but the VCF body is empty
The vcf files from GEN_CHROM_VCF do not have this issue
    bcftools view -V indels -m2 -M2 {input.vcf}
This command appears to work, so something is going wrong with the BAM Caller
Both outputs from BAM_CALLER are empty

Testing alternative way to get mask file, using mask code from Jaskaran's pipeline
    python /rds/project/rds-8b3VcZwY7rY/users/ag2427/hpc-work/msmc-tools/generate_multihetsep.py vcf_files/HydCol_CM068742.1.vcf.gz > test_CM068742.1_multihet.txt
    python /rds/project/rds-8b3VcZwY7rY/users/ag2427/hpc-work/20241011_jaskaran_script_msmc_fileprep/mask_file.py --input_file test_CM068742.1_multihet.txt --Chromosome_name CM068742.1 --output_BED Output_Mask/
This generates a negative mask of high variant hotspots to be avoided

Given that a mask is supposed to say which are good regions to call, 
I'm going to try generating multihepsep files with just ROH negative mask
This resulted in only some multihet files getting generated -- others were empty
Possibly need a mask to state what to include
Will create new code and use system to only look at regions which are aligned
Wrote 20250529_Mask_Maker.py to address this


Created 20250529_Plot_MSMC_Bootstrap.R to plot primary with bootstrapping separately from only the primary MSMC
Created separate rule to plot with bootstrapping


python {CALC_HET_PER_CHR_PY} {input.VAR_FILE} {input.CHROM_LENGTH_FILE} {params.CLADE} {params.SPEC_NAME} {params.TODAY_DATE} {params.WINDOW_LENGTH} {params.WINDOW_INTERVAL} {wildcards.CHROM} {input.ROH}


#### UPDATE ####
20250530 (May 30th 2025)

Due to syntax error in rm command, removed all old versions of codes from directory
Using this as an opportunity to re-organize the heterozygosity directory
    mkdir config_files
This will be the directory which contains the config files for all species I am analyzing
    mkdir profiles
This will contain the config.yaml file which specifies the slurm parameters:
    default-resources:
        slurm_account: "DURBIN-SL2-CPU"
        slurm_partition: "cclake"
        mem_mb: 100000
        runtime: "2h"
Moved MobBir_config.yml to config_files
Recreated config files for the other 6 test elasmobranch species
Worked on adding all files to github to make sure that data loss doesn't occur again

Resubmitted FASTGA for all 7 elasmobranch species

Went through old snakemake file and realized that the chrom list files were not required for any rules
Deleted them from the current config files

Looking at bony fish species to see what species have both reference and alternate haplotypes available
    Lycodes pacificus (blackbelly eelpout) -- LycPac
        Reference = GCA_028022725.1.fa.gz
        Alternate = GCA_028021495.1.fa.gz
    Acipenser ruthenus (sterlet)
        Reference = GCF_902713425.1.fa.gz
        Alternate = GCA_902713435.2.fa.gz
    Hoplias malabaricus (trahira)
        Reference = GCA_029633875.1.fa.gz
        Alternate = GCA_029633855.1.fa.gz
    Salminus brasiliensis (dorado)
        Reference = GCA_030463535.1.fa.gz
        Alternate - GCA_030448965.1.fa.gz
    Amia calva (bowfin)
        Reference = GCA_036373705.1.fa.gz
        Alternate - GCA_036365475.1.fa.gz
    Fundulus diaphanus (banded killifish)
        Reference = GCA_037039145.1.fa.gz
        Alternate = GCA_037038625.1.fa.gz
    Cyprinella venusta (blacktail shiner)
        Reference = GCA_038024135.1.fa.gz
        Alternate - GCA_038021265.1.fa.gz

Added in ALNCHAIN and ALNtoPAF and submitted for HydCol

Created LycPac_config.yml 
    Couldn't find generation time anywhere, so estimated it to be 2yrs based on estimation of lifespan of 5 years from fishbase
Submitted for LycPac

FILTER_PAF_CHR_ONLY failed for HepPer due to error:
    /usr/bin/bash: -c: line 1: syntax error near unexpected token `('


#### UPDATE ####
20250531 (May 31st, 2025)

Switched partition to cclake-himem from cclake
Submitting files to finish filtering/sorting the paf file, getting the ALN plot, and getting ALN and VAR files
Submitting HydCol, HepPer, LycPac, MobBir

FASTGA failed for LycPac and MobBir due to error:
    FATAL ERROR: failed to write temporary file /tmp/OneTextSchema-710457.schema errno 28
After looking into error, I resubmitted for LycPac and it appears to be running without error (at least 9min in)
Resubmitted for MobBir

Getting ALN and VAR files failed for HepPer and HydCol due to not finding paftools
Need to re-download script
Re-downloaded paftools.js from minimap2/misc on github
Resubmitted for HydCol and HepPer
Jobs completed apparently without error

Added in rules for ROH calculation and plotting
Took out Today_date variable in rules, and made respective modifications to remove today_date as a variable in 20250106_Plot_ROH.R and the ROH calculation python script
Submitted for HydCol and HepPer

Submitted for NarBan, HetFra, and HypSab

the rule for PLOT_ROH failed for HepPer due to the following error:
    Error in `geom_bar()`:
     48 ! Problem while computing aesthetics.
     49 ℹ Error occurred in the 1st layer.
     50 Caused by error in `check_aesthetics()`:
     51 ! Aesthetics must be either length 1 or the same as the data (47).
     52 ✖ Fix the following mappings: `x`.
     53 Backtrace:
Problem is because mt genome was also sequenced and recorded as a chromosome -- have to alter script ot make sure it is excluded
Modified line 40 in 20250106_Plot_ROH.R script to:
    Chrom <- as.data.frame(Chrom[1:num_all_chr,])
This way it will exclude at mtDNA
Resubmitted for HepPer

HydCol, LycPac, and NarBan finished successfully
HepPer finished successfully
Checked maps and realize they have issue of ROH plotting separately to chromosomes and not on top of them

Removed HepPer ROH Map pdf and resubmitted with corrections to 20250106_Plot_ROH.R to see if error is fixed
Fixed

Removed ROH maps for HydCol and NarBan
Resubmitted for both

Submitted for LycPac to run ROH codes

Submitted for MobBir and HemOce
HemOce failed with error:
    FATAL ERROR: failed to write temporary file /tmp/OneTextSchema-1122153.schema errno 28
Resubmitted

PLOT_ROH for HemOce failed with this error:
    47 Error in `levels<-`(`*tmp*`, value = as.character(levels)) : 
     48   factor level [54] is duplicated
     49 Calls: factor
     50 Execution halted

HypSab failed in PLOT_ROH with this error:
    Error in `levels<-`(`*tmp*`, value = as.character(levels)) : 
     48   factor level [34] is duplicated
     49 Calls: factor
     50 Execution halted

These two errors are because there are multiple sex chromosomes for these species
And that means that 'Sex' as a label is duplicated, so I can't use factor


#### UPDATE ####
20250602 (June 2nd, 2025)

Used interactive session to fix 20250106_Plot_ROH.R script to fix errors in HemOce and HypSab

Pushed changes to github
Process to push updates:
    git commit -a -m "INSERT MESSAGE FOR UPDATES HERE"
    git push

Added in commands for calculating FROH to Snakefile
Modified them to remove today_date

Submitted for HydCol, HepPer, HemOce, HypSab, NarBan, MobBir, and HetFra

Modified 20250520_find_het_per_chr_V4.py to remove the today_date and make it run for each individual chromosome rather than loop over all autosomal chromsomes
Modified 20250325_find_het_whole_genome_V3.py to remove today_date

Modified 20250123_Plot_het_per_chr.R to remove today_date from it, and to save outputs in png formation
Modified 20250325_Plot_het_whole_genome.R to remove today_date from it, and save output in png format

Added all outputs from heterozygosity rules to rule all
Submitted for HydCol, HepPer, HemOce, HypSab, MobBir, NarBan, and HetFra
Submitted for LycPac

Updated 20250529_Mask_Maker.py script and added rules for MSMC data prep, including masking, into snakefile
Updated RUN_PRIMARY_MSMC and PLOT_MSMC to get rid of today_date

Created input files for 
    Acipenser ruthenus (AciRut)
    Hoplias malabaricus (HopMal)
    Salminus brasiliensis (SalBra)
    Amia calva (AmiCal)
    Fundulus diaphanus (FunDia)
    Cyprinella venusta (CypVen)

Couldn't find any information for generation time for HopMal or any other species in its GEN_INPUT_MASK_VCF
Looked into the literature and found a paper by José Luís Costa Novaes & Edmir Daniel Carvalho
States most often HopMal is caught between 2.5-4.5yrs, but can live up to 9yrs
Estimating for now that their generation time is ~1.5 years


#### UPDATE ####
20250603 (June 3rd, 2025)

All jobs from yesterday finished
However whole genome FROH has not been generated for any species
PLOT_HET_PER_CHR also hasn't worked
Will re-run for HydCol to see what's going on

PLOT_HET_PER_CHR failed with error:
    sharks/HydCol/HydCol_CM068768.1_Het_Map.png (missing locally)

The whole genome map is being saved as a png, but the individual chromosome maps are being saved as svg
This is why png files per chr cannot be found
Due to syntax error in the 20250123_Plot_het_per_chr.R code
Fixed and resubmitted

Failed for whole genome FROH with error:
    Error in file(file, "rt") : cannot open the connection
     48 Calls: read.table -> file
     49 In addition: Warning message:
     50 In file(file, "rt") :
     51   cannot open file 'sharks/GCA_035084275.1/GCA_035084275.1_aligned.mm2.vcf': No such file or directory
Created 20250603_FROH_Calc_Whole_Genome_V2.R to clean up FROH code and remove unnecessary lines with vcf and some other lines
Modified WHOLE_FROH rule to use new script
Resubmitted
It worked!

Added all outputs for MSMC data prep to rule all
Submitted for HydCol
error on GEN_ROH_NEGATIVE_MASK
      38 ValueError: Unable to parse string "sharks/HydCol/CM068744.1_ROH_Results.txt"
     39 
     40 During handling of the above exception, another exception occurred:
     41 
     42 Traceback (most recent call last):
     43   File "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/20250528_ROH_Masker.py", line 27, in <module>
     44     chrom_length = pd.to_numeric(sys.argv[4])
     45                    ^^^^^^^^^^^^^^^^^^^^^^^^^^
     46   File "/home/ag2427/.conda/envs/snakemake/lib/python3.11/site-packages/pandas/core/tools/numeric.py", line 232, in to_numeric
     47     values, new_mask = lib.maybe_convert_numeric(  # type: ignore[call-overload]
     48                        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     49   File "lib.pyx", line 2433, in pandas._libs.lib.maybe_convert_numeric
     50 ValueError: Unable to parse string "sharks/HydCol/CM068744.1_ROH_Results.txt" at position 0
Syntax error in 20250528_ROH_Masker.py
Fixed

Error in MAKE_MASK:
     /usr/bin/bash: line 2: /rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/20250529_Mask_Maker.py: Permission denied
Realized this error was because I forgot to put 'python' before the script. Fixed

Resubmitted with both errors fixed
Error on GEN_ROH_NEGATIVE_MASK:
    Traceback (most recent call last):
     37   File "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/20250528_ROH_Masker.py", line 45, in <module>
     38     CREATE_NEGATIVE_BED_FILE(roh_dat, chrom, clade, spec_name)
     39   File "/rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/20250528_ROH_Masker.py", line 40, in CREATE_NEGATIVE_BED_FILE
     40     dat.iloc[:, 1] = int(dat.iloc[:, 1])
     41                      ^^^^^^^^^^^^^^^^^^^
     42   File "/home/ag2427/.conda/envs/snakemake/lib/python3.11/site-packages/pandas/core/series.py", line 248, in wrapper
     43     raise TypeError(f"cannot convert the series to {converter}")
     44 TypeError: cannot convert the series to <class 'int'>
Changed int() to pd.to_numeric()
Resubmitted

Appears to be working for HydCol, so also submitting for HepPer, HypSab, HemOce, HetFra, NarBan, MobBir
Also submitted for LycPac, AciRut, and HopMal

WHOLE_FROH failed for HepPer with error:
    head: cannot open 'sharks/HepPer/temp/CM068729.1_Var_Only.txt' for reading: No such file or directory
CM068729.1 is the mtDNA -- shouldn't have been recorded. Need to find a way to make sure this can be ignored
Manually removed CM068729.1 from the chrom_lists/HepPer_chroms.txt
For other species, modified GET_CHROM_LISTS:
    zcat < {input} | grep '>{params.CHROM_START_CHR}' | head -n {params.NUM_ALL_CHR} | sed 's/^>//' > {output}

For SalBra generation time, found paper by Tos et al. estimating age of sexual maturity at 1.82years
Estimated generation time as 2yrs
Submitted for SalBra, AmiCal, and CypVen


#### UPDATE ####
20250604 (June 4th, 2025)

Submitted for HydCol
Job finished, and was hypothetically successful, but the majority of multihet files came out empty
Chromosome paf files appear to have been generated succesfully
ROH negative mask files appear to have been made succesfully
Aln mask files appear to have been made succesfully
Multihet files were made for 17/40 chromosomes
Realized I hadn't modified the MAIN_MULTIHETSEP rule to include the Aln mask - it was only using the ROH negative mask
Updated rule
    python {GENERATE_MULTIHETSEP} {input.VCF} --mask={input.MASK} --negative_mask={input.ROH_negative_mask} > {output}
Resubmitted
First job ran for Chrom CM068768.1 -- multihet file still came out empty
Canceled rest of jobs

Possible that the problem is due to the VCF files listing all variants as homozygous (1/1)
Will change to 1/0
Modified GEN_CHROM_VCF, but will first see if this works on just Chrom CM068768.1
    k8 ../../../../paftools.js call -s HydCol -f /rds/project/rds-p67MZilb2eQ/projects/VGP/250430.VGP-Phase1/alignment/reference/sharks/GCA_035084275.1.fa.gz ../paf_files/HydCol_CM068768.1.paf > HydCol_CM068768.1.vcf
    gzip HydCol_CM068768.1.vcf
    zcat HydCol_CM068768.1.vcf.gz | sed  's/1\/1/0\/1/g' | bgzip -c > HydCol_CM068768.1_restructured.vcf.gz
    python /rds/project/rds-8b3VcZwY7rY/users/ag2427/hpc-work/msmc-tools/generate_multihetsep.py sharks/HydCol/MSMC/vcf_files/HydCol_CM068768.1_restructured.vcf.gz --mask=sharks/HydCol/MSMC/Aln_Mask/CM068768.1_Aln_Mask.bed --negative_mask=sharks/HydCol/MSMC/ROH_Negative_Mask/HydCol_CM068768.1_ROH_Mask.bed.gz > test_multihet.txt
test_multihet.txt still comes out empty

Wondering if it is coming out empty because the ROH negative mask script is empty
Will remove that and see if that fixes the issue
    python /rds/project/rds-8b3VcZwY7rY/users/ag2427/hpc-work/msmc-tools/generate_multihetsep.py sharks/HydCol/MSMC/vcf_files/HydCol_CM068768.1_restructured.vcf.gz --mask=sharks/HydCol/MSMC/Aln_Mask/CM068768.1_Aln_Mask.bed > test_multihet.txt
Success!
Modified MAIN_MULTIHETSEP rule
        if [ stat -c %c {input.ROH_negative_mask} -eq 51 ]; then
            python {GENERATE_MULTIHETSEP} {input.VCF} --mask={input.MASK} > {output}
        else
            python {GENERATE_MULTIHETSEP} {input.VCF} --mask={input.MASK} --negative_mask={input.ROH_negative_mask} > {output}
        fi
Resubmitted for HydCol
One job was successful, three put out empty files
Canceled job
The empty files were again with empty ROH_negative_mask files
Altered if statement
    if [ "$(stat -c %s {input.ROH_negative_mask})" -eq 51 ]; then
Resubmitted
It's working!
Submitted for HepPer, NarBan, MobBir, HetFra, HypSab, and HemOce

Downloaded msmc-tools from github using git clone
Pulled msmc2 software using cp
    cp -R /rds/project/rd109/rds-rd109-durbin-group/software/msmc2 msmc2
Getting error message when attempting to run primary MSMC2 on NarBan
    core.exception.ArrayIndexError@model/data.d(188): index [4] exceeds array of length 4
    ----------------
    ??:? _d_arraybounds_indexp [0x62efc5]
    ??:? bool model.data.has_missing_data(in char[][], ulong[2]) [0x5b63fe]
    ??:? model.data.SegSite_t[][] model.data.readSegSites(immutable(char)[], ulong[2][], bool) [0x5b5cf8]
    ??:? model.data.SegSite_t[][] msmc2.readDataFromFiles(immutable(char)[][], ulong[2][], bool) [0x5da79a]
    ??:? void msmc2.parseCommandLine(immutable(char)[][]) [0x5d8c9a]
    ??:? _Dmain [0x5d85b3]
Possible that this issue is due to indels being present instead of just SNPs
Testing this:
    bcftools view -V indels HydCol_CM068742.1.vcf.gz -Oz -o CM068742.1_no_indels.vcf.gz
    python /rds/project/rds-8b3VcZwY7rY/users/ag2427/hpc-work/msmc-tools/generate_multihetsep.py sharks/HydCol/MSMC/vcf_files/CM068742.1_no_indels.vcf.gz --mask=sharks/HydCol/MSMC/Aln_Mask/CM068742.1_Aln_Mask.bed --negative_mask=sharks/HydCol/MSMC/ROH_Negative_Mask/HydCol_CM068742.1_ROH_Mask.bed.gz > test.txt
test.txt came out empty
Tried this same command with the vcf file containing indels, using standard output, and it came out with results and was not empty
Not sure why this is the case, I checked both the mask and negative ROH mask, and the mask includes the positions of the first handful of SNPs while the ROH mask does not remove them
Will try to change 1/1 to 0/1 to see if this is the issue
    zcat CM068742.1_no_indels.vcf.gz | sed  's/1\/1/0\/1/g' | bgzip -c > HydCol_CM068742.1_restructured.vcf.gz
    python /rds/project/rds-8b3VcZwY7rY/users/ag2427/hpc-work/msmc-tools/generate_multihetsep.py sharks/HydCol/MSMC/vcf_files/HydCol_CM068742.1_restructured.vcf.gz --mask=sharks/HydCol/MSMC/Aln_Mask/CM068742.1_Aln_Mask.bed --negative_mask=sharks/HydCol/MSMC/ROH_Negative_Mask/HydCol_CM068742.1_ROH_Mask.bed.gz > test.txt
This appears to have worked! Now to test if it will work in MSMC
    msmc2/build/release/msmc2 -t 12 --fixedRecombination -o sharks/HydCol/MSMC/Primary_Results/HydCol.msmc2 test.txt
It worked!
Implementing the changes as a rule in snakemake, to then run on HydCol as a test
    rule FILTER_VCF:
    input:
        "{CLADE}/{SPEC_NAME}/MSMC/vcf_files/{SPEC_NAME}_{CHROM}.vcf.gz"
    output:
        "{CLADE}/{SPEC_NAME}/MSMC/filtered_vcf_files/{SPEC_NAME}_{CHROM}_filtered.vcf.gz"
    shell:
        """
        mkdir -p {wildcards.CLADE}/{wildcards.SPEC_NAME}/MSMC/filtered_vcf_files
        bcftools view -V indels {input} -Oz -o {wildcards.CLADE}/{wildcards.SPEC_NAME}/MSMC/{wildcards.CHROM}_no_indels.vcf.gz
        zcat {wildcards.CLADE}/{wildcards.SPEC_NAME}/MSMC/{wildcards.CHROM}_no_indels.vcf.gz | sed  's/1\/1/0\/1/g' | bgzip -c > {output}
        """
Submitted with updated rules for HydCol
Appears to have finished without issue
Will plot to look at results
    Rscript 20250508_Plot_MSMC.R sharks HydCol 1.25e-8 18.6 sharks/HydCol/HydCol_MSMC2.png sharks/HydCol/MSMC/Primary_Results/HydCol.msmc2.final.txt
It worked!
Submitted for all of HydCol
It worked!

Submitted for HepPer, HetFra, HemOce, NarBan, MobBir, and HypSab

Filter_vcf rule failed for HepPer with error:
    gzip: sharks/HepPer/MSMC/CM068643.1_no_indels.vcf.gz: No such file or directory
     37 RuleException:
     38 CalledProcessError in file /rds/project/rds-p67MZilb2eQ/projects/VGP/heterozygosity/Snakefile, line 574:
Fixed syntax error


#### UPDATE ####
20250605 (June 5th, 2025)

Resubmitted for HepPer, HetFra, HemOce


