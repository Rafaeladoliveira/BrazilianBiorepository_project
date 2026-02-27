#!/bin/bash
#06_Admixture_Analysis

MAIN_DIR=$(pwd) #Define the home base (project main folder)

OUTDIR="ADMIXTURE"
mkdir -p $OUTDIR

#Gather required files from previous steps
#We need the YTCMergeData_LD (Final output of script 05) and Correct_ReferenceRaceFile.txt (Final output of script 03)
cp YTC_RefPanel_Merged/YTCMergeData_LD ADMIXTURE/
cp FinalReferenceDataset/Correct_ReferenceRaceFile.txt ADMIXTURE/

cd $OUTDIR

#Create the Pop-File (Supervised/Label file) for ADMIXTURE
##ADMIXTURE uses a .pop file where 'OWN' indicates the target samples and EUR/AFR/etc. indicates the reference samples.
##Note: The .pop file must follow the order of the .fam file EXACTLY.
awk '{print $1}' Admixture_Ready_Dataset.fam > samples_order.txt
awk 'NR==FNR {pop[$1]=$2; next} {if($1 in pop) print pop[$1]; else print "OWN"}' Correct_ReferenceRaceFile.txt samples_order.txt > FinalRacefile.pop

#Run ADMIXTURE for K=4 (number of ancestral groups present on Reference Panel - EUR, AFR, ASN, AMR)
K=4
admixture32 --cv YTCMergeData_LD.bed $K | tee log${K}.out

#Run the R script "AncestryAnalysis_ADMIXTURE.R" (available on /scripts directory). This script allow us to plot the final ADMIXTURE result.
module load R/4.2.1 #use your R version. This step is important if you are running your analysis on a cluster.
echo "Generating ancestry plots..."
Rscript ../scripts/AncestryAnalysis_ADMIXTURE.R
echo "R analysis complete."

cd "$MAIN_DIR"


