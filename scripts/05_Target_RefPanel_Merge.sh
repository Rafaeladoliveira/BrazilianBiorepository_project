#!/bin/bash
#05_Target_RefPanel_Merge

MAIN_DIR=$(pwd) #Define the home base (project main folder)

OUTDIR="YTC_RefPanel_Merged"
mkdir -p $OUTDIR

#Ensure that you have the files: ReferenceDatasetFinal (Final output of script 03) and FinalYTC (Final output of script 04)
cp FinalReferenceDataset/ReferenceDatasetFinal YTC_RefPanel_Merged/
cp TargetCohort_QC/FinalYTC YTC_RefPanel_Merged/

cd $OUTDIR

#Extract the variants present in FinalYTC from the ReferenceDatasetFinal
awk '{print$2}' FinalYTC.bim > FinalYTC.txt
plink --bfile ReferenceDatasetFinal --extract FinalYTC.txt --make-bed --out ReferenceDatasetFinal1

#Extract the variants present in ReferenceDatasetFinal1 from the FinalYTC
awk '{print$2}' ReferenceDatasetFinal1.bim > ReferenceDatasetFinal1.txt
plink --bfile FinalYTC --extract ReferenceDatasetFinal1.txt --make-bed --out FinalYTC2

#Set the reference genome
awk '{print$2,$5}' ReferenceDatasetFinal1.bim > ReferenceDatasetFinal1.txt
plink --bfile FinalYTC2 --reference-allele ReferenceDatasetFinal1.txt --make-bed --out FinalYTC2-adj

#Error: Duplicate ID "10:96695774:CTG:C". 
#Remove this variant from the file from FINALYTC2
cut -f2 FinalYTC2.bim | sort | uniq -d > duplicated_ids.txt
plink --bfile FinalYTC2 --exclude duplicated_ids.txt --make-bed --out FinalYTC3

#Try to run the command for set reference genome again
plink --bfile FinalYTC3 --reference-allele ReferenceDatasetFinal1.txt --make-bed --out FinalYTC3-adj

#Merge Datasets
plink --bfile FinalYTC3-adj --bmerge ReferenceDatasetFinal1.bed ReferenceDatasetFinal1.bim ReferenceDatasetFinal1.fam --allow-no-sex --make-bed --out YTC_MergeData

#LD Prunning
plink --bfile YTC_MergeData --indep-pairwise 50 10 0.5
plink --bfile YTC_MergeData --extract plink.prune.in --make-bed --out YTCMergeData_LD


