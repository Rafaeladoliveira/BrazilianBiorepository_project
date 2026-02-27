#!/bin/bash
#03_Merge_Reference_Panel.sh

MAIN_DIR=$(pwd) #Define the home base (project main folder)

#Ensure that you have the files: 1000G_Phase3Merged_maf0.05 (Final output of script 01) and NAMRSamples_maf0.05 (Final output of script 02)
OUTDIR="FinalReferenceDataset"
mkdir -p $OUTDIR
cp 1000G_PLINK/1000G_QCed/1000G_Phase3Merged_maf0.05 FinalReferenceDataset/
cp nAMR_PLINK/nAMR_QCed/NAMRSamples_maf0.05 FinalReferenceDataset/
cp nAMR_PLINK/NativeAMR_Samples.txt FinalReferenceDataset/

cd $OUTDIR

#Set the reference genome
awk '{print$2,$5}' 1000G_Phase3Merged_maf0.05.bim > 1000G_Phase3Merged_maf0.05.txt
plink --bfile NAMRSamples_maf0.05 --reference-allele 1000G_Phase3Merged_maf0.05.txt --make-bed --out NAMRSamples_final-adj

#Merge Datasets
plink --bfile 1000G_Phase3Merged_maf0.05 --bmerge NAMRSamples_final-adj.bed NAMRSamples_final-adj.bim NAMRSamples_final-adj.fam --allow-no-sex --make-bed --out ReferenceDataset

#During merging appeard some warnings for duplicated positions "Warning: Variants 1:15274:A:T and 1:15274:A:G have the same position". 
#Remove duplicated variants
plink --bfile ReferenceDataset --list-duplicate-vars ids-only
if [ -f plink.dupvar ]; then
    plink --bfile ReferenceDataset --exclude plink.dupvar --make-bed --out ReferenceDataset_clean
else
    mv ReferenceDataset.bed ReferenceDataset_clean.bed
    mv ReferenceDataset.bim ReferenceDataset_clean.bim
    mv ReferenceDataset.fam ReferenceDataset_clean.fam
fi

#Remove Admixed American samples - PUR/CLM/MXL/PEL - from 1000G. The 1000GRacefile.txt is available on the /data directory.
awk '$3 == "AMR" {print $1, $1}' ../data/1000GRacefile.txt > admAMR_samples.txt
plink --bfile ReferenceDataset_clean --remove admAMR_samples.txt --make-bed --out ReferenceDatasetFinal

#Create a Final RaceFile (Population Map) for the Reference Panel - with Sample ID and Superpopulation code (IMPORTANT: NO HEADERS)
awk 'NR==FNR{remove[$1]; next} !($1 in remove)' admAMR_samples.txt 1000GRacefile.txt > New1000GRacefile.txt
tail -n +2 New1000GRacefile.txt | cut -f1,3 > New1000GRacefileFinal.txt 
awk '{print $1, $4}' NativeAMR_Samples.txt  > NativeAMR_SamplesFinal.txt  
sed -i 's/America/AMR/g' NativeAMR_SamplesFinal.txt
sed -i 's/EAS/ASN/g' New1000GRacefileFinal.txt
sed -i 's/SAS/ASN/g' New1000GRacefileFinal.txt
cat New1000GRacefileFinal.txt NativeAMR_SamplesFinal.txt > ReferenceDataset_Racefile.txt
awk '{print$1}' ReferenceDatasetFinal.fam > ReferenceDatasetFinal.txt
awk 'NR==FNR {pop[$1]=$2; next} {print $1, pop[$1]}' ReferenceDataset_Racefile.txt ReferenceDatasetFinal.txt > Correct_ReferenceRaceFile.txt

#Note: ReferenceDatasetFinal vcf has EAS, EUR and AFR samples from the 1000Genomes Project Dataset; and native AMR samples from the HGDP Dataset. In total, in this vcf we have 1712 samples.

#Return to Root
cd "$MAIN_DIR"
echo "Merge Complete. Final Panel: $OUTDIR/ReferenceDataset_Clean"