#!/bin/bash
#04_Target_Cohort_QC

MAIN_DIR=$(pwd) #Define the home base (project main folder)

OUTDIR="TargetCohort_QC"
mkdir -p $OUTDIR
cd $OUTDIR

#Convert target cohort to plink binary format (.bed/.bim/.fam)
plink --vcf your_target_cohort.vcf --double-id --make-bed --out your_target_cohort

#SNP and Sample-level Quality Control (QC)
#Remove SNPs and Samples with Missingness > 2%  (--geno/--mind 0.02)
plink --bfile your_target_cohort --geno 0.02 --mind 0.02 --make-bed --out your_target_cohort_genomind0.02

#Remove non-autosomal SNPs and variants with MAF < 0.05
 #Generate a file with autosomal SNPs
awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' your_target_cohort_genomind0.02.bim > your_target_cohort_Autosomalchr.txt
 #Extract and keep only the autosomal SNPs
plink --bfile your_target_cohort_genomind0.02 --extract your_target_cohort_Autosomalchr.txt --make-bed --out ytc_Autosomal
 #Exclude SNPs with MAF < 0.05
plink --bfile ytc_Autosomal --maf 0.05 --make-bed --out ytc_maf0.05

#Hardy-Weinberg Equilibrium (HWE)
#IMPORTANT: The HWE is only when controls are present. Therefore, if your target cohort has only samples of cancer cases skip this step.
plink --bfile ytc_Autosomal --hwe 1e-6 --make-bed --out ytc_HWE

#Heterozygosity 
 #Identify SNPs in high-LD (Linkage Desiquilibrium) - get independent SNPs
 #The inversion.txt (w/high inversion regions information) is available on /data directory.
plink --bfile ytc_HWE --exclude ../data/inversion.txt --range --indep-pairwise 50 5 0.2 --out indepSNP 
plink --bfile ytc_HWE --extract indepSNP.prune.in --het --out R_check 

 #Run the R script "heterozygosity_outliers_list.R" (available on /scripts directory)
module load R/4.2.1 #use your R version. This step is important if you are running your analysis on a cluster.
Rscript ../scripts/heterozygosity_outliers_list.R
echo "R analysis complete."

 #Transform the output file (fail-het-qc.txt) to be compatible with PLINK
 sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt

 #Remove samples which heterozygosity rate deviates more than 3SD from the mean
plink --bfile ytc_HWE --remove het_fail_ind.txt --make-bed --out ytc_HET

#Relatedness (IBD)
 #Check for relationships between individuals with a pihat > 0.2
plink --bfile ytc_HET --extract indepSNP.prune.in --genome --min 0.2 --out ytc_min0.2 #indepSNP.prune.in was generated on the previous code:"plink --bfile ytc_HWE --exclude inversion.txt --range --indep-pairwise 50  0.2 --out indepSNP "

 #Filter founders - individuals without parents in the cohort
plink --bfile ytc_HET --filter-founders --make-bed --out ytc_Founders
plink --bfile ytc_Founders --extract indepSNP.prune.in --genome --min 0.2 --out ytc_pihat_min0.2_in_founders
plink --bfile ytc_Founders --missing

 #Run the R script "Relatedness_LCR.R" (available on /scripts directory). This script removes the individual with the lowest call rate for each pair of 'related' individuals with a pihat > 0.2.
module load R/4.2.1 #use your R version. This step is important if you are running your analysis on a cluster.
Rscript ../scripts/Relatedness_LCR.R
echo "R analysis complete."

plink --bfile ytc_Founders --remove 0.2_low_call_rate_pihat.txt --make-bed --out ytc_Final

#Recode from Plink binary format to VCF
plink --bfile ytc_Final --recode vcf --out ytc_Final
bgzip ytc_Final.vcf
tabix -p vcf ytc_Final.vcf.gz

#Liftover (GRCh38 > GRCh37) and Normalization of Target Cohort VCF
 #Download the Reference.FASTA for GRCh37. This file is required for the liftover.
wget ftp://ftp.1000genomes.ebi.ac.uk//vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gunzip hs37d5.fa.gz
samtools faidx hs37d5.fa

 #Download the reference chain file
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz

 #Perform the liftover
CrossMap vcf hg38ToHg19.over.chain.gz ytc_Final.vcf.gz hs37d5.fa ytc_GRCh37.vcf
bgzip ytc_GRCh37.vcf
tabix -p vcf ytc_GRCh37.vcf.gz

 #Sort the lifted VCF file.
bcftools sort -Oz -o ytc_sGRCh37.vcf.gz ytc_GRCh37.vcf.gz
tabix -p vcf ytc_sGRCh37.vcf.gz

 #Perform Normalization via BCFtools
bcftools norm -m -both -f hs37d5.fa -Oz -o ytc_sGRCh37_norm.vcf.gz ytc_sGRCh37.vcf.gz
tabix -p vcf ytc_sGRCh37_norm.vcf.gz

#Unique ID to each Variant.
#The Reference Panel has the variants IDs as unique IDs - CHR:POS:REF:ALT. Therefore, our target cohort need to follow the same steps to be compatible.
bcftools annotate \
  -x ID \
  -I +'%CHROM:%POS:%REF:%ALT' \
  -Oz \
  -o ytc_uniqueID.vcf.gz \
  ytc_sGRCh37_norm.vcf.gz

bcftools index ytc_uniqueID.vcf.gz

#Remove exact duplicated variants
bcftools view -H ytc_uniqueID.vcf.gz | awk '{print $1, $2, $4, $5}' | sort | uniq -d #if this command has no output (meaning zero duplicates), it is not necessary to run the following 2 commands
bcftools norm -d exact -Oz -o ytc_nodups.vcf.gz ytc_uniqueID.vcf.gz
bcftools index ytc_nodups.vcf.gz

#Convert VCF file to Plink Binary format
plink2 \
  --vcf ytc_nodups.vcf.gz \
  --double-id \
  --set-missing-var-ids '@:#:$r:$a' \
  --make-bed \
  --out FinalYTC

cd "$MAIN_DIR"
echo "STEP 04 Complete!"