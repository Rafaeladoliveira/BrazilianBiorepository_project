#!/bin/bash
#01_reference_panel_prep.sh

MAIN_DIR=$(pwd) #Define the home base (project main folder)

#Download the 1000Genomes Phase 3 data for each chromossome. This data is on GRCh37 genomic build.

OUTDIR="1000G_PLINK"
mkdir -p $OUTDIR
cd $OUTDIR

echo "Starting download of 1000 Genomes Phase 3 (GRCh37)..."

for CHR in {1..22}; do
  FILE="ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
  URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/${FILE}"
  wget -c $URL
done

echo "Download complete. Files are located in $OUTDIR"

#Verify the Reference.FASTA used on 1000G data and download it. This file is required by bcftools norm to correctly check and align the REF/ALT alleles in the variants.

zless ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz | head -n 50
##reference= ftp://ftp.1000genomes.ebi.ac.uk//vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz

wget ftp://ftp.1000genomes.ebi.ac.uk//vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gunzip hs37d5.fa.gz
samtools faidx hs37d5.fa #creates an index file for the hs37d5.fa.

#Perform VCFs Normalization using BCFtools
REF="hs37d5.fa"

echo "Starting normalization and indexing..."
for CHR in {1..22}; do
    echo "--- Processing Chromosome: ${CHR} ---"
    
    IN="ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    OUT="1000G_CHR${CHR}.norm.vcf.gz"
    bcftools norm -m -both -f $REF -Oz -o $OUT $IN
    tabix -p vcf $OUT
done
echo "Normalization complete."

#Unique ID to each Variant.
#The 1000 Genomes Project dataset contains variants lacking standard rs-identifiers, denoted by a period ("."). Because PLINK (and other tools) may treat these missing identifiers as duplicate entries for the same variant.
#Therefore, it is important to assign unique identifiers based on genomic coordinates (%CHROM:%POS:%REF:%ALT) to ensure each SNP is uniquely and accurately represented.

mkdir -p 1000Genomes_UniqueID 

for CHR in {1..22}; do
    echo "--- Annotating Chromosome: ${CHR} ---"
    
    IN="1000G_CHR${CHR}.norm.vcf.gz"
    OUT="1000Genomes_UniqueID/chr${CHR}.unique.vcf.gz"
    bcftools annotate \
        -x ID \
        -I +'%CHROM:%POS:%REF:%ALT' \
        -Oz \
        -o $OUT \
        $IN
    bcftools index $OUT
done

#Remove exact duplicated variants
for CHR in {1..22}; do
    IN="1000Genomes_UniqueID/chr${CHR}.unique.vcf.gz"
    OUT="1000Genomes_UniqueID/chr${CHR}.nodups.vcf.gz"
    bcftools norm \
    -d exact \
    -Oz \
    -o $OUT \
    $IN
     bcftools index $OUT
done

#Convert all 1000G VCFs to PLINK Binary Format (.bed/.bim/.fam)
mkdir -p 1000G_Chr_BF

for CHR in {1..22}; do
    IN="1000Genomes_UniqueID/chr${CHR}.nodups.vcf.gz"
    OUT="1000G_Chr_BF/finalchr${CHR}"
    plink2 \
  --vcf $IN \
  --double-id \
  --set-missing-var-ids '@:#:$r:$a' \
  --make-bed \
  --out $OUT
done

#Merge all 1000G CHR binary files 
mkdir -p 1000G_FinalMerged

for CHR in {2..22}; do
  echo " 1000G_Chr_BF/finalchr${CHR}" >> merge_list.txt
done

plink --bfile 1000G_Chr_BF/finalchr1 --merge-list merge_list.txt --make-bed --out 1000G_FinalMerged/1000G_Phase3_merged

#SNP and Sample-level Quality Control (QC)

#Remove SNPs and samples with missingness > 2%  (--geno/--mind 0.02), and SNPs with MAF < 0.05
mkdir -p 1000G_QCed

plink --bfile 1000G_FinalMerged/1000G_Phase3_merged \
  --geno 0.02 \
  --mind 0.02 \
  --maf 0.05 \
  --allow-no-sex \
  --make-bed \
  --out 1000G_QCed/1000G_Phase3Merged_maf0.05

cd "$MAIN_DIR" #Return to the project main root (home base)
echo "Process complete. Current directory: $(pwd)"