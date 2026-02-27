#!/bin/bash
#02_HGDP_nAMR_prep.sh

MAIN_DIR=$(pwd) #Define the home base (project main folder)

#Filter the Native American (nAMR) samples from HGDPSamples.txt (available on /data directory).
#nAMR samples include samples from Maza, Pima, Surui, and Karitiana. 
#The Colombian samples on HGDP correspond to admixed Americans.
OUTDIR="nAMR_PLINK"
mkdir -p $OUTDIR
cd $OUTDIR

grep "America" ../data/HGDPSamples.txt > AMR_HGDP.txt
grep "Colombian" AMR_HGDP.txt > NonNative_Toexclude.txt
grep -v -F -f NonNative_Toexclude.txt AMR_HGDP.txt | sort > NativeAMR_Samples.txt

#Download the VCFs files of HGDP nAMR samples
#Removed 01006 and 01051 as they are unavailable on the Sanger FTP
SAMPLES=(00832 00837 00838 00843 00845 00849 00854 00856 00858 00859 \  #Definition of our samples list that will be used all the times that we run "${SAMPLES[@]}"
00860 00861 00862 00863 00864 00865 00868 00869 00870 00871 00872 00873 \
00875 00876 00877 00995 00999 01001 01009 01010 01013 01014 01019 01037 \
01041 01043 01050 01053 01055 01056 01057 01058 01059 01060)

for SAMPLE in "${SAMPLES[@]}"; do
  FILE="HGDP${SAMPLE}.hgdp_wgs.20190516.vcf.gz"
    URL="https://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/gVCFs/${FILE}"
  wget -c "$URL"
done 

#Download the index (.tbi) files of HGDP nAMR samples
for SAMPLE in "${SAMPLES[@]}"; do
  FILE="HGDP${SAMPLE}.hgdp_wgs.20190516.vcf.gz.tbi"
    URL="https://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/gVCFs/${FILE}"
  wget -c "$URL"
done

#Download the Reference.FASTA for GRCh37. This file is required for the liftover.
wget ftp://ftp.1000genomes.ebi.ac.uk//vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gunzip hs37d5.fa.gz
samtools faidx hs37d5.fa

#Download Chain file for the coordinate conversion during liftover.
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz

#Perform Liftover (GRCh38 > GRCh37) on nAMR samples data
mkdir -p nAMR_GRCh37

for SAMPLE in "${SAMPLES[@]}"; do
  FILE="HGDP${SAMPLE}.hgdp_wgs.20190516.vcf.gz"
  OUT="nAMR_GRCh37/hgdp${SAMPLE}_GRCh37.vcf"
  CrossMap vcf hg38ToHg19.over.chain.gz "$FILE" hs37d5.fa "$OUT"
  bgzip "$OUT"
  tabix -p vcf "${OUT}.gz"
  rm -f "${OUT}.unmap"
done

#Reheader GRCh37 VCF Files due to inconsistency between chr name in variants (chr1) and in the header (##contig=<ID=1)
cat > nAMR_GRCh37/Rename_chr.txt <<EOF  #this command creates the mapping file
chr1    1
chr2    2
chr3    3
chr4    4
chr5    5
chr6    6
chr7    7
chr8    8
chr9    9
chr10   10
chr11   11
chr12   12
chr13   13
chr14   14
chr15   15
chr16   16
chr17   17
chr18   18
chr19   19
chr20   20
chr21   21
chr22   22
chrX    X
chrY    Y
chrM    MT
EOF

for SAMPLE in "${SAMPLES[@]}"; do
  MAP="nAMR_GRCh37/Rename_chr.txt"
  IN="nAMR_GRCh37/hgdp${SAMPLE}_GRCh37.vcf.gz"
  OUT="nAMR_GRCh37/hgdp${SAMPLE}_GRCh37_nochr.vcf.gz"
  bcftools annotate --rename-chrs "$MAP" -Oz -o "$OUT" "$IN"
done

#Sort and Index GRCh37 nAMR VCF Files
for SAMPLE in "${SAMPLES[@]}"; do
  IN="nAMR_GRCh37/hgdp${SAMPLE}_GRCh37_nochr.vcf.gz"
  OUT="nAMR_GRCh37/hgdp${SAMPLE}_sGRCh37.vcf.gz"
  bcftools sort -Oz -o "$OUT" "$IN"
  tabix -p vcf "$OUT"
done

#Filter out non-standard REF alleles BEFORE normalizing and Perform VCFs Normalization using BCFtools
mkdir -p nAMR_Norm

for SAMPLE in "${SAMPLES[@]}"; do
  IN="nAMR_GRCh37/hgdp${SAMPLE}_sGRCh37.vcf.gz"
  OUT="nAMR_Norm/HGDP${SAMPLE}_sGRCh37.norm.vcf.gz"
  bcftools view -e 'REF!="A" & REF!="C" & REF!="G" & REF!="T" & REF!="N"' "$IN" | \
  bcftools norm -m -both -f hs37d5.fa  -Oz -o "$OUT"
  tabix -p vcf "$OUT"
done

#Unique ID Assingnment to each Variant.
#The 1000 Genomes Project dataset contains with unique identifiers based on genomic coordinates (%CHROM:%POS:%REF:%ALT). Therefore, HGDP needs to have variants ID in the same format to be compatible.
mkdir -p nAMR_UniqueID
for SAMPLE in "${SAMPLES[@]}"; do
  IN="nAMR_Norm/HGDP${SAMPLE}_sGRCh37.norm.vcf.gz"
  OUT="nAMR_UniqueID/HGDP${SAMPLE}_uniqueID.vcf.gz"
  bcftools annotate \
  -x ID \
  -I +'%CHROM:%POS:%REF:%ALT' \
  -Oz \
  -o "$OUT" \
"$IN"
  tabix -p vcf "$OUT"
done

#Remove non-variants positions (with ALT=".") from HGDP VCFs Files - not compatible with 1000G Data.
for SAMPLE in "${SAMPLES[@]}"; do
  IN="nAMR_UniqueID/HGDP${SAMPLE}_uniqueID.vcf.gz"
  OUT="nAMR_UniqueID/HGDP${SAMPLE}_variantsonly.vcf.gz"
  bcftools view -e  'ALT="."' -Oz -o "$OUT" "$IN"
  tabix -p vcf "$OUT" 
done

#Remove exact duplicated variants
for SAMPLE in "${SAMPLES[@]}"; do
  IN="nAMR_UniqueID/HGDP${SAMPLE}_variantsonly.vcf.gz"
  OUT="nAMR_UniqueID/HGDP${SAMPLE}_nodups.vcf.gz"
  bcftools norm -d exact -Oz -o "$OUT" "$IN"
  bcftools index "$OUT"
done

#At this phase, trying to merge all files output several errors that will be adressed one by one on the following commands.
#Error1: Incorrect number of FORMAT/PL values at [chr:pos], cannot merge. The tag is defined as Number=G, but found X values and Y alleles.
mkdir -p nAMR_Final
for SAMPLE in "${SAMPLES[@]}"; do
  IN="nAMR_UniqueID/HGDP${SAMPLE}_nodups.vcf.gz"
  OUT="nAMR_Final/HGDP${SAMPLE}.clean.vcf.gz"
  bcftools annotate -x FORMAT/PL,FORMAT/AD,FORMAT/PGT,FORMAT/PID -Oz -o "$OUT" "$IN"
  tabix -p vcf "$OUT"
done 

#Error2: Several errors like this [E::bcf_calc_ac] Incorrect allele ("31431904") in HGDP00865 at 1:54712
for file in nAMR_Final/*.clean.vcf.gz; do
  sample=$(bcftools query -l "$file")
  echo "Processing $sample ($file)"
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$file" | \
    awk -v s="$sample" '{print $1"\t"$2"\t"$3"\t"$4"\t"s}' >> nAMR_Final/all_variants.tsv
done
  #Run AWK logic to find mismatches
awk '{
  key = $1"\t"$2
  var = $3":"$4
  combo[key][var]++
  sample[key][var] = sample[key][var]","$5
}
END {
  for (pos in combo) {
    if (length(combo[pos]) > 1) {
      print "[WARNING] Inconsistent position: " pos
      for (v in combo[pos]) {
        print "\t" v " → Samples:" sample[pos][v]
      }
    }
  }
}' nAMR_Final/all_variants.tsv > nAMR_Final/inconsistent_positions.txt

  #Extract coordinates
awk '/^\[WARNING\] Inconsistent position:/ {split($0,a,"\t"); print a[1], a[2]}' nAMR_Final/inconsistent_positions.txt | awk '{print $4 "\t" $5}' | sort -u > nAMR_Final/bad_positions.txt

  #Create the final exclusion list
while read chr pos; do
  echo -e "$chr\t$pos"
done < nAMR_Final/bad_positions.txt > nAMR_Final/regions_to_exclude.txt

awk -v OFS="\t" '{print $1, $2-1, $2}' nAMR_Final/regions_to_exclude.txt > nAMR_Final/regions_to_exclude.bed #Important: Convert to BED Format for BCFtools (chr, start, end)

  #Fix Error2 by filtering out inconsistent positions
for SAMPLE in "${SAMPLES[@]}"; do
  IN="nAMR_Final/HGDP${SAMPLE}.clean.vcf.gz"
  OUT="nAMR_Final/HGDP${SAMPLE}.final2.vcf.gz"

  echo "Filtering inconsistencies for Sample: HGDP${SAMPLE}"

  bcftools view -T ^nAMR_Final/regions_to_exclude.bed "$IN" -Oz -o "$OUT"
  tabix -f "$OUT"
done

#Error3: Keep appearing the inconsistencies errors - Solve GT issues

echo -e "File\tCHROM\tPOS\tREF\tALT\tGT" > nAMR_Final/all_gt_problems.tsv

for file in nAMR_Final/*.final2.vcf.gz; do
  sample=$(bcftools query -l "$file")
  echo "Checking $file ($sample)..."

  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' "$file" | \
  awk -v f="$file" '
  {
    split($4, alt, ",")
    n_alt = length(alt)
    gsub(/[\/|]/, " ", $5)
    split($5, gt, " ")
    for (i in gt) {
      if (gt[i] != "." && gt[i] > n_alt) {
        print f"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5
        break
      }
    }
  }' >> nAMR_Final/all_gt_problems.tsv
done
echo "GT consistency check complete. Results in nAMR_Final/all_gt_problems.tsv"

  #Extract unique problematic positions per file for downstream filtering
cut -f1,2,3 nAMR_Final/all_gt_problems.tsv | tail -n +2 | sort | uniq > nAMR_Final/all_problems_cleaned.tsv

  #To each affected file:
mkdir -p nAMR_FINAL2
for file in $(cut -f1 nAMR_Final/all_problems_cleaned.tsv | sort | uniq); do
  
  echo "Filtering $file..."
  #Extract only this file’s positions:
  awk -v f="$file" '$1 == f { print $2"\t"$3 }' nAMR_Final/all_problems_cleaned.tsv > nAMR_Final/tmp_bad_positions.txt
  awk -v OFS="\t" '{print $1, $2-1, $2}' nAMR_Final/tmp_bad_positions.txt > nAMR_Final/tmp_bad_positions.bed #Important: Convert to BED Format for BCFtools (chr, start, end)

  #Apply the filter with bcftools to remove the inconsistent positions regarding GP:
  outname=$(basename "$file")
  bcftools view -T ^nAMR_Final/tmp_bad_positions.bed "$file" -Oz -o "nAMR_FINAL2/$outname"
  tabix -f "nAMR_FINAL2/$outname"
done

rm -f nAMR_Final/tmp_bad_positions.txt

  #Copy samples that had no errors into the final folder (nAMR_FINAL2)
for original_file in nAMR_Final/*.final2.vcf.gz; do
    base=$(basename "$original_file")
    if [ ! -f "nAMR_FINAL2/$base" ]; then
        echo "Copying healthy sample: $base"
        cp "$original_file" "nAMR_FINAL2/$base"
        tabix -p vcf "nAMR_FINAL2/$base"
    fi
done

#Merge all nAMR samples files 
# Generate a list of all final VCF files to be merged (only successfully processed samples present in nAMR_FINAL2 will be included; missing samples will be skipped)
mkdir -p nAMR_Merged
rm -f nAMR_Merged/mergelist.txt

for SAMPLE in "${SAMPLES[@]}"; do
  echo "nAMR_FINAL2/HGDP${SAMPLE}.final2.vcf.gz" >> nAMR_Merged/mergelist.txt
done

bcftools merge -l nAMR_Merged/mergelist.txt -Oz -o nAMR_Merged/HGDP_nAMR_Merge.vcf.gz
tabix -p vcf nAMR_Merged/HGDP_nAMR_Merge.vcf.gz

#Convert all VCF file to PLINK Binary Format (.bed/.bim/.fam)
plink --vcf nAMR_Merged/HGDP_nAMR_Merge.vcf.gz --make-bed --out nAMR_Merged/NAMRSamples

#SNP and Sample-level Quality Control (QC)
 #Remove SNPs and samples with missingness > 2% (--geno/--mind 0.02), and SNPs with MAF < 0.05
mkdir -p nAMR_QCed

plink --bfile nAMR_Merged/NAMRSamples \
  --geno 0.02 \
  --mind 0.02 \
  --maf 0.05 \
  --allow-no-sex \
  --keep-allele-order \
  --make-bed \
  --out nAMR_QCed/NAMRSamples_maf0.05

cd "$MAIN_DIR" #Return to the project main root (home base)
echo "Process complete. Current directory: $(pwd)"
