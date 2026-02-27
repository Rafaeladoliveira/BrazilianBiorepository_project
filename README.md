# ADMIXTURE Ancestry Analysis Pipeline
-----------------------------

**Project Name:** Brazilian Biorepository Ancestry Inference 

**Description:** This repository provides the computational pipeline for the global ancestry estimation of the Brazilian Biorepository cohort. This code supports the findings presented in the manuscript, specifically focusing on the population structure of individuals diagnosed with colorectal, breast, and cervical cancer.

**Citation:** Novaes, L. A. C., Oliveira, R. D., Tegami, I. L., Gonçalves, M. F. S., Ribeiro Junior, H. L., Reis, M. B., Moreno, D. A., Possati-Resende, J., Santos, F., Hashimoto, C., Antoniazzi, A., Baraldo, S., Romagnolo, L., dos Reis, R., da Silva, L. S., Leal, L. F., Guimarães, D. P., Marques, M. M. C., Ribeiro, A. G., and Reis, R. M. (2026). Brazilian biorepository to support genome-wide association studies of colorectal, breast, and cervical cancer.

* * * * *

### Required Softwares

| **Software** | **Version** | **Link/Citation** |
| --- | --- | --- |
| plink | 1.90b6.21 | [PLINK 1.9](https://www.cog-genomics.org/plink/) |
| plink2 | 2.00a3LM | [PLINK 2.0](https://www.cog-genomics.org/plink/2.0/) |
| R | 4.3.0 or above | [R Project](https://www.r-project.org/) |
| ADMIXTURE | 1.3.1 | [ADMIXTURE](https://dalexander.github.io/admixture/) |
| BCFtools | 1.20 or above | [BCFtools](https://samtools.github.io/bcftools/) |
| CrossMap | 0.6.0 or above | [CrossMap](https://pythonhosted.org/CrossMap/) |


### Required R Packages

| **R Package** | **Version** | **Link/Citation** |
| --- | --- | --- |
| tidyverse | 2.0.0 | [tidyverse](https://www.tidyverse.org/) |
| ggplot2 | 3.1.3.1 | [gplots](https://ggplot2.tidyverse.org/) |
| data.table | 1.14.10 | [data.table](https://github.com/Rdatatable/data.table) |

### Languages

| **Language** | **Function** | 
| --- | --- | 
| Bash| Pipeline Automatation & Data Parsing | 
| R | Data Visualization | 

### Data Availability
> [!NOTE]
> Due to the sensitive nature of human genomic information and the terms of the informed consent, the raw genotyping data for the samples analyzed in this manuscript are not publicly available within this repository.

To test this pipeline use your own genotyped data, replacing the placeholder filename **your_genodata.vcf**. All code was constructed for data in **GRCh38**. However, to ensure compatibility with the 1000 Genomes Phase 3 reference panel used for the ancestry analysis, our cohort was subsequently **lifted over to GRCh37** (see script X).


### Pipeline Workflow
* * * * *
To reproduce the results and figures presented in the manuscript, all scripts are provided in the /scripts directory. The scripts are numbered sequentially (01 to X) to indicate the required order of execution.

#### 1) Create a Reference Panel Suitable for Brazilian Cohorts

This initial step focuses on building a high-fidelity reference panel tailored for populations with complex admixture, such as Brazilian cohorts. By combining the 1000 Genomes Project with refined HGDP Native American lineages, we ensure a robust baseline for global ancestry inference. **This step was refined according to information described on the methodology described in the Materials and Methods section of this study: https://doi.org/10.3390/diagnostics15091098, in which ancestry evaluation was also conducted on a brazilian cohort**. 

List of scripts and their functions for **STEP1**:

| **Script** | **Function** |
| --- | --- |
| 01_1000G_panel_prep.sh | Downloads and standardizes the 1000G Phase 3 reference data |
| 02_HGDP_nAMR_prep.sh |  Extracts and prepares Native American samples from the HGDP dataset |
| 03_Merge_Reference_Panel.sh | Merges 1000G and HGDP-nAMR into a final unified Reference Panel |

 On the ***01_1000G_panel_prep.sh***, 1000Genomes Phase 3 data was downloaded, normalized, converted to PLINK Binary Format (.bed/.bim/.fam), and lastly, filtered using a missingness threshold of 2% for both samples and variants, and also a cut-off of 5% for variants. This QC protocol was applied to ensure high-fidelity ancestry inference.

To run ***01_1000G_panel_prep.sh*** script:

```bash
bash scripts/01_1000G_panel_prep.sh
```

> [!NOTE]
> **QC Strategy:** A missingness threshold of 2% (--geno and --mind) was used to exclude low-quality samples and variants. Additionally, a Minor Allele Frequency (MAF) cut-off of 5% (--maf 0.05) was applied to prioritize common variants that provide a stable signal for global ancestry inference.

 After preparing the 1000G Reference data, Whole genome Sequencing (WGS) data from HGDP was processed to extract high-quality Native American (nAMR) reference panel - important for Brazilian ancestry analysis. The script ***02_HGDP_nAMR_prep.sh*** performed nAMR samples selection, conversion of genomic coordinates from **GRCh38 to GRCh37 (hg19) via Liftover**, and variant normalization. Additionally, an advanced error-resolution protocol was implemented to solve inconsistencies identified during merging attempts. Lastly, the data was converted to plink binary format and filtered using the same QC cut-offs previously applied to 1000G Data.

 To run ***02_HGDP_nAMR_prep.sh*** script:

```bash
bash scripts/02_HGDP_nAMR_prep.sh
```
Lastly, using the script ***03_Merge_Reference_Panel.sh***, the 1000G QCed panel and the nAMR QCed panel (from HGDP) were merged to generate the final reference dataset. During this process, admixed American individuals from the 1000G data were excluded to ensure the Native American component remained unadmixed. This resulted in a unified, high-fidelity reference panel specifically appropriated for global ancestry analysis in Brazilian cohorts.

 To run ***03_Merge_Reference_Panel.sh*** script:

```bash
bash scripts/03_Merge_Reference_Panel.sh
```

> [!NOTE]
>**STEP 1 is a one-time procedure. The resulting unified Reference Panel is standardized and applied across all neoplasia cohorts described in this study, including the breast cancer, colorectal cancer, cervical cancer, and control cohorts.**

#### 2) Target Cohort Processing and Reference Integration

Following the creation of the Reference Panel (STEP1), the study cohorts (target samples) must undergo rigorous standardization and QC. **STEP 2 describes the pipeline used to refine the Brazilian cohorts and integrate them with the unified Reference Panel**. This process ensures that only high-quality, overlapping variants are used for the final ancestry estimation, maintaining consistency across the control dataset and also breast, colorectal, and cervical cancer datasets.

List of scripts and their functions for **STEP2**:

| **Script** | **Function** |
| --- | --- |
| 04_Target_Cohort_QC.sh | Performs Quality Control on the specific study cohort |
| 05_Target_RefPanel_Merge.sh |  Merges the QCed Target cohort with the Reference Panel and performs LD pruning|

The target cohorts were prepared for compatibility with the Reference Panel using the script ***04_Target_Cohort_QC.sh***. This script handles sample and variant filtering, relatedness checks, and coordinate liftover. To execute this script:

```bash
bash scripts/04_Target_Cohort_QC.sh
```

The final stage before ancestry estimation is the integration of the target cohort with the Reference Panel, followed by Linkage Disequilibrium (LD) pruning to ensure variant independence for ADMIXTURE. To run the ***05_Target_RefPanel_Merge.sh*** script:

```bash
bash scripts/05_Target_RefPanel_Merge.sh
```

> [!NOTE]
>**STEP 2 was repeated for each target Brazilian cohort described in the article, namely the controls cohort; and breast, colorectal and cervical cancer cohorts.**

#### 3) Global Ancestry Estimation using ADMIXTURE and Visualization using R

In the final stage, ADMIXTURE (v1.3.0) is utilized to estimate individual ancestry proportions for each of the study cohorts. We employ a supervised approach using the previously defined reference populations (AFR, EUR, ASN, and AMR(=Native American) to decompose the genomic background of the Brazilian samples into four ancestral components.

> [!IMPORTANT]
>**Supervised analysis ensures that the reference populations act as fixed anchors, providing more consistent ancestry estimates across different target cohorts.**


| **Script** | **Function** |
| --- | --- |
| 06_Admixture_Analysis.sh | Runs ADMIXTURE for K = 4. |

The results are processed in R to generate publication-quality ancestry barplots and to extract the mean ancestry proportions for each cancer cohort.

To run ***06_Admixture_Analysis.sh*** script:

```bash
bash scripts/06_Admixture_Analysis.sh
```

> [!NOTE]
>**STEP 3 was repeated for each target Brazilian cohort described in the article, namely the controls cohort; and breast, colorectal and cervical cancer cohorts.**




* * * * *
### Sources
**Some of the pipelines modified from:**



Marees, A. T., de Kluiver, H., Stringer, S., Vorspan, F., Curis, E., Marie‐Claire, C., & Derks, E. M. (2018). A tutorial on conducting genome‐wide association studies: Quality control and statistical analysis. International journal of methods in psychiatric research, 27(2), e1608. https://doi.org/10.1002/mpr.1608

D.H. Alexander, J. Novembre, and K. Lange. (2009). Fast model-based estimation of ancestry in unrelated individuals. Genome Research, 19:1655–1664. https://doi.org/10.1101/gr.094052.109

Liu, C.-C., Shringarpure, S., Lange, K., and Novembre, J. (2020). Exploring population structure with admixture models and principal component analysis. In J. Noyvert (Ed.), Statistical Population Genomics (pp. 67–86). Humana, New York, NY. https://doi.org/10.1007/978-1-0716-0199-0_4



* * * * *
### Final Remarks

> [!NOTE]
> This repository is created to enhance reproducibility. It will not be actively maintained.

> [!TIP]
> If you need assistance at any step, please feel free to [email](rafaeladiasoliveira98@gmail.com) me.

> [!IMPORTANT]
>  It is not possible to share the genotyped data used in this pipeline. Run this pipeline using your own data.

> [!CAUTION]
> If you use any part of this pipeline in your work, please make sure to cite the original work software/package and if possible our work.