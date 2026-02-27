
library(dplyr)

pihat_min0.2_in_founders <- read.csv("~/TargetCohort_QC/ytc_pihat_min0.2_in_founders.genome", sep="")

plink <- read.csv("~/TargetCohort_QC/plink.imiss", sep="")

select_fid <- function(fid1, fid2, iid1, iid2) {
  f_miss_fid1 <- plink %>%
    filter(FID == fid1 & IID == iid1) %>%
    select(F_MISS) %>%
    unlist()
  
  f_miss_fid2 <- plink %>%
    filter(FID == fid2 & IID == iid2) %>%
    select(F_MISS) %>%
    unlist()
  
  if (f_miss_fid1 > f_miss_fid2) {
    return(c(fid1, iid1))
  } else {
    return(c(fid2, iid2))
  }
}

# Apply the function to each line of pihat_min0.2_in_founders
result <- t(apply(pihat_min0.2_in_founders, 1, function(row) {
  select_fid(row["FID1"], row["FID2"], row["IID1"], row["IID2"])
}))

#Rename Columns
colnames(result) <- c("FID", "IID") #the file do not need header

# Salvar a variável result em um arquivo de texto
write.table(result, file = "~/TargetCohort_QC/0.2_low_call_rate_pihat.txt", sep = "\t", row.names = FALSE, quote = FALSE)

