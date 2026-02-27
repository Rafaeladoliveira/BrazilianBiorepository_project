#Load library "RColorBrewer"
if(!require(RColorBrewer)) install.packages("RColorBrewer", repos='http://cran.us.r-project.org')
library(RColorBrewer)

#Load your admixture output for merged dataset (Ref Panel + Target Cohort) - "YTCMergeData_LD.4.Q"
Merge_data <- read.table("YTCMergeData_LD.4.Q")

#Extract your target cohort samples
Only_YTC <- Merge_data[1:1111, ] #change "1111" and put the number of samples that correspond to your target cohort
colMeans(Only_YTC) #calculate the mean of each ancestry group (column) 

#Reorganize columns (Update these indices based on your manual check!)
#Goal: V1:EUR, V2:AFR, V3:AMR, V4:ASN
#Note: To conclude this step, you need to open the Merge_data table and look specifically to the samples from the reference panel and see which column exlains the ancestry group (have the highest value).
MergeData_Ordered <- Merge_data[, c(1,4,3,2)]
YTC_Ordered <- Only_YTC[, c(1,4,3,2)]
YTC_Ordered[1:1111, ]#Confirm if the means still the same for each ancestry
colMeans(YTC_Ordered) 

#Define colors and labels
my_colors <- brewer.pal(4, "Set1")
pop_labels <- c("EUR", "AFR", "AMR", "ASN")

#Plot: Merged Dataset
pdf("ADMIXTURE_Barplot_K4_YTCMergeData.pdf", width = 12, height = 5)
par(mar = c(1.5, 4, 2.5, 2), cex.lab = 0.75, cex.axis = 0.6)
barplot(t(as.matrix(MergeData_Ordered)),
        col = my_colors,
        ylab = "Ancestry Proportions",
        border = NA,
        space = 0)
legend("bottom", legend = pop_labels, fill = my_colors, horiz = TRUE, cex = 0.8, bty = "n")
dev.off()

#Plot: Target Cohort Only
pdf("ADMIXTURE_Barplot_K4_YTCData.pdf", width = 12, height = 5)
par(mar = c(1.5, 4, 2.5, 2), cex.lab = 0.75, cex.axis = 0.6)
barplot(t(as.matrix(YTC_Ordered)),
        col = my_colors,
        ylab = "Ancestry Proportions",
        border = NA,
        space = 0)
legend("bottom", legend = pop_labels, fill = my_colors, horiz = TRUE, cex = 0.8, bty = "n")
dev.off()

#Create new Ordered .Q file on the working directory
write.table(MergeData_Ordered, file = "YTCMergeData_LDP_ordered.4.Q", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)


    