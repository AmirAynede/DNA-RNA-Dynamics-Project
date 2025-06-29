###### STEP 0: BASIC PREPARATION ######
# 1. Create the project structure by running "00_setup_project_structure.R".
# 2. Set the working directory to the location of "DRD_2025_project".

# Install and load the necessary packages:
library(minfi)
library(minfiData)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(shinyMethyl)
library(AnnotationDbi)
library(sva)
library(gplots)
library(ggplot2)
library(qqman)
install.packages("lib/SummarizedExperiment_1.38.0.tar.gz", repos = NULL)

# STEP 1: DATA IMPORT AND OBJECT CREATION
# List files in the raw data folder
list.files("./data/raw/")

# Load the sample sheet (must have 'array' and 'slide' columns)
SampleSheet <- read.csv("./data/raw/SampleSheet_Report_II.csv", header = TRUE)
SampleSheet

# Load sample sheet using minfi's function
baseDir <- ("./data/raw")
targets <- read.metharray.sheet(baseDir)
targets

# Create RGChannelSet object
RGset <- read.metharray.exp(targets = targets)
save(RGset, file = "RGset_Report.RData")
load("RGset_Report.RData")

# Explore RGset object
RGset
str(RGset)

###### STEP 2: EXTRACT RED AND GREEN CHANNELS ######
# Extract fluorescence intensity information for Red and Green channels
Red <- data.frame(getRed(RGset))
dim(Red)
head(Red)
Green <- data.frame(getGreen(RGset))
dim(Green)
head(Green)


###### STEP 3: PROBE INFORMATION AND ADDRESS CHECKS ######
# Explore manifest and probe information
getManifest(RGset)
getManifestInfo(RGset)
getProbeInfo(RGset)
ProbeInfo_I <- data.frame(getProbeInfo(RGset))
dim(ProbeInfo_I)
head(getProbeInfo(RGset, type = "II"))
ProbeInfo_II <- data.frame(getProbeInfo(RGset, type = "II"))
dim(ProbeInfo_II)

# Check probe with address 18756452
ProbeInfo_I[ProbeInfo_I$AddressA == "18756452",]
ProbeInfo_I[ProbeInfo_I$AddressB == "18756452",]
ProbeInfo_II[ProbeInfo_II$AddressB == "18756452",]
ProbeInfo_II[ProbeInfo_II$AddressA == "18756452",]

#There is a Type II probe at address 18756452, identified by the name cg01523029.

Probe_Info_18756452 <- ProbeInfo_II[ProbeInfo_II$AddressA == "18756452",]

# Extract Red/Green fluorescence for address 18756452 and create a data frame to present it
Red_fluorescences <- Red[rownames(Red) == "18756452",]
Green_fluorescences <- Green[rownames(Green) == "18756452",]

# Use sample names from the column names of Red/Green (same order)
sample_names <- colnames(Red_fluorescences)

# Create summary dataframe
df_summary_probes <- data.frame(
  "Sample" = sample_names,
  "Red Fluorescence" = as.numeric (Red_fluorescences[1, ]),
  "Green Fluorescence" = as.numeric (Green_fluorescences[1,]),
  "Type" = rep("II", length(sample_names)),
  "Color" = rep("Both", length(sample_names))
)
print(df_summary_probes)

###### STEP 4: CREATE MSet.raw OBJECT ######
# Extract methylated and unmethylated signals
MSet.raw <- preprocessRaw(RGset)
MSet.raw

# Note that now the number of rows is 485512, exactly the number of probes according to the Manifest

save(MSet.raw, file = "./data/processed/MSet_raw.RData")
Meth <- as.matrix(getMeth(MSet.raw))
str(Meth)
head(Meth)
Unmeth <- as.matrix(getUnmeth(MSet.raw))
str(Unmeth)
head(Unmeth)

# Check probe cg01523029 in MethylSet
Unmeth[rownames(Unmeth) == "cg01523029",]
Meth[rownames(Meth) == "cg01523029",]

load("data/raw/Illumina450Manifest.RData")
Illumina450Manifest_clean <- Illumina450Manifest[Illumina450Manifest$CHR != "",]
Illumina450Manifest_clean <- droplevels(Illumina450Manifest_clean)
Probe_Info_cg01523029 <- Illumina450Manifest_clean[Illumina450Manifest_clean$IlmnID == "cg01523029",]

###### STEP 5: QUALITY CONTROL ######
# QC plot
pdf("QCplot.pdf")
qc <- getQC(MSet.raw)
qc
plotQC(qc)
dev.off()         

# Check negative controls
pdf("ControlStripPlot.pdf")
NegativeControl <- data.frame(getProbeInfo(RGset, type = "Control"))
table(NegativeControl$Type)
head(NegativeControl) 
controlStripPlot(RGset, controls = "NEGATIVE")
dev.off()

# Calculate detection p-values
Detection_pvalue <- detectionP(RGset)
save(Detection_pvalue, file = "Detection_pvalue.RData")
Over_Threshold <- Detection_pvalue > 0.05
head(Over_Threshold)
table(Over_Threshold)
sum(Over_Threshold)
sprintf("There are %s failed probes with a detection p-value higher than 0.05", sum(Over_Threshold))
Over_Threshold_Per_Sample <- colSums(Over_Threshold)
Over_Threshold_Per_Sample
summary(Over_Threshold)
failed <- data.frame(
  Sample = names(Over_Threshold_Per_Sample),
  n_Failed_Positions = Over_Threshold_Per_Sample
)
print(failed)

###### STEP 6: CALCULATE BETA AND M VALUES, PLOT DENSITIES ######
# Data preparation
beta <- getBeta(MSet.raw)
M <- getM(MSet.raw)
beta_df <- data.frame(beta)
M_df <- data.frame(M)

pheno <- read.csv("data/raw/SampleSheet_Report_II.csv", header = TRUE, stringsAsFactors = TRUE)

SampleSheet$Group
beta_ctrl <- beta_df[SampleSheet$Group == "CTRL",]
beta_dis <- beta_df[SampleSheet$Group == "DIS",]
M_ctrl <- M_df[SampleSheet$Group == "CTRL",]
M_dis <- M_df[SampleSheet$Group == "DIS",]

mean_of_beta_ctrl <- apply(beta_ctrl, 1, mean, na.rm = TRUE)
mean_of_beta_dis <- apply(beta_dis, 1, mean, na.rm = TRUE)
mean_of_M_ctrl <- apply(M_ctrl, 1, mean, na.rm = TRUE)
mean_of_M_dis <- apply(M_dis, 1, mean, na.rm = TRUE)

# Plot density of mean methylation values
pdf("Beta_M_Density_Plots.pdf")
plot_density <- function() {
  par(mfrow = c(1, 2), mar = c(5, 4, 4, 2))  
  plot(density(mean_of_beta_ctrl, na.rm = TRUE), main = "Density of Beta Values", col = "black", lwd = 2)
  lines(density(mean_of_beta_dis, na.rm = TRUE), col = "orange", lwd = 2)
  legend("topright", legend = c("CTL", "DIS"),
         col = c("black", "orange"), lwd = 2,
         cex = 0.9, bty = "n", inset = c(0, 0))
  plot(density(mean_of_M_ctrl, na.rm = TRUE), main = "Density of M Values", col = "black", lwd = 2)
  lines(density(mean_of_M_dis, na.rm = TRUE), col = "orange", lwd = 2)
  legend("topright", legend = c("CTL", "DIS"),
         col = c("black", "orange"), lwd = 2,
         cex = 0.9, bty = "n", inset = c(0, 0))
}
plot_density()
dev.off()

###### STEP 7: NORMALIZATION WITH preprocessNoob ######
# Subset beta values by probe type
dfI <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type == "I",]
dfI <- droplevels(dfI)
dfII <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type == "II",]
dfII <- droplevels(dfII)

beta_I <- beta[rownames(beta) %in% dfI$IlmnID,]
beta_II <- beta[rownames(beta) %in% dfII$IlmnID,]

mean_of_beta_I <- apply(beta_I, 1, mean)
mean_of_beta_II <- apply(beta_II, 1, mean)
d_mean_of_beta_I <- density(mean_of_beta_I, na.rm = TRUE)
d_mean_of_beta_II <- density(mean_of_beta_II, na.rm = TRUE)

# For Raw data, we already have the matrix of beta values, the d_mean_of_beta_I and d_mean_of_beta_II objects; we need to calculate the densities of the standard deviations, which can be calculated using the function sd():
sd_of_beta_I <- apply(beta_I, 1, sd, na.rm = TRUE)
sd_of_beta_II <- apply(beta_II, 1, sd, na.rm = TRUE)
d_sd_of_beta_I <- density(sd_of_beta_I,)
d_sd_of_beta_II <- density(na.omit(sd_of_beta_II))

preprocessNoob_results <- preprocessNoob(RGset)
beta_preprocessNoob <- getBeta(preprocessNoob_results)
save(beta_preprocessNoob, file = "beta_preprocessNoob.RData")

# Divide normalized beta matrix by probe type and calculate statistics
beta_preprocessNoob_I <- beta_preprocessNoob[rownames(beta_preprocessNoob) %in% dfI$IlmnID,]
beta_preprocessNoob_II <- beta_preprocessNoob[rownames(beta_preprocessNoob) %in% dfII$IlmnID,]
mean_of_beta_preprocessNoob_I <- apply(beta_preprocessNoob_I, 1, mean)
mean_of_beta_preprocessNoob_II <- apply(beta_preprocessNoob_II, 1, mean)
d_mean_of_beta_preprocessNoob_I <- density(mean_of_beta_preprocessNoob_I, na.rm = TRUE)
d_mean_of_beta_preprocessNoob_II <- density(mean_of_beta_preprocessNoob_II, na.rm = TRUE)
sd_of_beta_preprocessNoob_I <- apply(beta_preprocessNoob_I, 1, sd)
sd_of_beta_preprocessNoob_II <- apply(beta_preprocessNoob_II, 1, sd)
d_sd_of_beta_preprocessNoob_I <- density(sd_of_beta_preprocessNoob_I, na.rm = TRUE)
d_sd_of_beta_preprocessNoob_II <- density(sd_of_beta_preprocessNoob_II, na.rm = TRUE)

# Boxplot visualization by group
pdf("Plot_comparison_raw_preprocessNoob control vs disease.pdf", height = 7, width = 15)
par(mfrow = c(2, 3))
plot(d_mean_of_beta_I, col = "blue", main = "Raw beta")
lines(d_mean_of_beta_II, col = "red")
legend("topright",
       legend = c("Probe Type I", "Probe Type II"),
       col = c("blue", "red"),
       lwd = 2,
       cex = 0.9,
       inset = c(0.01, 0.01))
plot(d_sd_of_beta_I, col = "blue", main = "Raw sd")
lines(d_sd_of_beta_II, col = "red")
legend("topright",
       legend = c("Probe Type I", "Probe Type II"),
       col = c("blue", "red"),
       lwd = 2,
       cex = 0.9,
       inset = c(0.01, 0.01))
boxplot(beta,
        las=1,                      
        col=c("green", "magenta")[pheno$Group], 
        ylab='Beta values',
        xlab='Samples',
        main='Boxplot of Beta Values')
legend("bottom", 
       legend = c("Control", "Disease"),
       fill = c("green", "magenta"),
       horiz = TRUE,
       bty = 'n',
       inset = -0.21,       
       xpd = TRUE,          
       cex = 0.9)

plot(d_mean_of_beta_preprocessNoob_I, col = "blue", main = "preprocessNoob beta")
lines(d_mean_of_beta_preprocessNoob_II, col = "red")
legend("topright",
       legend = c("Probe Type I", "Probe Type II"),
       col = c("blue", "red"),
       lwd = 2,
       cex = 0.9,
       inset = c(0.01, 0.01))
plot(d_sd_of_beta_preprocessNoob_I, col = "blue", main = "preprocessNoob sd")
lines(d_sd_of_beta_preprocessNoob_II, col = "red")
legend("topright",
       legend = c("Probe Type I", "Probe Type II"),
       col = c("blue", "red"),
       lwd = 2,
       cex = 0.9,
       inset = c(0.01, 0.01))
boxplot(beta_preprocessNoob,
        las=1,                      
        col=c("green", "magenta")[pheno$Group], 
        ylab='Beta values',
        xlab='Samples',
        main='Boxplot of Normalised Beta Values')
legend("bottom", 
       legend = c("Control", "Disease"),
       fill = c("green", "magenta"),
       horiz = TRUE,
       bty = 'n',
       inset = -0.21,       
       xpd = TRUE,          
       cex = 0.9)
dev.off()

###### STEP 8: PCA ANALYSIS ######
# Perform PCA on normalized beta values
pca_results <- prcomp(t(beta_preprocessNoob), scale = TRUE)
pca_var <- pca_results$sdev^2
pca_var_explained <- pca_var / sum(pca_var)
barplot(pca_var_explained[1:10], 
        names.arg = paste0("PC", 1:10),
        col = "blue",
        main = "Scree Plot: Variance Explained by PCs",
        xlab = "Principal Components",
        ylab = "Proportion of Variance Explained",
        las = 2)

# Plot PCA by Group
par(mfrow = c(1, 1)) 
pdf("PCA_by_group.pdf")
palette(c('orange', 'green'))
plot(pca_results$x[, 1], pca_results$x[, 2], cex = 1, pch = 17, col = pheno$Group, main = 'PCA by group', xlab = "PC1", ylab = "PC2", xlim = c(-1000, 1000), ylim = c(-1000, 1000))
text(pca_results$x[, 1], pca_results$x[, 2], labels = (pheno$SampleID), cex = 0.7, pos = 1)
legend('bottomright', legend = levels(pheno$Group), col = c(1:nlevels(pheno$Group)), pch = 17)
dev.off()

# Plot PCA by Sex
pheno$Sex
pheno$Sex <- as.factor(pheno$Sex)
pdf('PCA_by_sex.pdf')
palette(c('pink', 'aquamarine'))
plot(pca_results$x[, 1], pca_results$x[, 2], cex = 2, pch = 17, col = pheno$Sex, main = 'PCA by sex', xlab = "PC1", ylab = "PC2", xlim = c(-1000, 1000), ylim = c(-1000, 1000))
text(pca_results$x[, 1], pca_results$x[, 2], labels = (pheno$SampleID), cex = 0.7, pos = 1)
legend('bottomright', legend = levels(pheno$Sex), col = c(1:nlevels(pheno$Sex)), pch = 17)
dev.off()

# Plot PCA by Sentrix_ID
pheno$Sentrix_ID
pheno$Sentrix_ID <- as.factor(pheno$Sentrix_ID)
pdf('PCA_by_Sentrix_ID.pdf')
palette(c('red', 'blue', 'green', 'yellow', 'purple'))
plot(pca_results$x[, 1], pca_results$x[, 2], cex = 2, pch = 17, col = pheno$Sentrix_ID, main = 'PCA by Sentrix_ID', xlab = "PC1", ylab = "PC2", xlim = c(-1000, 1000), ylim = c(-1000, 1000))
text(pca_results$x[, 1], pca_results$x[, 2], labels = (pheno$SampleID), cex = 0.7, pos = 1)
legend('bottomright', legend = levels(pheno$Sentrix_ID), col = c(1:nlevels(pheno$Sentrix_ID)), pch = 17)
dev.off()

###### STEP 8: PCA ANALYSIS (continued) ######
# Additional PCA visualizations
pca <- prcomp(beta_preprocessNoob, scale = TRUE)
par(mfrow = c(1, 1))
plot(pca$x[, 1], pca$x[, 2], cex = 2, pch = 18, xlab = "PC1", ylab = "PC2", xlim = c(-5, 7), main = "PCA")
text(pca$x[, 1], pca$x[, 2], labels = rownames(pca$x), pos = 4, cex = 0.5)

# PCA with ggplot2
pcr <- data.frame(pca$rotation[, 1:3], Group = targets$Group)
ggplot(pcr, aes(PC1, PC2, color = Group)) + geom_point(size = 4)
pcr <- data.frame(pca$rotation[, 1:3], Group = targets$Sex)
ggplot(pcr, aes(PC1, PC2, color = Group)) + geom_point(size = 4)

###### STEP 9: DIFFERENTIAL METHYLATION ANALYSIS ######
# Identify differentially methylated probes between CTRL and DIS groups using t-test and Wilcoxon test
# Example for a single probe:
t_test <- t.test(beta_preprocessNoob[1, ] ~ pheno$Group)
t_test
t_test$p.value

# Define functions to apply tests to all probes
My_ttest_function <- function(x) {
  t_test <- t.test(x ~ pheno$Group)
  return(t_test$p.value)
}

# Apply the functions to all probes
pValues_ttest <- apply(beta_preprocessNoob, 1, My_ttest_function)

# Create data.frames with beta values and p-values
final_ttest <- data.frame(beta_preprocessNoob, pValues_ttest)

# Order probes by p-value
final_ttest <- final_ttest[order(final_ttest$pValues_ttest), ]
head(final_ttest)

# Filter for significant probes (p < 0.05)
significant_ttest <- final_ttest[final_ttest$pValues_ttest <= 0.05, ]

###### STEP 10: MULTIPLE TEST CORRECTION ######
# Apply multiple test correction and set a significance threshold of 0.05
corrected_pValues_BH <- p.adjust(final_ttest$pValues_ttest, "BH")
corrected_pValues_Bonf <- p.adjust(final_ttest$pValues_ttest, "bonferroni")
final_ttest_corrected <- data.frame(final_ttest, corrected_pValues_BH, corrected_pValues_Bonf)
head(final_ttest_corrected)

# Visualize distributions of p-values and corrected p-values
colnames(final_ttest_corrected)
box_colors <- c("skyblue", "lightgreen", "salmon")
boxplot(final_ttest_corrected[, 9:11],
        main = "Distribution of Raw and Corrected p-values",
        ylab = "p-values",
        xlab = "Correction Method",
        col = box_colors,
        names = c("Raw", "BH", "Bonferroni"))
legend("bottomright",
       legend = c("Raw", "BH", "Bonferroni"),
       fill = box_colors,
       cex = 0.8)

# How many probes remain after the multiple test correction
dim(final_ttest_corrected[final_ttest_corrected$pValues_ttest <= 0.05, ])
dim(final_ttest_corrected[final_ttest_corrected$corrected_pValues_BH <= 0.05, ])
dim(final_ttest_corrected[final_ttest_corrected$corrected_pValues_Bonf <= 0.05, ])

###### STEP 11: VOLCANO AND MANHATTAN PLOTS ######
# Volcano plot
# Preparation of data
beta_first <- final_ttest_corrected[, 1:8]
beta_first_ctrl <- beta_first[, pheno$Group == "CTRL"]
beta_first_dis <- beta_first[, pheno$Group == "DIS"]
mean_beta_first_ctrl <- apply(beta_first_ctrl, 1, mean)
mean_beta_first_dis <- apply(beta_first_dis, 1, mean)
delta_first <- mean_beta_first_dis - mean_beta_first_ctrl
toVolcPlot <- data.frame(
  delta_beta = delta_first,
  negLog10p = -log10(final_ttest_corrected$pValues_ttest)
)

# Volcano plot with highlight of hypermethylated (Δβ > 0.1 & p < 0.05) and hypomethylated (Δβ < -0.1 & p < 0.05)
pdf("Volcano_plot.pdf")
plot(toVolcPlot$delta_beta, toVolcPlot$negLog10p,
     pch = 16, cex = 0.5,
     xlab = "Delta Beta (DIS - CTRL)",
     ylab = "-log10(p-value)",
     main = "Volcano Plot of Methylation Differences")
hyper <- which(toVolcPlot$delta_beta > 0.1 & final_ttest_corrected$pValues_ttest < 0.05)
points(toVolcPlot$delta_beta[hyper],
       toVolcPlot$negLog10p[hyper],
       col = "red", pch = 16, cex = 0.6)
hypo <- which(toVolcPlot$delta_beta < -0.1 & final_ttest_corrected$pValues_ttest < 0.05)
points(toVolcPlot$delta_beta[hypo],
       toVolcPlot$negLog10p[hypo],
       col = "blue", pch = 16, cex = 0.6)

abline(h = -log10(0.05), col = "darkgray", lty = 2)  # p = 0.05
abline(v = c(-0.1, 0.1), col = "darkgray", lty = 2)  # Δβ = ±0.1
legend("topright", legend = c("Hypermethylated", "Hypomethylated"),
       col = c("red", "blue"), pch = 16, bty = "n")
dev.off()

# Manhattan plot
# (During the lesson, the object final_ttest_corrected was prepared on a subset of probes)
# Annotation of the dataframe with genome annotation for each CpG probe
load('~/Downloads/L2/Illumina450Manifest_clean.RData')
final_ttest_corrected <- data.frame(rownames(final_ttest_corrected), final_ttest_corrected)
colnames(final_ttest_corrected)[1] <- "IlmnID"
final_ttest_corrected_annotated <- merge(final_ttest_corrected, Illumina450Manifest_clean, by = "IlmnID")
input_Manhattan <- final_ttest_corrected_annotated[colnames(final_ttest_corrected_annotated) %in% c("IlmnID", "CHR", "MAPINFO", "pValues_ttest")]
order_chr <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")
input_Manhattan$CHR <- factor(input_Manhattan$CHR, levels = order_chr)
input_Manhattan$CHR <- as.numeric(input_Manhattan$CHR)
manhattan(input_Manhattan, snp = "IlmnID", chr = "CHR", bp = "MAPINFO", p = "pValues_ttest")
pdf("Manhattan_plot.pdf")
manhattan(input_Manhattan, snp = "IlmnID", chr = "CHR", bp = "MAPINFO", p = "pValues_ttest", annotatePval = 0.00001, col = rainbow(24))
dev.off()

###### STEP 12: HEATMAP OF TOP 100 PROBES ######
# Define colorbar from phenotype group and custom color palette
colorbar <- pheno$Group
colorbar <- as.character(factor(colorbar, labels = c("royalblue", "orange")))
input_heatmap = as.matrix(final_ttest[1:100, 1:8])
col2 = colorRampPalette(c("green", "black", "red"))(100)

# Complete linkage 
pdf("heatmap_complete_linkage.pdf")
heatmap.2(input_heatmap,
          col = col2,
          Rowv = TRUE,
          Colv = TRUE,
          dendrogram = "both",
          key = TRUE,
          ColSideColors = colorbar,
          density.info = "none",
          trace = "none",
          scale = "none",
          symm = FALSE,
          main = "Complete Linkage")
legend("topright", legend = c("CTRL", "DIS"), fill = c("royalblue", "orange"), border = NA, bty = "n", cex = 0.8)
dev.off()

# Single linkage
pdf("heatmap_single_linkage.pdf")
heatmap.2(input_heatmap,
          col = col2,
          Rowv = TRUE,
          Colv = TRUE,
          hclustfun = function(x) hclust(x, method = "single"),
          dendrogram = "both",
          key = TRUE,
          ColSideColors = colorbar,
          density.info = "none",
          trace = "none",
          scale = "none",
          symm = FALSE,
          main = "Single Linkage")
legend("topright", legend = c("CTRL", "DIS"), fill = c("royalblue", "orange"), border = NA, bty = "n", cex = 0.8)
dev.off()

# Average linkage
pdf("heatmap_average_linkage.pdf")
heatmap.2(input_heatmap,
          col = col2,
          Rowv = TRUE,
          Colv = TRUE,
          hclustfun = function(x) hclust(x, method = "average"),
          dendrogram = "both",
          key = TRUE,
          ColSideColors = colorbar,
          density.info = "none",
          trace = "none",
          scale = "none",
          symm = FALSE,
          main = "Average Linkage")
legend("topright", legend = c("CTRL", "DIS"), fill = c("royalblue", "orange"), border = NA, bty = "n", cex = 0.8)
dev.off()