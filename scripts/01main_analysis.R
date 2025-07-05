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
install.packages("lib/SummarizedExperiment_1.38.0.tar.gz", repos = NULL)

# STEP 1: DATA IMPORT AND OBJECT CREATION
# List files in the raw data folder
list.files("./data/raw/")

# Load the sample sheet (must have 'array' and 'slide' columns)
sample_sheet <- read.csv("./data/raw/SampleSheet_Report_II.csv", header = TRUE)
sample_sheet

# Load sample sheet using minfi's function
base_dir <- "./data/raw"
targets <- read.metharray.sheet(base_dir)
targets

# Create RGChannelSet object
rgset <- read.metharray.exp(targets = targets)
save(rgset, file = "RGset_Report.RData")
load("RGset_Report.RData")

# Explore rgset object
rgset
str(rgset)

###### STEP 2: EXTRACT RED AND GREEN CHANNELS ######
# Extract fluorescence intensity info for Red and Green channels
red <- data.frame(getRed(rgset))
dim(red)
head(red)
green <- data.frame(getGreen(rgset))
dim(green)
head(green)

###### STEP 3: PROBE INFO AND ADDRESS CHECKS ######
# Explore manifest and probe info
getManifest(rgset)
getManifestInfo(rgset)
getProbeInfo(rgset)
probe_info_i <- data.frame(getProbeInfo(rgset))
dim(probe_info_i)
head(getProbeInfo(rgset, type = "II"))
probe_info_ii <- data.frame(getProbeInfo(rgset, type = "II"))
dim(probe_info_ii)

# Check probe with address 18756452
probe_info_i[probe_info_i$AddressA == "18756452", ]
probe_info_i[probe_info_i$AddressB == "18756452", ]
probe_info_ii[probe_info_ii$AddressB == "18756452", ]
probe_info_18756452 <- probe_info_ii[probe_info_ii$AddressA == "18756452", ]

# Extract Red/Green fluorescence for address 18756452
red_fluorescences <- red[rownames(red) == "18756452", ]
green_fluorescences <- green[rownames(green) == "18756452", ]

# Create summary table for probe 18756452
sample_names <- colnames(red_fluorescences)
df_summary_probes <- data.frame(
  sample = sample_names,
  red_fluorescence = as.numeric(red_fluorescences[1, ]),
  green_fluorescence = as.numeric(green_fluorescences[1, ]),
  type = rep("II", length(sample_names)),
  color = rep("Both", length(sample_names))
)
print(df_summary_probes)

###### STEP 4: CREATE mset_raw OBJECT ######
mset_raw <- preprocessRaw(rgset)
mset_raw
# Note that now the number of rows is 485512, exactly the number of probes according to the Manifest!
save(mset_raw, file = "./data/processed/MSet_raw.RData")
meth <- as.matrix(getMeth(mset_raw))
str(meth)
head(meth)
unmeth <- as.matrix(getUnmeth(mset_raw))
str(unmeth)
head(unmeth)

# Check probe cg01523029 in MethylSet
unmeth[rownames(unmeth) == "cg01523029", ]
meth[rownames(meth) == "cg01523029", ]

load("data/raw/Illumina450Manifest.RData")
illumina450_manifest_clean <- Illumina450Manifest[Illumina450Manifest$CHR != "", ]
illumina450_manifest_clean <- droplevels(illumina450_manifest_clean)
probe_info_cg01523029 <- illumina450_manifest_clean[illumina450_manifest_clean$IlmnID == "cg01523029", ]

# Out-of-band signals: compare rgset and methset for cg01523029
red[rownames(red) == "18756452", ]
green[rownames(green) == "18756452", ]
unmeth[rownames(unmeth) == "cg01523029", ]
meth[rownames(meth) == "cg01523029", ]

###### STEP 5: QUALITY CONTROL ######
# QC plot
pdf("QCplot.pdf")
qc <- getQC(mset_raw)
qc
plotQC(qc)
dev.off()

# Check negative controls
getProbeInfo(rgset, type = "Control")
negative_control <- data.frame(getProbeInfo(rgset, type = "Control"))
table(negative_control$Type)
#head(getProbeInfo(rgset, type = "Control")) - #ensuring rgset contains control probes
controlStripPlot(rgset, controls = "NEGATIVE")

# Calculate detection p-values
detection_pvalue <- detectionP(rgset)
save(detection_pvalue, file = "Detection_pvalue.RData")
over_threshold <- detection_pvalue > 0.05
head(over_threshold)
table(over_threshold)
sum(over_threshold)
failed_probe_message <- sprintf(
  "There are %s failed probes with a detection p-value higher than 0.05",
  sum(over_threshold)
)
over_threshold_per_sample <- colSums(over_threshold)
over_threshold_per_sample
summary(over_threshold)
failed <- data.frame(
  sample = names(over_threshold_per_sample),
  n_failed_positions = over_threshold_per_sample
)
print(failed)
print(failed_probe_message)

###### STEP 6: CALCULATE BETA AND M VALUES, PLOT DENSITIES ######
beta <- getBeta(mset_raw)
m <- getM(mset_raw)
beta_df <- data.frame(beta)
m_df <- data.frame(m)

pheno <- read.csv("data/raw/SampleSheet_Report_II.csv", header = TRUE, stringsAsFactors = TRUE)
pheno$sample_id <- paste(pheno$Sentrix_ID, pheno$Sentrix_Position, sep = "_")

wt_samples <- pheno$sample_id[pheno$Group == "CTRL"]
mut_samples <- pheno$sample_id[pheno$Group == "DIS"]

beta_ctrl <- beta_df[pheno$Group == "CTRL", ]
beta_dis <- beta_df[pheno$Group == "DIS", ]
m_ctrl <- m_df[pheno$Group == "CTRL", ]
m_dis <- m_df[pheno$Group == "DIS", ]

mean_of_beta_ctrl <- apply(beta_ctrl, 1, mean, na.rm = TRUE)
mean_of_beta_dis <- apply(beta_dis, 1, mean, na.rm = TRUE)
mean_of_m_ctrl <- apply(m_ctrl, 1, mean, na.rm = TRUE)
mean_of_m_dis <- apply(m_dis, 1, mean, na.rm = TRUE)

# Plot density of mean methylation values with function and legends
plot_density <- function() {
  par(mfrow = c(1, 2), mar = c(5, 4, 4, 2))
  plot(density(mean_of_beta_ctrl, na.rm = TRUE), main = "Density of Beta Values", col = "black", lwd = 2)
  lines(density(mean_of_beta_dis, na.rm = TRUE), col = "orange", lwd = 2)
  legend("top", legend = c("CTL", "DIS"), col = c("black", "orange"), lwd = 2, cex = 0.9, bty = "n", inset = c(-0.01, -0.01))
  plot(density(mean_of_m_ctrl, na.rm = TRUE), main = "Density of M Values", col = "black", lwd = 2)
  lines(density(mean_of_m_dis, na.rm = TRUE), col = "orange", lwd = 2)
  legend("topright", legend = c("CTL", "DIS"), col = c("black", "orange"), lwd = 2, cex = 0.9, bty = "n", inset = c(-0.01, -0.01))
}
plot_density()

###### STEP 7: NORMALIZATION WITH preprocessNoob ######
# Subset beta values by probe type
df_i <- illumina450_manifest_clean[illumina450_manifest_clean$Infinium_Design_Type == "I", ]
df_i <- droplevels(df_i)
df_ii <- illumina450_manifest_clean[illumina450_manifest_clean$Infinium_Design_Type == "II", ]
df_ii <- droplevels(df_ii)

beta_i <- beta[rownames(beta) %in% df_i$IlmnID, ]
beta_ii <- beta[rownames(beta) %in% df_ii$IlmnID, ]

mean_of_beta_i <- apply(beta_i, 1, mean)
mean_of_beta_ii <- apply(beta_ii, 1, mean)
d_mean_of_beta_i <- density(mean_of_beta_i, na.rm = TRUE)
d_mean_of_beta_ii <- density(mean_of_beta_ii, na.rm = TRUE)

sd_of_beta_i <- apply(beta_i, 1, sd, na.rm = TRUE)
sd_of_beta_ii <- apply(beta_ii, 1, sd, na.rm = TRUE)
d_sd_of_beta_i <- density(sd_of_beta_i)
d_sd_of_beta_ii <- density(na.omit(sd_of_beta_ii))

preprocess_noob_results <- preprocessNoob(rgset)
beta_preprocess_noob <- getBeta(preprocess_noob_results)
save(beta_preprocess_noob, file = "beta_preprocessNoob.RData")

# Divide normalized beta matrix by probe type and calculate stats
beta_preprocess_noob_i <- beta_preprocess_noob[rownames(beta_preprocess_noob) %in% df_i$IlmnID, ]
beta_preprocess_noob_ii <- beta_preprocess_noob[rownames(beta_preprocess_noob) %in% df_ii$IlmnID, ]
mean_of_beta_preprocess_noob_i <- apply(beta_preprocess_noob_i, 1, mean)
mean_of_beta_preprocess_noob_ii <- apply(beta_preprocess_noob_ii, 1, mean)
d_mean_of_beta_preprocess_noob_i <- density(mean_of_beta_preprocess_noob_i, na.rm = TRUE)
d_mean_of_beta_preprocess_noob_ii <- density(mean_of_beta_preprocess_noob_ii, na.rm = TRUE)
sd_of_beta_preprocess_noob_i <- apply(beta_preprocess_noob_i, 1, sd)
sd_of_beta_preprocess_noob_ii <- apply(beta_preprocess_noob_ii, 1, sd)
d_sd_of_beta_preprocess_noob_i <- density(sd_of_beta_preprocess_noob_i, na.rm = TRUE)
d_sd_of_beta_preprocess_noob_ii <- density(sd_of_beta_preprocess_noob_ii, na.rm = TRUE)

# now plotting
pdf("Plot_comparison_raw_preprocessNoob.pdf", height = 7, width = 15)
par(mfrow = c(2, 3))
plot(d_mean_of_beta_i, col = "blue", main = "raw beta")
lines(d_mean_of_beta_ii, col = "red")
plot(d_sd_of_beta_i, col = "blue", main = "raw sd")
lines(d_sd_of_beta_ii, col = "red")
boxplot(beta)

plot(d_mean_of_beta_preprocess_noob_i, col = "blue", main = "preprocessNoob beta")
lines(d_mean_of_beta_preprocess_noob_ii, col = "red")
plot(d_sd_of_beta_preprocess_noob_i, col = "blue", main = "preprocessNoob sd")
lines(d_sd_of_beta_preprocess_noob_ii, col = "red")
boxplot(beta_preprocess_noob)
dev.off()

# Improved boxplot visualization by group
pdf("Plot_comparison_raw_preprocessNoob control vs disease.pdf", height = 7, width = 15)
par(mfrow = c(2, 3))
plot(d_mean_of_beta_i, col = "blue", main = "raw beta")
lines(d_mean_of_beta_ii, col = "red")
legend("topright", legend = c("Probe Type I", "Probe Type II"), col = c("blue", "red"), lwd = 2, cex = 0.9, inset = c(0.01, 0.01))
plot(d_sd_of_beta_i, col = "blue", main = "raw sd")
lines(d_sd_of_beta_ii, col = "red")
legend("topright", legend = c("Probe Type I", "Probe Type II"), col = c("blue", "red"), lwd = 2, cex = 0.9, inset = c(0.01, 0.01))
boxplot(beta,
        las = 2,
        col = c("green", "magenta")[pheno$Group],
        ylab = "Beta values",
        xlab = "Samples",
        main = "Boxplot of Beta Values")
legend("bottom", legend = c("Control", "Disease"), fill = c("green", "magenta"), horiz = TRUE, bty = "n", inset = -0.25, xpd = TRUE, cex = 0.9)
plot(d_mean_of_beta_preprocess_noob_i, col = "blue", main = "preprocessNoob beta")
lines(d_mean_of_beta_preprocess_noob_ii, col = "red")
legend("topright", legend = c("Probe Type I", "Probe Type II"), col = c("blue", "red"), lwd = 2, cex = 0.9, inset = c(0.01, 0.01))
plot(d_sd_of_beta_preprocess_noob_i, col = "blue", main = "preprocessNoob sd")
lines(d_sd_of_beta_preprocess_noob_ii, col = "red")
legend("topright", legend = c("Probe Type I", "Probe Type II"), col = c("blue", "red"), lwd = 2, cex = 0.9, inset = c(0.01, 0.01))
boxplot(beta_preprocess_noob,
        las = 2,
        col = c("green", "magenta")[pheno$Group],
        ylab = "Beta values",
        xlab = "Samples",
        main = "Boxplot of Normalised Beta Values")
legend("bottom", legend = c("Control", "Disease"), fill = c("green", "magenta"), horiz = TRUE, bty = "n", inset = -0.25, xpd = TRUE, cex = 0.9)
dev.off()

###### STEP 8: PCA ANALYSIS ######
# Perform PCA on normalized beta values
pca_results <- prcomp(t(beta_preprocess_noob), scale = TRUE)
par(mfrow = c(1, 1))
plot(pca_results, main = "PCA results", ylim = c(0, 150000))

# Plot PCA by group, sex, and Sentrix_ID
palette(c("orange", "green"))
plot(pca_results$x[, 1], pca_results$x[, 2], cex = 1, pch = 17, col = pheno$Group, main = "PCA by group", xlab = "PC1", ylab = "PC2", xlim = c(-1000, 1000), ylim = c(-1000, 1000))
text(pca_results$x[, 1], pca_results$x[, 2], labels = (pheno$sample_id), cex = 0.7, pos = 1)
legend("bottomright", legend = levels(pheno$Group), col = c(1:nlevels(pheno$Group)), pch = 17)

pheno$Sex <- as.factor(pheno$Sex)
palette(c("pink", "aquamarine"))
plot(pca_results$x[, 1], pca_results$x[, 2], cex = 2, pch = 17, col = pheno$Sex, main = "PCA by sex", xlab = "PC1", ylab = "PC2", xlim = c(-1000, 1000), ylim = c(-1000, 1000))
text(pca_results$x[, 1], pca_results$x[, 2], labels = (pheno$sample_id), cex = 0.7, pos = 1)
legend("bottomright", legend = levels(pheno$Sex), col = c(1:nlevels(pheno$Sex)), pch = 17)

pheno$Sentrix_ID <- as.factor(pheno$Sentrix_ID)
palette(c("red", "blue", "green", "yellow", "purple"))
plot(pca_results$x[, 1], pca_results$x[, 2], cex = 2, pch = 17, col = pheno$Sentrix_ID, main = "PCA by Sentrix_ID", xlab = "PC1", ylab = "PC2", xlim = c(-1000, 1000), ylim = c(-1000, 1000))
text(pca_results$x[, 1], pca_results$x[, 2], labels = (pheno$sample_id), cex = 0.7, pos = 1)
legend("bottomright", legend = levels(pheno$Sentrix_ID), col = c(1:nlevels(pheno$Sentrix_ID)), pch = 17)

###### STEP 8: PCA ANALYSIS (continued) ######
# Additional PCA visualizations
pca <- prcomp(beta_preprocess_noob, scale = TRUE)
par(mfrow = c(1, 1))
plot(pca$x[, 1], pca$x[, 2], cex = 2, pch = 18, xlab = "PC1", ylab = "PC2", xlim = c(-500, 800), main = "PCA")
text(pca$x[, 1], pca$x[, 2], labels = rownames(pca$x), pos = 4, cex = 0.5)

# PCA with ggplot2
install.packages("ggplot2")
library(ggplot2)
pcr <- data.frame(pca$rotation[, 1:3], group = targets$Group)
ggplot(pcr, aes(PC1, PC2, color = group)) + geom_point(size = 4)
pcr <- data.frame(pca$rotation[, 1:3], group = targets$Sex)
ggplot(pcr, aes(PC1, PC2, color = group)) + geom_point(size = 4)

###### STEP 9: DIFFERENTIAL METHYLATION ANALYSIS ######
# Identify differentially methylated probes between CTRL and DIS groups using t-test and Wilcoxon test
# Example for a single probe:
t_test <- t.test(beta_preprocess_noob[1, ] ~ pheno$Group)
t_test
t_test$p.value
wilcox <- wilcox.test(beta_preprocess_noob[1, ] ~ pheno$Group)
wilcox
wilcox$p.value

# Define functions to apply tests to all probes
my_ttest_function <- function(x) {
  t_test <- t.test(x ~ pheno$Group)
  return(t_test$p.value)
}
my_mannwhitney_function <- function(x) {
  wilcox <- wilcox.test(x ~ pheno$Group)
  return(wilcox$p.value)
}

# Apply the functions to all probes
p_values_ttest <- apply(beta_preprocess_noob, 1, my_ttest_function)
p_values_wilcox <- apply(beta_preprocess_noob, 1, my_mannwhitney_function)

# Create data.frames with beta values and p-values
final_ttest <- data.frame(beta_preprocess_noob, p_values_ttest)
final_wilcox <- data.frame(beta_preprocess_noob, p_values_wilcox)

# Order probes by p-value
final_ttest <- final_ttest[order(final_ttest$p_values_ttest), ]
head(final_ttest)
final_wilcox <- final_wilcox[order(final_wilcox$p_values_wilcox), ]
head(final_wilcox)

# Filter for significant probes (p < 0.05)
significant_ttest <- final_ttest[final_ttest$p_values_ttest <= 0.05, ]
significant_wilcox <- final_wilcox[final_wilcox$p_values_wilcox <= 0.05, ]
intersection <- intersect(rownames(significant_ttest), rownames(significant_wilcox))
length(intersection)
# (Alternatively, dmpFinder can be used)

###### STEP 10: MULTIPLE TEST CORRECTION ######
# Apply multiple test correction and set a significance threshold of 0.05
corrected_p_values_bh <- p.adjust(final_ttest$p_values_ttest, "BH")
corrected_p_values_bonf <- p.adjust(final_ttest$p_values_ttest, "bonferroni")
final_ttest_corrected <- data.frame(final_ttest, corrected_p_values_bh, corrected_p_values_bonf)
head(final_ttest_corrected)

# Visualize distributions of p-values and corrected p-values
colnames(final_ttest_corrected)
boxplot(final_ttest_corrected[, 9:11]) # See BI's notes for the plot

# How many probes survive the multiple test correction?
dim(final_ttest_corrected[final_ttest_corrected$p_values_ttest <= 0.05, ])
dim(final_ttest_corrected[final_ttest_corrected$corrected_p_values_bh <= 0.05, ])
dim(final_ttest_corrected[final_ttest_corrected$corrected_p_values_bonf <= 0.05, ])
# Note: For both tests, there may be none that pass the test.

###### STEP 11: VOLCANO AND MANHATTAN PLOTS ######
# Volcano plot
beta_first <- final_ttest_corrected[, 1:8]
beta_first_ctrl <- beta_first[, pheno$Group == "CTRL"]
mean_beta_first_ctrl <- apply(beta_first_ctrl, 1, mean)
beta_first_dis <- beta_first[, pheno$Group == "DIS"]
mean_beta_first_dis <- apply(beta_first_dis, 1, mean)
delta_first <- mean_beta_first_dis - mean_beta_first_ctrl
to_volc_plot <- data.frame(delta_beta = delta_first, neg_log10p = -log10(final_ttest_corrected$p_values_ttest))
plot(to_volc_plot$delta_beta, to_volc_plot$neg_log10p,
     pch = 16, cex = 0.5,
     xlab = "Delta Beta (DIS - CTRL)",
     ylab = "-log10(p-value)",
     main = "Volcano Plot of Methylation Differences")
hyper <- which(to_volc_plot$delta_beta > 0.1 & final_ttest_corrected$p_values_ttest < 0.05)
points(to_volc_plot$delta_beta[hyper],
       to_volc_plot$neg_log10p[hyper],
       col = "red", pch = 16, cex = 0.6)
hypo <- which(to_volc_plot$delta_beta < -0.1 & final_ttest_corrected$p_values_ttest < 0.05)
points(to_volc_plot$delta_beta[hypo],
       to_volc_plot$neg_log10p[hypo],
       col = "blue", pch = 16, cex = 0.6)
abline(h = -log10(0.05), col = "darkgray", lty = 2)  # p = 0.05
abline(v = c(-0.1, 0.1), col = "darkgray", lty = 2)  # Delta Beta = ±0.1
legend("topright", legend = c("Hypermethylated", "Hypomethylated"),
       col = c("red", "blue"), pch = 16, bty = "n")

# Manhattan plot
install.packages("qqman")
library(qqman)
# (During the lesson, the object final_ttest_corrected was prepared on a subset of probes)
# Annotate dataframe with genome annotation for each CpG probe
final_ttest_corrected <- data.frame(rownames(final_ttest_corrected), final_ttest_corrected)
colnames(final_ttest_corrected)[1] <- "IlmnID"
final_ttest_corrected_clean <- final_ttest_corrected[, !duplicated(colnames(final_ttest_corrected))]
final_ttest_corrected_annotated <- merge(final_ttest_corrected_clean, illumina450_manifest_clean, by = "IlmnID")
input_manhattan <- final_ttest_corrected_annotated[colnames(final_ttest_corrected_annotated) %in% c("IlmnID", "CHR", "MAPINFO", "p_values_ttest")]
order_chr <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")
input_manhattan$CHR <- factor(input_manhattan$CHR, levels = order_chr)
input_manhattan$CHR <- as.numeric(input_manhattan$CHR)
manhattan(input_manhattan, 
          snp = "IlmnID", 
          chr = "CHR", 
          bp = "MAPINFO", 
          p = "p_values_ttest", 
          annotatePval = 0.00001, 
          col = rainbow(24))
pdf("Manhattan_plot.pdf")
manhattan(input_manhattan, snp = "IlmnID", chr = "CHR", bp = "MAPINFO", p = "p_values_ttest", annotatePval = 0.00001, col = rainbow(24))
dev.off()

###### STEP 12: HEATMAP OF TOP 100 PROBES ######
library(gplots)
colorbar <- as.character(factor(pheno$Group, labels = c("royalblue", "orange")))
input_heatmap <- as.matrix(final_ttest[1:100, 1:8])
col2 <- colorRampPalette(c("green", "black", "red"))(100)
# Complete linkage
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
legend("topright", legend = c("CTRL", "DIS"), fill = c("royalblue", "orange"), border = NA, bty = "n", cex = 0.8, xpd = TRUE, inset = c(0, -0.10))
# Single linkage
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
legend("topright", legend = c("CTRL", "DIS"), fill = c("royalblue", "orange"), border = NA, bty = "n", cex = 0.8, xpd = TRUE, inset = c(0, -0.10))
# Average linkage
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
legend("topright", legend = c("CTRL", "DIS"), fill = c("royalblue", "orange"), border = NA, bty = "n", cex = 0.8, xpd = TRUE, inset = c(0, -0.10))
# Custom color palette
col2 <- colorRampPalette(c("green", "black", "red"))(100)
heatmap.2(input_heatmap, col = col2, Rowv = TRUE, Colv = TRUE, dendrogram = "both", key = TRUE, ColSideColors = colorbar, density.info = "none", trace = "none", scale = "none", symm = FALSE)


