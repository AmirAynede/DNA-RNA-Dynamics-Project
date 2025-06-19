######               STEP 0              ###### PREPERATION(ENVIRONMENT SETTING)
'''
# As always, first set the working directory: 
# for the project you need the one in which you have the "InputData" folder
getwd()
setwd("~/Dropbox/DRD_2025/Final_Report-20250521/")
'''
# Install necessary packages:

library(minfi)
library(minfiData)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(shinyMethyl)
library(AnnotationDbi)
library(sva)
install.packages("./library/SummarizedExperiment_1.38.0.tar.gz", repos = NULL)

load('./data/raw/Illumina450Manifest_clean.RData')

######     STEP ONE    ######
###### Import raw data ######
###### List the files included in the folder ######

list.files("./data/raw/")

#SAMPLE SHEET FILE NEEDS TO BE DESIGNED CORRECTLY IN ORDER TO BE ABLE TO LOAD IT: have a look at the csv
#CORRECT SAMPLE SHEET FORMAT ==> REQUIRED COLUMNS: array and slide columns
SampleSheet <- read.csv("./data/raw/SampleSheet_Report_II.csv",header=T)

# sample of type of sample sheet that can actually be read by the loading function
SampleSheet 

# LOADING FUNCTION
# Set the directory in which the raw data are stored
# load the samplesheet using the function read.metharray.sheet
baseDir <- ("./data/raw")
# basedir arg combines array and sample name (slide) columns in the base name column
targets <- read.metharray.sheet(baseDir) 
# combines idat files that share name: 
# dataframe generated from read.metharray.sheet(baseDir) containing info necessary for loading function
targets  

# Create an object of class RGChannelSet using the function read.metharray.exp
RGset <- read.metharray.exp(targets = targets)
save(RGset,file="RGset_Report.RData")

# Let's explore the RGset object: complex R data 
# values are organised in specific slots identified by a specific name @name (@ access operator)
# to access do not use @ but specific function found in the help page ?RGChannelSet
RGset
str(RGset)

######               STEP TWO              ######
###### Create the Red and Green dataframes ######

# We extract the Green and Red Channels using the functions getGreen and getRed
# extract fluorescence intensity info associated to red
Red <- data.frame(getRed(RGset)) 
dim(Red)
head(Red)
# extract fluorescence intensity info associated to green
Green <- data.frame(getGreen(RGset))
dim(Green)
head(Green)

# why 622399 rows instead of 485k as expected? (we should mention this in out manuscript)
# We assign 2 ids for infinium I and 1 id for infinium II (2 + 1)

######               STEP THREE             ######
###### Find the Red and Green fluorescences ######
######  assigned to the address: '18756452' ######

# From the help page, you see that getManifest and getManifestInfo are two accessor functions for this class. Let's check their results:
getManifest(RGset)
getManifestInfo(RGset)
# getManifestInfo does not seem really useful! :-D

# What about the getProbeInfo() function? It returns a data.frame giving the type I, type II or control probes
getProbeInfo(RGset)
# Only the first and the last rows of the object are printed. We will create an object df containing the result of this function:
ProbeInfo_I <- data.frame(getProbeInfo(RGset))
# Strange...this df has XXXXXX rows, not 485512...why?
dim(ProbeInfo_I)
# It only stores Type I by default.
# The getProbeInfo function returns only Type I probes. However, this is not really specified in the help page. 
# Take home message: when you work with a package, try before you trust!
head(getProbeInfo(RGset, type = "II"))
ProbeInfo_II <- data.frame(getProbeInfo(RGset, type = "II"))
dim(ProbeInfo_II)
350036+135476 #485512 (we should mention this in out manuscript)

# CHECKING ADDRESS AND VERYFING IF IT IS A OR B

# From probes to Addresses
head(ProbeInfo_I)
head(ProbeInfo_II)

# I want to check the probe having address 18756452 in the Type I array: none
ProbeInfo_I[ProbeInfo_I$AddressA=="18756452",]
ProbeInfo_I[ProbeInfo_I$AddressB=="18756452",]

# there is not an AddressB_ID with code 18756452;
ProbeInfo_II[ProbeInfo_II$AddressB=="18756452",] 
# there is an AddressA with code 18756452 and it is associated to the probe cg25192902, type II
ProbeInfo_II[ProbeInfo_II$AddressA=="18756452",]

# Using this commands you can fill the table on the word document
Red[rownames(Red)=="18756452",]
Green[rownames(Green)=="18756452",]

# Are there out of band signals??? (optional:check it later)


#ask FO to explain

#green_probes <-Green['18756452',] # =Red[rownames...]
#green_probes
#red_probes <-Red['18756452',]
#red_probes
#the type 
#type <- Illumina450Manifest_clean[Illumina450Manifest_clean$Adress=="18756452",'Infinium_Design_Type']
#type <- droplevels(type) #to remove the levels


######          STEP FOUR         ######
###### Create the object MSet.raw ######
######                            ######

## CHUNK 4: Extract methylated and unmethylated signals using the function MSet.raw
MSet.raw <- preprocessRaw(RGset)
MSet.raw
# Note that now the number of rows is 485512, exactly the number of probes according to the Manifest!
save(MSet.raw,file="./data/processed/MSet_raw.RData")
Meth <- as.matrix(getMeth(MSet.raw))
str(Meth)
head(Meth)
Unmeth <- as.matrix(getUnmeth(MSet.raw))
str(Unmeth)
head(Unmeth)


# FO changed the probe (make sure if it is not FU)

# Let's check what happens to the probes that we considered before when we move from RGset to MethylSet
# cg25192902, type II (the first probe in the Illumina450Manifest_clean object)
Red[rownames(Red)=="18756452",]
Green[rownames(Green)=="18756452",]
Unmeth[rownames(Unmeth)=="cg01523029",]
Meth[rownames(Meth)=="cg01523029",] # changed the probe ID-FO

#### OUT OF BAND SIGNALS: differences bw RGset and Methset
# cg25192902, type II 

#(make sure if it is not FU)
Illumina450Manifest_clean[Illumina450Manifest_clean$IlmnID=="cg01523029",]
Red[rownames(Red)=="18756452",]
# But the two addresses are present also when I look at the Green object: these are out of band signals!
Green[rownames(Green)=="18756452",]
Unmeth[rownames(Unmeth)=="cg25192902",]
Meth[rownames(Meth)=="cg25192902",]


####### Step 5######
#QCplot
#Ask FO

pdf("QCplot.pdf")
qc <- getQC(MSet.raw)
qc
plotQC(qc)
dev.off()

#	check the intensity of negative controls using minfi
getProbeInfo(RGset, type = "Control")
NegativeControl <- data.frame(getProbeInfo(RGset, type = "Control"))
table(NegativeControl$Type)
#head(getProbeInfo(RGset, type = "Control")) - #ensuring RGset contains control probes
controlStripPlot(RGset, controls="NEGATIVE")

#	calculate detection pValues; for each sample, how many probes have a detection p-value higher than the threshold assigned to each group?
Detection_pvalue <- detectionP(RGset)
save(Detection_pvalue,file="Detection_pvalue.RData")

Over_Threshold <- Detection_pvalue>0.05
head(Over_Threshold)
table(Over_Threshold)
sum(Over_Threshold)
paste0(sprintf('There are %s failed probes with a detection p-value higher than 0.05', sum(failed)))
Over_Threshold_Per_Sample <- colSums(Over_Threshold)
Over_Threshold_Per_Sample
summary(Over_Threshold)#what the prof suggested

#####step 6######
#6.	Calculate raw beta and M values and plot the densities of mean methylation values,
beta <- getBeta(MSet.raw)
M <- getM(MSet.raw)
beta_df <- data.frame(beta)
M_df <- data.frame(M)

pheno <-SampleSheet
pheno$Group
beta_ctrl <- beta_df[pheno$Group=="CTRL",]
beta_dis <- beta_df[pheno$Group=="DIS",]
M_ctrl <- M_df[pheno$Group=="CTRL",]
M_dis <- M_df[pheno$Group=="DIS",]


mean_of_beta_ctrl <- apply(beta_ctrl,1,mean,na.rm=T)
mean_of_beta_dis <- apply(beta_dis,1,mean,na.rm=T)
mean_of_M_ctrl <- apply(M_ctrl,1,mean,na.rm=T)
mean_of_M_dis <- apply(M_dis,1,mean,na.rm=T)

#plotting the density of mean methylation values
##########need to fix how to put it all in both beta and m in one plot, also the names of the plots
pdf("Density_Mean_Methylation_Values.pdf")
par(mfrow=c(1,2))
plot(density(mean_of_beta_ctrl,na.rm=T), col="orange", lwd=2, main="Density of Mean Beta Values")
#lines(density(mean_of_M_ctrl, na.rm=T), col="purple",)
plot(density(mean_of_M_ctrl, na.rm=T), col="purple")
plot(density(mean_of_beta_dis,na.rm=T), col="orange")
#lines(density(mean_of_M_dis),col='purple')
plot(density(mean_of_M_dis, na.rm=T), col="purple")
dev.off()

########Step 7##########
#	Normalize the data using PreprocessNoob
#part1 - We want to subset the dataframe of beta values according to Type I and II
dfI <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="I",]
dfI <- droplevels(dfI)
dfII <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="II",]
dfII <- droplevels(dfII)

beta_I <- beta[rownames(beta) %in% dfI$IlmnID,]
beta_II <- beta[rownames(beta) %in% dfII$IlmnID,]

mean_of_beta_I <- apply(beta_I,1,mean)
mean_of_beta_II <- apply(beta_II,1,mean)
d_mean_of_beta_I <- density(mean_of_beta_I,na.rm=T)
d_mean_of_beta_II <- density(mean_of_beta_II,na.rm=T)
  # we could also skip calculation of density and just plot density like in previous step.

# For Raw data, we already have the matrix of beta values, the d_mean_of_beta_I and d_mean_of_beta_II objects; we need to calculate the densities of the standard deviations, which can be calculated using the function sd():
sd_of_beta_I <- apply(beta_I,1,sd,na.rm=T)
sd_of_beta_II <- apply(beta_II,1,sd,na.rm=T)
d_sd_of_beta_I <- density(sd_of_beta_I,)
d_sd_of_beta_II <- density(sd_of_beta_II) #doesn't work. We need to find out why- all na values should have been omited
d_sd_of_beta_II <- density(na.omit(sd_of_beta_II))

?preprocessNoob
# According to the help page, the input can be an RGset  or MehylSet object. Let's load the RGset object:
load("~/Dropbox/DRD_2025/3/RGset.RData")
RGset
?GenomicRatioSet

preprocessNoob_results <- preprocessNoob(RGset)
beta_preprocessNoob <- getBeta(preprocessNoob_results)
save(beta_preprocessNoob, file="beta_preprocessNoob.RData")

#divide beta_preprocessNoob matrix according to type I and type II probes, calculate the mean and the standartd deviation for each probe across the 8 samples and calculate the density distributions
beta_preprocessNoob_I <- beta_preprocessNoob[rownames(beta_preprocessNoob) %in% dfI$IlmnID,]
beta_preprocessNoob_II <- beta_preprocessNoob[rownames(beta_preprocessNoob) %in% dfII$IlmnID,]
mean_of_beta_preprocessNoob_I <- apply(beta_preprocessNoob_I,1,mean)
mean_of_beta_preprocessNoob_II <- apply(beta_preprocessNoob_II,1,mean)
d_mean_of_beta_preprocessNoob_I <- density(mean_of_beta_preprocessNoob_I,na.rm=T)
d_mean_of_beta_preprocessNoob_II <- density(mean_of_beta_preprocessNoob_II,na.rm=T)
sd_of_beta_preprocessNoob_I <- apply(beta_preprocessNoob_I,1,sd)
sd_of_beta_preprocessNoob_II <- apply(beta_preprocessNoob_II,1,sd)
d_sd_of_beta_preprocessNoob_I <- density(sd_of_beta_preprocessNoob_I,na.rm=T)
d_sd_of_beta_preprocessNoob_II <- density(sd_of_beta_preprocessNoob_II,na.rm=T)


#now plotting
pdf("Plot_comparison_raw_preprocessNoob.pdf",height=7,width=15)
par(mfrow=c(2,3))
plot(d_mean_of_beta_I,col="blue",main="raw beta")
lines(d_mean_of_beta_II,col="red")
plot(d_sd_of_beta_I,col="blue",main="raw sd")
lines(d_sd_of_beta_II,col="red")
boxplot(beta)

plot(d_mean_of_beta_preprocessNoob_I,col="blue",main="preprocessNoob beta")
lines(d_mean_of_beta_preprocessNoob_II,col="red")
plot(d_sd_of_beta_preprocessNoob_I,col="blue",main="preprocessNoob sd")
lines(d_sd_of_beta_preprocessNoob_II,col="red")
boxplot(beta_preprocessNoob)
dev.off()

#with better scale, but i'm not sure it helped
par(mfrow=c(2,3))

plot(d_mean_of_beta_I, col="blue", main="Raw Beta", ylim=c(0, 5))
lines(d_mean_of_beta_II, col="red")
plot(d_sd_of_beta_I, col="blue", main="Raw SD", ylim=c(0, 30))
lines(d_sd_of_beta_II, col="red")
boxplot(beta, main="Raw Beta Boxplot", ylim=c(0, 1))

plot(d_mean_of_beta_preprocessNoob_I, col="blue", main="PreprocessNoob Beta", ylim=c(0, 5))
lines(d_mean_of_beta_preprocessNoob_II, col="red")
plot(d_sd_of_beta_preprocessNoob_I, col="blue", main="PreprocessNoob SD", ylim=c(0, 60))
lines(d_sd_of_beta_preprocessNoob_II, col="red")
boxplot(beta_preprocessNoob, main="PreprocessNoob Beta Boxplot", ylim=c(0, 1))

#boxplot of normalised values
#attempt1
boxplot(beta_preprocessNoob, ylab='probe mean', xlab='probes mean', main='Boxplot of normalised beta values',col=pheno$Group,r..)
#attempt ai
pheno$Group <- as.factor(pheno$Group)
# Create boxplot for beta values, grouped by sample
par(mfrow=c(1,1))
boxplot(beta_preprocessNoob,
        las=2,                      # rotate x-axis labels
        col=c("skyblue", "green")[pheno$Group],  # color by group
        ylab='Beta values',
        xlab='Samples',
        main='Boxplot of Normalised Beta Values')
# legend
legend("topright", legend=levels(pheno$Group), fill=c("skyblue", "green"))

#now plotting with boxplot based on sample
pdf("Plot_comparison_raw_preprocessNoob control vs disease.pdf",height=7,width=15)
par(mfrow=c(2,3))
plot(d_mean_of_beta_I,col="blue",main="raw beta")
lines(d_mean_of_beta_II,col="red")
plot(d_sd_of_beta_I,col="blue",main="raw sd")
lines(d_sd_of_beta_II,col="red")
boxplot(beta,
        las=2,                      # rotate x-axis labels
        col=c("skyblue", "green")[pheno$Group],  # color by group
        ylab='Beta values',
        xlab='Samples',
        main='Boxplot of Beta Values')

plot(d_mean_of_beta_preprocessNoob_I,col="blue",main="preprocessNoob beta")
lines(d_mean_of_beta_preprocessNoob_II,col="red")
plot(d_sd_of_beta_preprocessNoob_I,col="blue",main="preprocessNoob sd")
lines(d_sd_of_beta_preprocessNoob_II,col="red")
boxplot(beta_preprocessNoob,
        las=2,                      # rotate x-axis labels
        col=c("skyblue", "green")[pheno$Group],  # color by group
        ylab='Beta values',
        xlab='Samples',
        main='Boxplot of Normalised Beta Values')
dev.off()
# we should fix the sample visualisation

###Step 8######
#8.	Perform a PCA on the matrix of normalized beta values 
pca_results <- prcomp(t(beta_preprocessNoob),scale=T)
par(mfrow=c(1,1))
plot(pca_results,main="PCA results",ylim=c(0,150000))

#the principal components are stored in the x slot of the pca_results object(list)
pca_results$x
#? Do the samples divide according to the group,sex of the samples or batch? 

plot(pca_results$x[,1],pca_results$x[,2],cex=1,pch=2,xlab="PC1",ylab="PC2",main='PCA', xlim=c(-1000,1000),ylim=c(-1000,1000))
text(pca_results$x[,1], pca_results$x[,2], labels=(pheno$SampleID), cex=0.7, pos=1)

pheno$Group
palette(c('orange','green'))
plot(pca_results$x[,1],pca_results$x[,2],cex=1,pch=17,col=pheno$Group,main='PCA by group',xlab="PC1",ylab="PC2",xlim=c(-1000,1000),ylim=c(-1000,1000))
text(pca_results$x[,1], pca_results$x[,2], labels=(pheno$SampleID), cex=0.7, pos=1)
legend('bottomright', legend=levels(pheno$Group), col=c(1:nlevels(pheno$Group)), pch=17)

pheno$Sex
pheno$Sex <- as.factor(pheno$Sex)
palette(c('pink','aquamarine'))
plot(pca_results$x[,1],pca_results$x[,2],cex=2,pch=17,col=pheno$Sex,main='PCA by sex',xlab="PC1",ylab="PC2",xlim=c(-1000,1000),ylim=c(-1000,1000))
text(pca_results$x[,1], pca_results$x[,2], labels=(pheno$SampleID), cex=0.7, pos=1)
legend('bottomright', legend=levels(pheno$Sex), col=c(1:nlevels(pheno$Sex)), pch=17)

pheno$Sentrix_ID
pheno$Sentrix_ID <- as.factor(pheno$Sentrix_ID)
palette(c('red','blue','green','yellow','purple'))
plot(pca_results$x[,1],pca_results$x[,2],cex=2,pch=17,col=pheno$Sentrix_ID,main='PCA by Sentrix_ID',xlab="PC1",ylab="PC2",xlim=c(-1000,1000),ylim=c(-1000,1000))
text(pca_results$x[,1], pca_results$x[,2], labels=(pheno$SampleID), cex=0.7, pos=1)
legend('bottomright', legend=levels(pheno$Sentrix_ID), col=c(1:nlevels(pheno$Sentrix_ID)), pch=17)

###Step 9######
#identify differentially methylated probes between CTRL and DIS groups using the function
#lets use the t-test and the wilcox test to compare the groups
t_test <- t.test(beta_preprocessNoob[1,] ~ pheno$Group)
t_test
wilcox <- wilcox.test(beta_preprocessNoob[1,]~ pheno$Group)
wilcox
# We are interested in the "p.value" element of these lists.
t_test$p.value
wilcox$p.value

# We will use the apply() function, that we have already met in Lesson 4 to calculate the mean of each row of the matrix. However, unlike the mean() function, t.test() and wilcox.test() functions do not return a value, but a list; therefore, we have to create an ad hoc function:
#helpful when looking at huge data saets, as it allows to apply a function to each row of a matrix or dataframe

My_ttest_function <- function(x) {
  t_test <- t.test(x~ pheno$Group)
  return(t_test$p.value)
} 
My_mannwhitney_function <- function(x) {
  wilcox <- wilcox.test(x~ pheno$Group)
  return(wilcox$p.value)
}

# Let's apply the function to few rows of beta_preprocessQuantile
first20k_beta_preprocessNoob <- beta_preprocessNoob[1:20000,]
pValues_ttest_20k <- apply(first20k_beta_preprocessNoob,1, My_ttest_function)
pValues_ttest_20k
pValues_wilcox_20k <- apply(first20k_beta_preprocessNoob, 1, My_mannwhitney_function)
pValues_wilcox_20k

# We can create a data.frame with all the beta values and the pValue column
final_ttest_20k <- data.frame(first20k_beta_preprocessNoob, pValues_ttest_20k)
final_wilcox_20k <- data.frame(first20k_beta_preprocessNoob, pValues_wilcox_20k)

# We can order the probes on the basis of the pValues column (from the smallest to the largest value)
final_ttest_20k <- final_ttest_20k[order(final_ttest_20k$pValues_ttest_20k),]
head(final_ttest_20k)

final_wilcox_20k <- final_wilcox_20k[order(final_wilcox_20k$pValues_wilcox_20k),]
head(final_wilcox_20k)

#filter for significant probes < 0.05
significant_ttest <- final_ttest_20k[final_ttest_20k$pValues_ttest_20k <= 0.05, ]
significant_wilcox <- final_wilcox_20k[final_wilcox_20k$pValues_wilcox_20k <= 0.05, ]

intersection <- intersect(rownames(significant_ttest),rownames(significant_wilcox))
length(intersection)
    ### we could use dmpFinder instead if you like

####Step 10######
#	Apply multiple test correction and set a significant threshold of 0.05.
corrected_pValues_BH <- p.adjust(final_ttest_20k$pValues_ttest_20k,"BH")
corrected_pValues_Bonf <- p.adjust(final_ttest_20k$pValues_ttest_20k,"bonferroni")
final_ttest_20k_corrected <- data.frame(final_ttest_20k, corrected_pValues_BH, corrected_pValues_Bonf)
head(finalttest_20k_corrected)

# We can visualize the distributions of the p-values and of the corrected p-values by boxplots:
colnames(final_ttest_20k_corrected)
boxplot(final_ttest_20k_corrected[,9:11])

# How many probes survive the multiple test correction?
dim(final_ttest_20k_corrected[final_ttest_20k_corrected$pValues_ttest<=0.05,])
dim(final_ttest_20k_corrected[final_ttest_20k_corrected$corrected_pValues_BH<=0.05,])
dim(final_ttest_20k_corrected[final_ttest_20k_corrected$corrected_pValues_Bonf<=0.05,])
  # for both test, I have none that pass the test. Is that a problem?

####Step 11######
#	Produce a volcano plot and a Manhattan plot
beta_first20k <- final_ttest_20k_corrected[,1:8]

beta_first20k_ctrl <- beta_first20k[,pheno$Group=="CTRL"]
mean_beta_first20k_ctrl <- apply(beta_first20k_ctrl,1,mean)
beta_first20k_dis <- beta_first20k[,pheno$Group=="DIS"]
mean_beta_first20k_dis <- apply(beta_first20k_dis,1,mean)

# Now we can calculate the difference between average values:
delta_first20k <- mean_beta_first20k_dis-mean_beta_first20k_ctrl
head(delta_first20k)

# Now we create a dataframe with two columns, one containing the delta values and the other with the -log10 of p-values
toVolcPlot <- data.frame(delta_first20k, -log10(final_ttest_20k_corrected$pValues_ttest))
head(toVolcPlot)
plot(toVolcPlot[,1], toVolcPlot[,2])

# Now I want to highlight the probes (that is, the points), that have a nominal pValue<0.01 and a delta > 0.1) - I set threshold to 0.05-FO
pdf("Volcano_plot.pdf")
plot(toVolcPlot[,1], toVolcPlot[,2],pch=16,cex=0.5,main='volcano plot')
toHighlight <- toVolcPlot[toVolcPlot[,1]>0.1 & toVolcPlot[,2]>(-log10(0.05)),]
head(toHighlight)
points(toHighlight[,1], toHighlight[,2],pch=16,cex=0.7,col="red")
dev.off()

# Now we can create a Manhattan plot
install.packages("qqman")
library(qqman)

# To calculate the Manhattan plot, during the lesson we will use the object final_ttest_reduced, that I prepared by applying the ttest function on a subset of the probes of the dataset. The object is in the Lesson 6 folder:
load("~/Dropbox/DRD_2025/6/final_ttest_reduced.RData") # I used the one the prof created but we should probably ask how it was made tao see if it applies to our case
# First we have to annotate our dataframe, that is add genome annotation information for each cpg probe. We will use the Illumina450Manifest_clean object that we previously created:
load('~/Dropbox/DRD_2025/2/Illumina450Manifest_clean.RData')

# We will use the merge() function to merge the final_ttest_corrected with the Illumina450Manifest_clean object
?merge
# The merge function performs the merging by using a column which is common to two dataframes and which has the same name in the two dataframes
head(Illumina450Manifest_clean)
head(final_ttest_reduced)

# We want to merge on the basis of the CpG probes, but unfortunately in the final_ttest_corrected object the CpG probes are stored in the rownames, not in a column. We can overcome this issue as follows:
final_ttest_reduced <- data.frame(rownames(final_ttest_reduced),final_ttest_reduced)
head(final_ttest_reduced)
colnames(final_ttest_reduced)
colnames(final_ttest_reduced)[1] <- "IlmnID"
colnames(final_ttest_reduced)

final_ttest_reduced_annotated <- merge(final_ttest_reduced, Illumina450Manifest_clean,by="IlmnID")
dim(final_ttest_reduced)
dim(Illumina450Manifest_clean)
dim(final_ttest_reduced_annotated)
str(final_ttest_reduced_annotated)

# Note that the dataframe is automathically reordered on the basis of alphabetical order of the column used for merging
head(final_ttest_reduced_annotated)

# Now we can create the input for the Manhattan plot analysis. The input object should contain 4 info: probe, chromosome, position on the chromosome and p-value. We will select these columns in the final_ttest_corrected_annotated object
input_Manhattan <- final_ttest_reduced_annotated[colnames(final_ttest_reduced_annotated) %in% c("IlmnID","CHR","MAPINFO","pValues_ttest")]
dim(input_Manhattan)
head(input_Manhattan)
str(input_Manhattan$CHR)
levels(input_Manhattan$CHR)
# It is better to reorder the levels of the CHR
order_chr <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
input_Manhattan$CHR <- factor(input_Manhattan$CHR,levels=order_chr )
levels(input_Manhattan$CHR)

#creation of manhattan plot
input_Manhattan$CHR <- as.numeric(input_Manhattan$CHR)
table(input_Manhattan$CHR)

# and finally we can produce our Manhattan plot
manhattan(input_Manhattan, snp="IlmnID",chr="CHR", bp="MAPINFO", p="pValues_ttest" )
# what is the blue line?
-log10(0.00001)
pdf("Manhattan_plot.pdf")
manhattan(input_Manhattan, snp="IlmnID",chr="CHR", bp="MAPINFO", p="pValues_ttest",annotatePval = 0.00001,col=rainbow(24) )
dev.off()

####Step 12######
#	Produce an heatmap of the top 100 significant, differentially methylated probes
library(gplots)
pheno$Group
colorbar <- c("green","orange","orange","green","orange","orange","green","green")

# In the following lines we will compare the results of hierchical clustering using different methods.
input_heatmap=as.matrix(final_ttest_20k[1:100,1:8])
# Complete (default options)
heatmap.2(input_heatmap,col=terrain.colors(100),Rowv=T,Colv=T,dendrogram="both",key=T,ColSideColors=colorbar,density.info="none",trace="none",scale="none",symm=F,main="Complete linkage")

# You can see that hierarchical clustering divides well control group and disease group; 

# Single
heatmap.2(input_heatmap,col=terrain.colors(100),Rowv=T,Colv=T,hclustfun = function(x) hclust(x,method = 'single'),dendrogram="both",key=T,ColSideColors=colorbar,density.info="none",trace="none",scale="none",symm=F,main="Single linkage")

# Average
heatmap.2(input_heatmap,col=terrain.colors(100),Rowv=T,Colv=T,hclustfun = function(x) hclust(x,method = 'average'),dendrogram="both",key=T,ColSideColors=colorbar,density.info="none",trace="none",scale="none",symm=F,main="Average linkage")

# You can set your palette of colors
col2=colorRampPalette(c("green","black","red"))(100)
heatmap.2(input_heatmap,col=col2,Rowv=T,Colv=T,dendrogram="both",key=T,ColSideColors=colorbar,density.info="none",trace="none",scale="none",symm=F)

#####Sorry guys. This is super sloppy