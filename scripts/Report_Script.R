# As always, first set the working directory: 
# for the project you need the one in which you have the "InputData" folder
getwd()
setwd('Final_Report-20250606')

# Install necessary packages:

library(minfi)
library(minfiData)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(shinyMethyl)
library(AnnotationDbi)
install.packages("../../../Practical/L3/SummarizedExperiment_1.38.0.tar.gz", repos = NULL)

######     STEP ONE    ######
###### Import raw data ######
######                 ######

###### List the files included in the folder ######

list.files("Input_Data/")

# SAMPLE SHEET FILE NEEDS TO BE DESIGNED CORRECTLY IN ORDER TO BE ABLE TO LOAD IT: have a look at the csv

#CORRECT SAMPLE SHEET FORMAT ==> REQUIRED COLUMNS: array and slide columns
SampleSheet <- read.csv("Input_Data/SampleSheet_Report_II.csv",header=T)

# sample of type of sample sheet that can actually be read by the loading function
SampleSheet 

# LOADING FUNCTION
# Set the directory in which the raw data are stored
# load the samplesheet using the function read.metharray.sheet
baseDir <- ("Input_data")
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
?RGChannelSet

######               STEP TWO              ######
###### Create the Red and Green dataframes ######
######                                     ######

# We extract the Green and Red Channels using the functions getGreen and getRed
# extract fluorescence intensity info associated to red
Red <- data.frame(getRed(RGset)) 
dim(Red)
head(Red)
# extract fluorescence intensity info associated to green
Green <- data.frame(getGreen(RGset))
dim(Green)
head(Green)

# why 622399 rows instead of 485k as expected? 
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
ProbeInfo <- data.frame(getProbeInfo(RGset))
# Strange...this df has XXXXXX rows, not 485512...why?
dim(ProbeInfo)
# It only stores Type I by default.
# The getProbeInfo function returns only Type I probes. However, this is not really specified in the help page. 
# Take home message: when you work with a package, try before you trust!
head(getProbeInfo(RGset, type = "II"))
ProbeInfo_II <- data.frame(getProbeInfo(RGset, type = "II"))
dim(ProbeInfo_II)
350036+135476 #485512

# CHECKING ADDRESS AND VERYFING IF IT IS A OR B

# From probes to Addresses
head(ProbeInfo)
head(ProbeInfo_II)

# I want to check the probe having address 18756452 in the Type I array: none
ProbeInfo[ProbeInfo$AddressA=="18756452",]
ProbeInfo[ProbeInfo$AddressB=="18756452",]

# there is not an AddressB_ID with code 18756452;
ProbeInfo_II[ProbeInfo_II$AddressB=="18756452",]
# there is an AddressA with code 18756452 and it is associated to the probe cg25192902, type II
ProbeInfo_II[ProbeInfo_II$AddressA=="18756452",]


# cg25192902, type II (the first probe in the Illumina450Manifest_clean object)
# Using this commands you can fill the table on the word document
Red[rownames(Red)=="18756452",]
Green[rownames(Green)=="18756452",]

# Are there out of band signals???

######          STEP FOUR         ######
###### Create the object MSet.raw ######
######                            ######

## CHUNK 4: Extract methylated and unmethylated signals using the function MSet.raw
MSet.raw <- preprocessRaw(RGset)
MSet.raw
# Note that now the number of rows is 485512, exactly the number of probes according to the Manifest!
save(MSet.raw,file="MSet_raw.RData")
Meth <- as.matrix(getMeth(MSet.raw))
str(Meth)
head(Meth)
Unmeth <- as.matrix(getUnmeth(MSet.raw))
str(Unmeth)
head(Unmeth)

# Let's check what happens to the probes that we considered before when we move from RGset to MethylSet
# cg25192902, type II (the first probe in the Illumina450Manifest_clean object)
Red[rownames(Red)=="18756452",]
Green[rownames(Green)=="18756452",]
Unmeth[rownames(Unmeth)=="cg25192902",]
Meth[rownames(Meth)=="cg25192902",]

#### OUT OF BAND SIGNALS: differences bw RGset and Methset
# cg25192902, type II 
Illumina450Manifest_clean[Illumina450Manifest_clean$IlmnID=="cg25192902",]
Red[rownames(Red)=="18756452",]
# But the two addresses are present also when I look at the Green object: these are out of band signals!
Green[rownames(Green)=="18756452",]
Unmeth[rownames(Unmeth)=="cg25192902",]
Meth[rownames(Meth)=="cg25192902",]


