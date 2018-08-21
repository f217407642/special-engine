#############################################
### Load Your Workspace
#############################################
source("https://bioconductor.org/biocLite.R") # works
source("http://www.bioconductor.org/biocLite.R") # does not work
library(GEOquery)

#############################################
### Access the data sets
#############################################
mydata <- getGEO("GSE13699")
mydata # list 2 items
mydata[[1]] # access first item in the list
# Alternative to immediate code above:
# acquire ftp link from right clicking the download link on the GEO page 
# (to the file: GSE13699-GPL6104_series_matrix.txt.gz)
setwd("your/favorite/location")

destfile <- "GSE data/GSE13699-GPL6104_series_matrix.txt.gz"
if(!file.exists(destfile)) download.file(remotefile, destfile) # checks to see if data is already downloaded and saved to directory
mydata<- getGEO(filename="your/favorite/location/GSE13699-GPL6104_series_matrix.txt.gz")
mydata # list
length(mydata) #  1 item in list 
dim(mydata) # 22,184 probes 126 samples in GPL6104 platform

#############################################
### Check out your data set
#############################################
exprs(mydata)[1:5,1:5] # 5 by 5 matrix 
summary(exprs(mydata)[,1:5]) # summary of 5 by 5 matrix 
boxplot(log2(exprs(mydata))) # check the quality of the probeset with boxplot, useless you have many samples!

#############################################
### Perform the array quality metric
#############################################
biocLite("arrayQualityMetrics")
library(arrayQualityMetrics) # Do a quality check with arrayQualityMetrics 
arrayQualityMetrics(mydata, force=TRUE) # force=TRUE argument will make sure that the output gets written if the specified directory name already exists
# Normalized! To have it raw, set do.logtransform=TRUE in the function 
# warned there were ALOT OF ISSUES
# by default saves at the location of your current directory 

#############################################
### Check out your expressio data set
#############################################
subset <- mydata[,1:4] # even though mydata is a list, you can select 
# still a list
dim(exprs(subset)) # dimensions of the subset: 22184 rows, 4 columns # expression data 
# OR
dim(subset) # 22184 probes, 4 samples

#############################################
### Check out your pheno data set
#############################################
dim(pData(subset)) # dimensions of the subset: 4 rows, 36 columns 
# data on the project and its people, 
colnames(pData(subset))

#############################################
### Check out your features data set
#############################################
dim(fData(subset)) # dimensions of the subset: 22184 rows, 26 columns 
# gene name information
colnames(fData(mydata)) # column names for fData
# or
colnames(features)
features <- fData(mydata) # save features data to a new variable
# subset <- mydata[1:100,] # filter by first 100 genes
all(rownames(fData(mydata)) == rownames(exprs(mydata))) # Should be true, because no changes made
#############################################
### Check for outlier arrays
#############################################
matrix <- exprs(mydata) # matrix!
outliers(matrix, method = c("KS"))# choose from either KS, sum, or upperquartile to find the outliers (outputs GSM)
outliers(matrix, method = c("sum"))
outliers(matrix, method = c("upperquartile"))



# We refer to the following github guide: https://github.com/bioinformatics-core-shared-training/microarray-analysis/blob/master/public-data.Rmd
# Noticed you must refer to GEO data set information to find the age of the patient's sample, sex, etc. fData and pData functions don't provide you with that!
# I was able to find outliers within the whole GSE13699 dataset though, do not know how to feasibly parse the outliers from the matrix w/o using python
# Have to decide how to definatively classify outliers! Say pass 3/3 outlier tests (KS,sum,upperquartile)
# Notice, this is just for the 