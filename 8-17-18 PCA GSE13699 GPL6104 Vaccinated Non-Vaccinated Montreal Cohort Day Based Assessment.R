###################################################
### Load all necessary packages
###################################################
source("https://bioconductor.org/biocLite.R") # works

biocLite("GEOquery")
biocLite("limma")
biocLite("Biobase")
biocLite("affy")
library(GEOquery)
library(Biobase)
library(limma)
library(affy)
library(knitr)
library(rmarkdown)
library(GEOquery)
library(stringr)
library(tibble)
library(stringr)

###################################################
### Clean up the datasets annotations
###################################################
# 1. Access the GPL6104 dataset
gset <- getGEO("GSE13699", GSEMatrix =TRUE)# 1 item length list
names(gset) # file name for now gset is a list with one element that was derived from file GSE46268_series_matrix.txt.gz.
if (length(gset) > 1) idx <- grep("GPL6104", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# 2. Assess the features portion of the datasets
fData(gset)[1:5,1:5]
colnames(fData(gset))
all(rownames(fData(gset)) == rownames(exprs(gset))) # SHOULD ALWAYS CHECK TO MAKE SURE ITS THE CASE!

# 3. Convert the phenotype data set into a data frame
gsms <- pData(gset)$geo_accession
pdata <- pData(gset) # access the phenotype/annotations of the datasets
pdata <- data.frame(pdata) # convert the data into a data frame. Makes is easier to manage!

pdata$title <- gsub(".*;","",pdata$title) # title column is now cohort column
colnames(pdata)[which(names(pdata) == "title")] <- "cohort"

pdata <- cbind(pdata,pdata[,10])  # duplicate characteristics_ch1 column
colnames(pdata)[which(names(pdata) == "pdata[, 10]")] <- "Age"

pdata <- cbind(pdata,pdata[,10])  # duplicate characteristics_ch1 column
colnames(pdata)[which(names(pdata) == "pdata[, 10]")] <- "Gender"

pdata <- cbind(pdata,pdata[,10])  # duplicate characteristics_ch1 column
colnames(pdata)[which(names(pdata) == "pdata[, 10]")] <- "Volunteer_Status"

pdata <- cbind(pdata,pdata[,10])  # duplicate characteristics_ch1 column
colnames(pdata)[which(names(pdata) == "pdata[, 10]")] <- "ID"

pdata$Age <- gsub(".*Age: ","",pdata$Age)
pdata$Age <- gsub(".*Age:","",pdata$Age)
pdata$Age <- gsub("^","Age",pdata$Age)

pdata$Volunteer_Status <- gsub("Vaccinated Volunteer.*","",pdata$Volunteer_Status)
pdata$Volunteer_Status <- NULL
pdata <- cbind(pdata,pdata[,10])  # duplicate characteristics_ch1 column
colnames(pdata)[which(names(pdata) == "pdata[, 10]")] <- "Volunteer_Status" # change column name given column name

pdata$Volunteer_Status <- gsub("YF.*","",pdata$Volunteer_Status)
pdata$Volunteer_Status <- sub("v","V",pdata$Volunteer_Status)
pdata$Volunteer_Status <- sub(" :","",pdata$Volunteer_Status)
pdata$Volunteer_Status <- sub(" ","",pdata$Volunteer_Status)
pdata$Volunteer_Status <- sub(" ","",pdata$Volunteer_Status)
pdata$Volunteer_Status <- sub("VaccinatedVolunteer","VV",pdata$Volunteer_Status)
pdata$Volunteer_Status <- sub("NonVV","NV",pdata$Volunteer_Status)

pdata$Gender <- NULL

pdata <- cbind(pdata,pdata[,10])  # duplicate characteristics_ch1 column
colnames(pdata)[which(names(pdata) == "pdata[, 10]")] <- "Gender"

pdata$Gender <- gsub("^.*?Gender","Gender",pdata$Gender) # remove before a string not including the string
pdata$Gender <- gsub("Age.*","",pdata$Gender) 
pdata$Gender <- gsub(":","",pdata$Gender) 
pdata$Gender <- gsub(";","",pdata$Gender) 
pdata$Gender <- gsub(".*Gender","",pdata$Gender) 
pdata$Gender <- gsub(" ","",pdata$Gender)

pdata$ID <- gsub("^.*?YF","YF",pdata$ID) # remove before a string not including the string
pdata$ID <- gsub(";.*","",pdata$ID) # remove all after a character
pdata$Volunteer_Status <- gsub(" ","",pdata$Volunteer_Status)

pdata <- cbind(pdata,pdata[,8])  # duplicate characteristics_ch1 column
colnames(pdata)[which(names(pdata) == "pdata[, 8]")] <- "Day" # change column name given column name

pdata$Day <- gsub("^.*?Day","Day",pdata$Day)
pdata$Day <- gsub(" ","",pdata$Day)


pdata <- cbind(pdata,pdata[,37])  # duplicate characteristics_ch1 column
colnames(pdata)[which(names(pdata) == "pdata[, 37]")] <- "Age_Group"
pdata$Age_Group <- gsub("^Age2","Twenties",pdata$Age_Group)
pdata$Age_Group <- gsub("^Age3","Thirties",pdata$Age_Group)
pdata$Age_Group <- gsub("^Age4","Forties",pdata$Age_Group)
pdata$Age_Group <- gsub("^Age1","Teens",pdata$Age_Group)
pdata$Age_Group <- gsub("s.*","",pdata$Age_Group)



pdata <- pdata[grep("Montreal Cohort",pdata$cohort),] #filter based on cohort

###################################################
### Filter Expression Dataset as necessary
###################################################

to.remove <- pdata[[2]] # contains GSMs of 100 vaccianted samples
x <- as.character(to.remove) # change the type to make it work with the function
ex <- exprs(gset) # matrix 
ex <- data.frame(ex) # data frame
new_ex <- ex[x] # data frame with 100 samples (+ vaccinated)
gf <- new_ex
gf <- na.omit(gf) # remove rows w/ NA in expression dataset 100 samples!
dim(gf) # 18563 probes 88 samples!

###################################################
### PCA Assessment
###################################################

fl <- pdata$Day # 
fl <- as.factor(fl)


design <- model.matrix(~ Day + 0, pdata) # accepts a formula, then a list (gset)
design 
# If row is 1, then it belongs to said group (Gn), if not, then row is 0
colnames(design) <- levels(fl)
fit <- lmFit(gf, design) # fitting a linear model to the expression matrix (rows as genes and samples as columns) given a PCA classification!
# list object 


# NEED TO CHANGE NAMES BELOW
cont.matrix <- makeContrasts(Day14-Day3,
                             Day14-Day7,
                             Day14-Day10,
                             Day14-Day28,
                             Day14-Day60,
                             Day14-Day0,
                             levels=design) # you can do three different contrasts and the tpes are within the design vector
# result is a matrix
# Construct the contrast matrix corresponding to specified contrasts of a set of parameters.



PC=prcomp(t(gf)) # class prcomp 
# principle components are character or numeric 
scores = predict(PC) # matrix
pc1 <- scores[ ,1] # first sample
pc2 <- scores[ ,2] # second sample!
shape <- as.numeric(fl) + 14 # numeric
plot(pc1, pc2, col=fl, pch=shape, cex=2)
legend("bottomright", pch=unique(shape), paste(unique(fl)))

# Created a PCA PLOT OF drawl days for the vaccinated samples! MUCH BETTER CORRELATIONS!
