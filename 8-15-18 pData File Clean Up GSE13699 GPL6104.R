#############################################
### Load your workspace
#############################################
source("https://bioconductor.org/biocLite.R") # works
#biocLite("GEOquery")
#biocLite("limma")
#biocLite("Biobase")
#biocLite("affy")
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
#############################################
### Clean up the annotation file
#############################################

# 1. Access the GPL6104 dataset
gset <- getGEO("GSE13699", GSEMatrix =TRUE)# 1 item length list
names(gset) # file name for now gset is a list with one element that was derived from file GSE46268_series_matrix.txt.gz.
if (length(gset) > 1) idx <- grep("GPL6104", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]] # gset is now annotations for a single platform from the GSE3699 study

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

pdata$Volunteer_Status <- gsub("YF.*","",pdata$Volunteer_Status) # remove characters after and including YF string
pdata$Volunteer_Status <- sub("v","V",pdata$Volunteer_Status) # substitute v character with V
pdata$Volunteer_Status <- sub(" :","",pdata$Volunteer_Status) # remove :
pdata$Volunteer_Status <- sub(" ","",pdata$Volunteer_Status) # remove single empty space
pdata$Volunteer_Status <- sub(" ","",pdata$Volunteer_Status) # remove single empty space (gsub function removes all empty strings)
pdata$Volunteer_Status <- sub("VaccinatedVolunteer","VV",pdata$Volunteer_Status) # make string substitution 
pdata$Volunteer_Status <- sub("NonVV","NV",pdata$Volunteer_Status) # make string substitution 

pdata$Gender <- NULL # delete the Gender column for now 

pdata <- cbind(pdata,pdata[,10])  # duplicate characteristics_ch1 column
colnames(pdata)[which(names(pdata) == "pdata[, 10]")] <- "Gender"

pdata$Gender <- gsub("^.*?Gender","Gender",pdata$Gender) # remove before a string not including the string
pdata$Gender <- gsub("Age.*","",pdata$Gender) 
pdata$Gender <- gsub(":","",pdata$Gender) 
pdata$Gender <- gsub(";","",pdata$Gender) 
pdata$Gender <- gsub(".*Gender","",pdata$Gender) # remove characters before and including the Gender string 
pdata$Gender <- gsub(" ","",pdata$Gender) # remove all empty spaces in the indicated column

pdata$ID <- gsub("^.*?YF","YF",pdata$ID) # remove before a string not including the string
pdata$ID <- gsub(";.*","",pdata$ID) # remove all characters after a given character
pdata$Volunteer_Status <- gsub(" ","",pdata$Volunteer_Status) # remove all spaces

pdata <- cbind(pdata,pdata[,8])  # duplicate characteristics_ch1 column
colnames(pdata)[which(names(pdata) == "pdata[, 8]")] <- "Day" # change column name given column name

pdata$Day <- gsub(".*Day","",pdata$Day)
pdata$Day <- gsub(" ","",pdata$Day)
pdata$Age <- gsub(".*e","",pdata$Age)


pdata <- cbind(pdata,pdata[,37])  # duplicate a column
colnames(pdata)[which(names(pdata) == "pdata[, 37]")] <- "Age_Group"

pdata$Age_Group <- gsub("^2","Twenties",pdata$Age_Group) # substitute a string for another string
pdata$Age_Group <- gsub("^3","Thirties",pdata$Age_Group)
pdata$Age_Group <- gsub("^4","Forties",pdata$Age_Group)
pdata$Age_Group <- gsub("^1","Teens",pdata$Age_Group)
pdata$Age_Group <- gsub("s.*","s",pdata$Age_Group)



pdata <- cbind(pdata,pdata[,19])  # duplicate 
colnames(pdata)[which(names(pdata) == "pdata[, 19]")] <- "Vaccine" # change column name given column name
pdata$Vaccine <- gsub(".*with","",pdata$Vaccine)
pdata$Vaccine <- gsub(".*the","",pdata$Vaccine)
pdata$Vaccine <- gsub(",.*","",pdata$Vaccine)
pdata$Vaccine <- gsub("Stamail","Stamaril",pdata$Vaccine)
pdata$Volunteer_Status <- gsub("VV","Vaccinated",pdata$Volunteer_Status)
pdata$Volunteer_Status <- gsub("NV","Non-Vaccinated",pdata$Volunteer_Status)


#############################################
### Create a data frame that includes 
### pertinent info
#############################################

humoral_response <- data.frame(
                      "ID" = c("YF001","YF002","YF003",
                      "YF004","YF005","YF006",
                      "YF007","YF008","YF009",
                      "YF010","YF011","YF012",
                      "YF013","YF014","YF015",
                      "YF016","YF017","YF018",
                      "YF019","YF020"),
                      "DPV_14" = c(1/640,1/640,1/80,
                                1/320,1/320,1/320,
                                1/160,0,1/640,
                                1/2560,1/1280,1/640,
                                1/640,1/1280,1/320,
                                1/640,1/160,1/640,
                                1/160,1/640),
                      "DPV_28" = c(1/640,1/640,1/320,
                                 1/320,1/160,1/320,
                                 1/320,1/1280,1/640,
                                 1/2560,1/1280,1/640,
                                 1/640,1/1280,1/320,
                                 1/640,1/160,1/640,
                                 1/160,1/640),
                      "DPV_60" = c(1/1280,1/320,1/80,
                                 1/320,1/160,1/80,
                                 1/320,1/1280,1/640,
                                 NA,1/1280,1/640,
                                 1/640,1/1280,1/320,
                                 1/640,1/160,1/640,
                                 1/160,1/640),
                      "DPV_90" = c(1280,320,80,
                                 1/640,1/320,1/80,
                                 1/320,1/1280,1/640,
                                 NA,1/1280,1/640,
                                 1/640,1/1280,1/320,
                                 1/640,1/160,1/640,
                                 1/160,1/640),
                      "DPV_180" = c(1/640,1/80,1/80,
                                  1/640,1/320,1/80,
                                  1/160,1/640,1/1280,
                                  1/NA,1/340,1/640,
                                  1/640,1/640,1/160,
                                  1/640,1/80,1/640,
                                  1/160,1/320),
                      "DPV_360" = c(1/640,1/80,1/80,
                                  1/640,1/160,1/640,
                                  1/160,1/640,1/1280,
                                  NA,1/320,1/320,
                                  1/640,1/640,1/160,
                                  1/640,1/40,1/320,
                                  1/160,1/320)
                      )
#############################################
### Merge the two datasets
#############################################

pdata <- merge(pdata,humoral_response,all = TRUE) # all set to true so that no rows are deleted from original pdata frame


#############################################
### Write the merged dataset to a file
#############################################

write.csv(pdata,file="GSE13699_GPL6104_pData.csv")

###################################################
### Summary
###################################################

# pdata has annotations for 126 different samples. Columnns of interest are named as below:
# Age_Group: one of 20s,30s,40s,or teens
# Day: one of 0,7,3,14,28,60, or 10
# Gender: male (m) or female (f)
# Volunteer_Status: vaccinated (VV) or non-vaccinated (NV)
# Age: numerical value
# ID: sample ID begins with YF
# Cohort: Montreal or Lausanne cohort
# Vaccine: One of Stamaril or YF-Vax
# DPV_14: Humoral immune response values for day post vaccination (DPV)
# DPV_28: Humoral immune response values for day post vaccination (DPV)
# DPV_60: Humoral immune response values for day post vaccination (DPV)
# DPV_90: Humoral immune response values for day post vaccination (DPV)
# DPV_180: Humoral immune response values for day post vaccination (DPV)
# DPV_360: Humoral immune response values for day post vaccination (DPV)