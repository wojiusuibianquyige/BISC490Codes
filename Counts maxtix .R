library(TCGAbiolinks)
library(R.utils)
#access TCGA data 
query <- GDCquery(project = "TCGA-GBM",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")
GDCdownload(query = query)
samplesDown <- getResults(query,cols=c("cases"))
options(stringsAsFactors = F)

#Move all files under "sampleFiles_GBM"
setwd("/Users/yiyang/Desktop/我太难了/BISC490/Data")
dir.create("sampleFiles_GBM")
filepath <- dir(path = "./RAWData", full.names = T)
for(wd in filepath){
    files <- dir(path = wd, pattern = "gz$")
    fromfilepath <- paste(wd, "/", files, sep = "")
    tofilepath <- paste("./sampleFiles_GBM/", files, sep = "")
    file.copy(fromfilepath, tofilepath)
}

#unzip all files and delete the originial ones
setwd("./sampleFiles_GBM")
countsFiles <- dir(path = "./", pattern = "gz$")
sapply(countsFiles, gunzip)


library(rjson)
library(dplyr)
library(limma)
library(stringr)
# json files
setwd("/Users/yiyang/Desktop/我太难了/BISC490/Data/sampleFiles_GBM")
jsonFile <- fromJSON(file = "/Users/yiyang/Desktop/我太难了/BISC490/TCGA-GBM/metadata.cart.2020-09-08.json")
filesNameToBarcode <- data.frame(filesName = c(), TCGA_Barcode = c())
for(i in 1:length(jsonFile)){
    TCGA_Barcode <- jsonFile[[i]][["associated_entities"]][[1]][["entity_submitter_id"]]
    file_name <- jsonFile[[i]][["file_name"]]
    filesNameToBarcode <- rbind(filesNameToBarcode, data.frame(filesName = file_name, TCGA_Barcode = TCGA_Barcode))
}
rownames(filesNameToBarcode) <- filesNameToBarcode[,1]

library(maftools)
#get counts matrix 
filesNameToBarcode <- filesNameToBarcode[-1]
setwd("/Users/yiyang/Desktop/我太难了/BISC490/Data/sampleFiles_GBM")
countsFileNames <- dir(pattern = "counts$")
allsampleRawCounts <- data.frame()
for(txtFile in countsFileNames){
    SampleCounts <- read.table(txtFile, header = F)
    rownames(SampleCounts) <- SampleCounts[,1]
    SampleCounts <- SampleCounts[-1]
    colnames(SampleCounts) == filesNameToBarcode$TCGA_Barcode
    if (dim(allsampleRawCounts)[1]==0){
        allsampleRawCounts <- SampleCounts
    }
    else{
        allsampleRawCounts <- cbind(allsampleRawCounts, SampleCounts)
    }
}
colnames(allsampleRawCounts) <- filesNameToBarcode$TCGA_Barcode
ensembl_id <- substr(row.names(allsampleRawCounts), 1, 15)
rownames(allsampleRawCounts) <- ensembl_id
write.csv(allsampleRawCounts, file="/Users/yiyang/Desktop/我太难了/BISC490/Data/Files/RawCounts.csv")
