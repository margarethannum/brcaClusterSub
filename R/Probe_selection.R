#####################################################
# Methylation Probe Selection using Identified DMRs #
# Author: Margaret Hannum                           #
#####################################################
setwd("C:/Users/hannumm/Documents/DMR Analysis")
load("tcga.meth.all.480.RData")

library(ChAMP)
library(wateRmelon)

tcga.meth.all.480 <- as.matrix(tcga.meth.all.480)
#Run filtering from ChAMP on all 480K probes in 10 samples (probe filtering does not require all samples and the file was too big)
tcga.meth.qc1 <- champ.filter(beta = tcga.meth.all.480[,1:10], pd = NULL, filterDetP = FALSE, filterBeads = FALSE, fixOutlier = FALSE)

#Apply filtered probes to full 885 samples
tcga.meth.qc <- subset(tcga.meth.all.480, rownames(tcga.meth.all.480) %in% rownames(tcga.meth.qc1$beta))
setwd("C:/Users/hannumm/Documents/DMR Analysis")
save(tcga.meth.qc, file = "tcga.meth.qc.RData") #This file is the 885 samples with 414206 cpgs filtered
load("tcga.meth.qc.RData")

#Get annotation for probe types
probetype <- subset(probe.features, 
                    rownames(probe.features) %in% rownames(tcga.meth.qc))
probetype$Type <- as.numeric(probetype$Type)
probetype <- probetype[ order(row.names(probetype)),]
setwd("C:/Users/hannumm/Documents/DMR Analysis")
save(probetype, file = "probetype.annot.RData") #This file is the cleaned probe types
load("probetype.annot.RData")
#tcga.meth.qc <- tcga.meth.qc[order(row.names(tcga.meth.qc)),]


begin.time = Sys.time() #28 seconds per sample
#Normalize: RUnning BMIQ Function on all 885 samples, 414K QC filtered probes
#norm.BetaVal=sapply(1:ncol(tcga.meth.qc), function(i){
 norm.BetaVal=sapply(1:3, function(i){
    
  print(i)
  BMIQ(tcga.meth.qc[,i], probetype$Type)$nbeta
})
row.names(norm.BetaVal)=row.names(tcga.meth.qc)
colnames(norm.BetaVal)=colnames(tcga.meth.qc)[1:3]
save(norm.BetaVal, file = "tcga.meth.qc.norm.RData") #File with 414K filtered, normalized probes in 885 samples
end.time = Sys.time()
end.time - begin.time

load("tcga.meth.qc.norm.RData")
tcga.meth.qc.df <- as.data.frame(tcga.meth.qc)
save(tcga.meth.qc.df, file = "tcga.meth.qc.df.RData")

library(data.table)
tcga.meth.qc.i <- fread("tcga.meth.100.RData", select = 1:2, header = TRUE)


#Trying to improve computing time
begin.time = Sys.time()
load("tcga.meth.100.RData")
load("tcga.meth.400.RData")
rownames(tcga.meth.f.500) <- tcga.meth.f.all[,1]
#tcga.meth.f.100 <- tcga.meth.f.all[,-1]
tcga.meth.f.500 <- subset(tcga.meth.f.500, rownames(tcga.meth.f.500) %in% rownames(probetype) )
end.time = Sys.time() #3.4 seconds

begin.time = Sys.time() #28 seconds per sample
#Normalize: RUnning BMIQ Function on all 885 samples, 414K QC filtered probes
#norm.BetaVal=sapply(1:ncol(tcga.meth.qc), function(i){
norm.BetaVal.5=sapply(1:100, function(i){
  
  print(i)
  BMIQ(tcga.meth.f.500[,i], probetype$Type)$nbeta
})
row.names(norm.BetaVal.5)=row.names(tcga.meth.f.500)
colnames(norm.BetaVal.5)=colnames(tcga.meth.f.500)[1:100]
end.time = Sys.time()
end.time - begin.time #27 seconds

#norm.BetaVal.200 <- cbind2(norm.BetaVal.2, norm.BetaVal)
save(norm.BetaVal.5, file = "tcga.meth.qc.norm.5.RData")

rm(list = ls())
load("probetype.annot.RData")
load("tcga.meth.100.RData")
load("tcga.meth.600.RData")
rownames(tcga.meth.f.600) <- tcga.meth.f.all[,1]
#tcga.meth.f.100 <- tcga.meth.f.all[,-1]
tcga.meth.f.600 <- subset(tcga.meth.f.600, rownames(tcga.meth.f.600) %in% rownames(probetype) )

#Normalize: RUnning BMIQ Function on all 885 samples, 414K QC filtered probes
#norm.BetaVal=sapply(1:ncol(tcga.meth.qc), function(i){
norm.BetaVal.6=sapply(1:100, function(i){
  
  print(i)
  BMIQ(tcga.meth.f.600[,i], probetype$Type)$nbeta
})
row.names(norm.BetaVal.6)=row.names(tcga.meth.f.600)
colnames(norm.BetaVal.6)=colnames(tcga.meth.f.600)[1:100]
save(norm.BetaVal.6, file = "tcga.meth.qc.norm.6.RData")

rm(list = ls())
load("probetype.annot.RData")
load("tcga.meth.100.RData")
load("tcga.meth.700.RData")
rownames(tcga.meth.f.700) <- tcga.meth.f.all[,1]
#tcga.meth.f.100 <- tcga.meth.f.all[,-1]
tcga.meth.f.700 <- subset(tcga.meth.f.700, rownames(tcga.meth.f.700) %in% rownames(probetype) )

#Normalize: RUnning BMIQ Function on all 885 samples, 414K QC filtered probes
#norm.BetaVal=sapply(1:ncol(tcga.meth.qc), function(i){
norm.BetaVal.7=sapply(1:100, function(i){
  
  print(i)
  BMIQ(tcga.meth.f.700[,i], probetype$Type)$nbeta
})
row.names(norm.BetaVal.7)=row.names(tcga.meth.f.700)
colnames(norm.BetaVal.7)=colnames(tcga.meth.f.700)[1:100]
save(norm.BetaVal.7, file = "tcga.meth.qc.norm.7.RData")

rm(list = ls())
load("probetype.annot.RData")
load("tcga.meth.100.RData")
load("tcga.meth.800.RData")
rownames(tcga.meth.f.800) <- tcga.meth.f.all[,1]
#tcga.meth.f.100 <- tcga.meth.f.all[,-1]
tcga.meth.f.800 <- subset(tcga.meth.f.800, rownames(tcga.meth.f.800) %in% rownames(probetype) )

#Normalize: RUnning BMIQ Function on all 885 samples, 414K QC filtered probes
#norm.BetaVal=sapply(1:ncol(tcga.meth.qc), function(i){
norm.BetaVal.8=sapply(1:100, function(i){
  
  print(i)
  BMIQ(tcga.meth.f.800[,i], probetype$Type)$nbeta
})
row.names(norm.BetaVal.8)=row.names(tcga.meth.f.800)
colnames(norm.BetaVal.8)=colnames(tcga.meth.f.800)[1:100]
save(norm.BetaVal.8, file = "tcga.meth.qc.norm.8.RData")

rm(list = ls())
load("probetype.annot.RData")
load("tcga.meth.100.RData")
load("tcga.meth.85.RData")
rownames(tcga.meth.f.85) <- tcga.meth.f.all[,1]
#tcga.meth.f.100 <- tcga.meth.f.all[,-1]
tcga.meth.f.85 <- subset(tcga.meth.f.85, rownames(tcga.meth.f.85) %in% rownames(probetype) )

#Normalize: RUnning BMIQ Function on all 885 samples, 414K QC filtered probes
#norm.BetaVal=sapply(1:ncol(tcga.meth.qc), function(i){
norm.BetaVal.9=sapply(1:85, function(i){
  
  print(i)
  BMIQ(tcga.meth.f.85[,i], probetype$Type)$nbeta
})
row.names(norm.BetaVal.9)=row.names(tcga.meth.f.85)
colnames(norm.BetaVal.9)=colnames(tcga.meth.f.85)[1:85]
save(norm.BetaVal.9, file = "tcga.meth.qc.norm.9.RData")

##Combine normalized set
#load("tcga.meth.qc.norm.9.RData") #Loaded all 9 files with separate 100 samples
norm.BetaVal.all <- cbind(norm.BetaVal.100, norm.BetaVal.2, norm.BetaVal.3, norm.BetaVal.4, norm.BetaVal.5, norm.BetaVal.6, norm.BetaVal.7, norm.BetaVal.8, norm.BetaVal.9)
save(norm.BetaVal.all, file = "tcga.meth.qc.norm.all.RData") #885 Samples with 414K clean, normalized probes

#List of probes that were identified in DMR process 
load("dmrcpg_with_gene.RData") #This file is the selected 69K probes with gene names and coordinates
out.full <- out.full[complete.cases(out.full[,2]),] #60677 that have gene symbol

#Subset of normalized probes that were in DMRs: 55035
norm.BetaVal.select <- subset(norm.BetaVal.all, rownames(norm.BetaVal.all) %in% out.full$cpg)
save(norm.BetaVal.select, file = "tcga.meth.qc.norm.DMR.select.RData")
load("tcga.meth.qc.norm.DMR.select.RData")

#Load TCGA gene expression data
setwd("H:/Biostatistics/Ronglai-Margaret/TCGABRCA")
load("H:/Biostatistics/Ronglai-Margaret/TCGABRCA/BRCA_all.Rdata")

#Annotation file removing extra genes after ";" symbol
out.full.split <- out.full
out.full.split$Gene_Symbol <- sub(";.*$", "", out.full.split$Gene_Symbol)

#Gene expression data
tcga.mrna <- t(BRCA_all$mrna)
#Split names since gene names are in front of | character
rownames(tcga.mrna) <- sapply(strsplit(rownames(tcga.mrna), "|", fixed = TRUE), '[[', 1)
#Subset gene expression and meth data with each other to match gene names
tcga.mrna.select <- subset(tcga.mrna, rownames(tcga.mrna) %in% out.full.split$Gene_Symbol) #4244 genes that are also in methylation data
out.full.select <- subset(out.full.split, out.full.split$Gene_Symbol %in% rownames(tcga.mrna.select)) #50K probes
tcga.meth.select <- subset(norm.BetaVal.select, rownames(norm.BetaVal.select) %in% out.full.select$cpg) #41K Probes
#Make methylation sample names match mrna sample names by removing end part of string
#colnames(tcga.meth.select) <- substr(colnames(tcga.meth.select),1,12)
#subset out.full.select based on cpgs in meth.select
out.full.select <- subset(out.full.select, out.full.select$cpg %in% rownames(tcga.meth.select))

#Subset methylation and expression based on sample names
#Have to transpose first since subset works on rows
tcga.mrna.select <- t(tcga.mrna.select)
tcga.meth.select <- t(tcga.meth.select)
tcga.mrna.select <- subset(tcga.mrna.select, rownames(tcga.mrna.select) %in% rownames(tcga.meth.select))
tcga.meth.select <- subset(tcga.meth.select, rownames(tcga.meth.select) %in% rownames(tcga.mrna.select))
###Matching sample names only left 663 samples??


tcga.meth.qc2 <- norm.BetaVal.select[-which(rowMeans(is.na(norm.BetaVal.select)) > 0.3),] #52094 cpg 885 sample
#Remove samples with sample-level coverage <= 0.95
tcga.meth.qc3 <- t(tcga.meth.qc2) #transpose to be able to do rowwise
tcga.meth.qc3 <-  tcga.meth.qc3[-which(is.na(tcga.meth.qc3) > 0.05),] #52094 cpg 884 sample (1 sample removed)
colnames(tcga.meth.qc3) <- substr(colnames(tcga.meth.qc3),1,16)
#Remove normal samples (those without "01A")
tcga.meth.qc4 <- c(grep(c("01A"), rownames(tcga.meth.qc3), fixed = TRUE, value = TRUE))
tcga.meth.select <- subset(tcga.meth.qc3, rownames(tcga.meth.qc3) %in% c(tcga.meth.qc4)) 
rownames(tcga.meth.select) <- substr(rownames(tcga.meth.select),1,12)
tcga.mrna.select <- subset(tcga.mrna.select, rownames(tcga.mrna.select) %in% rownames(tcga.meth.select))
tcga.meth.select <- subset(tcga.meth.select, rownames(tcga.meth.select) %in% rownames(tcga.mrna.select))
tcga.mrna.select <- t(tcga.mrna.select)
tcga.meth.select <- t(tcga.meth.select)
tcga.meth.select <- subset(tcga.meth.select, rownames(tcga.meth.select) %in% out.full.select$cpg)
out.full.select <- subset(out.full, out.full$cpg %in% rownames(tcga.meth.select)) #45K probes
tcga.mrna.select <- subset(tcga.mrna.select, rownames(tcga.mrna.select) %in% out.full.select$Gene_Symbol) #3816 genes that are also in methylation data
#3794 genes in 39290 probes in 646 samples

id <- NULL
rho <- NULL
sd <- NULL
for (i in 1:length(unique(out.full.select$Gene_Symbol))){
  #for (i in 1:100){
  
  gene <- unique(out.full.select$Gene_Symbol)[i]
  subanno <- subset(out.full.select, out.full.select$Gene_Symbol == gene)
  #CG PROBE ID FROM SUBSET use it to extract rows/cols methylation data
  meth.mat.i <- tcga.meth.select[subanno$cpg, , drop = F]
  exp.mat.i <- tcga.mrna.select[gene, , drop=F]
  rho.i <- cor(t(exp.mat.i), t(meth.mat.i), use="pairwise.complete.obs")
  sd.i <- sd(rho.i)
  sd <- c(sd, sd.i)
  rho <- c(rho, rho.i)
  id.i <- paste(colnames(rho.i), rownames(rho.i), sep=".")
  #id.i <- do.call(paste, c(subanno[1:2], sep = "."))
  id <- c(id, id.i)
}

names(rho) <- id
names(sd) <- unique(out.full.select$Gene_Symbol)
summary(sd)


setwd("H:/Biostatistics/Ronglai-Margaret/DMR Analysis/Output")
save(rho, file = "rho.RData")

rho.df <- as.data.frame(rho)

hist(rho.df$rho, main = "Histogram of Correlation")
summary(rho.df)

#Subset of cpgs with negative correlation less than -0.3
rho.df.1 <- subset(rho.df, rho.df$rho < -0.3) #4579 cpgs
