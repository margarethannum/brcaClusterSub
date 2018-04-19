################################################################
## TCGA BRCA Differentially Methylated Region Identification  ##
## Author: Margaret Hannum                                    ##
## Date: 6/21/2017                                            ##
################################################################

#Install ChAMP package from bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite("ChAMP")

source("https://bioconductor.org/biocLite.R")
biocLite("minfi")

source("https://bioconductor.org/biocLite.R")
biocLite("FDb.InfiniumMethylation.hg19")
biocLite("lumi")

library(utils)
library(ChAMP)
library(ChAMPdata)
library(tidyverse)
library(data.table)
library(minfi)
library(lumi)
library(wateRmelon)

### Using test data from ChAMP ###
testDir=system.file("extdata",package="ChAMPdata")
myLoad <- champ.load(testDir,arraytype="450K")
champ.process(directory = testDir)
#Testing Champ functions with in-package test data
champ.QC(beta = myLoad$beta, pheno = myLoad$pd$Sample_Group, 
         mdsPlot = TRUE, densityPlot = TRUE, dendrogram = TRUE, 
         PDFplot = TRUE, Rplot = TRUE, Feature.sel = "None", 
         resultsDir = "./CHAMP_QCimages/")
#Normalize
myNorm <- champ.norm(beta = myLoad$beta, method = "BMIQ", plotBMIQ = FALSE, arraytype = "450K")
detectCores()
#I kept getting error message below on champ.DMR package because the computer
#would not allow me to use more than one core. So I set cores = 1 and it ran fine.
myDMR <- champ.DMR(beta=myNorm,pheno=myLoad$pd$Sample_Group,method="DMRcate", cores = 1)
DMR.GUI(DMR = myDMR, beta = myNorm, pheno = myLoad$pd$Sample_Group, runDMP = TRUE, compare.group = NULL, arraytype = "450K")


#load pre-processed BRCA TCGA DNA methylation data
#90 tumor with matched normal-adjacent
setwd("H:/Biostatistics/Ronglai-Margaret/TCGABRCA")
load("H:/Biostatistics/Ronglai-Margaret/TCGABRCA/20160511_BRCA_raw_BetaVal_cover.RData")
#Load methylation data from Ronglai
#tmpdir <- tempdir()
#dir.create("DNAMeth.TCGA")
#untar("gdac.broadinstitute.org_BRCA.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2016012800.0.0.tar.gz", exdir = "DNAMeth.TCGA")
#system("tar -zxvf gdac.broadinstitute.org_BRCA.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2016012800.0.0.tar.gz")
#system("file gdac.broadinstitute.org_BRCA.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2016012800.0.0.tar.gz")

#where to save things later
setwd("C:/Users/hannumm/Documents/DMR Analysis")

#make phenotype vector which is necessary for champ.DMR function
sample.group <- c(replicate(90, "Tumor"),replicate(90,"Normal"))

#champ.DMR function applies bumphunter, DMRcate or ProbeLasso algorithms to detect DMRs
#beta "methylation beta valued dataset"- normalized. 
#pheno: categorical vector with phenotypes to be analyzed (cancer/normal)
#arraytype= "450K"
#method = "DMRcate"
#minProbes default is 7. Shuang used 3***
#adjPvalDmr default is 0.05
#Parameters specific for DMRcate algorithm:
  #rmSNPCH filters values by distance to SNP (default = TRUE) ?
  #fdr cutoff ?what should it be?
  #dist - max dist from CpG to SNP of probes to be filtered out (def = 2)?
  #lambda - kernel bandwidth (default = 1000) ?
  #C - scaling factor for bandwidth (default =2)?
#Value will be data.frame DMRcateDMR - suitable for champ.GSEA() analysis


## Using CHAMP pipeline with TCGA BRCA data provided by Ya, pre-processed ##
tcgaNorm <- as.matrix(raw.BetaVal.cover)
#DMRCate method
begin.time = Sys.time()
tcgaDMRcate <- champ.DMR(beta=tcgaNorm, pheno = sample.group, method="DMRcate", 
                     cores = 1, fdr = 0.05, minProbes = 7) #when minProbes=4, you get 13K+ regions!
#FDR cuttoff uses Benjamini-Hochberg procedure
#DMRcate method produced 6030 DMRs. 
end.time = Sys.time()
DMcatetime = end.time - begin.time #DMRcate method took 4.17 minutes
DMR.GUI(DMR = tcgaDMRcate, beta = tcgaNorm, pheno = sample.group, runDMP = TRUE, compare.group = NULL, arraytype = "450K")
#FDRCutoff 1%
tcgaDMRcate.1 <- champ.DMR(beta=tcgaNorm, pheno = sample.group, method="DMRcate", 
                         cores = 1, fdr = 0.01, minProbes = 7)
range(tcgaDMRcate.1$DMRcateDMR$maxbetafc)
DMR.GUI(DMR = tcgaDMRcate, beta = tcgaNorm, pheno = sample.group, runDMP = TRUE, compare.group = NULL, arraytype = "450K")
#FDRCutoff 1%

#Gene Set Enrichment Analysis with DMRcate results
myDMP <- champ.DMP(beta = tcgaNorm, pheno = sample.group)
dmcateGSEA <- champ.GSEA(beta = tcgaNorm, DMP = myDMP, DMR = tcgaDMRcate, arraytype = "450K", adjPval = 0.05)

#Bumphunter method
#begin.time = Sys.time()
#tcgaDMRBH <- champ.DMR(beta=tcgaNorm, pheno = sample.group,method="Bumphunter", 
 #                    cores = 1, maxGap = 1000, B = 250) #Default maxGap = 300, DMRcate paper uses 1000 to compare
#end.time = Sys.time() #10 resampling takes 3.09 min. 250 resampling takes 37.5 min.
#BHtime = end.time - begin.time
#bumphunter produced 1368 DMRs using 10 resample. 3113 DMRs using 250 resampling.
#DMR.GUI(DMR = tcgaDMR, beta = tcgaNorm, pheno = sample.group, runDMP = TRUE, compare.group = NULL, arraytype = "450K")

#Probe Lasso method
#begin.time = Sys.time()
#tcgaDMRPL <- champ.DMR(beta=tcgaNorm, pheno = sample.group,method="ProbeLasso", meanLassoRadius = 1000, cores = 1)
#end.time = Sys.time() #Method took 2.7 minutes with meanLasso = 375, 1.52 mins with meanLassoRadius = 1000
#PLtime = end.time - begin.time
#Probe Lasso discovered 327 DMRs with meanLassoRadius =375, 1401 DMRs with meanLassoRadius=1000.
#DMR.GUI(DMR = tcgaDMR, beta = tcgaNorm, pheno = sample.group, runDMP = TRUE, compare.group = NULL, arraytype = "450K")

save(tcgaDMR, file = "PL.tcgaDMR.RData")
save.image()
setwd("C:/Users/hannumm/Documents/DMR Analysis")
load("PL.tcgaDMR.RData")
mean(tcgaDMRPL$ProbeLassoDMR$width)
range(tcgaDMRPL$ProbeLassoDMR$width)
mean(tcgaDMRPL$ProbeLassoDMR$no.cpgs)
range(tcgaDMRcate$DMRcateDMR$no.cpgs)
load(".RData")

setwd("C:/Users/hannumm/Documents/R.Images")
save.image(file = "tcga.dmrID.RData")

load("tcga.dmrID.RData")

library(DMRcate)

#DMRcate results analysis
mean(tcgaDMRcate$DMRcateDMR$width)
range(tcgaDMRcate$DMRcateDMR$width)
mean(tcgaDMRcate$DMRcateDMR$no.cpgs)
range(tcgaDMRcate$DMRcateDMR$no.cpgs)

sum(tcgaDMRcate$DMRcateDMR$meanbetafc > 0)

#####################################
#Use full sample of data, apply DMRs identified using DMRcate
#####################################
setwd("H:/Biostatistics/Ronglai-Margaret/TCGABRCA")
setwd("C:/Users/hannumm/Documents/DMR Analysis") #Use file saved in local computer for faster loading times
#First 5 columns of full sample, which gives cg id, gene name, and coordinates
tcga.meth.genename <- fread("BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt", 
                     header = TRUE, skip = 1, select = 1:5, data.table = FALSE)



#Attempting to read in all columns, every fourth column, for one row. 0.309 seconds per cpg site on network data, 0.108 sec per cpg site on local computer file
begin.time = Sys.time()
tcga.meth.f <- fread("BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt", 
                     header = TRUE, skip = 1, nrows = 1, select = 1:3541, data.table = FALSE)[c(TRUE, rep(c(TRUE, FALSE, FALSE, FALSE),885))]
colnames(tcga.meth.f) <- fread("BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt", 
                               header = FALSE, skip = 0, nrows = 1, select = 1:3541, data.table = FALSE)[c(TRUE, rep(c(TRUE, FALSE, FALSE, FALSE),885))]
end.time = Sys.time()
end.time - begin.time


#Attempting to read in all 885 samples, every fourth column (Beta value), for all rows (485K probes) 
#Split it into 8 100 samples and 1 85 sample
setwd("C:/Users/hannumm/Documents/DMR Analysis")
begin.time = Sys.time()
tcga.meth.f.100 <- fread("BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt", 
                     header = TRUE, 
                     skip = 1, #this was done to preserve the correct classes, otherwise fread interprets everything as factors
                     #select = 1:3541, #total number of columns in the data set
                     select = 1:401,
                     #data.table = FALSE)[c(TRUE, rep(c(TRUE, FALSE, FALSE, FALSE),885))]  #To make subset of betas (First column is cpg names. Subsequent Column order was originally Beta, Gene Symbol, Chromosome, Gene Coordinate for each sample, so I made last three FALSE.)
                    data.table = FALSE)[c(TRUE, rep(c(TRUE, FALSE, FALSE, FALSE),100))]
colnames(tcga.meth.f.100) <- fread("BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt", 
                               header = FALSE, 
                               skip = 0, 
                               nrows = 1, #The first row is the true sample names. Above, I had to remove the row since the second row was also character but I wanted the columns to be interpreted as numeric
                             #  select = 1:3541, 
                               select = 1:401,
                             #  data.table = FALSE)[c(TRUE, rep(c(TRUE, FALSE, FALSE, FALSE),885))]
                             data.table = FALSE)[c(TRUE, rep(c(TRUE, FALSE, FALSE, FALSE),100))]
save(tcga.meth.f.100, file = "tcga.meth.100.RData")
end.time = Sys.time()
readtime <- end.time - begin.time #26 minutes to read in 100 samples, all 480K cpg sites
save(readtime, file = "readtime.RData")
load("readtime.RData")

#Read columns from subsequent ranges
tcga.meth.f.200 <- fread("BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt", 
                         header = TRUE, 
                         skip = 1, #this was done to preserve the correct classes, otherwise fread interprets everything as factors
                         #select = 1:3541, #total number of columns in the data set
                         select = 402:801,
                         #data.table = FALSE)[c(TRUE, rep(c(TRUE, FALSE, FALSE, FALSE),885))]  #To make subset of betas (First column is cpg names. Subsequent Column order was originally Beta, Gene Symbol, Chromosome, Gene Coordinate for each sample, so I made last three FALSE.)
                         data.table = FALSE)[c(rep(c(TRUE, FALSE, FALSE, FALSE),100))]
colnames(tcga.meth.f.200) <- fread("BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt", 
                                   header = FALSE, 
                                   skip = 0, 
                                   nrows = 1, #The first row is the true sample names. Above, I had to remove the row since the second row was also character but I wanted the columns to be interpreted as numeric
                                   #  select = 1:3541, 
                                   select = 402:801,
                                   #  data.table = FALSE)[c(TRUE, rep(c(TRUE, FALSE, FALSE, FALSE),885))]
                                   data.table = FALSE)[c(rep(c(TRUE, FALSE, FALSE, FALSE),100))]
save(tcga.meth.f.200, file = "tcga.meth.200.RData")
tcga.meth.f.300 <- fread("BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt", 
                         header = TRUE, 
                         skip = 1, #this was done to preserve the correct classes, otherwise fread interprets everything as factors
                         #select = 1:3541, #total number of columns in the data set
                         select = 802:1201,
                         #data.table = FALSE)[c(TRUE, rep(c(TRUE, FALSE, FALSE, FALSE),885))]  #To make subset of betas (First column is cpg names. Subsequent Column order was originally Beta, Gene Symbol, Chromosome, Gene Coordinate for each sample, so I made last three FALSE.)
                         data.table = FALSE)[c(rep(c(TRUE, FALSE, FALSE, FALSE),100))]
colnames(tcga.meth.f.300) <- fread("BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt", 
                                   header = FALSE, 
                                   skip = 0, 
                                   nrows = 1, #The first row is the true sample names. Above, I had to remove the row since the second row was also character but I wanted the columns to be interpreted as numeric
                                   #  select = 1:3541, 
                                   select = 802:1201,
                                   #  data.table = FALSE)[c(TRUE, rep(c(TRUE, FALSE, FALSE, FALSE),885))]
                                   data.table = FALSE)[c(rep(c(TRUE, FALSE, FALSE, FALSE),100))]
save(tcga.meth.f.300, file = "tcga.meth.300.RData")
tcga.meth.f.400 <- fread("BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt", 
                         header = TRUE, 
                         skip = 1, #this was done to preserve the correct classes, otherwise fread interprets everything as factors
                         #select = 1:3541, #total number of columns in the data set
                         select = 1202:1601,
                         #data.table = FALSE)[c(TRUE, rep(c(TRUE, FALSE, FALSE, FALSE),885))]  #To make subset of betas (First column is cpg names. Subsequent Column order was originally Beta, Gene Symbol, Chromosome, Gene Coordinate for each sample, so I made last three FALSE.)
                         data.table = FALSE)[c(rep(c(TRUE, FALSE, FALSE, FALSE),100))]
colnames(tcga.meth.f.400) <- fread("BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt", 
                                   header = FALSE, 
                                   skip = 0, 
                                   nrows = 1, #The first row is the true sample names. Above, I had to remove the row since the second row was also character but I wanted the columns to be interpreted as numeric
                                   #  select = 1:3541, 
                                   select = 1202:1601,
                                   #  data.table = FALSE)[c(TRUE, rep(c(TRUE, FALSE, FALSE, FALSE),885))]
                                   data.table = FALSE)[c(rep(c(TRUE, FALSE, FALSE, FALSE),100))]
save(tcga.meth.f.400, file = "tcga.meth.400.RData")
tcga.meth.f.500 <- fread("BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt", 
                         header = TRUE, 
                         skip = 1, #this was done to preserve the correct classes, otherwise fread interprets everything as factors
                         #select = 1:3541, #total number of columns in the data set
                         select = 1602:2001,
                         #data.table = FALSE)[c(TRUE, rep(c(TRUE, FALSE, FALSE, FALSE),885))]  #To make subset of betas (First column is cpg names. Subsequent Column order was originally Beta, Gene Symbol, Chromosome, Gene Coordinate for each sample, so I made last three FALSE.)
                         data.table = FALSE)[c(rep(c(TRUE, FALSE, FALSE, FALSE),100))]
colnames(tcga.meth.f.500) <- fread("BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt", 
                                   header = FALSE, 
                                   skip = 0, 
                                   nrows = 1, #The first row is the true sample names. Above, I had to remove the row since the second row was also character but I wanted the columns to be interpreted as numeric
                                   #  select = 1:3541, 
                                   select = 1602:2001,
                                   #  data.table = FALSE)[c(TRUE, rep(c(TRUE, FALSE, FALSE, FALSE),885))]
                                   data.table = FALSE)[c(rep(c(TRUE, FALSE, FALSE, FALSE),100))]
save(tcga.meth.f.500, file = "tcga.meth.500.RData")
tcga.meth.f.600 <- fread("BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt", 
                         header = TRUE, 
                         skip = 1, #this was done to preserve the correct classes, otherwise fread interprets everything as factors
                         #select = 1:3541, #total number of columns in the data set
                         select = 2002:2401,
                         #data.table = FALSE)[c(TRUE, rep(c(TRUE, FALSE, FALSE, FALSE),885))]  #To make subset of betas (First column is cpg names. Subsequent Column order was originally Beta, Gene Symbol, Chromosome, Gene Coordinate for each sample, so I made last three FALSE.)
                         data.table = FALSE)[c(rep(c(TRUE, FALSE, FALSE, FALSE),100))]
colnames(tcga.meth.f.600) <- fread("BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt", 
                                   header = FALSE, 
                                   skip = 0, 
                                   nrows = 1, #The first row is the true sample names. Above, I had to remove the row since the second row was also character but I wanted the columns to be interpreted as numeric
                                   #  select = 1:3541, 
                                   select = 2002:2401,
                                   #  data.table = FALSE)[c(TRUE, rep(c(TRUE, FALSE, FALSE, FALSE),885))]
                                   data.table = FALSE)[c(rep(c(TRUE, FALSE, FALSE, FALSE),100))]
save(tcga.meth.f.600, file = "tcga.meth.600.RData")
tcga.meth.f.700 <- fread("BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt", 
                         header = TRUE, 
                         skip = 1, #this was done to preserve the correct classes, otherwise fread interprets everything as factors
                         #select = 1:3541, #total number of columns in the data set
                         select = 2402:2801,
                         #data.table = FALSE)[c(TRUE, rep(c(TRUE, FALSE, FALSE, FALSE),885))]  #To make subset of betas (First column is cpg names. Subsequent Column order was originally Beta, Gene Symbol, Chromosome, Gene Coordinate for each sample, so I made last three FALSE.)
                         data.table = FALSE)[c(rep(c(TRUE, FALSE, FALSE, FALSE),100))]
colnames(tcga.meth.f.700) <- fread("BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt", 
                                   header = FALSE, 
                                   skip = 0, 
                                   nrows = 1, #The first row is the true sample names. Above, I had to remove the row since the second row was also character but I wanted the columns to be interpreted as numeric
                                   #  select = 1:3541, 
                                   select = 2402:2801,
                                   #  data.table = FALSE)[c(TRUE, rep(c(TRUE, FALSE, FALSE, FALSE),885))]
                                   data.table = FALSE)[c(rep(c(TRUE, FALSE, FALSE, FALSE),100))]
save(tcga.meth.f.700, file = "tcga.meth.700.RData")
tcga.meth.f.800 <- fread("BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt", 
                         header = TRUE, 
                         skip = 1, #this was done to preserve the correct classes, otherwise fread interprets everything as factors
                         #select = 1:3541, #total number of columns in the data set
                         select = 2802:3201,
                         #data.table = FALSE)[c(TRUE, rep(c(TRUE, FALSE, FALSE, FALSE),885))]  #To make subset of betas (First column is cpg names. Subsequent Column order was originally Beta, Gene Symbol, Chromosome, Gene Coordinate for each sample, so I made last three FALSE.)
                         data.table = FALSE)[c(rep(c(TRUE, FALSE, FALSE, FALSE),100))]
colnames(tcga.meth.f.800) <- fread("BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt", 
                                   header = FALSE, 
                                   skip = 0, 
                                   nrows = 1, #The first row is the true sample names. Above, I had to remove the row since the second row was also character but I wanted the columns to be interpreted as numeric
                                   #  select = 1:3541, 
                                   select = 2802:3201,
                                   #  data.table = FALSE)[c(TRUE, rep(c(TRUE, FALSE, FALSE, FALSE),885))]
                                   data.table = FALSE)[c(rep(c(TRUE, FALSE, FALSE, FALSE),100))]
save(tcga.meth.f.800, file = "tcga.meth.800.RData")
tcga.meth.f.85 <- fread("BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt", 
                         header = TRUE, 
                         skip = 1, #this was done to preserve the correct classes, otherwise fread interprets everything as factors
                         #select = 1:3541, #total number of columns in the data set
                         select = 3202:3541,
                         #data.table = FALSE)[c(TRUE, rep(c(TRUE, FALSE, FALSE, FALSE),885))]  #To make subset of betas (First column is cpg names. Subsequent Column order was originally Beta, Gene Symbol, Chromosome, Gene Coordinate for each sample, so I made last three FALSE.)
                         data.table = FALSE)[c(rep(c(TRUE, FALSE, FALSE, FALSE),85))]
colnames(tcga.meth.f.85) <- fread("BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt", 
                                   header = FALSE, 
                                   skip = 0, 
                                   nrows = 1, #The first row is the true sample names. Above, I had to remove the row since the second row was also character but I wanted the columns to be interpreted as numeric
                                   #  select = 1:3541, 
                                   select = 3202:3541,
                                   #  data.table = FALSE)[c(TRUE, rep(c(TRUE, FALSE, FALSE, FALSE),885))]
                                   data.table = FALSE)[c(rep(c(TRUE, FALSE, FALSE, FALSE),85))]
save(tcga.meth.f.85, file = "tcga.meth.85.RData")

#Quality Check/Control on new Data using CHAMP?
#Normalize?
#Run DMRCate 



getwd()
load("annotation_450K_v1.2.RData")
save(tcgaDMRcate, file = "tcgaDMRcate.RData")

###MAKE X 23 and Y 24 as convention
annotation.450K$CHR <- gsub("X", "23", annotation.450K$CHR)
annotation.450K$CHR <- gsub("Y", "23", annotation.450K$CHR)

#annotation file test
cpg.anot <- cbind(as.integer(annotation.450K$CHR), as.integer(annotation.450K$MAPINFO))
rownames(cpg.anot) <- c(annotation.450K$IlmnID)
colnames(cpg.anot) <- c("CHR", "MAPINFO")
cpg.anot <- as.data.frame(cpg.anot)
cpg.anot$IlmnID <- c(annotation.450K$IlmnID)
cpg.anot$address <- c(annotation.450K$AddressA_ID)
sum(tcgaDMRcate$DMRcateDMR$no.cpgs) #66554 cpgs
cpg.anot$genename <- c(annotation.450K$UCSC_RefGene_Name)




#MAKING dataframe for DMRs and removing "chr" from chromosome string so it will match with annotation
dmrs <- data.frame(tcgaDMRcate$DMRcateDMR$seqnames, tcgaDMRcate$DMRcateDMR$start, tcgaDMRcate$DMRcateDMR$end)
dmrs$tcgaDMRcate.DMRcateDMR.seqnames <- gsub("chr", "", dmrs$tcgaDMRcate.DMRcateDMR.seqnames)
dmrs$tcgaDMRcate.DMRcateDMR.seqnames <- as.integer(dmrs$tcgaDMRcate.DMRcateDMR.seqnames)

#Loop to identify cpgs in each DMR
begin.time <- Sys.time()
out <- NULL
for (i in 1:6030){
    chr <- dmrs[i, 1]
    start <- dmrs[i, 2]
    end <- dmrs[i, 3]
    outi <- subset(cpg.anot, CHR == chr & MAPINFO > start & MAPINFO < end )
    out <- rbind(out, outi)
}
end.time <- Sys.time()   
end.time - begin.time #2.26 minutes

tcga.meth.genename$Chromosome <- as.integer(tcga.meth.genename$Chromosome)

#Full set, before processing
begin.time <- Sys.time()
out.full <- NULL
for (i in 1:6030){
  chr <- dmrs[i, 1]
  start <- dmrs[i, 2]
  end <- dmrs[i, 3]
  out.fulli <- subset(tcga.meth.genename, Chromosome == chr & Genomic_Coordinate > start & Genomic_Coordinate < end )
  out.full <- rbind(out.full, out.fulli)
}
end.time <- Sys.time()   
end.time - begin.time #2.26 minutes
out.full <- out.full[,-2]
colnames(out.full)[1] <- "cpg"

   
setwd("C:/Users/hannumm/Documents/DMR Analysis")
save(out, file = "dmrcpgs.RData")
save(out.full, file = "dmrcpg_with_gene.RData")
save(tcga.meth.genename, file = "tcgamethgene.RData") #This file is the full 480K probes with gene names and coordinates
load("dmrcpg_with_gene.RData") #This file is the selected 69K probes with gene names and coordinates
load("tcgamethgene.RData")
out.full <- out.full[complete.cases(out.full[,2]),] #60677 that have gene symbol

length(unique(out.full$Gene_Symbol))

#Index where there are no multiples
length(grep(pattern = ";", x = out.full$Gene_Symbol, value = FALSE, invert = TRUE))

#index where there are multiples
#grep(pattern = ";", x = out.full$Gene_Symbol, value = FALSE)


#Splitting the gene symbols with more than one 
out.split <- out.full[grep(pattern = ";", x = out.full$Gene_Symbol, value = FALSE),] #6421
split <- strsplit(out.split$Gene_Symbol, ";")
#Repeat all elements of other rows by number of splits, make new df #15698
out.split.1 <- data.frame(cpg = rep(out.split$cpg, sapply(split, length)), Gene_Symbol = unlist(split), Chromosome = rep(out.split$Chromosome, sapply(split, length)), Genomic_Coordinate = rep(out.split$Genomic_Coordinate, sapply(split, length)))

#Select df without multiple gene symbols and remove ones with NA gene symbol #63070
out.split.sing <- out.full[grep(pattern = ";", x = out.full$Gene_Symbol, value = FALSE, invert = TRUE),] 
out.split.sing <- out.split.sing[complete.cases(out.split.sing[,2]),] #54256

#Combine cpgs with multiple genes and single genes
out.split.final <- rbind(out.split.1, out.split.sing)
setwd("C:/Users/hannumm/Documents/DMR Analysis")
save(out.split.final, file = "cpg.select.split.genename.RData") #This file is the selected cpgs with split gene names and genomic coordinates
load("cpg.select.split.genename.RData")

#Load in all 885 TCGA samples (100 at a time), then combine the subset that is in out.full
setwd("C:/Users/hannumm/Documents/DMR Analysis")
load("tcgamethgene.RData")
load("cpg.select.split.genename.RData")
load("tcga.meth.800.RData")
rownames(tcga.meth.f.800) <- tcga.meth.genename[,1]
#Test subset for tcga data #60677
tcga.meth.f.800 <- subset(tcga.meth.f.800, rownames(tcga.meth.f.800) %in% out.full$cpg)

load("tcga.meth.700.RData")
rownames(tcga.meth.f.700) <- tcga.meth.genename[,1]
#Test subset for tcga data #60677
tcga.meth.f.700 <- subset(tcga.meth.f.700, rownames(tcga.meth.f.700) %in% out.full$cpg)

load("tcga.meth.600.RData")
rownames(tcga.meth.f.600) <- tcga.meth.genename[,1]
#Test subset for tcga data #60677
tcga.meth.f.600 <- subset(tcga.meth.f.600, rownames(tcga.meth.f.600) %in% out.full$cpg)

load("tcga.meth.500.RData")
rownames(tcga.meth.f.500) <- tcga.meth.genename[,1]
#Test subset for tcga data #60677
tcga.meth.f.500 <- subset(tcga.meth.f.500, rownames(tcga.meth.f.500) %in% out.full$cpg)

load("tcga.meth.400.RData")
rownames(tcga.meth.f.400) <- tcga.meth.genename[,1]
#Test subset for tcga data #60677
tcga.meth.f.400 <- subset(tcga.meth.f.400, rownames(tcga.meth.f.400) %in% out.full$cpg)

load("tcga.meth.300.RData")
rownames(tcga.meth.f.300) <- tcga.meth.genename[,1]
#Test subset for tcga data #60677
tcga.meth.f.300 <- subset(tcga.meth.f.300, rownames(tcga.meth.f.300) %in% out.full$cpg)

load("tcga.meth.200.RData")

load("tcga.meth.100.RData")
rownames(tcga.meth.f.all) <- tcga.meth.f.all[,1]
tcga.meth.f.100 <- tcga.meth.f.all[,-1]
#Test subset for tcga data #60677
tcga.meth.f.100 <- subset(tcga.meth.f.all, rownames(tcga.meth.f.all) %in% out.full$cpg)
tcga.meth.f.100 <- tcga.meth.f.100[,-1]

load("tcga.meth.85.RData")
rownames(tcga.meth.f.85) <- tcga.meth.genename[,1]
#Test subset for tcga data #60677
tcga.meth.f.85 <- subset(tcga.meth.f.85, rownames(tcga.meth.f.85) %in% out.full$cpg)

tcga.meth.all.480 <- cbind(tcga.meth.f.100, tcga.meth.f.200, tcga.meth.f.300, tcga.meth.f.400, tcga.meth.f.500, tcga.meth.f.600, tcga.meth.f.700, tcga.meth.f.800, tcga.meth.f.85)
rownames(tcga.meth.all.480) <- tcga.meth.genename[,1]

tcga.meth.all <- cbind(tcga.meth.f.100, tcga.meth.f.200, tcga.meth.f.300, tcga.meth.f.400, tcga.meth.f.500, tcga.meth.f.600, tcga.meth.f.700, tcga.meth.f.800, tcga.meth.f.85)
setwd("C:/Users/hannumm/Documents/DMR Analysis")
save(tcga.meth.all, file = "tcga.meth.all.select.RData") #This file is the 885 samples with 60677 cpgs identified in DMRs
save(tcga.meth.all.480, file = "tcga.meth.all.480.RData") #This file is the 885 samples with 60677 cpgs identified in DMRs
load("tcga.meth.all.select.RData")
colnames(tcga.meth.all) <- substr(colnames(tcga.meth.all),1,12)

#tcga.meth.select.800 <- tcga.meth.select.800[complete.cases(tcga.meth.select.800),] #56747

####
# QC Measures on Methylation data
###
tcga.meth.qc <- as.matrix(tcga.meth.all)
tcga.meth.qc <- champ.filter(beta = tcga.meth.all, pd = NULL, filterDetP = FALSE, filterBeads = FALSE, fixOutlier = FALSE)
CpG.GUI(CpG = rownames(tcga.meth.qc), arraytype = "450K")

head(tcga.meth.qc1$beta[,1:10])

length(colMeans(is.na(tcga.meth.qc1$beta)))

tcga.meth.qc2 <- tcga.meth.qc$beta #55035 cpg 885 samples
#Remove probes with CpG-level coverage <= 0.7 (removed 2941 probes)
tcga.meth.qc2 <- tcga.meth.qc2[-which(rowMeans(is.na(tcga.meth.qc2)) > 0.3),] #52094 cpg 885 sample
#Remove samples with sample-level coverage <= 0.95
tcga.meth.qc4 <- t(tcga.meth.qc3) #transpose to be able to do rowwise
tcga.meth.qc4 <-  tcga.meth.qc4[-which(is.na(tcga.meth.qc4) > 0.05),] #52094 cpg 884 sample (1 sample removed)
tcga.meth.qc4 <- t(tcga.meth.qc4)
#impute missing values
tcga.meth.qc3 <- champ.impute(beta = tcga.meth.qc2, pd = pdtype$Design)
#Now use champ to normalize using the BMIQ method
tcga.meth.norm <- champ.norm(beta = tcga.meth.qc2, method = "BMIQ", cores = 1, plotBMIQ = TRUE,
                             resultsDir = "C:/Users/hannumm/Documents/DMR Analysis/Graphs/Champ.QC")
tcga.meth.norm <- BMIQ(tcga.meth.qc4)
pdtype <- as.data.frame(probeInfoALL.lv)
pdtype <-  subset(pdtype, pdtype$probeID %in% rownames(tcga.meth.qc4))

tcga.meth.f.800 <- as.matrix(t(tcga.meth.f.800))
#tcga.meth.qc <- t(tcga.meth.qc)
phen <- as.factor(c(rep("Tumor", 885)))
names(phen) <- rownames(tcga.meth.qc)
champ.QC(beta = tcga.meth.qc1, pheno = phen, resultsDir = "C:/Users/hannumm/Documents/DMR Analysis/Graphs/Champ.QC")
champ.QC(beta = tcga.meth.f.800, pheno = phen, resultsDir = "C:/Users/hannumm/Documents/DMR Analysis/Graphs/Champ.QC")

tcga.meth.f.100 <- as.matrix(tcga.meth.f.100)
tcga.meth.f.1001 <- champ.filter(beta = tcga.meth.f.100, pd = NULL, 
                                 filterDetP = FALSE, filterBeads = FALSE, fixOutlier = FALSE)
tcga.meth.norm <- champ.norm(beta = tcga.meth.f.1001$beta, method = "BMIQ", cores = 1, plotBMIQ = TRUE,
                             resultsDir = "C:/Users/hannumm/Documents/DMR Analysis/Graphs/Champ.QC")

probetype <- subset(probe.features, 
                  rownames(probe.features) %in% rownames(tcga.meth.qc2))
probetype$Type <- as.numeric(probetype$Type)
probetype <- probetype[ order(row.names(probetype)),]
tcga.meth.qc2 <- tcga.meth.qc2[order(row.names(tcga.meth.qc2)),]
#####
# Attempting to use wateRmelon package #
#####
norm.BetaVal=sapply(1:ncol(tcga.meth.qc2), function(i){
  print(i)
  BMIQ(tcga.meth.qc2[,i], probetype$Type)$nbeta
})
row.names(norm.BetaVal)=row.names(tcga.meth.qc2)
colnames(norm.BetaVal)=colnames(tcga.meth.qc2)

tcgabeta = NULL
for (i in 1:ncol(tcga.meth.qc2)){
  tcgabeta.i <- BMIQ(tcga.meth.qc2[,i], probetype$designType)$nbeta
}
  
probe.features$Type <- as.numeric(probe.features$Type)
test <- BMIQ(tcga.meth.qc2[,10], probetype$designType)
test <- BMIQ(tcga.meth.f.100[,10], probe.features$Type)

norm.BetaVal=sapply(1:ncol(tcga.meth.f.100), function(i){
  print(i)
  BMIQ(tcga.meth.f.100[,i], probe.features$Type)$nbeta
})
row.names(norm.BetaVal)=row.names(tcga.meth.f.100)
colnames(norm.BetaVal)=colnames(tcga.meth.f.100)


#######################################################
# NITTY GRITTY SUBSET METHYLATION AND EXPRESSION TIME #
#######################################################


#Load TCGA gene expression data
setwd("H:/Biostatistics/Ronglai-Margaret/TCGABRCA")
load("H:/Biostatistics/Ronglai-Margaret/TCGABRCA/BRCA_all.Rdata")

#Gene expression data
tcga.mrna <- t(BRCA_all$mrna)
#Split names since gene names are in front of | character
rownames(tcga.mrna) <- sapply(strsplit(rownames(tcga.mrna), "|", fixed = TRUE), '[[', 1)
#Subset gene expression and meth data with each other to match gene names
tcga.mrna.select <- subset(tcga.mrna, rownames(tcga.mrna) %in% out.full$Gene_Symbol) #3816 genes that are also in methylation data
out.full.select <- subset(out.full, out.full$Gene_Symbol %in% rownames(tcga.mrna.select)) #45K probes
tcga.meth.select <- subset(tcga.meth.all, rownames(tcga.meth.all) %in% out.full.select$cpg)
#Add columns for cpg and gene symbol (Decided I didn't need to do this)
#tcga.meth.select$cpg <- out.full.select$cpg
#tcga.meth.select$Gene_Symbol <- out.full.select$Gene_Symbol
#tcga.meth.select$id <- do.call(paste, c(tcga.meth.select[886:887], sep = "."))
#tcga.meth.select <- tcga.meth.select[c(888,887,886, 1:885)]

#Subset methylation and expression based on sample names
#Have to transpose first since subset works on rows
tcga.mrna.select <- t(tcga.mrna.select)
tcga.meth.select <- t(tcga.meth.select)
tcga.mrna.select <- subset(tcga.mrna.select, rownames(tcga.mrna.select) %in% rownames(tcga.meth.select))
tcga.meth.select <- subset(tcga.meth.select, rownames(tcga.meth.select) %in% rownames(tcga.mrna.select))
###Matching sample names only left 663 samples??

#Transpose back
tcga.mrna.select <- t(tcga.mrna.select)
tcga.meth.select <- t(tcga.meth.select)

tcga.mrna.select <- tcga.mrna.select[,order(colnames(tcga.mrna.select))]
tcga.meth.select <- tcga.meth.select[,order(colnames(tcga.meth.select))]

######
#Analysis with the correct split ##
###################################
tcga.mrna.select1 <- subset(tcga.mrna, rownames(tcga.mrna) %in% out.split.final$Gene_Symbol) #4409 genes
out.split.select1 <- subset(out.split.final, out.split.final$Gene_Symbol %in% rownames(tcga.mrna.select1)) #58K probes
tcga.meth.select1 <- subset(tcga.meth.all, rownames(tcga.meth.all) %in% out.split.select1$cpg)
#Subsetting based on gene symbols: 4409 genes in mrna select and 51729 cpgs in meth select
out.split.select1 <- subset(out.split.select1, out.split.select1$cpg %in% rownames(tcga.meth.select1))
tcga.mrna.select1 <- t(tcga.mrna.select1)
tcga.meth.select1 <- t(tcga.meth.select1)
tcga.mrna.select1 <- subset(tcga.mrna.select1, rownames(tcga.mrna.select1) %in% rownames(tcga.meth.select1))
tcga.meth.select1 <- subset(tcga.meth.select1, rownames(tcga.meth.select1) %in% rownames(tcga.mrna.select1))
###Matching sample names only left 663 samples??

#Transpose back
tcga.mrna.select1 <- t(tcga.mrna.select1)
tcga.meth.select1 <- t(tcga.meth.select1)

length(unique(out.split.select1$cpg)) #51729 unique genes

#Splitting the gene symbols in meth.select with more than one 
out.full1 <- subset(out.full, out.full$cpg %in% out.split.select1$cpg)
tcga.meth.select2 <- as.data.frame(tcga.meth.select1)
tcga.meth.select2$Gene_Symbol <- c(out.full1$Gene_Symbol)
tcga.meth.select2$cpg <- c(out.full1$cpg)

meth.split <- tcga.meth.select2

meth.split <- tcga.meth.select2[grep(pattern = ";", x = tcga.meth.select2$Gene_Symbol, value = FALSE),] #5862 cpgs with multiple genes
#split <- strsplit(meth.split$Gene_Symbol, ";")

#Repeat all elements of other rows by number of splits
meth.split <- meth.split %>%
  mutate(Gene_Symbol = strsplit(as.character(Gene_Symbol), ";")) %>%
  unnest(Gene_Symbol) #14420

#Repeat all elements of other rows by number of splits, make new df #15698

#Select df without multiple gene symbols and remove ones with NA gene symbol #63070
meth.split.sing <- tcga.meth.select2[grep(pattern = ";", x = tcga.meth.select2$Gene_Symbol, value = FALSE, invert = TRUE),] #45867
meth.split.sing <- meth.split.sing[complete.cases(meth.split.sing[,2]),] #54256

#Combine cpgs with multiple genes and single genes
meth.split.final <- rbind(meth.split, meth.split.sing) #60287
#rownames(meth.split.final) <- c(meth.split.final$cpg)
meth.split.final1 <- meth.split.final[,-665]
out.split.select1 <- subset(out.split.select1, out.split.select1$cpg %in% meth.split.final1$cpg)
meth.split.final1 <- subset(meth.split.final1, meth.split.final1$cpg %in% out.split.select1$cpg)


head(colnames(tcga.meth.select))

id <- NULL
rho <- NULL
for (i in 1:length(unique(out.full.select$Gene_Symbol))){
 #for (i in 1:100){
    
  gene <- unique(out.full.select$Gene_Symbol)[i]
  subanno <- subset(out.full.select, out.full.select$Gene_Symbol == gene)
  #CG PROBE ID FROM SUBSET use it to extract rows/cols methylation data
  meth.mat.i <- tcga.meth.select[subanno$cpg, , drop = F]
  exp.mat.i <- tcga.mrna.select[gene, , drop=F]
  rho.i <- cor(t(exp.mat.i), t(meth.mat.i), use="pairwise.complete.obs")
  rho <- c(rho, rho.i)
  id.i <- paste(colnames(rho.i), rownames(rho.i), sep=".")
  #id.i <- do.call(paste, c(subanno[1:2], sep = "."))
  id <- c(id, id.i)
  }

names(rho) <- id
setwd("C:/Users/hannumm/Documents/DMR Analysis")
save(rho, file = "rho.RData") #This file is the correlations for all 45867 cpgs with gene expression
load("rho.RData")

rho.df <- as.data.frame(rho)

hist(rho.df, main = "Histogram of Correlation")
summary(rho.df) #Min -0.4, 1Q -0.06, Med -0.0264, Mean -0.0265, 3Q 0.0139, Max 0.2146
#NA's: 2519

rho.df.3 <- subset(rho.df, rho.df$rho < -0.3) #Only 19 cpgs??

unique <- unique(out.full1$Gene_Symbol)
unique$Chromosome <- as.factor()


#Trying to make loop with summary information, using proper split data
id <- NULL
rho <- NULL
sd <- NULL
#for (i in 1:length(unique(out.split.select1$Gene_Symbol))){
  for (i in 1:39){
  
  gene <- unique(out.split.select1$Gene_Symbol)[i]
  subanno <- subset(out.split.select1, out.split.select1$Gene_Symbol == gene)
  #CG PROBE ID FROM SUBSET use it to extract rows/cols methylation data
  meth.mat.i <- meth.split.final1[subanno$cpg, , drop = F]
  exp.mat.i <- tcga.mrna.select1[gene, , drop=F]
  rho.i <- cor(t(exp.mat.i), t(meth.mat.i), use="pairwise.complete.obs")
  sd.i <- sd(rho.i)
  rho <- c(rho, rho.i)
  sd <- c(sd, sd.i)
  id.i <- paste(colnames(rho.i), rownames(rho.i), sep=".")
  #id.i <- do.call(paste, c(subanno[1:2], sep = "."))
  id <- c(id, id.i)
}
names(sd) <- unique(out.full.select$Gene_Symbol)
summary(sd)
#Rank vector from largest to smallest

####MAKE SUMMARY OF GENES and probes per gene/average





#cpg.test <- subset(cpg.anot, cpg.anot$CHR == 6)

#cpg.test[cpg.test$CHR == tcgaDMRcate$DMRcateDMR$seqnames | cpg.test$MAPINFO >= 33130696 | cpg.test$MAPINFO <= 33149777]
#cpg.test[which(cpg.test$MAPINFO >= 33130696 | cpg.test$MAPINFO <= 33149777)]
#cpg.test.1 <- subset(cpg.test, cpg.test$MAPINFO >= 33130696 & cpg.test$MAPINFO <= 33149777)
#cpg.test.2 <- subset(cpg.test, cpg.test$MAPINFO %in% 33130696:33149777)
#subset(probes, probes$XY.probes %in% cpg.test.1$IlmnID) #Even though above subsets produce 212 probes, this subset does not produce anything, indicating XYProbees output does not reflect what I thought it did?
#subset(probes.1, probes.1$crosshyb %in% cpg.test.1$IlmnID) #this one came up with 4 probes!

#cpg.all <- subset(cpg.anot, cpg.anot$MAPINFO %in% tcgaDMRcate$DMRcateDMR$start:tcgaDMRcate$DMRcateDMR$end)
#cpg.all <- filter(cpg.anot, filter(cpg.anot$MAPINFO >= tcgaDMRcate$DMRcateDMR$start & cpg.anot$MAPINFO <= tcgaDMRcate$DMRcateDMR$end))
#head(tcgaDMRcate$DMRcateDMR$start)


#####
#MATCHING CpGs in TCGA FULL DATA TO DMRS ##
#####

#probes <- as.data.frame(XY.probes) #Using XY Probes output from DMRcate
#probes.1 <- as.data.frame(crosshyb)

#tcga.meth.f.100.20K <- subset(tcga.meth.f.all, tcga.meth.f.all$`Hybridization REF` %in% probes$XY.probes)
#tcga.meth.f.100.20K1 <- subset(tcga.meth.f.all, tcga.meth.f.all$`Hybridization REF` %in% probes.1$crosshyb)

