################################################################
## TCGA BRCA Integrated Clustering with iCluster              ##
## Author: Margaret Hannum                                    ##
## Orig Date: 6/26/2017                                       ##
## Updated: 7/17/2017                                         ##
################################################################



#source("https://bioconductor.org/biocLite.R")
#biocLite("iClusterPlus")


library(iClusterPlus)
library(iC10)
library(gplots)
library(lattice)
library(fossil) #for rand index function
library(mclust)

setwd("H:/Biostatistics/Ronglai-Margaret/TCGABRCA")
load("H:/Biostatistics/Ronglai-Margaret/TCGABRCA/BRCA_all.Rdata")
head(BRCA_all$cn[,1:5])

#Gene expression data
tcga.mrna <- BRCA_all$mrna
#Split names since gene names are in front of | character
colnames(tcga.mrna) <- sapply(strsplit(colnames(tcga.mrna), "|", fixed = TRUE), '[[', 1)

####
#Match gene expression features
####
#Load Selected GE Features from top 1000 probes (715 features)
setwd("C:/Users/hannumm/Documents/CNA Analysis")
ge.feat <- read.csv("C:/Users/hannumm/Documents/CNA Analysis/Selected GE Features.csv")
#Subset of expression data from selected features- to get 579 features
tcga.mrna.df <- as.data.frame(t(tcga.mrna))
tcga.mrna.df <-subset(tcga.mrna.df, rownames(tcga.mrna.df) %in% ge.feat$Gene_symbol)

#Getting CNA Segments in matrix form
tcga.cna.seg <- read.table("H:/Biostatistics/Ronglai-Margaret/TCGABRCA/all_data_by_genes.txt", fill = TRUE, row.names = 1, header = TRUE)
#tcga.cna.seg.m <- as.matrix(tcga.cna.seg)
tcga.cna.seg.m <- data.matrix(tcga.cna.seg[5:1084])
#tcga.cna.seg.m <- tcga.cna.seg.m[,-(1:4)]
colnames(tcga.cna.seg.m) <- substr(colnames(tcga.cna.seg.m),1,12)
#tcga.cna.seg.m2 <- as.matrix(tcga.cna.seg.m)
#Have to sub - for . so we can match colnames with mrna data
colnames(tcga.cna.seg.m) <- gsub("\\.", "-", colnames(tcga.cna.seg.m))

#Match CNA features with ones selected on gene expression
tcga.cna.m.select <-subset(tcga.cna.seg.m, rownames(tcga.cna.seg.m) %in% rownames(tcga.mrna.df))
tcga.mrna.m.select <-as.matrix(subset(tcga.mrna.df, rownames(tcga.mrna.df) %in% rownames(tcga.cna.m.select)))

###
#Now we need to match the samples, since cna has 1080 but mrna has 960
#transpose matrices since subset works on rows
tcga.cna.m.select <- t(tcga.cna.m.select)
tcga.mrna.m.select <- t(tcga.mrna.m.select)
tcga.cna.m.select <-subset(tcga.cna.m.select, rownames(tcga.cna.m.select) %in% rownames(tcga.mrna.m.select))
tcga.cna.m.select <- tcga.cna.m.select[,order(colnames(tcga.cna.m.select))]
tcga.mrna.m.select <- tcga.mrna.m.select[,order(colnames(tcga.mrna.m.select))]

#Combine matrices into list after preparation, as required by iCluster
#tcga.list <- list(tcga.cna.m, tcga.mrna.m)
tcga.list <- list(tcga.cna.m.select, tcga.mrna.m.select)

#####################
# iCLUSTER ANALYSIS #
#####################
begin.time = Sys.time()
lambda <- alist()
lambda[[1]]=0
lambda[[2]]=0
fit.cn.ge <- iCluster2(tcga.list, K=10, lambda = lambda, method = c("lasso", "lasso"))
end.time = Sys.time()
end.time-begin.time

#compute proportion of deviance
plotiCluster(fit = fit.cn.ge)
pod.10 <- compute.pod(fit.cn.ge)
pod.10/10

#plot Heatmap
plotHeatmap(fit.cn.ge, tcga.list, type = c("gaussian", "gaussian"))

fit.list <- list()
for(k in 9:11){
  fit.cn[[k]] <- iCluster2(tcga.list, K=k, lambda = lambda, method = c("lasso", "lasso"))
  plotiCluster(fit = fit.cn[[k]])
  pod[[k]] <- compute.pod(fit.cn.[[k]])
  plot(pod[[k]]/k)
}

for(k in 9:11){
  plotiCluster(fit = iCluster2(tcga.list, K=k, lambda = lambda, method = c("lasso", "lasso")))
  plot(compute.pod(iCluster2(tcga.list, K=k, lambda = lambda, method = c("lasso", "lasso")))/k)
}

fit.cn3 <- iCluster2(tcga.list, K=3, lambda = lambda, method = c("lasso", "lasso"))
plotiCluster(fit = fit.cn3)
pod3 <- compute.pod(fit.cn3)

fit.cn4 <- iCluster2(tcga.list, K=4, lambda = lambda, method = c("lasso", "lasso"))
plotiCluster(fit = fit.cn4)
pod4 <- compute.pod(fit.cn4)

fit.cn5 <- iCluster2(tcga.list, K=5, lambda = lambda, method = c("lasso", "lasso"))
plotiCluster(fit = fit.cn5)
pod5 <- compute.pod(fit.cn5)

fit.cn6 <- iCluster2(tcga.list, K=6, lambda = lambda, method = c("lasso", "lasso"))
plotiCluster(fit = fit.cn6)
pod6 <- compute.pod(fit.cn6)

fit.cn7 <- iCluster2(tcga.list, K=7, lambda = lambda, method = c("lasso", "lasso"))
plotiCluster(fit = fit.cn7)
pod7 <- compute.pod(fit.cn7)

fit.cn8 <- iCluster2(tcga.list, K=8, lambda = lambda, method = c("lasso", "lasso"))
plotiCluster(fit = fit.cn8)
pod8 <- compute.pod(fit.cn8)

fit.cn9 <- iCluster2(tcga.list, K=9, lambda = lambda, method = c("lasso", "lasso"))
plotiCluster(fit = fit.cn9)
pod9 <- compute.pod(fit.cn9)

fit.cn10 <- iCluster2(tcga.list, K=10, lambda = lambda, method = c("lasso", "lasso"))
plotiCluster(fit = fit.cn10)
pod10 <- compute.pod(fit.cn10)

fit.cn11 <- iCluster2(tcga.list, K=11, lambda = lambda, method = c("lasso", "lasso"))
plotiCluster(fit = fit.cn11)
pod11 <- compute.pod(fit.cn11)

fit.cn12 <- iCluster2(tcga.list, K=12, lambda = lambda, method = c("lasso", "lasso"))
plotiCluster(fit = fit.cn12)
pod12 <- compute.pod(fit.cn12)

fit.cn13 <- iCluster2(tcga.list, K=13, lambda = lambda, method = c("lasso", "lasso"))
plotiCluster(fit = fit.cn13)
pod13 <- compute.pod(fit.cn13)

fit.cn14 <- iCluster2(tcga.list, K=14, lambda = lambda, method = c("lasso", "lasso"))
plotiCluster(fit = fit.cn14)
pod14 <- compute.pod(fit.cn14)

fit.cn15 <- iCluster2(tcga.list, K=15, lambda = lambda, method = c("lasso", "lasso"))
plotiCluster(fit = fit.cn15)
pod15 <- compute.pod(fit.cn15)

fit.cn16 <- iCluster2(tcga.list, K=16, lambda = lambda, method = c("lasso", "lasso"))
plotiCluster(fit = fit.cn16)
pod16 <- compute.pod(fit.cn16)

fit.cn17 <- iCluster2(tcga.list, K=17, lambda = lambda, method = c("lasso", "lasso"))
plotiCluster(fit = fit.cn17)
pod17 <- compute.pod(fit.cn17)

fit.cn18 <- iCluster2(tcga.list, K=18, lambda = lambda, method = c("lasso", "lasso"))
plotiCluster(fit = fit.cn18)
pod18 <- compute.pod(fit.cn18)

fit.cn19 <- iCluster2(tcga.list, K=19, lambda = lambda, method = c("lasso", "lasso"))
plotiCluster(fit = fit.cn19)
pod19 <- compute.pod(fit.cn19)

fit.cn20 <- iCluster2(tcga.list, K=20, lambda = lambda, method = c("lasso", "lasso"))
plotiCluster(fit = fit.cn20)
pod20 <- compute.pod(fit.cn20)

pod.all <- c(pod3, pod4, pod5, pod6, pod7, pod8, pod9, pod10, pod11, pod12, 
             pod13, pod14, pod15, pod16, pod17, pod18, pod19, pod20)
clust <- c(3:20)

pod.tab <- as.data.frame(cbind(clust, pod.all))

plot(x= pod.tab$clust, y=pod.tab$pod.all, 
     main = "Proportion of Deviation for Clusters 3-20",
     xlab = "K Clusters",
     ylab = "Proportion of Deviation")

###################
# Attempting Rand Index calculation comparing k=10,11,12 with ic10 assignments for TCGA data
######################
head(res2$class)
head(fit.cn10$clusters)
intclust.assign <- as.matrix(c(res2$class))
icluster.10.assign <- as.matrix(c(fit.cn10$clusters))
rand.index(intclust.assign, icluster.10.assign) #0.802
icluster.11.assign <-as.matrix(c(fit.cn11$clusters))
rand.index(intclust.assign, icluster.11.assign) #0.809
icluster.12.assign <-as.matrix(c(fit.cn12$clusters))
rand.index(intclust.assign, icluster.12.assign) #0.818



#Alternatively could use ic10 matchFeatures function, chose to just do the above
#tcga.mrna.select <- matchFeatures(Exp = t(tcga.mrna), Exp.by.feat = "gene", ref = "hg19")
#tcga.mrna.select <- normalizeFeatures(tcga.mrna.select)
#res2 <- iC10(tcga.mrna.select)
#summary(res2)

#tcga.cna.seg <- as.matrix(tcga.cna.seg)
#colnames(tcga.cna.seg) <- tcga.cna.seg[1,]
#tcga.cna.seg <- tcga.cna.seg[-1,]
#tcga.cna.seg.df <- as.data.frame(tcga.cna.seg)
#Need to extract 1000 genes for supervised clustering
#tcga.cna.seg.df <-subset(tcga.cna.seg.df, tcga.cna.seg.df$Gene %in% ge.feat$Gene_symbol)
#544 genes for 1085 pts

#Matching copy number and gene expression features with each other
#Able to match 542 features
#tcga.mrna.df <-subset(tcga.mrna.df, rownames(tcga.mrna.df) %in% tcga.cna.seg.df$Gene)
#tcga.cna.seg.df <-subset(tcga.cna.seg, rownames(tcga.cna.seg) %in% rownames(tcga.mrna.df))

#tcga.mrna.m <- as.matrix(tcga.mrna.df)
#tcga.mrna.m <- tcga.mrna.m[,order(colnames(tcga.mrna.m))]
#tcga.cna.seg.m <- as.matrix(tcga.cna.seg.df)
#rownames(tcga.cna.seg.m) <- tcga.cna.seg.m[,1]
#tcga.cna.seg.m <- tcga.cna.seg.m[,-(1:5)]
#Need to alter colnames of cna.seg.m so they match mrna.m
#colnames(tcga.cna.seg.m) <- sapply(strsplit(as.character(colnames(tcga.cna.seg.m)), split = "12"), '[[', 1)
#colnames(tcga.cna.seg.m) <-gsub(".{16}$", "", as.character(colnames(tcga.cna.seg.m)), fixed = FALSE)
#colnames(tcga.cna.seg.m) <- substr(colnames(tcga.cna.seg.m),1,12)
#tcga.cna.df <- as.data.frame(t(tcga.cna.seg.m))
#tcga.mrna.df <- as.data.frame(t(tcga.mrna.m))

#Subset only seems to work on rows, not cols
#tcga.cna.df <-subset(tcga.cna.df, as.character(rownames(tcga.cna.df)) %in% as.character(rownames(tcga.mrna.df)))
#tcga.cna.m <- as.matrix(tcga.cna.df)
#tcga.cna.m <- tcga.cna.m[,order(colnames(tcga.cna.m))]
##Trying to subset the cna.seg 1080 samples into the 960 samples from mrna.m
#tcga.cna.seg.m <-subset(tcga.cna.seg.m, as.character(colnames(tcga.cna.seg.m)) %in% as.character(colnames(tcga.mrna.m)))
#tcga.mrna.m <-subset(tcga.mrna.m, as.character(colnames(tcga.mrna.m)) %in% as.character(colnames(tcga.cna.seg.m)))

##icluster having an issue becuase cna.m is listed as chr class for some reason
indx <- sapply(tcga.cna.m, is.factor)
tcga.cna.m1 <- lapply(tcga.cna.m[indx], function(x) as.numeric(as.character(x)))
tcga.cna.m1 <- matrix(data=tcga.cna.m1, ncol=542, nrow = 960)
colnames(tcga.cna.m1) <- colnames(tcga.cna.m)
rownames(tcga.cna.m1) <- rownames(tcga.cna.m)



setwd("C:/Users/hannumm/Documents/R.Images")
save.image(file = "tcga.iCluster.RData")
save.image(file = "tcga.iCluster.new.RData")
load("tcga.iCluster.RData")
save(tcga.list, file = "tcga.list.icluster.RData" )
load("tcga.list.icluster.RData")


###
#All the stuff from ic10, DMR id and iCluster prior to analysis!
setwd("C:/Users/hannumm/Documents/R.Images")
save.image(file = "tcga.DMR.ic10.RData")
