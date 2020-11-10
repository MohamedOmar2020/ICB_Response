
############################################################################
### Clean working space ans set working directory
############################################################################

### Clean
rm(list=ls())

### Set working directory
setwd("/Volumes/Macintosh/Research/Projects/ICB_Response")

### Set global options
options(width=130)


############################################################################
### Load libraries
library(Biobase)
library(MergeMaid)
library(parallel)
library(mclust)
library(limma)
library(dplyr)
library(edgeR)
library(genefilter)
###########################################################################
### Useful functions
###########################################################################

### Functions
source("/Volumes/Macintosh/Research/Tutorials/plotEMdistr.R")

### Panel correlations
panel.cor <- function(x, y, digits = 3, prefix = "Correlation:\n", cex.cor=1)  {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  text(0.5, 0.5, txt, cex = cex.cor)
}

### scatterPlot  with regression line
panel.scatterFit <- function(x, y) {
  points(x, y, pch=16)
  abline(lsfit(x,y), col="blue", lty=2)
}

### Geometric mean
geoMean <- function(x, ...) exp(mean(log(x), na.rm=TRUE, ...))

### Nice histogram
niceHist <- function(x, probability = TRUE, nclass=25,
                     xlab="Overall Integrative Correlation", ylab="Density", ...) {
  hist(x, probability = probability, col="salmon",
       ylab=ylab, nclass=nclass, xlab=xlab, ...)
  lines(density(x), lwd=2, col="blue")
  abline(v=median(x), lwd=2, lty=2, col="blue")
}


###########################################################################
### Read/Load data
###########################################################################
load("./Objs/exprsAll.rda")

##################

### Quantile normalize
allExpr <- lapply(exprsAll, normalizeBetweenArrays)

### Check dimensions
sapply(exprsAll, dim)

names(exprsAll) <- c("NP", "HUGO", "TCGA", "RIAZ_Pre_aPD1", "RIAZ_On_aPD1", "MGH_On_aPD1", "MGH_On_aCTLA4", "VanAllen")
###########################################################################
### Prepapre object for MergeMaid analysis
###########################################################################

###########################################################################
### Uncollapsed data
### Get matrices out of the list (necessary to work with MergeMaid)
expr1 <- as.matrix(exprsAll$NP)
expr2 <- as.matrix(exprsAll$HUGO)
expr3 <- as.matrix(exprsAll$TCGA)
expr4 <- as.matrix(exprsAll$RIAZ_Pre_aPD1)
expr5 <- as.matrix(exprsAll$RIAZ_On_aPD1)
expr6 <- as.matrix(exprsAll$MGH_On_aPD1)
expr7 <- as.matrix(exprsAll$MGH_On_aCTLA4)
expr8 <- as.matrix(exprsAll$VanAllen)


  
### Merge data and clean object names
mgdExpr <- mergeExprs(expr1, expr2, expr3, expr4, expr5, expr6, expr7, expr8)



###########################################################################
### Now calculate the correlations
###########################################################################

### Calculate integrative correlation: original values
iCorP <- intCor(mgdExpr, method="pearson", exact = FALSE)

###########################################################################
### Explore reproducibility using plots
###########################################################################

###########################################################################
### ALL DATA

### Open device
png("./Figs/Correlation/iCorDensPearson.png", width=2000, height=2000, res=200)
### Set graphic layout
par(oma=rep(0.5,4), mar=c(2.75,2.75,0.75,0.25), mgp=c(1.75, 0.75, 0))
### Compute and plot distributions
set.seed(333)
iCorNullP <- intcorDens(mgdExpr, cex.legend=1.1, lwd=2)
### Close device
dev.off() 

### Open device
png("./Figs/Correlation/iCorHistPearson.png", width=2000, height=2000, res=200)
### Set graphic layout
par(oma=rep(0.25,4), mar=c(2.75,2.75,0.75,0.25), mgp=c(1.75, 0.75, 0))
### Compute and plot distributions
### hist(iCorP, col="salmon")
niceHist(integrative.cors(iCorP), main="Integrative Correlation")
### Close device
dev.off()


###########################################################################

###########################################################################
### EM-algorithm to split genes based on iCOR
###########################################################################
set.seed(333)

null <- unlist(iCorNullP)
iCorNullP_Modified <- na.omit(sample(null, length(integrative.cors(iCorP)), replace = TRUE))
### Create list of integrative correlations omit NAs
myDat <- data.frame("Corr_All"=na.omit(integrative.cors(iCorP)),
                    "Null_Distribution_All"= iCorNullP_Modified)


### ### Use Mclust with 3 normals using mclapply from the "parallel" package
datEM <- mclapply(myDat, Mclust, mc.cores=detectCores(), G=1:2)


###########################################################################
### Plot
png("./Figs/Correlation/iCorPearsonEMclass.png", height=3500, width=3500, res=350)

### Set layout
nc <- floor(sqrt(ncol(myDat)))
nr <- ceiling((ncol(myDat))/nc)
layout(matrix(1:(nc*nr), ncol=nc, nrow= nr, byrow=FALSE))
### Set par
par(list(mar=c(3, 3, 4, 1), mgp=c(1.5, 0.5, 0)))
### Plot the distributions
splitPoint <- mapply(exprDat=myDat, emOut=datEM, name=gsub("_mRNA", "", colnames(myDat)),
                     FUN=plotEMdistr, MoreArgs=list(quant = c(0.1, 0.9)),
                     forXlab="Integrative Correlation", forTitle="Integrative Correlation")
### Close screen and device
dev.off()


##########################################################################
###########################################################################
### Reproducible genes: Pearson
iCorPint <- integrative.cors(iCorP)
### More then selected quantile
thr <- 0.1
#thr <- quantile(iCorPint, 0.66)

table(iCorPint > thr)
keep <- which(iCorPint > thr)
iCorPint_T <- iCorPint[keep]

RGenes <- names(iCorPint_T)
save(iCorPint, RGenes, file = "./Objs/RGenes.rda")


###########################################################################
###########################################################################
### Find split point for each pair-wise correlation in comparison to null

myDatWithNull <- data.frame(
  ## Join real and null data to run EM-classification for ALL
  Corr_All = c(myDat$Corr_All, myDat$Null_Distribution_All))


### ### Use Mclust with 3 normals using mclapply from the "parallel" package
datEM <- mclapply(myDatWithNull, Mclust, mc.cores=detectCores(), G=1:2)


###########################################################################
### Plot
png("figs/iCorPearsonEMclassNULLvsREAL.png", height=2500, width=3500, res=300)
### Set layout
nc <- floor(sqrt(ncol(myDat)))
nr <- ceiling((ncol(myDat))/nc)
layout(matrix(1:(nc*nr), ncol=nc, nrow= nr, byrow=TRUE))
### Set par
par(list(mar=c(3, 3, 2.5, 1), mgp=c(1.5, 0.5, 0)))
### Plot the distributions
splitPointVSnull <- mapply(exprDat=myDatWithNull, emOut=datEM, name=gsub("_mRNA", "", colnames(myDat)),
                           FUN=plotEMdistr, MoreArgs=list(quant = c(0.1, 0.9),  upper=2.1),
                           forXlab="Integrative Correlation", forTitle="Integrative Correlation")
### Close screen and device
dev.off()


###########################################################################
### Save
save(eDat, eDatMed, eDatMean, eDatQuant, mgdExpr, mgdExprMed, mgdExprMean, mgdExprQuant,
     iCorP, iCorMedP, iCorMeanP,  iCorQuantP,  iCorNullP, iCorNullMedP, iCorNullMeanP, iCorNullQuantP,
     splitPoint, splitPointVSnull, file="objs/iCorData.rda")


###########################################################################
### Session information
sessionInfo()

### Clean working space and quit
rm(list=ls())
q("no")
