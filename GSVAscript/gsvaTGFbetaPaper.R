##################################################################################################################
# R script to run GSVA on Rody and BT set 
# By Radia Johnson, Ph.D
# Date: June 15th, 2017
# On UNIX LAPTOP
##################################################################################################################

library(GSEABase)
library(Biobase)
library(genefilter)
library(limma)
library(RColorBrewer)
library(GSVA)
library(gdata)
library(gplots)

# To convert entrez ids to gene symbols
library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)

source("convertGeneNames.R")
GeneConvTable <- readRDS("Gene.Conv.Table.rds")

getMissingSymbolIds <- function(Symbol.ID, Gene.Conv.Table=GeneConvTable){

    repIds <- match(names(Symbol.ID), as.character(Gene.Conv.Table$EntrezGeneID))
    inds <- which(!is.na(repIds))
    replaceGenes <- as.character(Gene.Conv.Table$GeneSymbol[repIds[inds]])
    names(replaceGenes) <- as.character(Gene.Conv.Table$EntrezGeneID[repIds[inds]])
    
    # look at the genes found
    print(head(replaceGenes))
    Symbol.ID[names(replaceGenes)] <- replaceGenes
    return(Symbol.ID)
    
}

# Set the working directory
setwd("~/Documents/TNBCimmune")

#-----------------------------------------------------------------
# GET THE MATRIX NEEDED TO RUN GSVA
#-----------------------------------------------------------------
load("ExpSetBT.RData") # loads edataBT

M1 <- edataBT

# Get the entrez ids for the gene names
geneSymbols.ID <- mapIds(org.Hs.eg.db, keys=rownames(M1), column="SYMBOL", keytype="ENTREZID", multiVals="first")

geneSymbols.ID1 <- getMissingSymbolIds(geneSymbols.ID)
dannot <- as.data.frame(geneSymbols.ID1)

M2 <- M1
rownames(M2) <- dannot[rownames(M1), 1]

#-----------------------------------------------------------------
# Load the clinical features for the dataset
#-----------------------------------------------------------------

clinical <- read.xls("Patient-ranking-Please-plot-PaperID.xlsx", row.names=1)
colnames(M2)[c(27, 33, 15)] <- c("86", "41", "85")
head(M2)

# Get the revised order for the patients to be plotted, stringsAsFactors=FALSE
paperOrder <- read.xls("/Users/radiajohnson/Documents/TNBCimmune/TGBF/TGFBeta-GSVA-order.xlsx", row.names=2, stringsAsFactors=FALSE)
colj <- paperOrder[, 9]
names(colj) <- rownames(paperOrder)

colk <- paperOrder[, 10]
names(colk) <- rownames(paperOrder)

signPaths<- read.xls("TGFBeta-pathway-significant-GSEA.xlsx", sheet=2, stringsAsFactors=FALSE)
signPATHS <- lapply(as.list(signPaths), function(x) x[x != ""])

#-----------------------------------------------------------------
# Load the TGFbeta genesets
#-----------------------------------------------------------------

signs2 <- read.xls("All-TGFB-up-and-down-list - only-C2BPnC6.xlsx", skip=1, header=TRUE, stringsAsFactors=FALSE)
signs2Names <- read.xls("All-TGFB-up-and-down-list - only-C2BPnC6.xlsx", nrows=2, header=TRUE, stringsAsFactors=FALSE)

GENESETS2 <- lapply(as.list(signs2), function(x) x[x != ""])
names(GENESETS2) <- colnames(signs2Names)

#-----------------------------------------------------------------
# Run GSVA analysis save for BT
#-----------------------------------------------------------------
ordSamps <- rownames(clinical)
rmSamp <- setdiff(ordSamps, colnames(M2))
ordSamps <- ordSamps[-which(ordSamps==rmSamp)]

geneCore2_ss2b <- gsva(M2, GENESETS2, method="gsva",  mx.diff=1, abs.ranking=FALSE, verbose=FALSE, parallel.sz=1) #gsva
write.csv(geneCore2_ss2b$es.obs[, ordSamps] , "GSVA_All-TGFB-up-and-down-list-only-C2BPnC6.csv")


#-----------------------------------------------------------------
# Generate the heatmaps
#-----------------------------------------------------------------
corrdist = function(x) as.dist(1-cor(t(x), method="pearson"))
hclust.avl = function(x) hclust(x, method="ward.D2")

pal <- c("green", "blue", "red")

mycol <- colorpanel(n=100,low="darkblue",mid="grey50",high="yellow")

grp <- factor(as.character(clinical[ordSamps, 6]), levels=c("green", "blue", "red"))
pathUsed <- signPATHS$BT
paths <- intersect(pathUsed, rownames(geneCore2_ss2b$es.obs))

# FOR PAPER FIGURE ----------------
pdf("GSVA_All-TGFB-AllorderedvF2.pdf", width=8, height=8)
heatmap.2(geneCore2_ss2b$es.obs[paths, ordSamps], scale="none", labCol=clinical[colnames(M2), 1], col=mycol, trace='none', dendrogram="row", Colv=NA, cexRow=0.8, cexCol=0.8, hclustfun=hclust.avl, distfun=corrdist, margins=c(18, 18),ColSideColors=pal[grp], key.title="", density.info="none")
legend("topright", legend=levels(grp), fill=pal, title="", cex=0.7)
dev.off()
# --------------------------


#-----------------------------------------------------------------
# Rerun for LCM Epi 
#-----------------------------------------------------------------
load("Venkata/LCMPARK2-Using38Samples/esetEpiGroups.RData") # loads edataBT

M1 <- edataTE # For epi

# Convert entrez ids to gene symbols
# Get the entrez ids for the gene names
geneSymbols.ID <- mapIds(org.Hs.eg.db, keys=rownames(M1), column="SYMBOL", keytype="ENTREZID", multiVals="first")

geneSymbols.ID1 <- getMissingSymbolIds(geneSymbols.ID)
dannot <- as.data.frame(geneSymbols.ID1)

M2 <- M1
rownames(M2) <- dannot[rownames(M1), 1]

#-----------------------------------------------------------------
# Run GSVA analysis and save results 
#-----------------------------------------------------------------
revSamNames <- sapply(strsplit(colnames(M2),"-"), function(x) x[2])
revSamNames[c(34, 37, 38)] <- c("85", "86", "41")
colnames(M2) <- revSamNames
ordSamps <- rownames(clinical)


geneCore2_ss2b <- gsva(M2, GENESETS2, method="gsva",  mx.diff=1, abs.ranking=FALSE, verbose=FALSE, parallel.sz=1) #gsva
write.csv(geneCore2_ss2b$es.obs[, ordSamps] , "GSVA_All-TGFB-up-and-down-list-only-C2BPnC6_EPI.csv")

# FOR PAPER FIGURE ----------------
grp <- factor(as.character(clinical[ordSamps, 6]), levels=c("green", "blue", "red"))
pathUsed <- signPATHS$EPI

pdf("GSVA_All-TGFB-Allordered_EPIvF2.pdf", width=8, height=8)
heatmap.2(geneCore2_ss2b$es.obs[pathUsed, ordSamps], scale="none", labCol=clinical[colnames(M2), 1], col=mycol, trace='none', dendrogram="row", Colv=NA, cexRow=0.8, cexCol=0.8, hclustfun=hclust.avl, distfun=corrdist, margins=c(18, 18),ColSideColors=pal[grp], key.title="", density.info="none")
legend("topright", legend=levels(grp), fill=pal, title="", cex=0.7)
dev.off()
# --------------------------


#-----------------------------------------------------------------
# Rerun for LCM Stroma
#-----------------------------------------------------------------
load("Venkata/LCMPARK2-Using38Samples/esetStrGroups.RData") # loads edataBT

M1 <- edataTS # For epi

# Convert entrez ids to gene symbols
# Get the entrez ids for the gene names
geneSymbols.ID <- mapIds(org.Hs.eg.db, keys=rownames(M1), column="SYMBOL", keytype="ENTREZID", multiVals="first")

geneSymbols.ID1 <- getMissingSymbolIds(geneSymbols.ID)
dannot <- as.data.frame(geneSymbols.ID1)

M2 <- M1
rownames(M2) <- dannot[rownames(M1), 1]

#-----------------------------------------------------------------
# Run GSVA analysis save results for Tina
#-----------------------------------------------------------------
revSamNames <- sapply(strsplit(colnames(M2),"-"), function(x) x[2])
revSamNames[c(34, 37, 38)] <- c("85", "86", "41")
colnames(M2) <- revSamNames
ordSamps <- rownames(clinical)

geneCore2_ss2b <- gsva(M2, GENESETS2, method="gsva",  mx.diff=1, abs.ranking=FALSE, verbose=FALSE, parallel.sz=1) #gsva
write.csv(geneCore2_ss2b$es.obs[, ordSamps] , "GSVA_All-TGFB-up-and-down-list-only-C2BPnC6_STRb.csv")

pathUsed <- signPATHS$STROMA
path <- intersect(pathUsed, rownames(geneCore2_ss2b$es.obs))

# FOR PAPER ----------------
grp <- factor(as.character(clinical[ordSamps, 6]), levels=c("green", "blue", "red"))

pdf("GSVA_All-TGFB-Allordered_STROMAvF.pdf", width=8, height=8)
heatmap.2(geneCore2_ss2b$es.obs[path, ordSamps], scale="none", labCol=clinical[colnames(M2), 1], col=mycol, trace='none', dendrogram="row", Colv=NA, cexRow=0.8, cexCol=0.8, hclustfun=hclust.avl, distfun=corrdist, margins=c(18, 18),ColSideColors=pal[grp], key.title="", density.info="none")
legend("topright", legend=levels(grp), fill=pal, title="", cex=0.7)
dev.off()
# --------------------------


