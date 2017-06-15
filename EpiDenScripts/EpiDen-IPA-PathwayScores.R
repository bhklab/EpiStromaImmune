###############################################################################
#
# This code is  implemented to obtain Metagene Signatures for EpiDensity
# 
# Step1: Leading edge genes are taken from each of the signficiant pathway
# Step2: Pathway score is obtained for each sample (22 high immune samples)
# using sig.score function in genefu R package
# Step3: A hierarchical clustering using the correlation measure as distance 
# and complete linkage is used to cluster the pathway scores for all patients. 
# Step4: To group highly correlated pathway-based signatures, we cut the 
# dendrogram at a height of 0.70 to define metagene signatures
# 
# Result: 6 Metagene signatures for Epi Density
#
###############################################################################

library(piano)
library(snow)
library(genefu)
library(xlsx)
library(pvclust)
library(gplots)

source("plotRunningSum.R")

# load the gsea
load('GSEA-IPAPathways-EpiDen.RData')
gsasummary <- GSAsummaryTable(gsares)
gsasummary <- gsasummary[order(gsasummary[,4]),]

gsea <- xlsx::read.xlsx('GSEA-EpiDen.xlsx',1) 
gsea[which(as.numeric(gsea[,4]) == "0"),4] <- 1/(1000+1) # collapse the p values to a single column
gsea[,5] <- p.adjust(gsea[,4],method="fdr",length(gsea[,4]))

pathways1 <- subset(gsea,gsea[,5] < 0.05)
## remove "Spermine Biosynthesis" as it has oly 2 genes
#a1 <- which(pathways1$Name == "Spermine Biosynthesis")
#pathways1 <- pathways1[-a1,]

leading.edgegenes <- list()
GeneSets <- list()
enrichmentScore <- as.numeric()

for(k in 1:nrow(pathways1)){
 t <- as.character(pathways1$Name[k])
 par(mar = rep(2, 4))
 aa <- plotRunningSum(gsaRes=gsares,geneSet=t,1) # Plot the enrichment score
 leading.edgegenes[[k]] <- aa$leadingEdge
 GeneSets[[k]] <- names(leading.edgegenes[[k]])
 enrichmentScore[k] <- aa$ES 
}
names(GeneSets) <- pathways1$Name

# Convert the gene list  to a dataframe and save it in an excel file
n <- sapply(GeneSets, length)
seq.max <- seq_len(max(n))
GeneSets.mat <- sapply(GeneSets, "[", i = seq.max)
write.xlsx2(GeneSets.mat, file="LeadingEdgeGenes-EpiDensity.xlsx", na="")

map <- read.xlsx2('../HeatMaps-IPAFullList-May16/EpiDen-MetaSignatures/BT-GeneEntID.xlsx',1) 
map$Gene <- as.character(map$Gene)
map$EntID <- as.character(map$EntID)

# Map Ent ID to Gene Symbol

# meta <- read.xlsx2('../HeatMaps-IPAFullList-May16/EpiDen-MetaSignatures/LeadingEdgeGenes-EpiDensity.xlsx',1) 
# Genes <- list()
# for(i in 1:ncol(meta)){
#   t <- which(as.character(meta[,i]) != "")
#   t1 <- as.character(meta[t,i])
#   Genes[[i]] <- map$EntID[match(t1,map$Gene)]
# }
# n <- sapply(Genes, length)
# seq.max <- seq_len(max(n))
# GeneSymbol.mat <- sapply(Genes, "[", i = seq.max)
# write.xlsx2(GeneSymbol.mat, file="LeadingEdgeGeneSymbol-EpiDensity.xlsx", na="")

#piano::geneSetSummary(res, "GUANOSINE NUCLEOTIDES <I>DE NOVO< I> BIOSYNTHESIS%HUMANCYC%PWY-7228")

# length of each group 
v <- sapply(GeneSets, length)
pdf("EpiDen-LeadingEdgeGenes.pdf")
par(mar = rep(4, 4))
barplot(v,ylab="Number of leading-edge genes",main="IPA Pathways-EpiDen Analysis",xlim=c(0,75),ylim=c(0,100),xaxt='n')
mtext(side = 1, text = "Pathways", line = 0.5)
dev.off()

# load the expression set and T cell Ratio
load("/Users/vmanem/Desktop/Project/EpiStroma-Immune/GeneSignature/HeatMaps-IPAFullList-May13/esetBT-LCMBTMapping-RedBlue.RData")
load("/Users/vmanem/Desktop/Project/EpiStroma-Immune/GeneSignature/HeatMaps-IPAFullList-May13//TCellRatio-RB.RData")

annot = data.frame(rownames(edataBT))
colnames(annot) <- "EntrezGene.ID"
rownames(annot) <- rownames(edataBT)
data <- edataBT

scores <- list()
mat <- matrix(NA,length(GeneSets),ncol(edataBT))
for(k in 1:length(GeneSets)){
  
  b <- GeneSets[[k]]
  b1 <- as.numeric(+1)
  
  x = data.frame(b,b,b1)
  colnames(x) <- c("probe","EntrezGene.ID","coefficient")
  x$probe <- as.character(x[,1])
  x$EntrezGene.ID <- as.character(x[,2])
  rownames(x) <- b
  res <- genefu::sig.score(x, t(data), annot,do.mapping=TRUE, signed=TRUE, verbose=TRUE)
  mat[k,] <- res$score
  scores[[k]] <- res$score
}
rownames(mat) <- pathways1$Name
colnames(mat) <- colnames(edataBT)
mat1 <- t(mat)
colnames(mat1) <- pathways1$Name

# CLustering Plot 3
pvcl <- pvclust(data=mat1, method.hclust="complete", method.dist="correlation",use.cor="pairwise.complete.obs", nboot=1000, r=1, store=FALSE)

# plot 2
pdf("EpiDen.pdf",width=200,height=80)
par(oma=c(10,2,10,50)) # bottom,left,top,right
plot(pvcl,cex = 5)
rect.hclust(pvcl$hclust,h=0.7,border=c("red","blue","green","yellow"))
dev.off()

####################################
#
# Cluster Analysis
#
####################################
hc = 0.70 # # Number of clusters 
grps = cutree(pvcl$hclust,h=hc) 
g <- data.frame(table(grps),stringsAsFactors=FALSE)
g1 <- data.frame(grps,stringsAsFactors=FALSE)

g1.genesets <- list()
g1.genesets.final <- list()
for(i in 1:nrow(g)){
  g1 <- data.frame(names(grps[which(grps == as.numeric(g$grps[i]))]))
  ind <- match(g1[,1],names(GeneSets))
  g1.genesets <- list()
  for(k in 1:length(ind)){
    t <- ind[k]
    g1.genesets[[k]] <- GeneSets[[t]]
  }
  g1.genesets.final[[i]] <- unique(Reduce(c,g1.genesets)) # unique genes in all gene sets
}
# length of each group, i.e. the number of unique genes 
l <- sapply(g1.genesets.final, length)
pdf("EpiDen-Cluster.pdf")
par(mar = c(5, 4, 4, 2) + 0.1)
barplot(l,names.arg=seq(1:4),las=0.75,xlab="Cluster",ylab="Number of Unique Genes",main="EpiDen:IPA",ylim=c(0,250))
dev.off()

# Color the clusters- For PAPER
pdf("EpiDen_ColorClusters-Paper.pdf",width=60,height=40)
par(oma=c(10,2,10,10)) # bottom,left,top,right
clusterCols <- c("red","blue","black","magenta","grey","yellow")
myClusterSideBar <- clusterCols[grps]
breaks=0:24
mycol <- colorpanel(n=length(breaks)-1,low="darkblue",mid="grey50",high="yellow")
heatmap.2(mat1,col=mycol,Colv = as.dendrogram(pvcl$hclust),dendrogram = "column", scale = "column",key = FALSE, Rowv = FALSE,
          density.info = "none",trace="none",cexRow=2.5,cexCol=2.0,margins=c(50,60),ColSideColors= myClusterSideBar)
legend("topright",legend = c("MetaSig1","MetaSig2","MetaSig3","MetaSig4","MetaSig5","MetaSig6"),col = unique(myClusterSideBar),lty= 1,lwd = 8,cex=4)
dev.off()

####################################

# Extract meta gene signatures= 6 meta signatures
# Extract all the pathways in that cluster

####################################

Sig1 <- g1.genesets.final[[1]] 
Sig1.pathways <- names(grps[which(grps == "1")])
Sig2 <- g1.genesets.final[[2]]
Sig2.pathways <- names(grps[which(grps == "2")])
Sig3 <- g1.genesets.final[[3]]
Sig3.pathways <- names(grps[which(grps == "3")])
Sig4 <- g1.genesets.final[[4]]
Sig4.pathways <- names(grps[which(grps == "4")])
Sig5 <- g1.genesets.final[[5]]
Sig5.pathways <- names(grps[which(grps == "5")])
Sig6 <- g1.genesets.final[[6]]
Sig6.pathways <- names(grps[which(grps == "6")])

save(Sig1,file="Sig1EpiDen.RData")
save(Sig2,file="Sig2EpiDen.RData")
save(Sig3,file="Sig3EpiDen.RData")
save(Sig4,file="Sig4EpiDen.RData")
save(Sig5,file="Sig5EpiDen.RData")
save(Sig6,file="Sig6EpiDen.RData")

save(Sig1.pathways,file="Sig1EpiDenPathways.RData")
save(Sig2.pathways,file="Sig2EpiDenPathways.RData")
save(Sig3.pathways,file="Sig3EpiDenPathways.RData")
save(Sig4.pathways,file="Sig4EpiDenPathways.RData")
save(Sig5.pathways,file="Sig5EpiDenPathways.RData")
save(Sig6.pathways,file="Sig6EpiDenPathways.RData")


