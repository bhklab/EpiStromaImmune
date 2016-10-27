###############################################################################
#
# This code is  implemented to obtain Metagene Signatures for SRI
# Supplementary Figure 9C in the manuscript
#
# Step1: Leading edge genes are taken from each of the signficiant pathway
# Step2: Pathway score is obtained for each sample (22 high immune samples)
# using sig.score function in genefu R package
# Step3: A hierarchical clustering using the correlation measure as distance 
# and complete linkage is used to cluster the pathway scores for all patients. 
# Step4: To group highly correlated pathway-based signatures, we cut the 
# dendrogram at a height of 0.70 to define metagene signatures
# 
# Result: 4 Metagene signatures for SRI
#
###############################################################################

library(Biobase)
library(GEOquery)
library(piano)
library(snow)
library(genefu)
library(xlsx)
library(pvclust)
library(gplots)
library(Hmisc)

source("plotRunningSum.R")

# load the gsea
load('GSEA-IPAPathways-Ratio.RData')
gsasummary <- GSAsummaryTable(gsares)
gsasummary <- gsasummary[order(gsasummary[,4]),]

gseares <- gsasummary
gseares1 <- gseares[,c("Name","Genes (tot)","Stat (dist.dir)","p (dist.dir.up)","p (dist.dir.dn)","Genes (up)","Genes (down)")]
gseares1[,"pval"] <- rowSums(cbind(gseares1[,4],gseares1[,5]), na.rm=TRUE) # collapse the p values of "up and down" to a single column

gseares1[which(as.numeric(gseares1[,"pval"]) == "0"),"pval"] <- 1/(1000+1)  # Replace p value of 0 with: 1/(# of Permutations+1)
gseares1[,"FDR"] <- p.adjust(gseares1[,"pval"],method="fdr",length(gseares1[,"pval"])) # fdr correction
gseares1 <- gseares1[,-(4:5)]
gseaSRI <- subset(gseares1,gseares1[,"FDR"] < 0.05) # 25 pathways are significantly enriched
gseaSRI <- gseaSRI[-which(gseaSRI$'Genes (tot)' < 5),] # remove the pathway that less than 5 genes

pathways1 <- gseaSRI

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

# length of each group 
v <- sapply(GeneSets, length)
pdf("SRI-LeadingEdgeGenes.pdf")
par(mar = rep(4, 4))
barplot(v,ylab="Number of leading-edge genes",main="IPA Pathways-Ratio Analysis",xlim=c(0,30),ylim=c(0,100),xaxt='n')
mtext(side = 1, text = "Pathways", line = 0.5)
dev.off()

# load the expression set of high immune patients (22 samples)
load("ExpBT22Samples.RData")

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

# CLustering
pvcl <- pvclust(data=mat1, method.hclust="complete", method.dist="correlation",use.cor="pairwise.complete.obs", nboot=1000, r=1, store=FALSE)

hc = 0.70 # height of the dendogram
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
pdf("SRIClusters1.pdf")
par(mar = c(5, 4, 4, 2) + 0.1)
barplot(l,names.arg=seq(1:4),las=0.75,xlab="Cluster",ylab="Number of Unique Genes",main="CoreDen:IPA",ylim=c(0,250))
dev.off()

# Heatmap of pathway scores
pdf("SRIClusters-Paper.pdf",width=60,height=40)
par(oma=c(10,2,10,10)) # bottom,left,top,right
clusterCols <- c("red","blue","black","magenta")
myClusterSideBar <- clusterCols[grps]
breaks=0:24
mycol <- colorpanel(n=length(breaks)-1,low="darkblue",mid="grey50",high="yellow")
heatmap.2(mat1,col=mycol,Colv = as.dendrogram(pvcl$hclust),dendrogram = "column", scale = "column",key = FALSE, Rowv = FALSE,
          density.info = "none",trace="none",cexRow=2.5,cexCol=2.0,margins=c(50,60),ColSideColors= myClusterSideBar)
legend("topright",legend = c("MetaSig1","MetaSig2","MetaSig3","MetaSig4"),col = unique(myClusterSideBar),lty= 1,lwd = 8,cex=4)
dev.off()

save(grps,file="SRI-4Clusters.RData")

# Extract meta gene signatures= 4 meta signatures
# Extract all the pathways in each metagene signature

Sig1 <- g1.genesets.final[[1]] 
Sig1.pathways <- names(grps[which(grps == "1")])
Sig2 <- g1.genesets.final[[2]]
Sig2.pathways <- names(grps[which(grps == "2")])
Sig3 <- g1.genesets.final[[3]]
Sig3.pathways <- names(grps[which(grps == "3")])
Sig4 <- g1.genesets.final[[4]]
Sig4.pathways <- names(grps[which(grps == "4")])

save(Sig1,file="Sig1Ratio.RData")
save(Sig2,file="Sig2Ratio.RData")
save(Sig3,file="Sig3Ratio.RData")
save(Sig4,file="Sig4Ratio.RData")

save(Sig1.pathways,file="Sig1RatioPathways.RData")
save(Sig2.pathways,file="Sig2RatioPathways.RData")
save(Sig3.pathways,file="Sig3RatioPathways.RData")
save(Sig4.pathways,file="Sig4RatioPathways.RData")


