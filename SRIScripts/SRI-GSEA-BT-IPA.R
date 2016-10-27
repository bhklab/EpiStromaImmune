########################################################################
#
# Gene Set Enrichment Analysis using IPA pathway annotations
# QIAGEN’s Ingenuity® Pathway Analysis (IPA, QIAGEN Redwood City, www.qiagen.com/ingenuity
# A total of 653 pathways are annotated in the IPA
#
# Number of samples = 22
# Number of genes = 19700
# Result: 24 significantly enriched pathways for an FDR < 5%
#
########################################################################

library(Biobase)
library(GEOquery)
library(piano)

# Download the expession set of Bulk Tumor from the GEO
gse <- getGEO('GSE88847',GSEMatrix=TRUE) 
show(gse)

PData <- pData(gse[[1]])
Ratio <- fData(gse[[1]])
edataBT <- exprs(gse[[1]])

rownames(Ratio) <- Ratio$Patient.ID

# select all high immune patients = 22 samples
highimmune <- c("TN1032","TN1005","TN1035","TN1028","TN1004","TN1017","TN1036","TN1030","TN1003","TN1020","TN1025",
                "TN1018","TN1024","TN1015","TN1026","TN1027","TN1007","TN1006","TN1029","TN1019","TN1022","TN1033")

edataBT <- edataBT[,highimmune] # Column names of exp set should be the same as rownames of Ratio
Ratio1 <- Ratio[highimmune,]

cor.val <- as.numeric()
p.val <- as.numeric()
for (i in 1:nrow(edataBT)) {
  tmp <- cor.test(edataBT[i,],Ratio1$Stromal.CD8.Ratio,method="spearman")
  cor.val[i] <- as.numeric(tmp$estimate)
  p.val[i] <- as.numeric(tmp$p.value)
}
Gene <- rownames(edataBT)
metric <- data.frame(Gene,cor.val,p.val)
metric$Gene <- as.character(metric$Gene)
colnames(metric) <- c("Gene","Correlation","Pvalue")

# load the IPA gene sets for the pathways
# Genesets have 2 columns = "Genes" and "IPAPathways"

load('GeneSets-IPA-FinalVersion-EntID.RData')
a <- loadGSC(gSets_IPA_EntID)

stats <- as.vector(metric$Correlation)
names(stats) <- as.vector(metric$Gene)
gsares <- runGSA(geneLevelStats=stats,gsc=a,nPerm=1000,geneSetStat="gsea",adjMethod="none")
save(gsares,file="GSEA-IPAPathways-SRI.RData")

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







