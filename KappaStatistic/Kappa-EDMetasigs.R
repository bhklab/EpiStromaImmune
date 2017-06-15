######################################################################
#
# This code is used to calculate the Kappa statistic
# for each signature along with the combined
# Immune+Fibrosis and Interferon+Cholesterol signatures
# 
######################################################################

library(genefu)
library(xlsx)

# load the meta signatures
EDmetasigs <- read.xlsx('CD-EDSigs.xlsx',2)
ED1 <- as.character(EDmetasigs$EDSig1[which(is.na(EDmetasigs$EDSig1) == F)])
ED2 <- as.character(EDmetasigs$EDSig2[which(is.na(EDmetasigs$EDSig2) == F)])
ED3 <- as.character(EDmetasigs$EDSig3[which(is.na(EDmetasigs$EDSig3) == F)])
ED4 <- as.character(EDmetasigs$EDSig4[which(is.na(EDmetasigs$EDSig4) == F)])
ED5 <- as.character(EDmetasigs$EDSig5[which(is.na(EDmetasigs$EDSig5) == F)])
ED6 <- as.character(EDmetasigs$EDSig6[which(is.na(EDmetasigs$EDSig6) == F)])

load('TCellRatio-RB.RData')

load('ExpSetBT.RData')
edataBT <- edataBT[,rownames(Ratio)]

data <- edataBT
annot = data.frame(rownames(data))
colnames(annot) <- "EntrezGene.ID"
rownames(annot) <- rownames(data)

#interferon scores
b <- ED2
b1 <- as.numeric(+1)
x = data.frame(b,b,b1)
colnames(x) <- c("probe","EntrezGene.ID","coefficient")
x$probe <- as.character(x[,1])
x$EntrezGene.ID <- as.character(x[,2])
rownames(x) <- b
scores.int <- sig.score(x, t(data), annot,do.mapping=TRUE, signed=TRUE, verbose=FALSE)

#cholesterol scores
b <- ED5
b1 <- as.numeric(+1)
x = data.frame(b,b,b1)
colnames(x) <- c("probe","EntrezGene.ID","coefficient")
x$probe <- as.character(x[,1])
x$EntrezGene.ID <- as.character(x[,2])
rownames(x) <- b
scores.chol <- sig.score(x, t(data), annot,do.mapping=TRUE, signed=TRUE, verbose=FALSE)

#meta score of interferon and cholesterol
IntScore <- scores.int$score
IntScore<- (IntScore-median(IntScore))/mad(IntScore)
IntPval <- pnorm(IntScore, lower.tail = FALSE)
cholScore <- scores.chol$score
cholScore <- -(cholScore-median(cholScore))/mad(cholScore)
cholPval <- pnorm(cholScore, lower.tail = FALSE)
PvalMat_IntHi <- cbind(IntPval, cholPval)
z1 <- apply(PvalMat_IntHi, 1, function(x){combine.test(x)})

intH.cholL <- cbind(rownames(Ratio),z1)
colnames(intH.cholL)[2] <- "Metavalue"
intH.cholL <- intH.cholL[order(intH.cholL[,2]),]
nf <- ceiling(0.454*22)
intH.cholL1 <- intH.cholL[1:nf,]

#####
IntScore <- scores.int$score
IntScore<- -(IntScore-median(IntScore))/mad(IntScore)
IntPval <- pnorm(IntScore, lower.tail = FALSE)
cholScore <- scores.chol$score
cholScore <- (cholScore-median(cholScore))/mad(cholScore)
cholPval <- pnorm(cholScore, lower.tail = FALSE)
PvalMat_Intlow <- cbind(IntPval, cholPval)
z2 <- apply(PvalMat_Intlow, 1, function(x){combine.test(x)})

intL.cholH <- cbind(rownames(Ratio),z2)
colnames(intL.cholH)[2] <- "Metavalue"
intL.cholH <- intL.cholH[order(intL.cholH[,2]),]
nf1 <- ceiling(0.545*22)
intL.cholH1 <- intL.cholH[1:nf1,]

# Check if there are any patients repeated in both the groups
y <- intersect(rownames(intH.cholL1),rownames(intL.cholH1))
intH.cholL1 <- intH.cholL1[-match(y,rownames(intH.cholL1)),]

##### Concordance as kappa statistic

# Kappa for ED2 metasig
int.pred <- data.frame(Ratio$Patient.ID,scores.int$score)
int.pred <- int.pred[order(-int.pred$scores.int.score),]
int.pred[,3] <- NA
colnames(int.pred)[3] <- "Predicted Group"
n1 <- ceiling(0.454*22)
int.pred[1:n1,3] <- "red"
int.pred[(n1+1):22,3] <- "blue"
int.pred[,4] <- NA
colnames(int.pred)[4] <- "True"
int.pred[,4] <- Ratio[rownames(int.pred),'CD8.based.grouping']
table(int.pred$`Predicted Group`,int.pred$True)
a = (Kappa(table(int.pred$`Predicted Group`,int.pred$True)))$Unweighted

# Kappa for ED5 metasig
chol.pred <- data.frame(Ratio$Patient.ID,scores.chol$score)
chol.pred <- chol.pred[order(-chol.pred$scores.chol.score),]
chol.pred[,3] <- NA
colnames(chol.pred)[3] <- "Predicted Group"
n1 <- ceiling(0.545*22)
chol.pred[1:n1,3] <- "blue"
chol.pred[(n1+1):22,3] <- "red"
chol.pred[,4] <- NA
colnames(chol.pred)[4] <- "True"
chol.pred[,4] <- Ratio[rownames(chol.pred),'CD8.based.grouping']
table(chol.pred$`Predicted Group`,chol.pred$True)
a = (Kappa(table(chol.pred$`Predicted Group`,chol.pred$True)))$Unweighted

# Kappa for meta score of ED2, ED5 metasig
meta.pred <- data.frame(rbind(intH.cholL1,intL.cholH1))
meta.pred[,3] <- NA
meta.pred[,4] <- NA
meta.pred[1:nrow(intH.cholL1),3] <- "red"
meta.pred[(nrow(intH.cholL1)+1):nrow(meta.pred),3] <- "blue"
meta.pred[,4] <- Ratio[rownames(meta.pred),'CD8.based.grouping']

colnames(meta.pred)[3] <- "predicted"
colnames(meta.pred)[4] <- "true"
table(meta.pred$predicted,meta.pred$true)
a = (Kappa(table(meta.pred$predicted,meta.pred$true)))$Unweighted

# Kappa for ED1
b <- ED1
b1 <- as.numeric(+1)
x = data.frame(b,b,b1)
colnames(x) <- c("probe","EntrezGene.ID","coefficient")
x$probe <- as.character(x[,1])
x$EntrezGene.ID <- as.character(x[,2])
rownames(x) <- b
scores.sig1 <- sig.score(x, t(data), annot,do.mapping=TRUE, signed=TRUE, verbose=FALSE)
sig1.pred <- data.frame(Ratio$Patient.ID,scores.sig1$score)
sig1.pred <- sig1.pred[order(-sig1.pred$scores.sig1.score),]
sig1.pred[,3] <- NA
colnames(sig1.pred)[3] <- "Predicted Group"
n1 <- ceiling(0.454*22)
sig1.pred[1:n1,3] <- "red"
sig1.pred[(n1+1):22,3] <- "blue"
sig1.pred[,4] <- NA
colnames(sig1.pred)[4] <- "True"
sig1.pred[,4] <- Ratio[rownames(sig1.pred),'CD8.based.grouping']
table(sig1.pred$`Predicted Group`,sig1.pred$True)
a = (Kappa(table(sig1.pred$`Predicted Group`,sig1.pred$True)))$Unweighted

# Kappa for ED3
b <- ED3
b1 <- as.numeric(+1)
x = data.frame(b,b,b1)
colnames(x) <- c("probe","EntrezGene.ID","coefficient")
x$probe <- as.character(x[,1])
x$EntrezGene.ID <- as.character(x[,2])
rownames(x) <- b
scores.sig3 <- sig.score(x, t(data), annot,do.mapping=TRUE, signed=TRUE, verbose=FALSE)
sig3.pred <- data.frame(Ratio$Patient.ID,scores.sig3$score)
sig3.pred <- sig3.pred[order(-sig3.pred$scores.sig3.score),]
sig3.pred[,3] <- NA
colnames(sig3.pred)[3] <- "Predicted Group"
n1 <- ceiling(0.454*22)
sig3.pred[1:n1,3] <- "red"
sig3.pred[(n1+1):22,3] <- "blue"
sig3.pred[,4] <- NA
colnames(sig3.pred)[4] <- "True"
sig3.pred[,4] <- Ratio[rownames(sig3.pred),'CD8.based.grouping']
table(sig3.pred$`Predicted Group`,sig3.pred$True)
a = (Kappa(table(sig3.pred$`Predicted Group`,sig3.pred$True)))$Unweighted

# Kappa for ED4
b <- ED4
b1 <- as.numeric(+1)
x = data.frame(b,b,b1)
colnames(x) <- c("probe","EntrezGene.ID","coefficient")
x$probe <- as.character(x[,1])
x$EntrezGene.ID <- as.character(x[,2])
rownames(x) <- b
scores.sig4 <- sig.score(x, t(data), annot,do.mapping=TRUE, signed=TRUE, verbose=FALSE)
sig4.pred <- data.frame(Ratio$Patient.ID,scores.sig4$score)
sig4.pred <- sig4.pred[order(-sig4.pred$scores.sig4.score),]
sig4.pred[,3] <- NA
colnames(sig4.pred)[3] <- "Predicted Group"
n1 <- ceiling(0.454*22)
sig4.pred[1:n1,3] <- "red"
sig4.pred[(n1+1):22,3] <- "blue"
sig4.pred[,4] <- NA
colnames(sig4.pred)[4] <- "True"
sig4.pred[,4] <- Ratio[rownames(sig4.pred),'CD8.based.grouping']
table(sig4.pred$`Predicted Group`,sig4.pred$True)
a = (Kappa(table(sig4.pred$`Predicted Group`,sig4.pred$True)))$Unweighted

# Kappa for ED6
b <- ED6
b1 <- as.numeric(+1)
x = data.frame(b,b,b1)
colnames(x) <- c("probe","EntrezGene.ID","coefficient")
x$probe <- as.character(x[,1])
x$EntrezGene.ID <- as.character(x[,2])
rownames(x) <- b
scores.sig6 <- sig.score(x, t(data), annot,do.mapping=TRUE, signed=TRUE, verbose=FALSE)
sig6.pred <- data.frame(Ratio$Patient.ID,scores.sig6$score)
sig6.pred <- sig6.pred[order(-sig6.pred$scores.sig6.score),]
sig6.pred[,3] <- NA
colnames(sig6.pred)[3] <- "Predicted Group"
n1 <- ceiling(0.454*22)
sig6.pred[1:n1,3] <- "red"
sig6.pred[(n1+1):22,3] <- "blue"
sig6.pred[,4] <- NA
colnames(sig6.pred)[4] <- "True"
sig6.pred[,4] <- Ratio[rownames(sig6.pred),'CD8.based.grouping']
table(sig6.pred$`Predicted Group`,sig6.pred$True)
a = (Kappa(table(sig6.pred$`Predicted Group`,sig6.pred$True)))$Unweighted



