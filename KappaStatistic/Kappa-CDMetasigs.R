######################################################################
#
# This code is used to calculate the Kappa statistic
# for each signature along with the combined
# Immune+Fibrosis and Interferon+Cholesterol signatures
# 
######################################################################

library(genefu)
library(xlsx)
require(vcd)

# load the meta signatures
CDmetasigs <- read.xlsx('CD-EDSigs.xlsx',1)
CD1 <- as.character(CDmetasigs$CDSig1[which(is.na(CDmetasigs$CDSig1) == F)])
CD2 <- as.character(CDmetasigs$CDSig2[which(is.na(CDmetasigs$CDSig2) == F)])
CD3 <- as.character(CDmetasigs$CDSig3[which(is.na(CDmetasigs$CDSig3) == F)])
CD4 <- as.character(CDmetasigs$CDSig4[which(is.na(CDmetasigs$CDSig4) == F)])

load('TCellRatio.RData')
Ratio$CD8.based.grouping[which(Ratio$CD8.based.grouping == "black")] <- "green"

load('ExpSetBT.RData')
data <- edataBT
annot = data.frame(rownames(data))
colnames(annot) <- "EntrezGene.ID"
rownames(annot) <- rownames(data)

#immune scores
b <- CD1
b1 <- as.numeric(+1)
x = data.frame(b,b,b1)
colnames(x) <- c("probe","EntrezGene.ID","coefficient")
x$probe <- as.character(x[,1])
x$EntrezGene.ID <- as.character(x[,2])
rownames(x) <- b
scores.imm <- sig.score(x, t(data), annot,do.mapping=TRUE, signed=TRUE, verbose=FALSE)

#fibrosis scores
b <- CD3
b1 <- as.numeric(+1)
x = data.frame(b,b,b1)
colnames(x) <- c("probe","EntrezGene.ID","coefficient")
x$probe <- as.character(x[,1])
x$EntrezGene.ID <- as.character(x[,2])
rownames(x) <- b
scores.fib <- sig.score(x, t(data), annot,do.mapping=TRUE, signed=TRUE, verbose=FALSE)

#meta score of immune and fibrosis
ImScore <- scores.imm$score
ImScore<- (ImScore-median(ImScore))/mad(ImScore)
ImPval <- pnorm(ImScore, lower.tail = FALSE)
FibScore <- scores.fib$score
FibScore <- -(FibScore-median(FibScore))/mad(FibScore)
fibPval <- pnorm(FibScore, lower.tail = FALSE)
PvalMat_ImHi <- cbind(ImPval, fibPval)
z1 <- apply(PvalMat_ImHi, 1, function(x){combine.test(x)})

ImmH.FibL <- cbind(rownames(Ratio),z1)
colnames(ImmH.FibL)[2] <- "Metavalue"
ImmH.FibL <- ImmH.FibL[order(ImmH.FibL[,2]),]
nf <- floor(0.6*37)
ImmH.FibL1 <- ImmH.FibL[1:nf,]

#####
ImScore <- scores.imm$score
ImScore<- -(ImScore-median(ImScore))/mad(ImScore)
ImPval <- pnorm(ImScore, lower.tail = FALSE)
FibScore <- scores.fib$score
FibScore <- (FibScore-median(FibScore))/mad(FibScore)
fibPval <- pnorm(FibScore, lower.tail = FALSE)
PvalMat_Imlow <- cbind(ImPval, fibPval)
z2 <- apply(PvalMat_Imlow, 1, function(x){combine.test(x)})

ImmL.FibH <- cbind(rownames(Ratio),z2)
colnames(ImmL.FibH)[2] <- "Metavalue"
ImmL.FibH <- ImmL.FibH[order(ImmL.FibH[,2]),]
nf1 <- ceiling(0.4*37)
ImmL.FibH1 <- ImmL.FibH[1:nf1,]

# Check if there are any patients repeated in both the groups
y <- intersect(rownames(ImmH.FibL1),rownames(ImmL.FibH1))
ImmH.FibL1 <- ImmH.FibL1[-match(y,rownames(ImmH.FibL1)),]

##### Concordance as kappa statistic

# Kappa for CDmetasig1
Imm.pred <- data.frame(Ratio$Patient.ID,scores.imm$score)
Imm.pred <- Imm.pred[order(-Imm.pred$scores.imm.score),]
Imm.pred[,3] <- NA
colnames(Imm.pred)[3] <- "Predicted Group"
n1 <- floor(0.6*37)
Imm.pred[1:n1,3] <- "red"
Imm.pred[(n1+1):37,3] <- "green"
Imm.pred[,4] <- NA
colnames(Imm.pred)[4] <- "True"
Imm.pred[,4] <- Ratio[rownames(Imm.pred),'CD8.based.grouping']
Imm.pred$True[Imm.pred$True == "blue"] <- "red"
table(Imm.pred$`Predicted Group`,Imm.pred$True)
a = (Kappa(table(Imm.pred$`Predicted Group`,Imm.pred$True)))$Unweighted

# Kappa for CDmetasig3
Fib.pred <- data.frame(Ratio$Patient.ID,scores.fib$score)
Fib.pred <- Fib.pred[order(-Fib.pred$scores.fib.score),]
Fib.pred[,3] <- NA
colnames(Fib.pred)[3] <- "Predicted Group"
n1 <- ceil(0.4*37)
Fib.pred[1:n1,3] <- "green"
Fib.pred[(n1+1):37,3] <- "red"
Fib.pred[,4] <- NA
colnames(Fib.pred)[4] <- "True"
Fib.pred[,4] <- Ratio[rownames(Fib.pred),'CD8.based.grouping']
Fib.pred$True[Fib.pred$True == "blue"] <- "red"
a <- (Kappa(table(Imm.pred$`Predicted Group`,Imm.pred$True)))$Unweighted
summary(a)

# kappa for meta score: CDmetasig1 and CDmetasig3
meta.pred <- data.frame(rbind(ImmH.FibL1,ImmL.FibH1))
meta.pred[,3] <- NA
meta.pred[,4] <- NA
meta.pred[1:nrow(ImmH.FibL1),3] <- "red"
meta.pred[22:nrow(meta.pred),3] <- "green"
meta.pred[,4] <- Ratio[rownames(meta.pred),'CD8.based.grouping']
meta.pred$V4[meta.pred$V4 == "blue"]  <- "red"
colnames(meta.pred)[3] <- "predicted"
colnames(meta.pred)[4] <- "true"
table(meta.pred$predicted,meta.pred$true)
a = (Kappa(table(meta.pred$predicted,meta.pred$true)))$Unweighted
summary(a)

#Kappa for metasig2 
b <- CD2
b1 <- as.numeric(+1)
x = data.frame(b,b,b1)
colnames(x) <- c("probe","EntrezGene.ID","coefficient")
x$probe <- as.character(x[,1])
x$EntrezGene.ID <- as.character(x[,2])
rownames(x) <- b
scores.sig2 <- sig.score(x, t(data), annot,do.mapping=TRUE, signed=TRUE, verbose=FALSE)

sig2.pred <- data.frame(Ratio$Patient.ID,scores.sig2$score)
sig2.pred <- sig2.pred[order(-sig2.pred$scores.sig2.score),]
sig2.pred[,3] <- NA
colnames(sig2.pred)[3] <- "Predicted Group"
n1 <- floor(0.6*37)
sig2.pred[1:n1,3] <- "red"
sig2.pred[(n1+1):37,3] <- "green"
sig2.pred[,4] <- NA
colnames(sig2.pred)[4] <- "True"
sig2.pred[,4] <- Ratio[rownames(sig2.pred),'CD8.based.grouping']
sig2.pred$True[sig2.pred$True == "blue"] <- "red"
table(sig2.pred$`Predicted Group`,sig2.pred$True)
a = (Kappa(table(sig2.pred$`Predicted Group`,sig2.pred$True)))$Unweighted

#Kappa for CDmetasig4
b <- CD4
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
n1 <- floor(0.6*37)
sig4.pred[1:n1,3] <- "red"
sig4.pred[(n1+1):37,3] <- "green"
sig4.pred[,4] <- NA
colnames(sig4.pred)[4] <- "True"
sig4.pred[,4] <- Ratio[rownames(sig4.pred),'CD8.based.grouping']
sig4.pred$True[sig4.pred$True == "blue"] <- "red"
table(sig4.pred$`Predicted Group`,sig4.pred$True)
a = (Kappa(table(sig4.pred$`Predicted Group`,sig4.pred$True)))$Unweighted



