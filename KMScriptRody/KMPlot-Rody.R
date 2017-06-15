######################################################################
# KM plot on the RODY data set
# Meta analysis
# Step1: Combine Immune and Fibrosis
# Split the patients by 60H and 40L
# Step2: Combine Interferon and Cholesterol
# Split by 50H and 50L
#
######################################################################

rm(list=ls())

library(survcomp)
library(genefu)
library(xlsx)

load('ExpSet-Rody.RData')
CDmetasigs <- read.xlsx('CD-EDSigs.xlsx',1)
CD1 <- as.character(CDmetasigs$CDSig1[which(is.na(CDmetasigs$CDSig1) == F)])
CD2 <- as.character(CDmetasigs$CDSig2[which(is.na(CDmetasigs$CDSig2) == F)])
CD3 <- as.character(CDmetasigs$CDSig3[which(is.na(CDmetasigs$CDSig3) == F)])
CD4 <- as.character(CDmetasigs$CDSig4[which(is.na(CDmetasigs$CDSig4) == F)])

EDmetasigs <- read.xlsx('CD-EDSigs.xlsx',2)
ED1 <- as.character(EDmetasigs$EDSig1[which(is.na(EDmetasigs$EDSig1) == F)])
ED2 <- as.character(EDmetasigs$EDSig2[which(is.na(EDmetasigs$EDSig2) == F)])
ED3 <- as.character(EDmetasigs$EDSig3[which(is.na(EDmetasigs$EDSig3) == F)])
ED4 <- as.character(EDmetasigs$EDSig4[which(is.na(EDmetasigs$EDSig4) == F)])
ED5 <- as.character(EDmetasigs$EDSig5[which(is.na(EDmetasigs$EDSig5) == F)])
ED6 <- as.character(EDmetasigs$EDSig6[which(is.na(EDmetasigs$EDSig6) == F)])

load('ExpSet-Rody.RData')
sampleinfo <- read.delim("GSE31519_TNBC_SampleInfo.txt", header=T)
rownames(sampleinfo) <- sampleinfo$SampleNames
tnbc.pdata <- sampleinfo
data <- edata2
annot = data.frame(rownames(data))
colnames(annot) <- "EntrezGene.ID"
rownames(annot) <- rownames(data)
# load the CD1 immune signature
b <- CD1
b1 <- as.numeric(+1)
x = data.frame(b,b,b1)
colnames(x) <- c("probe","EntrezGene.ID","coefficient")
x$probe <- as.character(x[,1])
x$EntrezGene.ID <- as.character(x[,2])
rownames(x) <- b
scores.imm <- sig.score(x, t(data), annot,do.mapping=TRUE, signed=TRUE, verbose=FALSE)

# load the CD3 Fibrosis signature
b <- CD3
b1 <- as.numeric(+1)
x = data.frame(b,b,b1)
colnames(x) <- c("probe","EntrezGene.ID","coefficient")
x$probe <- as.character(x[,1])
x$EntrezGene.ID <- as.character(x[,2])
rownames(x) <- b
scores.fib <- sig.score(x, t(data), annot,do.mapping=TRUE, signed=TRUE, verbose=FALSE)

scores.imm.fib <- data.frame(colnames(data),scores.imm$score,scores.fib$score)
tnbc.pdata[,13] <- scores.imm$score
colnames(tnbc.pdata)[13] <- "ImmPredictedScores"
tnbc.pdata[,14] <- scores.fib$score
colnames(tnbc.pdata)[14] <- "FibPredictedScores"

tnbc.surv <- tnbc.pdata[,c("ev120","fu_120","age","ImmPredictedScores","FibPredictedScores")]
tnbc.surv.complete <- tnbc.surv[complete.cases(tnbc.surv),]

##### Meta analysis of scores

ImScore <- tnbc.surv.complete$ImmPredictedScores
ImScore<- (ImScore-median(ImScore))/mad(ImScore)
ImPval <- pnorm(ImScore, lower.tail = FALSE)
FibScore <- tnbc.surv.complete$FibPredictedScores
FibScore <- -(FibScore-median(FibScore))/mad(FibScore)
fibPval <- pnorm(FibScore, lower.tail = FALSE)
PvalMat_ImHi <- cbind(ImPval, fibPval)
z1 <- apply(PvalMat_ImHi, 1, function(x){combine.test(x)})

ImmH.FibL <- cbind(tnbc.surv.complete,z1)
colnames(ImmH.FibL)[6] <- "Metavalue"
ImmH.FibL <- ImmH.FibL[order(ImmH.FibL$Metavalue),]
nf <- floor(0.6*nrow(tnbc.surv.complete))
ImmH.FibL1 <- ImmH.FibL[1:nf,]

#####
ImScore <- tnbc.surv.complete$ImmPredictedScores
ImScore<- -(ImScore-median(ImScore))/mad(ImScore)
ImPval <- pnorm(ImScore, lower.tail = FALSE)
FibScore <- tnbc.surv.complete$FibPredictedScores
FibScore <- (FibScore-median(FibScore))/mad(FibScore)
fibPval <- pnorm(FibScore, lower.tail = FALSE)
PvalMat_Imlow <- cbind(ImPval, fibPval)
z2 <- apply(PvalMat_Imlow, 1, function(x){combine.test(x)})

ImmL.FibH <- cbind(tnbc.surv.complete,z2)
colnames(ImmL.FibH)[6] <- "Metavalue"
ImmL.FibH <- ImmL.FibH[order(ImmL.FibH$Metavalue),]
nf1 <- ceiling(0.4*nrow(tnbc.surv.complete))
ImmL.FibH1 <- ImmL.FibH[1:nf1,]

# Check if there are any patients repeated in both the groups
y <- intersect(rownames(ImmH.FibL1),rownames(ImmL.FibH1))
ImmH.FibL1 <- ImmH.FibL1[-match(y,rownames(ImmH.FibL1)),]

# Construct the onject for KM

HighImm.LowFib <- ImmH.FibL1
HighImm.LowFib[,7] <- "ImmHFibL"
colnames(HighImm.LowFib)[7] <- "Group"

LowImm.HighFib <- ImmL.FibH1
LowImm.HighFib[,7] <- "ImmLFibH"
colnames(LowImm.HighFib)[7] <- "Group"

tnbc.surv1 <- rbind(HighImm.LowFib[,c(1:5,7)],LowImm.HighFib[,c(1:5,7)])

tnbc.surv1$ev120[which(tnbc.surv1$ev120 == "living" & !is.na(tnbc.surv1$ev120))] <- 0
tnbc.surv1$ev120[which(tnbc.surv1$ev120 == "deceased" & !is.na(tnbc.surv1$ev120))] <- 1
tnbc.surv1$ev120 <- as.numeric(tnbc.surv1$ev120)
tnbc.surv1$fu_120 <- ((tnbc.surv1$fu_120))/10
tnbc.surv1$Group <- factor(tnbc.surv1$Group,levels=c("ImmLFibH","ImmHFibL"))
tnbc.surv1[,7] <- NA
colnames(tnbc.surv1)[7] <- "Color"
tnbc.surv1[which(tnbc.surv1$Group == "ImmLFibH"),7] <- "Green"
tnbc.surv1[which(tnbc.surv1$Group == "ImmHFibL"),7] <- "Red"
tnbc.surv1$Color <- factor(tnbc.surv1$Color,levels=c("Green","Red"))

# Censor for 10 years
tnbc.surv1$fu_120 <- censor.time(as.numeric(tnbc.surv1$fu_120),tnbc.surv1$ev120,time.cens = 10)$surv.time.cens
surv.obj <- survfit(Surv(tnbc.surv1$fu_120,tnbc.surv1$ev120) ~ tnbc.surv1$Group)

# Logrank test
bb <- survdiff(Surv(tnbc.surv1$fu_120,tnbc.surv1$ev120) ~ tnbc.surv1$Group,rho=0)

pdf("MetaAnalysis_Imm&Fib.pdf")
plot(main = "Meta analysis (Immune and Fibrosis): 60%H+40%L", surv.obj,col =levels(tnbc.surv1$Color),lty = 1,lwd = 3,xlab = "Time (years)",
     ylab = "Probability of Relapse-free Survival",xlim=c(0,10))
legend("topright",fill = levels(tnbc.surv1$Color),legend = levels(tnbc.surv1$Group),bty = "n")
legend(0.5,0.4, legend=sprintf("p value (Log rank test)= %.9f",1 - pchisq(bb$chisq, length(bb$n) - 1)),bty="n",cex=0.9,horiz = TRUE)
dev.off()

##### Apply the interferon and cholesterol signature on High score patients

Hscore <- subset(tnbc.surv1,tnbc.surv1$Group == "ImmHFibL")
tnbc.rb <- data[,rownames(Hscore)]
tnbc.pdata1 <- tnbc.pdata[rownames(Hscore),]
data1 <- tnbc.rb
annot = data.frame(rownames(tnbc.rb))
colnames(annot) <- "EntrezGene.ID"
rownames(annot) <- rownames(tnbc.rb)

#interferon and cholesterol scores
b <- ED2
b1 <- as.numeric(+1)
x = data.frame(b,b,b1)
colnames(x) <- c("probe","EntrezGene.ID","coefficient")
x$probe <- as.character(x[,1])
x$EntrezGene.ID <- as.character(x[,2])
rownames(x) <- b
scores.int <- sig.score(x, t(data1), annot,do.mapping=TRUE, signed=TRUE, verbose=FALSE)
b <- ED5
b1 <- as.numeric(+1)
x = data.frame(b,b,b1)
colnames(x) <- c("probe","EntrezGene.ID","coefficient")
x$probe <- as.character(x[,1])
x$EntrezGene.ID <- as.character(x[,2])
rownames(x) <- b
scores.chol <- sig.score(x, t(data1), annot,do.mapping=TRUE, signed=TRUE, verbose=FALSE)

scores.int.chol <- data.frame(colnames(data1),scores.int$score,scores.chol$score)
tnbc.pdata1[,13] <- scores.int$score
colnames(tnbc.pdata1)[13] <- "IntPredictedScores"
tnbc.pdata1[,14] <- scores.chol$score
colnames(tnbc.pdata1)[14] <- "CholPredictedScores"

tnbc.surv <- tnbc.pdata1[,c("ev120","fu_120","age","IntPredictedScores","CholPredictedScores")]
tnbc.surv.complete <- tnbc.surv[complete.cases(tnbc.surv),]

#meta score of interferon and cholesterol
IntScore <- scores.int$score
IntScore<- (IntScore-median(IntScore))/mad(IntScore)
IntPval <- pnorm(IntScore, lower.tail = FALSE)
cholScore <- scores.chol$score
cholScore <- -(cholScore-median(cholScore))/mad(cholScore)
cholPval <- pnorm(cholScore, lower.tail = FALSE)
PvalMat_IntHi <- cbind(IntPval, cholPval)
z1 <- apply(PvalMat_IntHi, 1, function(x){combine.test(x)})

intH.cholL <- cbind(tnbc.surv.complete,z1)
colnames(intH.cholL)[6] <- "Metavalue"
intH.cholL <- intH.cholL[order(intH.cholL$Metavalue),]
nf <- ceiling(0.5*nrow(tnbc.surv.complete))
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

intL.cholH <- cbind(tnbc.surv.complete,z2)
colnames(intL.cholH)[6] <- "Metavalue"
intL.cholH <- intL.cholH[order(intL.cholH$Metavalue),]
nf1 <- ceiling(0.5*nrow(tnbc.surv.complete))
intL.cholH1 <- intL.cholH[1:nf1,]

# Check if there are any patients repeated in both the groups
y <- intersect(rownames(intH.cholL1),rownames(intL.cholH1))
intH.cholL1 <- intH.cholL1[-match(y,rownames(intH.cholL1)),]

# Construct the object for KM
HighInt.LowChol <- intH.cholL1
HighInt.LowChol[,7] <- "IntHCholL"
colnames(HighInt.LowChol)[7] <- "Group"

LowInt.HighChol <- intL.cholH1
LowInt.HighChol[,7] <- "IntLCholH"
colnames(LowInt.HighChol)[7] <- "Group"

tnbc.surv1 <- rbind(HighInt.LowChol[,c(1:5,7)],LowInt.HighChol[,c(1:5,7)])

tnbc.surv1$ev120[which(tnbc.surv1$ev120 == "living" & !is.na(tnbc.surv1$ev120))] <- 0
tnbc.surv1$ev120[which(tnbc.surv1$ev120 == "deceased" & !is.na(tnbc.surv1$ev120))] <- 1
tnbc.surv1$ev120 <- as.numeric(tnbc.surv1$ev120)
tnbc.surv1$fu_120 <- ((tnbc.surv1$fu_120))/10
tnbc.surv1$Group <- factor(tnbc.surv1$Group,levels=c("IntLCholH","IntHCholL"))
tnbc.surv1[,7] <- NA
colnames(tnbc.surv1)[7] <- "Color"
tnbc.surv1[which(tnbc.surv1$Group == "IntLCholH"),7] <- "Green"
tnbc.surv1[which(tnbc.surv1$Group == "IntHCholL"),7] <- "Red"
tnbc.surv1$Color <- factor(tnbc.surv1$Color,levels=c("Green","Red"))

# Censor for 10 years
tnbc.surv1$fu_120 <- censor.time(as.numeric(tnbc.surv1$fu_120),tnbc.surv1$ev120,time.cens = 10)$surv.time.cens
surv.obj <- survfit(Surv(tnbc.surv1$fu_120,tnbc.surv1$ev120) ~ tnbc.surv1$Group)

# Logrank test
bb <- survdiff(Surv(tnbc.surv1$fu_120,tnbc.surv1$ev120) ~ tnbc.surv1$Group,rho=0)

pdf("MetaAnalysis_Int&Chol.pdf")
plot(main = "Meta analysis (Interferon & Cholesterol): 50%H+50%L", surv.obj,col =levels(tnbc.surv1$Color),lty = 1,lwd = 3,xlab = "Time (years)",
     ylab = "Probability of Relapse-free Survival",xlim=c(0,10))
legend("topright",fill = levels(tnbc.surv1$Color),legend = levels(tnbc.surv1$Group),bty = "n")
legend(0.5,0.4, legend=sprintf("p value (Log rank test)= %.9f",1 - pchisq(bb$chisq, length(bb$n) - 1)),bty="n",cex=0.9,horiz = TRUE)
dev.off()


