###############################################################################
#
# This code is implemented to obtain the KM plot using CoreDen-Metagene1 
# Rody Data Set: GSE31519
# Survival Curve in Supplementary Figure 9d in the manuscript
# Script for Preprocessing of Rody dataset is also available in
# Rody, Achim, et al. "A clinically relevant gene signature in triple negative 
# and basal-like breast cancer." Breast Cancer Research 13.5 (2011): 1.
#
# Download the data from GSE31519
# 1. "GSE31519_complete_dataset.txt"
# 2. "GSE31519_TNBC_SampleInfo.txt" = 579 samples 
# 3. "GPL96-57554.txt" = Annotation file
#
# Step1: Use CoreDen-MetaSig1 and obtain score for each patient
# using sig.score function in genefu R package
# Step2: Split the patient groups as Low and high immune based on 1/3 and 2/3
# 
#
###############################################################################

library(genefu)
library(gplots)
library(Hmisc)
library(survcomp)

# Script from Rody's paper
samples <- read.delim("GSE31519_complete_dataset_new.txt", colClasses= "character", skip=1, nrows=1, header=FALSE)   # read the sample names from row-2
data <- read.delim("GSE31519_complete_dataset_new.txt", skip=4,header=FALSE)  # read data starting row-5
colnames(data)=t(samples)       # add sample names from row-2 as column names
data=data[,!samples==""]        # remove blank gap-columns with no sample-title in row-2 from expr data
rm(samples)   # delete used variable
rownames(data)=data[,1]     # use first column with Affy IDs as rownames
cd_compl_579=data[,-1]    # remove first column with Affy IDs
rm(data)  # delete used variable
cd_compl_579=cd_compl_579[,sort(colnames(cd_compl_579))]   # sort columns by sample name

# Renove the duplicate probes by selecting the most IQR

a <- as.numeric()
for (i in 1:nrow(cd_compl_579))
{  
  a[i] <- sum(is.na(cd_compl_579[i,]))
}

data <- cd_compl_579
data <- data[-which(a>0),]

data1 <- data
data1 <- as.matrix(data1)
iqr <- as.numeric()
for (i in 1:nrow(data1))
{ 
  iqr[i] <- IQR(data1[i,]) 
}

annot <- read.delim("GPL96-57554.txt",header=T)
bestprobe <- cbind(rownames(data1),iqr)
Ent <- as.character()
for (i in 1:nrow(data1))
{ 
  Ent[i] <- as.character(annot$ENTREZ_GENE_ID[match(bestprobe[i,1],annot$ID)])
}
bestprobe <- cbind(rownames(data1),iqr,Ent)

bestprobe <- data.frame(bestprobe)
colnames(bestprobe) <- c("Probe","IQR","Ent")

genes <- as.character(unique(bestprobe$Ent))
bestprobe_final <- as.character()
for (i in 1:length(genes))
{ 
  a <- which(bestprobe$Ent == genes[i])
  id <- which.max(bestprobe$IQR[a])
  bestprobe_final[i] <- as.character(bestprobe$Probe[a[id]])
}

edata <- data1[match(bestprobe_final,rownames(data1)),]
rnames <- as.character(bestprobe$Ent[match(rownames(edata),bestprobe$Probe)])
edata1 <- edata
rownames(edata1) <- rnames

# remove all those probes that have more than one entrez gene id
c <- grep("///",rownames(edata1))
edata2 <- edata1[-c,]

###############################################
#
# KM plot: High immune predicts good outcome
#
###############################################

load('ExpSet-Rody.RData')

sampleinfo <- read.delim("GSE31519_TNBC_SampleInfo.txt", header=T) 
rownames(sampleinfo) <- sampleinfo$SampleNames

tnbc.pdata <- sampleinfo
data <- edata2 # eset of Rody
annot = data.frame(rownames(edata2))
colnames(annot) <- "EntrezGene.ID"
rownames(annot) <- rownames(edata2)

# load the signature
load("../CoreDen-MetaSignatures/Sig1CoreDen.RData")
b <- Sig1
b1 <- as.numeric(+1)
x = data.frame(b,b,b1)
colnames(x) <- c("probe","EntrezGene.ID","coefficient")
x$probe <- as.character(x[,1])
x$EntrezGene.ID <- as.character(x[,2])
rownames(x) <- b

scores <- sig.score(x, t(data), annot,do.mapping=TRUE, signed=TRUE, verbose=TRUE)
scores.tnbc <- data.frame(colnames(edata2),scores$score)

tnbc.pdata[,13] <- scores$score
colnames(tnbc.pdata)[13] <- "PredictedScores"

chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 
test <- split(scores$score, cut2(scores$score, g=3))

lowscore.tnbc <- scores.tnbc[match(names(test[[1]]),rownames(scores.tnbc)),]
interscore.tnbc <- scores.tnbc[match(names(test[[2]]),rownames(scores.tnbc)),]
highscore.tnbc <- scores.tnbc[match(names(test[[3]]),rownames(scores.tnbc)),]

tnbc.surv <- matrix(NA,ncol(data),5)
tnbc.surv <- tnbc.pdata[,c("ev120","fu_120","age","PredictedScores")]
tnbc.surv$ev120[which(tnbc.surv$ev120 == "living" & !is.na(tnbc.surv$ev120))] <- 0
tnbc.surv$ev120[which(tnbc.surv$ev120 == "deceased" & !is.na(tnbc.surv$ev120))] <- 1
tnbc.surv$ev120 <- as.numeric(tnbc.surv$ev120)
tnbc.surv[rownames(lowscore.tnbc),5] <- "LowScore"
tnbc.surv[rownames(interscore.tnbc),5] <- "HighScore"
tnbc.surv[rownames(highscore.tnbc),5] <- "HighScore"
colnames(tnbc.surv)[5] <- "Class"
tnbc.surv.complete <- tnbc.surv[complete.cases(tnbc.surv),]
days.per.month <- 30.4368
days.per.year <- 365.242

tnbc.surv.complete$fu_120 <- ((tnbc.surv.complete$fu_120))/10
tnbc.surv.complete$Class <- factor(tnbc.surv.complete$Class,levels=c("LowScore","HighScore"))
tnbc.surv.complete[,6] <- NA 
colnames(tnbc.surv.complete)[6] <- "Color"
tnbc.surv.complete[which(tnbc.surv.complete$Class == "LowScore"),6] <- "Green"
tnbc.surv.complete[which(tnbc.surv.complete$Class == "HighScore"),6] <- "Red"
tnbc.surv.complete$Color <- factor(tnbc.surv.complete$Color,levels=c("Green","Red"))

surv.obj <- survfit(Surv(tnbc.surv.complete$fu_120,tnbc.surv.complete$ev120) ~ tnbc.surv.complete$Class)

bb <- survdiff(Surv(tnbc.surv.complete$fu_120,tnbc.surv.complete$ev120) ~ tnbc.surv.complete$Class,rho=0) # logrank test

######################################
#
# D.index and Concordance Index
#
######################################

xx <- scores.tnbc$scores.score
stime1 <- tnbc.surv$fu_120
sevent1 <- tnbc.surv$ev120
dind <- D.index(x=xx, surv.time=stime1, surv.event=sevent1,na.rm = TRUE)
#print(dind)
#dind$d.index
#dind$p.value

conindex <- concordance.index(x=xx, surv.time=stime1,surv.event=sevent1,method = "noether",na.rm=TRUE)
#print(conindex)
#conindex$c.index
#conindex$p.value

pdf("Immune-RodyDataSet.pdf")
plot(main = "CDMetaSig1-Rody", surv.obj,col =levels(tnbc.surv.complete$Color),lty = 1,lwd = 3,xlab = "Time (years)",ylab = "Probability of Relapse-free Survival",xlim=c(0,10))
legend("topright",fill = levels(tnbc.surv.complete$Color),legend = levels(tnbc.surv.complete$Class),bty = "n")
legend(0.5,0.4, legend=sprintf("p value (Log rank test)= %.4f",1 - pchisq(bb$chisq, length(bb$n) - 1)),bty="n",cex=0.9,horiz = TRUE)
legend(0.5,0.25, legend=sprintf("D Index= %.3f, p value= %.5f",dind$d.index,dind$p.value),bty="n",cex=0.9,horiz = TRUE)
legend(0.5,0.20, legend=sprintf("Concordace Index= %.3f, p value= %.5f",conindex$c.index,conindex$p.value),bty="n",cex=0.9,horiz = TRUE)
dev.off()

######################################
#
# D.index for all the 8 metasignatures
#
######################################

load('ExpSet-Rody.RData')

sampleinfo <- read.delim("GSE31519_TNBC_SampleInfo.txt", header=T) 
rownames(sampleinfo) <- sampleinfo$SampleNames

tnbc.pdata <- sampleinfo
data <- edata2 # eset of Rody
annot = data.frame(rownames(edata2))
colnames(annot) <- "EntrezGene.ID"
rownames(annot) <- rownames(edata2)

# load the signature
load("../CoreDen-MetaSignatures/Sig1CoreDen.RData")
load("../CoreDen-MetaSignatures/Sig2CoreDen.RData")
load("../CoreDen-MetaSignatures/Sig3CoreDen.RData")
load("../CoreDen-MetaSignatures/Sig4CoreDen.RData")
load("../SRI-MetaSignatures/Sig1SRI.RData")
load("../SRI-MetaSignatures/Sig2SRI.RData")
load("../SRI-MetaSignatures/Sig3SRI.RData")
load("../SRI-MetaSignatures/Sig4SRI.RData")

g <- list(CDSig1,CDSig2,CDSig3,CDSig4,SRISig1,SRISig2,SRISig3,SRISig4)

for (i in 1:length(g)) {
  
  b <- g[[i]]
  b1 <- as.numeric(+1)
  x = data.frame(b,b,b1)
  colnames(x) <- c("probe","EntrezGene.ID","coefficient")
  x$probe <- as.character(x[,1])
  x$EntrezGene.ID <- as.character(x[,2])
  rownames(x) <- b
  data <- tnbc.scaled.eset
  annot = data.frame(rownames(data))
  colnames(annot) <- "EntrezGene.ID"
  rownames(annot) <- rownames(data)
  scores <- sig.score(x, t(data), annot,do.mapping=TRUE, signed=TRUE, verbose=TRUE)
  
  xx <- scores$score
  stime1 <- tnbc.surv$fu_120
  sevent1 <- tnbc.surv$ev120
  dind <- D.index(x=xx, surv.time=stime1, surv.event=sevent1,na.rm = TRUE)
  
  }

