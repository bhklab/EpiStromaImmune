plotRunningSum <- function(gsaRes, geneSet, gseaParam=1) {
  
  statistics <- sort(gsaRes$geneLevelStats[,1],decreasing=TRUE)
  selectedStats <- which(names(statistics)%in%gsaRes$gsc[[geneSet]])
  
  S <- selectedStats
  r <- statistics
  p <- gseaParam
  
  m <- length(S)
  N <- length(r)
  NR <- (sum(abs(r[S])^p))
  
  Pprev <- 0
  P <- rep(NA,N)
  for(i in 1:N) {
    if(i%in%S) {
      if(r[i]==0) {
        P[i] <- Pprev + 0
      } else {
        P[i] <- Pprev + abs(r[i])^p/NR
      }
    } else {
      P[i] <- Pprev - 1/(N-m)
    }
    Pprev <- P[i]
  }
  
  if(max(P) > -min(P)) {
    geneSetStatistic <- max(P)
  } else {
    geneSetStatistic <- min(P)        
  }
  
  plot(0,col="white",axes=F,ylab="",xlab="")
  par(fig=c(0,1,0.65,1), new=TRUE)
  plot(P,type="l",ylab=" ",xlab=" ",col="white",axes=F,main=geneSet)
  box()
  abline(v=S)
  if(geneSetStatistic>0) {
    abline(v=S[S<=which(P==geneSetStatistic)],col="red")
  } else {
    abline(v=S[S>=which(P==geneSetStatistic)],col="red")
  }
  
  par(fig=c(0,1,0,0.85), new=TRUE)
  plot(P,type="l",ylab="Running Enrichment Score",xlab="Universal Gene Set (ranked)")
  abline(h=c(0,geneSetStatistic),v=which(P==geneSetStatistic))
  
  if(geneSetStatistic>0) {
    le <- r[S[S<=which(P==geneSetStatistic)]]
  } else {
    le <- r[S[S>=which(P==geneSetStatistic)]]
  }
  
  res <- list(ES=geneSetStatistic,runningSum=P,leadingEdge=le)
  
}
