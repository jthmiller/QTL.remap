#################################################################################################
##                                    QTL Mapping Program                                      ##
##                                                                                             ##
##    This R code of analyzing selective genotyping data is programmed for two methods, one is ##
##  the method by Lee et al. (2014), and the other is the method using mixture truncated       ##
##  model, for QTL detection. It allows the analysis of the selective genotyping data for QTL  ##
##  detection, including the estimation of positions and effects (additive, dominance and      ##
##  epistasis) of QTL, in the F2 population. The method by Lee et al. is based on the normal   ##
##  mixture model when both genotyped and ungenotyped individuals are fitted in the model for  ##
##  selective genotyping analysis. The method using mixture truncated model is to analyze      ##
##  selective genotyping data when only the genotyped individuals are considered in the        ##
##  analysis. This program is composed of five functions. Function 5 is written by Mr.         ##
##  Hsiang-An Ho for the mixture truncated model. The demonstration of using the code for      ##
##  analyzing a simulated data is given below.                                                 ##
##                                                                                             ##
##                                                   Author: Hsin-I Lee                        ##
##                                                          Institute of Statistical Science,  ##
##                                                          Academia Sinica, Taipei, Taiwan    ##
#################################################################################################

##### Functions #################################################################################
##                                                                                             ##

Design<-function(nQTL, domi=TRUE, epis=TRUE){
   da<-c(1,0,-1)
   dd<-c(-1/2,1/2,-1/2)
   DsA=DsD<-matrix(0,3^nQTL,nQTL)
   for(i in 1:nQTL){
     DsA[,i]<-rep(rep(da,3^(nQTL-i)),each=3^(i-1))
     DsD[,i]<-rep(rep(dd,3^(nQTL-i)),each=3^(i-1))
     }
   if(domi==TRUE)
     {Dsn<-cbind(DsA,DsD)}else{
      Dsn<-DsA}
   if(nQTL>1 & epis==TRUE)
     {DsIA=DsID<-matrix(0,3^nQTL,choose(nQTL,2))
      DsIAD<-matrix(0,3^nQTL,nQTL*(nQTL-1))
      i.a=i.ad<-0
      for(j in 1:nQTL){
        for(k in 1:nQTL){
          if(k>j)
            {i.a<-i.a+1
             DsIA[,i.a]<-DsA[,j]*DsA[,k]
             DsID[,i.a]<-DsD[,j]*DsD[,k]}
          if(j!=k)
            {i.ad<-i.ad+1
             DsIAD[,i.ad]<-DsA[,j]*DsD[,k]}
          }
        }
      if(domi==TRUE)
        {Dsn<-cbind(Dsn,DsIA,DsID,DsIAD)}else{
         Dsn<-cbind(Dsn,DsIA)}}
   return(Dsn)
   }

QTLfreq<-function(p, marker.map, marker.gty){
   mar.d<-marker.map[[p[1]]][p[2]]
   mgty<-marker.gty[[p[1]]][,p[2]:(p[2]+1)]
   Ns<-nrow(mgty)
   Q.freq<-matrix(0,Ns,3)
   if(p[3]==0)
     {for(n in 1:Ns){
        if(mgty[n,1]==2)
          {Q.freq[n,]<-c(1,0,0)}
        if(mgty[n,1]==1)
          {Q.freq[n,]<-c(0,1,0)}
        if(mgty[n,1]==0)
          {Q.freq[n,]<-c(0,0,1)}
        }}
   if(p[3]==mar.d)
     {for(n in 1:Ns){
        if(mgty[n,2]==2)
          {Q.freq[n,]<-c(1,0,0)}
        if(mgty[n,2]==1)
          {Q.freq[n,]<-c(0,1,0)}
        if(mgty[n,2]==0)
          {Q.freq[n,]<-c(0,0,1)}
        }}
   if(sum(p[3]==c(0,mar.d))==0)
     {rm<-(1-exp(-2*mar.d))/2
      rq<-(1-exp(-2*p[3]))/2
      cq<-(rm^2)/((rm^2)+((1-rm)^2))
      pq<-rq/rm
      for(n in 1:Ns){
        if(mgty[n,1]==2 & mgty[n,2]==2)
          {Q.freq[n,]<-c(1,0,0)}
        if(mgty[n,1]==2 & mgty[n,2]==1)
          {Q.freq[n,]<-c(1-pq,pq,0)}
        if(mgty[n,1]==2 & mgty[n,2]==0)
          {Q.freq[n,]<-c((1-pq)^2,2*pq*(1-pq),pq^2)}
        if(mgty[n,1]==1 & mgty[n,2]==2)
          {Q.freq[n,]<-c(pq,1-pq,0)}
        if(mgty[n,1]==1 & mgty[n,2]==1)
          {Q.freq[n,]<-c(cq*pq*(1-pq),1-2*cq*pq*(1-pq),cq*pq*(1-pq))}
        if(mgty[n,1]==1 & mgty[n,2]==0)
          {Q.freq[n,]<-c(0,1-pq,pq)}
        if(mgty[n,1]==0 & mgty[n,2]==2)
          {Q.freq[n,]<-c(pq^2,2*pq*(1-pq),(1-pq)^2)}
        if(mgty[n,1]==0 & mgty[n,2]==1)
          {Q.freq[n,]<-c(0,pq,1-pq)}
        if(mgty[n,1]==0 & mgty[n,2]==0)
          {Q.freq[n,]<-c(0,0,1)}
        }}
   return(Q.freq)
   }

Mixprop<-function(nQTL, nu=length(yu), posi, marker.map, marker.gty, model){
   QTL.freq<-QTLfreq(posi[[1]],marker.map,marker.gty)
   if(nQTL>=2)
     {for(m in 2:nQTL){
        QTL.f<-QTLfreq(posi[[m]],marker.map,marker.gty)
        QTL.freq<-cbind(QTL.f[,1]*QTL.freq,QTL.f[,2]*QTL.freq,QTL.f[,3]*QTL.freq)
        }}
   if(model==2)
     {mix.prop<-QTL.freq}
   if(model==1)
     {N<-nrow(QTL.freq)+nu
      freq.s<-colSums(QTL.freq)/N
      if(nQTL==1)
        {popu.freq<-c(0.25,0.5,0.25)}
      if(nQTL>=2)
        {M<-c(1,0)
         GAM<-matrix(0,2^nQTL,nQTL)
         Pg<-c(2,1,0)
         PG<-matrix(0,3^nQTL,nQTL)
         rmn<-rep(0,nQTL-1)
         for(ng in 1:nQTL){
           GAM[,ng]<-rep(rep(M,2^(nQTL-ng)),each=2^(ng-1))
           PG[,ng]<-rep(rep(Pg,3^(nQTL-ng)),each=3^(ng-1))
           if(ng<nQTL)
             {if(posi[[c(ng,1)]]==posi[[c(ng+1,1)]])
                {rmd<-sum(marker.map[[posi[[c(ng,1)]]]][posi[[c(ng,2)]]:(posi[[c(ng+1,2)]]-1)])+
                      posi[[c(ng+1,3)]]-posi[[c(ng,3)]]
                 rmn[ng]<-(1-exp(-2*rmd))/2}else{
                 rmn[ng]<-0.5}}
           }
         gamf=GAM.f<-matrix(0,2^nQTL,nQTL-1)
         GAM.freq<-rep(1,2^nQTL)
         for(nr in 1:(nQTL-1)){
           gamf[,nr]<-GAM[,nr]==GAM[,nr+1]
           for(ngam in 1:2^nQTL){
             GAM.f[ngam,nr]<-ifelse(gamf[ngam,nr]==1,1-rmn[nr],rmn[nr])
             }
           GAM.freq<-GAM.freq*GAM.f[,nr]
           }
         Gfreq<-rep(0,(2^nQTL)*(2^nQTL))
         Gtype<-matrix(0,(2^nQTL)*(2^nQTL),nQTL)
         wg<-0
         for(ga1 in 1:2^nQTL){
           for(ga2 in 1:2^nQTL){
             wg<-wg+1
             Gfreq[wg]<-(GAM.freq[ga1]/2)*(GAM.freq[ga2]/2)
             Gtype[wg,]<-GAM[ga1,]+GAM[ga2,]
             }
           }
         popu.freq<-rep(0,3^nQTL)
         ind<-rep(0,(2^nQTL)*(2^nQTL))
         for(pfi in 1:3^nQTL){
           for(gti in 1:((2^nQTL)*(2^nQTL))){
             if(sum(Gtype[gti,]==PG[pfi,])==nQTL)
               {ind[gti]<-pfi}
             }
           popu.freq[pfi]<-sum(Gfreq[ind==pfi])
           }}
      freq.u<-popu.freq-freq.s
      if(sum(freq.u<0)>0)
        {for(ifu in 1:(3^nQTL)){
           if(freq.u[ifu]<0)
             {freq.u[ifu]<-0}
           }}
      Freq.u<-matrix(rep(freq.u/sum(freq.u),each=nu),nu,3^nQTL)
      mix.prop<-rbind(QTL.freq,Freq.u)}
   return(mix.prop)
   }

QTL.estim.p<-function(nQTL, posi, ys, yu, D.matrix, marker.map, marker.gty, EM=1){
   Y<-c(ys,yu)
   N<-length(Y)
   nu<-length(yu)
   Freq<-Mixprop(nQTL, nu, posi, marker.map, marker.gty, model=1)
   n.para<-ncol(D.matrix)
   if(EM==1)
     {L0<-sum(log(dnorm(Y,mean=mean(Y),sd=var(Y)^0.5)))
      mut<-as.vector(mean(Y))
      Et<-rep(0,n.para)
      sigt<-var(Y)
      MUt<-c(rep(mut,3^nQTL))+D.matrix%*%matrix(Et,n.para,1)
      PIup<-matrix(0,N,3^nQTL)
      for(i.p in 1:3^nQTL){
        PIup[,i.p]<-Freq[,i.p]*dnorm(Y,mean=MUt[i.p],sd=sigt^0.5)
        }
      PI<-PIup/rowSums(PIup)
      R<-matrix(0,n.para,1)
      M=V<-matrix(0,n.para,n.para)
      one<-rep(1,N)
      for(m1 in 1:n.para){
        R[m1,]<-(t(Y-mut)%*%PI%*%D.matrix[,m1])/(one%*%PI%*%(D.matrix[,m1]*D.matrix[,m1]))
        for(m2 in 1:n.para){
          V[m1,m2]<-one%*%PI%*%(D.matrix[,m1]*D.matrix[,m2])
          if(m1!=m2)
            {M[m1,m2]<-(one%*%PI%*%(D.matrix[,m1]*D.matrix[,m2]))/
                       (one%*%PI%*%(D.matrix[,m1]*D.matrix[,m1]))}
          }
        }
      Et1<-R-M%*%Et
      mut1<-(one%*%(Y-(PI%*%D.matrix%*%Et1)))/N
      sigt1<-(t(Y-mut1)%*%(Y-mut1)-2*(t(Y-mut1)%*%PI%*%D.matrix%*%Et1)+(t(Et1)%*%V%*%Et1))/N
      MUt1<-c(rep(mut1,3^nQTL))+D.matrix%*%Et1
      Lih10=Lih1<-matrix(0,N,3^nQTL)
      for(i.L in 1:3^nQTL){
        Lih10[,i.L]<-Freq[,i.L]*dnorm(Y,mean=MUt[i.L],sd=sigt^0.5)
        Lih1[,i.L]<-Freq[,i.L]*dnorm(Y,mean=MUt1[i.L],sd=sigt1^0.5)
        }
      L10<-sum(log(rowSums(Lih10)))
      L1<-sum(log(rowSums(Lih1)))
      while(abs(L10-L1)>=0.00001)
        repeat{mut<-mut1
          Et<-Et1
          sigt<-sigt1
          MUt<-MUt1
          L10<-L1
          PIup<-matrix(0,N,3^nQTL)
          for(i.p in 1:3^nQTL){
            PIup[,i.p]<-Freq[,i.p]*dnorm(Y,mean=MUt[i.p],sd=sigt^0.5)
            }
          PI<-PIup/rowSums(PIup)
          for(m1 in 1:n.para){
            R[m1,]<-(t(Y-mut)%*%PI%*%D.matrix[,m1])/(one%*%PI%*%(D.matrix[,m1]*D.matrix[,m1]))
            for(m2 in 1:n.para){
              V[m1,m2]<-one%*%PI%*%(D.matrix[,m1]*D.matrix[,m2])
              if(m1!=m2)
                {M[m1,m2]<-(one%*%PI%*%(D.matrix[,m1]*D.matrix[,m2]))/
                           (one%*%PI%*%(D.matrix[,m1]*D.matrix[,m1]))}
              }
            }
          Et1<-R-M%*%Et
          mut1<-(one%*%(Y-(PI%*%D.matrix%*%Et1)))/N
          sigt1<-(t(Y-mut1)%*%(Y-mut1)-2*(t(Y-mut1)%*%PI%*%D.matrix%*%Et1)+(t(Et1)%*%V%*%Et1))/N
          MUt1<-c(rep(mut1,3^nQTL))+D.matrix%*%Et1
          Lih10=Lih1<-matrix(0,N,3^nQTL)
          for(i.L in 1:3^nQTL){
            Lih10[,i.L]<-Freq[,i.L]*dnorm(Y,mean=MUt[i.L],sd=sigt^0.5)
            Lih1[,i.L]<-Freq[,i.L]*dnorm(Y,mean=MUt1[i.L],sd=sigt1^0.5)
            }
          L10<-sum(log(rowSums(Lih10)))
          L1<-sum(log(rowSums(Lih1)))
          break(abs(L10-L1)<0.00001)}
      LRT<--2*(L0-L1)
      parameter<-c(mut1,Et1,sigt1)}
   if(EM==0)
     {lik0<-function(theta){
        mu0<-theta[1]
        sig0<-theta[2]
        likf0<--sum(log(dnorm(Y,mean=mu0,sd=sig0^0.5)))
        return(likf0)
        }
      L0<-optim(c(mean(Y),var(Y)),lik0)
      h0lrt<-L0$value
      lik1<-function(theta){
        mut<-theta[1]
        sigt<-theta[n.para+2]
        Et<-theta[-(c(1,n.para+2))]
        MUt<-c(rep(mut,3^nQTL))+D.matrix%*%matrix(Et,n.para,1)
        Lihood<-matrix(0,N,3^nQTL)
        for(i.L in 1:3^nQTL){
          Lihood[,i.L]<-Freq[,i.L]*dnorm(Y,mean=MUt[i.L],sd=sigt^0.5)
          }
        likf1<--sum(log(rowSums(Lihood)))
        return(likf1)
        }
      L1<-optim(c(mean(Y),rep(0,n.para),var(Y)),lik1,control=list(maxit=3000))
      parameter<-L1$par
      h1lrt<-L1$value
      LRT<-2*(h0lrt-h1lrt)}
   output<-vector("list",3)
   output[[1]]<-posi
   output[[2]]<-LRT
   output[[3]]<-parameter
   return(output)
   }

QTL.estim.t <- function(nQTL, posi, ys, tL, tR, D.matrix, marker.map, marker.gty, EM=1){
   maxit <- 3000
   tol <- 1e-5
   nS <- length(ys)
   nType <- nrow(D.matrix)
   nE <- ncol(D.matrix)
   Freq <- Mixprop(nQTL, nu=0, posi, marker.map, marker.gty, model=2)
   init <- c(mean(ys), rep(0, nE), sd(ys))
   negLL0 <- function(theta) {
     mu <- theta[1]
     sigma <- theta[2]
     tauL <- (tL - mu)/sigma
     tauR <- (tR - mu)/sigma
     U <- 1 + pnorm(tauL) - pnorm(tauR)
     LL <- sum(log(dnorm(ys, mu, sigma)/U))
     return(-LL)
   }
   temp <- optim(c(mean(ys), sd(ys)), negLL0, method = "Nelder-Mead", control=list(maxit = maxit))
   LL0 <- -1 * (temp$value)
   if (EM==1) {
     err <- 1
     mu0<-mean(ys)
     E<-c(rep(0, nE))
     sigma<-sd(ys)
     mu <- mu0 + matrix(D.matrix %*% E, nS, nType, byrow = TRUE)
     tauL <- (tL - mu0 - D.matrix %*% E)/sigma
     tauR <- (tR - mu0 - D.matrix %*% E)/sigma
     U <- matrix(1 + pnorm(tauL) - pnorm(tauR), nS, nType, byrow = TRUE)
     LL <- sum(log(rowSums(Freq * dnorm(ys, mu, sigma)/U)))
     for (iter in 1:maxit) {
       Pi <- Freq * dnorm(ys, mu, sigma)/U
       Pi <- Pi/rowSums(Pi)
       Amy <- matrix(dnorm(tauL) - dnorm(tauR), nS, nType, byrow = TRUE)/U
       Bob <- matrix(tauL * dnorm(tauL) - tauR * dnorm(tauR), nS, nType, byrow = TRUE)/U
       V <- matrix(0, nE, nE)
       for (i in 1:nE) {
         V[i, ] <- colSums(Pi) %*% (D.matrix[, i] * D.matrix)
       }
       M <- V/diag(V) - diag(1, nE)
       r <- (t(ys - mu0) %*% Pi + sigma * colSums(Pi * Amy)) %*%
            D.matrix/colSums(Pi %*% D.matrix^2)
       E <- t(r) - M %*% E
       mu0 <- mean(ys - Pi %*% D.matrix %*% E + sigma * rowSums(Pi * Amy))
       sigma <- (t(ys - mu0) %*% (ys - mu0) - 2 * t(ys - mu0) %*% Pi %*%
                D.matrix %*% E + t(E) %*% V %*% E)/(nS - sum(Pi * Bob))
       sigma <- sqrt(as.vector(sigma))
       mu <- mu0 + matrix(D.matrix %*% E, nS, nType, byrow = TRUE)
       tauL <- (tL - mu0 - D.matrix %*% E)/sigma
       tauR <- (tR - mu0 - D.matrix %*% E)/sigma
       U <- matrix(1 + pnorm(tauL) - pnorm(tauR), nS, nType, byrow = TRUE)
       LL1 <- sum(log(rowSums(Freq * dnorm(ys, mu, sigma)/U)))
       err <- abs(LL1 - LL)
       LL <- LL1
       if (err < tol) {
         break
       }
     }
     parameter <- c(mu0, E, sigma^2)
   }
   if (EM==0) {
     negLL1 <- function(theta) {
       mu0 <- theta[1]
       E <- theta[1 + 1:nE]
       sigma <- tail(theta, 1)
       mu <- mu0 + matrix(D.matrix %*% E, nS, nType, byrow = TRUE)
       tauL <- (tL - mu0 - D.matrix %*% E)/sigma
       tauR <- (tR - mu0 - D.matrix %*% E)/sigma
       U <- matrix(1 + pnorm(tauL) - pnorm(tauR), nS, nType, byrow = TRUE)
       LL <- sum(log(rowSums(Freq * (dnorm(ys, mu, sigma)/U))))
       return(-LL)
     }
     temp <- optim(init, negLL1, method = "Nelder-Mead", control=list(maxit=maxit, abstol=tol))
     LL1 <- -1 * (temp$value)
     parameter <- temp$par
     parameter[nE + 2] <- parameter[nE + 2]^2
   }
   LRT <- -2*(LL0 - LL1)
   output<-vector("list",3)
   output[[1]]<-posi
   output[[2]]<-LRT
   output[[3]]<-parameter
   return(output)
}
