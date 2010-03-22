plot.coxph <- function (x, term, ...){
  # x object returned by coxph
      if(missing(term)) {
          id<-1:length(coef(x))
            } else {
          nomeParent<-paste("\\(", term, "\\)",sep="")
          id<-grep(nomeParent,names(coef(x)))
          }
      #strsplit(strsplit(names(obj$pterms),", by.x = ")[[1]][2],"\\)")[[1]]
      b<-coef(x)[id]
      #b.data<- x$x[term==1,]
      nomeTempo<-strsplit(strsplit(names(x$pterms),c("psplineTime\\("))[[1]][2],",")[[1]][1]
      #se è stato chiamato come psplineTime(.., x=,..) non funziona allora
      evTime<-range(as.matrix(x$y)[,nomeTempo])
      deg<-eval(parse(text=strsplit(strsplit(names(x$pterms),c("deg ="))[[1]][2],"\\,")[[1]][1]))
      if(is.na(deg)) deg<-3
      nterm<-length(b)-deg
      xx<-seq(evTime[1],evTime[2],l=50)
      B<-bsplineSurv(xx, nterm=nterm, deg=deg)        
      y<- drop(B%*%b)
      se<-sqrt(diag(B%*%x$var[id,id]%*%t(B)))
      M<-cbind(y-1.96*se,y,y+1.96*se)
      ylab=if(missing(term)) "Log hazard ratio" else paste("Log hazard ratio for",term, "")
      matplot(xx, M, type="l", lty=c(2,1,2), lwd=2, xlab="Time",ylab=ylab, col=1,...)
      abline(h=0,lty=3)
 }
 