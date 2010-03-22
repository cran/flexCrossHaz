crosshaz <-function(obj, starting, term, it.max=20, eps=.00001 , display=FALSE, k=1){
      if(missing(term)) {
          id<-1:length(coef(obj))
            } else {
          nomeParent<-paste("\\(", term, "\\)",sep="")
          id<-grep(nomeParent,names(coef(obj)))
          }
      #strsplit(strsplit(names(obj$pterms),", by.x = ")[[1]][2],"\\)")[[1]]
      b<-coef(obj)[id]
      if(any(is.na(b))) stop("NA in coeffs .. incorrect `term'?")
      var.b<-obj$var[id,id]
      nomeTempo<-strsplit(strsplit(names(obj$pterms),c("psplineTime\\("))[[1]][2],",")[[1]][1]
      #se è stato chiamato come psplineTime(.., x=,..) non funziona allora
      evTime<-range(as.matrix(obj$y)[,nomeTempo])
      deg<-eval(parse(text=strsplit(strsplit(names(obj$pterms),c("deg ="))[[1]][2],"\\,")[[1]][1]))
      if(is.na(deg)) deg<-3
      nterm<-length(b)-deg
  if(display) message("starting value ", "\t", formatC(round(starting,2)))
  it.max0<-it.max
  toll<-10
  i<-1
  gamma0<-starting
  k<-min(1,abs(k))
    while(toll>eps){
      B.x0<-bsplineSurv(c(starting,range(evTime)),nterm=nterm, deg=deg)[1,]
      f<-drop(B.x0%*%b)
      B.x0.prime<-bsplineSurv(c(starting,range(evTime)),nterm=nterm,deriv=1)[1,]
      f.prime<-drop(B.x0.prime%*%b)
      old<-starting
      starting<-old-k*f/f.prime
      if(starting<min(evTime)||starting>max(evTime)) stop("Invalid estimate!")
      toll<-abs((old-starting)/starting)
      #if(display) message("iteration ", i, "\t", formatC(round(starting,2)))
      if(display) message("iteration ", i, "\t", formatC(starting,format="f",digits = 3))      
      if(i>=it.max) stop("max iterations attained")
      i<-i+1
      }
      v00<-crossprod(B.x0, var.b)%*%B.x0
      v11<-crossprod(B.x0.prime, var.b)%*%B.x0.prime
      v01<-crossprod(B.x0, var.b)%*%B.x0.prime
      # SE delta
      num<-old*f.prime-f 
      den<- f.prime
      var.den<- v11
      var.num<- v11*old^2 - 2*old*v01 + v00
      cov.num.den<- -v01+ old*v11
      se.num<-sqrt(var.num)
      se.den<-sqrt(var.den)
      t1<-num/se.num
      t2<-den/se.den
      r<-cov.num.den/(se.num*se.den)
      rho<-num/den
      s.e.rho<-sqrt((var.num+(rho^2)*var.den-2*rho*cov.num.den)/(den^2))
      ###  
      # B<-bsplineSurv(seq(range(evTime)[1],range(evTime)[2],l=40), nterm=nterm)
      r<-list(est=starting, starting=gamma0, alpha0=f, alpha1=f.prime,
          s.e.=s.e.rho, rangeTime=range(evTime), nterm=nterm, pterms=obj$pterms)
      #v00=v00, v11=v11, v01=v01, 
      class(r)<-"flexCrossHaz"
      r
  }
