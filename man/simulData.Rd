\name{simulData}
\alias{simulData}
\docType{data}
\title{
Simulated survival data having a step-function time varying log hazard ratio.
}
\description{
Survival data (two groups, each with 100 subjects) generated via a logistic regression 
setting, with log hazard ratio having a step-function equal to \code{-1.5+2.5*I(time>20)}. 
Data are in counting process format.
}
\usage{data(simulData)}
\format{
  A data frame with 6723 rows regarding 200 subjects and the following 5 variables.
  \describe{
    \item{\code{time}}{a numeric vector of survival or censoring time}
    \item{\code{status}}{censoring status}
    \item{\code{group}}{group variable}
    \item{\code{start}}{the starting time for the interval}
    \item{\code{id}}{subject identification number }
  }
}

\examples{
data(simulData)

#standard L2 penalty
fitsim<- coxph(Surv(start,time,status) ~ psplineTime(start, group,
         theta=0.7, pen="L2", nterm=15), data=simulData) 

#better: (approximate) L1 penalty        
fitsim1<-coxph(Surv(start,time,status) ~ psplineTime(start,group,
         theta=0.7, pen="L1", weights=coef(fitsim), nterm=15), 
         data=simulData)

#even better: (approximate) L1 penalty  with updated weights..   
fitsim2<-coxph(Surv(start,time,status) ~ psplineTime(start, group,
         theta=0.7, pen="L1", weights=coef(fitsim1), nterm=15), 
         data=simulData)
        
par(mfrow=c(2,2))
plot(fitsim)
plot(fitsim1)
plot(fitsim2)
}
\keyword{datasets}
