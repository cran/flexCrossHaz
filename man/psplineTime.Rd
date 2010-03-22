\name{psplineTime}
\alias{psplineTime}
\title{
Penalized splines using a B-basis and difference penalty
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Specifies a penalised spline basis for the predictor with time-varying effect in order to model 
a time-varying hazard ratio in the Cox model framework. It is a modification of 
\code{\link[survival]{pspline}} allowing to model a time-varying hazard ratio within a Cox model framework.
}
\usage{
psplineTime(x, df = 4, theta, nterm = 2.5 * df, by.x, pen.diff = 2,
           degree = 3, eps = 0.1, v, crs.point, Bs = FALSE, 
           method, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x }{
the time-variable.
}
  \item{df }{
	the desired degrees of freedom. One of the arguments \code{df} or \code{theta} must be 
  given, but not both.
}
  \item{theta }{
  roughness penalty for the fit. \code{theta=0} corresponds to an unconstrained fit 
  and \code{theta=1} provides a \code{pen.diff-1} polynomial fit.
}
  \item{nterm }{
 number of intervals in the range of survival times. See \code{\link{bsplineSurv}}. 
}
  \item{by.x }{
 the predictor with time-varying effect (both a factor or a continuous variable).
}
  \item{pen.diff }{
 order of penalty difference, usually \code{pen.diff=1} or \code{2} (default).
}
  \item{degree }{
degree of splines (defaults to 3). 
}
  \item{eps }{
 accuracy for \code{df}.
}
  \item{v }{
 vector of same length of the spline coefficient vector. Negative values mean that the corresponding 
 coefficient is strongly penalized. This argument is useful to enforce positiveness (negativeness) of 
 the time-varying effect. Usually \code{v} are the (minus) coefficients returned by an unconstrained 
 fit \code{coxph} object. See examples. 
}
  \item{crs.point }{
 the known crossing point. If provided the fitted log hazard ratio intersects the time axis at 
 just in \code{crs.point}. See example.
}
  \item{Bs }{
 if TRUE uses a modified B-spline basis as returned by \code{bs(x, df=nterm, degree=degree, intercept=TRUE)},
 otherwise an 'un-modified' B-spline basis where, unlike \code{pspline()} the first column is not modified.
}
  \item{method }{
the method for choosing the tuning parameter \code{theta}. If theta is given, then 'fixed' is 
assumed. If the degrees of freedom is given, then 'df' is assumed. If method='aic' then the degrees of freedom is chosen automatically using Akaike's information criterion.
}
  \item{\dots }{
optional arguments to the control function.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Object of class \code{coxph.penalty} containing the spline basis, with the appropriate 
attributes to be recognized as a penalized term by the coxph or survreg functions.
}
\section{Warning }{Currently \code{summary} and \code{print} methods don't work for
 \code{psplineTime}.}
\references{
Eilers, P.H.C. and Marx, B.D. (1996) Flexible smoothing with B-splines and penalties.
 \emph{Statistical Sciences} \bold{11}, 89--121.
 
 Muggeo, V.M.R. and Tagliavia, M (2010) A flexible approach to the crossing 
hazards problems. \emph{submitted} %\bold{xx}, xx--xx.
}
\author{
 Vito M. R. Muggeo, \email{vito.muggeo@unipa.it} 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[survival]{coxph}},  \code{\link[survival]{pspline}}, \code{\link{bs}}
}
\examples{
\dontrun{
data(gastric)
 # Plot the time-varying effect of treat  
plot(cox.zph(coxph(Surv(time,status) ~ treat, data=gastric))); abline(h=0)  
# Create a new dataset in counting process format,
#  with a start time, stop time, and event status for each record
eventTime<-sort(with(gastric,unique(time[status==1])))
gastric1<-survSplit(gastric, cut=eventTime, end="time", 
          event="status", start="start", id="id")
gastric1$id<-as.integer(gastric1$id)
gastric1<-gastric1[order(gastric1$id),]

#a (quasi) PH fit..
PHfit <- coxph(Surv(start,time,status) ~ psplineTime (start,theta=.9999,
        pen.diff=1, nterm=15, by.x=treat), data=gastric1, x=TRUE) 

fit1 <- coxph(Surv(start,time,status) ~ psplineTime (start,theta=0.7,
        pen.diff=1, nterm=15, by.x=treat), data=gastric1, x=TRUE) 

#a constrained fits.. 
#..fixed crossing points      
fit2 <- coxph(Surv(start,time,status) ~ psplineTime (start,theta=0.7,
         pen.diff=1, nterm=15, by.x=treat, crs.point=200), gastric1, x=TRUE) 
#..positive log hazard ratio
fit3 <- coxph(Surv(start,time,status) ~ psplineTime (start,theta=0.7,pen.diff=1, 
  nterm=15, by.x=treat, v=coef(fit)), gastric1, x=TRUE)

#..negative log hazard ratio
fit4<-coxph(Surv(start,time,status) ~ psplineTime (start,theta=0.7,pen.diff=1, 
  nterm=15, by.x=treat, v=-coef(fit)), gastric1, x=TRUE)

}
}
\keyword{regression}
%%\keyword{ nonparametric}
%%\keyword{ smooth}
%%\keyword{ survival}

