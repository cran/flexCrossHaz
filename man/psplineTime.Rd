\name{psplineTime}
\alias{psplineTime}
\title{
Penalized splines using a B-spline basis and difference penalty
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Specifies a B-spline basis of time for the predictor with a time-varying effect within the Cox model. 
It is a modification of \code{\link[survival]{pspline}} in package \code{survival}.
}
\usage{
psplineTime(x, by.x, weights, crs.point, pen = c("L2", "L1"), 
            b.sign = c("unspecified", "positive", "negative"), 
            monot = c("unspecified", "increasing", "decreasing"), 
            df = 4, theta, nterm = 2.5 * df, pen.diff = 1, 
            degree = 3, eps = 0.1, Bs = FALSE, method, ...)  
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
the time variable.
}
  \item{by.x }{
 the predictor with time-varying effect (both a factor or a continuous variable).
}
  \item{weights}{
 possible vector of the same length of the spline coefficient vector to be passed if a constrained
 time-varying effect is sought. It has to be provided when \code{pen="L1"}, and/or \code{b.sign="positive"} (or 
 \code{"negative"}), and/or \code{monot="increasing"} (or \code{"decreasing"}). Usually \code{weights} 
 are the coefficients returned by an unconstrained fit \code{coxph} object. See examples.
}
  \item{crs.point }{
 the known crossing point. If provided the fitted log hazard ratio intersects the time axis at 
 just in \code{crs.point}. See example.
}
  \item{pen}{
the penalty norm to apply to the coefficients differences. \code{pen="L2"} uses a standard quadratic 
  (\eqn{L_2}norm) penalty, while \code{pen="L1"} uses an approximate \eqn{L_1} norm given by \eqn{\sum\Delta\beta_j^2/|w|} where 
  \eqn{\Delta\beta_j}{diff(beta)_j} are the coefficient differences and \eqn{w} are weights, i.e. known coefficients, provided 
  in the \code{weights} argument. \code{pen="L1"} should be preferred when abrupt changes are expected 
  in the log hazard ratio (see examples in \code{?simulData}).
}
  \item{b.sign}{
a character string specifying if the time-varying effect sign is unconstrained (the 
default \code{"unspecified"}), or has to be enforced to be \code{"positive"} or  \code{"negative"}.
}
  \item{monot}{
a character string specifying if monotonicity has to be enforced (\code{"increasing"} or
 \code{"decreasing"}) or not (the default \code{"unspecified"}).
}
  \item{df}{
	the desired degrees of freedom. One of the arguments \code{df} or \code{theta} must be 
  given, but not both. 
}
  \item{theta}{
  roughness penalty for the fit. \code{theta=0} corresponds to an unconstrained fit 
  and \code{theta1} tending to 1 provides a \code{pen.diff-1} polynomial fit.
}
  \item{nterm }{
 number of intervals in the range of survival times. See \code{\link{bsplineSurv}}. 
}
  \item{pen.diff }{
 order of penalty difference, usually \code{1} or \code{2} (defaults to \code{pen.diff=1}).
}
  \item{degree }{
degree of splines (defaults to 3). 
}
  \item{eps }{
 accuracy for \code{df}.
}
  \item{Bs }{
 if TRUE uses a modified B-spline basis as returned by \code{bs(x, df=nterm, degree=degree, intercept=TRUE)},
 otherwise an 'un-modified' B-spline basis according to Eilers and Marx (1996). 
}
  \item{method }{
the method for choosing the tuning parameter \code{theta}. If theta is given, then 'fixed' is 
assumed. If the degrees of freedom is given, then 'df' is assumed. 
If method='aic' then the degrees of freedom is chosen automatically using Akaike's information criterion.
Currently only \code{method='fixed'} works.
}
  \item{\dots }{
optional arguments to the control function.
}
}

\details{
\code{psplineTime} allows straightforward modelling of time varying effect in Cox model when the 
\code{survival::coxph} is employed. The original \code{pspline()} does not allow ready inclusion of 
terms with time-varying effect. (for instance \code{pspline()} modifies the first column of the B-spline basis)
}
\value{
Object of class \code{coxph.penalty} containing the spline basis, with the appropriate 
attributes to be recognized as a penalized term by the coxph or survreg functions.
}
\details{
A single iteration could fail in order to enforce monotonicity or positiveness (negativeness); thus 
two or three iteration are recommended.
}
\section{Warning }{
Currently only \code{method="fixed"} works, i.e. the \code{theta} argument has to be specified.

Moreover \code{summary.coxph} and \code{print.coxph} methods do not 
work for object returned by \code{coxph} using \code{psplineTime}.

}
\references{
Eilers, P.H.C. and Marx, B.D. (1996) Flexible smoothing with B-splines and penalties.
 \emph{Statistical Sciences} \bold{11}, 89--121.
 
 Muggeo, V.M.R. and Tagliavia, M. (2010) A flexible approach to the crossing 
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
PHfit <-coxph(Surv(start,time,status) ~ psplineTime(start,by.x=treat,
        theta=.9999, pen.diff=1, nterm=15), gastric1) 

fit1 <-coxph(Surv(start,time,status) ~ psplineTime(start,by.x=treat, 
       theta=0.7, nterm=15), gastric1)

#a constrained fits.. 
#..fixed crossing points      
fit2 <-coxph(Surv(start,time,status) ~ psplineTime(start,by.x=treat,
       theta=0.7, nterm=15, crs.point=200), gastric1)

#..positive log hazard ratio
fit3 <-coxph(Surv(start,time,status) ~ psplineTime(start,by.x=treat, 
       theta=0.7, nterm=15, weights=coef(fit1), b.sign="p"),
       gastric1)

#..negative log hazard ratio
fit4 <-coxph(Surv(start,time,status) ~ psplineTime(start,by.x=treat,
       theta=0.7, nterm=15, weights=coef(fit1), b.sign="n"),
       gastric1)
}
}
\keyword{regression}
%%\keyword{ nonparametric}
%%\keyword{ smooth}
%%\keyword{ survival}

