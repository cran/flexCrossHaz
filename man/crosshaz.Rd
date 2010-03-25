\name{crosshaz}
\alias{crosshaz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Crossing-point in the smoothed time varying log hazard ratio
}
\description{
Estimates the crossing point through the standard Newton-Raphson method 
applied to the flexible estimate of the time varying log hazard ratio. 
The crossing-point estimate is returned along with its approximate standard error  
computed via the Delta method.
}
\usage{
crosshaz(obj, starting, term, it.max = 20, eps = 1e-05, 
    display = FALSE, k = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj }{
the Cox model as returned by the \code{coxph()} function with the linear predictor including 
a time varying smoothed effect via the \code{psplineTime()}.
}
  \item{starting }{
starting value for the crossing point (within the range of the survival times).
}
\item{term }{
the covariate having a time varying smoothed effect. It can be omitted when the Cox model includes 
only the spline coefficients of the fitted time-varying hazard ratio. 
}
  \item{it.max }{
integer giving the maximum number of iterations allowed in the Newton-Raphson algorithm.
}
  \item{eps }{
the threshold of tollerance before algorithm converges. Default is 1e-05.
}
  \item{display }{
logical indicating if the updated estimates of the crossing point should be printed at each iteration.
}
  \item{k }{
positive factor (less than 1) modifying the increments in crossing-point estimation.
}
}
\details{
The crossing-point estimate is attained via a standard Newton iterative method which searches for a 
root of the estimated smoothed log hazard ratio.
% If the crossing-point
% \eqn{gamma} exists, \eqn{beta(t)} intersect the time-axis just at \eqn{gamma}, 
% thus it is estimated solving for \eqn{gamma} the equation  \eqn{\hat{\beta(t)}=0}.
% Since \eqn{\hat{\beta(t)}=\sum_j b_j*B_j(t)} (where \eqn{B_j(t)} is the B-spline basis 
%constructed on time-variable) is smooth, the root of the above equation may be found 
% by the standard Newton-Raphson method, whose step is
% \deqn{\gamma_0-k*\frac{\hat{beta(gamma_0)}}{\hat{beta'(gamma_0)}}}
% where \eqn{gamma_0} is an appropriate solution (\code{starting}) and 
% \eqn{beta'(gamma_0)} is the first derivative of \eqn{beta(t)} evaluated in \eqn{gamma_0}.
At each iteration the updated estimate is usually given by \code{new.est<-old.est+increm}, 
where \code{increm} is given by the ratio of the objective function and its first derivative 
evaluated in \code{old.est}.

Sometimes to stabilize the procedure it can be useful to set \code{k<1} to 
modify the increments in the crossing-point estimation, thus the updated estimate 
is given by \code{new.est<-old.est+k*increm}.

The number of iterations is governed by \code{it.max} and each updated estimates 
is returned if \code{display= TRUE}.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{est }{the estimated crossing-point.}
  \item{starting }{initial value provided as argument.}
  \item{alpha0 }{the objective function evaluated in the next-to-last updated estimate.}
  \item{alpha1 }{the first derivative of the objective function evaluated in the next-to-last updated estimate.}
%%  \item{v00 }{variance of the objective function evaluated in the next-to-last updated estimate.}
%%  \item{v11 }{variance of the first derivative of the objective function evaluated in the next-to-last updated estimate.}
%%  \item{v01 }{covariance of the objective function and its first derivative evaluated in the next-to-last updated estimate.}   
  \item{s.e.}{crossing point standard error computed via the Delta method.}
  \item{rangeTime }{the range of survival times.} 
  \item{pterms}{the covariate with time varying smoothed effect.}
  \item{nterm}{the number of intervals in which survival times were divided. See \code{bsplineSurv}.} 
 
  %\item{B}{the B-spline basis evaluated in 40 points within the range of the survival times }     
}
\references{
Muggeo, V.M.R. and Tagliavia, M. (2010) A flexible approach to the crossing 
hazards problems. \emph{submitted} %\bold{xx}, xx--xx.
}
\author{
Vito M.R. Muggeo and Miriam Tagliavia 
}

\seealso{
\code{\link[survival]{coxph}}, \code{\link{psplineTime}}
}
\examples{
\dontrun{
%	or do  help(data=index)  for the standard data sets.
data(gastric)
%# Test for a time varying effect for variable treat
%cox.zph(coxph(Surv(time,status) ~ treat, data=gastric))

# Create a new dataset in counting process format,
eventTime<-sort(with(gastric,unique(time[status==1])))
gastric1<-survSplit(gastric, cut=eventTime, end="time", 
          event="status", start="start", id="id")
gastric1$id<-as.integer(gastric1$id)
gastric1<-gastric1[order(gastric1$id),]
 
# Fit a Cox model with a smooth effect for treat..
Fit1<-coxph(Surv(start,time,status) ~ psplineTime(start, by.x=treat,
      theta=0.7, nterm=15), data=gastric1) 

# .. and compute the crossing point
cross.p<-crosshaz(Fit1,starting=300) 
}
}
\keyword{ smooth }
