\name{bsplineSurv}
\alias{bsplineSurv}
\title{
B-spline basis for time varying effect
}
\description{
B-spline basis with equally-spaced knots. 
}
\usage{
bsplineSurv(x, nterm, deg = 3, deriv = 0)
}
\arguments{
  \item{x }{
a numeric vector of values at which to evaluate the B-spline functions or
 derivatives. In the crossing hazard context, the time variable. The function 
 does not apply to factor variables.
}
  \item{nterm }{
a positive integer giving the number of intervals in which the range of \code{x}
has to be divided.
}
  \item{deg }{
degree of the piecewise polynomial. Default is 3 for cubic splines.
}
  \item{deriv }{
an integer with values between 0 and deg. The derivative of the given order is
 evaluated at the x positions. Default is 0.
}
}
\value{
A matrix of dimension \code{length(x)*(nterm + deg)}.

\code{bsplineSurv()} is based on the function \code{\link[splines]{spline.des}}. 
It generates a B-spline basis matrix of rank `\code{nterm} + \code{deg}' of degree \code{deg}. The employed knots are equally spaced
outside the observed range of the variable \code{x}. This basis is used by \code{psplineTime} to estimate the time varying log hazard ratio in 
the Cox model.
}
\references{
Eilers, P.H.C. and Marx, B.D. (1996) Flexible smoothing with B-splines and penalties.
 \emph{Statistical Sciences} \bold{11}, 89--121.
}
\author{
 Vito M. R. Muggeo \email{vito.muggeo@unipa.it } 
}

\examples{
x<-1:100/100
base<-bsplineSurv(x,nterm=10)
matplot(x,base,type="l",main="A B-spline basis with equally spaced knots")
}
\keyword{ nonparametric }
\keyword{ smooth }

