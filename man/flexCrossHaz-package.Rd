\name{flexCrossHaz-package}
\alias{flexCrossHaz-package}
\alias{flexCrossHaz}
\docType{package}
\title{ 
Flexible crossing hazards in the Cox model
}
\description{
Estimation of Cox model with flexible time-varying effects via P-splines and possible crossing points.
 }
\details{
\tabular{ll}{
Package: \tab flexCrossHaz\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2010-03-04\cr
License: \tab GPL\cr
}
Given a covariate with a time-varying effect, the package \code{flexCrossHaz} allows to 
  include it in the `linear predictor' of the Cox model. Penalized splines ar employed 
  to obtain a smooth estimate of the time-varying effect, and the possible crossing point is estimated. 
  Several variables with time-varying effect are allowed.
}
\author{
Vito M.R. Muggeo \email{vito.muggeo@unipa.it} and Miriam Tagliavia \email{tagliavia@dssm.unipa.it}
% Maintainer: Who to complain to <yourfault@somewhere.net>
}
\references{
Verweij, P. and van Houwelingen, H. (1995) Time-dependent effects of fixed 
covariates in Cox regression. \emph{Biometrics} ; \bold{51}:1550--1556.

Mantel, N. and Stablein, D. (1988) The crossing hazard function problem. \emph{The
  Statistician} ; \bold{37}:59--64.
  
Liu, K., Qiu, P. and Sheng, J. (2007) Comparing two crossing hazard rates by Cox
 proportional hazards modelling. \emph{Statistics in Medicine} ; \bold{26}:375--391.  
 
Muggeo, V.M.R. and Tagliavia, M. (2010) A flexible approach to the crossing 
hazards problems. \emph{Submitted}.
}
\keyword{ package }
\keyword{regression}
\seealso{
\code{\link[survival]{coxph}}
}
