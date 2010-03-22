\name{plot.coxph}
\alias{plot.coxph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot time-varying effect for a coxph object
}
\description{
Plot the estimated time varying log hazard ratio.
}
\usage{
\method{plot}{coxph}(x, term, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
the Cox model as returned by the \code{coxph()} function. The linear predictor 
has to include time varying smoothed effect via the \code{psplineTime()} function.
}
  \item{term}{
the covariate having a time varying smoothed effect. It can be missing when 
the Cox model includes only the spline coefficients. 
}
  \item{\dots}{
other graphics parameters to pass on to plotting commands.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
None.
}
\references{
Muggeo, V.M.R. and Tagliavia, M (2010) A flexible approach to the crossing 
hazards problems. \emph{submitted} %\bold{xx}, xx--xx.
}
\author{
Miriam Tagliavia
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{cox.zph}}
}
\examples{
\dontrun{
#continues from ?psplineTime
par(mfrow=c(2,2))
plot(fit1)
plot(fit2)
plot(fit3)
plot(fit4)

#compare two plots..
par(mfrow=c(2,1))
plot(fit,term="treat",type="l")

plot(cox.zph(coxph(Surv(time,status) ~ treat, data=gastric),transform="identity"))
abline(h=0)
}
}

\keyword{ smooth }
\keyword{ survival }
