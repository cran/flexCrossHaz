\name{plot.coxph}
\alias{plot.coxph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of time-varying effect for a coxph object
}
\description{
Plot the estimated time varying log hazard ratio along with pointwise confidence intervals.
}
\usage{
\method{plot}{coxph}(x, term, add=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
the Cox model as returned by the \code{coxph()} function. The linear predictor 
has to include time-varying smooth effect via the \code{psplineTime()} function.
}
\item{term}{
the covariate having a time varying smooth effect. It can be missing when 
the Cox model includes only the spline coefficients. 
}
\item{add}{ logical, if \code{TRUE}, it adds the plot to current device (only the point estimate).
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
Muggeo, V.M.R. and Tagliavia, M. (2010) A flexible approach to the crossing 
hazards problems. \emph{submitted} %\bold{xx}, xx--xx.
}
\author{
Vito M.R. Muggeo and Miriam Tagliavia 
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
plot(fit1,term="treat")

plot(cox.zph(coxph(Surv(time,status) ~ treat, data=gastric),
    transform="identity"))
abline(h=0)  

#add the fitted line to the current plot..
plot(fit1)
plot(fit4,add=T,col="blue")
}
}



\keyword{ smooth }
\keyword{ survival }
