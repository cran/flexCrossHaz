\name{gastric}
\alias{gastric}
\docType{data}
\title{
Gastric Cancer Survival Data
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
Survival times in a clinical trial comparing two treatments of locally
non-resectable gastric carcinoma.
}
\usage{data(gastric)}
\format{
  A data frame with 90 observations on the following 3 variables.
  \describe{
    \item{\code{time}}{survival or censoring time}
    \item{\code{status}}{censoring status: 0 = alive, 1 = dead }
    \item{\code{treat}}{treatment variable: 0 = combination of chemotherapy and radiation, 1= only chemotherapy}
  }
}

\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
Stablein, D.M., Carter, W.H.,Jr and Novak, J.W. (1981) Analysis of survival data with 
nonproportional hazard function. \emph{Controlled Clinical Trials} \bold{2}, 149--159.
}
\references{
%%  ~~ possibly secondary sources and usages ~~
Orbe, J.,Ferreira, E. and Nunez-Anton, V. (2002) Comparing proportional 
hazards and accelerated failure time models for survival analysis.
 \emph{Statistics in Medicine} \bold{21}, 3493--3510.
}
\examples{
data(gastric)
%## maybe str(gastric) ; plot(gastric) ...
}
\keyword{datasets}
