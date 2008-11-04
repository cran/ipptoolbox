\name{bpafromprobs}
\alias{dsminmeanmax}
\alias{dsminmodemax}
\alias{dsmeanvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Obtain BPAs from statistical properties}
\description{
 Various functions to obtain conservative BPA estimates from statistical properties such as mean, variance, ...}
\usage{
dsminmeanmax(intervalnumber, min, mean, max)
dsminmodemax(intervalnumber, min, mode, max)
dsmeanvar(intervalnumber, mean, var)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{intervalnumber}{Number of focal elements of the discretization}
  \item{min}{Minimum of the BPA}
  \item{mean}{Mean of the BPA}
  \item{mode}{Mode of the BPA}
  \item{max}{Maximum of the BPA}
  \item{var}{Variance of the BPA}
}
\details{
\itemize{
\item \code{dsminmeanmax} generates a BPA from min, mean and max.
\item \code{dsminmodemax} generates a BPA from min, mode and max.
\item \code{dsmeanvar} generates a BPA from mean and variance.
}
A set of functions to obtain a BPA that bounds all CDFs with given statistical properties (mean, variance...). These functions are especially useful for quantification of expert estimates. dsminmeanmax generates a BPA from min, mean and max, dsminmodemax from min, mode and max and dsmeanvar from mean and variance. It is necessary to provide an amount of focal elements to be generated (discretization accuracy).

The example generates 3 BPAs. The first bounds all distributions with min=10, mode=60 and max=100. The second all with min=10, mean=30 and max=100. The third all with mean=30 and var=5.
}
\value{
BPA with \code{intervalnumber} focal elements bounding all possible distributions with the given parameters.
}
\author{Philipp Limbourg <p.limbourg@uni-due.de>}
\examples{
pmin=10
pmean=30
pmode=60
pmax=100
pvar=5
y1=dsminmodemax(100,pmin,pmode,pmax)
dscdf(y1);

y2=dsminmeanmax(100,pmin,pmean,pmax)
dscdf(y2);

y3=dsmeanvar(100,pmean,pvar)
dscdf(y3)
}