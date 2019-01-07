\name{bpafromdata}
\alias{dsecdffit}
\alias{dslapconf}
\alias{dsksconf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generating BPAs from point / interval data}
\description{
A set of functions to create empirical BPAs from data sets (either points or intervals).
}
\usage{
dsecdffit(x)
dslapconf(x, lims = c(-Inf, Inf))
dsksconf(x, conf = 0.95, lims = c(-Inf, Inf))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Array of points / matrix with interval data}
  \item{lims}{Optional: Limits of the BPA for cutting first and last focal elements.}
  \item{conf}{Confidence level of Kolmogorov-Smirnov fit in ]0,1[.}
}
\details{
\describe{
\item{dsecdffit}{Creates a BPA of the data set analoguous to probabilistic empirical CDF}
\item{dslapconf}{Creates a BPA of the data set according to the Laplace method}
\item{dsksconf}{Creates a BPA of the data set from the Kolmogorov-Smirnov bounds for a given confidence.}
}
The toolbox contains a set of functions to create empirical BPAs from data sets. These data sets may be either standard data values or intervals (e.g. caused by measurement imprecision). dsecdffit creates a BPA of the data set analogous to probabilistic empirical CDF. It is not conservative respective to distribution tails. dslapconf creates a BPA according to the Laplace method, dsksconf from the Kolmogorov-Smirnov bounds for a given confidence.

The example illustrates the three different methods on interval and point data. First, 20 random values and intervals are generated. Then, dsecdffit, dslapconf and dsksconf are used to generate fits for both data sets.
}
\value{
BPA according to the chosen method.
}
\references{Ferson, S., V. Kreinovich, et al. (2003). Constructing Probability Boxes and Dempster-Shafer Structures. Sandia Report. Albuquerque, USA, Sandia National Laboratories.}
\author{Philipp Limbourg <p.limbourg@uni-due.de>}
\examples{
print("Create random data sets")
setofpoints=rnorm(20,0,1)
setofintervals=cbind(setofpoints-runif(20,0,1),setofpoints+runif(20,0,1))
print("Plot ecdf for both points and interval data")
ecdfp=dsecdffit(setofpoints)
dscdf(ecdfp)

ecdfi=dsecdffit(setofintervals)
dscdf(ecdfi)

print("Plot Laplace BPA for both points and interval data")
lapp=dslapconf(setofpoints,lim=c(-5,5))
dscdf(lapp)

lapi=dslapconf(setofintervals,lim=c(-5,5))
dscdf(lapi)

print("Plot Kolmogorov-Smirnov 75 percent bound BPA for both points and interval data")

ksp=dsksconf(setofpoints,conf=0.75,lim=c(-5,5))
dscdf(ksp)

ksi=dsksconf(setofintervals,conf=0.75,lim=c(-5,5))
dscdf(ksi)
}
