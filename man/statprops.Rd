\name{statprops}
\alias{dsaggwidth}
\alias{dsconf}
\alias{dsexpect}
\alias{dsvariance}
\alias{dsbelpl}
\alias{dssummary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Statistical properties of a BPA}
\description{
  Functions for various different statistical properties of a BPA.}
\usage{
dsaggwidth(x)
dsconf(x, conf, confconf = NULL)
dsexpect(x)
dsvariance(x)
dsbelpl(x, a)
dssummary(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{BPA}
  \item{conf}{Confidence level in ]0,1[}
  \item{confconf}{Wilk's confidence bounds in ]0,1[ on the conf (if empty, ommitted).}
  \item{a}{Interval for Bel(x in a), Pl(x in a)}
}
\details{
Various different statistical properties of a BPA: Aggregated width (dsaggwidth), expected value (dsexpect), quantiles (dsconf), variance (dsvariance) and Belief/Plausibility (dsbelpl). The function dsconf is also able to return Wilks' bounds for a given confidence level. The function dssummary is the analogue to a summary statistics and computes a lot of statistics at once.

The example calculates a set of statistical properties of a pbox x.
}
\value{
Statistical property of a BPA (interval). dssummary returns a list of statistics.
}
\references{Kreinovich, V., G. Xiang, et al. (2006). "Computing mean and variance under Dempster-Shafer uncertainty: Towards faster algorithms." International Journal of Approximate Reasoning 42(3): 212-227.}
\author{Philipp Limbourg <p.limbourg@uni-due.de>}
\examples{
mu=dsstruct(c(10,12,1))
sigma=dsstruct(c(1,1.5,1))
x=dsadf('qnorm',1000,mu,sigma)
print("Sample of 1000 focal elements from a normal dist")
print("Mean:")
print(dsexpect(x))
print("Variance:")
print(dsvariance(x))
print("Median:")
print(dsconf(x,0.5))
print("Bel and Pl of x in [4,8]:")
print(dsbelpl(x,c(4,8)))
print("Aggregated width:")
print(dsaggwidth(x))
print("95 percent conf. level with 95 percent Wilks bounds")
print(dsconf(x,0.95,0.95))
print("Summary statistics")
print(dssummary(x))
}
