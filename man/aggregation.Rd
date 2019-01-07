\name{aggregation}
\alias{dsdempstersrule}
\alias{dsenvelope}
\alias{dsintersect}
\alias{dssaveintersect}
\alias{dsxaveraging}
\alias{dsxenvelope}
\alias{dswavg}
\alias{dswmix}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Aggregation}
\description{
  Various functions to aggregate a set of estimated BPAs to a single aggregation.
}
\usage{
dsdempstersrule(..., maxfocals = 1e+07)
dsenvelope(...)
dsintersect(...)
dssaveintersect(..., method = "averaging")
dsxaveraging(..., w = NULL, maxfocals = 1e+07)
dsxenvelope(..., maxfocals = 1e+07)
dswavg(..., w = NULL)
dswmix(..., w = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{1...n BPAs to be aggregated.}
  \item{w}{Vector of n weights for weighted aggregations}
  \item{maxfocals}{Maximal number of focal elements of the result. Warning: If not constrained, some rules use time and memory \eqn{O(|m|^|x|)}, |m| number of focals by structure, |x| number of structures. Meaning: 5 structures with 1000 focals: Maximal \eqn{10^15} bytes}
  \item{method}{method="averaging": In case of empty intersection, focal elements will be averaged. method="averaging": In case of empty intersection, focal elements will be enveloped.}
}
\details{
A set of functions to aggregate a set of estimated BPAs to a single aggregation. Some functions have an optional maxfocals parameter that limits the size of the resulting BPA (which theoretically is the product of all focal elements). The vector w contains weights in case of a weighted aggregation type. The dssaveintersect method has an optional parameter that specifies the alternative rule in case of empty intersections.

In the example, two p-boxes defined from exponential distributions are aggregated using dsdempstersrule and dswmix.
}
\value{
Aggregation of the BPAs according to the selected rule.
}
\references{Sentz, K. and S. Ferson (2002). Combination of Evidence in Dempster-Shafer Theory. Sandia Report. Albuquerque, USA, Sandia National Laboratories.}
\author{Philipp Limbourg <p.limbourg@uni-due.de>}
\examples{
lambda1=dsstruct(c(2,3,1))
dss1=dsodf('qexp',100,lambda1);
lambda2=dsstruct(c(5,6,1))
dss2=dsodf('qexp',100,lambda2);
y=dsdempstersrule(dss1,dss2);
y2=dswmix(dss1,dss2,w=c(2,0.5));
dscdf(y);

dscdf(y2);
}
