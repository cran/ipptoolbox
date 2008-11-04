\name{bpafromdist}
\alias{dskstest}
\alias{dsbelpltests}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Construct pboxes}
\description{
Kolmogorov-Smirnov-Fits for BPAs.
}
\usage{
dskstest(val, ds, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{val}{Data for K-S / C-M test}
  \item{ds}{BPA or handle of a pbox}
  \item{\dots}{Optional pbox parameters (points, BPAs)}
}
\details{
\item{dskstest}{Evaluates the fit of a pbox using a Kolmogorov-Smirnov test}

The method dskstest evaluates the fit of a pbox using a Kolmogorov-Smirnov test. The function works both on point values and interval data. The syntax is close to the R syntax for a probabilistic K-S test ks.test(). Both a BPA and a parametric pbox (in form of a CDF and a set of precise/imprecise parameters) can be tested.

The example generates a random set of points and a random set of intervals from a normal distribution (mean=0.3, sd=1.2). It tests, if a pbox with mean=[-0.5,0.5] and sd=1 fits the point data. In a second example, the user passes a normal CDF with parameters mean=0 and sd=[1,2] to the test. dskstest checks, if it fits the interval data.

}
\value{
BPA representing a pbox sampled from fhandle / K-S test results.
}
\author{Philipp Limbourg <p.limbourg@uni-due.de>}
\examples{
data=rnorm(100,0.3,1.2)
intervaldata=cbind(data,data+runif(100,0,0.2))

mu=0; sigma=1
mu2=dsstruct(c(-0.5,0.5,1)); sigma2=dsstruct(c(1,2,1))
pbox=dsodf(qnorm,1000,mu2,sigma)

print("K-S result, pbox on data:")
ks=dskstest(data,pbox)
print(ks)
print("K-S result, precise mu, imprecise sigma on interval data:")
ks2=dskstest(intervaldata,pnorm,list(mu,sigma2))
print(ks2)
}