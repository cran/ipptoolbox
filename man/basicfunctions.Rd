\name{basicfunctions}
\alias{dsstruct}
\alias{dsnorm}
\alias{dsred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Basic functions}
\description{
Helpful functions for creating and handling BPAs
}
\usage{
dsstruct(x)
dsnorm(y)
dsred(y, thres = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Matrix/BPA}
  \item{y}{BPA}
  \item{thres}{Minimal mass value of a focal element desired}
}
\details{
\itemize{
  \item \code{dsstruct} creates a new BPA from a matrix and normalizes it.
  \item \code{dsnorm} normalizes a BPA (i. e. masses sum up to 1, no lower bound larger than higher bound).
  \item \code{dsred} reduces the number of focal elements in a BPA by merging adjoint elements.
}
The IPP Toolbox represents BPAs in a discretized form. A BPA is stored as a matrix (3 columns) that lists the focal elements. The first and second columns contain lower and upper bounds, the third column the masses of the focal elements. The function dsstruct is a simple function to create a valid BPA. It takes a three-column matrix of the above type as an input. It performs some checks (upper bound >= lower bound, masses sum up to 1) and returns a valid BPA (possibly switching lower & upper bounds of focal elements, normalizing masses to 1) with the according warning messages.
}
\value{
}
\references{Ferson, S., V. Kreinovich, et al. (2003). Constructing Probability Boxes and Dempster-Shafer Structures. Sandia Report. Albuquerque, USA, Sandia National Laboratories.}
\author{Philipp Limbourg <p.limbourg@uni-due.de>}
\examples{
a=dsstruct(matrix(c(2,3,0.1,1.5,5,0.9),ncol=3,byrow=TRUE))
b=dsstruct(matrix(c(2,3,0.5,10,5,0.9),ncol=3,byrow=TRUE))
c=dsodf('qexp',10000,dsstruct(c(10,20,1)));
d=dsred(c,0.02)
dscdf(d);
}