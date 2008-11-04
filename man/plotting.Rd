\name{plotting}
\alias{dscdf}
\alias{dsecdf}
\alias{dsqqplot}
\alias{dsbel}
\alias{dspl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plotting BPAs}
\description{
Various routines to visualize BPAs
}
\usage{
dscdf(x, xrange = NULL, col = c(3, 4), ..., newplot = TRUE)
dsqqplot(ds, sample, points = FALSE)
dsbel(x)
dspl(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{BPA to plot}
  \item{xrange}{Optional: xrange=c(1,5) limits the x-axis to [1,5] }
  \item{col}{Optional: col = c(3,4) defines the colors of the Bel and the Pl line.}
  \item{\dots}{Optional: additional parameters passed to the plot routine (except linewidth, linetype and color)}
  \item{newplot}{Optional: Create a new plot or draw lines onto an existing plot.}
  \item{sample}{set of data values, e. g. c(1.5,2,3.5)}
  \item{ds}{BPA for qq plot resembling sample}
  \item{points}{Optional, points=TRUE: plot the sample values as points + line, otherwise only a straight line.}
}
\details{
\itemize{
\item \code{dscdf} plots Bel([-Inf,x]) and Pl([-Inf,x]). These functions bound all CDFs.
\item \code{dsqqplot} plots a quantile-quantile plot between data and a BPA.
\item \code{dsbel} returns Bel([-Inf,x]), e. g. for own plotting routines.
\item \code{dspl} returns Pl([-Inf,x]), e. g. for own plotting routines.
}
The function dscdf plots Bel([-Inf,x]) and Pl([-Inf,x]), the bounds on all CDFs. Arbitrary graphics parameters can be passed to dscdf. dsqqplot plots a quantile-quantile plot between data and a BPA. dsbel and dspl return the function Bel([-Inf,x])/ Pl([-Inf,x]) as a set of (x,y)-type points, e. g. for own plotting routines.

The example generates a BPA. It plots the BPA first without, then with additional parameters. Finally it uses the functions dsbel and dspl for plotting points into the diagram. Then it generates a BPA from a normal pbox with imprecise standard deviation. It tests with a QQ plot, if the function fits well to a random sample from N(0,1).
}
\value{
dscdf returns the x and y coordinates of the plotted lines for further use.
}
\author{Philipp Limbourg <p.limbourg@uni-due.de>}
\examples{
print("Plot new BPA")
a=dsstruct(rbind(c(1,2,0.3),c(1.5,5,0.5),c(4,6,0.2)))
dscdf(a)

dscdf(a,xrange=c(0,7),xlab="Important parameter")
print("Bel([-Inf,x])")
bel=dsbel(a)
print("Pl([-Inf,x])")
pl=dspl(a)
points(pl,lwd=5,col='red')
points(bel,lwd=5,col='pink')

print("See if a pbox from N(0,[0.5,1.5]) fits data from N(0,1)")
data=rnorm(100,0,1)
pbox=dsodf(qnorm,1000,0,dsstruct(c(0.5,1.5,1)))
print("Plot pbox")
dscdf(pbox)

print("qq plot")
dsqqplot(pbox,data)
}
