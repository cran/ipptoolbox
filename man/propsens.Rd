\name{propsens}
\alias{dsevalmc}
\alias{dssensitivity}
\alias{dsbound}
\alias{dsmonotonous}
\alias{dsmonotonous2}
\alias{dsopt}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Monte Carlo propagation and sensitivity analysis}
\description{
Routines for propagating BPAs and sensitivity analysis.
}
\usage{
dsevalmc(fhandle, x, mcIT, optimizer = dsbound, corr = NULL,
samples = NULL, fnoptions = NULL)

dssensitivity(x, parnums, fhandle, uncfn, mcIT, pinch_samples, 
pinch_type = NULL,optimizer = dsmonotonous, fnoptions = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fhandle}{Function handle of type f(x,...), where x is the vector of uncertain variables.}
  \item{x}{List of BPAs representing each variable}
  \item{mcIT}{Number of Monte Carlo iterations}
  \item{optimizer}{Optimizer for solving interval problems. Currently implemented: "dsbound" for monotonously increasing functions, "dsmonotonous" for monotonous functions (unknown direction) and "dsopt" for local optimization using L-BFGS.}
  \item{corr}{Optional: In case of dependent inputs: Vector of correlations}
  \item{samples}{Currently known samples of the function}
  \item{fnoptions}{options to pass to fhandle}

  \item{parnums}{Array of parameter numbers for which a sensitivity analysis should be performed.}
  \item{uncfn}{Uncertainty function handle (e.g. dsaggwidth). Needs to return a single numeric value for a BPA.}
  \item{pinch_samples}{Number of different pinches of an input.}
  \item{pinch_type}{Optional: Pinch a 'distribution' (default), 'interval' or 'point'.}
}
\details{
The Monte Carlo propagation engine is the core of the toolbox. It allows propagating BPAs through arbitrary functions. The user can define the amount of Monte Carlo samples and the optimization method. Currently, four methods are implemented (dsbound, dsmonotonous, dsmonotonous2, dsopt).
\describe{
\item{dsmonotonous}{Recommended for monotonous functions (increasing/decreasing/mixed).}
\item{dsmonotonous2}{Evaluates all corners of the focal element and is recommended, if the function is not monotonous but quite regular.}
\item{dsopt}{Local optimizer for nonmonotonous functions.}
\item{dsbound}{Applicable in the special case that the function is monotonously increasing.}
}
Dependency between focal elements is modelled by a Gaussian copula as proposed in (Ferson, Hajagos et al. 2004). The copula parameters for n inputs need to be given as a vector of size n*(n-1)/2 with the first component c1 in [0,1] representing the dependency between input 1&2, c2 between 1&3, ..., cn-1 between 1&n, cn between 2&3 and so on. This format is analogous to the format used in the R package "copula". The function returns a data structure with a result BPA and all function evaluations.
The method dssensitivity allows sensitivity analysis for a custom uncertainty measure. The function receives a list of BPAs x and a vector parnums giving all indices for which the sensitivity shall be computed. fhandle is the system function, uncfn the uncertainty measure. An arbitrary uncertainty measure can be used as long as it returns a single value for a BPA. A useful choice is "dsaggwidth", the aggregated width measure. The input BPAs are "pinch sample" times reduced to a structure with less uncertainty defined by "pinchtype" (either a point, a distribution or an interval). Each time, the BPA is propagated using mcIT Monte Carlo samples. The function returns a set of sensitivity indices.

The example illustrates the propagation process. Two functions are defined, f1 being monotonous and f2 nonmonotonous. The inputs dss1 and dss2 are defined as pboxes. The correlation is set to 0.5. Then, dss1 and dss2 are propagated through f1 under independence and dependence. In the next step, dsevalmc uses the optimizer dsopt to propagate x through f2. The method dssensitivity is called with the aggregated width measure to calculate AW sensitivity.
}
\value{
  \item{dsevalmc}{Returns list with two elements: [[1]] contains the result, [[2]] all calculated function points.}
}
\references{Ferson, S., J. Hajagos, et al. (2004). Dependence in Dempster-Shafer theory and probability bounds analysis. Sandia Report. Albuquerque, USA, Sandia National Laboratories.

Limbourg, P., R. Savic, et al. (2007). Fault Tree Analysis in an Early Design Stage using the Dempster-Shafer Theory of Evidence. European Conference on Safety and Reliability - ESREL 2007, Stavanger, Norway, Taylor and Francis.}
\author{Philipp Limbourg <p.limbourg@uni-due.de>}
\examples{
print("Define monotonous function f1 and nonmonotonous function f2");
f1=(function(x){ z <- sqrt(abs(x[,1]+ 2*x[,2]));});
f2=(function(x){ z <- sin(x[,1]+ x[,2]);});
print("Create input BPAs x");
mu1=dsstruct(c(2,3,1));
mu2=dsstruct(c(4,4,1));
dss1=dsodf('qnorm',200,mu1);
dss2=dsodf('qnorm',500,mu2);
x=list(dss1,dss2)
correlations=0.5;
print("Propagate through f1 under independence assumption / correlation 0.5");
y_f1_independent=dsevalmc(f1,x,1000)
y_f1_dependent=dsevalmc(f1,x,1000,dsbound,correlations)
print("Propagate through f2 using optimization");

y_f2=dsevalmc(f2,x,1000,dsopt)
print("Plot results:");
dscdf(y_f1_independent[[1]])

dscdf(y_f1_dependent[[1]])

dscdf(y_f2[[1]])

print("Sensitivity on AW measure, inputs 1 & 2");
s=dssensitivity(x,c(1,2),f1,dsaggwidth,20,100,'distribution');
print(s)
}
