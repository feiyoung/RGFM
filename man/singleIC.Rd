\name{singleIC}
\alias{singleIC}
\title{IC(PC) criteria for selecting number}
\description{
  IC(PC) criteria for selecting number of  factors in generalized factor models.
}
\usage{
  singleIC(X, group, type, q_set=1:10, dropout=0, dc_eps=1e-4, maxIter=10, omega=1/ncol(X), output=F, fast_version=T)
}
\arguments{
  \item{X}{a matrix with dimension of n*p(p=(p1+p2+..+p_d)),observational mixed data matrix, d is the types of variables, p_j is the dimension of j-th type variable.}
  \item{group}{a vector with length equal to p, specify each column of X belonging to which group.}
  \item{type}{a d-dimensional character vector, specify the type of variables in each group. For example, \code{type=c'poisson', 'binomial')}, and it is referred to the help file of \link{glm.fit} function for more details.}
  \item{q_set}{a positive integer vector, specify the candidates of factor number q, (optional) default as \code{c(1:10)} according to Bai,2013.}
  \item{q}{a positive integer or empty, specify the number of factors. If q is \code{NULL}, then IC criteria is used to dertemined $q$ automatically.}
  \item{dropout}{a proper subset of $[1, 2, ..., d]$,  specify which group to be dropped in obtaining the initial estimate of factor matrix $H$, and the aim is to ensure the convergence of algorithm leaded by weak signal variable types. Optional parameter with default as 0, no group dropping.}
  \item{dc_eps}{ positive real number, specify the tolerance of varing quantity of objective function in the algorithm. Optional parameter with default as \code{1e-4}.}
  \item{maxIter}{a positive integer, specify the times of iteration. Optional parameter with default as 50.}
  \item{output}{a logical value with TRUE or FALSE, specify whether ouput the mediate information in iteration process, (optional) default as FALSE.}
  \item{fast_version}{logical value with TRUE or FALSE, \code{fast_version = TRUE}: use the fast algorithm which omit the one-step updating, but it cannot ensure the estimation efficieny; \code{fast_version = FALSE}: use the original algorithm; (optional) default as \code{FALSE}; }
}
\details{
  This function also has the MATLAB version at \url{https://github.com/feiyoung/GFM/blob/master/singleIC.m}.
}
\value{
  return an integer, the estimated number of factors.
}
\references{
  Wei Liu, Huazhen Lin, Shurong Zheng & Jin Liu (2019) . Generalized factor model for ultra-high dimensional mixed data. Submitted.

  Bai, J. and Liao, Y. (2013). Statistical inferences using large esti-
    mated covariances for panel data and factor models.
}
\author{
  Liu Wei
}
\note{
  nothing
}
\seealso{
  nothing
}
\examples{
  ## Homogeneous  normal variables
  dat <- gendata(q = 2, n=100, p=100, rho=3)
  group <- rep(1,ncol(dat$X))
  type <- 'gaussian'
  # select q automatically
  singleIC(dat$X, group, type, q_set = 1:6, output = F)
}

\keyword{singleIC}
