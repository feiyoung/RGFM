\name{gfm}
\alias{gfm}
\title{Generalized Factor Model}
\description{
  This function is used to conduct the Generalized Factor Model.
}
\usage{
  gfm(X, group, type, q=NULL, parallel=T,dropout=0,
   dc_eps=1e-4, maxIter=50, q_set=1:10, output=T, fast_version=F)
}
\arguments{
  \item{X}{a matrix with dimension of n*p(p=(p1+p2+..+p_d)),observational mixed data matrix, d is the types of variables, p_j is the dimension of j-th type variable.}
  \item{group}{a vector with length equal to p, specify each column of X belonging to which group.}
  \item{type}{a d-dimensional character vector, specify the type of variables in each group. For example, \code{type=c'poisson', 'binomial')}, and it is referred to the help file of \link{glm.fit} function for more details.}
  \item{q}{a positive integer or empty, specify the number of factors. If q is \code{NULL}, then IC criteria is used to dertemined $q$ automatically.}
  \item{parallel}{a logical value with TRUE or FALSE, indicates wheter to use prallel computating. Optional parameter with default as TRUE.}
  \item{dropout}{a proper subset of $[1, 2, ..., d]$,  specify which group to be dropped in obtaining the initial estimate of factor matrix $H$, and the aim is to ensure the convergence of algorithm leaded by weak signal variable types. Optional parameter with default as 0, no group dropping.}
  \item{dc_eps}{ positive real number, specify the tolerance of varing quantity of objective function in the algorithm. Optional parameter with default as \code{1e-4}.}
  \item{maxIter}{a positive integer, specify the times of iteration. Optional parameter with default as 50.}
  \item{q_set}{a positive integer vector, specify the candidates of factor number q, (optional) default as \code{c(1:10)} according to Bai,2013.}
  \item{output}{a logical value with TRUE or FALSE, specify whether ouput the mediate information in iteration process, (optional) default as FALSE.}
  \item{fast_version}{logical value with TRUE or FALSE, \code{fast_version = TRUE}: use the fast algorithm which omit the one-step updating, but it cannot ensure the estimation efficieny; \code{fast_version = FALSE}: use the original algorithm; (optional) default as \code{FALSE}; }
}
\details{
  This function also has the MATLAB version at \url{https://github.com/feiyoung/GFM/blob/master/gfm.m}.
}
\value{
  return a list with class name 'gfm' and including following components,
  \item{hH}{a n*q matrix, the estimated factor matrix.}
  \item{hB}{a p*q matrix, the estimated loading matrix.}
  \item{hmu}{a p-dimensional vector, the estimated intercept terms.}
  \item{obj}{a real number, the value of objective function when the convergence achieves.}
  \item{q}{an integer, the used or estimated factor number.}
  \item{history}{a list including the following 7 components: (1)dB: the varied quantity of B in each iteration; (2)dH: the varied quantity of H in each iteration; (3)dc: the varied quantity of the objective function in each iteration; (4)c: the objective value in each iteration; (5) realIter: the real iterations to converge; (6)maxIter: the tolerance of maximum iterations; (7)elapsedTime: the elapsed time.}
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
  dropout <- 0
  eps2 <- 1e-3
  maxIter <- 30
  # specify q=2
  gfm1 <- gfm(dat$X, group, type, q=2, output = F)
  # select q automatically
  gfm2 <- gfm(dat$X, group, type, q=NULL, q_set = 1:6, output = F)
  # measure the performance of GFM estimators
  measurefun(gfm2$hH, dat$H0, type='ccor')
  measurefun(gfm2$hB, dat$B0, type='ccor')
## Heterogeous normal variables
  dat <- gendata(seed=1, n=100, p=100, type='heternorm', q=2, rho=4)
  group <- rep(1,ncol(dat$X))
  type <- 'gaussian'
  gfm3 <- gfm(dat$X, group, type, q=NULL, q_set = 1:4, output = F)
## Poisson variable
  q <- 3; p <- 200
  dat <- gendata(seed=1, n=200, p=p, type='pois', q=q, rho=4)
  group <- rep(1,ncol(dat$X))
  type <- 'poisson'
  system.time(
    gfm2 <- gfm(dat$X, group, type,parallel = F, q=NULL, q_set = 1:6, output = T, fast_version = T)
  )
  system.time(
    gfm2 <- gfm(dat$X, group, type,parallel = T, q=NULL, q_set = 1:6, output = T, fast_version = T)
  )
  measurefun(gfm2$hH, dat$H0, type='ccor')
  measurefun(gfm2$hB, dat$B0, type='ccor')

## mix of normal and Poisson
  dat <- gendata(seed=1, n=200, p=200, type='norm_pois', q=2, rho=2)
  group <- c(rep(1,ncol(dat$X)/2), rep(2,ncol(dat$X)/2))
  type <- c('gaussian','poisson')
  # user-specified q=2
  gfm2 <- gfm(dat$X, group, type, dropout = 2, q=2, output = F)
  measurefun(gfm2$hH, dat$H0, type='ccor')
  measurefun(gfm2$hB, dat$B0, type='ccor')
  #  select q automatically
  gfm2 <- gfm(dat$X, group, type, dropout = 2, q=NULL, q_set = 1:4,     output = F)
  measurefun(gfm2$hH, dat$H0, type='ccor')
  measurefun(gfm2$hB, dat$B0, type='ccor')

## mix of Poisson and Binomial
  dat <- gendata(seed=1, n=200, p=200, type='pois_bino', q=2, rho=3)
  group <- c(rep(1,ncol(dat$X)/2), rep(2,ncol(dat$X)/2))
  type <- c('poisson', 'binomial')
  #  select q automatically
  system.time(
    gfm3 <- gfm(dat$X, group, type, parallel = F,dropout = 2, q=NULL, q_set = 1:3, output = F, fast_version = T)
  )

  system.time(
    gfm3 <- gfm(dat$X, group, type, parallel = T,dropout = 2, q=NULL, q_set = 1:3, output = F, fast_version = T)
 )
}

\keyword{GFM}
