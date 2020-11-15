# RGFM
Generalized factor model for ultra-high dimensional variables with mixed types --R version

## Load R package GFM
```{r, eval=FALSE}
library(GFM)
example(gfm)
```

# Some Examples
```{r}
library(GFM)

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

```
