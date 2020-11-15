# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

## Homogeneous variance normal
q <- 2
dat <- gendata(q = q, n=500, p=500, rho=3)
X <- dat$X
t(dat$B0) %*% dat$B0
Fac <- Factorm(X)
group <- rep(1,ncol(X))
type <- 'gaussian'
dropout <- 0
eps2 <- 1e-3
maxIter <- 30
gfm1 <- gfm_eval_intercept_init(dat$X, group, type, q,
                                    dropout, eps2, maxIter=10,
                                output=0, parallel = F)

measurefun(gfm1$hH, dat$H0, type='ccor')
measurefun(gfm1$hB, dat$B0, type='ccor')
# cbind(gfm1$hmu, dat$mu0)
measurefun(gfm1$hH, dat$H0, type='fnorm')

gfm2final <-gfm_eval_intercept_osfinal(X, gfm1$hH, gfm1$hB,gfm1$hmu, group, type)
measurefun(gfm2final$hH, dat$H0, type='ccor')
measurefun(gfm2final$hB, dat$B0, type='ccor')
measurefun(gfm2final$hH, dat$H0, type='fnorm')

# select q automatically
gfm3 <- gfm(X, group, type, q=NULL, q_set = 1:6, output = F)




# Heterogeous variance normal
dat <- gendata(seed=1, n=100, p=100, type='heternorm', q=2, rho=4)
group <- rep(1,ncol(dat$X))
type <- 'gaussian'
dropout <- 0
eps2 <- 1e-3
maxIter <- 30
gfm2 <- gfm_eval_intercept_init(dat$X, group, type, q=2,
                                dropout, eps2, maxIter=10,
                                output=0, parallel = T)
gfm2
measurefun(gfm2$hH, dat$H0, type='ccor')
Xc <- scale(dat$X, scale = F)
measurefun(Factorm(Xc,q)$hH, dat$H0, type='ccor')
measurefun(gfm2$hH, dat$H0, type='fnorm')
measurefun(gfm2$hB, dat$B0, type='ccor')
measurefun(Factorm(Xc,q=2)$hB, dat$B0, type='ccor')

#  select q automatically
gfm3 <- gfm(dat$X, group, type, q=NULL, q_set = 1:4, output = F)

# Poisson variable
q <- 3; p <- 200
dat <- gendata(seed=1, n=200, p=p, type='pois', q=q, rho=4)
X <- dat$X
t(dat$B0) %*% dat$B0 /p
group <- rep(1,ncol(dat$X))
type <- 'poisson'
dropout <- 0
eps2 <- 1e-3
maxIter <- 30

system.time(
  gfm2 <- gfm_eval_intercept_init(dat$X, group, type, q=q,
                                  dropout, eps2, maxIter=30,
                                  output=T, parallel = F)
)
system.time(
  gfm2 <- gfm_eval_intercept_init(dat$X, group, type, q=q,
                                  dropout, eps2, maxIter=30,
                                  output=F, parallel = T)
)
gfm2
measurefun(gfm2$hH, dat$H0, type='ccor')
measurefun(gfm2$hB, dat$B0, type='ccor')
measurefun(gfm2$hH, dat$H0, type='fnorm')
Xc <- scale(dat$X, scale = F)
measurefun(Factorm(Xc,q=q)$hH, dat$H0, type='ccor')

#  select q automatically
gfm3 <- gfm(dat$X, group, type, q=NULL, q_set = 1:4, output = T, fast_version = T)


# mix of normal and Poisson
dat <- gendata(seed=1, n=200, p=200, type='norm_pois', q=2, rho=2)
t(dat$B0) %*% dat$B0
group <- c(rep(1,ncol(dat$X)/2), rep(2,ncol(dat$X)/2))
type <- c('gaussian','poisson')
dropout <- 1
eps2 <- 1e-3
maxIter <- 30
gfm2 <- gfm_eval_intercept_init(dat$X, group, type, q=2,
                                dropout=dropout, eps2=eps2, maxIter=30, output=T)
measurefun(gfm2$hH, dat$H0, type='ccor')
measurefun(gfm2$hH, dat$H0, type='fnorm')
measurefun(gfm2$hB, dat$B0, type='ccor')
gfm2final <- gfm_eval_intercept_osfinal(dat$X, gfm2$hH, gfm2$hB,gfm2$hmu, group, type)
measurefun(gfm2final$hH, dat$H0, type='ccor')
measurefun(gfm2final$hH, dat$H0, type='fnorm')
measurefun(gfm2final$hB, dat$B0, type='ccor')
#  select q automatically
gfm3 <- gfm(dat$X, group, type, dropout = 2, q=NULL, q_set = 1:4, output = F)

# mix of Poisson and Binomial
dat <- gendata(seed=1, n=200, p=200, type='pois_bino', q=2, rho=3)
X <- dat$X
group <- c(rep(1,ncol(dat$X)/2), rep(2,ncol(dat$X)/2))
type <- c('poisson', 'binomial')
dropout <- 1
eps2 <- 1e-3
gfm2 <- gfm_eval_intercept_init(dat$X, group, type, q=2,
                                dropout, eps2, maxIter=10, output=T)

measurefun(gfm2$hH, dat$H0, type='ccor')
measurefun(gfm2$hB, dat$B0, type='ccor')

gfm2final <- gfm_eval_intercept_osfinal(dat$X, gfm2$hH, gfm2$hB,gfm2$hmu, group, type)
measurefun(gfm2final$hH, dat$H0, type='ccor')
measurefun(gfm2final$hB, dat$B0, type='ccor')

#  select q automatically
system.time(
  gfm3 <- gfm(dat$X, group, type, parallel = F,dropout = 2, q=NULL, q_set = 1:3, output = F, fast_version = T)

)

# select q automatically in parallel
system.time(
  gfm3 <- gfm(dat$X, group, type, parallel = T,dropout = 2, q=NULL, q_set = 1:3, output = F, fast_version = T)

)



# 单变量并行计算test
myexp <- function(x) exp(1+x)
mymean <- function(x) mean(x+1)
mylog <- function(x, base=2, t=1,s=2){

  log(x, myexp(2+t*mymean(s)))
}
single_parallel(mylog,1:10000, varlist = c('myexp', 'mymean'), base=3, t=1, s=1)
# 多变量并行测试代码
multi_parallel(paste,c("a","b"),c("c","d"),c('E','F'),MoreArgs = list(sep="_"))

paraIC(4,X, group, type)
varlist <- c('X','gfm_eval_intercept_init', 'ICriteria', "Factorm",
             "localupdateB2","family2func","ortheB","Diag",
             "localupdateH2","ortheH", "objfunc")
single_parallel(paraIC,iterable =1:3, varlist= varlist, XX=X, group=group, type=type,
                dropout=0, eps2=1e-4, maxIter=10, omega=1/ncol(X), output=F, fast_version=T)

