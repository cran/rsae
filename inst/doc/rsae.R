## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "",
    prompt = TRUE,
    fig.align = "center"
)
library("rsae")

## -----------------------------------------------------------------------------
data("landsat")

## -----------------------------------------------------------------------------
landsat[32:34,]

## -----------------------------------------------------------------------------
bhfmodel <- saemodel(formula = HACorn ~ PixelsCorn + PixelsSoybeans,
                     area = ~ CountyName,
                     data = subset(landsat, subset = (outlier == FALSE)))

## -----------------------------------------------------------------------------
mlfit <- fitsaemodel(method = "ml", bhfmodel)

## -----------------------------------------------------------------------------
mlfit

## -----------------------------------------------------------------------------
summary(mlfit)

## ----echo = FALSE-------------------------------------------------------------
set.seed(12345)
n <- 200; beta <- c(1, 1)
cst <- rep(1, n)
x <- rnorm(n)
y <- as.matrix(cbind(cst, x)) %*% beta + rnorm(n)
areaid <- rep(1:10, each=10)
df <- data.frame(y=y, x=x, areaid=areaid)
m <- saemodel(y ~ x, area=~areaid, data=df)
fitsaemodel("ml", m)

## -----------------------------------------------------------------------------
huberfit <- fitsaemodel("huberm", bhfmodel, k = 1.5)

## -----------------------------------------------------------------------------
huberfit

## -----------------------------------------------------------------------------
summary(huberfit)

## ----eval = FALSE-------------------------------------------------------------
#  robpredict(fit, areameans = NULL, k = NULL, reps = NULL, progress_bar = TRUE)

## -----------------------------------------------------------------------------
dat <- unique(landsat[-33, c("MeanPixelsCorn", "MeanPixelsSoybeans", "CountyName")])
dat <- cbind(rep(1,12), dat)
rownames(dat) <- dat$CountyName
dat <- dat[,1:3]
dat

## -----------------------------------------------------------------------------
pred <- robpredict(mlfit, areameans = dat)
pred

## -----------------------------------------------------------------------------
pred <- robpredict(mlfit, areameans = dat, reps = 100,
                   progress_bar = FALSE)
pred

## -----------------------------------------------------------------------------
est <- fitsaemodel("ml", bhfmodel, niter = 3)

## -----------------------------------------------------------------------------
est

## -----------------------------------------------------------------------------
convergence(est)

