library(xts)
library(Quandl)
library(quantmod)
library(TTR)
EURUSD <- Quandl("QUANDL/EURUSD", start_date="2014-01-01",end_date="2014-07-01", type="xts")##需注册账号
USDEUR <- Quandl("QUANDL/USDEUR", start_date="2014-01-01",end_date="2014-07-01", type="xts")#需注册账号
szSymbols <- c("MSFT","ORCL","GOOG","INTL","AAPL","CSCO","SYMC","TSLA")
getSymbols(szSymbols,src="yahoo",from="2008-1-1",to=Sys.Date()) 
getSymbols("CJSY",from="1900-1-1",to=Sys.Date()) 
print(head(CJSY));print(tail(CJSY))
getSymbols("399001.sz",from = "2002-01-01",to = Sys.Date(),src = "yahoo",auto.assign=FALSE)##下载深成指数

att <- getSymbols("T", from = "1985-01-01", to = "2015-12-31", auto.assign = FALSE)
plot(density(att), main = "Distribution of AT&T Returns")
rug(jitter(attr))

##模拟1维的布朗运动


#模拟2维的布朗运动
D2_Wiener <- function(){
  windows(200, 75)
  par(mfrow = c(1, 3), oma = c(0, 0, 2, 0))
  for(i in 1:3){
    W1 <- cumsum(rnorm(100000))
    W2 <- cumsum(rnorm(100000))
    plot(W1,W2, type= "l", ylab = "", xlab = "")
  }
  mtext("2-dimensional Wiener-processes", outer = TRUE, cex = 1.5, line = -1)
}

D2_Wiener()

##模拟指定相关系数的布朗运动
Correlated_Wiener <- function(cor){
  windows(200, 75)
  par(mfrow = c(1, 3), oma = c(0, 0, 2, 0))
  for(i in 1:3){
    W1 <- cumsum(rnorm(100000))
    W2 <- cumsum(rnorm(100000))
    W3 <- cor*W1+sqrt(1-cor^2)*W2
    plot(W1,W3, type= "l", ylab = "", xlab = "")
  }
  mtext("Correlated Wiener-processes", outer = TRUE, cex = 1.5, line = -1)
}
Correlated_Wiener(0.6)


###1维几何布朗运动的模拟方法
  
simGBM <- function (SO,mu,sigma,T,numSteps,numRep1)
{
       dt <- T/numSteps
       nuT <- (mu-sigma^2/2)*dt
       sigmaT <- sqrt(dt)*sigma
       pathMatrix = matrix(nrow=numRep1,ncol=numSteps+1)
       pathMatrix[,1] <- SO
       for (i in 1:numRep1){
         for (j in 2:numSteps){
           
           pathMatrix[i,j] <- pathMatrix[i,j-1]*
             exp(rnorm(1,nuT,sigmaT))
           
         }
       }
return(pathMatrix)
}
     # 想要生成几条模拟将 numRep1 赋值为几
SO <- 50; mu <- 0.1; sigma <- 0.3; T <- 1
numSteps <- 1000; numRep1 <- 10
     
set.seed(5)
     
     path = matrix(nrow=numRep1,ncol=numSteps+1)
     path <- simGBM(SO,mu,sigma,T,numSteps,numRep1 )
     
     for(i in 1:numRep1) {
       if (i == 1)
         plot(0:numSteps, path[i,], col = i, type='l', ylim = c(0,100))
       else
         lines(0:numSteps,path[i,], col = i)
     }  
     
     dim(path)
###作业：如何产生具有相关性的几个股票价格的价格路径  
     
#######正态性检验######
 x <- rnorm(3000)
 ks.test(x,"pnorm")     
 shapiro.test(x)
 w <- c(75.0, 64.0, 47.4, 66.9, 62.2, 62.2, 58.7, 63.5,
        66.6, 64.0, 57.0, 69.0, 56.9, 50.0, 72.0)
 qqnorm(w); qqline(w)
## nortest包
 library(nortest)
 ad.test(x)
 cvm.test(x)
 pearson.test(x)
 sf.test(x)
#####copula####
 library(copula)
 ## Generate some data (just skip this step, let's pretend we don't see it)
 n <- 1000 # sample size
 set.seed(271) # for reproducibility
 U <- rCopula(n, copula = gumbelCopula(iTau(gumbelCopula(), tau = 0.5))) # geek zone
 X <- qnorm(U) # quantile transformation
 X. <- qexp(U) # quantile transformation
 U  <- pobs(X)
 U. <- pobs(X.)
 plot(U,  xlab = expression(U[1]),     ylab = expression(U[2]))
 plot(U., xlab = expression(U[1]*"'"), ylab = expression(U[2]*"'"))
 layout(matrix(1:4, ncol = 2))
 plot(U[,1],  ylab = expression(U[1]))
 plot(U[,2],  ylab = expression(U[2]))
 plot(U.[,1], ylab = expression(U[1]*"'"))
 plot(U.[,2], ylab = expression(U[2]*"'"))
 layout(1)
 
 ## Plots
 layout(rbind(1:2))
 plot(X,  xlab = expression(X[1]),     ylab = expression(X[2]))
 plot(X., xlab = expression(X[1]*"'"), ylab = expression(X[2]*"'"))
 ####sklar定理的可视化
 library(mvtnorm)
 
 set.seed(271)
 
 ## Sample from a t copula
 n <- 1000 # sample size
 d <- 2 # dimension
 rho <- 0.7 # off-diagonal entry in the correlation matrix P
 P <- matrix(rho, nrow = d, ncol = d) # build the correlation matrix P
 diag(P) <- 1
 nu <- 3.5 # degrees of freedom
 set.seed(271)
 X <- rmvt(n, sigma = P, df = nu) # n multiv. t observations
 U <- pt(X, df = nu) # n ind. realizations from the corresponding t copula
 Y <- qexp(U) # transform U (t copula) to Exp(1) margins
 
 ## Plot setup
 ind <- c(A = 119, B = 516, C = 53) # use 'plot(X); identify(X)' to identify these points
 cols <- c("royalblue3", "maroon3", "darkorange2")
 par(pty = "s")
 
 ## Scatter plot of a bivariate t distribution
 plot(X, xlab = quote(X[1]), ylab = quote(X[2]))
 points(X[ind["A"],1], X[ind["A"],2], pch = 21, bg = cols[1], col = cols[1])
 points(X[ind["B"],1], X[ind["B"],2], pch = 21, bg = cols[2], col = cols[2])
 points(X[ind["C"],1], X[ind["C"],2], pch = 21, bg = cols[3], col = cols[3])
 text(X[ind["A"],1], X[ind["A"],2], labels = "A", adj = c(0.5, -0.75), col = cols[1])
 text(X[ind["B"],1], X[ind["B"],2], labels = "B", adj = c(0.5, -0.75), col = cols[2])
 text(X[ind["C"],1], X[ind["C"],2], labels = "C", adj = c(0.5, -0.75), col = cols[3])
 
 ## Scatter plot of the corresponding t copula (1st part of Sklar's Theorem)
 plot(U, xlab = quote(U[1]), ylab = quote(U[2]))
 points(U[ind["A"],1], U[ind["A"],2], pch = 21, bg = cols[1], col = cols[1])
 points(U[ind["B"],1], U[ind["B"],2], pch = 21, bg = cols[2], col = cols[2])
 points(U[ind["C"],1], U[ind["C"],2], pch = 21, bg = cols[3], col = cols[3])
 text(U[ind["A"],1], U[ind["A"],2], labels = "A", adj = c(0.5, -0.75), col = cols[1])
 text(U[ind["B"],1], U[ind["B"],2], labels = "B", adj = c(0.5, -0.75), col = cols[2])
 text(U[ind["C"],1], U[ind["C"],2], labels = "C", adj = c(0.5, -0.75), col = cols[3])
 
 ## Scatter plot of the meta-t distribution with Exp(1) margins (2nd part of Sklar's Theorem)
 plot(Y, xlab = quote(Y[1]), ylab = quote(Y[2]))
 points(Y[ind["A"],1], Y[ind["A"],2], pch = 21, bg = cols[1], col = cols[1])
 points(Y[ind["B"],1], Y[ind["B"],2], pch = 21, bg = cols[2], col = cols[2])
 points(Y[ind["C"],1], Y[ind["C"],2], pch = 21, bg = cols[3], col = cols[3])
 text(Y[ind["A"],1], Y[ind["A"],2], labels = "A", adj = c(0.5, -0.75), col = cols[1])
 text(Y[ind["B"],1], Y[ind["B"],2], labels = "B", adj = c(0.5, -0.75), col = cols[2])
 text(Y[ind["C"],1], Y[ind["C"],2], labels = "C", adj = c(0.5, -0.75), col = cols[3])
 ## We see that the *relative* locations of all points remains the same.
 ## We thus only change the marginal distributions, but not the dependence
 ## (by applying the probability and quantile transformations). This also
 ## visually confirms the invariance principle. 
 
 
 ##正态copula######
 n <- 1000 # sample size
 d <- 5 # max. considered dimension
 tau <- 0.5 # Kendall's tau
 th <- iTau(normalCopula(), tau = tau)
 nc <- normalCopula(th)
 
 ## Copula (wire frame and level curves)
 wireframe2(nc, FUN = pCopula)
 contourplot2(nc, FUN = pCopula)
 
 ## Copula density (wire frame and level curves)
 wireframe2(nc, FUN = dCopula, delta = 0.02)
 contourplot2(nc, FUN = dCopula)
 
 ## Scatter plots
 nc. <- normalCopula(th, dim = 5)
 U <- rCopula(n, copula = nc.)
 plot(U[,1:2], xlab = expression(U[1]), ylab = expression(U[2])) # d = 2
 cloud2(U[,1:3], # d = 3
        xlab = expression(U[1]), ylab = expression(U[2]), zlab = expression(U[3]))
 pairs2(U, cex = 0.4, col = tblack(0.5)) # d = 5
 
 library(copula)
 library(gridExtra) # for grid.arrange()
 n <- 1000 # sample size
 set.seed(271)
 ## Define the copulas定义
 nc <- normalCopula(th.n) # Gauss copula
 gc <- gumbelCopula(th.g) # Gumbel copula
 cc <- claytonCopula(th.c) # Clayton copula
 tc <- tCopula(th.t) # t_4 copula
 ## Generate copula data
 U.nc <- rCopula(n, copula = nc)
 U.gc <- rCopula(n, copula = gc)
 U.cc <- rCopula(n, copula = cc)
 U.tc <- rCopula(n, copula = tc)
 
 ## 从Copula映射到N(0,1) Map to N(0,1) margins (meta-copula data)
 X.nc <- qnorm(U.nc)
 X.gc <- qnorm(U.gc)
 X.cc <- qnorm(U.cc)
 X.tc <- qnorm(U.tc)
 
 ##相关矩阵 Correlations
 cors <- vapply(list(X.nc, X.gc, X.cc, X.tc), function(x) cor(x)[1,2], NA_real_)
 stopifnot(all.equal(cors, rep(0.7, 4), tol = 0.015))
 
 
 
 ## Define the corresponding (meta-C) densities (via Sklar's Theorem)
 ## Note: Density f(x_1, x_2) = c(F_1(x_1), F_2(x_2)) * f_1(x_1) * f_2(x_2)
 ##                           = exp( log(c(F_1(x_1), F_2(x_2))) + log(f_1(x_1)) + log(f_2(x_2)) )
 dMetaCopulaN01 <- function(x, copula)
   exp(dCopula(pnorm(x), copula = copula, log = TRUE) + rowSums(dnorm(x, log = TRUE)))
 ## Alternatively, we could work with dMvdc() here
 
 ### 3 Explicit copulas #########################################################
 
 ### 3.1 Clayton copula #########################################################
 
 ## Define the Clayton copula object
 th <- iTau(claytonCopula(), tau = tau)
 cc <- claytonCopula(th)
 
 ## Copula (wire frame and level curves)
 wireframe2(cc, FUN = pCopula)
 contourplot2(cc, FUN = pCopula)
 
 ## Copula density (wire frame and level curves)
 wireframe2(cc, FUN = dCopula, delta = 0.02)
 contourplot2(cc, FUN = dCopula)
 
 ## Scatter plots
 cc. <- claytonCopula(th, dim = 5)
 U <- rCopula(n, copula = cc.)
 plot(U[,1:2], xlab = expression(U[1]), ylab = expression(U[2])) # d = 2
 cloud2(U[,1:3], # d = 3
        xlab = expression(U[1]), ylab = expression(U[2]), zlab = expression(U[3]))
 pairs2(U, cex = 0.4, col = tblack(0.5)) # d = 5
 
 
 ### 3.2 Gumbel copula ##########################################################
 
 ## Define the Gumbel copula object
 th <- iTau(gumbelCopula(), tau = tau)
 gc <- gumbelCopula(th)
 
 ## Copula (wire frame and level curves)
 wireframe2(gc, FUN = pCopula)
 contourplot2(gc, FUN = pCopula)
 
 ## Copula density (wire frame and level curves)
 wireframe2(gc, FUN = dCopula, delta = 0.02)
 contourplot2(gc, FUN = dCopula)
 
 ## Scatter plots
 gc. <- gumbelCopula(th, dim = 5)
 U <- rCopula(n, copula = gc.)
 plot(U[,1:2], xlab = expression(U[1]), ylab = expression(U[2])) # d = 2
 cloud2(U[,1:3], # d = 3
        xlab = expression(U[1]), ylab = expression(U[2]), zlab = expression(U[3]))
 pairs2(U, cex = 0.4, col = tblack(0.5)) # d = 5
 
 ## Scatter plot of a survival Gumbel copula
 pairs2(1-U, cex = 0.4, col = tblack(0.5))
 
 ## Only the first component flipped
 pairs2(U, cex = 0.4, col = tblack(0.5)) # original
 pairs2(cbind(1-U[,1], U[,2:5]), cex = 0.4, col = tblack(0.5)) # first flipped
 
 
 ### 3.3 Outer power Clayton copula #############################################
 
 ## Note: Outer power Clayton copulas can have both lower and upper tail dependence
 
 ## Define the outer power Clayton copula
 tauC <- 0.3 # Kendall's tau for the underlying Clayton copula
 thC <- copClayton@iTau(tauC) # choose Clayton's generator s.t. Kendall's tau is tauC
 opC <- opower(copClayton, thC) # define an outer power Clayton copula (its parameter theta is not specified yet)
 th <- opC@iTau(tau) # define the outer power Clayton copula s.t. has Kendall's tau is tau
 opcc <- onacopulaL(opC, list(th, 1:2)) # define the outer power Clayton copula
 
 ## Scatter plot
 U <- rCopula(n, copula = opcc)
 plot(U[,1:2], xlab = expression(U[1]), ylab = expression(U[2])) # d = 2
 
 
 ### 3.4 Marshall--Olkin copulas ################################################
 
 ## Note: Marshall--Olkin copulas have a singular component, i.e., a set of
 ##       Lebesgue measure 0 with positive probability mass assigned.
 
 ## Define the MO copula
 C <- function(u, alpha) {
   if(!is.matrix(u)) u <- rbind(u)
   stopifnot(0 <= alpha, alpha <= 1, length(alpha) == 2)
   pmin(u[,1]*u[,2]^(1-alpha[2]), u[,1]^(1-alpha[1])*u[,2])
 }
 
 ## Sampling
 rMO <- function(n, alpha) {
   stopifnot(n >= 1, 0 <= alpha, alpha <= 1, length(alpha) == 2)
   U. <- matrix(runif(n*3), ncol = 3)
   U <- cbind(pmax(U.[,1]^(1/(1-alpha[1])), U.[,3]^(1/alpha[1])),
              pmax(U.[,2]^(1/(1-alpha[2])), U.[,3]^(1/alpha[2])))
 }
 
 ## Define the singular component
 S.C <- function(u, alpha) {
   stopifnot(0 <= u, u <= 1,
             0 <= alpha, alpha <= 1, length(alpha) == 2)
   tau <- alpha[1] * alpha[2] / (alpha[1] + alpha[2] - alpha[1] * alpha[2])
   tau * pmin(u[,1]^alpha[1], u[,2]^alpha[2])^(1/tau)
 }
 
 ## Define the absolutely continuous component
 A.C <- function(u, alpha) C(u, alpha) - S.C(u, alpha)
 
 
 ## Scatter plot
 alpha <- c(0.2, 0.8)
 U <- rMO(n, alpha = alpha)
 plot(U, xlab = expression(U[1]), ylab = expression(U[2]))
 ## Interpretation: Given u_1 = 0.6, for example, u_2 lies with a non-zero
 ## probability p on the curve and with the remaining probability 1-p
 ## anywhere else (*not*: uniform) along the line u_1 = 0.6.
 
 ## Check the margins
 plot(U[,1], ylab = expression(U[1]))
 plot(U[,2], ylab = expression(U[2]))
 
 ## Evaluate the copula (and its singular and absolutely continuous components)
 grid <- expand.grid("u[1]" = u, "u[2]" = u) # build a grid
 val.C <- cbind(grid, "C(u[1],u[2])" = C(grid, alpha = alpha)) # append C
 val.S <- cbind(grid, "S[C](u[1],u[2])" = S.C(grid, alpha = alpha)) # append S.C
 val.A <- cbind(grid, "A[C](u[1],u[2])" = A.C(grid, alpha = alpha)) # append S.C
 
 ## Copula (wire frame and level curves)
 wireframe2(val.C) # wire frame plot
 contourplot2(val.C, xlim = 0:1, ylim = 0:1) # level curves
 ## The copula has a kink, that is, a so-called singular component.
 ## A singular component is a set of Lebesgue measure 0 where C puts
 ## mass at (here: a curve). This is better visible from the scatter plots below.
 
 ## Singular and absolutely continuous component (wire frame)
 wireframe2(val.S) # singular component
 wireframe2(val.A) # absolutely continuous component
 
 ######copula的估计##############
 library(xts)
 library(copula)
 library(qrmdata)
 library(qrmtools)
 
 
 ### 1 Working with the data ####################################################
 
 ## Load the time series
 data("SP500")
 data("FTSE")
 
 ## Build negative log-returns
 X.SP500 <- -returns(SP500)
 X.FTSE  <- -returns(FTSE)
 
 ## Merge
 X <- merge(X.SP500, X.FTSE, all = FALSE)
 time <- c("2003-01-01", "2012-12-31")
 X <- X[paste0(time, collapese = "/"),]
 
 ## Basic plot
 plot.zoo(X)
 
 ## Observe that there are zero values caused by market closures
 apply(X == 0, 2, sum)
 rexcl <- (X[,1] == 0) | (X[,2] == 0)
 X. <- X[!rexcl,]
 
 ## Aggregating by week
 dim(X.w <- apply.weekly(X., FUN = colSums))
 plot.zoo(X.w, type = "h")
 
 ## Compute pseudo-observations
 U <- as.matrix(pobs(X.w))
 plot(U, xlab = expression(U[1]), ylab = expression(U[2]))
 
 
 ### 2 Fitting copulas ##########################################################
 
 ## Compare various bivariate copulas
 fit.N <- fitCopula(normalCopula(),  data = U)
 fit.t <- fitCopula(tCopula(),       data = U) # df of freedom are estimated, too
 fit.C <- fitCopula(claytonCopula(), data = U)
 fit.G <- fitCopula(gumbelCopula(),  data = U)
 
 ## Comparing the likelihoods
 sort(c(N = fit.N@loglik, t = fit.t@loglik, C = fit.C@loglik, G = fit.G@loglik),
      decreasing = TRUE)
 
 