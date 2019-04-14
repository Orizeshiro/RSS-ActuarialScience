claim <- read.csv("/Users/timosie/Desktop/精算模型2018/claim.csv",header = T,stringsAsFactors = F)
#描述性统计画直方图
par(mfrow=c(1,1))
claims <- claim$Claim
max(claim$Claim)
hist(claims,freq = T,breaks = seq(0,13000000,by=1000000))

x<-claim$Claim
n<-length(x)

#指数分布检验
fexp<-function(p,x,n){
  a=-sum(x)/p-n*log(p)  #似然函数值
  return(a)
}
xmax <- optimize(fexp,c(1,100),tol = 0.0001, maximum = TRUE,x,n)$maximum
xmax

Fnx<-ecdf(x)
plot.ecdf(Fnx)
curve(pexp(x,1/xmax),xlim =c(0,14),ylab = "Fn(x)",add=TRUE) 
pexp(cumsum(table(x)/(length(x)+1)),1/xmax) #计算模型得到的分布函数的分布密度

cumsum(table(x)/length(x)) #计算原始的分布密度

plot(pexp(cumsum(table(x)/(length(x)+1)),1/xmax)-cumsum(table(x)/length(x)),ylab = "F(x)-Fn(x)",xlab = "x")
abline(h=0)

#pp图
a<-pexp(cumsum(table(x)/(length(x)+1)),1/xmax) #纵坐标
b<-cumsum(table(x)/(length(x)+1)) #横坐标
plot(b,a,xlim=c(0,1),ylim = c(0,0.8))
abline(0,1)  #对角线

#qq图
c<-qexp(b,1/xmax) #纵坐标，理论数据
plot(unique(x),c,xlim = c(0,2000),ylim = c(0,2000))
abline(0,1)

library(goftest)
library(MASS)
library(VGAM)
library(actuar)
#K-S检验
# Kolmogorov-Smirnov # the test statistic is "D"
ks.test(claim$Claim, "pexp")
#AD检验
ad.test(claim$Claim, "pexp")
#卡方检验
chisq.test(unique(claim$Claim),pexp(table(x)/(length(x)),xmax))

#帕累托分布检验
fpareto<-function(alpha,theta)
{
  ll <- n*log(alpha)+n*alpha*log(theta)+(alpha+1)*sum(log(x+theta))
  return(-ll)
}

pareto_max <- stats4::mle(minuslogl = fpareto,start=list(alpha = 1,theta = 1))
summary(pareto_max)
alpha <- 21009432
theta <- 5580958#求极大似然估计值

ppareto(cumsum(table(x)/(length(x)+1)),alpha,theta) #计算模型得到的分布函数的分布密度

plot(ppareto(cumsum(table(x)/(length(x)+1)),alpha,theta)-cumsum(table(x)/length(x)),ylab = "F(x)-Fn(x)",xlab = "x")
abline(h=0)
#pp图
a<-ppareto(cumsum(table(x)/(length(x)+1)),alpha,theta) #纵坐标
b<-cumsum(table(x)/(length(x)+1)) #横坐标
plot(b,a)
abline(0,1)  #对角线
#KS检验
ks.test(claim$Claim, "ppareto",shape=alpha,scale=theta)
#AD检验
ad.test(claim$Claim, "ppareto",shape=alpha,scale=theta)
#卡方检验
chisq.test(unique(claim$Claim),ppareto(table(x)/(length(x)),alpha,theta))

#伽马分布检验
fgamma<-function(param)
{
  ll <- (param[1]-1)*sum(log(x))-sum(x)/param[2]-n*param[1]*log(param[2])-n*log(prod(1:(param[1]-1)))
  return(-ll)
}

gamma_max <- optim(param <- c(10,5),fgamma,hessian = FALSE)
gamma_max$par
alpha <- gamma_max$par[1]
theta <- gamma_max$par[2]

pgamma(cumsum(table(x)/(length(x)+1)),alpha,theta)

plot(pgamma(cumsum(table(x)/(length(x)+1)),alpha,theta)-cumsum(table(x)/length(x)),ylab = "F(x)-Fn(x)",xlab = "x")
abline(h=0)

#pp图
a<-pgamma(unique(x),alpha,theta) #纵坐标
b<-cumsum(table(x)/(length(x)+1)) #横坐标
plot(b,a)
abline(0,1)  #对角线
#KS检验
ks.test(claim$Claim, "pgamma",shape=alpha,scale=theta)
#AD检验
ad.test(claim$Claim, "pgamma",shape=alpha,scale=theta)
#卡方检验
chisq.test(unique(claim$Claim),pgamma(table(x)/(length(x)),alpha,theta))

#广义beta分布
fgbeta<-function(param)
{
  ll <- n*log(prod(1:(param[1]+param[2]-1))/(prod(1:(param[1]-1))*prod(1:(param[2]-1))))
  +param[1]*param[4]*(sum(log(x))-n*log(param[3]))
  +(param[2]-1)*sum(log(1-(x/param[3])^param[4]))+n*log(param[4])-sum(log(x))
  return(-ll)
}

gbeta_max <- optim(param <- c(10,5,2,1),fgbeta,hessian = FALSE)
gbeta_max$par
a <- gbeta_max$par[1]
b <- gbeta_max$par[2]
theta <- gbeta_max$par[3]
t <- gbeta_max$par[4]

pgenbeta(cumsum(table(x)/(length(x)+1)),a,b,theta,t)

plot(pgenbeta(unique(x),a,b,theta,t)-cumsum(table(x)/length(x)),ylab = "F(x)-Fn(x)",xlab = "x")
abline(h=0)

#pp图
a<-pgenbeta(cumsum(table(x)/(length(x)+1)),a,b,theta,t) #纵坐标
b<-cumsum(table(x)/(length(x)+1)) #横坐标
plot(b,a)
abline(0,1)  #对角线
#KS检验
ks.test(claim$Claim, "pgenbeta",a,b,theta,t)
#AD检验
ad.test(claim$Claim, "pgenbeta",shape=alpha,scale=theta)
#卡方检验
chisq.test(unique(claim$Claim),pgenbeta(table(x)/(length(x)),a,b,theta,t))

