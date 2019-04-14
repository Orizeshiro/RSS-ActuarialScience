#用卷积法计算索赔额分布
#加载actuar包
library(actuar)
#索赔次数N的分布
pn = dbinom(0:5,5,0.5)
#索赔次数N的分布,由题意得a=2.25,ceta=2500
qpareto(0.99,2.25,2500)
#先计算帕累托分布得0.99分位数，以确定计算密度函数时得取值范围
fx = dpareto(seq(0,16856,by=10),2.25,2500)
#卷积法计算总索赔额复合分布
Fs1 = aggregateDist("convolution", model.freq = pn, model.sev = fx) 
num=seq(10,10000,by=10)
amount=c()
for (i in 1:1000)
{
  amount[i]=Fs1(10*i)
}
plot(num,amount,"l")

#用递推法计算索赔额分布
fx = dpareto(seq(0,5000,by=100),2.25,2500)
Fs2 = aggregateDist("recursive", model.freq = "binomial", model.sev = fx, prob = 0.5,size=5)
num=seq(10,10000,by=10)
amount=c()
for (i in 1:1000)
{
  amount[i]=Fs2(10*i)
}
plot(num,amount,"l")

#傅立叶变换法计算索赔额分布
x <- seq(0, 16000, by = 10)  #等间隔离散化，步长为10
F_par <- ppareto(x, shape = 2.25, scale = 2500) #产生pareto的分布函数
f_par <- diff(F_par)  #离散化以后伽马分布的概率
f_par <- c(f_par, 1-ppareto(10000, shape = 2, scale = 500))  #补足概率
sum(f_par)
f_par

plot(x, f_par, type = 's', col = 2)
f_par <- c(f_par, rep(0, 100))  #右尾补充足够的零
phi_x=fft(f_par)  #损失金额的特征函数
phi_s=exp(2.5*(fft(f_par)-1))  #累计损失的特征函数
fs <- Re(fft(phi_s, inverse = TRUE)/length(f_par) ) #累积索赔金额的分布
sum(fs)

Fs <- cumsum(fs) #求出s的分布
par(mfrow=c(1,2))
plot(10*(0:(length(fs)-1)), fs, type = 'h', col=2, xlab='s')
plot(10*(0:(length(Fs)-1)), Fs, type = 's', col=2, xlab='s')

#通过随机数求总索赔额分布
alpha <- 2.25
theta <- 2500
par(mfrow=c(1,2))

n <- 1:10
fn <- dbinom(1:10,2)

plot(n,fn,ylim=c(0,0.3),main="Frequency: Poisson")
abline(h=0,lty=2)

x <- seq(1,16000,10)
fx <- alpha*theta^alpha/(x+theta)^(alpha+1)

plot(x,fx,type="l",main="Severity: Pareto")

set.seed(123)
size <- 5000
S <- rep(NA,size)
N <- rpois(size,2)
for (i in 1:size){
  uu <- runif(N[i])
  X <- theta*((1-uu)^(-1/alpha)-1)
  S[i] <- sum(X)
}

par(mfrow=c(1,2))
hist(S,freq=F,breaks=100)
plot(ecdf(S),xlab="S")




