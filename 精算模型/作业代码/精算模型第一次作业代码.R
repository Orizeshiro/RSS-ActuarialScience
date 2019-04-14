#生成常见损失分布随机数
#正态分布
rnorm(n=10,mean=0,sd=1)
#均匀分布
runif(10,min=0,max=10)
#指数分布
rexp(10,rate=0.2)
#伽马分布
rgamma(10,rate=1,10)
#pareto分布
set.seed(10)
rpareto(10,3,2)
#对数正态分布
rlnorm(10,0,1)
#泊松分布
rpois(10,2)
#二项分布
rbinom(10,10,0.6)
#负二项分布
rnbinom(10,10,0.6)


#计算常见分布的矩、E(X^d)和E(X-d|X>d),以Pareto分布为例
library(actuar)
mpareto(shape=2,scale=1500,order=1)
b <- levpareto(1000,2,1500,order=1)
b
##设a=E(X-d|X>d),根据公式得
a <- (1500-b)/(0.4)^2
a
