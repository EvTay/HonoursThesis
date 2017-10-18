library(survival)
library(rms)
res_1<-res_0.1.df[res_0.1.df$ae=="AE_6",]
surv.fit2 <- survfit(Surv(res_1$stime,res_1$signal) ~ res_1$n,conf.type="log-log")
plot(surv.fit2,xlim = c(200,800),ylim = c(0,1.0),fun = function(x) 1-x,pch = "n",col=1:4,main = "Underlying rate 0.1, RR 2.5",ylab = "Cumulative Probability of Signal Detection",xlab = "Time(Weeks)")
abline(h=0.5,col = "red")
text(x=27,y=0.53,label="Median")
legend(x=700,y=0.2, legend = c("50/Day","100/Day","500/Day","1000/Day"), col=1:4, pch=1,bty="n")

