
theta_true <- 0.8
###Generate population
N = sample(seq(100000,400000),1)
A = round(theta_true*N)
B = N - A
Zpop <- sample(c(rep(1,A),rep(0,B)))
###Pull random sample from population
N_samp <- 200
Zsamp <- sample(Zpop,N_samp)
Y_samp <- sum(Zsamp)

###Functions to obtain the arguments representing the parameters of the Beta distributions
pa.fun<-function(m,n){
  pa=m*n
  return(pa)
}
pb.fun<-function(m,n){
  pb=n*(1-m)
  return(pb)
}
la.fun<-function(N_samp,Y_samp){
  la = Y_samp+1
  return(la)
}
lb.fun<-function(N_samp,Y_samp){
  lb=N_samp - Y_samp +1
  return(lb)
}
poafun<-function(m,n,N_samp,Y_samp){
  poa = Y_samp + (m*n)
  return(poa) 
}
pobfun<-function(m,n,N_samp,Y_samp){
  pob = N_samp - Y_samp +(n*(1-m)) 
  return(pob)
}

###Plot of the beta distributions
library(ggplot2)
betaplot_m <- function(pa,pb,la,lb,poa,pob){
  theta = seq(0,1,0.005)
  pr  = dbeta(theta,pa,pb)
  lk  = dbeta(theta,la,lb)
  po = dbeta(theta,poa,pob)##this should be a beta binom m = prob of success
  plot.dat<- data.frame(pi=rep(theta,3),pdf = c(pr,lk,po), 
                        dist=rep(c("Prior","Likelihood","Posterior"),each=201))
  p <- ggplot(data=plot.dat,aes(x=pi,y=pdf, colour = factor(dist) )) +
    geom_line(size=1.5)+
    ylab("PDF") + ylim(0,40) +
    xlab(expression(theta))+
    scale_color_discrete(name="Distribution")+
    theme(legend.position="top",
          legend.title=element_text( size=16, face="bold"),
          legend.text=element_text(size=16, face="bold"),
          axis.text.x = element_text(size=16, face="bold",colour="black"),
          axis.title.x = element_text(size=16, face="bold", colour="black"),
          axis.text.y = element_text(size=16, face = "bold", colour="black"),
          axis.title.y = element_text(size=16, face="bold",colour="black") )
  
  return(p)}


betaplot_m(200,800,157,45,355,843)##n=1000 m=0.2
betaplot_m(60,240,157,45,215,283)##n=300
betaplot_m(20,80,157,45,175,123)##n=100 m=0.2
