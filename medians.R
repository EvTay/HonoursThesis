
medianlist <- list()
sampsize <- c(50,100,500,1000)
                
for (p in 1:4) {
    ss <- sampsize[p]
    for (q in 1:11){
        median <- median(res_0.1.df$stime[which(res_0.1.df$ae ==eval(paste("AE_",q,sep = "")) & res_0.1.df$n == ss)])
        medianlist[(p-1)*11+q] <- median
      
       
    }
}
out.median <- data.frame(SampSize=rep(sampsize, each=11), RR=rep(sig,4), Mediandays=unlist(medianlist))
out.median
