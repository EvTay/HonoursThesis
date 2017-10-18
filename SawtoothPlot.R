library(ggplot2)

batch2prob <- data.frame(days = c(1:1000), bprob = NA)
batch2prob$bprob <- pnorm( batch2prob$days,mean=0.2*days,sd=25 ) + (1-pnorm( batch2prob$days,mean=0.5*days,sd=25 ))
plot(I(bprob-1)~days,data=batch2prob, ylab="Probability of Receiving Batch 2")

y_int <- 0.0012

lbplots <- function(ss) {
  
  avsim_time_summary1 <- aggregate(post_lb ~ Day + RR,data=res_update_0.001_df[which(res_update_0.001_df$N==ss),],FUN = "mean") 
  avsim_time_summary2 <- aggregate(post_mean ~ Day + RR,data=res_update_0.001_df[which(res_update_0.001_df$N==ss),],FUN = "mean") 
  avsim_time_summary <- merge(avsim_time_summary1, avsim_time_summary2, by=c("Day","RR"))
  
  p.ss <- ggplot()+
    geom_line(data=avsim_time_summary,aes(x=Day,y=post_lb, col=factor(RR)),lty=1) +
    geom_line(data=avsim_time_summary,aes(x=Day,y=post_mean, col=factor(RR)),lty=2) +
    geom_hline (aes (yintercept = y_int), lty=2) +
    ggtitle(paste("Vaccine uptake = ",ss," /week",sep="")) +
    scale_x_continuous(name = "Time(Weeks)") +
    scale_y_continuous(name = "Probability of AEFI", labels=function(n){format(n, scientific = FALSE)},
                       limits=c(0,0.0030))+
    scale_colour_discrete(name="Relative\nRisk")+
    scale_linetype_discrete(labels = c("95% Credible Interval LB","Posterior Mean"))+
    theme(axis.title.x = element_text(face="bold", size=16),
          axis.text.x  = element_text(size=14),
          axis.title.y = element_text(face="bold", size=14),
          axis.text.y  = element_text(size=12),
          legend.position=c(0.1,0.73),
          legend.text = element_text(size=12),
          legend.title = element_text(size=12),
          plot.title = element_text(face="bold",size=18)) +
          annotate("text", x = 95, y = y_int*1.05, label = "Threshold") +
          annotate("text", x = 15, y = 0.0009, label = "Posterior Mean") +
          annotate("text", x = 35, y = 0.0004, label = "Posterior Lower Bound")
  return(p.ss)
  
  
}

