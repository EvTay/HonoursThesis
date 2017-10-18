########AusVaxSafety Data
fever <- read.csv("Documents/FeverData.csv", header=TRUE, stringsAsFactors=FALSE)
fever[is.na(fever)] <- 0
fever$Date.f <- as.Date(fever$Date,"%d/%m/%Y")

plot(No.Records ~ Date.f,data=fever,pch=20,ylab = "Number of Records", xlab = "Date(year-month)")

library(ggplot2)
ggplot(fever, aes(x=Date.f, y=No.Records)) +
  ggtitle("Children Immunised With the Influenza Vaccine")+
  ylab("Number of Records")+
  xlab("Date(year-month)")+
  theme(axis.title.x = element_text(face="bold", size=16),
              axis.text.x  = element_text(size=14),
              axis.title.y = element_text(face="bold", size=14),
              axis.text.y  = element_text(size=16),
              plot.title = element_text(face="bold",size=18)) +
  geom_point()



plot(No.Events/No.Records~Date.f,data=fever, ylim=c(0,0.1),pch=20,
     ylab = "Proportion of Reported Numbers of Fever",
     xlab = "Date")
abline(h=0.031,col="red",lty=2)
legend("topleft", legend = "Mean Incidence of Fever Over the Three Years", col="red", lty=2, bty="n")

ggplot(fever, aes(x=Date.f, y=No.Events/No.Records)) +
  ggtitle("Incidence of Fever in Children Immunised with the Influenza Vaccine")+
  ylab("Proportion of Reported Numbers of Fever")+
  xlab("Date(year-month)")+
  ylim(0,0.1)+
  geom_hline (aes (yintercept = 0.031), color="red", lty=2,show.legend=TRUE)+
  scale_colour_manual(name="Legend",values=c("red"="Mean Incidence"))+
  theme(axis.title.x = element_text(face="bold", size=16),
        axis.text.x  = element_text(size=14),
        axis.title.y = element_text(face="bold", size=14),
        axis.text.y  = element_text(size=16),
        plot.title = element_text(face="bold",size=18))+
   annotate("text",y=0.033,x=min(fever$Date.f,na.rm=TRUE)+200,label="Mean Incidence of Fever")+
   geom_point()

#####Bayes_cal function
prior_a_init <- 2
prior_b_init <- (2/base)-prior_a_init
prior_a_init <- 1
prior_b_init <- 1


Bayes_cal <- function(resN) {
  
  resN$prior_a <- prior_a_init
  resN$prior_b <- prior_b_init
  resN$post_a <- resN$prior_a + resN$cAE
  resN$post_b <- resN$prior_b + resN$cvax - resN$cAE
  resN$post_mean <- resN$post_a/(resN$post_a+resN$post_b)
  resN$post_var <- (resN$post_a*resN$post_b)/((resN$post_a + resN$post_b +1)*(resN$post_a + resN$post_b)^2)
  resN$post_lb <- qbeta(p=0.025, shape1=resN$post_a, shape2=resN$post_b, lower.tail=TRUE)
  resN$post_ub <- qbeta(p=0.975, shape1=resN$post_a, shape2=resN$post_b, lower.tail=TRUE)
  return(resN)
} 

#######2016 season using a prior informed by the 2015 season

#######2015 season
fever2015 <- subset(fever, as.numeric(rownames(fever))<12)
Bayes_cal(fever2015)

#######2017 season using a prior informed by the 2016 season
fever2017 <- subset(fever, as.numeric(rownames(fever)) >=32)
rownames(fever2017) <- 1:nrow(fever2017)
fever2017$Week <- as.numeric(rownames(fever2017))
base <- 0.032 #the proportion of cAE/cvax at end of season 2016
basesig <- base*1.2
Bayes_cal(fever2017)

ggplot(Bayes_cal(fever2017), aes(x=Week, y=post_lb)) +
  ggtitle("Probability of AEFI by Week in 2017")+
  ylab("Probability of AEFI:95%CI Lower Bound")+
  xlab("Week")+
  ylim(0,0.1)+
  geom_hline(aes(yintercept=0.0384),lty=2,col="red")+
  theme(axis.title.x = element_text(face="bold", size=16),
        axis.text.x  = element_text(size=14),
        axis.title.y = element_text(face="bold", size=14),
        axis.text.y  = element_text(size=16),
        plot.title = element_text(face="bold",size=18)) +
    annotate("text",y = 0.04, x = 10,label="Threshold")+
  geom_line()



plot(post_lb ~ Week,data = Bayes_cal(fever2017), ylim = c(0,0.1),type='l',
     ylab="Probability of AEFI: 95% CI Lower Bound",xlab="Week")
abline(h=0.0384,lty=2,col="red")
legend("topleft", legend = "Threshold", col="red", lty=2, bty="n")

