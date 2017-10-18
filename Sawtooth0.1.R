library(plyr)

##############
##Initialise##
##############

base <- 0.1
basesig <- base*1.2
sampsize <- c(50,100,500,1000)
sig <- c(1,1.25,1.5,1.75,2,2.5,3,3.5,4,4.5,5)
days <- 104
mvpd <- 1000
sdvpd <- 10
sims <- 500
probvax <- c(pnorm(c(1:20),mean=10,sd=2),1-pnorm(c(21:26),mean=23,sd=0.75))

res_0.1.df <- data.frame(base = base, n = rep(sampsize,11*sims), 
                          sim = rep(c(1:sims),each=44), 
                          ae = rep(rep(c("AE_1","AE_2","AE_3","AE_4","AE_5",
                                         "AE_6","AE_7","AE_8","AE_9","AE_10","AE_11"),each=4),sims), 
                          signal = 0, stime = days) 


##############################
##Summarises the simulation ##
##############################

sim_res <- function(df, ae, N) {
  
  vpd50 <- as.integer(rnorm(days, mean=50,sd=sdvpd*50/100)) 
  vpd100 <- as.integer(rnorm(days, mean=100,sd=sdvpd*100/100))
  vpd500 <- as.integer(rnorm(days, mean=500,sd=sdvpd*500/100))
  df1000 <- df
  
  ID50 <- list()
  ID100 <- list()
  ID500 <- list()
  
  for (i in 1:days) {
    ID50[[i]] <- sample(df1000$ID[which(df1000$Day==i)], size = min(vpd50[i],vpd[i]))
    ID100[[i]] <- sample(df1000$ID[which(df1000$Day==i)], size = min(vpd100[i],vpd[i]))
    ID500[[i]] <- sample(df1000$ID[which(df1000$Day==i)], size = min(vpd500[i],vpd[i]))
  }
  
  ID500  <- c(as.numeric(unlist(ID500)))
  ID100  <- c(as.numeric(unlist(ID100)))
  ID50  <- c(as.numeric(unlist(ID50)))
  
  DF_50_1 <- aggregate(as.formula(paste(ae,"~ Day")), data=df1000[which(df1000$ID %in% ID50 & df1000$vax==1),], sum)  
  DF_50_2 <- aggregate(as.formula(paste(N, "~ Day")), data=df1000[which(df1000$ID %in% ID50 & df1000$vax==1),], sum) 
  DF_100_1 <- aggregate(as.formula(paste(ae,"~ Day")), data=df1000[which(df1000$ID %in% ID100 & df1000$vax==1),], sum)  
  DF_100_2 <- aggregate(as.formula(paste(N, "~ Day")), data=df1000[which(df1000$ID %in% ID100 & df1000$vax==1),], sum) 
  DF_500_1 <- aggregate(as.formula(paste(ae,"~ Day")), data=df1000[which(df1000$ID %in% ID500 & df1000$vax==1),], sum)  
  DF_500_2 <- aggregate(as.formula(paste(N, "~ Day")), data=df1000[which(df1000$ID %in% ID500 & df1000$vax==1),], sum) 
  DF_1000_1 <- aggregate(as.formula(paste(ae,"~ Day")), data=df1000[which(df1000$vax==1),], sum)  
  DF_1000_2 <- aggregate(as.formula(paste(N, "~ Day")), data=df1000[which(df1000$vax==1),], sum) 
  
  res50 <- merge(DF_50_1,DF_50_2,by="Day",all.x=TRUE,all.y=TRUE)
  names(res50) <- c("Day","AE","vax")
  res50$N <- 50
  res100 <- merge(DF_100_1,DF_100_2,by="Day",all.x=TRUE,all.y=TRUE)
  names(res100) <- c("Day","AE","vax")
  res100$N <- 100
  res500 <- merge(DF_500_1,DF_500_2,by="Day",all.x=TRUE,all.y=TRUE)
  names(res500) <- c("Day","AE","vax")
  res500$N <- 500
  res1000 <- merge(DF_1000_1,DF_1000_2,by="Day",all.x=TRUE,all.y=TRUE)
  names(res1000) <- c("Day","AE","vax")
  res1000$N <- 1000
  
  resN <- rbind(res50,res100,res500,res1000)
  return(resN)
}


##########################################################
## Bayesian data analysis using the Beta/Binomial Model ##
##########################################################

Bayes_cal <- function(sumdata) {
  
  resN <- sumdata
  resN$cAE <-unlist(by(INDICES=factor(resN$N), data=resN$AE, FUN=cumsum))
  resN$cvax <-unlist(by(INDICES=factor(resN$N), data=resN$vax, FUN=cumsum))
  tempday <- seq(1:104)
  tempday50 <- tempday[which(!(tempday %in% resN$Day[which(resN$N==50)]))]
  tempday100 <- tempday[which(!(tempday %in% resN$Day[which(resN$N==100)]))]
  tempday500 <- tempday[which(!(tempday %in% resN$Day[which(resN$N==500)]))]
  tempday1000 <- tempday[which(!(tempday %in% resN$Day[which(resN$N==1000)]))]
  x_resN <- data.frame(Day=c(tempday50,tempday100,tempday500,tempday1000),AE=0,vax=0,
                       N=c(rep(50,length(tempday50)),rep(100,length(tempday100)),rep(500,length(tempday500)),rep(1000,length(tempday1000))),
                       cAE=NA,cvax=NA)
  resN <- rbind(resN,x_resN)
  resN <- resN[order(resN$N,resN$Day),]
  for (p in 1:length(resN$Day)) {
    if (is.na(resN$cAE[p])) {
      if (p==1) {resN$cAE[p] <- 0;resN$cvax[p] <- 0} else {
        if (resN$N[p] != resN$N[p-1]) {resN$cAE[p] <- 0; resN$cvax[p] <- 0} 
        else {resN$cAE[p] <- resN$cAE[p-1];resN$cvax[p] <- resN$cvax[p-1]} 
      }
    } else {}
  }
  
  
  resN$prior_a <- 2
  resN$prior_b <- 18
  resN$post_a <- resN$prior_a + resN$cAE
  resN$post_b <- resN$prior_b + resN$cvax - resN$cAE
  resN$post_mean <- resN$post_a/(resN$post_a+resN$post_b)
  resN$post_var <- (resN$post_a*resN$post_b)/((resN$post_a + resN$post_b +1)*(resN$post_a + resN$post_b)^2)
  resN$post_lb <- qbeta(p=0.025, shape1=resN$post_a, shape2=resN$post_b, lower.tail=TRUE)
  resN$post_ub <- qbeta(p=0.975, shape1=resN$post_a, shape2=resN$post_b, lower.tail=TRUE)
  return(resN)
} 

###################
##Data Simulation##
###################

for (s in 1:sims) {
  vpd <- as.integer(rnorm(days, mean=mvpd,sd=sdvpd*mvpd/100))
  sim<-array(data=NA,dim=c(sum(vpd),15))
  colnames(sim) <- c('ID','Day','vpd','vax','AE_1','AE_2','AE_3','AE_4','AE_5','AE_6','AE_7','AE_8','AE_9','AE_10','AE_11')
  
  sim[,'ID'] <- 1:NROW(sim)
  idx <- c(0,cumsum(vpd))
  for (j in 1:days) {
    for (k in 1:vpd[j]){
      idx1 <- (idx[j]+1):idx[j+1]#creating a "chunk"
      i <- idx[j]+k
      jday <- j %% 26
      if (jday==0) {jday <- 26}
      sim[idx1,'Day'] <- j 
      sim[idx1,'vpd'] <- rep(vpd[j],idx[j+1]-idx[j])
      sim[i,'vax'] <- rbinom(1,1,probvax[jday]) 
      
      # assigning the AEs 
      for (l in 1:11) {
        signal<-base*sig[l]
        if (j >= 27 & j <= 52 & sim[i,'vax']==1) {sim[i,4+l] <- rbinom(1,1,signal)} else
        {sim[i,4+l] <- rbinom(1,1,base)}
        
      }
    }
    
  }

  
  
  data<-as.data.frame(sim)
  sim_df<-data[,1:15]
  
  
##########################################
## Derive the operating characteristics ##
##########################################
  

  for (m in 1:11){
    
    res  <- sim_res(df=sim_df, ae= eval(paste("AE_",m,sep = "")),  N="vax")
    res_update <- Bayes_cal(sumdata=res)
    if (s==1 & m==1) {res_update_0.1_df <- res_update} else
    {res_update_0.1_df <- rbind(res_update_0.1_df,res_update)}

    for (size in 1:4){

      temp.res <- res_update[which(res_update$N==sampsize[size]),]
      if (length(temp.res$Day[which(temp.res$post_lb>basesig)]) > 0) {
        res_0.1.df$signal[(s-1)*44 + (m-1)*4 + size]  <- 1 
        first_time <- min(temp.res$Day[which(temp.res$post_lb>basesig)])
        res_0.1.df$stime[(s-1)*44 + (m-1)*4 + size] <- first_time
      }
    }
    
    
  }
  
  
  }

res_update_0.1_df <- cbind(res_update_0.1_df, 
                            sim=c(rep(1:sims, each=days*4*11)), 
                            RR=c(rep(rep(sig,each=days*4),sims)))

save.image("~/base0.1.st.df.RData")

