library(plyr)

##############
##Initialise##
##############

base <- 0.1
basesig <- base*1.2
sampsize <- c(50,100,500,1000)
sig <- c(1,1.25,1.5,1.75,2,2.5,3,3.5,4,4.5,5)#relative risk ratios
vax <- 1.0 #everyone in the SmartVax data is vaccinated
days <- 1000
mvpd <- 1000
sdvpd <- 15
sims <- 500

res_0.1.df <- data.frame(base = base, n = rep(sampsize,5500), #500*11
                         sim = rep(c(1:sims),each=44), #4*11 - recording the outputs for each n and for each AE for each of 1 in 500 sims
                         ae = rep(rep(c("AE_1","AE_2","AE_3","AE_4","AE_5","AE_6","AE_7","AE_8","AE_9","AE_10","AE_11"),each=4),500), 
                         signal = 0, stime = days) 


##############
##Summarises the simulation ##
##############

sim_res <- function(df, ae, N ) {
  
  vpd50 <- as.integer(rnorm(days, mean=50,sd=sdvpd*50/100)) # the number sampled each day when it's 50/day
  vpd100 <- as.integer(rnorm(days, mean=100,sd=sdvpd*100/100))
  vpd500 <- as.integer(rnorm(days, mean=500,sd=sdvpd*500/100))
  df1000 <- df
  
  ID50 <- list()
  ID100 <- list()
  ID500 <- list()
  for (i in 1:days) {
    ID50[[i]] <- sample(df1000$ID[which(df1000$Day==i)], size = min(vpd50[i],vpd[i]))
    ID100[[i]] <- sample(df1000$ID[which(df1000$Day==i)], size = min(vpd100[i],vpd[i]))
    ID500[[i]] <- sample(df1000$ID[which(df1000$Day==i)], size = min(vpd500[i],vpd[i]))# which individuals to be sampled each day
  }
  
  ID500  <- c(as.numeric(unlist(ID500)))
  ID100  <- c(as.numeric(unlist(ID100)))
  ID50  <- c(as.numeric(unlist(ID50)))
  
  DF_50_1 <- aggregate(as.formula(paste(ae,"~ Day")), data=df1000[which(df1000$ID %in% ID50),], sum)  
  DF_50_2 <- aggregate(as.formula(paste(N, "~ Day")), data=df1000[which(df1000$ID %in% ID50),], sum) 
  DF_100_1 <- aggregate(as.formula(paste(ae,"~ Day")), data=df1000[which(df1000$ID %in% ID100),], sum)  
  DF_100_2 <- aggregate(as.formula(paste(N, "~ Day")), data=df1000[which(df1000$ID %in% ID100),], sum) 
  DF_500_1 <- aggregate(as.formula(paste(ae,"~ Day")), data=df1000[which(df1000$ID %in% ID500),], sum)  
  DF_500_2 <- aggregate(as.formula(paste(N, "~ Day")), data=df1000[which(df1000$ID %in% ID500),], sum) 
  DF_1000_1 <- aggregate(as.formula(paste(ae,"~ Day")), data=df1000, sum)  
  DF_1000_2 <- aggregate(as.formula(paste(N, "~ Day")), data=df1000, sum) 
  
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
  
  # x50 <- bayes_analysis(res50)
  # x100 <- bayes_analysis(res100)
  # x500 <- bayes_analysis(res500)
  # x1000 <- bayes_analysis(res1000)
  
  resN <- rbind(res50,res100,res500,res1000)
  return(resN)
}


##############
## Calculates the Bayes ##
##############

Bayes_cal <- function(sumdata) {
  resN <- sumdata
  resN$cAE <-unlist(by(INDICES=factor(resN$N), data=resN$AE, FUN=cumsum))
  resN$cvax <-unlist(by(INDICES=factor(resN$N), data=resN$vax, FUN=cumsum))
  resN$prior_a <- 1
  resN$prior_b <- 1
  resN$post_a <- resN$prior_a + resN$cAE
  resN$post_b <- resN$prior_b + resN$cvax - resN$cAE
  resN$post_mean <- resN$post_a/(resN$post_a+resN$post_b)
  resN$post_var <- (resN$post_a*resN$post_b)/((resN$post_a + resN$post_b +1)*(resN$post_a + resN$post_b)^2)
  resN$post_lb <- qbeta(p=0.025, shape1=resN$post_a, shape2=resN$post_b, lower.tail=TRUE)
  resN$post_ub <- qbeta(p=0.975, shape1=resN$post_a, shape2=resN$post_b, lower.tail=TRUE)
  return(resN)
} 



##############
##Simulate the data (aka "the big loop") ##
##############


for (s in 1:sims) {
  vpd <- as.integer(rnorm(days, mean=mvpd,sd=sdvpd*mvpd/100))
  sim<-array(data=NA,dim=c(sum(vpd),18))
  colnames(sim) <- c('ID','Day','vpd','Batch_a','Batch_b','vax','Batch','AE_1','AE_2','AE_3','AE_4','AE_5','AE_6','AE_7','AE_8','AE_9','AE_10','AE_11')
  
  sim[,'ID'] <- 1:NROW(sim)
  idx <- c(0,cumsum(vpd))
  for (j in 1:days) {
    for (k in 1:vpd[j]){
      idx1 <- (idx[j]+1):idx[j+1]#creating a "chunk"
      i <- idx[j]+k
      sim[idx1,'Day'] <- j
      sim[idx1,'vpd'] <- rep(vpd[j],idx[j+1]-idx[j])
      sim[i,'Batch_a'] <- rbinom(1,1,pnorm( j,mean=0.2*days,sd=25 ))+1#batch1 to 2 
      sim[i,'Batch_b'] <- rbinom(1,1,pnorm( j,mean=0.5*days,sd=25 ))+2#batch 2 to 3 
      sim[i,'vax'] <- rbinom(1,1,vax) 
      #assigning the Batch
      if(sim[i,'Batch_a']==1)  {sim[i,'Batch']<-1} else {
        sim[i,'Batch']<-sim[i,'Batch_b']*sim[i,'vax']}
      # assigning the AEs 
      for(l in 1:11){
        signal<-base*sig[l]
        if (sim[i,'Batch']==1) {sim[i,7+l] <- rbinom(1,1,base)} else {
          if (sim[i,'Batch']==2) {sim[i,7+l] <- rbinom(1,1,signal)} else {
            if (sim[i,'Batch']==3) {sim[i,7+l] <-rbinom(1,1,base)}
          }
        }
      }
      
    }
  }
  
  data<-as.data.frame(sim)
  sim_df<-data[,1:18]
  
  
  ##############
  ##  Derive the operating characteristics ##
  ##############
  
  
  for (m in 1:11){#the output is being sorted by the ae
    
    res  <- sim_res(df=sim_df, ae= eval(paste("AE_",m,sep = "")),  N="vax")
    res_update <- Bayes_cal(sumdata=res)
    
    if (s==1 & m==1) {res_update_0.1_df <- res_update} else
    {res_update_0.1_df <- rbind(res_update_0.1_df,res_update)}#500*4*104
    
    
    for (size in 1:4){
      
      temp.res <- res_update[which(res_update$N==sampsize[size]),]
      if (length(temp.res$Day[which(temp.res$post_lb>basesig)]) > 0) {
        res_0.1.df$signal[(s-1)*44 + (m-1)*4 + size]  <- 1 #a signal has been detected
        first_time <- min(temp.res$Day[which(temp.res$post_lb>basesig)])#time to signal detection
        res_0.1.df$stime[(s-1)*44 + (m-1)*4 + size] <- first_time
             }
    }
    
    
  }
  
}

res_update_0.1_df <- cbind(res_update_0.1_df, 
                           sim=c(rep(1:sims, each=days*4*11)), #1000*4*11=44000
                           RR=c(rep(rep(sig,each=days*4),sims)))#sims*4000


save.image("~/base0.1df.RData")

