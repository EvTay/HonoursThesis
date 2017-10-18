library(ggplot2)
sim_df$Batch <- as.factor(sim_df$Batch)
sumbatch_f <-data.frame(unlist(table(sim_df$Batch,sim_df$Day)), stringsAsFactors=FALSE)
sumbatch_n <-data.frame(unlist(table(sim_df$Day)),stringsAsFactors=FALSE)
names(sumbatch_f) <- c('Batch','Day','Freq')
names(sumbatch_n) <- c('Day','TotFreq')
sumbatch <- merge(sumbatch_f,sumbatch_n,by="Day",all.x=TRUE, all.y=TRUE)
sumbatch$Day <- as.numeric(levels(sumbatch$Day))[sumbatch$Day]
sumbatch$Batch <- as.factor(levels(sumbatch$Batch))[sumbatch$Batch]


p0 <- ggplot(data=sumbatch, aes(x=Day, y=Freq/TotFreq, group=Batch, colour=Batch)) +
  ylab("Proportion of Children Vaccinated") + geom_point() +
  theme(legend.title = element_text(size=16,face="bold"),#change title of legend
        #legend.title = element_blank(),
        legend.text=element_text(size=12, face="bold"),
        axis.title.x = element_text(size=14, face="bold", colour="black"),
        axis.text.y = element_text(size=12, face="bold",colour="black"),
        axis.title.y = element_text(size=14, face="bold",colour="black"),
        axis.text.x = element_text(size=12, face="bold",colour="black"))
