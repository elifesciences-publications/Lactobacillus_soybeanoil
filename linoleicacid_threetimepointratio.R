## R commands for Figures 4E, 4–4, 2, 2-1A,B
library("plyr")
library("lme4")
library("ggbeeswarm")

## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#select the data file
data<-read.delim(file="linoleicacidgrowthdata.txt",header=TRUE)


#blank the no-oil and oil values for each time point
nooilcolumnA<-which(names(data)=="timeA_nooil")
nooilcolumnE<-which(names(data)=="timeE_nooil")
oilcolumnA<-which(names(data)=="timeA_oil")
oilcolumnE<-which(names(data)=="timeE_oil")
data["timeA_nooilblanked"]<-"NA"
data["timeB_nooilblanked"]<-"NA"
data["timeC_nooilblanked"]<-"NA"
data["timeD_nooilblanked"]<-"NA"
data["timeE_nooilblanked"]<-"NA"
data["timeA_oilblanked"]<-"NA"
data["timeB_oilblanked"]<-"NA"
data["timeC_oilblanked"]<-"NA"
data["timeD_oilblanked"]<-"NA"
data["timeE_oilblanked"]<-"NA"


for(samp in 1:dim(data)[1]) {
  #blanking the no-oil samples
  ODvaluesnoil<-data[samp,nooilcolumnA:nooilcolumnE] 
  ODvaluesnoil<-as.matrix(ODvaluesnoil) #ready for plotting
  
  blanknames<-strsplit(as.character(data$blank_nooil[samp]),split=", ")[[1]]
  plate<-as.character(data$plate[samp])  
  plategrep<-paste("\\b",plate,"\\b",sep="")
  platespecific<-data[grep(plategrep,data$plate),]
  blanklocation1<-which(platespecific$location_nooil == blanknames[1])
  blanklocation2<-which(platespecific$location_nooil == blanknames[2]) 
  #check blank values
 # par(mfrow=c(1,2))
  #plot(timemins, as.matrix(platespecific[blanklocation1,nooilcolumnA:nooilcolumnE]),pch=16)
  #plot(timemins, as.matrix(platespecific[blanklocation2,nooilcolumnA:nooilcolumnE]))
  
  blankvalues<-as.matrix(platespecific[blanklocation1,nooilcolumnA:nooilcolumnE]+platespecific[blanklocation2,nooilcolumnA:nooilcolumnE])/2
  
  blanked_nooil<-ODvaluesnoil-blankvalues
  timeA_nooilblankedlocation<-which(names(data)=="timeA_nooilblanked")
  timeE_nooilblankedlocation<-which(names(data)=="timeE_nooilblanked")
  names(blanked_nooil)<-paste(names(blanked_nooil),"_blanked", sep="")
  data[samp,timeA_nooilblankedlocation:timeE_nooilblankedlocation]<-as.matrix(blanked_nooil)
  #plot(timemins,(blankedvaluesnoil))
  
  #blank the oil samples
  ODvaluesoil<-data[samp,oilcolumnA:oilcolumnE] 
  ODvaluesoil<-as.matrix(ODvaluesoil) #ready for plotting
  
  blanknames<-strsplit(as.character(data$blank_oil[samp]),split=", ")[[1]]
  plate<-as.character(data$plate[samp])  
  plategrep<-paste("\\b",plate,"\\b",sep="")
  platespecific<-data[grep(plategrep,data$plate),]
  
  blanklocation1<-which(platespecific$location_oil == blanknames[1])
  blanklocation2<-which(platespecific$location_oil == blanknames[2]) 
  #check blank values
  #par(mfrow=c(1,2))
  #plot(timemins, as.matrix(platespecific[blanklocation1,oilcolumnA:oilcolumnE]),pch=16)
  #plot(timemins, as.matrix(platespecific[blanklocation2,oilcolumnA:oilcolumnE]))
  
  blankvalues<-as.matrix(platespecific[blanklocation1,oilcolumnA:oilcolumnE]+platespecific[blanklocation2,oilcolumnA:oilcolumnE])/2
  
  blanked_oil<-ODvaluesoil-blankvalues
  timeA_oilblankedlocation<-which(names(data)=="timeA_oilblanked")
  timeE_oilblankedlocation<-which(names(data)=="timeE_oilblanked")
  names(blanked_oil)<-paste(names(blanked_oil),"_blanked", sep="")
  data[samp,timeA_oilblankedlocation:timeE_oilblankedlocation]<-as.matrix(blanked_oil)
  
  rm(plategrep,blanked_nooil,blankvalues,ODvaluesnoil,ODvaluesoil,platespecific,plate,blanklocation1,blanklocation2,blanknames,timeA_nooilblankedlocation,timeA_oilblankedlocation,timeE_nooilblankedlocation,timeE_oilblankedlocation)
  }

#the calculations for the ratio of growth in oil to that in no oil for each of the last three timepoints:
for (samp2 in 1:dim(data)[1]) {
  if (as.numeric(data[samp2,]$timeC_nooilblanked)>0) {
    data[samp2,"Coil_Cnoil"]<-as.numeric(data[samp2,]$timeC_oilblanked)/as.numeric(data[samp2,]$timeC_nooilblanked) } else {data[samp2,"Coil_Cnoil"]<-"NA" }
  if (is.na(data[samp2,]$timeD_nooilblanked) == "FALSE") {if (as.numeric(data[samp2,]$timeD_nooilblanked)>0) {
    data[samp2,"Doil_Dnoil"]<-as.numeric(data[samp2,]$timeD_oilblanked)/as.numeric(data[samp2,]$timeD_nooilblanked) } else {data[samp2,"Doil_Dnoil"]<-"NA" }}
  if (as.numeric(data[samp2,]$timeE_nooilblanked)>0) {
    data[samp2,"Eoil_Enoil"]<-as.numeric(data[samp2,]$timeE_oilblanked)/as.numeric(data[samp2,]$timeE_nooilblanked) } else {data[samp2,"Eoil_Enoil"]<-"NA" }
} 
  
data["mouse_colonytype"]<-paste(data$mouse,".",data$colonytype, sep="")
data["cage_mouse_colonytype"]<-paste(data$cage,".",data$mouse,".",data$colonytype, sep="")
data["diet_colonytype"]<-paste(data$diet,".",data$colonytype, sep="")
data["ratiovariance"]<-NA
data["meanCDEratio"]<-NA

#summarize the growth ratios for the last three timepoints as mean and variance
ratiovarcol<-which(names(data) == "ratiovariance")
meanCDEratiocol<-which(names(data) == "meanCDEratio")
ratioEcol<-which(names(data) == "Eoil_Enoil")
ratioDcol<-which(names(data) == "Doil_Dnoil")
ratioCcol<-which(names(data) == "Coil_Cnoil")
for(samp in 1:dim(data)[1]) {
  data[samp,ratiovarcol]<-var(as.numeric(c(data[samp,ratioEcol],data[samp,ratioDcol],data[samp,ratioCcol])), na.rm=TRUE)
  data[samp,meanCDEratiocol]<-mean(as.numeric(c(data[samp,ratioEcol],data[samp,ratioDcol],data[samp,ratioCcol])), na.rm=TRUE)
  }

mean(data$ratiovariance,na.rm=TRUE)

#subset just the experimental samples and ignore the controls
nocontrols<-subset(data, data$isolation.source !="control" &data$isolation.source !="blank" &data$study !="blank")
nocontrols<-droplevels(nocontrols)

#these samples did not grow well in the no oil treatment and should be removed (see below)
which(nocontrols$timeE_nooilblanked<.1)  ###important!!!

grew<-subset(nocontrols, nocontrols$timeE_nooilblanked >=.1 )  #exclude samples that did not grow
grew<-droplevels(grew)
mean(grew$ratiovariance,na.rm=TRUE)

#subset the data into the mouse diet study and the naturally isolated strains from Jens Walter
AP4<-subset(grew,grew$study == "AP4" & grew$colonytype != "" &grew$colony != "")
AP4<-droplevels(AP4)
JW<-subset(grew,grew$study == "JW")

#summarize the mouse diet study data by mouse diet and Lactobacillus type; white colonies are L. reuteri; non-white colonies are L. johnsonii; these observations were confirmed by full length 16S rRNA gene sequencing
AP4_mouse_colonytype_summary<-summarySE(AP4, measurevar="meanCDEratio", groupvars=c("diet_colonytype","mouse", "colonytype", "colony", "diet", "body.part", "cage") )
AP4_mouse_colonytype_summarysorted<-arrange(AP4_mouse_colonytype_summary,xtfrm(diet_colonytype),xtfrm(colonytype),xtfrm(mouse))
AP4_mouse_colonytype_summarysorted$diet_colonytype <- factor(AP4_mouse_colonytype_summarysorted$diet_colonytype, levels=c("LF.white","HF.white","LF.nonwhite","HF.nonwhite") )

my_mean2=aggregate(AP4_mouse_colonytype_summarysorted$meanCDEratio, by=list(AP4_mouse_colonytype_summarysorted$diet_colonytype) , mean) ; colnames(my_mean2)=c("diet_colonytype" , "mean")
my_sd2=aggregate(AP4_mouse_colonytype_summarysorted$meanCDEratio , by=list(AP4_mouse_colonytype_summarysorted$diet_colonytype) , sd) ; colnames(my_sd2)=c("diet_colonytype" , "sd")
my_info2=merge(my_mean2 , my_sd2 , by.x=1 , by.y=1)


AP4plotting<-AP4_mouse_colonytype_summary
AP4plotting<-arrange(AP4plotting,xtfrm(diet_colonytype),xtfrm(colonytype),xtfrm(mouse))
AP4plotting$diet_colonytype <- factor(AP4plotting$diet_colonytype, levels=c("LF.white","HF.white", "LF.nonwhite","HF.nonwhite") )


#here is one that uses the mean of the summarized data:

#Figure 4E
ggplot(AP4plotting) +     
  geom_point(aes(x=diet_colonytype , y = meanCDEratio,colour = diet, shape=colonytype),size=3  ) +
  geom_errorbar(data = AP4plotting, aes(x = diet_colonytype, y = sd, ymin = meanCDEratio - sd, ymax = meanCDEratio + sd, colour=diet), width=.01) +
  geom_point(data = my_info2, aes(x=diet_colonytype , y = mean),  shape =95, size = 10, stroke=5,colour="black") +
  geom_errorbar(data = my_info2, aes(x = diet_colonytype, y = sd, ymin = mean - sd, ymax = mean + sd), width=.2, size=0.8) +
  xlab("") +
  ylab("18:2 growth ratio") +
   scale_colour_manual(values=c("orange", "green")) + scale_fill_manual(values=c("orange", "green")) +
  theme_classic() +
  theme(
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.x= element_text(size=20),
    axis.title.y= element_text(size=20),
    axis.text.x= element_text(size=15),
    axis.text.y= element_text(size=15)
  )


colorpalatte<-c("mediumspringgreen", "blue","cornflowerblue","springgreen4", "tan","plum","firebrick1","orange","skyblue1","purple4","yellow","magenta")

#Figure 4 – Figure supplement 4 
ggplot(AP4plotting) +     
  geom_point(aes(x=diet_colonytype , y = meanCDEratio, fill=cage, colour=cage, shape=colonytype), size=3  ) +
  geom_errorbar(data = AP4plotting, aes(x = diet_colonytype, y = sd, ymin = meanCDEratio - sd, ymax = meanCDEratio + sd), width=.01, colour="black") +
  geom_point(data = my_info2, aes(x=diet_colonytype , y = mean),  shape =95, size = 10, stroke=5,colour="black") +
  geom_errorbar(data = my_info2, aes(x = diet_colonytype, y = sd, ymin = mean - sd, ymax = mean + sd), width=.2, size=0.8) +
  xlab("") +
  ylab("18:2 growth ratio") +
  scale_colour_manual(values=colorpalatte, name="Cage")+
  scale_fill_manual(values=colorpalatte,name="Cage") +
  guides(shape=FALSE) +
  scale_shape_manual(values=c(21,24))+
  #scale_colour_manual(values=c("orange", "green","orange", "green","orange", "green","orange", "green","orange", "green","orange", "green")) + #scale_fill_manual(values=c("black", "black")) +
  theme_classic() +
  theme(
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.x= element_text(size=20),
    axis.title.y= element_text(size=20),
    axis.text.x= element_text(size=15),
    axis.text.y= element_text(size=15)
  )


kruskal.test(meanCDEratio~diet_colonytype, data=subset(AP4plotting, AP4plotting$colonytype=="white"))
# Kruskal-Wallis rank sum test
# 
# data:  meanCDEratio by diet_colonytype
# Kruskal-Wallis chi-squared = 4.2529, df = 1, p-value = 0.03918
# 
t.test(meanCDEratio~diet_colonytype, data=subset(AP4plotting, AP4plotting$colonytype=="white"))
# 
# 
# Welch Two Sample t-test
# 
# data:  meanCDEratio by diet_colonytype
# t = -2.1685, df = 125.94, p-value = 0.032
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.085507310 -0.003906846
# sample estimates:
#   mean in group LF.white mean in group HF.white 
# 0.04364834             0.08835542 


kruskal.test(meanCDEratio~diet_colonytype, data=subset(AP4plotting, AP4plotting$colonytype=="nonwhite"))
# Kruskal-Wallis rank sum test
# 
# data:  meanCDEratio by diet_colonytype
# Kruskal-Wallis chi-squared = 13.971, df = 1, p-value = 0.0001857
# 
t.test(meanCDEratio~diet_colonytype, data=subset(AP4plotting, AP4plotting$colonytype=="nonwhite"))

# Welch Two Sample t-test
# 
# data:  meanCDEratio by diet_colonytype
# t = 3.2216, df = 158.14, p-value = 0.001548
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.02723796 0.11355121
# sample estimates:
#   mean in group LF.nonwhite mean in group HF.nonwhite 
# 0.14728707                0.07689248 



#naturally isolated strains from Jens Walter

JW<-droplevels(JW)


#summarize the data by host
JWsummarybystrain<-summarySE(JW, measurevar = "meanCDEratio", groupvars=c("Strain","isolation.source", "Species") )
JWsummarybystrain_sorted<-arrange(JWsummarybystrain,xtfrm(isolation.source),xtfrm(Strain))
JWsummarybystrain_sorted$isolation.source <- factor(JWsummarybystrain_sorted$isolation.source, levels=c("mouse","rat","pig","human","chicken","turkey", "sourdough") )


#to plot for each host separately
ggplot(subset(JWsummarybystrain_sorted, JWsummarybystrain_sorted$Species =="L. reuteri"), aes(x=isolation.source, y=meanCDEratio)) + 
  geom_point(size=4,stat="identity") +
  geom_errorbar(aes(ymin=meanCDEratio-sd, ymax=meanCDEratio+sd), width=.05,colour="black") +
  xlab("source") +
  ylab("ratio") +
  #scale_colour_manual(values=c("orange", "green")) + scale_fill_manual(values=c("orange", "green")) +
  theme_classic() +
  theme(
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.x= element_text(size=20),
    axis.title.y= element_text(size=20),
    axis.text.x= element_text(size=15),
    axis.text.y= element_text(size=15)
  )

#need to group the hosts by rodent and nonrodent 
JWsummarybystrain_sorted_LR<-subset(JWsummarybystrain_sorted, JWsummarybystrain_sorted$Species =="L. reuteri")
JWsummarybystrain_sorted_LR["hostgroup"]<-NA
JWsummarybystrain_sorted_LR["rodentgroup"]<-NA


for (samp3 in 1:dim(JWsummarybystrain_sorted_LR)[1]) {
  if (JWsummarybystrain_sorted_LR[samp3,]$isolation.source=="mouse" | JWsummarybystrain_sorted_LR[samp3,]$isolation.source=="rat") {
    JWsummarybystrain_sorted_LR[samp3,]$hostgroup<-"rodents";
    JWsummarybystrain_sorted_LR[samp3,]$rodentgroup<-"yes"
  } else if (JWsummarybystrain_sorted_LR[samp3,]$isolation.source=="pig" ) {
    JWsummarybystrain_sorted_LR[samp3,]$hostgroup<-"porcine";
    JWsummarybystrain_sorted_LR[samp3,]$rodentgroup<-"no"
  } else if (JWsummarybystrain_sorted_LR[samp3,]$isolation.source=="human" ) {
    JWsummarybystrain_sorted_LR[samp3,]$hostgroup<-"human";
    JWsummarybystrain_sorted_LR[samp3,]$rodentgroup<-"no"
  } else if (JWsummarybystrain_sorted_LR[samp3,]$isolation.source=="chicken" | JWsummarybystrain_sorted_LR[samp3,]$isolation.source=="turkey" ) {
    JWsummarybystrain_sorted_LR[samp3,]$hostgroup<-"poultry";
    JWsummarybystrain_sorted_LR[samp3,]$rodentgroup<-"no"
  } else if (JWsummarybystrain_sorted_LR[samp3,]$isolation.source=="sourdough" ) {
   JWsummarybystrain_sorted_LR[samp3,]$hostgroup<-"sourdough";
    JWsummarybystrain_sorted_LR[samp3,]$rodentgroup<-"no"
  }
}


#reorganize the data
JWsummarybystrain_sorted_LR<-arrange(JWsummarybystrain_sorted_LR,xtfrm(hostgroup))
JWsummarybystrain_sorted_LR$hostgroup <- factor(JWsummarybystrain_sorted_LR$hostgroup, levels=c("rodents","porcine","human","poultry", "sourdough") )


is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}

#add additional information about the naturally isolated strains
JWinfo<-read.delim(file="JWstraininfo.txt",header=TRUE, na.strings = "NA")

JWsummarybystrain_sorted_LR["Clade"]<-NA
JWsummarybystrain_sorted_LR["Humansite"]<-NA
JWsummarybystrain_sorted_LR["Country"]<-NA


for(sample in 1:dim(JWsummarybystrain_sorted_LR)[1]) {
  samplename<-JWsummarybystrain_sorted_LR[sample,1]
  samplerow<-which(JWinfo$Strain==as.character(samplename))
  if (!is.integer0(samplerow)) {
    clade<-as.character(JWinfo[samplerow,]$Clade)
    humansite<-as.character(JWinfo[samplerow,]$Human.site)
    country<-as.character(JWinfo[samplerow,]$Country)
    JWsummarybystrain_sorted_LR[sample,]$Clade<-clade
    JWsummarybystrain_sorted_LR[sample,]$Humansite<-humansite
    JWsummarybystrain_sorted_LR[sample,]$Country<-country
  }
}


#Figure 2
ggplot(data=JWsummarybystrain_sorted_LR, aes(x=hostgroup, y=meanCDEratio)) + 
  geom_point(size=4,stat="identity") +
  geom_errorbar(aes(ymin=meanCDEratio-sd, ymax=meanCDEratio+sd), width=.05,colour="black") +
  xlab("host source") +
  ylab("18:2 grwoth ratio") +
  #scale_colour_manual(values=c("orange", "green")) + scale_fill_manual(values=c("orange", "green")) +
  theme_classic() +
  theme(
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.x= element_text(size=20),
    axis.title.y= element_text(size=20),
    axis.text.x= element_text(size=15),
    axis.text.y= element_text(size=15)
  )

JWsummarybystrain_sorted_LR$rodentgroup<-as.factor(JWsummarybystrain_sorted_LR$rodentgroup)
kruskal.test(meanCDEratio ~ rodentgroup , data=JWsummarybystrain_sorted_LR)
# Kruskal-Wallis rank sum test
# 
# data:  meanCDEratio by rodentgroup
# Kruskal-Wallis chi-squared = 15.313, df = 1, p-value = 9.11e-05.

t.test(meanCDEratio ~ rodentgroup , data=JWsummarybystrain_sorted_LR)
# Welch Two Sample t-test
# 
# data:  meanCDEratio by rodentgroup
# t = 4.0327, df = 28.607, p-value = 0.0003733
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.1162317 0.3557453
# sample estimates:
#   mean in group no mean in group yes 
# 0.4541600         0.2181715



humansubset<-subset(JWsummarybystrain_sorted_LR, JWsummarybystrain_sorted_LR$isolation.source=="human")
humansubset<-droplevels(humansubset)

#Figure 2 – Figure supplement 1A. 
ggplot(data=humansubset, aes(x=Humansite, y=meanCDEratio) ) + 
  geom_point(size=4,stat="identity") +
  geom_errorbar(aes(ymin=meanCDEratio-sd, ymax=meanCDEratio+sd), width=.05) +
  xlab("human site") +
  ylab("18:2 grwoth ratio") +
  #scale_colour_manual(values=c("orange", "green")) + scale_fill_manual(values=c("orange", "green")) +
  theme_classic() +
  theme(
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.x= element_text(size=20),
    axis.title.y= element_text(size=20),
    axis.text.x= element_text(size=15),
    axis.text.y= element_text(size=15)
  )


#Figure 2 – Figure supplement 1B. 
ggplot(data=JWsummarybystrain_sorted_LR, aes(x=Clade, y=meanCDEratio, colour=hostgroup)) + 
  geom_point(size=4,stat="identity") +
  geom_errorbar(aes(ymin=meanCDEratio-sd, ymax=meanCDEratio+sd), width=.05) +
  xlab("Clade") +
  ylab("18:2 grwoth ratio") +
  theme_classic() +
  theme(
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.x= element_text(size=20),
    axis.title.y= element_text(size=20),
    axis.text.x= element_text(size=15),
    axis.text.y= element_text(size=15)
  )


