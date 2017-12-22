library("plyr")
library("reshape2")
library("beeswarm")


is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}
#open mapping file
setwd("~/Box Sync/FAMI/16S/FAMIAP1234")
mappingfile<-read.delim(file="mappingfile_FAMIAP1234_all.txt",sep="\t",header = TRUE)
#mappingfile<-read.delim(file="mappingfile.txt",sep="\t",header = TRUE)
#open OTU Table
#to read it into R, first remove the line"# Constructed from biom file" and remove the "#" from the entry "#OTU ID"
OTUtable<-read.delim(file="output/open_ref_otus/wf_bdiv_40000/otu_table_mc2_w_tax_even40000.biom.txt", sep="\t",header = TRUE)
#OTUtable<-read.delim(file="outputwithcontrols/wf_bdiv_40000/otu_table_mc2_w_tax_even40000.biomnoheader.txt", sep="\t",header = FALSE,check.names=FALSE)
OTUtablenonrare<-read.delim(file="output/open_ref_otus/otu_table_mc2_w_tax.biom.biomnoheader.txt", sep="\t",header = FALSE,check.names=FALSE)

#flip the OTU table to make it easier to work with
#note the last column is taxonomy
OTUtableflip<-t(as.matrix(OTUtable))
colnames(OTUtableflip) <- c(as.character(unlist(OTUtableflip[1,])))
names(OTUtableflip)<-colnames(OTUtableflip)
OTUtableflip = OTUtableflip[-1, ]
OTUtableflipdf<-as.data.frame(OTUtableflip)
#same for the table not rarefied for the poop data
OTUtableflipnonrare<-t(as.matrix(OTUtablenonrare))
colnames(OTUtableflipnonrare) <- c(as.character(unlist(OTUtableflipnonrare[1,])))
names(OTUtableflipnonrare)<-colnames(OTUtableflipnonrare)
OTUtableflipnonrare = OTUtableflipnonrare[-1, ]
OTUtableflipdfnonrare<-as.data.frame(OTUtableflipdfnonrare)

#set up the mapping file to add the OTU data
numOTUs<-2  #number of OTUs to add to the table
nextcolumn<-dim(mappingfile)[2]+1
lastcolumn<-nextcolumn+numOTUs-1
mappingfile[,nextcolumn:lastcolumn]<-NA
mappingfile["timepoint"]<-"NA"  #for the poop data

LRcol<-which(names(OTUtableflipdf)=="692154")
LJcol<-which(names(OTUtableflipdf)=="592160")
OTUnames<-c("692154","592160")
names(mappingfile)[nextcolumn:lastcolumn]<-as.character(OTUnames)

#make another identical table told hold the data rarefied at 10000
mappingfileforpoop<-mappingfile



#put the OTU counts of interest into the mappingfile for ease of use
for(sample in 1:dim(mappingfile)[1]) {  #go through the subset file
  samplename<-mappingfile[sample,1] #get the sample name
  samplerow<-which(OTUtableflipdf$`OTU ID` == as.character(samplename)) #get the row in the biomfile that has the sample
  if (!is.integer0(samplerow)) {
    OTUcounts<-c(as.numeric(as.character(OTUtableflipdf[samplerow,LRcol])),as.numeric(as.character(OTUtableflipdf[samplerow,LJcol])))
    mappingfile[sample,nextcolumn:lastcolumn]<-OTUcounts
  }
}

#for plotting the poop data
#use a nonrarefied OTU table for the poops
#same as before but also adds a column for plotting the poop samples, the poops were collected at times 0, 4,5,6,10 and 11 weeks
for(sample in 1:dim(mappingfileforpoop)[1]) {  #go through the subset file
  samplename<-mappingfileforpoop[sample,1] #get the sample name
  samplerow<-which(OTUtableflipdfnonrare$`OTU ID` == as.character(samplename)) #get the row in the biomfile that has the sample
  if (!is.integer0(samplerow)) {
    OTUcounts<-c(as.numeric(as.character(OTUtableflipdfnonrare[samplerow,LRcol])),as.numeric(as.character(OTUtableflipdfnonrare[samplerow,LJcol])))
    mappingfileforpoop[sample,nextcolumn:lastcolumn]<-OTUcounts
    if (mappingfileforpoop$SampleType[sample] == "Poop") { 
      if (as.numeric(as.character(mappingfileforpoop$SampleWeek[sample]))<4) {
        mappingfileforpoop[sample,]$timepoint<-0
      }
      else if (as.numeric(as.character(mappingfileforpoop$SampleWeek[sample]))>3 & as.numeric(mappingfileforpoop$SampleWeek[sample])<7) {
        mappingfileforpoop[sample,]$timepoint<-5
      }
      else if (as.numeric(as.character(mappingfileforpoop$SampleWeek[sample]))>8) {
        mappingfileforpoop[sample,]$timepoint<-10
      }
    }
  }
}


#subset AP3
AP3subset<-subset(mappingfile,mappingfile$FAMI=="FAMIAP3"& mappingfile$SequenceCounts>=40000) 
AP3subset<-droplevels(AP3subset)
AP3subsetsort<-arrange(AP3subset,xtfrm(Diet))
AP3subsetsort$Diet<-factor(AP3subsetsort$Diet, levels=c("LF", "HF"))

#Figure 4 – Figure supplement 3A
beeswarm(`692154`/40000 ~ Diet, data = AP3subsetsort, cex = 2, cex.lab=1.4, pch=21, col=c( "black","black"),bg =c( "green","orange"), ylab="Relative Abundance", main="OTU 692154")
bxplot(`692154`/40000 ~ Diet, data = AP3subsetsort, width = .5, add = TRUE)


#Figure 4 – Figure supplement 3B
beeswarm(`592160`/40000 ~ Diet, data = AP3subsetsort, cex = 2, cex.lab=1.4, pch=21, col=c( "black","black"),bg =c( "green","orange"), ylab="Relative Abundance", main="OTU 592160")
bxplot(`592160`/40000 ~ Diet, data = AP3subsetsort, width = .5, add = TRUE)


#subset AP1
AP1subset<-subset(mappingfile,mappingfile$FAMI=="FAMIAP1"& mappingfile$SequenceCounts>=40000) 
AP1subset<-droplevels(AP1subset)
AP1subsetsort<-arrange(AP1subset,xtfrm(Diet))
AP1subsetsort$Diet<-factor(AP1subsetsort$Diet, levels=c("LF", "HF"))

#Figure 4 – Figure supplement 3C
beeswarm(`592160`/40000 ~ Diet, data = AP1subsetsort, cex = 2, cex.lab=1.4, pch=21, col=c( "black","black"),bg =c( "green","orange"), ylab="Relative Abundance", main="OTU 592160")
bxplot(`592160`/40000 ~ Diet, data = AP1subsetsort, width = .5, add = TRUE)

#Figure 4 – Figure supplement 3D
beeswarm(`692154`/40000 ~ Diet, data = AP1subsetsort, cex = 2, cex.lab=1.4, pch=21, col=c( "black","black"),bg =c( "green","orange"), ylab="Relative Abundance", main="OTU 692154")
bxplot(`692154`/40000 ~ Diet, data = AP1subsetsort, width = .5, add = TRUE)

AP1poop10000<-subset(mappingfileforpoop,mappingfileforpoop$FAMI=="FAMIAP1_poop" &mappingfileforpoop$SequenceCounts>=10000)
AP1poop10000<-droplevels(AP1poop10000)
AP1poop10000["Diet_Timepoint"]<-paste(AP1poop10000$Diet, ".W",AP1poop10000$timepoint,sep="")
AP1poop10000sort<-arrange(AP1poop10000,xtfrm(Diet_Timepoint))
AP1poop10000sort$Diet_Timepoint<-factor(AP1poop10000sort$Diet_Timepoint, levels=c("LF.W0", "HF.W0","LF.W5","HF.W5","LF.W10","HF.W10"))

#Figure 4 – Figure supplement 3E
beeswarm(`692154`/SequenceCounts ~ Diet_Timepoint, data = AP1poop10000sort, cex = 2,spacing = .8, cex.lab=1.4, pch=21, col=c( "black","black"),bg =c( "green","orange"), ylab="Unrarified Relative Abundance", xlab="Diet.Timepoint in Weeks", main="OTU 692154")
bxplot(`692154`/SequenceCounts ~ Diet_Timepoint, data = AP1poop10000sort, width = .5, add = TRUE)


#Figures 4C and D and Figure 4 - Figure supplement 2
library("lme4")
library(ggplot2)
#subset AP4
AP4subset<-subset(mappingfile,mappingfile$FAMI=="FAMIAP4" & mappingfile$SequenceCounts>=40000)
AP4subset<-droplevels(AP4subset)
#add a column describing Diet and PMA status
AP4subset["DietPMA"]<-"NA"
for(samp in 1:dim(AP4subset)[1]) {
  if (AP4subset[samp,]$Diet == "LF" &  AP4subset[samp,]$PMA == "no" ) {
    AP4subset[samp,]$DietPMA <- "LFno"
  } else if (AP4subset[samp,]$Diet == "LF" &  AP4subset[samp,]$PMA == "yes" ) {
    AP4subset[samp,]$DietPMA <- "LFPMA"
  } else if (AP4subset[samp,]$Diet == "HF" &  AP4subset[samp,]$PMA == "no" ) {
    AP4subset[samp,]$DietPMA <- "HFno"
  } else if (AP4subset[samp,]$Diet == "HF" &  AP4subset[samp,]$PMA == "yes" ) {
    AP4subset[samp,]$DietPMA <- "HFPMA"
  }
}
AP4subsetsort<-arrange(AP4subset,xtfrm(DietPMA))
AP4subsetsort$DietPMA<-factor(AP4subsetsort$DietPMA, levels=c("LFno", "HFno","LFPMA","HFPMA"))
AP4subsetsortnoPMAsaline<-subset(AP4subsetsort,AP4subsetsort$PMA=="no" | (AP4subsetsort$PMA=="yes" & AP4subsetsort$Gavage=="18_2") )

AP4noPMA<-subset(AP4subsetsort,AP4subsetsort$PMA=="no")
AP4yesPMA<-subset(AP4subsetsort,AP4subsetsort$PMA=="yes")
AP4noPMA["Diet_Cage"]<-paste(AP4noPMA$Diet,".",AP4noPMA$Cage,sep="")
AP4noPMAsort<-arrange(AP4noPMA,xtfrm(Diet_Cage))
AP4noPMAsort$Diet_Cage<-factor(AP4noPMAsort$Diet_Cage, levels=c("LF.239556","LF.239557","LF.239560", "LF.239562","LF.239564","LF.239567","HF.239558","HF.239559", "HF.239561" ,"HF.239563","HF.239565","HF.239566"))


#Figures 4C
beeswarm(`692154`/40000 ~ DietPMA, data = AP4subsetsortnoPMAsaline, cex = 2, cex.lab=1.4, pch=21, col=c( "black","black"),bg =c( "green","orange"), ylab="Relative Abundance", main="OTU 692154", corral="wrap", spacing=.6)
bxplot(`692154`/40000 ~ DietPMA, data = AP4subsetsortnoPMAsaline, width = .5, add = TRUE)

kruskal.test(`692154`/40000 ~ DietPMA,data=subset(AP4subsetsortnoPMAsaline,AP4subsetsortnoPMAsaline$Diet=="LF"))
Kruskal-Wallis rank sum test

#data:  `692154`/40000 by DietPMA
#Kruskal-Wallis chi-squared = 1.7567, df = 1, p-value = 0.185

t.test(`692154`/40000 ~ DietPMA,data=subset(AP4subsetsortnoPMAsaline,AP4subsetsortnoPMAsaline$Diet=="LF"))
#Welch Two Sample t-test

#data:  `692154`/40000 by DietPMA
#t = -1.0231, df = 13.599, p-value = 0.3241
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.025599513  0.009094967
#sample estimates:
#  mean in group LFno mean in group LFPMA 
#0.007134091         0.015386364 

reutLF<-lmer(`692154` ~ DietPMA + (1|Cage), data=subset(AP4subsetsortnoPMAsaline,AP4subsetsortnoPMAsaline$Diet=="LF"), REML = FALSE)
reutLF.null<-lmer(`692154` ~ (1|Cage), data=subset(AP4subsetsortnoPMAsaline,AP4subsetsortnoPMAsaline$Diet=="LF"), REML = FALSE)
anova(reutLF,reutLF.null)
#Data: subset(AP4subsetsortnoPMAsaline, AP4subsetsortnoPMAsaline$Diet ==  ...
 #            Models:
  #             reutLF.null: `692154` ~ (1 | Cage)
  #           reutLF: `692154` ~ DietPMA + (1 | Cage)
#             Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
#             reutLF.null  3 517.44 521.93 -255.72   511.44                         
#             reutLF       4 517.00 522.98 -254.50   509.00 2.4425      1     0.1181
             
kruskal.test(`692154`/40000 ~ DietPMA,data=subset(AP4subsetsortnoPMAsaline,AP4subsetsortnoPMAsaline$Diet=="HF"))
# Kruskal-Wallis rank sum test
# 
# data:  `692154`/40000 by DietPMA
# Kruskal-Wallis chi-squared = 4.2471, df = 1, p-value = 0.03932

t.test(`692154`/40000 ~ DietPMA,data=subset(AP4subsetsortnoPMAsaline,AP4subsetsortnoPMAsaline$Diet=="HF"))
# Welch Two Sample t-test
# 
# data:  `692154`/40000 by DietPMA
# t = -2.9785, df = 9.8984, p-value = 0.01399
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.08541395 -0.01225147
# sample estimates:
#   mean in group HFno mean in group HFPMA 
# 0.01218229          0.06101500 

reutHF<-lmer(`692154` ~ DietPMA + (1|Cage), data=subset(AP4subsetsortnoPMAsaline,AP4subsetsortnoPMAsaline$Diet=="HF"), REML = FALSE)
reutHF.null<-lmer(`692154` ~ (1|Cage), data=subset(AP4subsetsortnoPMAsaline,AP4subsetsortnoPMAsaline$Diet=="HF"), REML = FALSE)
anova(reutHF,reutHF.null)
# Data: subset(AP4subsetsortnoPMAsaline, AP4subsetsortnoPMAsaline$Diet ==  ...
#              Models:
#                reutHF.null: `692154` ~ (1 | Cage)
#              reutHF: `692154` ~ DietPMA + (1 | Cage)
#              Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
#              reutHF.null  3 597.70 602.28 -295.85   591.70                             
#              reutHF       4 581.66 587.77 -286.83   573.66 18.036      1  2.167e-05 ***
#                ---
#                Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
             

reutmodel<-lmer(`692154` ~ Diet + (1|Cage), data = AP4noPMAsort, REML = FALSE)
reutmodel.null<-lmer(`692154` ~ (1|Cage), data = AP4noPMAsort, REML = FALSE)
anova(reutmodel.null,reutmodel)
# Data: AP4noPMAsort
# Models:
#   reutmodel.null: `692154` ~ (1 | Cage)
# reutmodel: `692154` ~ Diet + (1 | Cage)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# reutmodel.null  3 728.18 733.66 -361.09   722.18                         
# reutmodel       4 729.57 736.88 -360.78   721.57 0.6116      1     0.4342

reutmodel<-lmer(`692154` ~ Diet + (1|Cage), data = AP4yesPMA, REML = FALSE)
reutmodel.null<-lmer(`692154` ~ (1|Cage), data = AP4yesPMA, REML = FALSE)
anova(reutmodel.null,reutmodel)
# Data: AP4yesPMA
# Models:
#   reutmodel.null: `692154` ~ (1 | Cage)
# reutmodel: `692154` ~ Diet + (1 | Cage)
# Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)
# reutmodel.null  3 766.52 771.87 -380.26   760.52                        
# reutmodel       4 766.23 773.36 -379.11   758.23 2.292      1       0.13


#Figures 4D
beeswarm(`592160`/40000 ~ DietPMA, data = AP4subsetsortnoPMAsaline, cex = 2, cex.lab=1.4, pch=21, col=c( "black","black"),bg =c( "green","orange"), ylab="Relative Abundance", main="OTU 692154", corral="wrap", spacing=.6)
bxplot(`592160`/40000 ~ DietPMA, data = AP4subsetsortnoPMAsaline, width = .5, add = TRUE)

kruskal.test(`592160`/40000 ~ DietPMA,data=subset(AP4subsetsortnoPMAsaline,AP4subsetsortnoPMAsaline$Diet=="LF"))
# Kruskal-Wallis rank sum test
# 
# data:  `592160`/40000 by DietPMA
# Kruskal-Wallis chi-squared = 4.5737, df = 1, p-value = 0.03247

t.test(`592160`/40000 ~ DietPMA,data=subset(AP4subsetsortnoPMAsaline,AP4subsetsortnoPMAsaline$Diet=="LF"))
# Welch Two Sample t-test
# 
# data:  `592160`/40000 by DietPMA
# t = -1.7123, df = 15.086, p-value = 0.1073
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.36314711  0.03950848
# sample estimates:
#   mean in group LFno mean in group LFPMA 
# 0.1403148           0.3021341 

johnLF<-lmer(`592160` ~ DietPMA + (1|Cage), data=subset(AP4subsetsortnoPMAsaline,AP4subsetsortnoPMAsaline$Diet=="LF"), REML = FALSE)
johnLF.null<-lmer(`592160` ~ (1|Cage), data=subset(AP4subsetsortnoPMAsaline,AP4subsetsortnoPMAsaline$Diet=="LF"), REML = FALSE)
anova(johnLF,johnLF.null)
# Data: subset(AP4subsetsortnoPMAsaline, AP4subsetsortnoPMAsaline$Diet ==  ...
#              Models:
#                johnLF.null: `592160` ~ (1 | Cage)
#              johnLF: `592160` ~ DietPMA + (1 | Cage)
#              Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)  
#              johnLF.null  3 699.76 704.25 -346.88   693.76                           
#              johnLF       4 697.07 703.06 -344.53   689.07 4.6862      1     0.0304 *
#                ---
#                Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
             
kruskal.test(`592160`/40000 ~ DietPMA,data=subset(AP4subsetsortnoPMAsaline,AP4subsetsortnoPMAsaline$Diet=="HF"))
# Kruskal-Wallis rank sum test
# 
# data:  `592160`/40000 by DietPMA
# Kruskal-Wallis chi-squared = 12.893, df = 1, p-value = 0.0003298

t.test(`592160`/40000 ~ DietPMA,data=subset(AP4subsetsortnoPMAsaline,AP4subsetsortnoPMAsaline$Diet=="HF"))
# Welch Two Sample t-test
# 
# data:  `592160`/40000 by DietPMA
# t = -3.6613, df = 12.157, p-value = 0.003192
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.4508169 -0.1147377
# sample estimates:
#   mean in group HFno mean in group HFPMA 
# 0.1212802           0.4040575 

johnHF<-lmer(`592160` ~ DietPMA + (1|Cage), data=subset(AP4subsetsortnoPMAsaline,AP4subsetsortnoPMAsaline$Diet=="HF"), REML = FALSE)
johnHF.null<-lmer(`592160` ~ (1|Cage), data=subset(AP4subsetsortnoPMAsaline,AP4subsetsortnoPMAsaline$Diet=="HF"), REML = FALSE)
anova(johnHF,johnHF.null)
# Data: subset(AP4subsetsortnoPMAsaline, AP4subsetsortnoPMAsaline$Diet ==  ...
#              Models:
#                johnHF.null: `592160` ~ (1 | Cage)
#              johnHF: `592160` ~ DietPMA + (1 | Cage)
#              Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
#              johnHF.null  3 716.13 720.71 -355.07   710.13                             
#              johnHF       4 700.55 706.65 -346.27   692.55 17.585      1  2.748e-05 ***
#                ---
#                Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#              


johnmodel<-lmer(`592160` ~ Diet + (1|Cage), data = AP4noPMAsort, REML = FALSE)
johnmodel.null<-lmer(`592160` ~ (1|Cage), data = AP4noPMAsort, REML = FALSE)
anova(johnmodel.null,johnmodel)
# Data: AP4noPMAsort
# Models:
#   johnmodel.null: `592160` ~ (1 | Cage)
# johnmodel: `592160` ~ Diet + (1 | Cage)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# johnmodel.null  3 947.12 952.61 -470.56   941.12                         
# johnmodel       4 948.98 956.29 -470.49   940.98 0.1445      1     0.7039


johnmodel<-lmer(`592160` ~ Diet + (1|Cage), data = AP4yesPMA, REML = FALSE)
johnmodel.null<-lmer(`592160` ~ (1|Cage), data = AP4yesPMA, REML = FALSE)
anova(johnmodel.null,johnmodel)
# Data: AP4yesPMA
# Models:
#   reutmodel.null: `592160` ~ (1 | Cage)
# reutmodel: `592160` ~ Diet + (1 | Cage)
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# reutmodel.null  3 934.05 939.40 -464.02   928.05                         
# reutmodel       4 935.59 942.73 -463.80   927.59 0.4557      1     0.4996

#Figure 

#plot by cage
AP4noPMA["Diet_Cage"]<-paste(AP4noPMA$Diet,".",AP4noPMA$Cage,sep="")
AP4noPMAsort<-arrange(AP4noPMA,xtfrm(Diet_Cage))
AP4noPMAsort$Diet_Cage<-factor(AP4noPMAsort$Diet_Cage, levels=c("LF.239556","LF.239557","LF.239560", "LF.239562","LF.239564","LF.239567","HF.239558","HF.239559", "HF.239561" ,"HF.239563","HF.239565","HF.239566"))

#Figure 4 – Figure supplement 2C 

ggplot(data=AP4noPMAsort) +
  geom_boxplot(aes(x=Diet_Cage, y=`692154`/40000, fill=Diet)) +
  xlab("Cage") +
  ylab("Relative Abundance") +
  scale_colour_manual(values=c("green","orange")) + scale_fill_manual(values=c("green","orange")) +
  theme_classic() +
  theme(
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.x= element_text(size=20),
    axis.title.y= element_text(size=20),
    axis.text.x= element_text(size=5),
    axis.text.y= element_text(size=15)
  )
