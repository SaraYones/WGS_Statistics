
#Convert from xlsx to CSV
#convert("AML_WGS_Strelka_Filtered.xlsx", "AML_WGS_Strelka_Filtered.csv")
mutations=fread("20190507_AML_WGS_Strelka_proteincoding.csv",header=TRUE)
mutationsAll=fread("AML_WGS_Strelka_Filtered.csv",header=TRUE)
cases_A_P=read.xlsx("Adult_and_pediatric_cases.xlsx",header=TRUE,sheetName="Sheet1")
cases_A=cases_A_P$Adult.cases
cases_P=cases_A_P$Pediatriccases
mutations=as.data.frame(mutations)
mutations2=as.data.frame(mutations2)
mutationsAll=as.data.frame(mutationsAll)
meta_data=fread("AML_Filtered_Metadata.csv")
View(meta_data[order(meta_data$"Age at AML onset (Years)"),])
#Sort by Age and replace all occurences of _ with - in order to appear on titv

meta_data=meta_data[order(meta_data$"Age at AML onset (Years)"),]
meta_data$`Sample ID`=str_replace(meta_data$`Sample ID`, "_", "-")
meta_data=as.data.frame(meta_data)


#Some trials to change the values in a col but it worked with within function
#apply(as.data.frame(mutations[which(mutations[,"Variant_Classification"]=="frameshift_variant"),]),1, function(x) if(x["Variant_Type"]=="INS") {x["Variant_Classification"]="Frame_Shift_Ins"} else if(x["Variant_Type"]=="DEL") {x["Variant_Classification"]="Frame_Shift_Del"})
#hello=apply(mutations,1, function(x) if(x["Variant_Type"]=="INS") {x["Variant_Classification"]="Frame_Shift_Ins"} else if(x["Variant_Type"]=="DEL") {x["Variant_Classification"]="Frame_Shift_Del"})

#Changing a the value of frameshift_Variant in Variant_Classification col to Frame_Shift_Ins based on the Variant Type
mutations<-within(mutations, Variant_Classification[Variant_Type=='INS' & Variant_Classification=='frameshift_variant'] <- 'Frame_Shift_Ins')
mutations<-within(mutations, Variant_Classification[Variant_Type=='DEL' & Variant_Classification=='frameshift_variant'] <- 'Frame_Shift_Del')
#Comparison between Adults (D vs R)
mutations_Adults=mutationsAll[mutationsAll$Tumor_Sample_Barcode %in% cases_A,]
#Comparison between Pediatric (D vs R)
mutations_Pediatric=mutationsAll[mutationsAll$Tumor_Sample_Barcode %in% cases_P,]
#Comparison between Adults /Pediatric (D vs D)

#in this comparison i wont do the filteration here but in the part of the script that defines mutationDiagnosis and mutationRelapse
#mutations_A_P_Diagnosis=mutationsAll[mutationsAll$Tumor_Sample_Barcode %in% append(levels(cases_A_P$Adult.cases),levels(cases_A_P$Pediatriccases)),]
#mutations_A_P_Diagnosis=mutations_A_P_Diagnosis[mutations_A_P_Diagnosis$sample=="D",]

mutationsAlltemp=mutationsAll

mutationsAll=mutations_Pediatric
mutationsAll=mutations_Adults
#mutationsAll=mutations_A_P_Diagnosis

#return back to mutationsAll
mutationsAll=mutationsAlltemp
#Combine Names of tumor sample barcode to also include if its D or R to appear on titv 
mutations$Tumor_Sample_Barcode_temp=mutations$Tumor_Sample_Barcode
mutations$Tumor_Sample_Barcode=paste(mutations$Tumor_Sample_Barcode,"-",mutations$sample,sep="")
mutationsAll$Tumor_Sample_Barcode_temp=mutationsAll$Tumor_Sample_Barcode
mutationsAll$Tumor_Sample_Barcode=paste(mutationsAll$Tumor_Sample_Barcode,"-",mutationsAll$sample,sep="")

#From mutationAll we can divide into Diagnosis and Relapse
#mutationDiagnosis=mutationsAll[mutationsAll[,"sample"]=="D",]
mutationDiagnosisPC=mutations[mutations$Tumor_Sample_Barcode %in% Unpaired_LatestD,]
mutationRelapsePC=mutations[mutations$Tumor_Sample_Barcode %in% Unpaired_LatestR,]

#Purest Relapse Unmatched
mutationRelapse=mutationsAll[mutationsAll$Tumor_Sample_Barcode %in% purest_R$PuresetRelapse,]
mutationDiagnosis=mutationsAll[mutationsAll$Tumor_Sample_Barcode %in% purest_R$DiagnosisCases,]

#Purest Relapse Matched
mutationsAll$Variant_Classification=rep("Missense_Mutation",dim(mutationsAll)[1])
mutationRelapse=mutationsAll[mutationsAll$Tumor_Sample_Barcode %in% purest_R_matched$PuresetRelapse,]
mutationDiagnosis=mutationsAll[mutationsAll$Tumor_Sample_Barcode %in% purest_R_matched$DiagnosisCases,]
mutationRelapse=filterVariants(mutationDiagnosis,mutationRelapse)
#Comparison (D VS D)
#mutationsAll=mutationsAll[mutationsAll$sample=="D",]
mutationsAll=mutationsAll[mutationsAll$Tumor_Sample_Barcode %in% purest_R$DiagnosisCases,]
#Comparison (R VS R)
mutationsAll=mutationsAll[mutationsAll$Tumor_Sample_Barcode %in% purest_R$PuresetRelapse,]
#mutationsAll=mutationsAll[mutationsAll$sample!="D",]
#Comparison (D VS D)
mutationDiagnosisAdult=mutationsAll[mutationsAll$Tumor_Sample_Barcode_temp %in% cases_A ,]
mutationDiagnosisPediatric=mutationsAll[mutationsAll$Tumor_Sample_Barcode_temp %in% cases_P ,]
#mutationDiagnosisPediatric=filterVariants(mutationDiagnosisAdult,mutationDiagnosisPediatric)

#I have to have mutationDiagnosisAdult and mutationDiagnosisPediatric in order to get mutationRelapseAdult and Pediatric
#Tumor_Sample_Barcode_temp is really important 
#Comparison (R VS R)
mutationRelapseAdult=mutationsAll[mutationsAll$Tumor_Sample_Barcode_temp %in% cases_A ,]
mutationRelapseAdult=filterVariants(mutationDiagnosisAdult,mutationRelapseAdult)
mutationRelapsePediatric=mutationsAll[mutationsAll$Tumor_Sample_Barcode_temp %in% cases_P ,]
mutationRelapsePediatric=filterVariants(mutationDiagnosisPediatric,mutationRelapsePediatric)


filterVariants=function(mutationDiagnosis,mutationRelapse)
{
  
cases_D=mutationDiagnosis$Tumor_Sample_Barcode
cases_D=unlist(strsplit(cases_D,"-D"))
mutationDiagnosis=as.data.frame(cbind(mutationDiagnosis,cases_D))
cases_R=mutationRelapse$Tumor_Sample_Barcode
cases_R=gsub("(AML[0-9][0-9][0-9])-.*","\\1",mutationRelapse$Tumor_Sample_Barcode)
mutationRelapse=as.data.frame(cbind(mutationRelapse,cases_R))
mutationRelapse=as.data.table(mutationRelapse)
#remove variants that are in the diagnosis cases (this function is in data.table package)
setkey(mutationRelapse,cases_R,Chromosome,Hugo_Symbol,Start_Position,End_Position,Reference_Allele,Tumor_Seq_Allele2)

index=mutationRelapse[.(mutationDiagnosis$cases_D,mutationDiagnosis$Chromosome,mutationDiagnosis$Hugo_Symbol,mutationDiagnosis$Start_Position,mutationDiagnosis$End_Position,mutationDiagnosis$Reference_Allele,mutationDiagnosis$Tumor_Seq_Allele2),which=TRUE]
notNA=!is.na(mutationRelapse[.(mutationDiagnosis$cases_D,mutationDiagnosis$Chromosome,mutationDiagnosis$Hugo_Symbol,mutationDiagnosis$Start_Position,mutationDiagnosis$End_Position,mutationDiagnosis$Reference_Allele,mutationDiagnosis$Tumor_Seq_Allele2),which=TRUE])
index=unique(index[which(notNA)])
mutationRelapse=mutationRelapse[-index,]

return(mutationRelapse)
}

calculateTiTvPvalues=function(TiTvMatrix,flagPaired,orderNames)
{
  TiTv=unique(as.character(TiTvMatrix$variable))
  names(TiTv)=TiTv
  TiTv=TiTv[order(factor(names(TiTv), levels = orderNames))]
  classes=unique(as.character(TiTvMatrix$classType))
  pvalues=NULL  
  for (i in 1:length(TiTv))
  {
  temp=TiTvMatrix[TiTvMatrix$variable==TiTv[i],]
  #Normality test
  result1=shapiro.test(temp[temp$classType==classes[1],"value"])
  result2=shapiro.test(temp[temp$classType==classes[2],"value"])
  #t-test
  if(result1$p.value>0.05 && result2$p.value>0.05)
  { print(paste("using t-test for",TiTv[i]))
    #print(resultDiagnosis)
    if(flagPaired=="Paired")
    {
      result=compare_means(value ~ classType, data = temp, paired = TRUE,method = "t.test")
      pvalues=append(pvalues,result$p.adj)
    }else{
      result=compare_means(value ~ classType, data = temp, paired = FALSE,method = "t.test")
      pvalues=append(pvalues,result$p.adj)
      
    }
  }else{#Wilxcon any of them is significantly differnt than normal 
    print(paste("using Wilcoxon-test for",TiTv[i]))
    if(flagPaired=="Paired")
    {
      result=compare_means(value ~ classType, data = temp, paired = TRUE,method = "wilcox.test")
      pvalues=append(pvalues,result$p.adj)
      
    }else{

      result=compare_means(value ~ classType, data = temp, paired = FALSE,method = "wilcox.test")
      pvalues=append(pvalues,result$p.adj)
    }
    
  }
  }
  names(pvalues)<-TiTv
  return(pvalues)
}

laml = read.maf(maf = mutations,useAll = TRUE)
lamlAll = read.maf(maf = mutationsAll,useAll = TRUE)
lamlDiagnosis=read.maf(maf = mutationDiagnosis,useAll = TRUE)
lamlRelapse=read.maf(maf = mutationRelapse,useAll = TRUE)
#Comparison (D VS D)
lamlDiagnosisAdults=read.maf(maf = mutationDiagnosisAdult,useAll = TRUE)
lamlDiagnosisPediatric=read.maf(maf = mutationDiagnosisPediatric,useAll = TRUE)

#Comparison (R VS R)
lamlRelapseAdults=read.maf(maf = mutationRelapseAdult,useAll = TRUE)
lamlRelapsePediatric=read.maf(maf = mutationRelapsePediatric,useAll = TRUE)



lamlDiagnosisPC=read.maf(maf = mutationDiagnosisPC,useAll = TRUE)
lamlRelapsePC=read.maf(maf = mutationRelapsePC,useAll = TRUE)

laml2=read.maf(maf=mutations2)
temp= gsub("(AML[0-9][0-9][0-9])(-.*)","\\1",gsub("(AML[0-9][0-9][0-9])(_.*)","\\1", meta_data$`Sample ID`))
match(temp,mutations$Tumor_Sample_Barcode)
#If i want to change the name of a col
#names(metadata)[names(metadata)=="patient"]="Tumor_Sample_Barcode"

getSampleSummary(laml)
pdf( "SummaryMAF.pdf", onefile=TRUE)
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
laml.titv2 = titv(maf = laml2, plot = FALSE, useSyn = TRUE)
laml.titvAll = titv(maf = lamlAll, plot = FALSE, useSyn = TRUE)
laml.titvAll = titv(maf = lamlAll, plot = FALSE, useSyn = TRUE)
#Comparisons(D V R)
laml.titvDiagnosis=titv(maf = lamlDiagnosis, plot = FALSE, useSyn = TRUE)
laml.titvRelapse=titv(maf = lamlRelapse, plot = FALSE, useSyn = TRUE)
#Comparison (D vs D)
laml.titvDiagnosisA=titv(maf = lamlDiagnosisAdults, plot = FALSE, useSyn = TRUE)
laml.titvDiagnosisP=titv(maf = lamlDiagnosisPediatric, plot = FALSE, useSyn = TRUE)
#Comparison (R vs R)
laml.titvRelapseA=titv(maf = lamlRelapseAdults, plot = FALSE, useSyn = TRUE)
laml.titvRelapseP=titv(maf = lamlRelapsePediatric, plot = FALSE, useSyn = TRUE)


#plot titv summary

pdf( "TiVSTvAll.pdf", onefile=TRUE)
plotTiTv(res = laml.titv,showBarcodes=TRUE,sampleOrder=meta_data$`Sample ID`)
col = c('coral4', 'lightcyan4', 'cornflowerblue', 'lightsalmon1', 'yellow', 'deeppink3')
names(col) = c('C>T', 'C>G', 'C>A', 'T>A', 'T>C','T>G')
plotTiTv(res = laml.titvAll,showBarcodes=TRUE,sampleOrder=meta_data$`Sample ID`,color = col)
plotTiTv(res = laml.titv2,showBarcodes=TRUE)
dev.off()
#-------------------------------------------------------------------------------------------

pdf( "TiVSTvDiagnosis.PDF", onefile=TRUE)
plotTiTv(res = laml.titvDiagnosis,showBarcodes=TRUE,sampleOrder=meta_dataDiagnosis$`Sample ID`)
col = c('coral4', 'lightcyan4', 'cornflowerblue', 'lightsalmon1', 'yellow', 'deeppink3')
names(col) = c('C>T', 'C>G', 'C>A', 'T>A', 'T>C','T>G')

meta_dataDiagnosis=meta_data[meta_data[,"Stage"]=="D",]
plotTiTv(res = laml.titvDiagnosis,showBarcodes=TRUE,sampleOrder=meta_dataDiagnosis$`Sample ID`,color = col)
dev.off()


#--------------------------------------------------------------------------------------------

pdf( "TiVSTvRelapse.PDF", onefile=TRUE)
#plotTiTv(res = laml.titvRelapse,showBarcodes=TRUE,sampleOrder=meta_dataRelapse$`Sample ID`)
col = c('coral4', 'lightcyan4', 'cornflowerblue', 'lightsalmon1', 'yellow', 'deeppink3')
names(col) = c('C>T', 'C>G', 'C>A', 'T>A', 'T>C','T>G')

meta_dataRelapse=meta_data[meta_data[,"Sample ID"] %in% purest_R$PuresetRelapse,]
plotTiTv(res = laml.titvRelapse,showBarcodes=TRUE,sampleOrder=meta_dataRelapse$`Sample ID`,color = col)
dev.off()

#--------------------------------------------------------------------------------------------
laml.titvlist=list(laml.titvDiagnosis,laml.titvRelapse,laml.titvDiagnosis,laml.titvRelapse)
plotType=list("box","box","bar","bar")
meta_datalist=list(meta_dataDiagnosis,meta_dataRelapse,meta_dataDiagnosis,meta_dataRelapse)
#---------------------------------------------------------------------------------------------

myplots <- list()  # new empty list
for (i in 1:length(laml.titvlist))
{
  local({
    i <- i
    # theme_set(theme_cowplot(font_size=4)) # reduce default font size
    # The first is for plotting 1 color group for all meta data variables
    # p1 <- qplot(ordered_metadata_for_exploratory[which(ordered_metadata_for_exploratory[,"clusters"]=="red"),meta[i]], as.factor(rownames(ordered_metadata_for_exploratory)[which(ordered_metadata_for_exploratory[,"clusters"]=="red")]),colour = I("red"),xlab = meta[i], ylab = "samples")+theme(plot.background=element_rect(fill="white", colour=NA))#+theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
    #The second for plotting all color group for just 1 meta data variable at a time
    print("hello")
    if(plotType[i]=="box")
    {print("hi")
      #  p1 <- plotTiTv(res = laml.titvlist[[i]],plotType=as.character(plotType[i]),color = col)
      plotTiTv(res = laml.titvlist[[i]],plotType=as.character(plotType[i]),color = col,textSize=0.5,baseFontSize =0.5,axisTextSize =c(0.5,0.5))
      p1<-recordPlot()
      
      print("I am here")
    }else{
      plotTiTv(res = laml.titvlist[[i]],plotType=as.character(plotType[i]),showBarcodes=TRUE,sampleOrder=meta_datalist[[i]]$`Sample ID`,color = col,textSize=0.5,baseFontSize =0.5,axisTextSize =c(0.5,0.5))
      p1<-recordPlot()
      #p1 <-plotTiTv(res = laml.titvlist[[i]],showBarcodes=TRUE,sampleOrder=meta_datalist[[i]]$`Sample ID`,color = col)
    }
    print(i)
    print(p1)
    myplots[[i]] <<- p1  # add each plot into plot list
  })
}
#multiplot(plotlist = myplots ,ncol=2,nrow=6)
# ggarrange(plotlist =  myplots ,ncol=2,nrow=5)
# pdf(paste(getwd(),"red",".pdf",sep=""), onefile=TRUE)

# The solution to the overlaying problem that used to take place in multiple grid was setting the base_aspect_ratio


plot2by2=plot_grid(plotlist = myplots ,labels="auto",nrow=2,ncol=2, scale = 0.9) #label_size=1,align="l")
plot_grid(plotlist = myplots,align="tblr", nrow = 4,ncol=6, rel_widths = c(4,4),rel_heights = c(3,3))


save_plot("TiTVGrid.png", plot2by2)
#,base_aspect_ratio=1.5
#  ncol = 2, # we're saving a grid plot of 2 columns
#  nrow = 2, # and 2 rows
# each individual subplot should have an aspect ratio of 1.3
# base_aspect_ratio = 1.3



library("gridExtra")


ggarrange(myplots[[1]],myplots[[2]] , 
          labels = c("A", "B"),
          ncol = 2, nrow = 2)

#This is the main script that we run for transitions and transversions
#-----------Plot Diagnosis and Relapse Transitions and transversion infront of one another  
pdf( "AllPatients/TiTV/TiVSTvBoth.pdf", onefile=TRUE)
classType=append(rep("D",dim(as.data.frame(laml.titvDiagnosis))[1]),rep("R",dim(as.data.frame(laml.titvRelapse))[1]))
figureTitv=rbind(as.data.frame(laml.titvDiagnosis),as.data.frame(laml.titvRelapse))
#Comparison (D vs D)
classType=append(rep("D-Adults",dim(as.data.frame(laml.titvDiagnosisA))[1]),rep("D-Pediatric",dim(as.data.frame(laml.titvDiagnosisP))[1]))
figureTitv=rbind(as.data.frame(laml.titvDiagnosisA),as.data.frame(laml.titvDiagnosisP))
#Comparison(R vs R)
classType=append(rep("R-Adults",dim(as.data.frame(laml.titvRelapseA))[1]),rep("R-Pediatric",dim(as.data.frame(laml.titvRelapseP))[1]))
figureTitv=rbind(as.data.frame(laml.titvRelapseA),as.data.frame(laml.titvRelapseP))
#Comparisons (D vs R)



figureTitv=as.data.frame(figureTitv)
figureTitv$classType=classType
figureTitv=figureTitv[,c("classType","fraction.contribution.C.A","fraction.contribution.C.G","fraction.contribution.C.T","fraction.contribution.T.C","fraction.contribution.T.A","fraction.contribution.T.G")]
colnames(figureTitv)<-c("classType","C:G>A:T","C:G>G:C","C:G>T:A","T:A>C:G", "T:A>A:T", "T:A>G:C")
df.m <- melt(figureTitv, id.var = "classType")
tgc <- summarySE(df.m, measurevar= "value", groupvars=c("classType","variable"))

tgc2 <- tgc
tgc2$variable <- factor(tgc$variable)
ggplot(tgc2, aes(x=variable, y=value, fill=classType)) + 
  geom_boxplot(position=position_dodge(), stat="identity",ymin=value-se, ymax=value+se) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))

require(ggplot2)
#First Trial

# p1<-ggplot(data = df.m, aes(x=variable, y=value))+geom_boxplot(aes(fill=classType),alpha=0.3)+geom_jitter(aes(fill=classType),position=position_jitter(0.2))+theme(axis.text=element_text(size=9),legend.title=element_text(size=10), 
#                                                                                          legend.text=element_text(size=9))+theme(legend.title = element_blank())+labs(x = "")+labs(y = "")#scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
#Trials

#Good 
#  p1<-ggplot(data = df.m, aes(x=variable, y=value))+geom_boxplot(aes(fill=classType),alpha=0.3)+geom_point(aes(y=value, group=classType,shape=as.factor(classType)),alpha=1,size=1.5, position = position_dodge(width=0.75)) +scale_fill_grey() +theme(axis.text=element_text(size=9),legend.title=element_text(size=10),                                                                                                                                                                                                                                                     legend.text=element_text(size=9))+theme(legend.title = element_blank())+labs(x = "")+labs(y = "")


#Reordered descendingly based on medians and make the fill color grey scale and make the shape of the jitter different

#Best Fit   

#To add whisker bars we use stat_boxplot and make size so small to see the dots from behind   
#Format p-values to include only 3 significant figures.
p1<-ggplot(data = df.m,aes(x = fct_reorder(variable, value,.desc=TRUE), y = value,fill=classType))+stat_boxplot( geom='errorbar', linetype=1,size =0.1, width=0.5,position = position_dodge(width=0.75))+theme_bw()+ geom_point(aes(y=value, group=classType),alpha=1,size=0.2, position = position_dodge(width=0.75))+scale_y_continuous(breaks = seq(0, 60, by=15), limits=c(0,60))
p1<-p1+geom_boxplot(inherit.aes = TRUE,aes(fill=classType),alpha=0.3,outlier.size=0,lwd=0.1,stat = "boxplot")+
  scale_fill_grey() +theme(axis.text=element_text(size=8), axis.title=element_text(size=8),legend.title=element_text(size=10),legend.text=element_text(size=9))+theme(legend.title = element_blank())+labs(x = "")+labs(y = "Frequency of substitutions")+stat_compare_means(aes(group = classType,label=sprintf("p = %5.4f", as.numeric(..p.format..))),label.y = 105,paired = FALSE ,label.x.npc="middle",size = 2)

                   


#dont format the pvalues 
p1<-ggplot(data = df.m,aes(x = fct_reorder(variable, value,.desc=TRUE), y = value,fill=classType))+stat_boxplot( geom='errorbar', linetype=1,size =0.1, width=0.5,position = position_dodge(width=0.75))+theme_bw()+ geom_point(aes(y=value, group=classType,legend=FALSE),alpha=1,size=0.2, position = position_dodge(width=0.75))+scale_y_continuous(breaks = seq(0, 60, by=15), limits=c(0,60))
p1<-p1+geom_boxplot(inherit.aes = TRUE,aes(fill=classType),alpha=0.3,outlier.size=0,lwd=0.1,stat = "boxplot")+scale_fill_grey() +theme(axis.text=element_text(size=8), axis.title=element_text(size=12),legend.title=element_text(size=10), 
                           legend.text=element_text(size=9))+theme(legend.title = element_blank())+labs(x = "")+labs(y = "Frequency of substitutions")+stat_compare_means(aes(group = classType),label.y = 105,paired=FALSE,label.x.npc="middle",size = 2,  label = "p.format")

#p1<-recordPlot()
#Using the package of the box plots
#p1<-ggplot2.boxplot(data=df.m, xName='variable',yName='value', groupName='classType', 
#                   position=position_dodge(0.8),
#                  backgroundColor="white", groupColors=c('#999999','#E69F00'),
#                 addDot=TRUE, dotSize=0.01,legendPosition="top") 
# p1<- ggplot(data = df.m, aes(x=variable, y=value))+geom_boxplot(aes(fill=classType),alpha=0.3)+geom_jitter(aes(y=value, group=classType),size=0.1) +scale_fill_grey() +theme(axis.text=element_text(size=9),legend.title=element_text(size=10), 
#                                                                                                                                                                                                                                                     legend.text=element_text(size=9))+theme(legend.title = element_blank())+labs(x = "")+labs(y = "")

#------------------------------------


#p1<-ggplot(data = df.m, aes(x=variable, y=value))+geom_point(aes(y=value, group=classType,colour=as.factor(classType)),size=5, position = position_dodge(width=0.75)) + geom_boxplot(aes(fill=classType),alpha=0.3,coef=1)+scale_colour_grey()+theme(axis.text=element_text(size=9),legend.title=element_text(size=10), 
#                                                                                                                                                                                                  legend.text=element_text(size=9))+theme(legend.title = element_blank())+labs(x = "")+labs(y = "")




#p1<-ggplot(data = df.m, aes(x=variable, y=value))+geom_point(aes(y=value, group=classType,colour=as.factor(classType)),size=0.1, position = position_dodge(width=0.75)) + scale_colour_grey()+theme(axis.text=element_text(size=9),legend.title=element_text(size=10), 
#                                                                                                                                                                                                                                                legend.text=element_text(size=9))+theme(legend.title = element_blank())+labs(x = "")+labs(y = "")




#---------------------------------------------------------------------------------------------


#Previous
p1<-p1+theme(legend.title = element_blank())+labs(x = "")+labs(y = "")+scale_fill_grey()+stat_compare_means(aes(group = classType), label = "p.format",size = 2) 

# +scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+stat_compare_means(aes(group = classType), label = "p.format",size = 2)


#------------Plot Fractions-------------------------------------------------

figureTitvFraction=rbind(as.data.frame(laml.titvDiagnosis)[,c("TiTv.fractions.Ti","TiTv.fractions.Tv")],as.data.frame(laml.titvRelapse)[,c("TiTv.fractions.Ti","TiTv.fractions.Tv")])
#Comparison (D vs D)
figureTitvFraction=rbind(as.data.frame(laml.titvDiagnosisA)[,c("TiTv.fractions.Ti","TiTv.fractions.Tv")],as.data.frame(laml.titvDiagnosisP)[,c("TiTv.fractions.Ti","TiTv.fractions.Tv")])
#Comparison (R vs R)
figureTitvFraction=rbind(as.data.frame(laml.titvRelapseA)[,c("TiTv.fractions.Ti","TiTv.fractions.Tv")],as.data.frame(laml.titvRelapseP)[,c("TiTv.fractions.Ti","TiTv.fractions.Tv")])



figureTitvFraction$classType=classType
colnames(figureTitvFraction)<-c("Ti","Tv", "classType")
df.mFraction <- melt(figureTitvFraction, id.var = "classType")

shapiro.test(my_data$len)
require(ggplot2)
# p2<-ggplot(data = df.mFraction, aes(x=variable, y=value))+geom_boxplot(aes(fill=classType),outlier.size=0) +theme(axis.text=element_text(size=9),legend.title=element_text(size=10), 
#                                                                                          legend.text=element_text(size=9))+theme(legend.title = element_blank())+labs(x = "")+labs(y = "")+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

# p2+theme(legend.title = element_blank())+labs(x = "")+labs(y = "")+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

# p2<-p2+stat_compare_means(aes(group = classType), label = "p.format",size = 2)


p2<-ggplot(data =  df.mFraction,aes(x=variable, y=value,fill=classType))+stat_boxplot( geom='errorbar', linetype=1,size =0.1, width=0.5,position = position_dodge(width=0.75))+scale_y_continuous(breaks = seq(0, 80, by=15), limits=c(0,80))+geom_boxplot(aes(fill=classType),alpha=0.3,outlier.size = 0,lwd=0.1)+theme_bw() +geom_point(aes(y=value, group=classType),alpha=1,size=0.2, position = position_dodge(width=0.75)) +scale_fill_grey() +theme(axis.text=element_text(size=9), axis.title=element_text(size=12),legend.title=element_text(size=10), 
                                                                                                                                                                                                                                                                                                                                                                                            legend.text=element_text(size=9))+theme(legend.title = element_blank())+labs(x = "Relapse/adult vs pediatric")+labs(y = "Transversion frequency")+stat_compare_means(aes(group = classType), label = "p.format",label.y=85,label.x.npc="middle",size = 2)
#Show only transversions
df.mFraction=df.mFraction[df.mFraction$variable=="Tv",]
df.mFraction$variable=as.character(df.mFraction$variable)
#Change it according to the plot
df.mFraction=df.mFraction %>% 
  mutate(variable= replace(variable, variable  == 'Tv', "Relapse/adult vs pediatric"))

df.mFraction=df.mFraction %>% 
  mutate(variable= replace(variable, variable  == 'Relapse/adult vs pediatric', "Tv"))

#FacetWrap helped me implement the last suggestion of Svea to defina the box plot for adult and pediatric for the transversions part

p2<-ggplot(data =  df.mFraction,aes(x=variable, y=value,fill=classType))+stat_boxplot( geom='errorbar', linetype=1,size =0.1, width=0.5,position = position_dodge(width=0.75))+scale_y_continuous(breaks = seq(0, 80, by=15), limits=c(0,80))+geom_boxplot(aes(fill=classType),alpha=0.3,outlier.size = 0,lwd=0.1)+theme_bw() +geom_point(aes(y=value, group=classType),alpha=1,size=0.2, position = position_dodge(width=0.75)) +scale_fill_grey() +theme(axis.text=element_text(size=8), axis.title=element_text(size=8),legend.title=element_text(size=10), 
                                                                                                                                                                                                                                                                                                                                                                                            legend.text=element_text(size=9))+theme(legend.title = element_blank())+labs(x = "")+labs(y = "Transversion frequency")+facet_wrap(~ classType)+stat_compare_means(aes(group = classType), label = "p.format",label.y=85,label.x.npc="middle",size = 2)


#p2<-recordPlot()
dev.off()

#In case i use tiff it should be called first and then plot_grid is called without assigning it to a variable then afterwards dev.off()
#If i want to use save_plot then dev off first then pcombined then save plot 
tiff(paste("AllPatients/TiTV/TiVSTvBothgrey_DRAdultTrial1","-",Sys.Date(),".tiff"), units="in", width=10, height=5, res=400)
pcombined<-plot_grid( as.grob(as.ggplot(p1)), as.grob(as.ggplot(p2)), labels="AUTO",label_size = 10)

#dev.off()
save_plot(paste("AllPatients/TiTV/TiVSTvBothgrey","-",Sys.Date(),".pdf",sep=""),pcombined, ncol = 2,nrow=1,base_aspect_ratio=1.1)
#For only transversions plot without aspect ratio

#bitmap(paste("AllPatients/TiTV/TiVSTvBothgrey_DiagnosisTrial2","-",Sys.Date(),".tiff"), height = 4, width = 4, units = 'in', type="tifflzw", res=300)
pcombined <- plot_grid( as.grob(as.ggplot(p1)), as.grob(as.ggplot(p2)), labels="AUTO",label_size = 10)

dev.off()


save_plot(paste("AllPatients/TiTV/TiVSTvBothgrey_DiagnosisTrial","-",Sys.Date(),".pdf",sep=""),pcombined ,resolution(600))


#  ggarrange(p1,p2 , 
#           labels = c("A", "B"),
#          ncol = 2, nrow = 1)

# save_plot("AllPatients/TiTV/TiVSTvBoth.pdf", plot2by2,base_aspect_ratio=1.5)
#  ncol = 2, # we're saving a grid plot of 2 columns

# dev.off()
