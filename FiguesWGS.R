#FiguresWGS
library(maftools)
library(xlsx)
library("rio")
library(data.table)
library(stringr)
library(cowplot)
library(ggplot2)
library("RColorBrewer")
#For ggarrange
library("gridExtra")
library("ggpubr")
library(forcats)
library("easyGgplot2")
library(xlsx)
#Used with frequency and difference between cohorts plot
library(grid)
library(gridExtra)
#For legends
library(lemon)

#Convert from xlsx to CSV
#convert("AML_WGS_Strelka_Filtered.xlsx", "AML_WGS_Strelka_Filtered.csv")
mutations=fread("20190507_AML_WGS_Strelka_proteincoding.csv",header=TRUE)
mutationsAll=fread("AML_WGS_Strelka_Filtered.csv",header=TRUE)
mutations2=fread("genome.wustl.edu_LAML.IlluminaGA_DNASeq.Level_2.2.13.0.somatic.maf",header=TRUE)
purest_R=fread("PurestR.csv",header=TRUE)
Unpaired=fread("MEAndCooccurence.csv",header=TRUE)
Unpaired_LatestD=Unpaired$`Unpaired: All D (WGS and WES)`
Unpaired_LatestR=Unpaired$`Unpaired: All latest  R (WGS and WES)`
recurrent_genes=c("ARID1A","ASXL1","ASXL2","BCOR","BCORL1","CEBPA","CSF3R","DNMT3A","EZH2","FAM5C","FLT3","GATA1","GATA2","H3F3A","HNRNPK","IDH1","IDH2","IKZF1","KIT","KMT2A","KRAS","MGA","MYC","NF1","NPM1","NRAS","NRNX3","PHF6","PTPN11","RAD21","RUNX1","SETD2","SMC1A","SMC3","SRSF1","SRSF2","STAG2","TET2","TP53","U2AF1","UBTF","WT1")


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

#Combine Names of tumor sample barcode to also include if its D or R to appear on titv 
mutations$Tumor_Sample_Barcode_temp=mutations$Tumor_Sample_Barcode
mutations$Tumor_Sample_Barcode=paste(mutations$Tumor_Sample_Barcode,"-",mutations$sample,sep="")
mutationsAll$Tumor_Sample_Barcode_temp=mutationsAll$Tumor_Sample_Barcode
mutationsAll$Tumor_Sample_Barcode=paste(mutationsAll$Tumor_Sample_Barcode,"-",mutationsAll$sample,sep="")

#From mutationAll we can divide into Diagnosis and Relapse
#mutationDiagnosis=mutationsAll[mutationsAll[,"sample"]=="D",]
mutationDiagnosisPC=mutations[mutations$Tumor_Sample_Barcode %in% Unpaired_LatestD,]
mutationRelapsePC=mutations[mutations$Tumor_Sample_Barcode %in% Unpaired_LatestR,]
mutationRelapse=mutationsAll[mutationsAll$Tumor_Sample_Barcode %in% purest_R$PuresetRelapse,]
mutationDiagnosis=mutationsAll[mutationsAll$Tumor_Sample_Barcode %in% purest_R$DiagnosisCases,]


laml = read.maf(maf = mutations,useAll = TRUE)
mutationsAll$Variant_Classification=rep("Missense_Mutation",dim(mutationsAll)[1])
lamlAll = read.maf(maf = mutationsAll,useAll = TRUE)
lamlDiagnosis=read.maf(maf = mutationDiagnosis,useAll = TRUE)
lamlRelapse=read.maf(maf = mutationRelapse,useAll = TRUE)
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
laml.titvDiagnosis=titv(maf = lamlDiagnosis, plot = FALSE, useSyn = TRUE)
laml.titvRelapse=titv(maf = lamlRelapse, plot = FALSE, useSyn = TRUE)
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
   
  #-----------Plot Diagnosis and Relapse Transitions and transversion infront of one another  
   pdf( "AllPatients/TiTV/TiVSTvBoth.pdf", onefile=TRUE)
   classType=append(rep("D",dim(as.data.frame(laml.titvDiagnosis))[1]),rep("R",dim(as.data.frame(laml.titvRelapse))[1]))
   figureTitv=rbind(as.data.frame(laml.titvDiagnosis),as.data.frame(laml.titvRelapse))
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
p1<-ggplot(data = df.m,aes(x = fct_reorder(variable, value,.desc=TRUE), y = value,fill=classType))+stat_boxplot( geom='errorbar', linetype=1,size =0.1, width=0.5,position = position_dodge(width=0.75))+theme_bw()+ geom_point(aes(y=value, group=classType,legend=FALSE),alpha=1,size=0.2, position = position_dodge(width=0.75))
p1<-p1+geom_boxplot(inherit.aes = TRUE,aes(fill=classType),alpha=0.3,outlier.size=0,lwd=0.1,stat = "boxplot")+
   scale_fill_grey() +theme(axis.text=element_text(size=8), axis.title=element_text(size=12),legend.title=element_text(size=10), 
 legend.text=element_text(size=9))+theme(legend.title = element_blank())+labs(x = "")+labs(y = "% mutation")+stat_compare_means(aes(group = classType),label.y = 85,label.x.npc="middle",size = 2,  label = "p.format")
   
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
   figureTitvFraction$classType=classType
   colnames(figureTitvFraction)<-c("Ti","Tv", "classType")
   df.mFraction <- melt(figureTitvFraction, id.var = "classType")
   
   shapiro.test(my_data$len)
   require(ggplot2)
  # p2<-ggplot(data = df.mFraction, aes(x=variable, y=value))+geom_boxplot(aes(fill=classType),outlier.size=0) +theme(axis.text=element_text(size=9),legend.title=element_text(size=10), 
   #                                                                                          legend.text=element_text(size=9))+theme(legend.title = element_blank())+labs(x = "")+labs(y = "")+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
   
  # p2+theme(legend.title = element_blank())+labs(x = "")+labs(y = "")+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
   
  # p2<-p2+stat_compare_means(aes(group = classType), label = "p.format",size = 2)
   
   
   p2<-ggplot(data =  df.mFraction,aes(x=variable, y=value,fill=classType))+stat_boxplot( geom='errorbar', linetype=1,size =0.1, width=0.5,position = position_dodge(width=0.75))+geom_boxplot(aes(fill=classType),alpha=0.3,outlier.size = 0,lwd=0.1)+theme_bw() +geom_point(aes(y=value, group=classType),alpha=1,size=0.2, position = position_dodge(width=0.75)) +scale_fill_grey() +theme(axis.text=element_text(size=9), axis.title=element_text(size=12),legend.title=element_text(size=10), 
 legend.text=element_text(size=9))+theme(legend.title = element_blank())+labs(x = "")+labs(y = "% mutation")+stat_compare_means(aes(group = classType), label = "p.format",label.y=85,label.x.npc="middle",size = 2)
   
   
 #p2<-recordPlot()
   dev.off()
   
   
  
   pcombined <- plot_grid( as.grob(as.ggplot(p1)), as.grob(as.ggplot(p2)), labels="AUTO",label_size = 10)
   save_plot(paste("AllPatients/TiTV/TiVSTvBothgrey","-",Sys.Date(),".pdf",sep=""),pcombined, ncol = 2,nrow=1,base_aspect_ratio=1.1)
   
   
   
 #  ggarrange(p1,p2 , 
  #           labels = c("A", "B"),
   #          ncol = 2, nrow = 1)
   
  # save_plot("AllPatients/TiTV/TiVSTvBoth.pdf", plot2by2,base_aspect_ratio=1.5)
             #  ncol = 2, # we're saving a grid plot of 2 columns
             
  # dev.off()
   
#-------Lollipop plots------------------------------------------------------
lollipopPlot(maf = laml, gene = 'DNMT3A', AACol = 'Amino_acids', showMutationRate = TRUE)
lollipopPlot(maf = laml2, gene = 'DNMT3A', AACol = 'amino_acid_change', showMutationRate = TRUE)

laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools') 
#------------Mutual Exclusion and co-occuring Plots---------------------------------------
pdf( "AllPatients/MEAndCooccurrence/ME_DiagosisGeneList.pdf", onefile=TRUE)
somaticInteractions(maf = lamlDiagnosisPC,genes=recurrent_genes, pvalue = c(0.05, 0.),returnAll = TRUE)
dev.off()

pdf( "AllPatients/MEAndCooccurrence/ME_RelapseGeneList.pdf", onefile=TRUE)
somaticInteractions(maf = lamlRelapsePC, genes=recurrent_genes, pvalue = c(0.05, 0.1),returnAll = TRUE)
dev.off()


#------------------------------Plot  Bar plots for comparing Cohorts---------------------------------------------
cohorts=read.xlsx("Gene_mutation_frequency_TCGA_TARGET_Linda_cohort_draft_190624.xlsx",sheetName = "Sheet1",header=TRUE)
#To reorder the plot based on TCGA frequency
cohorts=cohorts[order(cohorts$TCGA,decreasing = TRUE),]
colnames(cohorts)<-c("Gene.Name","TCGA", "Linda.Adult.Diagnosis","Linda.Adult.Relapse","TARGET","Linda.Pediatric.Diagnosis","Linda.Pediatric.Relapse")
cohortsFigureU=cohorts[,c("Gene.Name","TCGA", "Linda.Adult.Diagnosis","Linda.Adult.Relapse")]
cohortsFigureU=cohortsFigureU[order(cohortsFigureU$TCGA,decreasing = TRUE),]
cohortsFigureL=cohorts[,c("Gene.Name","TARGET","Linda.Pediatric.Diagnosis","Linda.Pediatric.Relapse")]
cohortsFigureL=cohortsFigureU[order(cohortsFigureL$TCGA,decreasing = TRUE),]
cohortsFigureU=melt(cohortsFigureU)
cohortsFigureL=melt(cohortsFigureL)
pdf( "AllPatients/FrequencyCohorts/Comparisons.pdf", onefile=TRUE)

View(cohortsFigureU)

#Old
cohorts1<-ggplot(cohortsFigureU, aes(Gene.Name,value)) +
  geom_bar(aes(fill = variable), stat = "identity",
           position = position_dodge(0.8), width = 0.7)+scale_fill_grey()+theme(legend.position = "bottom",axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(1,-1,1,0), "mm")) +scale_y_reverse() + coord_flip()
#New trial
cohorts1<-ggplot(cohortsFigureU, aes(Gene.Name,value)) +
  geom_bar(aes(fill = variable),position = position_dodge(1.0), width = 0.7, stat = "identity")+scale_fill_grey()+theme_bw()+theme(legend.position = "left",axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(3,0.5,1,0), "mm"))+labs(fill='Adults Cohort')  +scale_y_reverse() + coord_flip()



#New trial but with modifications
#--------------------------------
#plot.background = element_rect(),    # Background of the entire plot

#panel.background = element_rect(),   # Background of plotting area
#panel.border = element_rect(), 
#panel.grid.major.x = element_line()
#remove the y line using axis.line.y=element_b
#reverse the labels to match with the graph guides(fill = guide_legend(reverse = TRUE)) after coordi flip
#axis .text.x to increase the font of the xaxis after coordin flip
#guides(fill = guide_legend(reverse = TRUE)) to reverse the legend order

#cohortsFigureU$variable <- factor(cohortsFigureU$variable, levels = c("TCGA","Linda.Adult.Diagnosis","Linda.Adult.Relapse"))


#cohortsFigureU$variable <- factor(cohortsFigureU$variable, levels = rev(levels(cohortsFigureU$variable)))

#To reorder the plot based on TCGA frequency

cohortsFigureU$Gene.Name <- factor(cohortsFigureU$Gene.Name, levels = rev(unique(cohortsFigureU$Gene.Name)))
#------------------------------------------------
cohorts1<-ggplot(cohortsFigureU, aes(Gene.Name,value) ) +
  geom_bar(aes(fill = variable),position = position_dodge(0.7), width = 0.7, stat = "identity")+scale_fill_grey()+theme(plot.background = element_rect(),panel.border = element_rect(),panel.background = element_rect(),panel.grid.major.x = element_line(size=0.001,colour='lightgrey'),
 legend.position = "left",legend.title = element_text(size=8.0),legend.text=element_text(size=8.0),axis.line.y = element_blank(),axis.text.x =element_text(size=10.0),axis.text.y = element_blank(),axis.title.x = element_blank(), axis.title.y = element_blank(),  axis.ticks.x = element_blank() ,axis.ticks.y = element_blank(), plot.margin = unit(c(3,0.5,1,0), "mm"))+labs(fill='Adults Cohort')  +scale_y_reverse() + coord_flip()+ guides(fill = guide_legend(reverse = TRUE))

cohorts1<-reposition_legend(cohorts1, 'bottom left')



#Colorful
cohorts1<-ggplot(data = cohortsFigureU, aes(Gene.Name, value)) +
  geom_bar(aes(fill = variable),stat = "identity") + ggtitle("Number of sales staff") +
  theme(legend.position = "left",axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        plot.margin = unit(c(1,-1,1,0), "mm")) +
  scale_y_reverse() + coord_flip()

#Colorful
cohorts2 <- ggplot(data = cohortsFigureL, aes(Gene.Name, value)) +xlab(NULL)+
  geom_bar(aes(fill = variable),stat = "identity") + ggtitle("Sales (x $1000)") +
  theme(legend.position = "right",axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = unit(c(1,0,1,-1), "mm")) +
   coord_flip()

#Old
cohorts2<-ggplot(cohortsFigureL, aes(Gene.Name,value)) +
  geom_bar(aes(fill = variable), stat = "identity",
           position = position_dodge(0.8), width = 0.7)+scale_fill_grey()+
  theme(legend.position = "bottom",axis.title.x = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(1,0,1,-1), "mm")) +
  coord_flip()



#New trial
cohorts2<-ggplot(cohortsFigureL, aes(Gene.Name,value)) +
  geom_bar(aes(fill = variable), stat = "identity",position = position_dodge(1.0), width = 0.7)+scale_fill_grey()+theme_bw()+
  theme(legend.position = "right",axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(3,0,1,-1), "mm"))+labs(fill='Pediatric Cohort') +
  coord_flip()
#above,empty space,below,


#New trial but with modifications
#----------------------------------
#plot.background = element_rect(),    # Background of the entire plot

#panel.background = element_rect(),   # Background of plotting area
#panel.border = element_rect(), 

cohortsFigureL$Gene.Name <- factor(cohortsFigureL$Gene.Name, levels = rev(unique(cohortsFigureU$Gene.Name)))

cohorts2<-ggplot(cohortsFigureL, aes(Gene.Name,value)) +
  geom_bar(aes(fill = variable), stat = "identity",position = position_dodge(0.7), width = 0.7)+scale_fill_grey()+
  theme(plot.background = element_rect(),panel.border = element_rect(),panel.background = element_rect(),panel.grid.major.x = element_line(size=0.001,colour='lightgrey'),
legend.position = "right",legend.title = element_text(size=8.0),legend.text=element_text(size=8.0),axis.text.x =element_text(size=8.0),axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = unit(c(3,0,1,-1), "mm"))+labs(fill='Pediatric Cohort') +
  coord_flip()+guides(fill = guide_legend(reverse = TRUE))
cohorts2<-reposition_legend(cohorts2, 'bottom right')

cohorts$Gene.Name <- factor(cohorts$Gene.Name, levels = rev(unique(cohortsFigureU$Gene.Name)))

g.mid<-ggplot(cohorts,aes(x=1,y=Gene.Name))+geom_text(aes(label=Gene.Name))+
  geom_segment(aes(x=0.94,xend=0.96,yend=Gene.Name))+
  geom_segment(aes(x=1.04,xend=1.06,yend=Gene.Name))+
 # geom_segment(aes(x=0.94,xend=0.96,yend=Gene.Name))+
  ggtitle("")+
  ylab(NULL)+
  scale_x_continuous(expand=c(0,0),limits=c(0.94,1.065))+
  theme(axis.title=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_blank(),
       axis.ticks.y=element_blank(),
        panel.background=element_blank(),
        axis.text.x=element_text(color=NA),
        axis.ticks.x=element_line(color=NA),
        plot.margin =  unit(c(1,-1,2.5,-1), "mm"))


gg1 <- ggplot_gtable(ggplot_build(as.ggplot(cohorts1)))
gg2 <- ggplot_gtable(ggplot_build(as.ggplot(cohorts2)))
gg.mid <- ggplot_gtable(ggplot_build(g.mid))

cohortsFigure<-grid.arrange(gg1,gg.mid,gg2,ncol=3,widths=c(2/9,1/9,2/9))
save_plot(paste("AllPatients/FrequencyCohorts/Comparisons","-",Sys.Date(),".pdf",sep=""),cohortsFigure, ncol = 2,nrow=2,base_aspect_ratio=1)


dev.off()

cohorts3<-multiplot(cohorts1, cohorts2, cols = 2)


ggbarplot(cohortsFigureU, x = "Gene.Name", y = "value",
          fill = "variable",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "desc",          # Sort the value in descending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "MPG z-score",
          legend.title = "group",
          rotate = TRUE,
          ggtheme = theme_minimal()
)


 