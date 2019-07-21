mutationsAll=fread("AML_WGS_Strelka_Filtered.csv",header=TRUE)
Unpaired=fread("MEAndCooccurence.csv",header=TRUE)
Unpaired_LatestD=Unpaired$`Unpaired: All D (WGS and WES)`
Unpaired_LatestR=Unpaired$`Unpaired: All latest  R (WGS and WES)`
recurrent_genes=c("ARID1A","ASXL1","ASXL2","BCOR","BCORL1","CEBPA","CSF3R","DNMT3A","EZH2","FAM5C","FLT3","GATA1","GATA2","H3F3A","HNRNPK","IDH1","IDH2","IKZF1","KIT","KMT2A","KRAS","MGA","MYC","NF1","NPM1","NRAS","NRNX3","PHF6","PTPN11","RAD21","RUNX1","SETD2","SMC1A","SMC3","SRSF1","SRSF2","STAG2","TET2","TP53","U2AF1","UBTF","WT1")
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



#New trial but with latest modifications
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
#cohorts1<-ggplot(data = cohortsFigureU, aes(Gene.Name, value)) +
 # geom_bar(aes(fill = variable),stat = "identity") + ggtitle("Number of sales staff") +
#  theme(legend.position = "left",axis.title.x = element_blank(), 
 #       axis.title.y = element_blank(), 
#        axis.text.y = element_blank(), 
#        axis.ticks.y = element_blank(), 
#        plot.margin = unit(c(1,-1,1,0), "mm")) +
#  scale_y_reverse() + coord_flip()

#Colorful
#cohorts2 <- ggplot(data = cohortsFigureL, aes(Gene.Name, value)) +xlab(NULL)+
#  geom_bar(aes(fill = variable),stat = "identity") + ggtitle("Sales (x $1000)") +
#  theme(legend.position = "right",axis.title.x = element_blank(), axis.title.y = element_blank(), 
#        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
#      plot.margin = unit(c(1,0,1,-1), "mm")) +
 # coord_flip()

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


#New trial but with latest modifications
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

#After migration

