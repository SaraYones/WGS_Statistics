#mutationsAll=fread("AML_WGS_Strelka_Filtered.csv",header=TRUE)
mutationsMECO=fread("Genes_IDs_ME_and_Co_occurrence_190717.csv",header=TRUE)
Unpaired=fread("MEAndCooccurence.csv",header=TRUE)

PediatricMECO=read.xlsx("MEAndCooccurence.xlsx",sheetName = "Pediatric",header=TRUE)
AdultsMECO=read.xlsx("MEAndCooccurence.xlsx",sheetName = "Adults",header=TRUE)

#Old (Unpaired)
Unpaired_LatestD=Unpaired$`Unpaired: All D (WGS and WES)`
Unpaired_LatestR=Unpaired$`Unpaired: All latest  R (WGS and WES)`
recurrent_genes=c("ARID1A","ASXL1","ASXL2","BCOR","BCORL1","CEBPA","CSF3R","DNMT3A","EZH2","FAM5C","FLT3","GATA1","GATA2","H3F3A","HNRNPK","IDH1","IDH2","IKZF1","KIT","KMT2A","KRAS","MGA","MYC","NF1","NPM1","NRAS","NRNX3","PHF6","PTPN11","RAD21","RUNX1","SETD2","SMC1A","SMC3","SRSF1","SRSF2","STAG2","TET2","TP53","U2AF1","UBTF","WT1")
mutationsAll=as.data.frame(mutationsAll)
#---------------------------------------------------------
meta_data=fread("AML_Filtered_Metadata.csv")
View(meta_data[order(meta_data$"Age at AML onset (Years)"),])
#Sort by Age and replace all occurences of _ with - in order to appear on titv

meta_data=meta_data[order(meta_data$"Age at AML onset (Years)"),]
meta_data$`Sample ID`=str_replace(meta_data$`Sample ID`, "_", "-")
meta_data=as.data.frame(meta_data)


#This is all old script that was used before decided on the matched pairs
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
#-------------------------------------------------------------------------------
mutationDiagnosisPediatric=mutationsMECO[mutationsMECO$Tumor_Sample_Barcode %in% PediatricMECO$PediatricDiagnosis,]
mutationRelapsePediatric=mutationsMECO[mutationsMECO$Tumor_Sample_Barcode %in% PediatricMECO$PediatricRelapse,]
mutationDiagnosisAdult=mutationsMECO[mutationsMECO$Tumor_Sample_Barcode %in% AdultsMECO$AdultDiagnosis,]
mutationRelapseAdult=mutationsMECO[mutationsMECO$Tumor_Sample_Barcode %in% AdultsMECO$AdultRelapse,]

mutationDiagnosisAll=mutationsMECO[mutationsMECO$Tumor_Sample_Barcode %in% append(levels(PediatricMECO$PediatricDiagnosis),levels(AdultsMECO$AdultDiagnosis)),]
mutationRelapseAll=mutationsMECO[mutationsMECO$Tumor_Sample_Barcode %in% append(levels(PediatricMECO$PediatricRelapse),levels(AdultsMECO$AdultRelapse)),]

lamlDiagnosisPediatric=read.maf(maf = mutationDiagnosisPediatric,useAll = TRUE)
lamlRelapsePediatric=read.maf(maf = mutationRelapsePediatric,useAll = TRUE)
lamlDiagnosisAdult=read.maf(maf = mutationDiagnosisAdult,useAll = TRUE)
lamlRelapseAdult=read.maf(maf = mutationRelapseAdult,useAll = TRUE)
lamlDiagnosisAll=read.maf(maf = mutationDiagnosisAll,useAll = TRUE)
lamlRelapseAll=read.maf(maf = mutationRelapseAll,useAll = TRUE)

recurrent_genes=unique(mutationsMECO$Hugo_Symbol)


pdf( paste("AllPatients/MEAndCooccurrence/ME_DiagosisPediatric","-",Sys.Date(),".pdf",sep=""), onefile=TRUE)
somaticInteractions(maf = lamlDiagnosisPediatric,genes=recurrent_genes, pvalue = c(0.05, 0.1),returnAll = TRUE)
dev.off()


pdf( paste("AllPatients/MEAndCooccurrence/ME_RelapsePediatric","-",Sys.Date(),".pdf",sep=""), onefile=TRUE)

somaticInteractions(maf = lamlRelapsePediatric, genes=recurrent_genes, pvalue = c(0.05, 0.1),returnAll = TRUE)
dev.off()



pdf( paste("AllPatients/MEAndCooccurrence/ME_DiagosisAdult","-",Sys.Date(),".pdf",sep=""), onefile=TRUE)
somaticInteractions(maf = lamlDiagnosisAdult,genes=recurrent_genes, pvalue = c(0.05, 0.1),returnAll = TRUE)
dev.off()


pdf( paste("AllPatients/MEAndCooccurrence/ME_RelapseAdult","-",Sys.Date(),".pdf",sep=""), onefile=TRUE)

somaticInteractions(maf = lamlRelapseAdult, genes=recurrent_genes, pvalue = c(0.05, 0.1),returnAll = TRUE)
dev.off()


pdf( paste("AllPatients/MEAndCooccurrence/ME_DiagosisAll","-",Sys.Date(),".pdf",sep=""), onefile=TRUE)
somaticInteractions(maf = lamlDiagnosisAll,genes=recurrent_genes, pvalue = c(0.05, 0.1),returnAll = TRUE)
dev.off()


pdf( paste("AllPatients/MEAndCooccurrence/ME_RelapseAll","-",Sys.Date(),".pdf",sep=""), onefile=TRUE)

somaticInteractions(maf = lamlRelapseAll, genes=recurrent_genes, pvalue = c(0.05, 0.1),returnAll = TRUE)
dev.off()

