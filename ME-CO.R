mutationsAll=fread("AML_WGS_Strelka_Filtered.csv",header=TRUE)
Unpaired=fread("MEAndCooccurence.csv",header=TRUE)

PediatricMECO=read.xlsx("MEAndCooccurence.xlsx",sheetName = "Pediatric",header=TRUE)
AdultsMECO=read.xlsx("MEAndCooccurence.xlsx",sheetName = "Adults",header=TRUE)
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


pdf( "AllPatients/MEAndCooccurrence/ME_DiagosisGeneList.pdf", onefile=TRUE)
somaticInteractions(maf = lamlDiagnosisPC,genes=recurrent_genes, pvalue = c(0.05, 0.),returnAll = TRUE)
dev.off()

pdf( "AllPatients/MEAndCooccurrence/ME_RelapseGeneList.pdf", onefile=TRUE)
somaticInteractions(maf = lamlRelapsePC, genes=recurrent_genes, pvalue = c(0.05, 0.1),returnAll = TRUE)
dev.off()
