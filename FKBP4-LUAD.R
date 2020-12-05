# Run in R

# Load packages
library("survival")
library("survminer")
library("dplyr")

# Load data
rna.seq <- as.matrix(read.table("RNA_Seq_HiSeqV2", sep="\t",header=TRUE,row.names = 1))
survival <- read.table("Survival_LUAD_survival.txt", sep="\t",header=TRUE, row.names = 1 )
clin.data <- read.table("Clinical_LUAD_clinicalMatrix", sep= "\t", header=T, row.names=1)
prot.exp <- as.matrix(read.table("Protein_RPPA_RBN", sep="\t", header=T, row.names=1))
dna.meth <- read.table("HumanMethylation450", sep="\t", header=T, row.names=1)
CNV <- as.matrix(read.table("Gistic2_CopyNumber_Gistic2_all_data_by_genes",sep="\t",head=T,row.names=1))

# Loading DNA Methylation annotation
annot.dna.meth <- readRDS("annot450k.rds") # annotation was taken from Imperial College BRC Server

# Fixing PX ID format
rownames(clin.data)<-gsub(rownames(clin.data), pattern="-", replace=".")
rownames(survival)<-gsub(rownames(survival), pattern="-", replace=".")

# RPPA and survival
OS.Time.prot <- survival[colnames(prot.exp),"OS.time"]
OS.Event.prot <- as.numeric(survival[colnames(prot.exp),"OS"])
OS.prot <- Surv(OS.Time.prot,OS.Event.prot)

# RRPA cox regression
Results.OS_prot.exp<- array(NA, c(nrow(prot.exp),4))
colnames(Results.OS_prot.exp)<-c("HR","LCI","UCI","PVAL")
rownames(Results.OS_prot.exp)<-rownames(prot.exp)
Results.OS_prot.exp <- as.data.frame(Results.OS_prot.exp)

for(i in 1:nrow(prot.exp)){
  coxphmodel2 <- coxph(OS.prot~ as.numeric(prot.exp[i,]))
  Results.OS_prot.exp$HR[i] <- summary(coxphmodel2)$coef[1,2]
  Results.OS_prot.exp$LCI[i] <- summary(coxphmodel2)$conf.int[1,3]
  Results.OS_prot.exp$UCI[i] <- summary(coxphmodel2)$conf.int[1,4]
  Results.OS_prot.exp$PVAL[i] <- summary(coxphmodel2)$coef[1,5]
}

# Adjusting for multiple testing using FDR method
Results.OS_prot.exp$FDR <- p.adjust(Results.OS_prot.exp$PVAL,method="fdr") 
Results.OS_prot.exp<-Results.OS_prot.exp[order(Results.OS_prot.exp$FDR, decreasing=F),]

# Check RPPA data with FDR < 0.05
Results.OS_prot.exp   # no protein were statistically significant after adjusting the p-values


###


# RNAseq and survival
OS.Time.rna <- survival[colnames(rna.seq),"OS.time"]
OS.Event.rna <- as.numeric(survival[colnames(rna.seq),"OS"])
OS.rna <- Surv(OS.Time.rna,OS.Event.rna)

# RNAseq cox regression
Results.OS_rna.seq<- array(NA, c(nrow(rna.seq),4))
colnames(Results.OS_rna.seq)<-c("HR","LCI","UCI","PVAL")
rownames(Results.OS_rna.seq)<-rownames(rna.seq)
Results.OS_rna.seq <- as.data.frame(Results.OS_rna.seq)

for(i in 1:nrow(rna.seq)){
  coxphmodel <- coxph(OS.rna~ as.numeric(rna.seq[i,]))
  Results.OS_rna.seq$HR[i] <- summary(coxphmodel)$coef[1,2]
  Results.OS_rna.seq$LCI[i] <- summary(coxphmodel)$conf.int[1,3]
  Results.OS_rna.seq$UCI[i] <- summary(coxphmodel)$conf.int[1,4]
  Results.OS_rna.seq$PVAL[i] <- summary(coxphmodel)$coef[1,5]
}

# Adjusting for multiple testing using FDR method
Results.OS_rna.seq$FDR <- p.adjust(Results.OS_rna.seq$PVAL,method="fdr")
Results.OS_rna.seq<-Results.OS_rna.seq[order(Results.OS_rna.seq$FDR, decreasing=F),]

# Check gene data with FDR < 0.05
Results.OS_rna.seq # statiscally significant genes were detected - proceed to univariate analysis

# Percentage of genes that have FDR < 0.05
AA <- length(which(Results.OS_rna.seq[,5] < 0.05))
Aa <- length(Results.OS_rna.seq[,5])
((AA/Aa)*100) # 6.02%

# RNASeq multivariate cox regression
clin.data <- clin.data[colnames(rna.seq),]

# confounding factors from clinical data
# Gender
gender <- rep(NA, nrow(clin.data))
gender[clin.data$gender=="FEMALE"] <- 0
gender[clin.data$gender=="MALE"] <- 1

# age
age<-as.numeric(clin.data$age_at_initial_pathologic_diagnosis)

# histological type
histology <- as.factor(clin.data$histological_type)

# anatomic neoplasm subdivisioon
neoplasm <- as.factor(clin.data$anatomic_neoplasm_subdivision)

# location
location <- as.factor(clin.data$location_in_lung_parenchyma)

# Pathologic stage
stage.III <- grep("III", clin.data$pathologic_stage)
stage.IV <- grep("IV", clin.data$pathologic_stage)
stage.high <- rep(0, nrow((clin.data)))
stage.high [c(stage.III,stage.IV)] <- 1

# Smoking
# 0 = non-smokers, 1 = active smokers and Current reformed smoker for <, > or = 15 years
smoke.pack <-clin.data$number_pack_years_smoked
smoke.pack[which(clin.data$tobacco_smoking_history==1)]<-0
smoke <-rep(NA,nrow(clin.data))
smoke[which(clin.data$tobacco_smoking_history==1)]<-0
smoke[which(clin.data$tobacco_smoking_history>1)]<-1

# new tumour event prior initial treatment (NTE)
NTE <- rep(NA, nrow(clin.data))
NTE[clin.data$new_tumor_event_after_initial_treatment=="YES"] <- 1
NTE[clin.data$new_tumor_event_after_initial_treatment=="NO"] <- 0

# PX received target molecular therapy
targeted.therapy <-rep(NA,nrow(clin.data))
targeted.therapy[which(clin.data$targeted_molecular_therapy=="NO")]<-0
targeted.therapy[which(clin.data$targeted_molecular_therapy=="YES")]<-1

# Summary of confounding factors
summary(coxph(OS.rna ~ gender))$coef                # no statistical difference
summary(coxph(OS.rna ~ age))$coef                   # no statistical difference
summary(coxph(OS.rna ~ histology))$coef             # no statistical difference
summary(coxph(OS.rna ~ neoplasm))$coef              # no statistical difference
summary(coxph(OS.rna ~ stage.high))$coef            # SIGNIFICANT (p < 3.378704e-11)
summary(coxph(OS.rna ~ smoke))$coef                 # no statistical difference
summary(coxph(OS.rna ~ smoke.pack))$coef            # no statistical difference
summary(coxph(OS.rna ~ NTE))$coef                   # SIGNIFICANT (p < 2.614033e-10)
summary(coxph(OS.rna ~ targeted.therapy))$coef      # no statistical difference
summary(coxph(OS.rna ~ location))$coef              # no statistical difference

# NTE cox regression
coxph(OS.rna ~ NTE)
Results.OS_NTE<-array(NA, c(nrow(rna.seq),4))
colnames(Results.OS_NTE)<-c("HR","LCI","UCI","PVAL")
rownames(Results.OS_NTE)<-rownames(rna.seq)
Results.OS_NTE<-as.data.frame(Results.OS_NTE)

for(i in 1:nrow(rna.seq)){
  coxphmodel_NTE <- coxph(OS.rna ~ NTE)
  Results.OS_NTE$HR[i] <- summary(coxphmodel_NTE)$coef[1,2]
  Results.OS_NTE$LCI[i] <- summary(coxphmodel_NTE)$conf.int[1,3]
  Results.OS_NTE$UCI[i] <- summary(coxphmodel_NTE)$conf.int[1,4]
  Results.OS_NTE$PVAL[i] <- summary(coxphmodel_NTE)$coef[1,5]
}
Results.OS_NTE$FDR <- p.adjust(Results.OS_NTE$PVAL, method =  "fdr")

summary(coxphmodel_NTE)$coef[1,2]       # HR = 2.666519
summary(coxphmodel_NTE)$conf.int[1,3]   # LCI = 1.967217
summary(coxphmodel_NTE)$conf.int[1,4]   # HCI= 3.614407
summary(coxphmodel_NTE)$coef[1,5]       # p = 2.614033e-10; FDR = 2.614033e-10


# Stage cox regression
coxph(OS.rna ~ stage.high)
Results.OS_stage<-array(NA, c(nrow(rna.seq),4))
colnames(Results.OS_stage)<-c("HR","LCI","UCI","PVAL")
rownames(Results.OS_stage)<-rownames(rna.seq)
Results.OS_stage<-as.data.frame(Results.OS_stage)

for(i in 1:nrow(rna.seq)){
  coxphmodel_stage <- coxph(OS.rna ~ stage.high)
  Results.OS_stage$HR[i] <- summary(coxphmodel_stage)$coef[1,2]
  Results.OS_stage$LCI[i] <- summary(coxphmodel_stage)$conf.int[1,3]
  Results.OS_stage$UCI[i] <- summary(coxphmodel_stage)$conf.int[1,4]
  Results.OS_stage$PVAL[i] <- summary(coxphmodel_stage)$coef[1,5]
}
Results.OS_stage$FDR <- p.adjust(Results.OS_stage$PVAL, method =  "fdr")

summary(coxphmodel_stage)$coef[1,2]       # HR = 2.643894
summary(coxphmodel_stage)$conf.int[1,3]   # LCI = 1.983362
summary(coxphmodel_stage)$conf.int[1,4]   # HCI =3.524406
summary(coxphmodel_stage)$coef[1,5]       # p = 3.378704e-11; FDR = 3.378704e-11


# RNASeq multivariate cox regression
# Only accounted for stage.high and new.tumour.pre
# Did not adjust for other factors
Results.OS_rna.seq_fctrs<-array(NA, c(nrow(rna.seq),4))
colnames(Results.OS_rna.seq_fctrs)<-c("HR","LCI","UCI","PVAL")
rownames(Results.OS_rna.seq_fctrs)<-rownames(rna.seq)
Results.OS_rna.seq_fctrs<-as.data.frame(Results.OS_rna.seq_fctrs)

for(i in 1:nrow(rna.seq)){
  coxphmodel3 <- coxph(OS.rna ~ as.numeric(rna.seq[i,]+stage.high+NTE))
  Results.OS_rna.seq_fctrs$HR[i] <- summary(coxphmodel3)$coef[1,2]
  Results.OS_rna.seq_fctrs$LCI[i] <- summary(coxphmodel3)$conf.int[1,3]
  Results.OS_rna.seq_fctrs$UCI[i] <- summary(coxphmodel3)$conf.int[1,4]
  Results.OS_rna.seq_fctrs$PVAL[i] <- summary(coxphmodel3)$coef[1,5]
}

Results.OS_rna.seq_fctrs <- as.data.frame(Results.OS_rna.seq_fctrs) # Multivariate analysis

# Adjust for multiple testing using FDR
Results.OS_rna.seq_fctrs$FDR <- p.adjust(Results.OS_rna.seq_fctrs$PVAL, method =  "fdr")
Results.OS_rna.seq_fctrs<-Results.OS_rna.seq_fctrs[order(Results.OS_rna.seq_fctrs$FDR, decreasing=F),]

# Identify potential candidates
Results.OS_rna.seq_fctrs

# Identified FKBP4 to be a good novel candidate based on adj.p.vals and Pubmed Search
# Stratifying PX based on the expression of FKBP4 (median)
summary(rna.seq["FKBP4",])
Results.OS_rna.seq["FKBP4",]        # Univariate cox regression
Results.OS_rna.seq_fctrs["FKBP4",]  # Multivariate cox regression
FKBP4.high <- as.numeric(rna.seq["FKBP4",]>median(rna.seq["FKBP4",]))

# KM plot
png("KM_plots.png", width=9,height=9,units='in',res=300)
grid <- matrix(c(1,1,2,3), nrow = 2, ncol = 2, byrow = T)
layout(grid)
plot(survfit(OS.rna ~ FKBP4.high), col=c("#E7B800", "#2E9FDF"),lwd=2,mark.time=TRUE, xlab="OS Time (days)", ylab="Survival Probability")
  legend("topright",legend=c("High FKBP4 expression","Low FKBP4 expression"),col=c("#2E9FDF","#E7B800"),lwd=2)
  text(4000,0.6,"HR=1.63 (95% CI, 1.43-1.85)")
  text(5100,0.6, "FDR")
  text(5500,0.6, "< 0.001")
plot(survfit(OS.rna ~ NTE), col=c("#E7B800", "#2E9FDF"),lwd=2,mark.time=TRUE, xlab="OS Time (days)", ylab="Survival Probability")
  legend("topright",legend=c("NTE+", "NTE-"),col=c("#2E9FDF","#E7B800"),lwd=2)
  text(4800,0.7,"HR=2.67 (95% CI, 1.97-3.61)")
  text(4300,0.65, "FDR")
  text(5200,0.65, "< 0.001")
plot(survfit(OS.rna ~ stage.high), col=c("#E7B800", "#2E9FDF"),lwd=2,mark.time=TRUE, xlab="OS Time (days)", ylab="Survival Probability")
  legend("topright",legend=c("High Stage", "Low Stage"),col=c("#2E9FDF","#E7B800"),lwd=2)
  text(4800,0.7,"HR=2.64 (95% CI, 1.98-3.52)")
  text(4300,0.65, "FDR")
  text(5200,0.65, "< 0.001")
dev.off()

# FKBP4 KM-plot: dependencies" survminer" and "dpylr"
png("KMplot_RiskTable.png")
ggsurvplot(survfit(Surv(OS.Time.rna, OS.Event.rna) ~ FKBP4.high, data = survival), conf.int = F, risk.table = "absolute", risk.table.y.text.col = TRUE, palette = c("#E7B800", "#2E9FDF"),
           xlab = "OS Time (days)", legend = "bottom", xlim = c(0,7500), break.time.by = 1000, legend.labs = c("High FKBP4", "Low FKBP4"))
dev.off()


##


# Higher FKBP4 expression results to poorer prognosis
# DNA Methylation of FKBP4
rna.seq.meth<-rna.seq[,which(is.element(colnames(rna.seq),colnames(dna.meth)))]
dna.meth2<- dna.meth[,which(is.element(colnames(dna.meth),colnames(rna.seq.meth)))]

# Align PX IDs
rna.seq.meth <- as.matrix(rna.seq.meth[,order(colnames(rna.seq.meth))])
dna.meth2 <- as.matrix(dna.meth2[,order(colnames(dna.meth2))])

# Check for methylated CpG site on FKBP4 gene
FKBP4.meth <- rownames(annot.dna.meth[which(annot.dna.meth$UCSC_RefGene_Name=="FKBP4"),])
FKBP4.meth

# Check annotation
annot.dna.meth[FKBP4.meth,]

# Filtering methylated CpG islands in FKBP4 gene
meth.data.FKBP4 <- dna.meth2[FKBP4.meth,]

# Exclusion CpG islands that were missing in 50% of the PX
NA.Count_FKBP4.meth<-apply(meth.data.FKBP4,1,function(x) sum(as.numeric(is.na(x))))
Exclude<-as.numeric(NA.Count_FKBP4.meth>0.5*ncol(meth.data.FKBP4))
meth.data.FKBP4<-meth.data.FKBP4[which(Exclude==0),]

# Correlation test between DNA methylation (beta-values) and RNAseq (log2(x+1), RBN): Spearman correlation
Result.FKBP4.meth <- array(NA,c(nrow(meth.data.FKBP4),4))
rownames(Result.FKBP4.meth)<-rownames(meth.data.FKBP4)
colnames(Result.FKBP4.meth)<-c("Cor.FKBP4","Cor.test.FKBP4","Mean.high.FKBP4","Mean.low.FKBP4")
FKBP4.high.meth <- as.numeric(rna.seq.meth["FKBP4",]>median(rna.seq.meth["FKBP4",]))

for (i in 1:nrow(meth.data.FKBP4)){
  Result.FKBP4.meth [i,1]<-cor.test(as.numeric(rna.seq.meth["FKBP4",]),as.numeric(meth.data.FKBP4[i,]),method="spearman", use="c")$est
  Result.FKBP4.meth [i,2]<-cor.test(as.numeric(rna.seq.meth["FKBP4",]),as.numeric(meth.data.FKBP4[i,]),method="spearman", use="c")$p.value
}
Result.FKBP4.meth[,3]<-apply(meth.data.FKBP4[,which(FKBP4.high.meth==1)],1,mean,na.rm=T)
Result.FKBP4.meth[,4]<-apply(meth.data.FKBP4[,which(FKBP4.high.meth==0)],1,mean,na.rm=T)
Result.FKBP4.meth

# Beta-value - more intuitive in biological interpretion
# M-value - more statistically valid for the differential analysis of methylation levels.

# Produce an object for DNA Methylation M-values
# Convert beta-values to m-values. Formula: mvalues = log2(x/(1-x))
meth.data.FKBP4.Mvals <- apply(meth.data.FKBP4, MARGIN = 2, function(x) log2(x/(1-x)))

# Correlation test between DNA methylation (M-values) and RNAseq (log2(x+1), RBN): Spearman correlation
Result.FKBP4.meth.MVals <- array(NA,c(nrow(meth.data.FKBP4.Mvals),4))
rownames(Result.FKBP4.meth.MVals)<-rownames(meth.data.FKBP4.Mvals)
colnames(Result.FKBP4.meth.MVals)<-c("Cor.FKBP4","Cor.test.FKBP4","Mean.high.FKBP4","Mean.low.FKBP4")
FKBP4.high.meth.Mvals <- as.numeric(rna.seq.meth["FKBP4",]>median(rna.seq.meth["FKBP4",]))

for (i in 1:nrow(meth.data.FKBP4.Mvals)){
  Result.FKBP4.meth.MVals[i,1]<-cor.test(as.numeric(rna.seq.meth["FKBP4",]),as.numeric(meth.data.FKBP4.Mvals[i,]),method="spearman", use="c")$est
  Result.FKBP4.meth.MVals[i,2]<-cor.test(as.numeric(rna.seq.meth["FKBP4",]),as.numeric(meth.data.FKBP4.Mvals[i,]),method="spearman", use="c")$p.value
}
Result.FKBP4.meth.MVals[,3]<-apply(meth.data.FKBP4.Mvals[,which(FKBP4.high.meth.Mvals==1)],1,mean,na.rm=T)
Result.FKBP4.meth.MVals[,4]<-apply(meth.data.FKBP4.Mvals[,which(FKBP4.high.meth.Mvals==0)],1,mean,na.rm=T)


# Check for CpG island
Result.FKBP4.meth.MVals
Result.FKBP4.meth

# Identified "cg04611395" to be correlated with largest change in Beta-value abd M-values between groups of high and low FKBP4 expression 
#               Cor.FKBP4 Cor.test.FKBP4 Mean.high.FKBP4 Mean.low.FKBP4
# cg04611395  0.104384099    0.022605706      0.41217710     0.37383766
# Difference in Beta-value  =  0.038
# Difference in M-value     =  0.22

# Most significant CpG Island =cg04611395 
# Plots RRA Expression and DNA Methylation for FKBP4 gene
png("GeneExpVsMeth.png",width=6,height=6,units='in',res=300)
plot(as.numeric(meth.data.FKBP4["cg04611395",]),as.numeric(rna.seq.meth["FKBP4",]), xlab="cg04611395 (M-values)", ylab="FKBP4 expression")
abline(lm(rna.seq.meth["FKBP4",]~dna.meth2["cg04611395",]), col="red")
text(0.7, 8.5, "Cor = 0.104 (p < 0.05)")
dev.off()


##

# CNV and FKBP4
# Visually check the raw data
rna.seq.CNV<-rna.seq[,which(is.element(colnames(rna.seq),colnames(CNV)))]
CNV<- CNV[,which(is.element(colnames(CNV),colnames(rna.seq.CNV)))]

FKBP4.CNV <- CNV["FKBP4",]
summary(FKBP4.CNV)

# CNV plot
png("hist_CNV.png",width=6,height=6,units='in',res=300)
hist(FKBP4.CNV,xlab = "GISTIC2 gene-level estimate score", ylab = "Frequency",ylim = c(0,100), breaks=50 , col="blue", border=F,main="")
dev.off()



