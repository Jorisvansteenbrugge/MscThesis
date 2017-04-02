set.seed(1396558)
library("plyr")
library("reshape")
library("igraph")
library('knitr')
library("hexbin")
library("RColorBrewer")
library("ctr")

RNAseq_Annotated_Matrix <- read.table("/home/dev/tools/ctr/sample_data/EBPR_RNAseq_Annotation_Matrix",
                      sep = "\t", quote = "", row.names = NULL, stringsAsFactors = F,
                      fill=T)

KO_pathways_table <- read.table("/home/dev/tools/ctr/sample_data/ko00001.keg_tab_delimited",
                   sep="\t", quote = "", row.names = NULL, stringsAsFactors = F)

#Reformatting
Clade_IIA_bins<-c(39,61,71,92,96,99)
RNAseq_Annotated_Matrix$Bin<-gsub(".*\\.(.*)\\..*", "\\1", RNAseq_Annotated_Matrix[,1])
RNAseq_Annotated_Matrix_full<-RNAseq_Annotated_Matrix
RNAseq_Annotated_Matrix$Bin[which(RNAseq_Annotated_Matrix$Bin %in% Clade_IIA_bins)]<-39
high_quality_bins<-c(8,28,25,7,46,39,22,38,54,53,48,45,31,42,16,33,26,40,36,21,27,17,19,32,14,11,30,43,35,29,23,58,41,20,15,37,49,50)
RNAseq_Annotated_Matrix<-RNAseq_Annotated_Matrix[which(RNAseq_Annotated_Matrix$Bin %in% high_quality_bins),]

# Give column names
sample_names<-rep(NA,dim(RNAseq_Annotated_Matrix)[2]-3)
for (i in 1:dim(RNAseq_Annotated_Matrix)[2]-3) {
  sample_names[i]<-paste("Sample",i,sep="")
}

colnames(RNAseq_Annotated_Matrix)<-c("Locus_ID",sample_names,"KO","Bin")
high_quality_bins<-names(table(RNAseq_Annotated_Matrix$Bin))
kable(RNAseq_Annotated_Matrix[1:5,], caption = "Raw Data")

All_KOs<-names(which(table(RNAseq_Annotated_Matrix$KO)>=5))[-1]
PHA_module <- c("K01895","K00925","K00625","K00626","K00023", "K03821")
polyP_module <- c("K00937","K02040","K02037","K02038","K02036","K07636","K07657","K03306")
Glycogen_module<- c("K00700","K00688","K00703","K01214","K15778","K00975","K00845")

no_feature<- c(9159700,4459877,9826273,8171512,9542765,10522313)
ambiguous<- c(3940698,2023389,4675033,3308789,6446272,5966543)
not_aligned<- c(0,19317660,0,0,0,0)

RNAseq_Annotated_Matrix<-RNAseq_Normalize(RNAseq_Annotated_Matrix,no_feature,ambiguous,not_aligned)
kable(RNAseq_Annotated_Matrix[1:5,], caption = "Normalized Data")







RNAseq_Normalize(RNAseq_Annotated_Matrix,no_feature,ambiguous,not_aligned)
test <- RNAseq_Normalize(RNAseq_Annotated_Matrix,no_feature,ambiguous,not_aligned)
      

library(rbenchmark)
benchmark(replications=c(1),RNAseq_NormalizEe(RNAseq_Annotated_Matrix,no_feature,ambiguous,not_aligned), 
          RNAseq_Normalize(RNAseq_Annotated_Matrix,no_feature,ambiguous,not_aligned))




#############################################################################################
RNAseq_Annotated_Matrix<-Create_Rank_Columns(RNAseq_Annotated_Matrix)
kable(RNAseq_Annotated_Matrix[1:5,], caption = "Data with Rank Column added")

RNAseq_Annotation_Matrix_no_sd_of_zero<-which_rows_with_no_sd(RNAseq_Annotated_Matrix)

##############
library(rbenchmark)
benchmark(replications=rep(10,1), RNAseq_Normalize(RNAseq_Annotated_Matrix,no_feature,ambiguous,not_aligned))
##############

#Slow
I_KOs_Background <- Individual_KOs_Background(RNAseq_Annotation_Matrix_no_sd_of_zero,100)

#Calculating significant distances (person and euclidean)
t.test_KO_random_pearson<-t.test(I_KOs_Background$KO_pairwise_gene_correlation,I_KOs_Background$random_pairwise_gene_correlation, alternative="greater") # x > y (NULL)
t.test_H_KO_H_random_pearson<-t.test(I_KOs_Background$H_KO_pairwise_gene_correlation,I_KOs_Background$H_random_pairwise_gene_correlation, alternative="greater")
t.test_H_KO_KO_pearson<-t.test(I_KOs_Background$H_KO_pairwise_gene_correlation,I_KOs_Background$KO_pairwise_gene_correlation, alternative="greater")
t.test_KO_random_euclidean<-t.test(I_KOs_Background$KO_pairwise_gene_euclidean,I_KOs_Background$random_pairwise_gene_euclidean, alternative="less") # x > y (NULL)
t.test_H_KO_H_random_euclidean<-t.test(I_KOs_Background$H_KO_pairwise_gene_euclidean,I_KOs_Background$H_random_pairwise_gene_euclidean, alternative="less")
t.test_H_KO_KO_euclidean<-t.test(I_KOs_Background$H_KO_pairwise_gene_euclidean,I_KOs_Background$KO_pairwise_gene_euclidean, alternative="less")

#printing P-values
t.test_KO_random_pearson$p.value
t.test_H_KO_H_random_pearson$p.value
t.test_KO_random_euclidean$p.value
t.test_H_KO_H_random_euclidean$p.value


#Comparing the two KO distributions###################################################################################333
par(mfrow=c(2,2),mar=c(3,3,3,1))
# plot 1
plot(density(I_KOs_Background$random_pairwise_gene_correlation,adjust = 2,na.rm=TRUE),ylim=c(0,1),xlab="",ylab="",main="")
points(density(I_KOs_Background$KO_pairwise_gene_correlation,adjust = 2),typ="l",col="blue")
mtext(paste("p-value = ",signif(t.test_KO_random_pearson$p.value,2)),side=3,col="blue",padj=2,cex=.75)
title(ylab="Density", line=2, cex.lab=1)
title(xlab="PC", line=2, cex.lab=1)

# plot 2
plot(density(I_KOs_Background$H_random_pairwise_gene_correlation,adjust = 2),ylim=c(0,1),xlab="",ylab="",main=" ")
points(density(I_KOs_Background$H_KO_pairwise_gene_correlation,adjust = 2),typ="l",col="red")
mtext(paste("p-value = ",signif(t.test_H_KO_H_random_pearson$p.value,2)),side=3,col="red",padj=2,cex=.75)
title(ylab="Density", line=2, cex.lab=1)
title(xlab="PC", line=2, cex.lab=1)

# plot 3
plot(density(I_KOs_Background$random_pairwise_gene_euclidean,adjust = 2),typ="l" ,ylim=c(0,1.25),xlab="",ylab="",main="")
points(density(I_KOs_Background$KO_pairwise_gene_euclidean,adjust = 2),typ="l",col="blue")
title(ylab="Density", line=2, cex.lab=1)
title(xlab="NRED", line=2, cex.lab=1)
mtext(paste("p-value = ",signif(t.test_KO_random_euclidean$p.value,2)),side=3,col="blue",padj=2,cex=.75)

# plot 4
plot(density(I_KOs_Background$H_random_pairwise_gene_euclidean,adjust = 2),typ="l" ,ylim=c(0,1.25),xlab="",ylab="",main="")
points(density(I_KOs_Background$H_KO_pairwise_gene_euclidean,adjust = 2),typ="l",col="red")
title(ylab="Density", line=2, cex.lab=1)
title(xlab="NRED", line=2, cex.lab=1)
title(" \n\nComparison of random & functional \n pairwise comparisons", outer=TRUE)
mtext(paste("p-value = ",signif(t.test_H_KO_H_random_euclidean$p.value,2)),side=3,col="red",padj=2,cex=.75)
####################################################################################################################

mu_pearson<-mean(I_KOs_Background$random_pairwise_gene_correlation)[1]
sd_pearson<-sd(I_KOs_Background$random_pairwise_gene_correlation)[1]
# A KO distribution
mu_KO_pearson<-mean(I_KOs_Background$KO_pairwise_gene_correlation)[1]
sd_KO_pearson<-sd(I_KOs_Background$KO_pairwise_gene_correlation)[1]
# A maximized background KO distribution
mu_H_KO_pearson<-mean(I_KOs_Background$H_KO_pairwise_gene_correlation)[1]
sd_H_KO_pearson<-sd(I_KOs_Background$H_KO_pairwise_gene_correlation)[1]

# means & standard deviations for calculating Z scores (NRED)
# A completely random background distribution
mu_euclidean<-mean(I_KOs_Background$random_pairwise_gene_euclidean)[1]
sd_euclidean<-sd(I_KOs_Background$random_pairwise_gene_euclidean)[1]
# A KO distribution
mu_KO_euclidean<-mean(I_KOs_Background$KO_pairwise_gene_euclidean)[1]
sd_KO_euclidean<-sd(I_KOs_Background$KO_pairwise_gene_euclidean)[1]
# A maximized background KO distribution
mu_H_KO_euclidean<-mean(I_KOs_Background$H_KO_pairwise_gene_euclidean)[1]
sd_H_KO_euclidean<-sd(I_KOs_Background$H_KO_pairwise_gene_euclidean)[1]

Z_random_pairwise_gene_correlation<-((I_KOs_Background$random_pairwise_gene_correlation)-mu_pearson)/sd_pearson
Z_random_pairwise_gene_euclidean<-((I_KOs_Background$random_pairwise_gene_euclidean)-mu_euclidean)/sd_euclidean

Z_H_KO_pairwise_gene_correlation<-((I_KOs_Background$H_KO_pairwise_gene_correlation)-mu_pearson)/sd_pearson
Z_H_KO_pairwise_gene_euclidean<-((I_KOs_Background$H_KO_pairwise_gene_euclidean)-mu_euclidean)/sd_euclidean

# df_random<-as.data.frame(cbind((Z_random_pairwise_gene_correlation),(Z_random_pairwise_gene_euclidean)))
# df_KO<-as.data.frame(cbind((KO_pairwise_gene_correlation),(KO_pairwise_gene_euclidean)))
# df_H_KO<-as.data.frame(cbind((Z_H_KO_pairwise_gene_euclidean),(H_KO_pairwise_gene_euclidean)))

Z_df_random<-as.data.frame(cbind(-(Z_random_pairwise_gene_correlation),(Z_random_pairwise_gene_euclidean)))
Z_df_H_KO<-as.data.frame(cbind(-(Z_H_KO_pairwise_gene_correlation),(Z_H_KO_pairwise_gene_euclidean)))

rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
plot(hexbin(Z_df_random), colramp=rf, mincnt=1, maxcnt=75,ylab="NRED",xlab="PC")

plot(hexbin(Z_df_H_KO), colramp=rf, mincnt=1, maxcnt=75,ylab="NRED",xlab="PC")
