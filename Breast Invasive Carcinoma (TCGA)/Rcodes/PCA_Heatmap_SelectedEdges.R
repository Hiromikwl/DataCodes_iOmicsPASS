############################################################
### R-codes for outputs from Subnetwork Discovery module ###
############################################################
## The codes takes a subset of the output data to generate the following results:
## 1) PCA plot of selected edges, separating the subtypes
## 2) Heatmap of the co-expression data of selected edges across the subtypes


X = read.delim("../Output/Expressiondata_edges.txt",as.is=T,check.names=F)
survivingN = read.delim("../Output/EdgesSelected_minThres.txt",as.is=T)
sampleProb = read.delim("../Output/SampleClass_Probabilities.txt",as.is=T)
subtype = read.delim("../TCGA_BRCA_subtypes.txt",as.is=T)

basal_edge = survivingN$Edge[which(survivingN$Basal.like_absdiknew>0)]
her2_edge = survivingN$Edge[which(survivingN$HER2.enriched_absdiknew>0)]
lumA_edge = survivingN$Edge[which(survivingN$Luminal.A_absdiknew>0)]
lumB_edge = survivingN$Edge[which(survivingN$Luminal.B_absdiknew>0)]
Xsurv = X[which(X$Edge %in% survivingN$Edge),]
row.names(Xsurv) = Xsurv$Edge
Xsurv = Xsurv[,-1]

basal_sub = survivingN[which(survivingN$Basal.like_absdiknew>0),c(1:8)]
her2_sub = survivingN[which(survivingN$HER2.enriched_absdiknew>0),c(1:4,9:12)]
lumA_sub = survivingN[which(survivingN$Luminal.A_absdiknew>0),c(1:4,13:16)]
lumB_sub = survivingN[which(survivingN$Luminal.B_absdiknew>0),c(1:4,17:20)]


y = subtype$PAM50[match(colnames(Xsurv), subtype$SubjectID)]
ycol = as.numeric(factor(y))
yCOL = rep("purple",length(ycol))
yCOL[which(ycol==1)] = "blue"
yCOL[which(ycol==2)] = "maroon"
yCOL[which(ycol==3)] = "green"


library("gplots")

pdf("heatmap_EdgesSurv.pdf",height=16,width=10, useDingbats = F)
HM= heatmap.2(as.matrix(Xsurv), trace="n", Rowv = T, Colv=T,col=bluered(20),breaks=seq(-5,5,by=0.5),ColSideColors  = yCOL,
              labRow=F,labCol =F, mar=c(6,18), keysize=1)

par(lend = 1)
legend("right",      # location of the legend on the heatmap plot
       legend = c("Basal","HER2E","Lum A","Lum B"), # category labels
       col =unique(yCOL),  # color key
       lty= 1,             # line style
       lwd = 10   )         # line width
dev.off()

######################################################################
## creating heatmap for selected edges separately for each subtype ##

Xsurv_basal = Xsurv[rev(HM$rowInd),which(y=="Basal-like")]
Xsurv_her2 = Xsurv[rev(HM$rowInd),which(y=="HER2-enriched")]
Xsurv_lumA = Xsurv[rev(HM$rowInd),which(y=="Luminal A")]
Xsurv_lumB = Xsurv[rev(HM$rowInd),which(y=="Luminal B")]

tmpB = sampleProb[match(colnames(Xsurv_basal), sampleProb$Subject),]
ClassP_basal = as.numeric(tmpB$TrueClass==tmpB$PredictedClass)  
yCOL = rep("grey",nrow(tmpB))
yCOL[which(ClassP_basal==0)] = "brown"

xCOL=rep("grey",nrow(Xsurv_basal))
xCOL[which(row.names(Xsurv_basal) %in% basal_edge)] = "cyan"

CC = tmpB$TrueClass
CC[which(ClassP_basal==0)] = tmpB$PredictedClass[which(ClassP_basal==0)]

pdf("heatmap_EdgesSurv_Basal.pdf",height=16,width=8, useDingbats = F)
HM_basal= heatmap.2(as.matrix(Xsurv_basal), trace="n", Rowv = F, Colv=T,col=bluered(20),breaks=seq(-5,5,by=0.5),RowSideColors  = xCOL,ColSideColors  = yCOL,
              labRow=F,labCol =CC, mar=c(6,6), keysize=1)
dev.off()

tmpB = sampleProb[match(colnames(Xsurv_her2), sampleProb$Subject),]
ClassP_her2 = as.numeric(tmpB$TrueClass==tmpB$PredictedClass)  
yCOL = rep("grey",nrow(tmpB))
yCOL[which(ClassP_her2==0)] = "brown"

xCOL=rep("grey",nrow(Xsurv_her2))
xCOL[which(row.names(Xsurv_her2) %in% her2_edge)] = "cyan"
CC = tmpB$TrueClass
CC[which(ClassP_her2==0)] = tmpB$PredictedClass[which(ClassP_her2==0)]


pdf("heatmap_EdgesSurv_HER2E.pdf",height=16,width=8, useDingbats = F)
HM_her2 = heatmap.2(as.matrix(Xsurv_her2), trace="n", Rowv = F, Colv=T,col=bluered(20),breaks=seq(-5,5,by=0.5),RowSideColors  = xCOL,ColSideColors  = yCOL,
              labRow=F,labCol =CC, mar=c(6,6), keysize=1)
dev.off()


tmpB = sampleProb[match(colnames(Xsurv_lumA), sampleProb$Subject),]
ClassP_lumA = as.numeric(tmpB$TrueClass==tmpB$PredictedClass)  
yCOL = rep("grey",nrow(tmpB))
yCOL[which(ClassP_lumA==0)] = "brown"

xCOL=rep("grey",nrow(Xsurv_lumA))
xCOL[which(row.names(Xsurv_lumA) %in% lumA_edge)] = "cyan"
CC = tmpB$TrueClass
CC[which(ClassP_lumA==0)] = tmpB$PredictedClass[which(ClassP_lumA==0)]

pdf("heatmap_EdgesSurv_LumA.pdf",height=16,width=8, useDingbats = F)
HM_lumA= heatmap.2(as.matrix(Xsurv_lumA), trace="n", Rowv = F, Colv=T,col=bluered(20),breaks=seq(-5,5,by=0.5),RowSideColors  = xCOL,ColSideColors  = yCOL,
              labRow=F,labCol =CC, mar=c(6,6), keysize=1)
dev.off()


tmpB = sampleProb[match(colnames(Xsurv_lumB), sampleProb$Subject),]
ClassP_lumB = as.numeric(tmpB$TrueClass==tmpB$PredictedClass)  
yCOL = rep("grey",nrow(tmpB))
yCOL[which(ClassP_lumB==0)] = "brown"

xCOL=rep("grey",nrow(Xsurv_lumB))
xCOL[which(row.names(Xsurv_lumB) %in% lumB_edge)] = "cyan"
CC = tmpB$TrueClass
CC[which(ClassP_lumB==0)] = tmpB$PredictedClass[which(ClassP_lumB==0)]

pdf("heatmap_EdgesSurv_LumB.pdf",height=16,width=8, useDingbats = F)
HM_lumB = heatmap.2(as.matrix(Xsurv_lumB), trace="n", Rowv = F, Colv=T,col=bluered(20),breaks=seq(-5,5,by=0.5),RowSideColors  = xCOL,ColSideColors  = yCOL,
              labRow=F,labCol =CC, mar=c(6,6), keysize=1)
dev.off()

#######################################################################################
### performing Principal component analysis (PCA) on selected edges ##

row.names(X) = X$Edge
X = X[,-1]
Y = subtype$PAM50[match(colnames(X),subtype$SubjectID)]

Xnew = X[which(row.names(X) %in% survivingN$Edge),]


outPCA = prcomp(t(as.matrix(Xnew)),center=T,scale =T) ## first PC 55.3%, 6.7%, 4.9%

scores = as.data.frame(outPCA$x)
PC1 = scores[,1]
PC2 = scores[,2]
PC3 = scores[,3]
ColLab = Y


tmp_col =c("royalblue","maroon","green","purple")
samples = Y
uniqS = unique(Y)
colvec =rep(tmp_col[1],ncol(X))
colvec[which(samples==uniqS[2])] = tmp_col[2]
colvec[which(samples==uniqS[3])] = tmp_col[3]
colvec[which(samples==uniqS[4])] = tmp_col[4]

pchCol =rep(16,ncol(X))
pchCol[which(samples==uniqS[2])] = 17
pchCol[which(samples==uniqS[3])] = 3
pchCol[which(samples==uniqS[4])] = 4


sampleProb = sampleProb[match(colnames(Xsurv),sampleProb$Subject),]

ClassProb = rep(NA,nrow(sampleProb))
for(i in 1:nrow(sampleProb)){
  if(sampleProb$TrueClass[i]==uniqS[1]) ClassProb[i] = sampleProb$Prob_Basal.like[i]
  if(sampleProb$TrueClass[i]==uniqS[2]) ClassProb[i] = sampleProb$Prob_HER2.enriched[i]
  if(sampleProb$TrueClass[i]==uniqS[3]) ClassProb[i] = sampleProb$Prob_Luminal.A[i]
  if(sampleProb$TrueClass[i]==uniqS[4]) ClassProb[i] = sampleProb$Prob_Luminal.B[i]
}


### highlighting subjects who has class probabilities < 0.4 in red color ###
oddS = colnames(Xsurv)[which(ClassProb<0.4)]


pdf("PCAplot_BRCAclustering.pdf",height=6,width=12, useDingbats = F)
par(oma=c(2,1,2,1),mar=c(4,4,2,1),mfrow=c(1,2))
plot(PC1,PC2,cex =0.7,type="n",ylim=c(-50,50),xlim=c(-100,100),xlab="PC1 (Variance explained=55.3%)",ylab="PC2 (Variance explained=6.7%)",cex.lab = 0.8)
points(PC1,PC2,pch=pchCol,col=colvec)
abline(h=0,lty=2)
abline(v=0,lty=2)
legend(-100,65, uniqS, xpd = TRUE, horiz = TRUE, inset = c(0, 0),
       bty = "n", pch = c(16,17,3,4),col=tmp_col,cex=0.8)
plot(PC1,PC3,cex =0.7,type="n",ylim=c(-50,50),xlim=c(-100,100),xlab="PC1 (Variance explained=55.3%)",ylab="PC3 (Variance explained=4.9%)",cex.lab = 0.8)
points(PC1,PC3,pch=pchCol,col=colvec)
abline(h=0,lty=2)
abline(v=0,lty=2)
legend(-100,65, uniqS, xpd = TRUE, horiz = TRUE, inset = c(0, 0),
       bty = "n", pch = c(16,17,3,4),col=tmp_col,cex=0.8)
dev.off()

pdf("PCAplot_BRCAclustering_withClassProb.pdf",height=6,width=12, useDingbats = F)
par(oma=c(2,1,2,1),mar=c(4,4,2,1),mfrow=c(1,2))
plot(PC1,PC2,cex =0.7,type="n",ylim=c(-50,50),xlim=c(-100,100),xlab="PC1 (Variance explained=55.3%)",ylab="PC2 (Variance explained=6.7%)",cex.lab = 0.8)
points(PC1,PC2,pch=pchCol,col=colvec)
points(PC1[match(oddS,colnames(Xsurv))],PC2[match(oddS,colnames(Xsurv))]  ,pch=pchCol[match(oddS,colnames(Xsurv))],col=2)
abline(h=0,lty=2)
abline(v=0,lty=2)
text(PC1[match(oddS,colnames(Xsurv))],PC2[match(oddS,colnames(Xsurv))], labels =sampleProb$PredictedClass[match(oddS,colnames(Xsurv))],cex=0.5)
legend(-100,65, uniqS, xpd = TRUE, horiz = TRUE, inset = c(0, 0),
       bty = "n", pch = c(16,17,3,4),col=tmp_col,cex=0.8)
plot(PC1,PC3,cex =0.7,type="n",ylim=c(-50,50),xlim=c(-100,100),xlab="PC1 (Variance explained=55.3%)",ylab="PC3 (Variance explained=4.9%)",cex.lab = 0.8)
points(PC1,PC3,pch=pchCol,col=colvec)
points(PC1[match(oddS,colnames(Xsurv))],PC3[match(oddS,colnames(Xsurv))] ,pch=pchCol[match(oddS,colnames(Xsurv))],col=2)
abline(h=0,lty=2)
abline(v=0,lty=2)
legend(-100,65, uniqS, xpd = TRUE, horiz = TRUE, inset = c(0, 0),
       bty = "n", pch = c(16,17,3,4),col=tmp_col,cex=0.8)
text(PC1[match(oddS,colnames(Xsurv))],PC3[match(oddS,colnames(Xsurv))], labels =sampleProb$PredictedClass[match(oddS,colnames(Xsurv))],cex=0.5)
dev.off()





