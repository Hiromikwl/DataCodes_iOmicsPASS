############################################################
### R-codes for outputs from Subnetwork Discovery module ###
############################################################
## The codes takes a subset of the output data to generate the following results:
## 1) PCA plot of selected edges, separating the subtypes
## 2) Heatmap of the co-expression data of selected edges across the subtypes

X = read.delim("../Output/Expressiondata_edges.txt",as.is=T,check.names=F)
survivingN = read.delim("../Output/EdgesSelected_minThres.txt",as.is=T)
sampleProb = read.delim("../Output/SampleClass_Probabilities.txt",as.is=T)
subtype = read.delim("../SubtypeInfo_colorectal.txt",as.is=T)

type1_edge = survivingN$Edge[which(survivingN$CIN_absdiknew >0)]
type2_edge = survivingN$Edge[which(survivingN$Invasive_absdiknew >0)]
type3_edge = survivingN$Edge[which(survivingN$MSIorCIMP_absdiknew >0)]
Xsurv = X[which(X$Edge %in% survivingN$Edge),]
row.names(Xsurv) = Xsurv$Edge
Xsurv = Xsurv[,-1]  # 1991

y = subtype$expression_subtype[match(colnames(Xsurv), subtype$patient)]
ycol = as.numeric(factor(y))
yCOL = rep("maroon",length(ycol))
yCOL[which(ycol==1)] = "green"
yCOL[which(ycol==2)] = "purple"

rowMn = apply(Xsurv,1,mean)
XsurvN = sweep(Xsurv,1,rowMn)

library("gplots")

pdf("heatmap_EdgesSurv.pdf",height=16,width=14, useDingbats = F)
HM= heatmap.2(as.matrix(XsurvN), trace="n", Rowv = T, Colv=T,col=bluered(20),breaks=seq(-5,5,by=0.5),ColSideColors  = yCOL,
              labRow=F,labCol =F, mar=c(6,18), keysize=1)

par(lend = 1)
legend("right",      # location of the legend on the heatmap plot
       legend = c("CIN","Invasive","MSI/CIMP"), # category labels
       col =unique(yCOL),  # color key
       lty= 1,             # line style
       lwd = 10   )         # line width
dev.off()

######################################################################
## creating heatmap for selected edges separately for each subtype ##
Xsurv_type1 = XsurvN[rev(HM$rowInd),which(y=="CIN")]
Xsurv_type2 = XsurvN[rev(HM$rowInd),which(y=="Invasive")]
Xsurv_type3 = XsurvN[rev(HM$rowInd),which(y=="MSI/CIMP")]

Rowbar1 = Rowbar2 =Rowbar3 =rep("grey",nrow(Xsurv_type1))
Rowbar1[which(row.names(Xsurv_type1) %in% type1_edge)] = "cyan"
Rowbar2[which(row.names(Xsurv_type2) %in% type2_edge)] = "cyan"
Rowbar3[which(row.names(Xsurv_type3) %in% type3_edge)] = "cyan"

tmpB = sampleProb[match(colnames(Xsurv_type1), sampleProb$Subject),]
ClassP_type1 = as.numeric(tmpB$TrueClass==tmpB$PredictedClass)  
yCOL = rep("grey",nrow(tmpB))
yCOL[which(ClassP_type1==0)] = "brown"

CC = tmpB$TrueClass
CC[which(ClassP_type1==0)] = tmpB$PredictedClass[which(ClassP_type1==0)]

pdf("heatmap_EdgesSurv_CIN.pdf",height=15,width=6, useDingbats = F)
HM= heatmap.2(as.matrix(Xsurv_type1), trace="n", Rowv = F, Colv=T,col=bluered(20),breaks=seq(-5,5,by=0.5),RowSideColors  = Rowbar1,ColSideColors  = yCOL,
              labRow=F,labCol =CC, mar=c(6,6), keysize=1)

dev.off()

tmpB = sampleProb[match(colnames(Xsurv_type2), sampleProb$Subject),]
ClassP_type2 = as.numeric(tmpB$TrueClass==tmpB$PredictedClass)  
yCOL = rep("grey",nrow(tmpB))
yCOL[which(ClassP_type2==0)] = "brown"

CC = tmpB$TrueClass
CC[which(ClassP_type2==0)] = tmpB$PredictedClass[which(ClassP_type2==0)]

pdf("heatmap_EdgesSurv_Invasive.pdf",height=15,width=6, useDingbats = F)
HM= heatmap.2(as.matrix(Xsurv_type2), trace="n", Rowv = F, Colv=T,col=bluered(20),breaks=seq(-5,5,by=0.5),RowSideColors  = Rowbar2,ColSideColors  = yCOL,
              labRow=F,labCol =CC, mar=c(6,6), keysize=1)

dev.off()

tmpB = sampleProb[match(colnames(Xsurv_type3), sampleProb$Subject),]
ClassP_type3 = as.numeric(tmpB$TrueClass==tmpB$PredictedClass)  
yCOL = rep("grey",nrow(tmpB))
yCOL[which(ClassP_type3==0)] = "brown"

CC = tmpB$TrueClass
CC[which(ClassP_type3==0)] = tmpB$PredictedClass[which(ClassP_type3==0)]
pdf("heatmap_EdgesSurv_MSICIMP.pdf",height=15,width=6, useDingbats = F)
HM= heatmap.2(as.matrix(Xsurv_type3), trace="n", Rowv = F, Colv=T,col=bluered(20),breaks=seq(-5,5,by=0.5),RowSideColors  = Rowbar3,ColSideColors  = yCOL,
              labRow=F,labCol =CC, mar=c(6,6), keysize=1)

dev.off()


#######################################################################################
### performing Principal component analysis (PCA) on selected edges ##

Y = y
Xnew = Xsurv

outPCA = prcomp(t(as.matrix(Xnew)),center=T,scale =T)   ## first PC 24.6%, 11.8%, 8.8%

scores = as.data.frame(outPCA$x)
PC1 = scores[,1]
PC2 = scores[,2]
PC3 = scores[,3]
ColLab = Y

tmp_col =c("cyan","royalblue","maroon")
samples = Y
uniqS = unique(Y)
colvec =rep(tmp_col[1],ncol(X))
colvec[which(samples==uniqS[2])] = tmp_col[2]
colvec[which(samples==uniqS[3])] = tmp_col[3]


pchCol =rep(16,ncol(X))
pchCol[which(samples==uniqS[2])] = 17
pchCol[which(samples==uniqS[3])] = 3

sampleProb = sampleProb[match(colnames(Xsurv),sampleProb$Subject),]
ClassProb = rep(NA,nrow(sampleProb))

for(i in 1:nrow(sampleProb)){
  if(sampleProb$TrueClass[i]=="CIN") ClassProb[i] = sampleProb$Prob_CIN[i]
  if(sampleProb$TrueClass[i]=="Invasive") ClassProb[i] = sampleProb$Prob_Invasive[i]
  if(sampleProb$TrueClass[i]=="MSIorCIMP") ClassProb[i] = sampleProb$Prob_MSIorCIMP[i]
}

### highlighting subjects who has class probabilities < 0.4 in red color ###
oddS = colnames(Xsurv)[which(ClassProb<0.4)]


pdf("PCAplot_CRCclustering.pdf",height=6,width=6, useDingbats = F)
par(oma=c(2,1,2,1),mar=c(4,4,2,1))
plot(PC1,PC2,cex =0.7,type="n",xlab="PC1 (Variance explained=24.6%)",ylab="PC2 (Variance explained=11.8%)",cex.lab = 0.8)
points(PC1,PC2,pch=pchCol,col=colvec)
abline(h=0,lty=2)
abline(v=0,lty=2)
par(oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0),new=TRUE)
legend(-15,35, uniqS, xpd = TRUE, horiz = TRUE, inset = c(0, 0),
       bty = "n", pch = c(16,17,3),col=tmp_col,cex=0.8)

dev.off()


pdf("PCAplot_CRCclustering_WithClassProb.pdf",height=6,width=6, useDingbats = F)
par(oma=c(2,1,2,1),mar=c(4,4,2,1))
plot(PC1,PC2,cex =0.7,type="n",xlab="PC1 (Variance explained=24.6%)",ylab="PC2 (Variance explained=11.8%)",cex.lab = 0.8)
points(PC1,PC2,pch=pchCol,col=colvec)
points(PC1[match(oddS,colnames(Xsurv))],PC2[match(oddS,colnames(Xsurv))]  ,pch=pchCol[match(oddS,colnames(Xsurv))],col=2)
abline(h=0,lty=2)
abline(v=0,lty=2)
par(oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0),new=TRUE)
legend(-15,35, uniqS, xpd = TRUE, horiz = TRUE, inset = c(0, 0),
       bty = "n", pch = c(16,17,3),col=tmp_col,cex=0.8)
text(PC1[match(oddS,colnames(Xsurv))],PC2[match(oddS,colnames(Xsurv))], labels =sampleProb$PredictedClass[match(oddS,colnames(Xsurv))],cex=0.7)
text(PC1[match(oddS,colnames(Xsurv))],PC2[match(oddS,colnames(Xsurv))], labels =oddS,cex=0.6,pos = 3)

dev.off()

########################################## end #########################################


