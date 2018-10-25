##########################################################
### R-codes for outputs from Pathway Enrichment Module ###
##########################################################
## The codes takes a subset of the output data to generate the following results:
## 1) Individual table of pathways that are functionally enriched in the CRC subtypes
## 2) Heatmap of the level of significance of biological pathways up- and down-regulated in each subtype


type1up = read.delim("../Output/CIN_Enrichment_up.txt",as.is=T)
type1down = read.delim("../Output/CIN_Enrichment_down.txt",as.is=T)
type2up = read.delim("../Output/Invasive_Enrichment_up.txt",as.is=T)
type2down = read.delim("../Output/Invasive_Enrichment_down.txt",as.is=T)
type3up = read.delim("../Output/MSIorCIMP_Enrichment_up.txt",as.is=T)
type3down = read.delim("../Output/MSIorCIMP_Enrichment_down.txt",as.is=T)
pathAnnot =read.delim("../../Network Data/Pathways_CPDB_GO_human.txt",as.is=T)


## performing Benjamini-Hochberg method to correct for multiple testing corrections on the p-value ##

ALLPval1 = c(type1up$HypergeoPval, type2up$HypergeoPval, type3up$HypergeoPval)
ALLPval2 = c(type1down$HypergeoPval, type2down$HypergeoPval, type3down$HypergeoPval)
AdjPval1 = p.adjust(ALLPval1,method = "BH")
AdjPval2 = p.adjust(ALLPval2,method = "BH")

a = nrow(type1up)
b = nrow(type1up)+nrow(type2up)
c = nrow(type1up)+nrow(type2up) + nrow(type3up)
type1up$AdjPval = AdjPval1[1:a]
type2up$AdjPval = AdjPval1[(a+1):b]
type3up$AdjPval = AdjPval1[(b+1):c]

a = nrow(type1down)
b = nrow(type1down)+nrow(type2down)
c = nrow(type1down)+nrow(type2down) + nrow(type3down)
type1down$AdjPval = AdjPval2[1:a]
type2down$AdjPval = AdjPval2[(a+1):b]
type3down$AdjPval = AdjPval2[(b+1):c]

type1up$sig=type2up$sig=type3up$sig=type1down$sig=type2down$sig=type3down$sig=""
type1up$sig[which(type1up$AdjPval<0.05& type1up$Enriched_Edgesize>=3)]="*"
type2up$sig[which(type2up$AdjPval<0.05 & type2up$Enriched_Edgesize>=3)]="*"
type3up$sig[which(type3up$AdjPval<0.05& type3up$Enriched_Edgesize>=3)]="*"
type1down$sig[which(type1down$AdjPval<0.05 & type1down$Enriched_Edgesize>=3)]="*"
type2down$sig[which(type2down$AdjPval<0.05 & type2down$Enriched_Edgesize>=3)]="*"
type3down$sig[which(type3down$AdjPval<0.05 & type3down$Enriched_Edgesize>=3)]="*"

type1up = type1up[order(type1up$AdjPval,decreasing=F),]
type2up = type2up[order(type2up$AdjPval,decreasing=F),]
type3up = type3up[order(type3up$AdjPval,decreasing=F),]
type1down = type1down[order(type1down$AdjPval,decreasing=F),]
type2down = type2down[order(type2down$AdjPval,decreasing=F),]
type3down = type3down[order(type3down$AdjPval,decreasing=F),]

sigP_type1up = type1up$Pathway[which(type1up$sig=="*")]
sigP_type2up = type2up$Pathway[which(type2up$sig=="*")]
sigP_type3up = type3up$Pathway[which(type3up$sig=="*")]

sigP_type1down = type1down$Pathway[which(type1down$sig=="*")]
sigP_type2down = type2down$Pathway[which(type2down$sig=="*")]
sigP_type3down = type3down$Pathway[which(type3down$sig=="*")]


ALLP1 = union(sigP_type1up,union(sigP_type2up, sigP_type3up))
ALLP2 = union(sigP_type1down,union(sigP_type2down, sigP_type3down))
ALLP = union(ALLP1,ALLP2)

pval_type1up = type1up$AdjPval[match(ALLP1,type1up$Pathway)]
pval_type2up = type2up$AdjPval[match(ALLP1,type2up$Pathway)]
pval_type3up = type3up$AdjPval[match(ALLP1,type3up$Pathway)]
combup = data.frame(ALLP1,pval_type1up,pval_type2up,pval_type3up)
combup$Annot = pathAnnot$Function[match(ALLP1,pathAnnot$Pathwayid)]
combup = combup[,c(1,5,2:4)]

pval_type1down = type1down$AdjPval[match(ALLP2,type1down$Pathway)]
pval_type2down = type2down$AdjPval[match(ALLP2,type2down$Pathway)]
pval_type3down = type3down$AdjPval[match(ALLP2,type3down$Pathway)]
combdown = data.frame(ALLP2,pval_type1down,pval_type2down,pval_type3down)
combdown$Annot = pathAnnot$Function[match(ALLP2,pathAnnot$Pathwayid)]
combdown = combdown[,c(1,5,2:4)]

colnames(combup)[1] =colnames(combdown)[1] = "Pathway"
comb = merge(combup,combdown,all=T)

colnames(comb)[1] = "Pathway"
type1up= type1up[,c(1:2,9,3:8,10)]
type2up= type2up[,c(1:2,9,3:8,10)]
type3up= type3up[,c(1:2,9,3:8,10)]
type1down= type1down[,c(1:2,9,3:8,10)]
type2down= type2down[,c(1:2,9,3:8,10)]
type3down= type3down[,c(1:2,9,3:8,10)]


write.table(type1up,"CIN_upEnrichmentsummary.txt",sep="\t",row.names=F,quote=F)
write.table(type2up,"Invasive_upEnrichmentsummary.txt",sep="\t",row.names=F,quote=F)
write.table(type3up,"MSIorCIMP_upEnrichmentsummary.txt",sep="\t",row.names=F,quote=F)
write.table(type1down,"CIN_downEnrichmentsummary.txt",sep="\t",row.names=F,quote=F)
write.table(type2down,"Invasive_downEnrichmentsummary.txt",sep="\t",row.names=F,quote=F)
write.table(type3down,"MSIorCIMP_downEnrichmentsummary.txt",sep="\t",row.names=F,quote=F)


## 34 pathways were picked to be represented on heatmap ##

selected = read.csv("../SelectedPathways_CRC.csv",as.is=T,header=T)
combN =comb[which(comb$Pathway%in%selected$Pathway),]

## heatmap ##

PVAL = as.matrix(combN[,c(3:8)])
for(i in 1:ncol(PVAL)){
  cid = which(PVAL[,i]==0)
  if(length(cid)>0) PVAL[cid,i] =0.00001
}


PVALnew =PVAL
PVALnew[,c(1:3)] = -1*log10(PVAL[,c(1:3)])
PVALnew[,c(4:6)] =  log10(PVAL[,c(4:6)])
for(i in 1:ncol(PVALnew)){
  cid = which(is.na(PVALnew[,i]))
  if(length(cid)>0) PVALnew[cid,i] =0
}


library("RColorBrewer")
library("gplots")


pdf("heatmap_CRCPathway_Selected.pdf",height=18,width=10, useDingbats = F)
HM2= heatmap.2(PVALnew, trace="n", Rowv = T, Colv=F,col=bluered(20),breaks=seq(-3,3,by=0.3),
               labRow=combN$Annot,labCol =rep(c("CIN","Invasive","MSI/CIMP"),2), mar=c(6,23), keysize=1,cexCol = 1,cexRow =1.2,
               lhei = c(0.5,4),sepwidth=c(0.01,0.01),sepcolor="black",colsep=1:ncol(PVALnew),rowsep=1:nrow(PVALnew))

dev.off()


combN = combN[rev(HM2$rowInd),]
write.table(combN,"SigPathways_Ord.txt",sep="\t",row.names=F,quote=F,na="")

########################################## end #########################################
