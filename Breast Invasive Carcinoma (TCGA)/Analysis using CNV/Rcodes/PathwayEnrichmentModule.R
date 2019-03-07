##########################################################
### R-codes for outputs from Pathway Enrichment Module ###
##########################################################
## The codes takes a subset of the output data to generate the following results:
## 1) Individual table of pathways that are functionally enriched in the BRCA subtypes
## 2) Heatmap of the level of significance of biological pathways up- and down-regulated in each subtype


basalup = read.delim("../Output/Basal-like_Enrichment_up.txt",as.is=T)
basaldown = read.delim("../Output/Basal-like_Enrichment_down.txt",as.is=T)
her2up = read.delim("../Output/HER2-enriched_Enrichment_up.txt",as.is=T)
her2down = read.delim("../Output/HER2-enriched_Enrichment_down.txt",as.is=T)
lumAup = read.delim("../Output/Luminal A_Enrichment_up.txt",as.is=T)
lumAdown = read.delim("../Output/Luminal A_Enrichment_down.txt",as.is=T)
lumBup = read.delim("../Output/Luminal B_Enrichment_up.txt",as.is=T)
lumBdown = read.delim("../Output/Luminal B_Enrichment_down.txt",as.is=T)
pathAnnot =read.delim("../../../Network Data/Pathways_CPDB_GO_human.txt",as.is=T)

## no pathways were down-regulated in Her2 subtype, we will exclude it for subsequent analysis ##

## performing Benjamini-Hochberg method to correct for multiple testing corrections on the p-value ##
ALLPval1 = c(basalup$HypergeoPval, her2up$HypergeoPval, lumAup$HypergeoPval, lumBup$HypergeoPval)
ALLPval2 = c(basaldown$HypergeoPval, lumAdown$HypergeoPval, lumBdown$HypergeoPval)
AdjPval1 = p.adjust(ALLPval1,method = "BH")
AdjPval2 = p.adjust(ALLPval2,method = "BH")

a = nrow(basalup)
b = nrow(basalup)+nrow(her2up)
c = nrow(basalup)+nrow(her2up) + nrow(lumAup)
d = nrow(basalup)+nrow(her2up) + nrow(lumAup) +nrow(lumBup)
basalup$AdjPval = AdjPval1[1:a]
her2up$AdjPval = AdjPval1[(a+1):b]
lumAup$AdjPval = AdjPval1[(b+1):c]
lumBup$AdjPval = AdjPval1[(c+1):d]

a = nrow(basaldown)
b = nrow(basaldown)+nrow(her2down)
c = nrow(basaldown)+nrow(her2down) + nrow(lumAdown)
d = nrow(basaldown)+nrow(her2down) + nrow(lumAdown) +nrow(lumBdown)
basaldown$AdjPval = AdjPval2[1:a]
#her2down$AdjPval = AdjPval2[(a+1):b]
lumAdown$AdjPval = AdjPval2[(b+1):c]
lumBdown$AdjPval = AdjPval2[(c+1):d]

basalup$sig=her2up$sig=lumAup$sig=lumBup$sig=basaldown$sig=lumAdown$sig=lumBdown$sig=""
basalup$sig[which(basalup$AdjPval<0.05 & basalup$Enriched_Edgesize>=3)]="*"
her2up$sig[which(her2up$AdjPval<0.05 & her2up$Enriched_Edgesize>=3)]="*"
lumAup$sig[which(lumAup$AdjPval<0.05 & lumAup$Enriched_Edgesize>=3)]="*"
lumBup$sig[which(lumBup$AdjPval<0.05 & lumBup$Enriched_Edgesize>=3)]="*"
basaldown$sig[which(basaldown$AdjPval<0.05 & basaldown$Enriched_Edgesize>=3)]="*"
#her2down$sig[which(her2down$AdjPval<0.05)]="*"
lumAdown$sig[which(lumAdown$AdjPval<0.05 & lumAdown$Enriched_Edgesize>=3)]="*"
lumBdown$sig[which(lumBdown$AdjPval<0.05 & lumBdown$Enriched_Edgesize>=3)]="*"

basalup = basalup[order(basalup$AdjPval,decreasing=F),]
her2up = her2up[order(her2up$AdjPval,decreasing=F),]
lumAup = lumAup[order(lumAup$AdjPval,decreasing=F),]
lumBup = lumBup[order(lumBup$AdjPval,decreasing=F),]
basaldown = basaldown[order(basaldown$AdjPval,decreasing=F),]
#her2down = her2down[order(her2down$AdjPval,decreasing=F),]
lumAdown = lumAdown[order(lumAdown$AdjPval,decreasing=F),]
lumBdown = lumBdown[order(lumBdown$AdjPval,decreasing=F),]


sigP_basalup = basalup$Pathway[which(basalup$sig=="*")]
sigP_her2up = her2up$Pathway[which(her2up$sig=="*")]
sigP_lumAup = lumAup$Pathway[which(lumAup$sig=="*")]
sigP_lumBup = lumBup$Pathway[which(lumBup$sig=="*")]
sigP_basaldown = basaldown$Pathway[which(basaldown$sig=="*")]
#sigP_her2down = her2down$Pathway[which(her2down$sig=="*")]
sigP_lumAdown = lumAdown$Pathway[which(lumAdown$sig=="*")]
sigP_lumBdown = lumBdown$Pathway[which(lumBdown$sig=="*")]

ALLPup = union(sigP_basalup,union(sigP_her2up, union(sigP_lumAup,sigP_lumBup)))
ALLPdown = union(sigP_basaldown, union(sigP_lumAdown,sigP_lumBdown))

pval_basalup = basalup$AdjPval[match(ALLPup,basalup$Pathway)]
pval_her2up = her2up$AdjPval[match(ALLPup,her2up$Pathway)]
pval_lumAup = lumAup$AdjPval[match(ALLPup,lumAup$Pathway)]
pval_lumBup = lumBup$AdjPval[match(ALLPup,lumBup$Pathway)]

pval_basaldown = basaldown$AdjPval[match(ALLPdown,basaldown$Pathway)]
#pval_her2down = her2down$AdjPval[match(ALLPdown,her2down$Pathway)]
pval_lumAdown = lumAdown$AdjPval[match(ALLPdown,lumAdown$Pathway)]
pval_lumBdown = lumBdown$AdjPval[match(ALLPdown,lumBdown$Pathway)]

basalup= basalup[,c(1,2,9,3:8,10)]
her2up= her2up[,c(1,2,9,3:8,10)]
lumAup= lumAup[,c(1,2,9,3:8,10)]
lumBup= lumBup[,c(1,2,9,3:8,10)]
basaldown= basaldown[,c(1,2,9,3:8,10)]
#her2down= her2down[,c(1,8,2,9,3:7,10)]
lumAdown= lumAdown[,c(1,2,9,3:8,10)]
lumBdown= lumBdown[,c(1,2,9,3:8,10)]

## output enrichment tables for each subtype ##

write.table(basalup,"Basal_upEnrichmentsummary.txt",sep="\t",row.names=F,quote=F)
write.table(her2up,"Her2_upEnrichmentsummary.txt",sep="\t",row.names=F,quote=F)
write.table(lumAup,"LumA_upEnrichmentsummary.txt",sep="\t",row.names=F,quote=F)
write.table(lumBup,"LumB_upEnrichmentsummary.txt",sep="\t",row.names=F,quote=F)
write.table(basaldown,"Basal_downEnrichmentsummary.txt",sep="\t",row.names=F,quote=F)
#write.table(her2down,"Enrichment/Her2_downEnrichmentsummary.txt",sep="\t",row.names=F,quote=F)
write.table(lumAdown,"LumA_downEnrichmentsummary.txt",sep="\t",row.names=F,quote=F)
write.table(lumBdown,"LumB_downEnrichmentsummary.txt",sep="\t",row.names=F,quote=F)


## creating heatmap figure ##
ALLPaths = union(ALLPup, ALLPdown)

pval_basalup = basalup$AdjPval[match(ALLPaths,basalup$Pathway)]
pval_her2up = her2up$AdjPval[match(ALLPaths,her2up$Pathway)]
pval_lumAup = lumAup$AdjPval[match(ALLPaths,lumAup$Pathway)]
pval_lumBup = lumBup$AdjPval[match(ALLPaths,lumBup$Pathway)]

pval_basaldown = basaldown$AdjPval[match(ALLPaths,basaldown$Pathway)]
#pval_her2down = her2down$AdjPval[match(ALLPdown,her2down$Pathway)]
pval_lumAdown = lumAdown$AdjPval[match(ALLPaths,lumAdown$Pathway)]
pval_lumBdown = lumBdown$AdjPval[match(ALLPaths,lumBdown$Pathway)]
COMB = data.frame(ALLPaths,pval_basalup,pval_her2up,pval_lumAup,pval_lumBup, pval_basaldown, pval_lumAdown, pval_lumBdown)
COMB$Annot = pathAnnot$Function[match(ALLPaths,pathAnnot$Pathwayid)]
COMB = COMB[,c(1,9,2:8)]

comb_new = COMB

## 44 pathways were picked to be represented on heatmap ##
list = read.csv("../SelectedPathways_BRCA.csv",as.is=T,header=T)

ppp = unique(list$Pathway)
comb_new = COMB[which(COMB$ALLPaths %in%ppp),]

## remove those sig p-values in subtuypes with too few sig edges ##
for(i in 3:9){
  #cid = which(comb_new[,i]<0.05)
  #if(length(cid)>0){
  check = TRUE
  if(i==3) dd = basalup
  if(i==4) dd = her2up
  if(i==5) dd = lumAup
  if(i==6) dd = lumBup
  if(i==7) dd = basaldown
  if(i==8) dd = lumAdown
  if(i==9) dd = lumBdown
  
  tmp = dd[match(comb_new$ALLPaths, dd$Pathway),]
  if(!all(tmp$sig=="*")) check = FALSE
  if(!check){
    ind = which(tmp$sig!="*")
    comb_new[ind,i] = NA
  }
  #}
}




comb_new$pval_her2down = NA
comb_new = comb_new[,c(1:7,10,8:9)]


PVAL = as.matrix(comb_new[,c(3:10)])
for(i in 1:ncol(PVAL)){
  cid = which(PVAL[,i]==0)
  if(length(cid)>0) PVAL[cid,i] =0.00001
}

PVALnew =PVAL
PVALnew[,c(1:4)] = -1*log10(PVAL[,c(1:4)])
PVALnew[,c(5:8)] =  log10(PVAL[,c(5:8)])
for(i in 1:ncol(PVALnew)){
  cid = which(is.na(PVALnew[,i]))
  if(length(cid)>0) PVALnew[cid,i] =0
}

library("RColorBrewer")
library("gplots")

pdf("heatmap_BRCAPathway_selectedPathways.pdf",height=15,width=8, useDingbats = F)
HM2= heatmap.2(PVALnew, trace="n", Rowv = T, Colv=F,col=bluered(20),breaks=seq(-3,3,by=0.3),
              labRow=comb_new$Annot,labCol =rep(c("Basal-like","HER2E","Luminal A","Luminal B"),2), mar=c(6,24), keysize=0.8,cexCol = 1,cexRow =1.2,
              lhei = c(0.5,4),sepwidth=c(0.01,0.01),sepcolor="black",colsep=1:ncol(PVALnew),rowsep=1:nrow(PVALnew))
dev.off()

combN = comb_new[rev(HM2$rowInd),]
write.table(combN,"SigPathways_BRCA_Ord.txt",sep="\t",row.names=F,quote=F,na="")

########################################## end #########################################
