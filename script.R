library(survival)
library(survminer)
library(ggplot2)
head(data)
fit <- survfit(Surv(time, event) ~ group, data = data)
print(fit)
fit_cox <- coxph(Surv(time, event) ~ group, data = data)
print(fit_cox)
cox.zph(fit_cox)
ggsurvplot(fit = fit, data = data, fun = "pct",
           palette = c("#00468B", "#ED0000", "#42B540", "#0099B4", "#925E9F"),
           linetype = 1, pval = TRUE, 
           censor = TRUE, censor.size = 7,
           risk.table = FALSE, conf.int = FALSE)



library(DESeq2)
library(limma)
library(ggplot2)
library(ggpubr)
library(biomaRt)
library(ggrepel)
library(ggsci)
cleanCount <- function(rawdat,
                       lowlimit=5) {
  rawdat <- data.frame(rawdat)
  rawdat <- avereps(rawdat[,-c(1)],ID=rawdat[,1])
  n0 <- apply(rawdat<=lowlimit,1,sum)
  i0 <- which(n0 > ncol(rawdat)/2)
  rawdat <- rawdat[-i0,]
  return(rawdat)
}
##DESeq2
DEA <- function(rawdat,
                org='org.Hs.eg.db',
                type,
                group,
                N,M,
                volcanoPlot=TRUE,
                batch_compare=FALSE,
                anno) {
  colDate=data.frame(row.names=colnames(rawdat),
                     condition=factor(group,levels=c(N,M)))
  rawdat <- round(rawdat,0)
  dds<-DESeqDataSetFromMatrix(countData =rawdat,colData = colDate,design = ~condition)
  dds<-DESeq(dds)
  sizeFactors(dds)
  ##export results
  res<-results(dds,contrast = c('condition',M,N))
  res<-as.data.frame(res)
  res<-cbind(rownames(res),res) 
  res <- na.omit(res)
  colnames(res)<- c(type,"baseMean" ,"log2FoldChange","lfcSE" ,"stat","pvalue","padj" )
  res$direciton=ifelse(res$pvalue < 0.05 & abs(res$log2FoldChange) > 1,
                       ifelse(res$log2FoldChange>1, 'Up', 'Down'),'NOT')
  ##boxplot for normalization
  if(batch_compare==TRUE){
    rld <- rlogTransformation(dds) 
    exprSet_new=assay(rld)
    par(cex = 0.7)
    n.sample=ncol(rawdat)
    if(n.sample>40) par(cex = 0.5)
    cols <- rainbow(n.sample*1.2)
    par(mfrow=c(2,2))
    boxplot(rawdat, col = cols,main="expression value",las=2)
    boxplot(exprSet_new, col = cols,main="expression value",las=2)
    hist(unlist(rawdat))
    hist(exprSet_new)
  }
  ##transfom ID
  if(type!='SYMBOL'){
    gene.df <- bitr(res[,1],
                    fromType = type,
                    toType = c("SYMBOL"),
                    OrgDb = org)
    res <- merge(gene.df,res,by=type)
  }
  ##volcanoplot
  if(volcanoPlot==TRUE){
    p1=ggplot(data = res, aes(x = log2FoldChange, y = -log10(pvalue),colour = direciton,fill=direciton))+
      geom_point(alpha = 0.5,size=5)+
      scale_colour_manual(values = c("#4DBBD5", "grey", "#E64B35"))+ 
      xlab('log2(Fold change)')+ylab('-log10(P Value)')+ 
      geom_hline(yintercept = -log10(0.05),lty=2,lwd=1,alpha=0.8)+ 
      geom_vline(xintercept = c(1,-1),lty=2,lwd=1,alpha=0.8)+ 
      theme_bw()+
      theme(panel.border=element_rect(size = 1),
            plot.background = element_rect(size=1),
            axis.line = element_line(colour = "black"),
            axis.text = element_text(size = 15, colour = "black"), 
            axis.title = element_text(size = 18),
            axis.ticks = element_line(size=1),
            legend.position = "none")
  }
  if(is.na(anno)==FALSE){
    res$anno_name=res$SYMBOL
    res$anno_name[which(res$anno_name!=anno)]=NA
    anno_name=data.frame(res[which(res$SYMBOL==anno),])
    p1=p1+geom_point(data = anno_name,aes(x = log2FoldChange, y = -log10(pvalue)),
                     shape=1,colour='black',size=5,stroke=1,na.rm=TRUE)+
      geom_text_repel(data=res,aes(x = log2FoldChange, y = -log10(pvalue),label = anno_name),size = 4,
                      box.padding = unit(0.6,"lines"),
                      segment.color = "black",segment.size=0.7,
                      show.legend = FALSE,
                      na.rm = TRUE,max.overlaps = 1000000000,colour = 'black')
  }
  ##export
  if(batch_compare==TRUE){  
    return(list(deres=res,
                normalized_martrix=data.frame(exprSet_new),
                volcano=p1,
                dds=dds))
  }else{
    return(list(deres=res,volcano=p1,dds=dds))
  }
}


exp_count <- fread("/RR/count.csv",header = T)
exp_count <- data.frame(exp_count)
exp_count <- cleanCount(exp_count,lowlimit = 0)
group <- ifelse(exp_count$SEC14L3>=median(exp1$DLD),'H','L')
identical(colnames(exp_count),names(group))
res <- DEA(rawdat = exp_count,type = 'SYMBOL',group=group,N='L',M='H',anno = NA)


library(data.table)
library(GseaVis)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
list_ge <- read.gmt('D:/RR/c5.go.bp.v2023.2.Hs.symbols.gmt')
list_ke <- read.gmt('D:/RR/c2.cp.kegg_medicus.v2023.2.Hs.symbols.gmt')
list_ke2 <- read.gmt('D:/RR/c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt')
list_copper <- list_ge[grep("COPPER",list_ge$term),]
delist <- fread("sec.txt")
delist <- delist[!duplicated(delist$gene_name),]
genelist <- delist$log2FoldChange
names(genelist) <- delist$gene_name
genelist <- sort(genelist,decreasing = T)

genelist <- res$deres$log2FoldChange
names(genelist) <- res$deres$SYMBOL
genelist <- sort(genelist,decreasing = T)
gsea.re1<- GSEA(genelist,  
                TERM2GENE = list_ge,  
                pvalueCutoff = 1,  
                pAdjustMethod = 'BH',seed = 123456) 

gsea.re2<- GSEA(genelist,  
                TERM2GENE = list_ke2,  
                pvalueCutoff = 1,  
                pAdjustMethod = 'BH',seed = 123456) 

gseaNb(object = gsea.re1,subPlot = 2,base_size = 20,lineSize = 1,
       newHtCol = c('#E71F19','white','#5B96D0'),newCurveCol = c('#5B96D0','white','#E71F19'),
       geneSetID = c('GOBP_DETOXIFICATION_OF_COPPER_ION'),addPval = T,
       pvalX = 0.9,pvalY = 0.9,
       pCol = 'black',newGsea = T)