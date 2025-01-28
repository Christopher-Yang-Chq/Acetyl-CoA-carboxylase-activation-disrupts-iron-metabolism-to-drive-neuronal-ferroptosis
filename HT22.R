setwd('~/Downloads/code/')
library(tidyverse)

da <- read.csv('countsHT22.csv')
da <- da[which(da$coding_type=='protein_coding'),]
data <- da[,c(8:10,2:7)]
rownames(data) <- da$id
colnames(data) <- c('NC1','NC2','NC3','IH1','IH2','IH3','IHKO1','IHKO2','IHKO3')

group_list <- c(rep('NC',3),rep('IH',3),rep('IHKO',3))
group <- factor(group_list, levels = c('NC','IH','IHKO'), ordered = F)
annot <- data.frame(row.names = colnames(data),
                    group = group)

library(FactoMineR)
library(factoextra);library(ggpubr)
dat.pca <- PCA(t(data), graph = F)
pca <- fviz_pca_ind(dat.pca,
                    geom.ind = c('point', 'text'),
                    pointsize = 4.5,
                    labelsize = 4,
                    label = "ind",
                    repel = T,
                    #invisible = 'ind',
                    col.ind = group,
                    legend.title="group", 
                    palette = c('#E64B35FF','#4DBBD5FF','#7E6148FF'),
                    addEllipses = F,
                    ellipse.type="confidence",
                    ellipse.level=0.9,
                    axes.linetype = NA,
                    mean.point = F)+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth =1, linetype="solid"))
pca
pdf(file = '3group.pdf',width = 6,height = 6)
pca
dev.off()

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = annot,
                              design = ~ group)
rld <- rlog(dds, blind = F)
datt <- as.data.frame(assay(rld))
#dds$group<- relevel(dds$group, ref="NC")
keep <- rowSums(counts(dds)) >= 1.5*ncol(datt)
dds_filt <- dds[keep, ]
dds2 <- DESeq(dds_filt, quiet = F)
res <- results(dds2, contrast = c('group', 'IHKO', 'IH'))
DEG <- as.data.frame(res)
write.csv(DEG,file = 'DEGht22IHKOvsIH.csv',row.names = T)

library(viridis)
degihnc <- read.csv('DEGht22IHvsNC.csv',row.names = 1)
degihkoih <- read.csv('DEGht22IHKOvsIH.csv',row.names = 1)
degihnc <- degihnc[which(degihnc$padj<0.05 & abs(degihnc$log2FoldChange)>1),]
degihkoih <- degihkoih[which(degihkoih$padj<0.05 & abs(degihkoih$log2FoldChange)>1),]
gene <- intersect(rownames(degihkoih),rownames(degihnc))
final <- data[gene,]
# final <- data[rownames(DEG)[which(DEG$padj<0.05 & abs(DEG$log2FoldChange)>3)],]
x <- t(scale(t(as.matrix(final))))
ann_cols <- list(group = c(NC = '#E64B35FF', IH = '#4DBBD5FF', IHKO = '#7E6148FF'))
p <- pheatmap::pheatmap(na.omit(x), show_colnames = T, show_rownames = F,
                        annotation_colors = ann_cols,
                        fontsize = 10,
                        fontsize_row = 7,
                        fontsize_col = 7,
                        legend_breaks = -3:3,
                        angle_col = 45,
                        cluster_cols = F,
                        cluster_rows = T,
                        annotation_col = annot,
                        cellwidth = 30,
                        cellheight = 2.5,
                        #cluster_rows = F,
                        scale = 'row',
                        color = viridis(8))

a=x[p$tree_row[['order']],]
pheatmap::pheatmap(a, show_colnames = T, show_rownames = F,
                   annotation_colors = ann_cols,
                   fontsize = 10,
                   fontsize_row = 7,
                   fontsize_col = 7,
                   legend_breaks = -3:3,
                   angle_col = 45,
                   cluster_cols = F,
                   cluster_rows = F,
                   annotation_col = annot,
                   cellwidth = 30,
                   cellheight = 2.5,
                   #cluster_rows = F,
                   #scale = 'row',
                   color = viridis(8))



pdf(file = 'heatHT22.pdf',height = 9)
p
dev.off()


library(clusterProfiler)
library(org.Mm.eg.db)
DEG <- read.csv('DEGht22IHKOvsIH.csv',row.names = 1)
dUp <- DEG[which(DEG$log2FoldChange>0 & DEG$padj<0.05),]
convert <- bitr(geneID = rownames(dUp),
                fromType = 'SYMBOL',
                toType = 'ENTREZID',
                OrgDb = org.Mm.eg.db)
ego <- enrichGO(gene = convert$ENTREZID,
                OrgDb = org.Mm.eg.db,
                keyType = 'ENTREZID',
                pvalueCutoff =1,
                ont = 'ALL',
                qvalueCutoff =0.2,
                minGSSize = 5,
                maxGSSize = 5000,
                readable = T)
bbb=ego@result
bbb <- bbb[which(bbb$p.adjust<0.05),]
lipid <- bbb[str_detect(bbb$Description,'lipid'),]
oxida <- bbb[str_detect(bbb$Description,'oxida'),]
death <- bbb[str_detect(bbb$Description,'death'),]
iron <- bbb[str_detect(bbb$Description,'iron|ferr'),]
lyso <- bbb[str_detect(bbb$Description,'lyso|TOR'),]

for (i in c('lipid','oxida','death','iron','lyso','bbb')) {
  enrichmentScore=apply(get(i),1,function(x){
    gr=eval(parse(text = x['GeneRatio']))
    br=eval(parse(text = x['BgRatio']))
    efs=round(gr/br,2)
    efs
  })
  assign(paste0(i,'EF'),enrichmentScore)
}
for (j in c('lipid','oxida','death','iron','lyso','bbb')) {
  a=get(j)
  a$FoldEnrichment <- get(paste0(j,'EF'))
  assign(j,a)
}

lipid <- lipid[str_detect(lipid$Description,'biosyn'),]
bbb <- bbb[order(bbb$FoldEnrichment,decreasing = T),]
ihnc <- rbind(lyso[c(6,7),],bbb[c(188:191,213,214),])
ihnc <- rbind(lipid[c(1,4),],oxida[c(2,3),],death[c(3,6),],iron[c(6,7),])

  
lll <- letters[1:8]
for (i in 1:8) {
  filler=str_split(ihkoih$geneID[i],'/',simplify = T)
  assign(lll[i],t(filler)[,1])
}
gene <- unique(c(a,b,c,d,e,f,g,h))
intre <- data.frame(g=gene)
write.csv(intre,file = 'ihkoihDownHt22.csv')

ihkoih <- rbind(lipid[c(8,9),],oxida[c(1,4),],death[c(5,8),],iron[c(1,2),])
lipid <- lipid[str_detect(lipid$Description,'metab'),]
bbb <- bbb[order(bbb$FoldEnrichment,decreasing = T),]
ihkoih <- rbind(lyso,lipid[c(2,3),],bbb[c(640,642:643,647),])

library(viridis)
ihnc <- read.csv('ht22_ihnc_down.csv',row.names = 1)
pdf(file = 'GO_ihncDown_ht22.pdf',height = 6,width = 8)
ggplot(ihnc, aes(x=FoldEnrichment, y=Description)) +
  geom_point(aes(y=reorder(Description,FoldEnrichment),color=pvalue, size=Count))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=0,hjust = 1,vjust=0.5))+
  scale_color_viridis()+
  labs(x='FoldEnrichment',y='Description')+guides(size=guide_legend(order=3))
dev.off()

ihkoih <- read.csv('ht22_ihkoih_up.csv',row.names = 1)
pdf(file = 'GO_ihkoihUp_ht22.pdf',height = 6,width = 8)
ggplot(ihkoih, aes(x=FoldEnrichment, y=Description)) +
  geom_point(aes(y=reorder(Description,FoldEnrichment),color=pvalue, size=Count))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=0,hjust = 1,vjust=0.5))+
  scale_color_viridis()+
  labs(x='FoldEnrichment',y='Description')+guides(size=guide_legend(order=3))
dev.off()




