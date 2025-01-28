setwd('~/Downloads/code/')
library(tidyverse)

da <- read.csv('countsBV2.csv')
da <- da[which(da$coding_type=='protein_coding'),]
data <- da[,c(8:10,2:7)]
rownames(data) <- da$id
colnames(data) <- c('NC1','NC2','NC3','IH1','IH2','IH3','IHKO1','IHKO2','IHKO3')

group_list <- c(rep('NC',3),rep('IH',3),rep('IHKO',3))
group <- factor(group_list, levels = c('NC','IH','IHKO'),ordered = T)
annot <- data.frame(row.names = colnames(data),
                    group = group)

library(FactoMineR)
library(factoextra);library(ggpubr)
dat.pca <- PCA(t(data[,c(1:6)]), graph = F)
group_list <- c(rep('NC',3),rep('IH',3))
group <- factor(group_list, levels = c('NC','IH'),ordered = T)
pca <- fviz_pca_ind(dat.pca,
                    geom.ind = c('point', 'text'),
                    pointsize = 4.5,
                    labelsize = 4,
                    label = "ind",
                    repel = T,
                    #invisible = 'ind',
                    col.ind = group,
                    legend.title="group", 
                    palette = c('#E64B35FF','#4DBBD5FF'),
                    addEllipses = F,
                    ellipse.type="confidence",
                    ellipse.level=0.9,
                    axes.linetype = NA,
                    mean.point = F)+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth =1, linetype="solid"))
pca
pdf(file = '2group.pdf',width = 6,height = 6)
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
res <- results(dds2, contrast = c('group', 'IH', 'NC'))
DEG <- as.data.frame(res)
write.csv(DEG,file = 'DEGbv2IHvsNC.csv',row.names = T)

library(viridis)
DEG <- read.csv('DEGbv2IHvsNC.csv',row.names = 1)
final <- data[rownames(DEG)[which(DEG$padj<0.05 & abs(DEG$log2FoldChange)>3)],c(1:6)]
x <- t(scale(t(as.matrix(final))))
ann_cols <- list(group = c(NC = '#E64B35FF', IH = '#4DBBD5FF'))
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

pdf(file = 'heatBV22.pdf',height = 9)
p
dev.off()


library(clusterProfiler)
library(org.Mm.eg.db)
DEG <- read.csv('DEGbv2IHvsNC.csv',row.names = 1)
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

for (i in c('bbb')) {
  enrichmentScore=apply(get(i),1,function(x){
    gr=eval(parse(text = x['GeneRatio']))
    br=eval(parse(text = x['BgRatio']))
    efs=round(gr/br,2)
    efs
  })
  assign(paste0(i,'EF'),enrichmentScore)
}
for (j in c('bbb')) {
  a=get(j)
  a$FoldEnrichment <- get(paste0(j,'EF'))
  assign(j,a)
}

bbb=dplyr::arrange(bbb, Count, FoldEnrichment)
ihnc=bbb[c(646:648,588:589,715:717),]

inflam <- bbb[str_detect(bbb$Description,'inflam|C1|complement|IL|(i|I)nterleukin|TNF|tumor necrosis factor'),]

ekegg <- enrichKEGG(gene = convert$ENTREZID,
                    organism = 'mmu',
                    keyType = 'kegg',
                    pvalueCutoff = 1,
                    minGSSize = 5,
                    maxGSSize = 500)
kkk <- ekegg@result
kkk <- kkk[which(kkk$pvalue<0.05),]
inflamk <- kkk[str_detect(kkk$Description,'(i|I)nflam|C1|(c|C)omplement'),]

 
lll <- letters[1:3]
for (i in 1:3) {
  filler=str_split(inflam$geneID[i],'/',simplify = T)
  assign(lll[i],t(filler)[,1])
}
gene <- unique(c(a,b,c))


ihnc <- inflam[c(2,11,70,20,33,42,58,31),]

library(viridis)
ihnc <- read.csv('bv2_ihncUp.csv',row.names = 1)
pdf(file = 'GO_ihncUp_bv2.pdf',height = 6,width = 8)
ggplot(ihnc, aes(x=FoldEnrichment, y=Description)) +
  geom_point(aes(y=reorder(Description,FoldEnrichment),color=pvalue, size=Count))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=0,hjust = 1,vjust=0.5))+
  scale_color_viridis()+
  labs(x='FoldEnrichment',y='Description')+guides(size=guide_legend(order=3))
dev.off()
