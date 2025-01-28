setwd('~/Downloads/code/')
library(readxl)
library(tidyverse)

count <- read.csv('count.csv',row.names = 1)
count <- count[which(count$gene_biotype=='protein_coding'),]
data <- count
colnames(data) <- c('NC1','NC2','NC3','NC4','NC5','IH1','IH2','IH3','IH4','IH5')

group_list <- c(rep('NC',5),rep('IH',5))
group <- factor(group_list, levels = c('NC', 'IH'), ordered = F)
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
write.csv(DEG,file = 'DEG.csv',row.names = T)

library(clusterProfiler)
library(org.Mm.eg.db)
DEG <- read.csv('DEG.csv',row.names = 1)
dUp <- DEG[which(DEG$log2FoldChange>0&DEG$padj<0.05),]
convert <- bitr(geneID = rownames(dUp),
                fromType = 'ENSEMBL',
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
lipid <- bbb[str_detect(bbb$Description,'lipid'),]
oxida <- bbb[str_detect(bbb$Description,'oxida'),]
death <- bbb[str_detect(bbb$Description,'death'),]
iron <- bbb[str_detect(bbb$Description,'iron|Ferr'),]
lyso <- bbb[str_detect(bbb$Description,'lyso'),]

save(lipid,oxida,death,iron,lyso,file = 'ncih_Down.Rdata')

ekegg <- enrichKEGG(gene = convert$ENTREZID,
                    organism = 'mmu',
                    keyType = 'kegg',
                    pvalueCutoff = 1,
                    minGSSize = 5,
                    maxGSSize = 500)
kkk <- ekegg@result
kkk <- kkk[which(kkk$p.adjust<0.05),]
lipid <- kkk[str_detect(kkk$Description,'(l|L)ipid|(f|F)atty'),]
oxida <- kkk[str_detect(kkk$Description,'[o,O]xi'),]
death <- kkk[str_detect(kkk$subcategory,'death|osis'),]
iron <- kkk[str_detect(kkk$Description,'iron|fer'),]
lyso <- kkk[str_detect(kkk$Description,'TOR'),]
kegg <- rbind(lipid,oxida,death,iron,lyso)
kegg <- kegg[!is.na(kegg$Description),]
enrichmentScore=apply(kegg,1,function(x){
  gr=eval(parse(text = x['GeneRatio']))
  br=eval(parse(text = x['BgRatio']))
  efs=round(gr/br,2)
  efs
})
kegg$FoldEnrichment <- enrichmentScore
save(kegg,file = 'ncih-kegg.Rdata')

b <- DEG[which(DEG$padj<0.05),]
convert <- bitr(geneID = rownames(b),
                fromType = 'ENSEMBL',
                toType = 'ENTREZID',
                OrgDb = org.Mm.eg.db)
convert$logFC <- b$log2FoldChange[match(convert$ENSEMBL,rownames(b))]
gene_list <- convert$logFC
names(gene_list) <- convert$ENTREZID
gene_list <- sort(gene_list, decreasing = T)
Gmt <- read.gmt('~/Desktop/R_Project/msigdb/msigdb.v2024.1.Mm.entrez.gmt')
goGsea <- GSEA(
  gene_list,
  TERM2GENE = Gmt,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05
)
ggg <- goGsea@result
posi <- ggg[ggg$enrichmentScore>0,]
library(GseaVis)
gseaNb(goGsea, geneSetID = rownames(ggg)[66])


library(viridis)
lipid <- lipid[str_detect(lipid$Description,'biosyn'),]
ihnc <- rbind(lipid,death[c(2,3),],bbb[c(830),])
bbb <- dplyr::arrange(bbb, Count, desc(FoldEnrichment))
ihnc <- rbind(bbb[c(76,159,123:125,113:114),])

ihnc <- read.csv('GO_ihnc_up.csv',row.names = 1)
pdf(file = 'GO_ihnc_up.pdf',height = 7,width = 8)
ggplot(ihnc, aes(x=FoldEnrichment, y=Description)) +
  geom_point(aes(y=reorder(Description,FoldEnrichment),color=pvalue, size=Count))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=0,hjust = 1,vjust=0.5))+
  scale_color_viridis()+
  labs(x='FoldEnrichment',y='Description')+guides(size=guide_legend(order=3))
dev.off()

load('ncih-kegg.Rdata')
kegg$Description <- str_split(kegg$Description,' - ',simplify = T)[,1]
pdf(file = 'KEGG_ncih.pdf')
ggplot(data=kegg, aes(x=reorder(Description,FoldEnrichment), y=FoldEnrichment)) + 
  geom_bar(aes(fill=pvalue),stat = 'identity') + 
  scale_fill_viridis()+
  coord_flip() + 
  theme_bw() +  ylab("FoldEnrichment") + xlab("Description") 
dev.off()


final <- dataShow[rownames(DEG)[which(DEG$padj<0.05)],]
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

pdf(file = 'heat.pdf',height = 18)
p
dev.off()
 


