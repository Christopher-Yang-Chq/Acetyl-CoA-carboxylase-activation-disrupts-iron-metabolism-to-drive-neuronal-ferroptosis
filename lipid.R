setwd('~/Downloads/code/')
library(tidyverse)

raw <- read.csv('result/00.Raw_files/PDX-2024-0006.NEG.all_sample.csv',row.names = 1)
subPL <- raw[str_detect(raw$Abbrev, 'PC|PE|PS|PI|PG|PA'),]
subPL$type <- str_extract(subPL$Abbrev, 'PC|PE|PS|PI|PG|PA')

raw$a <- str_split(raw$Description,'/',simplify = T)[,1]
raw$b <- str_split(raw$Description,'/',simplify = T)[,2]

raw$a <- str_extract(raw$a, '(:[0-9]+)')
raw$b <- str_extract(raw$b, '(:[0-9]+)')

raw$a <- str_split(raw$a,':',simplify = T)[,2]
raw$b <- str_split(raw$b,':',simplify = T)[,2]

posi <- which(raw$a=='' & raw$b=='')
raw$a[posi] <- str_split(str_extract(raw$Abbrev[posi],'(:[0-9]+)'),':',simplify = T)[,2]
raw[,c(43,44)] <- lapply(raw[,c(43,44)], function(x)as.numeric(x))
raw <- raw[which(!is.na(raw$a)),]

cordi <- which(raw$b != '')
corddd <- which(is.na(raw$b))
raw$FA <- NA
raw$FA[cordi] <- ifelse(raw$a[cordi] > 0 | raw$b[cordi] > 0,
                        ifelse(raw$a[cordi] >1 | raw$b[cordi] >1, 'PUFA', 'MUFA'), 'SFA')
raw$FA[corddd] <- ifelse(raw$a[corddd] > 0,
                         ifelse(raw$a[corddd] >1, 'PUFA', 'MUFA'), 'SFA') 
# write.csv(raw, 'data/subFA.csv')


dat <- raw[,c(1:18)]
boxplot(dat,outline=FALSE, col=c(rep('green', 6), rep('blue', 6),rep('red',6)), notch=T, las=2)
ddd <- gather(dat, sam, exp, NC1:IH_E6)
ggplot(ddd,aes(x=exp))+
  geom_density(aes(color=sam))

data1 <- dat[,7:18]  


s1 <- colnames(data1)[7:12] 
s2 <- colnames(data1)[1:6] 
pvalue = padj = log2FoldChange = matrix(0, nrow(data1), 1)
for(i in 1:nrow(data1)){
  
  pvalue[i, 1] = p.value = t.test(c(t(data1[i, s1])), c(t(data1[i, s2])), exact = F)$p.value
  log2FoldChange[i, 1] = log2(mean(c(t(data1[i, s1]))) / mean(c(t(data1[i, s2]))))
}

padj = p.adjust(as.vector(pvalue), "fdr", n = length(pvalue))
rTable = data.frame(log2FoldChange, pvalue, padj, row.names = rownames(data1))
treatment_log10 <- signif(apply(data1[rownames(rTable), s1], 1, mean), 4)
control_log10 <- signif(apply(data1[rownames(rTable), s2], 1, mean), 4)

cat("mark DGE\n") 
DGE <- rep("NC", nrow(data1))
DGE[((rTable$pvalue) < 0.05) & (rTable$log2FoldChange > 0)] = "UP"
DGE[((rTable$pvalue) < 0.05) & (rTable$log2FoldChange < 0)] = "DN"
table(DGE)
compound = rownames(data1)
rTable = data.frame(treatment_log10, control_log10, rTable[, c("log2FoldChange", "pvalue", "padj")], DGE)
head(rTable)

rTable$Description <- raw$Description
rTable$Abbrev <- raw$Abbrev
write.csv(rTable,file = 'iheihDce.csv',row.names = T)

ihnc <- read.csv('ihncDce.csv',row.names = 1)
iheih <- read.csv('iheihDce.csv',row.names = 1)
compound <- unique(c(rownames(ihnc)[which(ihnc$DGE != 'NC')],rownames(iheih)[which(iheih$DGE != 'NC')]))
dat <- raw[,c(1:12)]
final <- dat[compound,]
x <- t(scale(t(as.matrix(final))))
group_list <- c(rep('NC',6),rep('IH',6),rep('IH_E',6))
group <- factor(group_list, levels = c('NC','IH','IH_E'), ordered = T)
annot <- data.frame(row.names = colnames(dat),
                    group = group)
ann_cols <- list(group = c(NC = '#E64B35FF', IH = '#4DBBD5FF', IH_E = '#00A087FF'))
library(viridis)
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
                        cellheight = 3.5,
                        #cluster_rows = F,
                        scale = 'row',
                        color = viridis(8),
                        border = F)

pdf(file = '3heatUnion.pdf',width = 12, height = 8)
p
dev.off()

library(ggplot2)
ihnc$type <- str_extract(raw$Abbrev, 'PC|PE|PS|PI|PG|PA')
ihnc$type <- ifelse(ihnc$DGE=='NC',NA,ihnc$type)

iheih$type <- str_extract(raw$Abbrev, 'PC|PE|PS|PI|PG|PA')
iheih$type <- ifelse(iheih$DGE=='NC',NA,iheih$type)

volc_plot <- ggplot(ihnc, aes(log2FoldChange, -log10(pvalue))) +
  geom_point(size = 1.5) + 
  #geom_vline(xintercept = c(-1,1), linetype='dotted')+
  geom_hline(yintercept = 1.30103, linetype='dotted')+
  xlab(expression("log"[2]*" fold change")) +
  ylab(expression("-log"[10]*" p-value")) + 
  #scale_x_continuous(limits = c(-2.5, 2.5)) + 
  #scale_color_manual(values = c("#00468BFF", "#ED0000FF", "#42B540FF","#925E9FFF","#FDAF91FF"))+
  theme_bw()+
  theme(panel.grid=element_blank())
volc_plot

ihnccol <- ihnc[!is.na(ihnc$type),]

iheihcol <- iheih[!is.na(iheih$type),]

volc_plot1 <- volc_plot+
  geom_point(data = ihnccol,aes(color=type),size = 1.5)+
  scale_color_manual(values = c("#00468BFF", "#ED0000FF", "#42B540FF",#"#0099B4FF",
                                "#925E9FFF","#FDAF91FF"))  
volc_plot1

showihnc <- ihnc[which(ihnc$DGE!='NC'),]
showihnc <- showihnc[c(2,10,11,30,31,35,36,38,40),]

showiheih <- iheih[which(iheih$DGE!='NC'),]
showiheih <- showiheih[c(6,12,13,18,19,38,52,53,54,55,58,64,66,69,72,73,77,80,82,90,95,98,99,101,102),]

library(ggrepel)
volc_plot2 <- volc_plot1 +
  geom_text_repel(data = subset(showihnc, log2FoldChange>0),
                   force = 20,color="grey20",size=3,point.padding = 0.5,hjust = 0.5,
                   aes(log2FoldChange, -log10(pvalue), label = Description),
                   #arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
                   segment.color="grey20",segment.size=0.2,segment.alpha=0.8,nudge_y=1,nudge_x = 1.5)+
  geom_text_repel(data = subset(showihnc, log2FoldChange<0),
                  force = 20,color="grey20",size=3,point.padding = 0.5,hjust = 0.5,
                  aes(log2FoldChange, -log10(pvalue), label = Description),
                  #arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
                  segment.color="grey20",segment.size=0.2,segment.alpha=0.8,nudge_y=1,nudge_x = -1.5)
volc_plot2

pdf(file = 'ihncVolcano.pdf',width = 8,height = 6)
volc_plot2
dev.off()


subPL <- read.csv('data/subPL.csv',row.names = 1)

showiheih$FA <- subPL$FA[match(rownames(showiheih),rownames(subPL))]
showihnc$FA <- subPL$FA[match(rownames(showihnc),rownames(subPL))]

showiheih <- showiheih[order(showiheih$DGE, showiheih$FA),]
showihnc <- showihnc[order(showihnc$DGE, showihnc$FA),]

dat <- raw[rownames(showiheih),7:18]
dat <- raw[rownames(showihnc),1:12]

rownames(dat) <- showiheih$Description
rownames(dat) <- showihnc$Description

x <- t(scale(t(as.matrix(dat))))
group_list <- c(rep('IH',6),rep('IH_E',6))
group_list <- c(rep('NC',6),rep('IH',6))

group <- factor(group_list, levels = c('IH','IH_E'), ordered = T)
group <- factor(group_list, levels = c('NC','IH'), ordered = T)

annot <- data.frame(row.names = colnames(dat),
                    group = group)
annotr <- data.frame(row.names = showiheih$Description,
                     FA = showiheih$FA)
annotr <- data.frame(row.names = showihnc$Description,
                     FA = showihnc$FA)

ann_cols <- list(group = c(IH = '#4DBBD5FF', IH_E = '#00A087FF'), 
                 FA = c(SFA = '#EE4C97FF', MUFA = '#FFDC91FF', PUFA = '#6F99ADFF'))
ann_cols <- list(group = c(NC = '#E64B35FF', IH = '#4DBBD5FF'),
                 FA = c(SFA = '#EE4C97FF', MUFA = '#FFDC91FF', PUFA = '#6F99ADFF'))
library(viridis)
p <- pheatmap::pheatmap(x, show_colnames = T, show_rownames = T,
                        annotation_colors = ann_cols,
                        fontsize = 10,
                        fontsize_row = 7,
                        fontsize_col = 7,
                        legend_breaks = -3:3,
                        angle_col = 45,
                        cluster_cols = F,
                        cluster_rows = F,
                        annotation_col = annot,
                        annotation_row = annotr,
                        cellwidth = 30,
                        cellheight = 8,
                        #cluster_rows = F,
                        scale = 'row',
                        color = viridis(8),
                        border = F)

pdf(file = 'ihncHeatSelectFA.pdf',width = 10, height = 7.5)
p
dev.off()


library(ropls);library(ggrepel)
dat <- t(raw[,1:18])
group <- c(rep('NC',6),rep('IH',6),rep('IH_E',6))
plsda <-  opls(dat, group)
sample.score = plsda@scoreMN[,1] %>%  
  as.data.frame() %>%
  mutate(group = group,
         o1 = plsda@scoreMN[,2]) 
head(sample.score)
colnames(sample.score) <- c('p1','group','o1')
sample.score$group <- factor(sample.score$group, levels = c('NC','IH','IH_E'), ordered = T)
p <- ggplot(sample.score, aes(p1, o1, color = group, shape= group)) +
  geom_hline(yintercept = 0, linetype = 'dashed', linewidth = 0.1) + 
  geom_vline(xintercept = 0, linetype = 'dashed', linewidth = 0.1) +
  geom_point(size=4.5) +
  geom_text_repel(aes(label = rownames(sample.score)), size = 4, vjust = 2)+
  #geom_point(aes(-10,-10), color = 'white') +
  labs(x = 'P1 (18%)',y = 'P2 (14%)') +
  stat_ellipse(level = 0.95, linetype = 1, geom = "polygon",
               aes(fill=group), alpha=0.1,
               size = 0.8, show.legend = T) + 
  scale_fill_manual(values = c('#E64B35FF','#4DBBD5FF','#00A087FF'))+
  scale_color_manual(values = c('#E64B35FF','#4DBBD5FF','#00A087FF')) +
  theme_bw() +
  theme(legend.position = c(0.1,0.85),
        legend.title = element_blank(),
        legend.text = element_text(color = 'black',size = 12),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black',size = 15),
        axis.title = element_text(color = 'black',size = 15),
        axis.ticks = element_line(color = 'black'))
p
ggsave(filename = 'plsda.pdf',width = 7, height = 6.5)




