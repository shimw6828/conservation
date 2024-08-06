res <- read.csv("/home/shimw/project/conservation/human_enhancer_DE.csv", row.names = 1)
res <- tibble::rownames_to_column(res, var = "gene")
df_data = read.csv("/home/panxl/CRC2/all_cluater")
cluster_type = read.csv("/home/panxl/CRC2/all_cluater",check.names = F)%>%
  dplyr::select("Gene.stable.ID","Human_logFC1",'cluster')
cluster_type$cluster = stringr::str_sub(df_data$cluster,0,9)
enhancer_logf <- dplyr::left_join(res,cluster_type,by = c("gene" = "Gene.stable.ID"))
enhancer_logf <- tidyr::drop_na(enhancer_logf)
nrow(enhancer_logf)
nrow(enhancer_logf[enhancer_logf$padj<0.05,])
dim(PeaktoGene)
head(enhancer_logf)

cor.test(tmp[tmp$cluster%in%c("Cluster 1"),]$log2FoldChange,
         tmp[tmp$cluster%in%c("Cluster 1"),]$Human_logFC1)

##点的分布很不好看
ggplot(data = enhancer_logf, aes(x = log2FoldChange, y = Human_logFC1, colour=cluster))+
  geom_point()+facet_wrap(~ cluster, ncol=3)





#############################
###上面的分析表明，即使合并了信号强度，相关性也没能提升多少，不过一致性14比23高





#################################################
###画log2foldchange的点图
head(enhancer_logf)
table(enhancer_logf$cluster)
enhancer_logf$Human_logFC1

plot_logf <- function(cluster_select){
  tmp_peak <- enhancer_logf[enhancer_logf$cluster==cluster_select,]
  p1 = ggplot(data = tmp_peak, aes(x = log2FoldChange, y = -log10(padj), colour=Human_logFC1)) +
    geom_point(aes(size = abs(Human_logFC1)),alpha=0.5) +
    geom_vline(xintercept=c(0),lty=4,col="grey",lwd=0.7) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.7) +
    scale_size(range = c(0.1,3))+
    theme_classic()+scale_color_gradient2(midpoint=0, low="blue", mid="white",
                                          high="red", space ="Lab", limits = c(0, 4),na.value = "red" )
  return(p1)
}

plot_logf("Cluster 1")


enhancer_logf$logf_color = case_when(enhancer_logf$Human_logFC1>3 ~ 3,
                                     enhancer_logf$Human_logFC1< -3 ~ -3,
                                          TRUE ~ enhancer_logf$Human_logFC1)
enhancer_logf$logf_size = case_when(enhancer_logf$Human_logFC1>4 ~ 4,
                                    enhancer_logf$Human_logFC1< -4 ~ -4,
                                         TRUE ~ enhancer_logf$Human_logFC1)
p = ggplot(data = enhancer_logf, aes(x = log2FoldChange, y = -log10(padj), colour=logf_color)) +
  geom_point(aes(size = abs(logf_size)),alpha=0.5,shape = 19,stroke = 0) +
  geom_vline(xintercept=c(0),lty=4,col="grey",lwd=0.7) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.7) +
  scale_size(range = c(0,4))+
  labs( x='-log10 (P_value)',y = "Enhancer intensity log2 (FoldChange)", color = "",size = "") +
  theme_classic(base_family = "sans",base_size = 16,base_line_size = 0.9)+
  scale_color_gradient2(midpoint=0, low="blue", mid="white",
                        high="red", space ="Lab", limits = c(-3, 3),breaks = c(-2,0,2))+
  facet_wrap(~ cluster, ncol=3)+
  theme(legend.title = element_text(size=14),legend.position = c(0.9, 0),
        legend.justification = c(1, 0),legend.title.align = 0.5,
        axis.title=element_text(size=18),
        axis.text=element_text(size=16),
        strip.background = element_blank(),legend.box = "horizontal",
        strip.text.x = element_text(size = 18)) +
  guides(    color = guide_colorbar(order = 1),
             fill = guide_legend(order = 0))
library(cowplot)


Cairo::CairoPDF("/home/shimw/project/conservation/enhancer_logf2_with_expression.pdf",width = 10,height = 5.6)
ggdraw(p) + draw_label("Gene expression log2 (FoldChange)", x = 0.83, y = 0.40)##取消了标签，自己手动标了个lable当图例
dev.off()


cor.test(cluster_PeaktoGene[cluster_PeaktoGene$cluster%in%c("Cluster 4"),]$log2FoldChange,
         cluster_PeaktoGene[cluster_PeaktoGene$cluster%in%c("Cluster 4"),]$Human_logFC1)

###############
##数字多，更好看些
res <- read.csv("/home/shimw/project/conservation/human_enhancer_DE.csv",row.names = 1)
res$chr = stringr::str_split_fixed(row.names(res),"_",3)[,1]
res$start = stringr::str_split_fixed(row.names(res),"_",3)[,2]
res$end = stringr::str_split_fixed(row.names(res),"_",3)[,3]
res = res[,c('chr','start','end','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj')]
readr::write_tsv(res,"/home/shimw/project/conservation/human_enhancer_DE.bed",col_names =F)
library(ChIPseeker)
library(GenomicFeatures)
hg19 <-  makeTxDbFromGFF("/NAS/luozh/CRC_conservation/reference/gencode.v32lift37.annotation.gtf", format = "gtf")
peak <- readPeakFile("/home/shimw/project/conservation/human_enhancer_DE.bed")
peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000), TxDb=hg19)
# plotAnnoPie(peakAnno)
peakAnno <- as.GRanges(peakAnno)
PeaktoGene <- data.frame("gene" = gsub("\\..*", "", peakAnno$geneId),"baseMean"=peakAnno$V4,
                         'log2FoldChange'=peakAnno$V5,'lfcSE'=peakAnno$V6,
                         'stat'=peakAnno$V7,'pvalue' = peakAnno$V8,'padj' = peakAnno$V9, 
                         "h_peak_id" = stringr::str_c(as.vector(seqnames(peakAnno)),start(peakAnno),end(peakAnno),sep = "_"))
dim(PeaktoGene[PeaktoGene$padj<0.05,])
# PeaktoGene = PeaktoGene[PeaktoGene$padj<0.05,]
# PeaktoGene = PeaktoGene[PeaktoGene$padj<0.05&abs(PeaktoGene$log2FoldChange)>1,]
dim(PeaktoGene)
df_data = read.csv("/home/panxl/CRC2/all_cluater")
df_data$cluster = stringr::str_sub(df_data$cluster,0,9)
df_data = df_data[,c("Gene.stable.ID","Human_logFC1","cluster")]
names(df_data)[1] = "gene"
table(df_data$gene%in%PeaktoGene$gene)
dim(PeaktoGene)
table(PeaktoGene$gene%in%df_data$gene)
cluster_PeaktoGene = dplyr::left_join(df_data,PeaktoGene)%>%tidyr::drop_na()
cluster_PeaktoGene <- cluster_PeaktoGene[,c("h_peak_id","gene","log2FoldChange","padj","Human_logFC1","cluster")]
names(cluster_PeaktoGene) = c("h_peak_id","h_gene","h_log2FoldChange","h_padj","Human_logFC1","cluster")
# readr::write_csv(cluster_PeaktoGene,"/home/shimw/project/conservation/enhancer_count/human_enhancer_DE_annotated.csv")
# 计算靠近上调enhancer的基因数
up_enhancer = cluster_PeaktoGene$h_padj<0.05&cluster_PeaktoGene$h_log2FoldChange>0
nrow(df_data[df_data$cluster=="Cluster 1",])
length(unique(cluster_PeaktoGene[cluster_PeaktoGene$cluster=="Cluster 1"&up_enhancer,]$h_gene))
495/1529
nrow(df_data[df_data$cluster=="Cluster 4",])
length(unique(cluster_PeaktoGene[cluster_PeaktoGene$cluster=="Cluster 4"&up_enhancer,]$h_gene))
214/560
# 计算靠近下调enhancer的基因数
down_enhancer = cluster_PeaktoGene$h_padj<0.05&cluster_PeaktoGene$h_log2FoldChange<0
nrow(df_data[df_data$cluster=="Cluster 2",])
length(unique(cluster_PeaktoGene[cluster_PeaktoGene$cluster=="Cluster 2"&down_enhancer,]$h_gene))
256/1092
nrow(df_data[df_data$cluster=="Cluster 3",])
length(unique(cluster_PeaktoGene[cluster_PeaktoGene$cluster=="Cluster 3"&down_enhancer,]$h_gene))
613/2049



plot_logf <- function(cluster_select){
  tmp_peak <- cluster_PeaktoGene[cluster_PeaktoGene$cluster==cluster_select,]
  p1 = ggplot(data = tmp_peak, aes(x = log2FoldChange, y = -log10(padj), colour=Human_logFC1)) +
    geom_point(aes(size = abs(Human_logFC1)),alpha=0.5) +
    geom_vline(xintercept=c(0),lty=4,col="grey",lwd=0.7) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.7) +
    scale_size(range = c(0.1,3))+
    theme_classic()+scale_color_gradient2(midpoint=0, low="blue", mid="white",
                                          high="red", space ="Lab", limits = c(0, 3),na.value = "red" )
  return(p1)
}
###为了画图一致，header变回去
names(cluster_PeaktoGene) <- c("h_peak_id", "gene","log2FoldChange","padj","Human_logFC1","cluster")
cluster_PeaktoGene$logf_color = case_when(cluster_PeaktoGene$Human_logFC1>3 ~ 3,
                                          cluster_PeaktoGene$Human_logFC1< -3 ~ -3,
                                          TRUE ~ cluster_PeaktoGene$Human_logFC1)
cluster_PeaktoGene$logf_size = case_when(cluster_PeaktoGene$Human_logFC1>4 ~ 4,
                                          cluster_PeaktoGene$Human_logFC1< -4 ~ -4,
                                          TRUE ~ cluster_PeaktoGene$Human_logFC1)
p1 = ggplot(data = cluster_PeaktoGene, aes(x = log2FoldChange, y = -log10(padj), colour=logf_color)) +
  geom_point(aes(size = abs(logf_size)),alpha=0.5,shape = 19,stroke = 0) +
  geom_vline(xintercept=c(0),lty=4,col="grey",lwd=0.7) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.7) +
  scale_size(range = c(0,4))+
  labs( y='-log10 (P-value)',x = "H3K27ac log2 (FoldChange)", color = "",size = "") +
  theme_classic(base_family = "sans",base_size = 16,base_line_size = 0.9)+
  scale_color_gradient2(midpoint=0, low="blue", mid="white",
                        high="red", space ="Lab", limits = c(-3, 3),breaks = c(-2,0,2))+
  facet_wrap(~ cluster, ncol=3)+
  theme(legend.title = element_text(size=14),legend.position = c(0.9, 0),
        legend.justification = c(1, 0),legend.title.align = 0.5,
        axis.title=element_text(size=18),
        axis.text=element_text(size=16),
        strip.background = element_blank(),legend.box = "horizontal",
        strip.text.x = element_text(size = 18)) +
  guides(    color = guide_colorbar(order = 1),
             fill = guide_legend(order = 0))

library(cowplot)

Cairo::CairoPDF("/home/shimw/project/conservation/human_enhancer_logf2_with_expression.pdf",width = 10,height = 5.6)
ggdraw(p1) + draw_label("Gene expression log2 (FoldChange)", x = 0.83, y = 0.40)
dev.off()



plot_logf("Cluster 1")

cluster_select="Cluster 1"

############################################
##判断同一个基因的peak是否都大于0或者小于0
mm = purrr::map(unique(cluster_PeaktoGene$gene),function(gg){
  nnn = cluster_PeaktoGene[cluster_PeaktoGene$gene==gg,]
  return(sum(nnn$log2FoldChange>0)>length(nnn$log2FoldChange)*0.6|sum(nnn$log2FoldChange<0)>length(nnn$log2FoldChange)*0.7)
})%>%unlist()
table(mm)


##############################################################
####小鼠
###使用了罗师兄直接提供的文件计算。

# mouse_peak <- readr::read_tsv("/home/zhihl/Project/CRC/Chip_analysis/dff/version0821/final_result/all_diff_data_H3K27ac_add_symbol.txt")
# mouse_peak <- mouse_peak[,c("peak_id","ensemblID","log2FoldChange.tumorVSctrl", "padj.tumorVSctrl")]
# 
# cluster_type = read.csv("/home/panxl/CRC2/all_cluater",check.names = F)%>%
#   dplyr::select("Mouse.gene.stable.ID","Mouse_logFC1",'cluster')%>%
#   dplyr::mutate(cluster = stringr::str_sub(cluster,0,9))
# 
# mouse_peak_to_gene <- left_join(mouse_peak,cluster_type,by = c("ensemblID"="Mouse.gene.stable.ID"))%>%
#   tidyr::drop_na()
# head(mouse_peak_to_gene)
# table(mouse_peak_to_gene$cluster)
# # cor.test(mouse_peak_to_gene[mouse_peak_to_gene$cluster%in%c("Cluster 2"),]$log2FoldChange.tumorVSctrl,
# #          mouse_peak_to_gene[mouse_peak_to_gene$cluster%in%c("Cluster 2"),]$Mouse_logFC1)
# 
# ggplot(data = mouse_peak_to_gene, aes(x = log2FoldChange.tumorVSctrl, y = Mouse_logFC1, colour=cluster))+
#   geom_point()+facet_wrap(~ cluster, ncol=3)


####保持一致
# res <- read.csv("/home/shimw/project/conservation/mouse_enhancer_DE.csv",row.names = 1)
# res$chr = stringr::str_split_fixed(row.names(res),"_",3)[,1]
# res$start = stringr::str_split_fixed(row.names(res),"_",3)[,2]
# res$end = stringr::str_split_fixed(row.names(res),"_",3)[,3]
# res = res[,c('chr','start','end','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj')]
# readr::write_tsv(res,"/home/shimw/project/conservation/mouse_enhancer_DE.bed",col_names =F)
mm10 <-  makeTxDbFromGFF("/home/shimw/project/conservation/gencode.vM20.basic.annotation.gtf", format = "gtf")
peak <- readPeakFile("/home/shimw/project/conservation/mouse_enhancer_DE.bed")
peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000), TxDb=mm10)
# plotAnnoPie(peakAnno)
peakAnno <- as.GRanges(peakAnno)
PeaktoGene <- data.frame("gene" = gsub("\\..*", "", peakAnno$geneId),"baseMean"=peakAnno$V4,
                         'log2FoldChange'=peakAnno$V5,'lfcSE'=peakAnno$V6,
                         'stat'=peakAnno$V7,'pvalue' = peakAnno$V8,'padj' = peakAnno$V9)
dim(PeaktoGene[PeaktoGene$padj<0.05,])
# PeaktoGene = PeaktoGene[PeaktoGene$padj<0.05,]
# PeaktoGene = PeaktoGene[PeaktoGene$padj<0.05&abs(PeaktoGene$log2FoldChange)>1,]
dim(PeaktoGene)
df_data = read.csv("/home/panxl/CRC2/all_cluater")
df_data$cluster = stringr::str_sub(df_data$cluster,0,9)
df_data = df_data[,c("Gene.stable.ID","Human_logFC1",'cluster')]
names(df_data)[1] = "gene"
table(df_data$gene%in%PeaktoGene$gene)
dim(PeaktoGene)
table(PeaktoGene$gene%in%df_data$gene)
cluster_PeaktoGene = dplyr::left_join(df_data,PeaktoGene)%>%tidyr::drop_na()
table(cluster_PeaktoGene[cluster_PeaktoGene$log2FoldChange>0,]$cluster)

up_enhancer = cluster_PeaktoGene$h_padj<0.05&cluster_PeaktoGene$h_log2FoldChange>0
nrow(df_data[df_data$cluster=="Cluster 1",])
length(unique(cluster_PeaktoGene[cluster_PeaktoGene$cluster=="Cluster 1"&up_enhancer,]$h_gene))
495/1529
nrow(df_data[df_data$cluster=="Cluster 4",])
length(unique(cluster_PeaktoGene[cluster_PeaktoGene$cluster=="Cluster 4"&up_enhancer,]$h_gene))
214/560
# 计算靠近下调enhancer的基因数
down_enhancer = cluster_PeaktoGene$h_padj<0.05&cluster_PeaktoGene$h_log2FoldChange<0
nrow(df_data[df_data$cluster=="Cluster 2",])
length(unique(cluster_PeaktoGene[cluster_PeaktoGene$cluster=="Cluster 2"&down_enhancer,]$h_gene))
256/1092
nrow(df_data[df_data$cluster=="Cluster 3",])
length(unique(cluster_PeaktoGene[cluster_PeaktoGene$cluster=="Cluster 3"&down_enhancer,]$h_gene))
613/2049




# nrow(mouse_peak_to_gene[mouse_peak_to_gene$padj.tumorVSctrl<0.05,])
nrow(mouse_peak_to_gene)
mouse_peak_to_gene$logf_color = case_when(mouse_peak_to_gene$Mouse_logFC1>3 ~ 3,
                                          mouse_peak_to_gene$Mouse_logFC1< -3 ~ -3,
                                          TRUE ~ mouse_peak_to_gene$Mouse_logFC1)
mouse_peak_to_gene$logf_size = case_when(mouse_peak_to_gene$Mouse_logFC1>4 ~ 4,
                                         mouse_peak_to_gene$Mouse_logFC1< -4 ~ -4,
                                         TRUE ~ mouse_peak_to_gene$Mouse_logFC1)
p1 = ggplot(data = mouse_peak_to_gene, aes(x = log2FoldChange.tumorVSctrl, y = -log10(padj.tumorVSctrl), colour=logf_color)) +
  geom_point(aes(size = abs(logf_size)),alpha=0.5,shape = 19,stroke = 0) +
  geom_vline(xintercept=c(0),lty=4,col="grey",lwd=0.7) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.7) +
  scale_size(range = c(0,4))+
  labs( y='-log10 (P-value)',x = "H3K27ac log2 (FoldChange)", color = "",size = "") +
  theme_classic(base_family = "sans",base_size = 16,base_line_size = 0.9)+
  scale_color_gradient2(midpoint=0, low="blue", mid="white",
                        high="red", space ="Lab", limits = c(-3, 3),breaks = c(-2,0,2))+
  facet_wrap(~ cluster, ncol=3)+
  theme(legend.title = element_text(size=14),legend.position = c(0.9, 0),
        legend.justification = c(1, 0),legend.title.align = 0.5,
        axis.title=element_text(size=18),
        axis.text=element_text(size=16),
        strip.background = element_blank(),legend.box = "horizontal",
        strip.text.x = element_text(size = 18)) +
  guides(    color = guide_colorbar(order = 1),
             fill = guide_legend(order = 0))

Cairo::CairoPDF("/home/shimw/project/conservation/mouse_enhancer_logf2_with_expression.pdf",width = 10,height = 5.6)
ggdraw(p1) + draw_label("Gene expression log2 (FoldChange)", x = 0.83, y = 0.40)
dev.off()



# /home/shimw/project/conservation/enhancer_count/human_enhancer_DE_annotated.csv
# /home/shimw/project/conservation/enhancer_count/mouse_enhancer_DE_annotated.csv
## 人类靠近enhancer的基因数量
cluster_PeaktoGene <- readr::read_csv("/home/shimw/project/conservation/enhancer_count/human_enhancer_DE_annotated.csv")
up_enhancer = cluster_PeaktoGene$h_padj<0.05&cluster_PeaktoGene$h_log2FoldChange>0
nrow(df_data[df_data$cluster=="Cluster 1",])
length(unique(cluster_PeaktoGene[cluster_PeaktoGene$cluster=="Cluster 1"&up_enhancer,]$h_gene))
495/1529
nrow(df_data[df_data$cluster=="Cluster 4",])
length(unique(cluster_PeaktoGene[cluster_PeaktoGene$cluster=="Cluster 4"&up_enhancer,]$h_gene))
214/560
# 计算靠近下调enhancer的基因数
down_enhancer = cluster_PeaktoGene$h_padj<0.05&cluster_PeaktoGene$h_log2FoldChange<0
nrow(df_data[df_data$cluster=="Cluster 2",])
length(unique(cluster_PeaktoGene[cluster_PeaktoGene$cluster=="Cluster 2"&down_enhancer,]$h_gene))
256/1092
nrow(df_data[df_data$cluster=="Cluster 3",])
length(unique(cluster_PeaktoGene[cluster_PeaktoGene$cluster=="Cluster 3"&down_enhancer,]$h_gene))
613/2049

# 小鼠靠近enhancer的基因数量
cluster_PeaktoGene <- readr::read_csv("/home/shimw/project/conservation/enhancer_count/mouse_enhancer_DE_annotated.csv")
up_enhancer = cluster_PeaktoGene$m_padj<0.05&cluster_PeaktoGene$m_log2FoldChange>0
length(unique(df_data[df_data$cluster=="Cluster 1",]$Mouse.gene.stable.ID))
length(unique(cluster_PeaktoGene[cluster_PeaktoGene$cluster=="Cluster 1"&up_enhancer,]$m_gene))
259/1529
length(unique(df_data[df_data$cluster=="Cluster 2",]$Mouse.gene.stable.ID))
length(unique(cluster_PeaktoGene[cluster_PeaktoGene$cluster=="Cluster 2"&up_enhancer,]$m_gene))
299/1092
# 计算靠近下调enhancer的基因数
down_enhancer = cluster_PeaktoGene$m_padj<0.05&cluster_PeaktoGene$m_log2FoldChange<0
length(unique(df_data[df_data$cluster=="Cluster 3",]$Mouse.gene.stable.ID))
length(unique(cluster_PeaktoGene[cluster_PeaktoGene$cluster=="Cluster 3"&down_enhancer,]$m_gene))
918/2049
length(unique(df_data[df_data$cluster=="Cluster 4",]$Mouse.gene.stable.ID))
length(unique(cluster_PeaktoGene[cluster_PeaktoGene$cluster=="Cluster 4"&down_enhancer,]$m_gene))
224/560





