library('GenomicRanges') 
library('rtracklayer') 
library('dplyr') 
library(GenomicFeatures)
library(ChIPseeker)
bedfile = "/home/zhluo/Project/CRC_conservation/step13_chip_seq_DEG/all_narrow_peak_human/human_H3K27ac_master.sort.merge.bed"
data<-read.table(bedfile,header=F)
bigwig_dir <- "/home/zhluo/Project/CRC_conservation/step13_chip_seq_DEG/all_narrow_peak_human/bigwig_enhancer"
list.files(bigwig_dir,full.names = F)
pr_files = list.files(bigwig_dir,full.names = T)
file_name = list.files(bigwig_dir,full.names = F)%>%stringr::str_remove(".bigwig")
names(pr_files)=file_name
outputfname="/home/mwshi/project/conservation/human_enhancer_tss_counts.csv"
colnames(data)<-c("chr","start","end")
bed<-with(data,GRanges(chr,IRanges(start,end)))
hg19 <-  makeTxDbFromGFF("/home/mwshi/project/conservation/gencode.v32lift37.annotation.gtf", format = "gtf")
peakAnno <- annotatePeak(bed, tssRegion=c(-3000, 3000), TxDb=hg19)
peakAnno <- as.GRanges(peakAnno)
peakAnno = peakAnno[startsWith(peakAnno$annotation,"Promoter"),]
peakAnno$geneId = stringr::str_split_fixed(peakAnno$geneId,"\\.",2)[,1]
bw=BigWigFileList(pr_files) 
counts <- matrix(NA, nrow = length(unique(peakAnno$geneId)), ncol = length(bw))
colnames(counts) <- names(bw)
chromnames=levels(seqnames(bed))
for(i in seq_len(length(bw))) 
{
  last=1 	 
  coverage <- import(bw[[i]], as = 'RleList')
  for(j in unique(peakAnno$geneId))
  {
    tmp_peak = peakAnno[peakAnno$geneId==j]
    tmp_chr = as.vector(seqnames(tmp_peak)[1])
    range_vals=ranges(tmp_peak)
    cur_coverage=coverage[[tmp_chr]] 
    newvals=sum(Views(cur_coverage,range_vals))
    counts[last, i] <-sum(newvals) 
    last=last+1
  }
}
counts <- round(counts / 36, 0)
row.names(counts) <- unique(peakAnno$geneId)
write.csv(counts,file=outputfname,quote = F)
counts <- read.csv(outputfname,row.names = 1)
library(DESeq2)
cData <- data.frame(sample = colnames(counts),
                    "group" = factor(ifelse(startsWith(colnames(counts),"N"),"Normal","Tumor"),
                                     levels = c("Tumor","Normal")))
rownames(cData) <- colnames(counts)
d.deseq <- DESeqDataSetFromMatrix(countData = counts, colData = cData,design = ~group)
d.deseq <- DESeq(d.deseq)
res <- results(d.deseq, contrast=c("group","Tumor","Normal"))
write.csv(res,"/home/mwshi/project/conservation/human_enhancer_tss_DE.csv",quote = F)
# res <- read.csv("/home/mwshi/project/conservation/human_enhancer_DE.csv",row.names = 1)
# res$chr = stringr::str_split_fixed(row.names(res),"_",3)[,1]
# res$start = stringr::str_split_fixed(row.names(res),"_",3)[,2]
# res$end = stringr::str_split_fixed(row.names(res),"_",3)[,3]
# res = res[,c('chr','start','end','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj')]
# readr::write_tsv(res,"/home/mwshi/project/conservation/human_enhancer_DE.bed",col_names =F)

####################################################################
###信号合并后进行画图
res <- read.csv("/home/shimw/project/conservation/human_enhancer_DE.csv",row.names = 1)
res <- tibble::rownames_to_column(res,var = "gene")
df_data = read.csv("/home/panxl/CRC2/all_cluater")
df_data$cluster = stringr::str_sub(df_data$cluster,0,9)
df_data = df_data[,c("Gene.stable.ID","cluster")]
names(df_data)[1] = "gene"
table(res$gene%in%df_data$gene)
cluster_enhancer = dplyr::left_join(df_data,res)%>%tidyr::drop_na()
dim(cluster_enhancer)
my_comparisons <- list(c("Cluster 1", "Cluster 2"), c("Cluster 1", "Cluster 3"), c("Cluster 1", "Cluster 4"), c("Cluster 1", "Cluster 5"))



ggplot(data = na.omit(cluster_enhancer),aes(x = cluster,y = log2FoldChange,fill=cluster)) +
  geom_violin() +geom_boxplot(width=0.1, fill="white")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif", tip.length=0.02,method.args = list(alternative = "greater"), method = "t.test") +
  theme_classic(base_family = "sans",base_size = 18,base_line_size = 1.1)+
  scale_fill_manual(values=c("Cluster 1" = "#CB9C7A", "Cluster 2" = "#8696a7", "Cluster 3" = "#CDB97D",
                             "Cluster 4" = "#7b8b6f", "Cluster 5" = "#A59B95"))+
  labs(title = "peak intensity", y="Log2foldchange", x=NULL, fill=NULL) +
  scale_x_discrete(labels = c("Cluster 1" = "Cluster1", "Cluster 2" = "Cluster2", "Cluster 3" = "Cluster3",
                              "Cluster 4" = "Cluster4", "Cluster 5" = "Cluster5"))+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5,size=16),
        axis.title=element_text(size=18),
        axis.text=element_text(size=15))



######################################
###以下为合并之前的代码
####下载到本地服务器后进行注释
library(ChIPseeker)
hg19 <-  makeTxDbFromGFF("/NAS/luozh/CRC_conservation/reference/gencode.v32lift37.annotation.gtf", format = "gtf")
peak <- readPeakFile("/home/shimw/project/conservation/human_enhancer_DE.bed")

####bed为没有筛选promoter的，csv为筛选
peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000), TxDb=hg19)
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
# Cluster 1 Cluster 2 Cluster 3 Cluster 4 Cluster 5 
# 2602       818      1500      1020      6376 
length(distinct(cluster_PeaktoGene[,c(1,2)])$cluster)


my_comparisons <- list(c("Cluster 1", "Cluster 2"), c("Cluster 1", "Cluster 3"), c("Cluster 1", "Cluster 4"), c("Cluster 1", "Cluster 5"))
ggplot(data = na.omit(cluster_PeaktoGene),aes(x = cluster,y = log2FoldChange,fill=cluster)) +
  geom_violin(adjust = .5) +geom_boxplot(width=0.1, fill="white")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif", tip.length=0.02,method.args = list(alternative = "greater"), method = "t.test") +
  theme_classic(base_family = "sans",base_size = 18,base_line_size = 1.1)+
  scale_fill_manual(values=c("Cluster 1" = "#CB9C7A", "Cluster 2" = "#8696a7", "Cluster 3" = "#CDB97D",
                             "Cluster 4" = "#7b8b6f", "Cluster 5" = "#A59B95"))+
  labs(title = "peak intensity", y="Log2foldchange", x=NULL, fill=NULL) +
  scale_x_discrete(labels = c("Cluster 1" = "Cluster1", "Cluster 2" = "Cluster2", "Cluster 3" = "Cluster3",
                              "Cluster 4" = "Cluster4", "Cluster 5" = "Cluster5"))+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5,size=16),
        axis.title=element_text(size=18),
        axis.text=element_text(size=15))


peakAnno[peakAnno$geneId%in%df_data[df_data$cluster=="Cluster 1",]$gene,]
peakAnno_kk <- annotatePeak(peakAnno[peakAnno$geneId%in%df_data[df_data$cluster=="Cluster 1",]$gene&peakAnno$V9<0.05,], tssRegion=c(-3000, 3000), TxDb=hg19)
plotAnnoPie(peakAnno_kk)

plot_cluster_anno <- function(cluter_type){
  peakAnno_kk <- annotatePeak(peakAnno[peakAnno$geneId%in%df_data[df_data$cluster==cluter_type,]$gene&peakAnno$V9<0.05,], tssRegion=c(-3000, 3000), TxDb=hg19)
  plotAnnoPie(peakAnno_kk)
}
plot_cluster_anno("Cluster 1")
plot_cluster_anno("Cluster 2")
plot_cluster_anno("Cluster 3")
plot_cluster_anno("Cluster 4")






