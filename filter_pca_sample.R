query.exp.hg19 <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  platform = "Illumina HiSeq", 
  file.type  = "normalized_results",
  experimental.strategy = "RNA-Seq",
  legacy = TRUE
)
GDCdownload(query.exp.hg19)

GDCdownload(query.exp.hg19 ,files.per.chunk = 30)


######################################################################

#readr::write_csv(metadata, "/home/shimw/project/Deconvolution/tcga_metadata.csv")
metadata = readr::read_csv( "/home/shimw/project/Deconvolution/tcga_metadata.csv")

cancer_count <- readr::read_csv("/home/shimw/project/conservation/tcgaCount.csv")
cancer_count <- tibble::column_to_rownames(cancer_count, var = "gene")
cancer_type = "COAD"
metadata <- metadata%>%
  dplyr::filter(abbreviation==cancer_type)%>%
  tidyr::drop_na()
metadata = metadata[metadata$sample%in%names(cancer_count),]
metadata = metadata[metadata$sample_type%in%c("Primary Tumor","Solid Tissue Normal"),]
table(metadata$sample_type)
COAD_count <- cancer_count[,metadata$sample]
rm(cancer_count)

write.csv(COAD_count,"/home/shimw/project/conservation/COAD_hg19_count.csv",quote = F)
readr::write_csv(metadata,"/home/shimw/project/conservation/COAD_hg19_metadata.csv")

#########################################################################
#oursdata
Gene_exp=read.csv("/home/panxl/CRC/RNA_seq/wuda_Gene_exp.csv")
rownames(Gene_exp)=Gene_exp$X
Gene_exp = Gene_exp[, -1]
#样本信息表
group1 <- c(rep("tumor",72),rep("native",72))
sample_info1 <- data.frame(Group = factor(group1, levels = c( "native","tumor")))
rownames(sample_info1) <- colnames(Gene_exp)
Group<-sample_info1$Group
sample_info1$DataSet="ourdata"
sample_info1$sample=rownames(sample_info1)
COAD_count <- read.csv("/home/shimw/project/conservation/COAD_hg19_count.csv",row.names = 1,check.names = F)
metadata = readr::read_csv("/home/shimw/project/conservation/COAD_hg19_metadata.csv")
metadata$Group <- ifelse(metadata$sample_type=="Primary Tumor","tumor","native")
metadata$DataSet = "TCGA"
tmpmeta <- metadata[,c("Group","DataSet","sample")]
sample_info <- rbind(sample_info1,tmpmeta)
##修一下基因名
row.names(COAD_count) <- stringr::str_split_fixed(row.names(COAD_count),pattern = "\\.",2)[,1]
gene=intersect(rownames(Gene_exp),rownames(COAD_count))
length(gene)#57000
ourdata=Gene_exp[gene,]
tcga=COAD_count[gene,]
all_data=cbind(ourdata,tcga)#652


###########################################################################
#相关性筛选


filter_cor <- function(cancer_group){
  tumor_cortest <- rcorr(as.matrix(all_data[sample_info[sample_info$Group==cancer_group,"sample"]]), type = "pearson")
  tumor_cortest=tumor_cortest$r
  tumor_cortest_sum <- apply(tumor_cortest, 1, sum)%>%
    as.data.frame()%>%
    set_colnames("r")
  tumor_cortest_sum$r = tumor_cortest_sum$r-1
  kk =mean(tumor_cortest_sum$r) - 1.5*sd(tumor_cortest_sum$r)
  tumor_cortest_sum[tumor_cortest_sum$r>kk,]
  row.names(sample_info) <- sample_info$sample
  return(row.names(tumor_cortest_sum)[tumor_cortest_sum$r>kk])
}




length(filter_cor("native"))
length(filter_cor("tumor"))

table(sample_info$Group)




barcode <- c(filter_cor("native"),filter_cor("tumor"))
row.names(sample_info) <- sample_info$sample
tmp_data <- all_data[,barcode]
tmp_info <- sample_info[barcode,]


dds <- DESeqDataSetFromMatrix(countData = tmp_data,
                              colData = tmp_info,
                              design = ~ Group)#57000
##筛选
dds<-dds[rowSums(counts(dds))>50,]#33559
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds)
vsd <- vst(dds, blind=FALSE)
human_wgcna_data<-as.data.frame(t(assay(vsd)))

##
dge <- DGEList(counts = tmp_data)
dge <- calcNormFactors(dge)
design <- model.matrix(~tmp_info$Group)
colnames(design) <- levels(tmp_info$Group)
rownames(design) <- colnames(tmp_data)
v <- voom(dge, design)
vomm_count = v$E

tmp_info$DataSet
# data.ComBat = ComBat(t(human_wgcna_data), batch = vsd$DataSet, mod = model.matrix(~vsd$Group))
data.ComBat = ComBat(vomm_count, batch = tmp_info$DataSet, mod = model.matrix(~tmp_info$Group))
##################################
data.ComBat <- removeBatchEffect(vomm_count,
batch = tmp_info$DataSet,
design =  model.matrix(~tmp_info$Group))



pca <- prcomp(t(data.ComBat))



library("FactoMineR")
library("factoextra")
fviz_pca_ind(pca,
             #axes = c(2,3),
             geom.ind = "point",
             col.ind = tmp_info$Group ,
             addEllipses = TRUE,
             legend.title = "Groups",
             invisible = c("ind.sup")
) +
  scale_color_manual(values=c("#930708", "#073c68")) +
  xlab(paste0(PCA, ": ", round(percentVar[A] * 100),"% variance")) +
  ylab(paste0(PCB, ": ", round(percentVar[B] * 100),"% variance")) +
  theme_classic() +
  theme(axis.title = element_text(size=28),
        axis.text = element_text(size=24),
        legend.title = element_blank(),
        legend.text = element_text(size=24),
        legend.key.size = unit(0.6,"cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
  )

#####################################################
##迭代筛选
pcaData = pca$x
pcaData = merge(pcaData[,1:10], tmp_info, by="row.names")
rownames(pcaData) = pcaData$Row.names
pcaData = pcaData[,-c(1)]
#pcaData = merge(pcaData, sampleinfo, by="row.names")
pcaData


filter_sd<- function(pData){
  #pData=pcaData
  normal_pd=pData[pData$Group=="native",]
  sd_num <- apply(normal_pd[,1:5], 2, sd)
  mean_num <- apply(normal_pd[,1:5], 2, mean)
  kk = sweep(data.matrix(normal_pd[,1:5]), 2, 3*sd_num, FUN = "-")
  mm = sweep(data.matrix(normal_pd[,1:5]), 2, 3*sd_num, FUN = "+")
  tt = sweep(mm, 2, mean_num, FUN = ">")&sweep(kk, 2, mean_num, FUN = "<")
  normal_keep = normal_pd$sample[apply(tt,1,all)]
  
  tumor_pd=pData[pData$Group=="tumor",]
  sd_num <- apply(tumor_pd[,1:5], 2, sd)
  mean_num <- apply(tumor_pd[,1:5], 2, mean)
  kk = sweep(data.matrix(tumor_pd[,1:5]), 2, 3*sd_num, FUN = "-")
  mm = sweep(data.matrix(tumor_pd[,1:5]), 2, 3*sd_num, FUN = "+")
  tt = sweep(mm, 2, mean_num, FUN = ">")&sweep(kk, 2, mean_num, FUN = "<")
  tumor_keep = tumor_pd$sample[apply(tt,1,all)]
  return(c(normal_keep,tumor_keep))
}

unfilter_df <- data.ComBat
unfilter_tmp_info=tmp_info
for (i in seq(1,5)) {
  print(i)
  pca <- prcomp(t(unfilter_df))
  pcaData = pca$x
  pcaData = merge(pcaData[,1:10], unfilter_tmp_info, by="row.names")
  rownames(pcaData) = pcaData$Row.names
  pcaData = pcaData[,-c(1)]
  #pcaData = merge(pcaData, sampleinfo, by="row.names")
  keep_sample <- filter_sd(pcaData)
  unfilter_df=unfilter_df[,keep_sample]
  unfilter_tmp_info = unfilter_tmp_info[keep_sample,]
  print(nrow(pcaData)-length(keep_sample))
  if (nrow(pcaData)-length(keep_sample)==0) {
    break
  }
}

pca <- prcomp(t(unfilter_df))

fviz_pca_ind(pca,
             #axes = c(2,3),
             geom.ind = "point",
             col.ind = unfilter_tmp_info$Group ,
             addEllipses = TRUE,
             legend.title = "Groups",
             invisible = c("ind.sup")
) +
  scale_color_manual(values=c("#930708", "#073c68")) +
  xlab(paste0(PCA, ": ", round(percentVar[A] * 100),"% variance")) +
  ylab(paste0(PCB, ": ", round(percentVar[B] * 100),"% variance")) +
  theme_classic() +
  theme(axis.title = element_text(size=28),
        axis.text = element_text(size=24),
        legend.title = element_blank(),
        legend.text = element_text(size=24),
        legend.key.size = unit(0.6,"cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
  )

length(filter_cor("native"))
length(filter_cor("tumor"))
table(unfilter_tmp_info$Group[unfilter_tmp_info$DataSet=="ourdata"])
table(sample_info$Group[sample_info$DataSet=="ourdata"])
table(unfilter_tmp_info$Group[unfilter_tmp_info$DataSet=="TCGA"])
table(sample_info$Group[sample_info$DataSet=="TCGA"])



#########################################################
library(readxl)

##########
total_data = read.csv("/home/panxl/CRC/all_mouse_data/summary_count.txt",row.names = 1)
#sample_info
library(readxl)
change_DataSet <- read_excel("/home/panxl/CRC/all_mouse_data/change_DataSet.xlsx")
change_DataSet=na.omit(change_DataSet)
colnames(change_DataSet)=c("SRR_id","DataSet")
group=c(rep("Normal",31,),rep("Tumor",41),rep("Normal",3,),rep("Tumor",6))#normal-7 tumor-5  #72  84
sample_info <- data.frame(Group = factor(group, levels = c("Normal", "Tumor")))
sample_info$SRR_id=colnames(total_data)

sample_info=left_join(sample_info,change_DataSet,by = "SRR_id")
rownames(sample_info)=colnames(total_data)
sample_info$DataSet[c(73:81)]="ourdata"
sample_info$DataSet=factor(sample_info$DataSet,levels = c("GSE102648","GSE111615","GSE116374","GSE124029",
                                                          "GSE132814","GSE136029","GSE142154",
                                                          "GSE152032","GSE155777","GSE158276","GSE164960",
                                                          "GSE166277","GSE57533","GSE60837",
                                                          "GSE97013","GSE98496","ourdata"))
sample_info=sample_info[order(sample_info$DataSet),]

cor_normal <- cor(total_data[sample_info[sample_info$Group==cancer_group,"SRR_id"]], method = 'pearson')
library(ggcorrplot)
ggcorrplot(tumor_cortest,hc.order=TRUE,outline.color="white",
           type="lower",colors = c("#6D9EC1", "white", "#E46726"),
           ggtheme = ggplot2::theme_void())



################################################
filter_cor <- function(cancer_group){
  tumor_cortest <- rcorr(as.matrix(total_data[sample_info[sample_info$Group==cancer_group,"SRR_id"]]), type = "pearson")
  tumor_cortest=tumor_cortest$r
  tumor_cortest_sum <- apply(tumor_cortest, 1, sum)%>%
    as.data.frame()%>%
    set_colnames("r")
  tumor_cortest_sum$r = tumor_cortest_sum$r-1
  kk =mean(tumor_cortest_sum$r) - 1.5*sd(tumor_cortest_sum$r)
  tumor_cortest_sum[tumor_cortest_sum$r>kk,]
  row.names(sample_info) <- sample_info$sample
  return(row.names(tumor_cortest_sum)[tumor_cortest_sum$r>kk])
}
length(filter_cor("Normal"))
length(filter_cor("Tumor"))
table(sample_info$Group)
barcode = c(filter_cor("Normal"), filter_cor("Tumor"))

tmp_data <- total_data[,barcode]
tmp_info <- sample_info[barcode,]
tmp_info$DataSet <- factor(as.character(tmp_info$DataSet))


dge <- DGEList(counts = tmp_data)
dge <- calcNormFactors(dge)
design <- model.matrix(~tmp_info$Group)
colnames(design) <- levels(tmp_info$Group)
rownames(design) <- colnames(tmp_data)
v <- voom(dge, design)
mouse_vomm_count = v$E
mouse_data.ComBat = ComBat(mouse_vomm_count, batch = tmp_info$DataSet, mod = model.matrix(~tmp_info$Group))

pca <- prcomp(t(mouse_data.ComBat))


library("FactoMineR")
library("factoextra")
fviz_pca_ind(pca,
             #axes = c(2,3),
             geom.ind = "point",
             col.ind = tmp_info$Group ,
             addEllipses = TRUE,
             legend.title = "Groups",
             invisible = c("ind.sup")
) +
  scale_color_manual(values=c("#930708", "#073c68")) +
  xlab(paste0(PCA, ": ", round(percentVar[A] * 100),"% variance")) +
  ylab(paste0(PCB, ": ", round(percentVar[B] * 100),"% variance")) +
  theme_classic() +
  theme(axis.title = element_text(size=28),
        axis.text = element_text(size=24),
        legend.title = element_blank(),
        legend.text = element_text(size=24),
        legend.key.size = unit(0.6,"cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
  )



filter_sd<- function(pData){
  #pData=pcaData
  normal_pd=pData[pData$Group=="Normal",]
  sd_num <- apply(normal_pd[,1:5], 2, sd)
  mean_num <- apply(normal_pd[,1:5], 2, mean)
  kk = sweep(data.matrix(normal_pd[,1:5]), 2, 3*sd_num, FUN = "-")
  mm = sweep(data.matrix(normal_pd[,1:5]), 2, 3*sd_num, FUN = "+")
  tt = sweep(mm, 2, mean_num, FUN = ">")&sweep(kk, 2, mean_num, FUN = "<")
  normal_keep = normal_pd$SRR_id[apply(tt,1,all)]
  
  tumor_pd=pData[pData$Group=="Tumor",]
  sd_num <- apply(tumor_pd[,1:5], 2, sd)
  mean_num <- apply(tumor_pd[,1:5], 2, mean)
  kk = sweep(data.matrix(tumor_pd[,1:5]), 2, 3*sd_num, FUN = "-")
  mm = sweep(data.matrix(tumor_pd[,1:5]), 2, 3*sd_num, FUN = "+")
  tt = sweep(mm, 2, mean_num, FUN = ">")&sweep(kk, 2, mean_num, FUN = "<")
  tumor_keep = tumor_pd$SRR_id[apply(tt,1,all)]
  return(c(normal_keep,tumor_keep))
}

unfilter_df <- mouse_data.ComBat
unfilter_tmp_info=tmp_info
for (i in seq(1,5)) {
  print(i)
  pca <- prcomp(t(unfilter_df))
  pcaData = pca$x
  pcaData = merge(pcaData[,1:10], unfilter_tmp_info, by="row.names")
  rownames(pcaData) = pcaData$Row.names
  pcaData = pcaData[,-c(1)]
  #pcaData = merge(pcaData, sampleinfo, by="row.names")
  keep_sample <- filter_sd(pcaData)
  unfilter_df=unfilter_df[,keep_sample]
  unfilter_tmp_info = unfilter_tmp_info[keep_sample,]
  print(nrow(pcaData)-length(keep_sample))
  if (nrow(pcaData)-length(keep_sample)==0) {
    break
  }
}

pca <- prcomp(t(unfilter_df))

fviz_pca_ind(pca,
             #axes = c(2,3),
             geom.ind = "point",
             col.ind = unfilter_tmp_info$Group ,
             addEllipses = TRUE,
             legend.title = "Groups",
             invisible = c("ind.sup")
) +
  scale_color_manual(values=c("#930708", "#073c68")) +
  xlab(paste0(PCA, ": ", round(percentVar[A] * 100),"% variance")) +
  ylab(paste0(PCB, ": ", round(percentVar[B] * 100),"% variance")) +
  theme_classic() +
  theme(axis.title = element_text(size=28),
        axis.text = element_text(size=24),
        legend.title = element_blank(),
        legend.text = element_text(size=24),
        legend.key.size = unit(0.6,"cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
  )

table(unfilter_tmp_info$Group)
table(tmp_info$Group)



heatmap_data = cor(mouse_data.ComBat)
annotation_col = data.frame(
  dd = factor(sample_info[,c("Group")]))
rownames(annotation_col) = rownames(sample_info)
library(pheatmap)
out=pheatmap(heatmap_data, annotation_col = annotation_col, show_rownames = F, show_colnames = T)



######################################################################
##汇总

library(psych)
filter_cor <- function(countdata,metadata,cancer_group){
# tumor_cortest <- rcorr(as.matrix(countdata[,row.names(metadata[metadata$Group==cancer_group,])]), type = "spearman")
# normal_cortest <- rcorr(as.matrix(countdata[,row.names(metadata[metadata$Group!=cancer_group,])]), type = "spearman")

tumor_cortest = rcorr(as.matrix(countdata[,row.names(metadata[metadata$Group==cancer_group,])]), type = "spearman")
normal_cortest = corr.test(countdata[,row.names(metadata[metadata$Group==cancer_group,])], countdata[,row.names(metadata[metadata$Group!=cancer_group,])], method= "spearman")
normal_cortest=normal_cortest$r
tumor_cortest=tumor_cortest$r

tumor_cortest_sum <- apply(tumor_cortest, 1, sum)%>%
  as.data.frame()%>%
  set_colnames("r")
tumor_cortest_mean <- apply(tumor_cortest, 1, mean)%>%
as.data.frame()%>%
set_colnames("r")
normal_cortest_mean <- apply(normal_cortest, 1, mean)%>%
  as.data.frame()%>%
  set_colnames("r")
kk = tumor_cortest_mean$r - normal_cortest_mean$r
mm =mean(kk) - 3*sd(kk)
ss = mean(tumor_cortest_sum$r) - 3*sd(tumor_cortest_sum$r)

return(row.names(tumor_cortest_sum)[kk>mm&tumor_cortest_sum$r>ss])
}





filter_sd<- function(pData,normal,tumor){
  #pData=pcaData
  normal_pd=pData[pData$Group==normal,]
  sd_num <- apply(normal_pd[,1:5], 2, sd)
  mean_num <- apply(normal_pd[,1:5], 2, mean)
  kk = sweep(data.matrix(normal_pd[,1:5]), 2, 3*sd_num, FUN = "-")
  mm = sweep(data.matrix(normal_pd[,1:5]), 2, 3*sd_num, FUN = "+")
  tt = sweep(mm, 2, mean_num, FUN = ">")&sweep(kk, 2, mean_num, FUN = "<")
  normal_keep = row.names(normal_pd)[apply(tt,1,all)]
  
  tumor_pd=pData[pData$Group==tumor,]
  sd_num <- apply(tumor_pd[,1:5], 2, sd)
  mean_num <- apply(tumor_pd[,1:5], 2, mean)
  kk = sweep(data.matrix(tumor_pd[,1:5]), 2, 3*sd_num, FUN = "-")
  mm = sweep(data.matrix(tumor_pd[,1:5]), 2, 3*sd_num, FUN = "+")
  tt = sweep(mm, 2, mean_num, FUN = ">")&sweep(kk, 2, mean_num, FUN = "<")
  tumor_keep = row.names(tumor_pd)[apply(tt,1,all)]
  return(c(normal_keep,tumor_keep))
}

library(Hmisc)
library(magrittr)

filter_sample <- function(countdata, mdata, normal = "normal",tumor = "tumor"){
  tmp_data <- countdata
  # tmp_data = tmp_data[rowSums(tmp_data)>ncol(tmp_data)*0.8,]
  dim(tmp_data)
  tmp_info <- mdata
  tmp_info$DataSet <- factor(as.character(tmp_info$DataSet))
  dge <- DGEList(counts = tmp_data)
  dge <- calcNormFactors(dge)
  cut = which(apply(cpm(dge), 1, max)<1)
  dge <- dge[-cut,]
  design <- model.matrix(~tmp_info$Group)
  colnames(design) <- levels(tmp_info$Group)
  rownames(design) <- colnames(tmp_data)
  v <- voom(dge,design)
  vomm_count = v$E
  unfilter_df = ComBat(vomm_count, batch = tmp_info$DataSet, mod = model.matrix(~tmp_info$Group))
  unfilter_tmp_info=tmp_info
  for (i in seq(1,5)) {
    print(i)
    old_sample <- ncol(unfilter_df)
    tumor_keep_sample <- filter_cor(unfilter_df,unfilter_tmp_info,tumor)
    normal_keep_sample <- filter_cor(unfilter_df,unfilter_tmp_info,normal)
    keep_sample <- c(normal_keep_sample,tumor_keep_sample)
    unfilter_tmp_info = unfilter_tmp_info[keep_sample,]
    unfilter_df = ComBat(vomm_count[,keep_sample], batch = unfilter_tmp_info$DataSet, mod = model.matrix(~unfilter_tmp_info$Group))
    print(old_sample-length(keep_sample))
    if (old_sample-length(keep_sample)==0) {
      break
    }
  }
  # annotation_col = data.frame(
  #   dd = factor(unfilter_tmp_info[,c("Group")]))
  # rownames(annotation_col) = rownames(unfilter_tmp_info)
  # hh=cor(unfilter_df, method = "spearman")
  # pheatmap::pheatmap(hh, annotation_col = annotation_col, show_rownames = F, show_colnames = F, clustering_method="ward.D")
  # 
  # 
  return(unfilter_df)
}

##############################human
###human
##整理数据
COAD_count <- read.csv("/home/shimw/project/conservation/COAD_hg19_count.csv",row.names = 1,check.names = F)
metadata = readr::read_csv("/home/shimw/project/conservation/COAD_hg19_metadata.csv")
Gene_exp=read.csv("/home/panxl/CRC/RNA_seq/wuda_Gene_exp.csv")
rownames(Gene_exp)=Gene_exp$X
Gene_exp = Gene_exp[, -1]
#样本信息表
group1 <- c(rep("tumor",72),rep("native",72))
sample_info1 <- data.frame(Group = factor(group1, levels = c( "native","tumor")))
rownames(sample_info1) <- colnames(Gene_exp)
Group<-sample_info1$Group
sample_info1$DataSet="ourdata"
sample_info1$sample=rownames(sample_info1)
metadata$Group <- ifelse(metadata$sample_type=="Primary Tumor","tumor","native")
metadata$DataSet = "TCGA"
tmpmeta <- metadata[,c("Group","DataSet","sample")]
sample_info <- rbind(sample_info1,tmpmeta)
row.names(sample_info) = sample_info$sample
##修一下基因名
row.names(COAD_count) <- stringr::str_split_fixed(row.names(COAD_count),pattern = "\\.",2)[,1]
gene=intersect(rownames(Gene_exp),rownames(COAD_count))
length(gene)#57000
ourdata=Gene_exp[gene,]
tcga=COAD_count[gene,]
all_data=cbind(ourdata,tcga)#652
# countdata=all_data
# mdata=sample_info
# normal = "native"
# tumor = "tumor"
human_data <- filter_sample(all_data,sample_info,normal = "native",tumor = "tumor")

human_info <- sample_info[colnames(human_data),]
# write.csv(as.data.frame(human_data),"/home/shimw/project/conservation/filter_human_data.csv",quote = F)
# write.csv(human_info,"/home/shimw/project/conservation/filter_human_info.csv",quote = F)
human_data <- read.csv("/home/shimw/project/conservation/filter_human_data.csv",row.names = 1, check.names = F)
human_info <- read.csv("/home/shimw/project/conservation/filter_human_info.csv",row.names = 1, check.names = F)


annotation_col = data.frame(
  dd = factor(human_info[,c("Group")]),
  dataset = factor(human_info[,c("DataSet")]))
rownames(annotation_col) = rownames(human_info)
hh=cor(human_data, method = "spearman")
pheatmap::pheatmap(hh, annotation_col = annotation_col, show_rownames = F, show_colnames = F, clustering_method="ward.D")



pca <- prcomp(t(human_data))
percentVar <- pca$sdev^2 / sum(pca$sdev^2)
PCA = "PC1"
PCB = "PC2"
A=as.numeric(gsub("PC", "", PCA))
B=as.numeric(gsub("PC", "", PCB))


fviz_pca_ind(pca,
             #axes = c(2,3),
             geom.ind = "point",
             col.ind = human_info$Group ,
             addEllipses = TRUE,
             legend.title = "Groups",
             invisible = c("ind.sup")
) +
  scale_color_manual(values=c("#930708", "#073c68")) +
  xlab(paste0(PCA, ": ", round(percentVar[A] * 100),"% variance")) +
  ylab(paste0(PCB, ": ", round(percentVar[B] * 100),"% variance")) +
  theme_classic() +
  theme(axis.title = element_text(size=28),
        axis.text = element_text(size=24),
        legend.title = element_blank(),
        legend.text = element_text(size=24),
        legend.key.size = unit(0.6,"cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
  )
##########################################
##mouse
total_data = read.csv("/home/panxl/CRC/all_mouse_data/summary_count.txt",row.names = 1)
#sample_info
dim(total_data)
library(readxl)
change_DataSet <- read_excel("/home/panxl/CRC/all_mouse_data/change_DataSet.xlsx")
change_DataSet=na.omit(change_DataSet)
colnames(change_DataSet)=c("SRR_id","DataSet")
group=c(rep("Normal",31,),rep("Tumor",41),rep("Normal",3,),rep("Tumor",6))#normal-7 tumor-5  #72  84
sample_info <- data.frame(Group = factor(group, levels = c("Normal", "Tumor")))
sample_info$SRR_id=colnames(total_data)

sample_info=left_join(sample_info,change_DataSet,by = "SRR_id")
rownames(sample_info)=colnames(total_data)
sample_info$DataSet[c(73:81)]="ourdata"
sample_info$DataSet=factor(sample_info$DataSet,levels = c("GSE102648","GSE111615","GSE116374","GSE124029",
                                                          "GSE132814","GSE136029","GSE142154",
                                                          "GSE152032","GSE155777","GSE158276","GSE164960",
                                                          "GSE166277","GSE57533","GSE60837",
                                                          "GSE97013","GSE98496","ourdata"))
# sample_info=sample_info[order(sample_info$DataSet),]
# SRR7447026
# SRR7447028
# SRR7447027
###时间不对，去掉
sample_info = sample_info[!sample_info$SRR_id%in%c("SRR7447026","SRR7447027","SRR7447028"),]
total_data = total_data[,sample_info$SRR_id]
dim(total_data)
row.names(sample_info) = sample_info$SRR_id
# countdata = countdata[rowSums(countdata)>ncol(countdata)*0.5,]
sample_info$DataSet <- factor(as.character(sample_info$DataSet))
# countdata=total_data
# mdata=sample_info
# normal = "Normal"
# tumor = "Tumor"
mouse_data <- filter_sample(total_data,sample_info,normal = "Normal",tumor = "Tumor")

mouse_info <- sample_info[colnames(mouse_data),]
# write.csv(as.data.frame(mouse_data),"/home/shimw/project/conservation/filter_mouse_data.csv",quote = F)
# write.csv(mouse_info,"/home/shimw/project/conservation/filter_mouse_info.csv",quote = F)
mouse_data <- read.csv("/home/shimw/project/conservation/filter_mouse_data.csv",row.names = 1, check.names = F)
mouse_info <- read.csv("/home/shimw/project/conservation/filter_mouse_info.csv",row.names = 1)

annotation_col = data.frame(
  dd = factor(mouse_info[,c("Group")]))
rownames(annotation_col) = rownames(mouse_info)
hh=cor(mouse_data, method = "spearman")
pheatmap::pheatmap(hh, annotation_col = annotation_col, show_rownames = F, clustering_method="ward.D")


pca <- prcomp(t(mouse_data))
percentVar <- pca$sdev^2 / sum(pca$sdev^2)
PCA = "PC1"
PCB = "PC2"
A=as.numeric(gsub("PC", "", PCA))
B=as.numeric(gsub("PC", "", PCB))

fviz_pca_ind(pca,
    #axes = c(2,3),
             geom.ind = "point",
             col.ind = mouse_info$Group ,
             addEllipses = TRUE,
             legend.title = "Groups",
             invisible = c("ind.sup")
) +
  scale_color_manual(values=c("#930708", "#073c68")) +
  xlab(paste0(PCA, ": ", round(percentVar[A] * 100),"% variance")) +
  ylab(paste0(PCB, ": ", round(percentVar[B] * 100),"% variance")) +
  theme_classic() +
  theme(axis.title = element_text(size=28),
        axis.text = element_text(size=24),
        legend.title = element_blank(),
        legend.text = element_text(size=24),
        legend.key.size = unit(0.6,"cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
  )












HOM_one2one=read.csv("/home/panxl/CRC/RNA_seq/HOM_ONE2ONE.txt",sep=",")
HOM_one2one=HOM_one2one[HOM_one2one$Mouse.homology.type=="ortholog_one2one",]
table(HOM_one2one$Gene.stable.ID%in%row.names(human_data)&HOM_one2one$Mouse.gene.stable.ID%in%row.names(mouse_data))





library(WGCNA)
data(BrainLists);
listGenes = unique(as.character(BrainLists[,1]))
set.seed(100)
geneR = sort(sample(listGenes,2000))
categories = sort(rep(standardColors(10),200))
categories[sample(1:2000,200)] = "grey"
file1 = tempfile();
file2 = tempfile();
write(c("TESTLIST1",geneR[300:400], sep="\n"), file1)
write(c("TESTLIST2",geneR[800:1000],sep="\n"), file2)

testResults = userListEnrichment(
  geneR, labelR=categories, 
  fnIn=c(file1, file2),
  catNmIn=c("TEST1","TEST2"), 
  nameOut = NULL, useBrainLists=TRUE, omitCategories ="grey")










