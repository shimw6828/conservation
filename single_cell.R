###############
# 人和小鼠单细胞聚类图以及marker标记热图Fig1
############ 导入人的数据#####################
library(tidyverse)
library(clusterProfiler)
library(Seurat)
library("data.table")
library("readxl")
# raw_count<- fread("/home/zhangy/project/SMC_colon/data/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz")
# ## 病人信息和细胞信息
# annotation<- fread("/home/zhangy/project/SMC_colon/data/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz")
# 
# load("/home/zhangy/project/SMC_colon/output/0100/0120_cluster.RData")
# seu$Cell_type
# 
# 
# Idents(seu)<-seu$Cell_type
# cell_markers<-lapply(unique(seu$Cell_type), function(x){
#   markers <- FindMarkers(seu, logfc.threshold = 0.01,
#                          group.by = 'Class', 
#                          ident.1 = "Tumor",ident.2 = "Normal",
#                          subset.ident = x)
#   markers <- dplyr::arrange(markers,desc(avg_log2FC))
#   markers$cell_type = x
#   markers = tibble::rownames_to_column(markers,var = "gene")
#   return(markers)
# })
# 
# filter_marker <- function(no_filter_marker){purrr::map(no_filter_marker,function(x){
#   x = x[x$p_val_adj<0.05&abs(x$avg_log2FC)>0.25,]
#   x=tidyr::drop_na(x,p_val_adj)
#   return(x)
# })}
# 
# filter_markers <- filter_marker(cell_markers)
# names(filter_markers) <- unique(seu$Cell_type)
# table(seu$Cell_type[seu$Class=="Tumor"])
# table(seu$Cell_type[seu$Class=="Normal"])
# 
# df_data = read.csv("/home/panxl/CRC2/all_cluater")
# df_data = df_data[,c("Mouse.gene.stable.ID","cluster")]
# names(df_data)[1] = "gene"
# 
# ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "grch37.ensembl.org")
# symboltoid = getBM(attributes=c('ensembl_gene_id', 'external_gene_name'), 
#                    mart = ensembl)
# 
# head(symboltoid)
# 
# seu$orig.ident
# seu@assays[["RNA"]]@counts[1:10,1:10]

load("/home/zhangy/project/conservation/output/single/result/human_mouse_seob.RData")


hh_gene_data<-read.csv("/home/panxl/CRC2/all_cluater")
hh_gene_data$cluster = stringr::str_sub(hh_gene_data$cluster,0,9)
library(biomaRt)
ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "grch37.ensembl.org")
idtosymbol = getBM(attributes=c('ensembl_gene_id', 'external_gene_name'), 
                   filters = 'ensembl_gene_id', 
                   values = hh_gene_data$Gene.stable.ID, 
                   mart = ensembl)
hh_gene_data = hh_gene_data[hh_gene_data$Gene.stable.ID%in%idtosymbol$ensembl_gene_id,]
hh_gene_data$SYMBOL<-idtosymbol[match(hh_gene_data$Gene.stable.ID,idtosymbol$ensembl_gene_id),]$external_gene_name
hh_gene_data <- na.omit(hh_gene_data)
dim(hh_gene_data)
hh_gene_data = hh_gene_data[hh_gene_data$SYMBOL%in%rownames(seu@assays[["RNA"]]@counts),]
cluster_type = "Cluster 1"
calc_cbscore <- function(cluster_type){
  c_score <- apply(seu@assays[["RNA"]]@counts[
    hh_gene_data$SYMBOL[hh_gene_data$cluster==cluster_type],],
    2,
    function(x){
      table(x>0)["TRUE"]/length(hh_gene_data$SYMBOL[hh_gene_data$cluster==cluster_type])
    })
  
  e_score = apply(seu@assays[["RNA"]]@counts,2,function(x){
    sum(x[hh_gene_data$SYMBOL[hh_gene_data$cluster==cluster_type]])/sum(x)
  })
  combind_score = c_score * e_score
  kk = 1/(-log(combind_score))
  return(kk)
}

load("/home/zhangy/project/conservation/output/single/result/human_score.Rdata")  
seu$cluster1_score =calc_cbscore("Cluster 1")
seu$cluster2_score =calc_cbscore("Cluster 2")
seu$cluster3_score =calc_cbscore("Cluster 3")
seu$cluster4_score =calc_cbscore("Cluster 4")
seu$cluster5_score =calc_cbscore("Cluster 5")
meta_df <- seu@meta.data
# rm(seu)
# unique(meta_df[meta_df$Class=="Tumor"&meta_df$Cell_type=="Epithelial cells",]$Cell_subtype)
# 全是CMS分型内
# meta_df[meta_df$Class=="Tumor"&meta_df$Cell_type=="Epithelial cells",]$Cell_type = "Malignant epithelial cell"
meta_df <- meta_df[,c("Cell_type","cluster1_score","cluster2_score",
                      "cluster3_score","cluster4_score","cluster5_score",
                      "Sex","Stage","seurat_clusters","Patient","Class")]


write.csv(meta_df,"/home/shimw/project/conservation/single_cell/meta_df.csv",quote = F)
meta_df<-read.csv("/home/shimw/project/conservation/single_cell/meta_df.csv",header = T,row.names = 1)
meta_df_gather <- tidyr::gather(meta_df,score_type,score,cluster1_score:cluster5_score)
readr::write_csv(meta_df_gather,"/home/shimw/project/conservation/single_cell/meta_df_gather.csv")
meta_df_gather <- readr::read_csv("/home/shimw/project/conservation/single_cell/meta_df_gather.csv")
meta_df_gather$Cell_type <- factor(meta_df_gather$Cell_type,levels = rev(c("Malignant epithelial cells","Epithelial cells","Stromal cells",
                                                                           "B cells","Mast cells","Myeloids",
                                                                           "T cells")))


label_map <- c(
  `cluster1_score`="Cluster 1",
  `cluster2_score`="Cluster 2",
  `cluster3_score`="Cluster 3",
  `cluster4_score`="Cluster 4",
  `cluster5_score`="Cluster 5"
)

p1 <- ggplot(data = na.omit(meta_df_gather),aes(x = Cell_type,y = score,fill=Cell_type)) +
  # geom_violin(adjust = .5) +
  geom_boxplot()+coord_flip()+facet_wrap(~score_type, nrow = 1, scales = "free_x",labeller = as_labeller(label_map))+
  # stat_compare_means(comparisons = my_comparisons,label = "p.signif", tip.length=0.02,method.args = list(alternative = "greater"), method = "t.test") +
  scale_y_continuous(breaks = scales::pretty_breaks(3)) +
  theme_classic(base_family = "sans",base_size = 18,base_line_size = 0.6)+
  scale_fill_manual(values=c('Malignant epithelial cells'="#24325F", 
                             "Adj epithelial cells"="#B8E4F9", #深蓝色
                             'Epithelial cells'="#94643F", # 棕色
                             'Stromal cells'="#D2691E", #黄色
                             'Myeloids'="#FB696B", # 红色
                             "Mast cells"="#7472AF", #紫色
                             'T cells'="#536F2E", # 绿色
                             'B cells'="#E78CCB"))+
  labs(title = "Human", y="CCB score", x=NULL, fill=NULL) +
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5,size=16),
        axis.title=element_text(size=18),
        axis.text=element_text(size=15),
        strip.background = element_blank(),
        panel.spacing = unit(1.5, "lines"))
ggsave("/home/shimw/project/conservation/single_cell/human_crb_box.pdf",p1,width = 13,height = 3)


#####################################################################
##小鼠的同样计算CRB score
rm(list = ls())
load("/home/zhangy/project/conservation/output/single/result/mouse_score.Rdata")  
mm_gene_data<-read.csv("/home/panxl/CRC2/all_cluater")
mm_gene_data$cluster = stringr::str_sub(mm_gene_data$cluster,0,9)
# library(biomaRt)
# ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "feb2014.archive.ensembl.org",ensemblRedirect = FALSE)
# kk = listAttributes(ensembl)
# idtosymbol = getBM(attributes=c('ensembl_gene_id', 'external_gene_name'), 
#                    mart = ensembl)
#不行，只能取e75自己手动下，我把header也顺便改了
idtosymbol <-  readr::read_csv("/home/shimw/project/conservation/mm10_idtosymbol.csv")
mm_gene_data = mm_gene_data[mm_gene_data$Mouse.gene.stable.ID%in%idtosymbol$ensembl_gene_id,]
mm_gene_data$SYMBOL<-idtosymbol[match(mm_gene_data$Mouse.gene.stable.ID,idtosymbol$ensembl_gene_id),]$external_gene_name

mm_gene_data <- na.omit(mm_gene_data)
##一个基因重复了，Mtcp1，cluster 5的基因，都不显著。所以直接去掉不影响
mm_gene_data = mm_gene_data[mm_gene_data$SYMBOL!="Mtcp1",]
mm_gene_data = mm_gene_data[mm_gene_data$SYMBOL%in%rownames(seob@assays[["RNA"]]@counts),]
cluster_type = "Cluster 1"
calc_cbscore <- function(cluster_type){
  c_score <- apply(seob@assays[["RNA"]]@counts[
    mm_gene_data$SYMBOL[mm_gene_data$cluster==cluster_type],],
    2,
    function(x){
      table(x>0)["TRUE"]/length(mm_gene_data$SYMBOL[mm_gene_data$cluster==cluster_type])
    })
  
  e_score = apply(seob@assays[["RNA"]]@counts,2,function(x){
    sum(x[mm_gene_data$SYMBOL[mm_gene_data$cluster==cluster_type]])/sum(x)
  })
  combind_score = c_score * e_score
  kk = 1/(-log(combind_score))
  return(kk)
}
seob$cluster1_score =calc_cbscore("Cluster 1")
seob$cluster2_score =calc_cbscore("Cluster 2")
seob$cluster3_score =calc_cbscore("Cluster 3")
seob$cluster4_score =calc_cbscore("Cluster 4")
seob$cluster5_score =calc_cbscore("Cluster 5")
load("/home/zhangy/project/conservation/output/single/result/mouse_score.Rdata")  
##zhangy已经计算过了直接读入
meta_df <- seob@meta.data

# meta_df[meta_df$tissue=="tumor"&meta_df$Cell_type=="Epithelial cells",]$Cell_type = "Malignant epithelial cell"
# meta_df[meta_df$tissue=="adj"&meta_df$Cell_type=="Malignant epithelial cell",]$Cell_type = "Epithelial cells"
# meta_df[meta_df$tissue=="wt"&meta_df$Cell_type=="Malignant epithelial cell",]$Cell_type = "Epithelial cells"
meta_df[meta_df$Cell_type=="Adj epithelial cell",]$Cell_type = "Epithelial cells"
table(meta_df$tissue,meta_df$Cell_type)

meta_df_gather <- tidyr::gather(meta_df,score_type,score,cluster1_score:cluster5_score)
meta_df_gather <- meta_df_gather[,c("sample","tissue","Cell_type","score_type","score")]
# readr::write_csv(meta_df_gather,"/home/shimw/project/conservation/single_cell/meta_df_mm_gather.csv")
# meta_df_gather <- readr::read_csv("/home/shimw/project/conservation/single_cell/meta_df_mm_gather.csv")
# meta_df_gather = meta_df_gather[meta_df_gather$Cell_type!="Stem cells",]
meta_df_gather$Cell_type <- factor(meta_df_gather$Cell_type,levels = rev(c("Malignant epithelial cells",
                                                                           "Adj epithelial cells",
                                                                           "Epithelial cells","Stromal cells",
                                                                           "B cells","Myeloids",
                                                                           "T cells")))

label_map <- c(
  `cluster1_score`="Cluster 1",
  `cluster2_score`="Cluster 2",
  `cluster3_score`="Cluster 3",
  `cluster4_score`="Cluster 4",
  `cluster5_score`="Cluster 5"
)

p1 <- ggplot(data = na.omit(meta_df_gather),aes(x = Cell_type,y = score,fill=Cell_type)) +
  # geom_violin(adjust = .5) +
  geom_boxplot()+coord_flip()+facet_wrap(~score_type, nrow = 1, scales = "free_x",labeller = as_labeller(label_map))+
  # stat_compare_means(comparisons = my_comparisons,label = "p.signif", tip.length=0.02,method.args = list(alternative = "greater"), method = "t.test") +
  scale_y_continuous(breaks = scales::pretty_breaks(3)) +
  theme_classic(base_family = "sans",base_size = 18,base_line_size = 0.6)+
  scale_fill_manual(values=c('Malignant epithelial cells'="#24325F", 
                             "Adj epithelial cells"="#B8E4F9", #深蓝色
                             'Epithelial cells'="#94643F", # 棕色
                             'Stromal cells'="#D2691E", #黄色
                             'Myeloids'="#FB696B", # 红色
                             "Mast cells"="#7472AF", #紫色
                             'T cells'="#536F2E", # 绿色
                             'B cells'="#E78CCB"))+
  labs(title = "Mouse", y="CCB score", x=NULL, fill=NULL) +
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5,size=16),
        axis.title=element_text(size=18),
        axis.text=element_text(size=15),
        strip.background = element_blank(),
        panel.spacing = unit(1.5, "lines"))
ggsave("/home/shimw/project/conservation/single_cell/mouse_crb_box.pdf",p1,width = 13,height = 3)


table(meta_df_gather$Cell_type)


mm = apply(seob@assays[["RNA"]]@counts,
  2,
  function(x){
    return(table(x>0)["TRUE"])
  })

seob@meta.data$gene_num <- mm
meta_df <- seob@meta.data
ggplot(meta_df,aes(x = Cell_type,y = gene_num)) +geom_boxplot()





#################################
load("/home/zhangy/project/conservation/output/single/human_score.Rdata")
load("/home/zhangy/project/conservation/output/single/mouse_score.Rdata")

human_meta <- seu@meta.data
mouse_meta <- seob@meta.data


human_meta_gather <- tidyr::gather(human_meta,score_type,score,cluster1_score:cluster5_score)
mouse_meta_gather <- tidyr::gather(mouse_meta,score_type,score,cluster1_score:cluster5_score)
mouse_meta_gather[mouse_meta_gather$Cell_type!="adj",]
mouse_meta_gather$Class = if_else(mouse_meta_gather$tissue=="tumor","Tumor","Normal")
human_meta_gather = human_meta_gather[,c("Class","Cell_type","score_type","score")]
mouse_meta_gather = mouse_meta_gather[,c("Class","Cell_type","score_type","score")]
human_meta_gather$species = "human"
mouse_meta_gather$species = "mouse"
meta_gather <- rbind(human_meta_gather, mouse_meta_gather)
meta_gather



# human_meta_gather$Cell_type
# human_meta_gather$Cell_type <- factor(meta_df_gather$Cell_type,levels = rev(c("Malignant epithelial cells","Epithelial cells","Stromal cells",
#                                                                            "B cells","Myeloids",
#                                                                            "T cells")))
or = c("Malignant epithelial cells","Epithelial cells","Stromal cells",
  "B cells","Myeloids",
  "T cells")
cl = "cluster1_score"
sp = "human"
meta_gather[meta_gather$Cell_type=="Malignant epithelial cell",]$Cell_type="Malignant epithelial cells"
meta_gather = meta_gather[meta_gather$Cell_type!="Mast cells",]
library(Ipaper)

is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

plot_box = function(sp, cl, or1,or2){
  or = c(or1,or2)
  tmp_meta <- meta_gather[meta_gather$species==sp&meta_gather$score_type==cl,]
  tmp_meta$Cell_type <- factor(tmp_meta$Cell_type,levels = or)
  tmp_meta = na.omit(tmp_meta)
  kk = t.test(tmp_meta[tmp_meta$Cell_type%in%or1,]$score,tmp_meta[tmp_meta$Cell_type%in%or2,]$score,alternative = "greater")
  tmp_meta <- tmp_meta %>%
    group_by(Cell_type) %>%
    mutate(outlier = ifelse(is_outlier(score), 0, score))
  or1_ymax <-  max(tmp_meta[tmp_meta$Cell_type%in%or1,]$outlier)
  or2_ymax <-  max(tmp_meta[tmp_meta$Cell_type%in%or2,]$outlier)
  or1_n = length(or1)
  or2_n = length(or2)
  p_sig = case_when(kk$p.value>0.05 ~ "ns",
            kk$p.value>0.01 ~ "*",
            kk$p.value>0.001 ~ "**",
            kk$p.value>0.0001 ~ "***",
            TRUE ~ "****")
  hhh = or1_ymax/20
  p = ggplot(data = na.omit(tmp_meta),aes(x = Cell_type,y = score,fill=Cell_type)) +
    # geom_violin(adjust = .5) +
    geom_boxplot2(width.errorbar = 0.5)+
    # geom_boxplot()+
    geom_segment(aes(x=1, y=or1_ymax+hhh, xend=or1_n, yend=or1_ymax+hhh),size=0.1)+
    geom_segment(aes(x=1+or1_n, y=or2_ymax+hhh, xend=or1_n+or2_n, yend=or2_ymax+hhh),size=0.1)+
    geom_segment(aes(x=((1+or1_n)/2), y=or1_ymax+hhh, xend=((1+or1_n)/2), yend=or1_ymax+(hhh*2)),size=0.1)+
    geom_segment(aes(x=((1+or2_n)/2)+or1_n, y=or2_ymax+hhh, xend=((1+or2_n)/2)+or1_n, yend=or1_ymax+(hhh*2)),size=0.1)+
    geom_segment(aes(x=((1+or1_n)/2), y=or1_ymax+(hhh*2), xend=((1+or2_n)/2)+or1_n, yend=or1_ymax+(hhh*2)),size=0.1)+
    annotate("text", x=(1+(or1_n+or2_n)/2+or1_n)/2, y=or1_ymax+(hhh*3), label=p_sig,size=5)+
    # stat_compare_means(comparisons = my_comparisons,label = "p.signif", tip.length=0.02,method.args = list(alternative = "greater"), method = "t.test") +
    scale_y_continuous(breaks = scales::pretty_breaks(3)) +
    theme_classic(base_family = "sans",base_size = 18,base_line_size = 0.6)+
    scale_fill_manual(values=c('Malignant epithelial cells'="#24325F", 
                               # "Adj epithelial cells"="#B8E4F9", #深蓝色
                               'Epithelial cells'="#94643F", # 棕色
                               'Stromal cells'="#D2691E", #黄色
                               'Myeloids'="#FB696B", # 红色
                               'T cells'="#536F2E", # 绿色
                               'B cells'="#E78CCB"))+    
    scale_x_discrete(labels = c('Malignant epithelial cells' = "Mc", 
                                "Epithelial cells" = "Ec", 
                                "Stromal cells" = "Sc",
                                "Myeloids" = "Ms", 
                                "T cells" = "Tc",
                                "B cells" = "Bc"))+
    labs(title = sp, y=NULL, x=NULL, fill=NULL) +
    theme(legend.position = "none",
          axis.title=element_text(size=12),
          axis.text=element_text(size=12))
  return(p)
}


p_1_h <- plot_box("human","cluster1_score",
                  c("Malignant epithelial cells","Stromal cells"),
                  c("Epithelial cells","B cells","Myeloids","T cells"))+
  labs(title = NULL, y=NULL, x=NULL, fill=NULL) 
p_1_m <- plot_box("mouse","cluster1_score",
                  c("Malignant epithelial cells","Stromal cells"),
                  c("Epithelial cells","B cells","Myeloids","T cells"))+
  labs(title = NULL, y=NULL, x=NULL, fill=NULL) 


p_2_h <- plot_box("human","cluster2_score",
                  c("B cells","Myeloids","T cells"),
                  c("Malignant epithelial cells","Epithelial cells","Stromal cells")
                  )+
  labs(title = NULL, y=NULL, x=NULL, fill=NULL) 
p_2_m <- plot_box("mouse","cluster2_score",
                  c("B cells","Myeloids","T cells"),
                  c("Malignant epithelial cells","Epithelial cells","Stromal cells")
                  )+
  labs(title = NULL, y=NULL, x=NULL, fill=NULL)

p_3_h <- plot_box("human","cluster3_score",
                  c("Epithelial cells"),
                  c("Malignant epithelial cells","Stromal cells","B cells","Myeloids","T cells")
)+
  labs(title = NULL, y=NULL, x=NULL, fill=NULL) 
p_3_m <- plot_box("mouse","cluster3_score",
                  c("Epithelial cells"),
                  c("Malignant epithelial cells","Stromal cells","B cells","Myeloids","T cells")
)+
  labs(title = NULL, y=NULL, x=NULL, fill=NULL)

p_4_h <- plot_box("human","cluster4_score",
                  c("Malignant epithelial cells","Epithelial cells","Stromal cells"),
                  c("B cells","Myeloids","T cells")
)+
  labs(title = NULL, y=NULL, x=NULL, fill=NULL) 
p_4_m <- plot_box("mouse","cluster4_score",
                  c("Malignant epithelial cells","Epithelial cells","Stromal cells"),
                  c("B cells","Myeloids","T cells")
)+
  labs(title = NULL, y=NULL, x=NULL, fill=NULL)

p_5_h <- plot_box("human","cluster5_score",
                  c("Malignant epithelial cells","Stromal cells"),
                  c("Epithelial cells","B cells","Myeloids","T cells")
)+
  labs(title = NULL, y=NULL, x=NULL, fill=NULL) 
p_5_m <- plot_box("mouse","cluster5_score",
                  c("Malignant epithelial cells","Stromal cells"),
                  c("Epithelial cells","B cells","Myeloids","T cells")
)+
  labs(title = NULL, y=NULL, x=NULL, fill=NULL)

library("ggpubr")
ggarrange(p_1_h, p_1_m,
             p_2_h, p_2_m,
             p_3_h, p_3_m,
             p_4_h, p_4_m,
             p_5_h, p_5_m,
             ncol = 2, nrow = 5)


