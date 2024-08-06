rm(list = ls()); gc()
library(GenomicFeatures)
library(GenomicRanges)
library(ChIPseeker)
library(dplyr)
hg19 <- makeTxDbFromGFF("/NAS/luozh/CRC_conservation/reference/gencode.v32lift37.annotation.gtf", format = "gtf")
mm10 <- makeTxDbFromGFF("/home/shimw/project/conservation/gencode.vM20.basic.annotation.gtf", format = "gtf")
df_data = read.csv("/home/panxl/CRC2/all_cluater")
cluster1_genes = df_data[df_data$cluster == "Cluster 1: Up-regulated in both species", "Gene.stable.ID"] %>% as.vector()
cluster2_genes = df_data[df_data$cluster == "Cluster 2: Up-regulated in mouse reverse in human", "Gene.stable.ID"] %>% as.vector()
cluster3_genes = df_data[df_data$cluster == "Cluster 3: Down-regulated in both species", "Gene.stable.ID"] %>% as.vector()
cluster4_genes = df_data[df_data$cluster == "Cluster 4: Up-regulated in human reverse in mouse", "Gene.stable.ID"] %>% as.vector()
cluster5_genes = df_data[df_data$cluster == "Cluster 5: Not-significant in both species", "Gene.stable.ID"] %>% as.vector()

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x <= (qnt[1] - H)] <- NA
  y[x >= (qnt[2] + H)] <- NA
  y
}

####
elementFT <- c("Quies -> Enh","ReprPC -> Enh","Tss -> Enh","Quies -> Het","Enh -> Quies",
               "Het -> Quies", "ReprPC -> Quies","Quies -> ReprPC","ReprPC -> TssBiv") 

#####处理与画图
print_change_plot <- function(locifile_path,my_comparisons,fig_title){
  df_active_a = read.csv(locifile_path, sep="\t")
  gr_active_a <- makeGRangesFromDataFrame(df_active_a,keep.extra.columns = T)
  peakAnno_a <- annotatePeak(gr_active_a, tssRegion=c(-3000, 3000), TxDb=hg19)
  peakAnno_a <- as.GRanges(peakAnno_a)
  peakAnno_a <- data.frame("gene" = gsub("\\..*", "", peakAnno_a$geneId),
                           "chr" = as.vector(seqnames(peakAnno_a)),
                           "start" = start(peakAnno_a),"end" = end(peakAnno_a),
                           "From" = peakAnno_a$From, "To" = peakAnno_a$To)
  peakAnno_a$FT  <-stringr::str_c(peakAnno_a$From," -> ",peakAnno_a$To)
  
  ####内部分转变类型
  data <- purrr::map(elementFT, function(e){
    tmp_genset = peakAnno_a[peakAnno_a$FT==e,]
    n1_a = as.vector(table(tmp_genset$gene[tmp_genset$gene %in% cluster1_genes]))
    n2_a = as.vector(table(tmp_genset$gene[tmp_genset$gene %in% cluster2_genes]))
    n3_a = as.vector(table(tmp_genset$gene[tmp_genset$gene %in% cluster3_genes]))
    n4_a = as.vector(table(tmp_genset$gene[tmp_genset$gene %in% cluster4_genes]))
    n5_a = as.vector(table(tmp_genset$gene[tmp_genset$gene %in% cluster5_genes]))
    gene_num=c(n1_a, n2_a, n3_a, n4_a, n5_a)
    group = c(rep("cluster 1", length(n1_a)), rep("cluster 2", length(n2_a)), rep("cluster 3", length(n3_a)), rep("cluster 4", length(n4_a)), rep("cluster 5", length(n5_a)) )
    data = data.frame(region_num=gene_num, group=group)
    data <- data %>%
      group_by(group)%>%
      mutate(value = remove_outliers(region_num))
    data <- na.omit(data)
    data$FT = e
    return(data)
  })%>%dplyr::bind_rows()
  
  print(table(data$group,data$FT))
  # n1_a = as.vector(table(peakAnno_a$gene[peakAnno_a$gene %in% cluster1_genes]))
  # n2_a = as.vector(table(peakAnno_a$gene[peakAnno_a$gene %in% cluster2_genes]))
  # n3_a = as.vector(table(peakAnno_a$gene[peakAnno_a$gene %in% cluster3_genes]))
  # n4_a = as.vector(table(peakAnno_a$gene[peakAnno_a$gene %in% cluster4_genes]))
  # n5_a = as.vector(table(peakAnno_a$gene[peakAnno_a$gene %in% cluster5_genes]))
  # gene_num=c(n1_a, n2_a, n3_a, n4_a, n5_a)
  # group = c(rep("cluster 1", length(n1_a)), rep("cluster 2", length(n2_a)), rep("cluster 3", length(n3_a)), rep("cluster 4", length(n4_a)), rep("cluster 5", length(n5_a)) )
  # data = data.frame(region_num=gene_num, group=group)
  # data <- data %>%
  #   group_by(group)%>%
  #   mutate(value = remove_outliers(region_num))
  # data <- na.omit(data)
  p1 <- ggplot(data = data,aes(x = group,y = region_num,fill=group)) +
    geom_violin(adjust = .5) +geom_boxplot(width=0.1, fill="white")+
    facet_wrap(~ FT, ncol=3, scales = "free_y")+
    stat_compare_means(comparisons = my_comparisons,label = "p.signif", tip.length=0.02,method.args = list(alternative = "greater"), method = "t.test") +
    theme_classic(base_family = "sans",base_size = 18,base_line_size = 1.1)+
    scale_fill_manual(values=c("cluster 1" = "#CB9C7A", "cluster 2" = "#8696a7", "cluster 3" = "#CDB97D",
                               "cluster 4" = "#7b8b6f", "cluster 5" = "#A59B95"))+
    labs(title = fig_title, y="Chromatin state density in gene", x=NULL, fill=NULL) +
    scale_x_discrete(labels = c("cluster 1" = "Cluster1", "cluster 2" = "Cluster2", "cluster 3" = "Cluster3",
                                "cluster 4" = "Cluster4", "cluster 5" = "Cluster5"))+
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5,size=16),##0.5居中
          axis.title=element_text(size=18),
          axis.text=element_text(size=15),
          strip.background = element_blank())
  return(p1)
}

fig_title = "Conserved orthologous block distribution"
locifile_path = "/NAS/luozh/CRC_conservation/step10_chromatin_state/figure_data/all_conserved_state_change.txt"
my_comparisons <- list(c("cluster 1", "cluster 5"), c("cluster 2", "cluster 5"), c("cluster 3", "cluster 5"), c("cluster 4", "cluster 5"))
p = print_change_plot(locifile_path,my_comparisons,fig_title)
ggsave("/home/shimw/project/conservation/state_change/conserved_block_distribution.pdf",p,width=14,height = 11)

fig_title = "Human specific orthologous block distribution"
locifile_path = "/NAS/luozh/CRC_conservation/step10_chromatin_state/figure_data/human_specific_state_change.txt"
my_comparisons <- list(c("cluster 1", "cluster 5"), c("cluster 2", "cluster 5"), c("cluster 3", "cluster 5"), c("cluster 4", "cluster 5"))
p = print_change_plot(locifile_path,my_comparisons,fig_title)
ggsave("/home/shimw/project/conservation/state_change/human_specific_block_distribution.pdf",p,width=14,height = 11)

fig_title = "Mouse specific orthologous block distribution"
locifile_path = "/NAS/luozh/CRC_conservation/step10_chromatin_state/figure_data/mouse_specific_state_change.txt"
my_comparisons <- list(c("cluster 1", "cluster 5"), c("cluster 2", "cluster 5"), c("cluster 3", "cluster 5"), c("cluster 4", "cluster 5"))
p = print_change_plot(locifile_path,my_comparisons,fig_title)
ggsave("/home/shimw/project/conservation/state_change/mouse_specific_block_distribution.pdf",p,width=14,height = 11)


## Human track图
locifile_path = "/NAS/luozh/CRC_conservation/step10_chromatin_state/figure_data/human_specific_state_change.txt"
df_active_a = read.csv(locifile_path, sep="\t")
gr_active_a <- makeGRangesFromDataFrame(df_active_a,keep.extra.columns = T)
peakAnno_a <- annotatePeak(gr_active_a, tssRegion=c(-3000, 3000), TxDb=hg19)
peakAnno_a <- as.GRanges(peakAnno_a)
peakAnno_a <- data.frame("gene" = gsub("\\..*", "", peakAnno_a$geneId),
                         "chr" = as.vector(seqnames(peakAnno_a)),
                         "start" = start(peakAnno_a),"end" = end(peakAnno_a),
                         "From" = peakAnno_a$From, "To" = peakAnno_a$To)
peakAnno_a$FT  <-stringr::str_c(peakAnno_a$From," -> ",peakAnno_a$To)
track_df = peakAnno_a[peakAnno_a$FT=="Tss -> Enh",]
track_df = track_df %>%
  mutate(cluster = case_when(
    gene%in%cluster1_genes ~ "cluster1",
    gene%in%cluster2_genes ~ "cluster2",
    gene%in%cluster3_genes ~ "cluster3",
    gene%in%cluster4_genes ~ "cluster4",
    gene%in%cluster5_genes ~ "cluster5"
  ))%>%
  tidyr::drop_na()%>%
  select(gene,FT,cluster)%>%
  distinct()
readr::write_csv(track_df, "~/project/conservation/track/state_change_gene.csv")

# ENSG00000028137 TNFRSF1B
# ENSG00000130766 SESN2


## Mouse track图########
locifile_path = "/NAS/luozh/CRC_conservation/step10_chromatin_state/figure_data/mouse_specific_state_change.txt"
df_active_a = read.csv(locifile_path, sep="\t")
gr_active_a <- makeGRangesFromDataFrame(df_active_a,keep.extra.columns = T)
peakAnno_a <- annotatePeak(gr_active_a, tssRegion=c(-3000, 3000), TxDb=hg19)
peakAnno_a <- as.GRanges(peakAnno_a)
peakAnno_a <- data.frame("gene" = gsub("\\..*", "", peakAnno_a$geneId),
                         "chr" = as.vector(seqnames(peakAnno_a)),
                         "start" = start(peakAnno_a),"end" = end(peakAnno_a),
                         "From" = peakAnno_a$From, "To" = peakAnno_a$To)
peakAnno_a$FT  <-stringr::str_c(peakAnno_a$From," -> ",peakAnno_a$To)
track_df = peakAnno_a[peakAnno_a$FT=="Tss -> Enh",]
track_df = track_df %>%
  mutate(cluster = case_when(
    gene%in%cluster1_genes ~ "cluster1",
    gene%in%cluster2_genes ~ "cluster2",
    gene%in%cluster3_genes ~ "cluster3",
    gene%in%cluster4_genes ~ "cluster4",
    gene%in%cluster5_genes ~ "cluster5"
  ))%>%
  tidyr::drop_na()%>%
  select(gene,FT,cluster)%>%
  distinct()
## 转换成mouse gene id
mouse_df = left_join(track_df, df_data[,c("Gene.stable.ID", "Mouse.gene.stable.ID")],
          by = c("gene"= "Gene.stable.ID"))%>%
  select(Mouse.gene.stable.ID, FT, cluster)%>%
  rename(gene = Mouse.gene.stable.ID)

readr::write_csv(mouse_df, "~/project/conservation/track/mouse_state_change_gene.csv")



get_change_df <- function(locifile_path){
  df_active_a = read.csv(locifile_path, sep="\t")
  gr_active_a <- makeGRangesFromDataFrame(df_active_a,keep.extra.columns = T)
  peakAnno_a <- annotatePeak(gr_active_a, tssRegion=c(-3000, 3000), TxDb=hg19)
  peakAnno_a <- as.GRanges(peakAnno_a)
  peakAnno_a <- data.frame("gene" = gsub("\\..*", "", peakAnno_a$geneId),
                           "chr" = as.vector(seqnames(peakAnno_a)),
                           "start" = start(peakAnno_a),"end" = end(peakAnno_a),
                           "From" = peakAnno_a$From, "To" = peakAnno_a$To)
  peakAnno_a$FT  <-stringr::str_c(peakAnno_a$From," -> ",peakAnno_a$To)
  
  ####内部分转变类型
  data <- purrr::map(elementFT, function(e){
    tmp_genset = peakAnno_a[peakAnno_a$FT==e,]
    n1_a = as.vector(table(tmp_genset$gene[tmp_genset$gene %in% cluster1_genes]))
    n2_a = as.vector(table(tmp_genset$gene[tmp_genset$gene %in% cluster2_genes]))
    n3_a = as.vector(table(tmp_genset$gene[tmp_genset$gene %in% cluster3_genes]))
    n4_a = as.vector(table(tmp_genset$gene[tmp_genset$gene %in% cluster4_genes]))
    n5_a = as.vector(table(tmp_genset$gene[tmp_genset$gene %in% cluster5_genes]))
    gene_num=c(n1_a, n2_a, n3_a, n4_a, n5_a)
    group = c(rep("cluster 1", length(n1_a)), rep("cluster 2", length(n2_a)), rep("cluster 3", length(n3_a)), rep("cluster 4", length(n4_a)), rep("cluster 5", length(n5_a)) )
    data = data.frame(region_num=gene_num, group=group)
    data <- data %>%
      group_by(group)%>%
      mutate(value = remove_outliers(region_num))
    data <- na.omit(data)
    data$FT = e
    return(data)
  })%>%dplyr::bind_rows()
  return(data)
}

locifile_path = "/NAS/luozh/CRC_conservation/step10_chromatin_state/figure_data/all_conserved_state_change.txt"
all_conserved_df = get_change_df(locifile_path)
all_conserved_df$type = "Conserved"
locifile_path = "/NAS/luozh/CRC_conservation/step10_chromatin_state/figure_data/human_specific_state_change.txt"
human_specific_df = get_change_df(locifile_path)
human_specific_df$type = "Human specific"
locifile_path = "/NAS/luozh/CRC_conservation/step10_chromatin_state/figure_data/mouse_specific_state_change.txt"
mouse_specific_df = get_change_df(locifile_path)
mouse_specific_df$type = "Mouse specific"
all_df = dplyr::bind_rows(list(all_conserved_df,human_specific_df,mouse_specific_df))
my_comparisons <- list(c("cluster 1", "cluster 5"), c("cluster 2", "cluster 5"), c("cluster 3", "cluster 5"), c("cluster 4", "cluster 5"))
p1 <- ggplot(data = all_df[all_df$FT=="Tss -> Enh",],aes(x = group,y = region_num,fill=group)) +
  geom_violin(adjust = .5) +geom_boxplot(width=0.1, fill="white")+
  facet_wrap(~ type, ncol=3, scales = "free_y")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif", tip.length=0.02,method.args = list(alternative = "greater"), method = "t.test") +
  theme_classic(base_family = "sans",base_size = 18,base_line_size = 1.1)+
  scale_fill_manual(values=c("cluster 1" = "#CB9C7A", "cluster 2" = "#8696a7", "cluster 3" = "#CDB97D",
                             "cluster 4" = "#7b8b6f", "cluster 5" = "#A59B95"))+
  labs(title = "Tss -> Enh", y="Chromatin state density in gene", x=NULL, fill=NULL) +
  scale_x_discrete(labels = c("cluster 1" = "Cluster1", "cluster 2" = "Cluster2", "cluster 3" = "Cluster3",
                              "cluster 4" = "Cluster4", "cluster 5" = "Cluster5"))+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5,size=16),##0.5居中
        axis.title=element_text(size=18),
        axis.text=element_text(size=14),
        strip.background = element_blank())

ggsave("/home/shimw/project/conservation/state_change/Tss_Enh_block_distribution.pdf",p1,width=13.5,height = 4.5)



#######
###算相关性
df_data = read.csv("/home/panxl/CRC2/all_cluater")
get_peak_number <- function(locifile_path){
  df_active_a = read.csv(locifile_path, sep="\t")
  gr_active_a <- makeGRangesFromDataFrame(df_active_a,keep.extra.columns = T)
  peakAnno_a <- annotatePeak(gr_active_a, tssRegion=c(-3000, 3000), TxDb=hg19)
  peakAnno_a <- as.GRanges(peakAnno_a)
  peakAnno_a <- data.frame("gene" = gsub("\\..*", "", peakAnno_a$geneId),
                           "chr" = as.vector(seqnames(peakAnno_a)),
                           "start" = start(peakAnno_a),"end" = end(peakAnno_a),
                           "From" = peakAnno_a$From, "To" = peakAnno_a$To)
  peakAnno_a$FT  <-stringr::str_c(peakAnno_a$From," -> ",peakAnno_a$To)
  data <- purrr::map(elementFT, function(e){
    tmp_genset = peakAnno_a[peakAnno_a$FT==e,]
    tmp_genset <- tmp_genset%>%dplyr::group_by(gene)%>%summarise(peak_num = n())
    tmp_genset <- dplyr::left_join(tmp_genset, df_data,by = c("gene" = "Gene.stable.ID"))%>%tidyr::drop_na()
    # cor.test(tmp_genset$peak_num,tmp_genset$Human_logFC1)
    tmp_genset$FT = e
    # tmp_genset <- tmp_genset%>%
    #   mutate(value = remove_outliers(peak_num))
    return(tmp_genset)
  })%>%dplyr::bind_rows()
  data <- tidyr::drop_na(data)
  return(data)
}


c("Quies -> Enh","ReprPC -> Enh","Tss -> Enh",
  "Quies -> Het","Enh -> Quies","Het -> Quies",
  "ReprPC -> Quies","Quies -> ReprPC","ReprPC -> TssBiv")
######conserved
all_trans = c("Quies -> Enh","ReprPC -> Enh","Tss -> Enh",
              "Quies -> Het","Enh -> Quies","Het -> Quies",
              "ReprPC -> Quies","Quies -> ReprPC","ReprPC -> TssBiv")
sub_trans = c("Quies -> Enh","Tss -> Enh",
              "Enh -> Quies")
######这里没用，下面用了


locifile_path = "/NAS/luozh/CRC_conservation/step10_chromatin_state/figure_data/all_conserved_state_change.txt"
conserved_df <- get_peak_number(locifile_path)
conserved_df = conserved_df[conserved_df$FT%in% c("Quies -> Enh","Tss -> Enh",
                                                  "Enh -> Quies"),]
conserved_df$FT = factor(conserved_df$FT,
                         levels = c("Quies -> Enh","Tss -> Enh",
                                    "Enh -> Quies"))
conserved_df<-conserved_df%>%tidyr::drop_na()
mm = purrr::map(unique(conserved_df$FT),function(e){
  tmp_df = conserved_df[conserved_df$FT==e,]
  kk = cor.test(tmp_df$peak_num,tmp_df$Human_logFC1)
  print(e)
  print(paste("r =",round(kk$estimate,2)))
  print(paste("p-value =",kk$p.value))
})

conserved_df = conserved_df%>%group_by(FT)%>%mutate(cols = densCols(log2(peak_num),
                                                     Human_logFC1, 
                                                     colramp=colorRampPalette(c("#809fff","green","yellow", "red")), nbin = 100))



f_labels <- data.frame(FT =c("Quies -> Enh","Tss -> Enh",
                               "Enh -> Quies"),
                       label = c("r = 0.04, p-value = 0.07",
                                 # "r = 0.16, p-value = 0.12",
                                 "r = -0.23, p-value = 2.3e-17",
                                 # "r = -0.14, p-value = 0.05",
                                 "r = -0.23, p-value = 1.2e-10"
                                 # "r = 0.14, p-value = 0.59",
                                 # "r = 0.01, p-value = 0.89",
                                 # "r = -0.11, p-value = 0.004",
                                 # "r = -0.04, p-value = 0.42"
                                 ),
                       xa = c(4,3,4),
                       ya = c(7,3,3))
f_labels$FT = factor(f_labels$FT,
                       levels =c("Quies -> Enh","Tss -> Enh",
                                 "Enh -> Quies"))
p1 <- ggplot(data = conserved_df,aes(x=log2(peak_num),y=Human_logFC1)) +
  geom_point(size=0.1,color = conserved_df$cols) +
  # geom_smooth(method=lm,se = F,color="red")+
  geom_text( aes(x = xa,y = ya,label = label), data = f_labels)+
  theme_classic(base_family = "sans",base_size = 18,base_line_size = 1)+
  scale_y_continuous(expand = expansion(mult=c(0,0.1)))+
  facet_wrap(~ FT, ncol=3, scales = "free")+
  labs(y="Gene expression fold change",x="log2(the density of orthologous block)",title="Conserved state change")+
  theme_classic(base_line_size = 1) +
  theme(axis.text = element_text(size=12, color="black"),
        title = element_text(size = 14,
                             color = "black",
                             hjust = 0.5),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.text.y = element_text(size = 12,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.ticks.length.x = unit(1,"mm"),
        strip.background = element_blank(),
        plot.margin = margin(t = 0.1, r = 0.6,b = 0.1,l = 0.1 ,unit = "cm"),
        panel.spacing = unit(1, "lines"),
        strip.text = element_text(size = 12))

ggsave("/home/shimw/project/conservation/state_density/conserved.pdf",p1,width = 10,height = 3.6)


########################################################################################################
#####human specific 
locifile_path = "/NAS/luozh/CRC_conservation/step10_chromatin_state/figure_data/human_specific_state_change.txt"
human_df <- get_peak_number(locifile_path)
tmp_trans = sub_trans
human_df = human_df[human_df$FT%in% tmp_trans,]
human_df$FT = factor(human_df$FT,
                         levels = tmp_trans)

mm = purrr::map(unique(human_df$FT),function(e){ 
  tmp_df = human_df[human_df$FT==e,]
  kk = cor.test(tmp_df$peak_num,tmp_df$Human_logFC1)
  print(as.vector(e))
  print(paste("r =",round(kk$estimate,2)))
  print(paste("p-value =",kk$p.value))
})


human_df = human_df%>%group_by(FT)%>%mutate(cols = densCols(log2(peak_num),
                                            Human_logFC1, 
                                            colramp=colorRampPalette(c("#809fff","green","yellow", "red")), nbin = 100))

f_labels <- data.frame(
  FT =tmp_trans,
  label = c("r = 0.04, p-value = 0.02",
            # "r = 0.06, p-value = 0.19",
            "r = -0.25, p-value = 2.1e-53",
            # "r = -0.09, p-value = 0.0005"
            "r = -0.11, p-value = 5.0e-09"
            # "r = -0.1, p-value = 0.40",
            # "r = -0.01, p-value =  0.71",
            # "r = -0.11, p-value = 4.6e-05",
            # "r = -0.05, p-value = 0.048"
             ),
  xa = c(4,2.5,4),
  ya = c(6,5,5))
f_labels$FT = factor(f_labels$FT,
                     levels =tmp_trans)
p2 <- ggplot(data = human_df,aes(x=log2(peak_num),y=Human_logFC1)) +
  geom_point(size=0.1,color = human_df$cols) +
  # geom_smooth(method=lm,se = F,color="red")+
  geom_text( aes(x = xa,y = ya,label = label), data = f_labels)+
  theme_classic(base_family = "sans",base_size = 18,base_line_size = 1)+
  scale_y_continuous(expand = expansion(mult=c(0,0.1)))+
  facet_wrap(~ FT, ncol=3, scales = "free")+
  labs(y="Gene expression fold change",x="log2(the density of orthologous block)",title="Human specific state change")+
  theme_classic(base_line_size = 1) +
  theme(axis.text = element_text(size=12, color="black"),
        title = element_text(size = 14,
                             color = "black",
                             hjust = 0.5),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.text.y = element_text(size = 12,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.ticks.length.x = unit(1,"mm"),
        strip.background = element_blank(),
        plot.margin = margin(t = 0.1, r = 0.6,b = 0.1,l = 0.1 ,unit = "cm"),
        panel.spacing = unit(1, "lines"),
        strip.text = element_text(size = 12))
ggsave("/home/shimw/project/conservation/state_density/human_specific.pdf",p2,width = 10,height = 3.6)

########################################################################################################
#####mouse specific 
locifile_path = "/NAS/luozh/CRC_conservation/step10_chromatin_state/figure_data/mouse_specific_state_change.txt"


# tmp_trans = all_trans
tmp_trans = sub_trans
mouse_df <- get_peak_number(locifile_path)
mouse_df = mouse_df[mouse_df$FT%in% tmp_trans,]
mouse_df$FT = factor(mouse_df$FT,
                     levels = tmp_trans)

mm = purrr::map(unique(mouse_df$FT),function(e){
  tmp_df = mouse_df[mouse_df$FT==e,]
  kk = cor.test(tmp_df$peak_num,tmp_df$Mouse_logFC1)
  print(as.vector(e))
  print(paste("r =",round(kk$estimate,2)))
  print(paste("p-value =",kk$p.value))
})


mouse_df = mouse_df%>%group_by(FT)%>%mutate(cols = densCols(log2(peak_num),
                                                            Mouse_logFC1, 
                                                            colramp=colorRampPalette(c("#809fff","green","yellow", "red")), nbin = 100))

f_labels <- data.frame(
  FT =tmp_trans,
  label = c("r = 0.1, p-value = 4.1e-13",###
            # "r = 0.31, p-value = 1.4e-06",
            "r = -0.08, p-value = 0.002",###
            # "r = -0.2, p-value = 3.0e-13",
            "r = -0.1, p-value = 9.8e-05"###
            # "r = 0.21, p-value = 0.22",
            # "r = 0.02, p-value =  0.58",
            # "r = -0.1, p-value = 7.1e-05",
            # "r = 0.11, p-value = 0.02"
  ),
  xa = c(4.5,2.5,4),
  ya = c(7,4.5,4.5))
f_labels$FT = factor(f_labels$FT,
                     levels =tmp_trans)
p3 <- ggplot(data = mouse_df,aes(x=log2(peak_num),y=Mouse_logFC1)) +
  geom_point(size=0.1,color = mouse_df$cols) +
  # geom_smooth(method=lm,se = F,color="red")+
  geom_text( aes(x = xa,y = ya,label = label), data = f_labels)+
  theme_classic(base_family = "sans",base_size = 18,base_line_size = 1)+
  scale_y_continuous(expand = expansion(mult=c(0,0.1)))+
  facet_wrap(~ FT, ncol=3, scales = "free")+
  labs(y="Gene expression fold change",x="log2(the density of orthologous block)",title="Mouse specific state change")+
  theme_classic(base_line_size = 1) +
  theme(axis.text = element_text(size=12, color="black"),
        title = element_text(size = 14,
                             color = "black",
                             hjust = 0.5),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.text.y = element_text(size = 12,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.ticks.length.x = unit(1,"mm"),
        strip.background = element_blank(),
        plot.margin = margin(t = 0.1, r = 0.6,b = 0.1,l = 0.1 ,unit = "cm"),
        panel.spacing = unit(1, "lines"),
        strip.text = element_text(size = 12))

ggsave("/home/shimw/project/conservation/state_density/mouse_specific.pdf",p3,width = 10,height = 3.6)


###############################
###最终合并
conserved_df <- conserved_df[c("gene","peak_num","Human_logFC1","cluster","FT","cols")]
names(conserved_df)[3] = "logFC1"
conserved_df$group = "Conserved"
human_df <- human_df[c("gene","peak_num","Human_logFC1","cluster","FT","cols")]
names(human_df)[3] = "logFC1"
human_df$group = "Human specific"
mouse_df <-  mouse_df[c("gene","peak_num","Mouse_logFC1","cluster","FT","cols")]
names(mouse_df)[3] = "logFC1"
mouse_df$group = "Mouse specific"

all_df = dplyr::bind_rows(list(conserved_df,human_df,mouse_df))

q_labels=data.frame(
  group =factor(c("Conserved","Human specific","Mouse specific"),levels = c("Conserved","Human specific","Mouse specific")),
  label = c("r = 0.04, p-value = 0.07",
            "r = 0.04, p-value = 0.02",
            "r = 0.1, p-value = 4.1e-13"
  ),
  xa = c(4,4,4.5),
  ya = c(7,6,7))

p4 <- ggplot(data = all_df[all_df$FT == "Quies -> Enh",],aes(x=log2(peak_num),y=logFC1)) +
  geom_point(size=0.1,color = all_df[all_df$FT == "Quies -> Enh",]$cols) +
  # geom_smooth(method=lm,se = F,color="red")+
  geom_text( aes(x = xa,y = ya,label = label), data = q_labels)+
  theme_classic(base_family = "sans",base_size = 18,base_line_size = 1)+
  scale_y_continuous(expand = expansion(mult=c(0,0.1)))+
  facet_wrap( ~ group, ncol = 3, scales = "free")+
  labs(y="Gene expression fold change",x=NULL,title="Quies -> Enh")+
  theme_classic(base_line_size = 1) +
  theme(axis.text = element_text(size=12, color="black"),
        title = element_text(size = 14,
                             color = "black",
                             hjust = 0.5),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.text.y = element_text(size = 12,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.ticks.length.x = unit(1,"mm"),
        strip.background = element_blank(),
        plot.margin = margin(t = 0.1, r = 0.6,b = 0.1,l = 0.1 ,unit = "cm"),
        # panel.spacing = unit(1, "lines"),
        strip.text = element_text(size = 12))

e_labels=data.frame(
  group =factor(c("Conserved","Human specific","Mouse specific"),levels = c("Conserved","Human specific","Mouse specific")),
  label = c("r = -0.23, p-value = 1.2e-10",
            "r = -0.11, p-value = 5.0e-09",
            "r = -0.1, p-value = 9.8e-05"###
  ),
  xa = c(4,5,4.5),
  ya = c(4.5,6,5))
p5 = ggplot(data = all_df[all_df$FT == "Enh -> Quies",],aes(x=log2(peak_num),y=logFC1)) +
  geom_point(size=0.1,color = all_df[all_df$FT == "Enh -> Quies",]$cols) +
  # geom_smooth(method=lm,se = F,color="red")+
  geom_text( aes(x = xa,y = ya,label = label), data = e_labels)+
  theme_classic(base_family = "sans",base_size = 18,base_line_size = 1)+
  scale_y_continuous(expand = expansion(mult=c(0,0.1)))+
  facet_wrap( ~ group, ncol = 3, scales = "free")+
  labs(y="Gene expression fold change",x="log2(the density of orthologous blocks)",title="Enh -> Quies")+
  theme_classic(base_line_size = 1) +
  theme(axis.text = element_text(size=12, color="black"),
        title = element_text(size = 14,
                             color = "black",
                             hjust = 0.5),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.text.y = element_text(size = 12,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.ticks.length.x = unit(1,"mm"),
        strip.background = element_blank(),
        plot.margin = margin(t = 0.1, r = 0.6,b = 0.1,l = 0.1 ,unit = "cm"),
        # panel.spacing = unit(1, "lines"),
        strip.text = element_text(size = 12))

library(ggpubr)

p7 = ggarrange(p4,p5,ncol=1,nrow=2)
ggsave("/home/shimw/project/conservation/state_density/orther.pdf",p7,width = 10,height = 7)




t_labels=data.frame(
  group =factor(c("Conserved","Human specific","Mouse specific"),levels = c("Conserved","Human specific","Mouse specific")),
  label = c("r = -0.23, p-value = 2.3e-17",
            "r = -0.25, p-value = 2.1e-53",
            "r = -0.08, p-value = 0.002"###
  ),
  xa = c(3,3,2.8),
  ya = c(4.5,6,5))
p6 = ggplot(data = all_df[all_df$FT == "Tss -> Enh",],aes(x=log2(peak_num),y=logFC1)) +
  geom_point(size=0.1,color = all_df[all_df$FT == "Tss -> Enh",]$cols) +
  # geom_smooth(method=lm,se = F,color="red")+
  geom_text( aes(x = xa,y = ya,label = label), data = t_labels)+
  theme_classic(base_family = "sans",base_size = 18,base_line_size = 1)+
  scale_y_continuous(expand = expansion(mult=c(0,0.1)))+
  facet_wrap( ~ group, ncol = 3, scales = "free")+
  labs(y="Gene expression fold change",x="log2(the density of orthologous blocks)",title="Tss -> Enh")+
  theme_classic(base_line_size = 1) +
  theme(axis.text = element_text(size=12, color="black"),
        title = element_text(size = 14,
                             color = "black",
                             hjust = 0.5),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.text.y = element_text(size = 12,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.ticks.length.x = unit(1,"mm"),
        strip.background = element_blank(),
        plot.margin = margin(t = 0.1, r = 0.6,b = 0.1,l = 0.1 ,unit = "cm"),
        # panel.spacing = unit(1, "lines"),
        strip.text = element_text(size = 12))
ggsave("/home/shimw/project/conservation/state_density/Tss_Enh.pdf",p6,width = 10,height = 3.6)




