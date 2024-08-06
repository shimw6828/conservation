library(ggpubr)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
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



df_active = read.csv("/NAS/luozh/CRC_conservation/step10_chromatin_state/state_compare/ortholog_state/mouse_quiescent_ActiveEnhancer_loci.txt", sep="\t")
#1. quiscent_active_enhancer_loci.txt
#2. promoter_active_enhancer_loci.txt
#3. human_quiescent_ActiveEnhancer_loci.txt
#4. mouse_quiescent_ActiveEnhancer_loci.txt



print_change_plot <- function(locifile_path,my_comparisons,fig_title){
  df_active = read.csv(locifile_path, sep="\t")
  gr_active=GRanges(seqnames=df_active$chr,
                    ranges=IRanges(start=df_active$start,
                                   end=df_active$end),
                    strand = rep("+", nrow(df_active)),
                    p_value = df_active$id
  )
  peakAnno_active <- annotatePeak(gr_active, tssRegion=c(-3000, 3000),
                                  TxDb=txdb, annoDb="org.Hs.eg.db")
  data_active = peakAnno_active@anno
  n1_a = as.vector(table(data_active$ENSEMBL[data_active$ENSEMBL %in% cluster1_genes]))
  n2_a = as.vector(table(data_active$ENSEMBL[data_active$ENSEMBL %in% cluster2_genes]))
  n3_a = as.vector(table(data_active$ENSEMBL[data_active$ENSEMBL %in% cluster3_genes]))
  n4_a = as.vector(table(data_active$ENSEMBL[data_active$ENSEMBL %in% cluster4_genes]))
  n5_a = as.vector(table(data_active$ENSEMBL[data_active$ENSEMBL %in% cluster5_genes]))
  gene_num=c(n1_a, n2_a, n3_a, n4_a, n5_a)
  group = c(rep("cluster 1", length(n1_a)), rep("cluster 2", length(n2_a)), rep("cluster 3", length(n3_a)), rep("cluster 4", length(n4_a)), rep("cluster 5", length(n5_a)) )
  data = data.frame(region_num=gene_num, group=group)
  data <- data %>%
    group_by(group) %>%
    mutate(value = remove_outliers(region_num))
  data <- na.omit(data)
  p1 <- ggplot(data = data,aes(x = group,y = region_num,fill=group)) +
    geom_violin(adjust = .5) +geom_boxplot(width=0.1, fill="white")+
    stat_compare_means(comparisons = my_comparisons,label = "p.signif", tip.length=0.02,method.args = list(alternative = "greater"), method = "t.test") +
    theme_classic(base_family = "sans",base_size = 18,base_line_size = 1.1)+
    scale_fill_manual(values=c("cluster 1" = "#CB9C7A", "cluster 2" = "#8696a7", "cluster 3" = "#CDB97D",
                               "cluster 4" = "#7b8b6f", "cluster 5" = "#A59B95"))+
    labs(title = fig_title, y="Chromatin state density in gene", x=NULL, fill=NULL) +
    scale_x_discrete(labels = c("cluster 1" = "Cluster1", "cluster 2" = "Cluster2", "cluster 3" = "Cluster3",
                                "cluster 4" = "Cluster4", "cluster 5" = "Cluster5"))+
    theme(legend.position = "none",plot.title = element_text(hjust = 0,size=16),
          axis.title=element_text(size=18),
          axis.text=element_text(size=15))
  return(p1)
  
}

#1. quiscent_active_enhancer_loci.txt
fig_title = "Quiescent to active enhancer in human and mouse"
locifile_path = "/NAS/luozh/CRC_conservation/step10_chromatin_state/state_compare/ortholog_state/quiscent_active_enhancer_loci.txt"
my_comparisons <- list(c("cluster 1", "cluster 5"), c("cluster 2", "cluster 5"), c("cluster 3", "cluster 5"), c("cluster 4", "cluster 5"))
p = print_change_plot(locifile_path,my_comparisons,fig_title)
ggsave("/home/shimw/project/conservation/state_change/quiscent_active_enhancer_both.pdf",p,width=6,height = 5)


#2. promoter_active_enhancer_loci.txt

fig_title = "Promoter to active enhancer in human and mouse"
locifile_path = "/NAS/luozh/CRC_conservation/step10_chromatin_state/state_compare/ortholog_state/promoter_active_enhancer_loci.txt"
my_comparisons <- list(c("cluster 1", "cluster 5"), c("cluster 2", "cluster 5"), c("cluster 3", "cluster 5"), c("cluster 4", "cluster 5"))
p = print_change_plot(locifile_path,my_comparisons,fig_title)
ggsave("/home/shimw/project/conservation/state_change/promoter_active_enhancer_both.pdf",p,width=6,height = 5)


#3. human_quiescent_ActiveEnhancer_loci.txt
fig_title = "Quiescent to active enhancer only in human"
locifile_path = "/NAS/luozh/CRC_conservation/step10_chromatin_state/state_compare/ortholog_state/human_quiescent_ActiveEnhancer_loci.txt"
my_comparisons <- list(c("cluster 1", "cluster 5"), c("cluster 2", "cluster 5"), c("cluster 3", "cluster 5"), c("cluster 4", "cluster 5"))
p = print_change_plot(locifile_path,my_comparisons,fig_title)
ggsave("/home/shimw/project/conservation/state_change/quiscent_active_enhancer_human.pdf",p,width=6,height = 5)


#4. mouse_quiescent_ActiveEnhancer_loci.txt
fig_title = "Quiescent to active enhancer only in mouse"
locifile_path = "/NAS/luozh/CRC_conservation/step10_chromatin_state/state_compare/ortholog_state/mouse_quiescent_ActiveEnhancer_loci.txt"
my_comparisons <- list(c("cluster 1", "cluster 5"), c("cluster 2", "cluster 5"), c("cluster 3", "cluster 5"), c("cluster 4", "cluster 5"))
p = print_change_plot(locifile_path,my_comparisons,fig_title)
ggsave("/home/shimw/project/conservation/state_change/quiscent_active_enhancer_mouse.pdf",p,width=6,height = 5)


#5. mouse_promoter_ActiveEnhancer_loci.txt
fig_title = "Promoter to active enhancer only in mouse"
locifile_path = "/NAS/luozh/CRC_conservation/step10_chromatin_state/state_compare/ortholog_state/mouse_promoter_ActiveEnhancer_loci.txt"
my_comparisons <- list(c("cluster 1", "cluster 5"), c("cluster 2", "cluster 5"), c("cluster 3", "cluster 5"), c("cluster 4", "cluster 5"))
p = print_change_plot(locifile_path,my_comparisons,fig_title)
ggsave("/home/shimw/project/conservation/state_change/promoter_active_enhancer_mouse.pdf",p,width=6,height = 5)


#6. human_promoter_ActiveEnhancer_loci.txt
fig_title = "Promoter to active enhancer only in human"
locifile_path = "/NAS/luozh/CRC_conservation/step10_chromatin_state/state_compare/ortholog_state/human_promoter_ActiveEnhancer_loci.txt"
my_comparisons <- list(c("cluster 1", "cluster 5"), c("cluster 2", "cluster 5"), c("cluster 3", "cluster 5"), c("cluster 4", "cluster 5"))
p = print_change_plot(locifile_path,my_comparisons,fig_title)
ggsave("/home/shimw/project/conservation/state_change/promoter_active_enhancer_human.pdf",p,width=6,height = 5)

