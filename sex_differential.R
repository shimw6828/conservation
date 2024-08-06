###sex differential analysis
####normal的性别差异
###########################先试试ERE表格的基因
ERE_gene <- readr::read_tsv("/home/shimw/project/conservation/Sex/sex_ER_genes.txt",col_names = c("human","mouse"))%>%
  tidyr::drop_na()

library("biomaRt")
ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "grch37.ensembl.org")
EREtoid = getBM(attributes=c('ensembl_gene_id', 'external_gene_name'), 
                   filters = 'external_gene_name', 
                   values = ERE_gene$human, 
                   mart = ensembl)

df_data = read.csv("/home/panxl/CRC2/all_cluater")
table(df_data$Mouse.gene.stable.ID%in%kk$Ensembl.Gene.ID)
df_data$cluster = stringr::str_sub(df_data$cluster,0,9)
ere_df <- df_data[df_data$Gene.stable.ID%in%EREtoid$ensembl_gene_id,]

table(df_data[df_data$Gene.stable.ID%in%EREtoid$ensembl_gene_id,]$cluster)

cgroups =unique(df_data$cluster)
#cgroups ="Cluster 1"
odds_ere_ratios <- purrr::map(cgroups,function(cgroups){
  genes_fisher_cluster = matrix(c(nrow(df_data[df_data$cluster==cgroups&df_data$Gene.stable.ID%in%EREtoid$ensembl_gene_id,]),
                                  nrow(df_data[df_data$cluster==cgroups&!(df_data$Gene.stable.ID%in%EREtoid$ensembl_gene_id),]),
                                  nrow(df_data[df_data$cluster!=cgroups&df_data$Gene.stable.ID%in%EREtoid$ensembl_gene_id,]),
                                  nrow(df_data[df_data$cluster!=cgroups&!(df_data$Gene.stable.ID%in%EREtoid$ensembl_gene_id),])),
                                nrow = 2,
                                dimnames =list(c("housekeeping", "non_housekeeping"),
                                               c("cluster1", "Not cluster1")))
  
  
  OR_1 = fisher.test(genes_fisher_cluster, conf.level = 0.95)
  #return(tibble("boxLabels"=stringr::str_replace(cgroups," ",""),"boxOdds"=OR_1$estimate,"boxCILow"=OR_1$conf.int[1],"boxCIHigh"=OR_1$conf.int[2]))
  return(tibble("boxLabels"=cgroups,"boxOdds"=OR_1$estimate,"boxCILow"=OR_1$conf.int[1],"boxCIHigh"=OR_1$conf.int[2]))
})%>%dplyr::bind_rows()



AR_gene <- readr::read_csv("/home/shimw/project/conservation/Sex/AR.csv")%>%
  tidyr::drop_na()
ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "grch37.ensembl.org")
ARtoid = getBM(attributes=c('ensembl_gene_id', 'external_gene_name'), 
                   filters = 'external_gene_name', 
                   values = AR_gene$Gene, 
                   mart = ensembl)

# df_data = read.csv("/home/panxl/CRC2/all_cluater")
table(df_data$Mouse.gene.stable.ID%in%kk$Ensembl.Gene.ID)
# df_data$cluster = stringr::str_sub(df_data$cluster,0,9)
ar_df <- df_data[df_data$Gene.stable.ID%in%ARtoid$ensembl_gene_id,]

table(df_data[df_data$Gene.stable.ID%in%ARtoid$ensembl_gene_id,]$cluster)

cgroups =unique(df_data$cluster)
#cgroups ="Cluster 1"
odds_ar_ratios <- purrr::map(cgroups,function(cgroups){
  genes_fisher_cluster = matrix(c(nrow(df_data[df_data$cluster==cgroups&df_data$Gene.stable.ID%in%ARtoid$ensembl_gene_id,]),
                                  nrow(df_data[df_data$cluster==cgroups&!(df_data$Gene.stable.ID%in%ARtoid$ensembl_gene_id),]),
                                  nrow(df_data[df_data$cluster!=cgroups&df_data$Gene.stable.ID%in%ARtoid$ensembl_gene_id,]),
                                  nrow(df_data[df_data$cluster!=cgroups&!(df_data$Gene.stable.ID%in%ARtoid$ensembl_gene_id),])),
                                nrow = 2,
                                dimnames =list(c("housekeeping", "non_housekeeping"),
                                               c("cluster1", "Not cluster1")))
  
  
  OR_1 = fisher.test(genes_fisher_cluster, conf.level = 0.95)
  #return(tibble("boxLabels"=stringr::str_replace(cgroups," ",""),"boxOdds"=OR_1$estimate,"boxCILow"=OR_1$conf.int[1],"boxCIHigh"=OR_1$conf.int[2]))
  return(tibble("boxLabels"=cgroups,"boxOdds"=OR_1$estimate,"boxCILow"=OR_1$conf.int[1],"boxCIHigh"=OR_1$conf.int[2]))
})%>%dplyr::bind_rows()


odds_ere_ratios$type = "Estrogen-responsive gene"
odds_ar_ratios$type = "Androgen-responsive gene"
scaleFUN <- function(x) sprintf("%.0f", x)
all_odds = dplyr::bind_rows(list(odds_ere_ratios,odds_ar_ratios))
p1 = ggplot(all_odds, aes(x = boxOdds, y = boxLabels,color=boxLabels)) + 
  geom_vline(aes(xintercept = 1.0), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = 1, height = 
                   .15 ) +
  geom_point(size = 2) +
  facet_wrap(~ type, ncol=3,scales = 'free_x') +
  labs(title = NULL, y=NULL, x="Odds ratio", fill=NULL) +
  scale_x_continuous(labels=scaleFUN)+
  scale_color_manual(values=c("Cluster 1" = "#CB9C7A", "Cluster 2" = "#8696a7", "Cluster 3" = "#CDB97D",
                              "Cluster 4" = "#7b8b6f", "Cluster 5" = "#A59B95"))+
  theme_bw(base_family = "sans",base_size = 16,base_line_size = 0.5)+
  theme(axis.text.y = element_text(size = 14,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        strip.text.x = element_text(size = 13),
        panel.spacing.x = unit(2, "mm"),
        legend.position = "none"
  )



#########################



#####gsva score
COAD_count <- read.csv("/home/shimw/project/conservation/COAD_hg19_count.csv",row.names = 1,check.names = F)
row.names(COAD_count) <- stringr::str_split_fixed(row.names(COAD_count),pattern = "\\.",2)[,1]
metadata = readr::read_csv("/home/shimw/project/conservation/COAD_hg19_metadata.csv")
metadata$Group <- ifelse(metadata$sample_type=="Primary Tumor","tumor","native")
clin <- GDCquery_clinic("TCGA-COAD", type = "clinical", save.csv = F)

table(stringr::str_sub(metadata$sample,1,12)%in%clin$submitter_id)
clin = clin[match(stringr::str_sub(metadata$sample,1,12),clin$submitter_id),]
metadata$Sex = clin$gender
metadata <- tidyr::drop_na(metadata)
COAD_count <- COAD_count[,metadata$sample]





# wuda_sex = readr::read_csv("/home/shimw/project/conservation/Sex/wuda_sex.csv")
# 
# stringr::str_c("A",wuda_sex$Sample)
# wuda_sex = data.frame("sample" = c(stringr::str_c("A",wuda_sex$Sample),stringr::str_c("B",wuda_sex$Sample)),
#                       "Sex" = c(wuda_sex$Gender,wuda_sex$Gender))
# wuda_sex = wuda_sex[match(human_info[human_info$DataSet=="ourdata",]$sample,wuda_sex$sample),]
# human_info[human_info$DataSet=="ourdata",]$Sex = tolower(wuda_sex$Sex)
# human_info <- tidyr::drop_na(human_info)
# 
# human_data_tumor_male <- human_data[,human_info[human_info$Group=="tumor"&human_info$Sex=="male",]$sample]
# human_data_tumor_female <- human_data[,human_info[human_info$Group=="tumor"&human_info$Sex=="female",]$sample]
# human_data_normal_male <- human_data[,human_info[human_info$Group=="native"&human_info$Sex=="male",]$sample]
# human_data_normal_female <- human_data[,human_info[human_info$Group=="native"&human_info$Sex=="female",]$sample]

df_intestine = read.csv("/home/shimw/project/enhancer_map/conserved/colon-specific.csv")
df_intestine_human = df_data[df_data$Mouse.gene.stable.ID%in%df_intestine$Gene,]
ar_df
ere_df


gene_list <- list("intestine"=df_intestine_human$Gene.stable.ID,
                  "ar" = ar_df$Gene.stable.ID,
                  "ere" = ere_df$Gene.stable.ID)


library(GSVA)
# g.tumor.male <- gsva(as.matrix(human_data_tumor_male), gene_list, verbose=FALSE, parallel.sz=1)
# g.tumor.female <- gsva(as.matrix(human_data_tumor_female), gene_list, verbose=FALSE, parallel.sz=1)
# g.normal.male <- gsva(as.matrix(human_data_normal_male), gene_list, verbose=FALSE, parallel.sz=1)
# g.normal.female <- gsva(as.matrix(human_data_normal_female), gene_list, verbose=FALSE, parallel.sz=1)
# 
# wilcox.test(g.tumor.male,g.tumor.female)
# t.test(g.tumor.male,g.tumor.female)
# wilcox.test(g.normal.male,g.normal.female)
# t.test(g.normal.male,g.normal.female)

COAD_count[,metadata[metadata$Sex=="male"&metadata$Group=="tumor",]$sample]
COAD_count[,metadata[metadata$Sex=="female"&metadata$Group=="tumor",]$sample]

COAD_count[,metadata[metadata$Sex=="male"&metadata$Group=="native",]$sample]
COAD_count[,metadata[metadata$Sex=="female"&metadata$Group=="native",]$sample]
g.tumor.male <- gsva(as.matrix(COAD_count[,metadata[metadata$Sex=="male"&metadata$Group=="tumor",]$sample]), gene_list, verbose=FALSE, parallel.sz=1,kcdf = "Poisson")
g.tumor.female <- gsva(as.matrix(COAD_count[,metadata[metadata$Sex=="female"&metadata$Group=="tumor",]$sample]), gene_list, verbose=FALSE, parallel.sz=1,kcdf = "Poisson")
g.normal.male <- gsva(as.matrix(COAD_count[,metadata[metadata$Sex=="male"&metadata$Group=="native",]$sample]), gene_list, verbose=FALSE, parallel.sz=1,kcdf = "Poisson")
g.normal.female <- gsva(as.matrix(COAD_count[,metadata[metadata$Sex=="female"&metadata$Group=="native",]$sample]), gene_list, verbose=FALSE, parallel.sz=1,kcdf = "Poisson")
wilcox.test(g.tumor.male,g.tumor.female)
t.test(g.tumor.male,g.tumor.female)
wilcox.test(g.normal.male,g.normal.female)
t.test(g.normal.male,g.normal.female)


###生存分析
##
survival_info <- readr::read_tsv("/home/shimw/project/Deconvolution/Survival_SupplementalTable_S1_20171025_xena_sp")
survival_info <- tibble::column_to_rownames(survival_info,var = "sample")

#table(metadata[metadata$Group=="tumor",]$sample%in%survival_info$sample)  提前检查过，都在里面
plot_surv <- function(gsetn,figtitle){
  g.tumor <- data.frame("sample" = c(colnames(g.tumor.male),colnames(g.tumor.female)),
                        "score" = c(g.tumor.male[gsetn,],g.tumor.female[gsetn,]))
  tmp_survival_info <- survival_info[g.tumor$sample,c("gender","OS","OS.time")]
  tmp_survival_info$score <- g.tumor$score
  data = tmp_survival_info %>%
    survminer::surv_cutpoint(
      time = "OS.time", event = "OS",
      variables = c("score"),
      minprop = 0.25, progressbar = TRUE
    ) %>%
    survminer::surv_categorize(labels = c("Low", "High")) %>%
    data.frame()
  
  tmp_survival_info$type = data$score
  tmp_survival_info$type <- factor(tmp_survival_info$type, levels = c("Low","High"))
  fit <- survfit(Surv(OS.time, OS) ~ type, data = tmp_survival_info)
  
  library(survminer)
  pp = ggsurvplot(fit, data = tmp_survival_info, pval = TRUE, pval.method = TRUE,palette = "aaas",
                  size = 1.2,
                  font.legend = c(14, "black"),
                  font.x = c(14, "bold", "black"),
                  font.y = c(14, "bold", "black"),
                  font.tickslab = c(12, "bold", "black"),
                  xlab = "Duration overall survival (days)",
                  ncensor.plot = FALSE, # plot the number of censored subjects at time t
                  surv.plot.height = 0.7,
                  risk.table.height = 0.15,
                  ncensor.plot.height = 0.15,
                  legend.labs = c("Low","High"),
                  ggtheme = ggplot2::theme_classic(),
                  legend.title=figtitle
  )
  return(pp)
}


names(gene_list)
p_intestine = plot_surv("intestine","Intestine gene expression")

p_ar<-plot_surv("ar","Androgen-responsive gene expression")
p_ere<-plot_surv("ere","Estrogen-responsive gene expression")


############################
###分cluster


# geneL = ar_df[ar_df$cluster=="Cluster 3",]$Gene.stable.ID
COAD_tumor_count <- COAD_count[,metadata[metadata$Group=="tumor",]$sample]
plot_surv_new <- function(geneL,figtitle){
  g.tumor <- gsva(as.matrix(COAD_tumor_count), list("tmp"=geneL), verbose=FALSE, parallel.sz=1,kcdf = "Poisson")
  samplelist = names(g.tumor[1,])
  tmp_survival_info <- survival_info[samplelist,c("gender","OS","OS.time")]
  tmp_survival_info$score <- g.tumor[1,]
  data = tmp_survival_info %>%
    survminer::surv_cutpoint(
      time = "OS.time", event = "OS",
      variables = c("score"),
      minprop = 0.25, progressbar = TRUE
    ) %>%
    survminer::surv_categorize(labels = c("Low", "High")) %>%
    data.frame()
  
  tmp_survival_info$type = data$score
  tmp_survival_info$type <- factor(tmp_survival_info$type, levels = c("Low","High"))
  fit <- survfit(Surv(OS.time, OS) ~ type, data = tmp_survival_info)
  
  # library(survminer)
  pp = ggsurvplot(fit, data = tmp_survival_info, pval = TRUE, pval.method = TRUE,palette = "aaas",
                  size = 1.2,
                  font.legend = c(14, "black"),
                  font.x = c(14, "bold", "black"),
                  font.y = c(14, "bold", "black"),
                  font.tickslab = c(12, "bold", "black"),
                  xlab = "Duration overall survival (days)",
                  ncensor.plot = FALSE, # plot the number of censored subjects at time t
                  surv.plot.height = 0.7,
                  risk.table.height = 0.15,
                  ncensor.plot.height = 0.15,
                  legend.labs = c("Low","High"),
                  ggtheme = ggplot2::theme_classic(),
                  legend.title=figtitle
  )
  return(pp)
}

geneL1 = ar_df[ar_df$cluster=="Cluster 1",]$Gene.stable.ID
p_cluster1_ar <- plot_surv_new(geneL1,"Androgen-responsive gene expression in cluster1")
geneL2 = ar_df[ar_df$cluster=="Cluster 2",]$Gene.stable.ID
####cluster2显著
p_cluster2_ar <- plot_surv_new(geneL2,"Androgen-responsive gene expression in cluster2")

geneL3 = ar_df[ar_df$cluster=="Cluster 3",]$Gene.stable.ID
p_cluster3_ar <- plot_surv_new(geneL3,"Androgen-responsive gene expression in cluster3")
geneL4 = ar_df[ar_df$cluster=="Cluster 4",]$Gene.stable.ID
p_cluster4_ar <- plot_surv_new(geneL4,"Androgen-responsive gene expression in cluster4")
geneL5 = ar_df[ar_df$cluster=="Cluster 5",]$Gene.stable.ID
p_cluster5_ar <- plot_surv_new(geneL5,"Androgen-responsive gene expression in cluster5")


egeneL1 = ere_df[ere_df$cluster=="Cluster 1",]$Gene.stable.ID
###cluster1 显著
p_cluster1_ere <- plot_surv_new(egeneL1,"Estrogen-responsive gene expression in cluster1")
egeneL2 = ere_df[ere_df$cluster=="Cluster 2",]$Gene.stable.ID
p_cluster2_ere <- plot_surv_new(egeneL2,"Estrogen-responsive gene expression in cluster2")
egeneL3 = ere_df[ere_df$cluster=="Cluster 3",]$Gene.stable.ID
#cluster3显著
p_cluster3_ere <- plot_surv_new(egeneL3,"Estrogen-responsive gene expression in cluster3")
egeneL4 = ere_df[ere_df$cluster=="Cluster 4",]$Gene.stable.ID
##cluster4 非常显著
p_cluster4_ere <- plot_surv_new(egeneL4,"Estrogen-responsive gene expression in cluster4")
egeneL5 = ere_df[ere_df$cluster=="Cluster 5",]$Gene.stable.ID
p_cluster5_ere <- plot_surv_new(egeneL5,"Estrogen-responsive gene expression in cluster5")



