



# kk = read.csv("/home/shimw/project/conservation/DN_DS.csv",header=T)
kk = read.csv("/home/shimw/project/conservation/MouseToHuman.csv",header=T)

kk = kk%>%tidyr::drop_na()
kk$dn_ds = kk$dN/kk$dS

kk = kk[,c("Ensembl.Gene.ID","dn_ds")]


dim(kk)
df_data = read.csv("/home/panxl/CRC2/all_cluater")
table(df_data$Mouse.gene.stable.ID%in%kk$Ensembl.Gene.ID)
df_data$cluster = stringr::str_sub(df_data$cluster,0,9)





# nn = left_join(df_data,kk,by = c("Gene.stable.ID" = "Ensembl.Gene.ID"))
nn = left_join(df_data,kk,by = c("Mouse.gene.stable.ID" = "Ensembl.Gene.ID"))





nn = tidyr::drop_na(nn)
cor.test(nn$Human_logFC1,nn$dn_ds)

nn[nn$Human_logFC1<1,]$dn_ds
nn[nn$Human_logFC1>1,]$dn_ds
t.test(nn[nn$Human_logFC1<1,]$dn_ds,nn[nn$Human_logFC1>1,]$dn_ds)


dim(nn)
nn$dn_ds[nn$cluster=="Cluster 5"]
summary(nn$dn_ds[nn$cluster=="Cluster 1"])
summary(nn$dn_ds[nn$cluster=="Cluster 2"])
summary(nn$dn_ds[nn$cluster=="Cluster 3"])
summary(nn$dn_ds[nn$cluster=="Cluster 4"])
summary(nn$dn_ds[nn$cluster=="Cluster 5"])

nn = filter(nn,cluster!="Cluster 5")
nn = dplyr::mutate(nn,"NewCluster" = ifelse(nn$cluster%in%c("Cluster 1","Cluster 3"),"cluster13","cluster24"))
nn= filter(nn,dn_ds<1)

my_comparisons <- list(c("cluster13", "cluster24"))
# my_comparisons <- list(c("Cluster 1", "Cluster 2"), c("Cluster 1", "Cluster 3"), c("Cluster 1", "Cluster 4"), c("Cluster 1", "Cluster 5"))
p1 = ggplot(data = nn,aes(x = NewCluster,y = dn_ds,fill=NewCluster)) +
  geom_violin(adjust = .5) +geom_boxplot(width=0.1, fill="white")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif", tip.length=0.02,method.args = list(alternative = "less"), method = "t.test") +
  theme_classic(base_family = "sans",base_size = 18,base_line_size = 1.1)+
  scale_fill_manual(values=c("cluster13" = "#CB9C7A", "cluster24" = "#8696a7"))+
  labs(title = NULL, y="dN/dS", x=NULL, fill=NULL) +
  scale_x_discrete(labels = c("cluster13" = "Conservation", "cluster24" = "Divergence"))+
  theme(legend.position = "none",plot.title = element_text(hjust = 0,size=16),
        axis.title=element_text(size=18),
        axis.text=element_text(size=15))


# 
# ggplot(data = nn,aes(x = NewCluster,y = dn_ds,fill=NewCluster)) +
#   geom_violin(adjust = .5) +geom_boxplot(width=0.1, fill="white")+
#   stat_compare_means(comparisons = my_comparisons,label = "p.signif", tip.length=0.02,method.args = list(alternative = "less"), method = "t.test") +
#   theme_classic(base_family = "sans",base_size = 18,base_line_size = 1.1)+
#   scale_fill_manual(values=c("Cluster 1" = "#CB9C7A", "Cluster 2" = "#8696a7", "Cluster 3" = "#CDB97D",
#                              "Cluster 4" = "#7b8b6f", "Cluster 5" = "#A59B95"))+
#   labs(title = "dN/dS", y="dN/dS", x=NULL, fill=NULL) +
#   # scale_x_discrete(labels = c("cluster 1" = "Cluster1", "cluster 2" = "Cluster2", "cluster 3" = "Cluster3",
#   #                             "cluster 4" = "Cluster4", "cluster 5" = "Cluster5"))+
#   theme(legend.position = "none",plot.title = element_text(hjust = 0,size=16),
#         axis.title=element_text(size=18),
#         axis.text=element_text(size=15))






#####################################################
###house keeping gene的富集

##转id
library("biomaRt")
HK_genes <- readr::read_tsv("/home/shimw/project/conservation/HK_genes.txt",col_names = c("symbol","ncbi_id"))
ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "grch37.ensembl.org")
symboltoid = getBM(attributes=c('ensembl_gene_id', 'external_gene_name'), 
      filters = 'external_gene_name', 
      values = HK_genes$symbol, 
      mart = ensembl)
HK_df = df_data[df_data$Gene.stable.ID%in%symboltoid$ensembl_gene_id,]
table(HK_df$cluster)["Cluster 1"]
table(HK_df$cluster)["Cluster 5"]
####odds ratio
cgroup="Cluster 1"
cgroups = c("Cluster 1","Cluster 2","Cluster 3","Cluster 4","Cluster 5")
odds_ratios <- purrr::map(cgroups,function(cgroup){
  genes_fisher_cluster = matrix(c(nrow(df_data[df_data$cluster==cgroup&df_data$Gene.stable.ID%in%symboltoid$ensembl_gene_id,]),
                                  nrow(df_data[df_data$cluster==cgroup&!(df_data$Gene.stable.ID%in%symboltoid$ensembl_gene_id),]),
                                  nrow(df_data[df_data$cluster!=cgroup&df_data$Gene.stable.ID%in%symboltoid$ensembl_gene_id,]),
                                  nrow(df_data[df_data$cluster!=cgroup&!(df_data$Gene.stable.ID%in%symboltoid$ensembl_gene_id),])),
                                nrow = 2,
                                dimnames =list(c("housekeeping", "non_housekeeping"),
                                               c("cluster5", "Not cluster5"))
                                
  )
  OR_1 = fisher.test(genes_fisher_cluster, conf.level = 0.95)
  return(tibble("boxLabels"=stringr::str_replace(cgroup," ",""),"boxOdds"=OR_1$estimate,"boxCILow"=OR_1$conf.int[1],"boxCIHigh"=OR_1$conf.int[2]))
})%>%dplyr::bind_rows()
c("Cluster 1" = "#CB9C7A", "Cluster 2" = "#8696a7", "Cluster 3" = "#CDB97D",
  "Cluster 4" = "#7b8b6f", "Cluster 5" = "#A59B95")
p <- ggplot(odds_ratios, aes(x = boxOdds, y = boxLabels)) + 
  geom_vline(aes(xintercept = 1.0), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = 1, height = 
                   .15, color = c("#CB9C7A", "#8696a7", "#CDB97D", "#7b8b6f", "#A59B95")) +
  geom_point(size = 2, color = c("#CB9C7A", "#8696a7", "#CDB97D", "#7b8b6f", "#A59B95")) +
  scale_color_uchicago() +labs(title = "House keeping", y=NULL, x="Odds ratio", fill=NULL) +
  theme(axis.title = element_text(size=18),
        axis.title.y = element_blank(),plot.title = element_text(size=16,hjust = 0.5),
        #axis.text = element_text(size=24, color="black"),
        axis.text.x = element_text(size = 16,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.text.y = element_text(size = 16,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        legend.title = element_blank(),
        legend.text = element_text(size=12), legend.key.size = unit(0.6,"cm"),
        axis.ticks.length.x = unit(2,"mm"),
        axis.ticks.length.y = unit(2,"mm")
  )










