#####################################################
###house keeping gene的富集
df_data = read.csv("/home/panxl/CRC2/all_cluater")
df_data$cluster = stringr::str_sub(df_data$cluster,0,9)
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
odds_house_ratios <- purrr::map(cgroups,function(cgroup){
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
library(ggsci)
p <- ggplot(odds_house_ratios, aes(x = boxOdds, y = boxLabels)) + 
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


############################################################
#####合并一起画
###肠道特异性计算
df_intestine = read.csv("/home/shimw/project/enhancer_map/conserved/colon-specific.csv")
H0M_logFC1= read.csv("/home/panxl/CRC2/all_cluater")
H0M_logFC1$cluster=gsub("Cluster 1: Up-regulated in both species","Cluster 1",H0M_logFC1$cluster)
H0M_logFC1$cluster=gsub("Cluster 2: Up-regulated in mouse reverse in human","Cluster 2",H0M_logFC1$cluster)
H0M_logFC1$cluster=gsub("Cluster 3: Down-regulated in both species","Cluster 3",H0M_logFC1$cluster)
H0M_logFC1$cluster=gsub("Cluster 4: Up-regulated in human reverse in mouse","Cluster 4",H0M_logFC1$cluster)
H0M_logFC1$cluster=gsub("Cluster 5: Not-significant in both species","Cluster 5",H0M_logFC1$cluster)
df_data=H0M_logFC1
colnames(df_data)
cgroups =unique(df_data$cluster)
odds_intestine_ratios <- purrr::map(cgroups,function(cgroups){
  genes_fisher_cluster = matrix(c(nrow(df_data[df_data$cluster==cgroups&df_data$Mouse.gene.stable.ID%in%df_intestine$Gene,]),
                                  nrow(df_data[df_data$cluster==cgroups&!(df_data$Mouse.gene.stable.ID%in%df_intestine$Gene),]),
                                  nrow(df_data[df_data$cluster!=cgroups&df_data$Mouse.gene.stable.ID%in%df_intestine$Gene,]),
                                  nrow(df_data[df_data$cluster!=cgroups&!(df_data$Mouse.gene.stable.ID%in%df_intestine$Gene),])),
                                nrow = 2,
                                dimnames =list(c("housekeeping", "non_housekeeping"),
                                               c("cluster1", "Not cluster1")))
  
  
  OR_1 = fisher.test(genes_fisher_cluster, conf.level = 0.95)
  #return(tibble("boxLabels"=stringr::str_replace(cgroups," ",""),"boxOdds"=OR_1$estimate,"boxCILow"=OR_1$conf.int[1],"boxCIHigh"=OR_1$conf.int[2]))
  return(tibble("boxLabels"=cgroups,"boxOdds"=OR_1$estimate,"boxCILow"=OR_1$conf.int[1],"boxCIHigh"=OR_1$conf.int[2]))
})%>%dplyr::bind_rows()
######################
###物种特异
df_species = read.csv("/home/shimw/project/enhancer_map/conserved/species-specific.csv")#744
ensembl_ID_change=read.csv("/home/panxl/ensembl_37gene.txt")#28255
H0M_logFC2=left_join(H0M_logFC1,ensembl_ID_change,c("Gene.stable.ID"="Gene.stable.ID"))
H0M_logFC3=na.omit(H0M_logFC2)
df_data=H0M_logFC3
colnames(df_species)

cgroups =unique(df_data$cluster)
#cgroups ="Cluster 1"
odds_specific_ratios <- purrr::map(cgroups,function(cgroups){
  genes_fisher_cluster = matrix(c(nrow(df_data[df_data$cluster==cgroups&df_data$Gene.name%in%df_species$gene_id,]),
                                  nrow(df_data[df_data$cluster==cgroups&!(df_data$Gene.name%in%df_species$gene_id),]),
                                  nrow(df_data[df_data$cluster!=cgroups&df_data$Gene.name%in%df_species$gene_id,]),
                                  nrow(df_data[df_data$cluster!=cgroups&!(df_data$Gene.name%in%df_species$gene_id),])),
                                nrow = 2,
                                dimnames =list(c("housekeeping", "non_housekeeping"),
                                               c("cluster1", "Not cluster1")))
  
  
  OR_1 = fisher.test(genes_fisher_cluster, conf.level = 0.95)
  #return(tibble("boxLabels"=stringr::str_replace(cgroups," ",""),"boxOdds"=OR_1$estimate,"boxCILow"=OR_1$conf.int[1],"boxCIHigh"=OR_1$conf.int[2]))
  return(tibble("boxLabels"=cgroups,"boxOdds"=OR_1$estimate,"boxCILow"=OR_1$conf.int[1],"boxCIHigh"=OR_1$conf.int[2]))
})%>%dplyr::bind_rows()

###########
###house keeping gene上面已经算了
#合并
odds_specific_ratios$type = "specific"
odds_intestine_ratios$type = "intestine"
odds_house_ratios$type = "house"



label_map <- c(
  `specific`="Species specific",
  `intestine`="Intestine specific",
  `house`="House keeping"
)
all_odds = dplyr::bind_rows(list(odds_specific_ratios,odds_intestine_ratios,odds_house_ratios))
all_odds$boxLabels = stringr::str_replace(all_odds$boxLabels," ","")
all_odds$type = factor(all_odds$type,levels = c("specific","intestine","house"))


scaleFUN <- function(x) sprintf("%.0f", x)
# ggplot(all_odds, aes(x = boxOdds, y = boxLabels,color=boxLabels)) + 
#   geom_vline(aes(xintercept = 1.0), size = .25, linetype = "dashed") + 
#   geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = 1, height = 
#                    .15 ) +
#   geom_point(size = 2) +
#   facet_wrap(~ type, ncol=3,labeller = as_labeller(label_map),scales = 'free_x') +
#   labs(title = NULL, y=NULL, x="Odds ratio", fill=NULL) +
#   scale_x_continuous(labels=scaleFUN)+
#   scale_color_manual(values=c("Cluster1" = "#CB9C7A", "Cluster2" = "#8696a7", "Cluster3" = "#CDB97D",
#                              "Cluster4" = "#7b8b6f", "Cluster5" = "#A59B95"))+
#   theme_bw(base_family = "sans",base_size = 16,base_line_size = 0.5)+
#   theme(axis.title = element_text(size=18),
#         panel.spacing = unit(0.3, "lines"),
#         axis.title.y = element_blank(),plot.title = element_text(size=16,hjust = 0.5),
#         #axis.text = element_text(size=24, color="black"),
#         axis.text.x = element_text(size = 16,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
#         axis.text.y = element_text(size = 16,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
#         legend.title = element_blank(),
#         legend.text = element_text(size=12), legend.key.size = unit(0.6,"cm"),
#         axis.ticks.length.x = unit(2,"mm"),
#         axis.ticks.length.y = unit(2,"mm"),
#         # strip.background = element_blank(),
#         strip.text = element_text(size = 18),
#         panel.spacing.x = unit(1, "mm"),
#         legend.position = "none"
#   )
p1 = ggplot(all_odds, aes(x = boxOdds, y = boxLabels,color=boxLabels)) + 
  geom_vline(aes(xintercept = 1.0), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = 1, height = 
                   .15 ) +
  geom_point(size = 2) +
  facet_wrap(~ type, ncol=3,labeller = as_labeller(label_map),scales = 'free_x') +
  labs(title = NULL, y=NULL, x="Odds ratio", fill=NULL) +
  scale_x_continuous(labels=scaleFUN)+
  scale_color_manual(values=c("Cluster1" = "#CB9C7A", "Cluster2" = "#8696a7", "Cluster3" = "#CDB97D",
                              "Cluster4" = "#7b8b6f", "Cluster5" = "#A59B95"))+
  theme_bw(base_family = "sans",base_size = 16,base_line_size = 0.5)+
  theme(axis.text.y = element_text(size = 14,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        strip.text.x = element_text(size = 13),
        panel.spacing.x = unit(2, "mm"),
        legend.position = "none"
  )

cairo_pdf("/home/shimw/project/conservation/specific_housekeeping_odds.pdf",width=9.5,height = 3.5)
p1
dev.off()






