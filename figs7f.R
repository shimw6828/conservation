cms_type = read.csv("/home/panxl/CRC2/human_cmstypes.csv",check.names = F)
human_data <- read.csv("/home/shimw/project/conservation/filter_human_data.csv",row.names = 1, check.names = F)
human_info <- read.csv("/home/panxl/CRC2/filter_human_info.csv",row.names = 1, check.names = F)
normal_sample = human_info[human_info$Group=="Normal",]$sample
CMS1_sample = cms_type[cms_type$type=="CMS1",]$barcode
CMS2_sample = cms_type[cms_type$type=="CMS2",]$barcode
CMS3_sample = cms_type[cms_type$type=="CMS3",]$barcode
CMS4_sample = cms_type[cms_type$type=="CMS4",]$barcode
cms_types = data.frame(sample = c(normal_sample, CMS1_sample,CMS2_sample,CMS3_sample,CMS4_sample),
                       group=c(rep("Normal", length(normal_sample)), rep("CMS1", length(CMS1_sample)),
                               rep("CMS2", length(CMS2_sample)), rep("CMS3", length(CMS3_sample)),
                               rep("CMS4", length(CMS4_sample))))


cluster_type = read.csv("/home/panxl/CRC2/all_cluater",check.names = F)
cluster1_gene = cluster_type$Gene.stable.ID[cluster_type$cluster=="Cluster 1: Up-regulated in both species"]
cluster2_gene = cluster_type$Gene.stable.ID[cluster_type$cluster=="Cluster 2: Up-regulated in mouse reverse in human"]
cluster3_gene = cluster_type$Gene.stable.ID[cluster_type$cluster=="Cluster 3: Down-regulated in both species"]
cluster4_gene = cluster_type$Gene.stable.ID[cluster_type$cluster=="Cluster 4: Up-regulated in human reverse in mouse"]
cluster5_gene = cluster_type$Gene.stable.ID[cluster_type$cluster=="Cluster 5: Not-significant in both species"]





#remotes::install_github('rpkgs/gg.layers')
library(Ipaper)
library(ggplot2)

plot_box = function(sp, sp_title, y1, y2){
  data = human_data[sp, c(normal_sample, CMS1_sample,CMS2_sample,CMS3_sample,CMS4_sample)]%>%
    tibble::rownames_to_column("gene")%>%
    tidyr::gather(sample, exp, -gene)%>%
    left_join(cms_types)%>%
    group_by(gene, group)%>%
    summarize(exp = mean(exp, na.rm = TRUE))
kk = t.test(data[data$group%in%c("CMS1","CMS2","CMS3"),]$exp,data[data$group%in%c("CMS4"),]$exp)
p_sig = case_when(kk$p.value>0.05 ~ "ns",
                  kk$p.value>0.01 ~ "*",
                  kk$p.value>0.001 ~ "**",
                  kk$p.value>0.0001 ~ "***",
                  TRUE ~ "****")

p = ggplot(data = na.omit(data),aes(x = group,y = exp,fill=group)) +
  # geom_violin(adjust = .5) +
  geom_boxplot(linetype = 2, outlier.shape = NA)+
  stat_boxplot(aes(ymin=..lower..,ymax=..upper..), outlier.shape = NA) +
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..),width=0.3,color="black")+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..),width=0.3,color="black")+
  geom_segment(aes(x=0.8, y=y1, xend=3.2, yend=y1),size=0.3)+
  geom_segment(aes(x=3.8, y=y1, xend=4.2, yend=y1),size=0.3)+
  geom_segment(aes(x=2, y=y1, xend=2, yend=y2),size=0.3)+
  geom_segment(aes(x=4, y=y1, xend=4, yend=y2),size=0.3)+
  geom_segment(aes(x=2, y=y2, xend=4, yend=y2),size=0.3)+
  annotate("text", x=3, y=y2+1, label=p_sig,size=5)+
  theme_classic(base_family = "sans",base_size = 18,base_line_size = 0.6)+
  scale_fill_manual(values=c('CMS1'="#b97573", 
                             'CMS2'="#686854",
                             'CMS3'="#5a5e66", 
                             'CMS4'="#816147", 
                             'Normal'="#536F2E"))+
  labs(y="Expression", x=NULL, fill=NULL, title=sp_title) +
  theme(title = element_text(size = 14,
                             color = "black"),
        legend.position = "none",
        axis.title = element_text(size=14),
        #axis.text = element_text(size=24, color="black"),
        axis.text.x = element_text(size = 14,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.text.y = element_text(size = 14,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.ticks.length.x = unit(2,"mm"),
        axis.ticks.length.y = unit(2,"mm"))
  return(p)
}

p1 = plot_box(cluster1_gene, "Gene expression in Cluster 1", 12, 19)
ggsave("/home/shimw/project/conservation/CMS_expression/cluster1.png",p1, width = 4.5,height = 4.5)
p2 = plot_box(cluster2_gene, "Gene expression in Cluster 2",12, 19)
ggsave("/home/shimw/project/conservation/CMS_expression/cluster2.png",p2, width = 4.5,height = 4.5)
p3 = plot_box(cluster3_gene, "Gene expression in Cluster 3",15.5, 19)
ggsave("/home/shimw/project/conservation/CMS_expression/cluster3.png",p3, width = 4.5,height = 4.5)
p4 = plot_box(cluster4_gene, "Gene expression in Cluster 4",13, 19)
ggsave("/home/shimw/project/conservation/CMS_expression/cluster4.png",p4, width = 4.5,height = 4.5)
p5 = plot_box(cluster5_gene, "Gene expression in Cluster 5",12, 19)
ggsave("/home/shimw/project/conservation/CMS_expression/cluster5.png",p5, width = 4.5,height = 4.5)
