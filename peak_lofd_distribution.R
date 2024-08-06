



####人与小鼠peak的logf
library(ChIPseeker)
####enhancer_logf2_with_expression.R脚本中已经计算
"/home/shimw/project/conservation/enhancer_count/human_enhancer_DE_annotated.csv"
#############mouse
mouse_peak <- readr::read_tsv("/home/zhihl/Project/CRC/Chip_analysis/dff/version0821/final_result/all_diff_data_H3K27ac_add_symbol.txt")
mouse_peak <- mouse_peak[,c("peak_id","ensemblID","log2FoldChange.tumorVSctrl", "padj.tumorVSctrl")]
names(mouse_peak) = c("m_peak_id","m_gene","m_log2FoldChange","m_padj")
df_data = read.csv("/home/panxl/CRC2/all_cluater")
df_data$cluster = stringr::str_sub(df_data$cluster,0,9)
df_data = df_data[,c("Mouse.gene.stable.ID","Mouse_logFC1",'cluster')]
mouse_peak <- left_join(mouse_peak,df_data,by = c("m_gene"="Mouse.gene.stable.ID"))%>%
  tidyr::drop_na()
readr::write_csv(mouse_peak,"/home/shimw/project/conservation/enhancer_count/mouse_enhancer_DE_annotated.csv")




###########################################################################################################################
####开始
human_peak <- readr::read_csv("/home/shimw/project/conservation/enhancer_count/human_enhancer_DE_annotated.csv")
mouse_peak <- readr::read_csv("/home/shimw/project/conservation/enhancer_count/mouse_enhancer_DE_annotated.csv")
mouse_position <- readr::read_tsv("/home/zhihl/Project/CRC/Chip_analysis/peak_dir_0821/enhancer/enhancer.peak.unique.ID.bed",
                                  col_names = c("chr","start","end","m_peak_id"))
table(mouse_peak$m_peak_id%in%mouse_position$m_peak_id)
mouse_peak <- left_join(mouse_peak, mouse_position)

peak_all <- tibble::tibble(
  "species" = c(rep("human", nrow(human_peak)), rep("mouse", nrow(mouse_peak))),
  "chr"=c(stringr::str_split_fixed(human_peak$h_peak_id, "_",3)[,1], mouse_peak$chr),
  "start" =c(stringr::str_split_fixed(human_peak$h_peak_id, "_",3)[,2],mouse_peak$start),
  "end" =c(stringr::str_split_fixed(human_peak$h_peak_id, "_",3)[,3],mouse_peak$end),
  "closest gene" = c(human_peak$h_gene, mouse_peak$m_gene),
  "closest symbol" = c(),
  "log2FoldChange" = c(human_peak$h_log2FoldChange, mouse_peak$m_log2FoldChange),
  "padj" = c(human_peak$h_padj, mouse_peak$m_padj),
  "cluster" = c(human_peak$cluster, mouse_peak$cluster)
  
)
library(biomaRt)
ensembl37 <- useMart(host="https://grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
attributes <- c("ensembl_gene_id", "external_gene_name")
human_symbol <- getBM(attributes=attributes, filters="ensembl_gene_id", values=peak_all[peak_all$species=="human",]$`closest gene`, mart=ensembl37)

ensembl10 <- useMart(host="https://nov2020.archive.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
mouse_symbol <- getBM(attributes=attributes, filters="ensembl_gene_id", values=peak_all[peak_all$species=="mouse",]$`closest gene`, mart=ensembl10)
symbol_all = rbind(human_symbol, mouse_symbol)
peak_all = left_join(peak_all, symbol_all, by = c("closest gene"="ensembl_gene_id"))
peak_all <- dplyr::rename(peak_all, "symbol" = external_gene_name)

readr::write_csv(peak_all, "/home/shimw/project/conservation/enhancer_count/expression_dif_result.csv")
peak_all <- readr::read_csv("/home/shimw/project/conservation/enhancer_count/expression_dif_result.csv")
filter_peak_all <- peak_all%>%
  dplyr::filter(padj < 0.05)
readr::write_csv(filter_peak_all, "/home/shimw/project/conservation/enhancer_count/expression_dif_result_padj.csv")


####同源peak在peak_NACC.R中已经经过整理了
orth_peak <- readr::read_csv("/home/shimw/project/conservation/NACC/peak_orth_map.csv")
orth_peak[orth_peak$h_peak_id=="chr11_9.6e+07_96000845",]$h_peak_id = "chr11_96000000_96000845"
orth_peak <- orth_peak[,c("h_peak_id", "m_peak_id")]
#####之前已经看过了，同源的peak，基因一定同源，human_peak筛选过必须在cluster内，orth必须有同源，所以peak并没有完全覆盖
human_peak <- human_peak[human_peak$h_peak_id%in%orth_peak$h_peak_id,]


table(orth_peak$m_peak_id%in%mouse_peak$m_peak_id)
orth_peak <- left_join(orth_peak, mouse_peak)%>%tidyr::drop_na()
orth_peak <- left_join(orth_peak, human_peak)%>%tidyr::drop_na()
orth_peak <- left_join(orth_peak, mouse_peak)%>%tidyr::drop_na()

dim(human_peak[human_peak$h_log2FoldChange>1,])
dim(human_peak)

orth_peak <- dplyr::mutate(orth_peak,type = case_when(h_log2FoldChange>1&m_log2FoldChange>1 ~ "both upregulate",
                                  h_log2FoldChange>1&m_log2FoldChange<=1 ~ "only human upregulate",
                                  h_log2FoldChange<=1&m_log2FoldChange>1 ~ "only mouse upregulate",
                                  h_log2FoldChange< -1&m_log2FoldChange< -1 ~ "both donwregulate",
                                  h_log2FoldChange< -1&m_log2FoldChange>= -1 ~ "only human  donwregulate",
                                  h_log2FoldChange>=-1&m_log2FoldChange< -1 ~ "only mouse  donwregulate",
                                  TRUE ~"not change"))
orth_peak <- orth_peak%>%dplyr::filter(type!="not change")

res_sum <- summarise(group_by(orth_peak,cluster,type),peak_num=n())


# res_sum$peak_num <-ifelse(grepl("donwregulate", res_sum$type),-res_sum$peak_num,res_sum$peak_num)


mm = res_sum%>%
  summarise(perc = (peak_num/sum(peak_num))*100)
res_sum$peak_perc <- mm$perc
# res_sum$peak_perc <-ifelse(grepl("donwregulate", res_sum$type),-res_sum$peak_perc,res_sum$peak_perc)

# ggplot(res_sum,)


res_sum$type <- factor(res_sum$type,levels = c("both upregulate","only human upregulate",
                                               "only mouse upregulate","both donwregulate" ,
                                               "only human  donwregulate","only mouse  donwregulate"))
ggplot(data = res_sum, mapping = aes(x = cluster, y = peak_perc)) + geom_bar(stat = 'identity')+
  facet_wrap(~ type, ncol=3)



########使用富集方式计算显著性

calc_ratio <- function(cgroup,ctype){
  genes_fisher_cluster = matrix(c(nrow(orth_peak[orth_peak$cluster==cgroup&orth_peak$type==ctype,]),
                                  nrow(orth_peak[orth_peak$cluster==cgroup&!(orth_peak$type==ctype),]),
                                  nrow(orth_peak[orth_peak$cluster=="Cluster 5"&orth_peak$type==ctype,]),
                                  nrow(orth_peak[orth_peak$cluster=="Cluster 5"&!(orth_peak$type==ctype),])),
                                nrow = 2,
                                dimnames =list(c("condition", "non_condition"),
                                               c("cluster5", "Not cluster5"))
                                
  )
  OR_1 = fisher.test(genes_fisher_cluster, conf.level = 0.95,alternative = "greater")
  return(tibble("boxLabels"=stringr::str_replace(cgroup," ",""),"boxTypes"=ctype,
                "boxOdds"=OR_1$estimate,"boxCILow"=OR_1$conf.int[1],
                "boxCIHigh"=OR_1$conf.int[2],"pval"=OR_1$p.value))
}
cgroups = unique(orth_peak$cluster)
odds_peak_ratios <- purrr::map(cgroups,function(cgroup){
  
  purrr::map(c("both upregulate","only human upregulate",
               "only mouse upregulate","both donwregulate" ,
               "only human  donwregulate","only mouse  donwregulate"),function(ctype){
                 calc_ratio(cgroup,ctype)
               })%>%dplyr::bind_rows()
  
})%>%dplyr::bind_rows()



# ggplot(data = res_sum, mapping = aes(x = cluster, y = peak_perc)) + geom_bar(stat = 'identity')+
#   facet_wrap(~ type, ncol=3)


odds_peak_ratios$boxTypes <- factor(odds_peak_ratios$boxTypes,levels = c("both upregulate","only human upregulate",
                                               "only mouse upregulate","both donwregulate" ,
                                               "only human  donwregulate","only mouse  donwregulate"))



label_map <- c(
  `both upregulate`=sprintf('Human \u2191 Mouse \u2191'),
  `only human upregulate`=sprintf('Human \u2191 Mouse \u2192'),
  `only mouse upregulate`=sprintf('Human \u2192 Mouse \u2191'),
  `both donwregulate`=sprintf('Human \u2193 Mouse \u2193'),
  `only human  donwregulate`=sprintf('Human \u2193 Mouse \u2192'),
  `only mouse  donwregulate`=sprintf('Human \u2192 Mouse \u2193')
)


odds_peak_ratios$star=ifelse(odds_peak_ratios$pval>0.1,"",ifelse(odds_peak_ratios$pval>0.05,"",ifelse(odds_peak_ratios$pval>0.01,"*",ifelse(odds_peak_ratios$pval>0.001,"**","***"))))
# geom_point(size = 2, color = c("#CB9C7A", "#8696a7", "#CDB97D", "#7b8b6f", "#A59B95")) +
p1 <-ggplot(odds_peak_ratios, aes(y = boxOdds, x = boxLabels,fill=boxLabels)) + 
  geom_bar(stat = "identity", width = 0.8)+
  geom_text(aes(y = boxOdds +  1, label = star))+
  facet_wrap(~ boxTypes, ncol=3,labeller = as_labeller(label_map))+
  labs(title = "Enrichment of differential peaks over clusters", y="Odds ratio", x=NULL, fill=NULL) +
  theme_bw(base_family = "sans",base_size = 16,base_line_size = 0.5)+
  scale_fill_manual(values=c("Cluster1" = "#CB9C7A", "Cluster2" = "#8696a7", "Cluster3" = "#CDB97D",
                             "Cluster4" = "#7b8b6f", "Cluster5" = "#A59B95"))+
  theme(axis.title = element_text(size=14),panel.spacing = unit(0.3, "lines"),
        plot.title = element_text(size=14,hjust = 0.5),
        #axis.text = element_text(size=24, color="black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        legend.title = element_blank(),
        legend.text = element_text(size=12), legend.key.size = unit(0.6,"cm"),
        # axis.ticks.length.x = unit(2,"mm"),
        # axis.ticks.length.y = unit(2,"mm"),
        # strip.background = element_blank(),
        legend.box = "horizontal",
        strip.text.x = element_text(size = 12)
  )

# ggsave("/home/shimw/project/conservation/peaks_enrichment_clusters.pdf",p1,width=9,units = "in",height = 5)

cairo_pdf("/home/shimw/project/conservation/peaks_enrichment_clusters.pdf",width=9,height = 5)
p1
dev.off()




























