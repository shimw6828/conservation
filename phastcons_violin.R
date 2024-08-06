library(GenomicFeatures)
library(GenomicRanges)
library(GenomicScores)
library(dplyr)
# txdb <-  makeTxDbFromGFF("/home/shimw/project/conservation/gencode.v32lift37.annotation.gtf", format = "gtf")
txdb = makeTxDbFromGFF("/home/mwshi/project/conservation/gencode.vM20.basic.annotation.gtf", format="gtf")
# txdb = makeTxDbFromGFF("/home/mwshi/project/conservation/gencode.v32lift37.annotation.gtf", format="gtf")

#longest transcript
all_transcript = transcriptsBy(txdb,by = "gene")
genetotx = purrr::map(names(all_transcript),function(x){
  y = all_transcript[[x]]
  tx_id = y[order(width(y),decreasing = TRUE)[1],]$tx_name
  return(tibble::tibble("gene"=x,"tx_id" = tx_id))
})%>%dplyr::bind_rows()
genetotx$gene = gsub("\\..*", "", genetotx$gene)
# readr::write_csv(genetotx,"/home/shimw/project/conservation/phastcons/genetotx.csv")
readr::write_csv(genetotx,"/home/mwshi/project/conservation/phastcons/genetotx.csv")
# genetotx <- readr::read_csv("/home/shimw/project/conservation/phastcons/genetotx.csv")
genetotx <- readr::read_csv("/home/mwshi/project/conservation/phastcons/genetotx.csv")
rm(all_transcript)
# df_data = read.csv("/home/panxl/CRC2/all_cluater")
df_data = read.csv("/home/mwshi/project/conservation/all_cluater")
df_data = df_data[,c("Mouse.gene.stable.ID","Gene.stable.ID","cluster")]
names(df_data)[1] = "gene"
df_data = dplyr::left_join(df_data,genetotx)
df_data$cluster = stringr::str_sub(df_data$cluster,0,9)

##################################################################################
#准备函数与phastCons对象
# library(phastCons100way.UCSC.hg19)
# phast <- phastCons100way.UCSC.hg19
phast <- getGScores("phastCons60way.UCSC.mm10")
#根据Grange对象计算保守性值
calc_score <- function(ob_grange){
  kk = c()
  for (x in names(ob_grange)) {
    # print(x)
    y=ob_grange[[x]]
    if (length(y)==0) {
      kk = c(kk,NA)
      next
    }
    tmp_y = gscores(phast,y,summaryFun = sum)
    tmp_y = tmp_y[!is.na(tmp_y$default),]
    tmp_score <- sum(tmp_y$default)/sum(width(tmp_y))
    kk = c(kk,tmp_score)
  }
  kk = unlist(kk)
  return(tibble::tibble("tx_id"=names(ob_grange),"score"=kk))
}
#exon
all_exon = exonsBy(txdb, by="tx", use.names = T)
all_exon = all_exon[df_data$tx_id]
length(all_exon)
exon_scores <- calc_score(all_exon)
# readr::write_csv(exon_scores,"/home/shimw/project/conservation/phastcons/exon_scores.csv")
readr::write_csv(exon_scores,"/home/mwshi/project/conservation/phastcons/exon_scores.csv")

#promoter TSS
all_tss <- promoters(genes(txdb), upstream=1500)
all_tss$gene_id =  gsub("\\..*", "", all_tss$gene_id )
all_tss = all_tss[all_tss$gene_id%in%df_data$gene,]
tss_scores = gscores(phast,all_tss,summaryFun = mean)
tss_scores = data.frame("gene"=tss_scores$gene_id,"score" = tss_scores$default,stringsAsFactors = F)
readr::write_csv(tss_scores,"/home/mwshi/project/conservation/phastcons/tss_scores.csv")
#fiveUTR
all_fiveUTRs = fiveUTRsByTranscript(txdb,use.names = T)
#部分转录本没有UTR
#ENSMUST00000080797.7
all_fiveUTRs = all_fiveUTRs[df_data$tx_id[df_data$tx_id%in%names(all_fiveUTRs)]]
fiveUTRs_scores <- calc_score(all_fiveUTRs)
# readr::write_csv(fiveUTRs_scores,"/home/shimw/project/conservation/phastcons/fiveUTRs_scores.csv")
readr::write_csv(fiveUTRs_scores,"/home/mwshi/project/conservation/phastcons/fiveUTRs_scores.csv")
#threeUTR
all_threeUTRs = threeUTRsByTranscript(txdb,use.names = T)
all_threeUTRs = all_threeUTRs[df_data$tx_id[df_data$tx_id%in%names(all_threeUTRs)]]
threeUTRs_scores <- calc_score(all_threeUTRs)
# readr::write_csv(threeUTRs_scores,"/home/shimw/project/conservation/phastcons/threeUTRs_scores.csv")
readr::write_csv(threeUTRs_scores,"/home/mwshi/project/conservation/phastcons/threeUTRs_scores.csv")



#introns
all_intron = intronsByTranscript(txdb,use.names = T)
all_intron = all_intron[df_data$tx_id[df_data$tx_id%in%names(all_intron)]]
intron_scores <- calc_score(all_intron)
# readr::write_csv(intron_scores,"/home/shimw/project/conservation/phastcons/intron_scores.csv")
readr::write_csv(intron_scores,"/home/mwshi/project/conservation/phastcons/intron_scores.csv")

#splicing site
#就是intron两端
kk = c()
for (x in names(all_intron)) {
  y=all_intron[[x]]
  if (length(y)==0) {
    kk = c(kk,NA)
    next
  }
  y = GRanges(seqnames=seqnames(y)[1],IRanges(start=c(start(y),end(y))))
  tmp_y = gscores(phast,y,summaryFun = mean)
  tmp_y = tmp_y[!is.na(tmp_y$default),]
  tmp_score <- sum(tmp_y$default)/sum(width(tmp_y))
  kk = c(kk,tmp_score)
}
splicing_score <- tibble::tibble("tx_id"=names(all_intron),"score"=kk)
readr::write_csv(splicing_score,"/home/mwshi/project/conservation/phastcons/splicing_score.csv")


#######################################################################
#####画图
library(ggpubr)
my_comparisons <- list(c("Cluster 1", "Cluster 2"), c("Cluster 1", "Cluster 3"), c("Cluster 1", "Cluster 4"), c("Cluster 1", "Cluster 5"))
plot_phastcons <- function(phastdf,fig_title){
  p1 <- ggplot(data = na.omit(phastdf),aes(x = cluster,y = score,fill=cluster)) +
    geom_violin(adjust = .5) +geom_boxplot(width=0.1, fill="white")+
    stat_compare_means(comparisons = my_comparisons,label = "p.signif", tip.length=0.02,method.args = list(alternative = "greater"), method = "t.test") +
    theme_classic(base_family = "sans",base_size = 18,base_line_size = 1.1)+
    scale_fill_manual(values=c("Cluster 1" = "#CB9C7A", "Cluster 2" = "#8696a7", "Cluster 3" = "#CDB97D",
                               "Cluster 4" = "#7b8b6f", "Cluster 5" = "#A59B95"))+
    labs(title = fig_title, y="Phastcons score", x=NULL, fill=NULL) +
    scale_x_discrete(labels = c("Cluster 1" = "Cluster1", "Cluster 2" = "Cluster2", "Cluster 3" = "Cluster3",
                                "Cluster 4" = "Cluster4", "Cluster 5" = "Cluster5"))+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5,size=16),
          axis.title=element_text(size=18),
          axis.text=element_text(size=15))
  return(p1)
}
genetotx <- readr::read_csv("/home/shimw/project/conservation/phastcons/genetotx.csv")
df_data = read.csv("/home/panxl/CRC2/all_cluater")
# df_data = read.csv("/home/mwshi/project/conservation/all_cluater")
df_data = df_data[,c("Mouse.gene.stable.ID","Gene.stable.ID","cluster")]
names(df_data)[1] = "gene"
df_data = dplyr::left_join(df_data,genetotx)
df_data$cluster = stringr::str_sub(df_data$cluster,0,9)




#####exon
exon_scores <- readr::read_csv("/home/shimw/project/conservation/phastcons/exon_scores.csv")
exon_scores = left_join(exon_scores,df_data)
p1=plot_phastcons(exon_scores,"Exons")
ggsave("/home/shimw/project/conservation/phastcons/exon_phastcons.pdf",p1,width=6,height = 5)

#####promoter TSS
tss_scores <- readr::read_csv("/home/shimw/project/conservation/phastcons/tss_scores.csv")
tss_scores = left_join(tss_scores,df_data)
p1=plot_phastcons(tss_scores,"Promoter (2kb upstream TSS)")
ggsave("/home/shimw/project/conservation/phastcons/promoter_phastcons.pdf",p1,width=6,height = 5)

tss_scores = na.omit(tss_scores)
t.test(tss_scores[tss_scores$cluster=="Cluster 1",]$score,tss_scores[tss_scores$cluster=="Cluster 4",]$score)


#####fiveUTR
fiveUTRs_scores <- readr::read_csv("/home/shimw/project/conservation/phastcons/fiveUTRs_scores.csv")
fiveUTRs_scores = left_join(fiveUTRs_scores,df_data)
p1=plot_phastcons(fiveUTRs_scores,"5′ UTR")
# ggsave("/home/shimw/project/conservation/phastcons/fiveUTR_phastcons.png",p1,width=6,height = 5)
Cairo::CairoPDF("/home/shimw/project/conservation/phastcons/fiveUTR_phastcons.pdf",width=6,height = 5)
p1
dev.off()


#####threeUTRs
threeUTRs_scores <- readr::read_csv("/home/shimw/project/conservation/phastcons/threeUTRs_scores.csv")
threeUTRs_scores = left_join(threeUTRs_scores,df_data)
p1=plot_phastcons(threeUTRs_scores,"3′ UTR")
# ggsave("/home/shimw/project/conservation/phastcons/threeUTR_phastcons.pdf",p1,width=6,height = 5)
Cairo::CairoPDF("/home/shimw/project/conservation/phastcons/threeUTR_phastcons.pdf",width=6,height = 5)
p1
dev.off()
#####introns
intron_scores <- readr::read_csv("/home/shimw/project/conservation/phastcons/intron_scores.csv")
intron_scores = left_join(intron_scores,df_data)
p1=plot_phastcons(intron_scores,"Introns")
ggsave("/home/shimw/project/conservation/phastcons/intron_phastcons.pdf",p1,width=6,height = 5)

#####splicing site
splicing_score <- readr::read_csv("/home/shimw/project/conservation/phastcons/splicing_score.csv")
splicing_score = left_join(splicing_score,df_data)
# splicing_score = splicing_score[,c("tx_id","score","gene","cluster")]
p1=plot_phastcons(splicing_score,"Splice Sites")
ggsave("/home/shimw/project/conservation/phastcons/splicing_phastcons.pdf",p1,width=6,height = 5)


#####all 
exon_scores$type = "exon"
tss_scores$type = "tss"
fiveUTRs_scores$type = "fiveUTRs"
threeUTRs_scores$type = "threeUTRs"
intron_scores$type = "intron"
splicing_score$type = "splicing"

all_score <- dplyr::bind_rows(list(exon_scores,tss_scores,
                                   fiveUTRs_scores,threeUTRs_scores,
                                   intron_scores,splicing_score))%>%
  tidyr::drop_na()

all_score$type = factor(all_score$type,
                        levels = c("exon","tss","fiveUTRs",
                                   "threeUTRs","intron","splicing"))


label_map <- c(
  `exon`="Exons",
  `tss`="Promoter (2kb upstream TSS)",
  `fiveUTRs`="5′ UTR",
  `threeUTRs`="3′ UTR",
  `intron`="Introns",
  `splicing`="Splice Sites"
)


p1 <- ggplot(data = na.omit(all_score),aes(x = cluster,y = score,fill=cluster)) +
  geom_violin(adjust = .5) +geom_boxplot(width=0.1, fill="white")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif", tip.length=0.02,method.args = list(alternative = "greater"), method = "t.test") +
  theme_classic(base_family = "sans",base_size = 18,base_line_size = 1.1)+
  facet_wrap(~ type, ncol=3,labeller = as_labeller(label_map))+
  scale_fill_manual(values=c("Cluster 1" = "#CB9C7A", "Cluster 2" = "#8696a7", "Cluster 3" = "#CDB97D",
                             "Cluster 4" = "#7b8b6f", "Cluster 5" = "#A59B95"))+
  labs( y="Phastcons score", x=NULL, fill=NULL) +
  scale_x_discrete(labels = c("Cluster 1" = "Cluster1", "Cluster 2" = "Cluster2", "Cluster 3" = "Cluster3",
                              "Cluster 4" = "Cluster4", "Cluster 5" = "Cluster5"))+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5,size=16),
        axis.title=element_text(size=18),
        axis.text=element_text(size=12),
        strip.background = element_blank())
# ggsave("/home/shimw/project/conservation/phastcons/all_phastcons.pdf",p1,width=12,height = 8)


cairo_pdf("/home/shimw/project/conservation/phastcons/all_phastcons.pdf",width=12,height = 8)
p1
dev.off()



######################################################################################
#####state change情况与不变情况的score比较
library(phastCons100way.UCSC.hg19)
phast <- phastCons100way.UCSC.hg19
library('rtracklayer') 




"/NAS/luozh/CRC_conservation/step10_chromatin_state/state_compare/ortholog_state/phastcons/human_mouse_state_not_conservation_loci.txt"
human_mouse_state_not_conservation <- import("/NAS/luozh/CRC_conservation/step10_chromatin_state/state_compare/ortholog_state/phastcons/human_mouse_state_not_conservation_loci.txt", format="bed")
human_mouse_state_not_conservation <- annotatePeak(human_mouse_state_not_conservation, tssRegion=c(-3000, 3000), TxDb=hg19)
human_mouse_state_not_conservation <- as.GRanges(human_mouse_state_not_conservation)
not_conservation_score <- gscores(phast,human_mouse_state_not_conservation)


human_mouse_state_conservation <- import("/NAS/luozh/CRC_conservation/step10_chromatin_state/state_compare/ortholog_state/phastcons/human_mouse_state_conservation_loci.txt", format="bed")
conservation_score <- gscores(phast,human_mouse_state_conservation)
not_conservation_score$default
t.test(not_conservation_score$default,conservation_score$default)
all_score <- data.frame("score" = c(conservation_score$default,not_conservation_score$default),
                        "type" = c(rep("Conservation",length(conservation_score$default)),
                                   rep("Differentiation",length(not_conservation_score$default)))
                        )

my_comparisons <- list(c("Conservation", "Differentiation"))
p1 <-ggplot(data = all_score,aes(x = type,y = score,fill=type)) +
  geom_violin(adjust = .5) +geom_boxplot(width=0.1, fill="white")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif", tip.length=0.02,method.args = list(alternative = "greater"), method = "t.test") +
  theme_classic(base_family = "sans",base_size = 18,base_line_size = 1.1)+
  scale_fill_manual(values=c("Conservation" = "#CB9C7A", "Differentiation" = "#8696a7"))+
  labs(title = "state region", y="Phastcons score", x=NULL, fill=NULL) +
  scale_x_discrete(labels = c("Conservation" = "Conservation", "Differentiation" = "Differentiation"))+
  theme(legend.position = "none",plot.title = element_text(hjust = 0,size=16),
        axis.title=element_text(size=18),
        axis.text=element_text(size=15))
p1
summary(all_score)
ggsave("/home/shimw/project/conservation/phastcons/state_change_phastcons.pdf",p1,width=12,height = 8)



#############################################################################################
###state change phastcons
library(phastCons100way.UCSC.hg19)
phast <- phastCons100way.UCSC.hg19
df_active = read.csv("/NAS/luozh/CRC_conservation/step10_chromatin_state/state_compare/ortholog_state/mouse_quiescent_ActiveEnhancer_loci.txt", sep="\t")
"/NAS/luozh/CRC_conservation/step10_chromatin_state/state_compare/ortholog_state/quiscent_active_enhancer_loci.txt"
"/NAS/luozh/CRC_conservation/step10_chromatin_state/state_compare/ortholog_state/human_quiescent_ActiveEnhancer_loci.txt"
mouse_quiescent_ActiveEnhancer = makeGRangesFromDataFrame(read.csv("/NAS/luozh/CRC_conservation/step10_chromatin_state/state_compare/ortholog_state/mouse_quiescent_ActiveEnhancer_loci.txt", sep="\t"), 
                              keep.extra.columns = F)
score_mouse_quiescent_ActiveEnhancer <- gscores(phast,mouse_quiescent_ActiveEnhancer)
calc_active_score <- function(txt_path){
  df_active = read.csv(txt_path, sep="\t")
  grange_active <- makeGRangesFromDataFrame(df_active, keep.extra.columns = F)
  score_active <- gscores(phast,grange_active)
  return(score_active)
}
mouse_Q_A <- calc_active_score("/NAS/luozh/CRC_conservation/step10_chromatin_state/state_compare/ortholog_state/mouse_quiescent_ActiveEnhancer_loci.txt")
human_Q_A <- calc_active_score("/NAS/luozh/CRC_conservation/step10_chromatin_state/state_compare/ortholog_state/human_quiescent_ActiveEnhancer_loci.txt")
all_Q_A <- calc_active_score("/NAS/luozh/CRC_conservation/step10_chromatin_state/state_compare/ortholog_state/quiscent_active_enhancer_loci.txt")
summary(mouse_Q_A$default)
summary(human_Q_A$default)
summary(all_Q_A$default)





