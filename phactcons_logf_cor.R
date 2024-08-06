rm(list = ls())
library(dplyr)

genetotx <- readr::read_csv("/home/shimw/project/conservation/phastcons/genetotx.csv")
df_data = read.csv("/home/panxl/CRC2/all_cluater")
df_data$cluster = stringr::str_sub(df_data$cluster,0,9)
df_data$abs_logF <- abs(df_data$Human_logFC1-df_data$Mouse_logFC1)
# df_data = read.csv("/home/mwshi/project/conservation/all_cluater")
df_data = df_data[,c("Mouse.gene.stable.ID","Gene.stable.ID","abs_logF","cluster")]
names(df_data)[1] = "gene"
df_data = dplyr::left_join(df_data,genetotx)


####phastcons score
exon_scores <- readr::read_csv("/home/shimw/project/conservation/phastcons/exon_scores.csv")
exon_scores = left_join(exon_scores,df_data)

####nacc score
human_pairs_top20_cor = readr::read_csv("/home/shimw/project/conservation/NACC/human_pairs_distance.csv")
mouse_pairs_top20_cor = readr::read_csv("/home/shimw/project/conservation/NACC/mouse_pairs_distance.csv")
mouse_pairs_top20_cor$NACC =  (mouse_pairs_top20_cor$d+human_pairs_top20_cor$d)/2

exon_scores_nacc = left_join(exon_scores,mouse_pairs_top20_cor)

library(ggpointdensity)
cor.test(exon_scores_nacc$NACC, exon_scores_nacc$score)#-0.06552183 p-value 1.882e-14
cols <- densCols(exon_scores_nacc$NACC,exon_scores_nacc$score, colramp=colorRampPalette(c("blue","green","yellow", "red")), nbin = 100)

p1 = ggplot(exon_scores_nacc, aes(x=NACC,y=score)) +
  geom_point(size=0.01,color = cols) +
  # geom_pointdensity()+
  scale_colour_viridis(option = "magma",direction = -1) +
  annotate("text", x = 3, y = 1.1, label = "r= -0.07, p-value 1.9e-14",
           color="#350E20FF",size = 7 )+# x = -3, y = -0.3
  # geom_smooth(method=lm,se = F,color="red")+
  labs(x="NACC",y="Phastcons score",title="")+
  theme_classic(base_line_size = 1) +
  theme(axis.title = element_text(size=28),
        axis.text = element_text(size=24, color="black"),
        title = element_text(size = 30,
                             color = "black",
                             hjust = 0.5),
        axis.text.x = element_text(size = 24,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.text.y = element_text(size = 24,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        legend.title = element_blank(),
        legend.text = element_text(size=12), legend.key.size = unit(0.6,"cm"),
        axis.ticks.length.x = unit(2,"mm"))
ggsave("/home/shimw/project/conservation/NACC/NACC_phastcons.pdf",p1,width=6.5,height = 5.5)






cor.test(exon_scores_nacc$score, exon_scores_nacc$abs_logF)#-0.1613364 p-value < 2.2e-16
p1 = ggplot(exon_scores_nacc, aes(x=score,y=log10(abs_logF))) +
  geom_point(size=1,color = exon_scores_nacc$cols) +
  # scale_color_uchicago() +
  annotate("text", x = 0.5, y = 1.8, label = "r= -0.1613, p-value < 2.2e-16",
           color="#350E20FF",size = 7 )+# x = -3, y = -0.3
  # geom_smooth(method=lm,se = F,color="black",size = 0.5,alpha = 0.5)+
  labs(y="Absolute deviation of fold change",x="Phastcons score",title="")+
  theme_classic(base_line_size = 1) +
  theme(axis.title = element_text(size=28),
        axis.text = element_text(size=24, color="black"),
        title = element_text(size = 30,
                             color = "black",
                             hjust = 0.5),
        axis.text.x = element_text(size = 24,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.text.y = element_text(size = 24,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        legend.title = element_blank(),
        legend.text = element_text(size=12), legend.key.size = unit(0.6,"cm"),
        axis.ticks.length.x = unit(2,"mm"))
ggsave("/home/shimw/project/conservation/NACC/logf_phastcons.pdf",p1,width=7,height = 6.3)
###########################################################################
#################################其他元件的相关性
exon_scores <- readr::read_csv("/home/shimw/project/conservation/phastcons/exon_scores.csv")
exon_scores = left_join(exon_scores,df_data)
exon_scores_nacc = left_join(exon_scores,mouse_pairs_top20_cor)
exon_scores_nacc = exon_scores_nacc[,c("gene","score","abs_logF","cluster","NACC")]
exon_scores_nacc$type = "Exons"
exon_scores_nacc$cols <- densCols(exon_scores_nacc$score,log10(exon_scores_nacc$abs_logF), colramp=colorRampPalette(c("#809fff","green","yellow", "red")), nbin = 100)

tss_scores <- readr::read_csv("/home/shimw/project/conservation/phastcons/tss_scores.csv")
tss_scores = left_join(tss_scores,df_data)
tss_scores_nacc = left_join(tss_scores,mouse_pairs_top20_cor)
tss_scores_nacc = tss_scores_nacc[,c("gene","score","abs_logF","cluster","NACC")]
tss_scores_nacc$type = "Promoter"
cor.test(tss_scores_nacc$score, tss_scores_nacc$abs_logF)#-0.08863388 p-value < 2.2e-16 
tss_scores_nacc$cols <- densCols(tss_scores_nacc$score,log10(tss_scores_nacc$abs_logF), colramp=colorRampPalette(c("#809fff","green","yellow", "red")), nbin = 100)

fiveUTRs_scores <- readr::read_csv("/home/shimw/project/conservation/phastcons/fiveUTRs_scores.csv")
fiveUTRs_scores = left_join(fiveUTRs_scores,df_data)
fiveUTRs_scores_nacc = left_join(fiveUTRs_scores,mouse_pairs_top20_cor)
fiveUTRs_scores_nacc = fiveUTRs_scores_nacc[,c("gene","score","abs_logF","cluster","NACC")]
fiveUTRs_scores_nacc$type = "5′ UTR"
cor.test(fiveUTRs_scores_nacc$score, fiveUTRs_scores_nacc$abs_logF)#-0.08115137 p-value < 2.2e-16
fiveUTRs_scores_nacc$cols <- densCols(fiveUTRs_scores_nacc$score,log10(fiveUTRs_scores_nacc$abs_logF), colramp=colorRampPalette(c("#809fff","green","yellow", "red")), nbin = 100)

threeUTRs_scores <- readr::read_csv("/home/shimw/project/conservation/phastcons/threeUTRs_scores.csv")
threeUTRs_scores = left_join(threeUTRs_scores,df_data)
threeUTRs_scores_nacc = left_join(threeUTRs_scores,mouse_pairs_top20_cor)
threeUTRs_scores_nacc = threeUTRs_scores_nacc[,c("gene","score","abs_logF","cluster","NACC")]
threeUTRs_scores_nacc$type = "3′ UTR"
cor.test(threeUTRs_scores_nacc$score, threeUTRs_scores_nacc$abs_logF)#-0.1249493 p-value < 2.2e-16
threeUTRs_scores_nacc$cols <- densCols(threeUTRs_scores_nacc$score,log10(threeUTRs_scores_nacc$abs_logF), colramp=colorRampPalette(c("#809fff","green","yellow", "red")), nbin = 100)

intron_scores <- readr::read_csv("/home/shimw/project/conservation/phastcons/intron_scores.csv")
intron_scores = left_join(intron_scores,df_data)
intron_scores_nacc = left_join(intron_scores,mouse_pairs_top20_cor)
intron_scores_nacc = intron_scores_nacc[,c("gene","score","abs_logF","cluster","NACC")]
intron_scores_nacc$type = "Introns"
cor.test(intron_scores_nacc$score, intron_scores_nacc$abs_logF)#-0.06667509 p-value = 7.797e-14
intron_scores_nacc$cols <- densCols(intron_scores_nacc$score,log10(intron_scores_nacc$abs_logF), colramp=colorRampPalette(c("#809fff","green","yellow", "red")), nbin = 100)

splicing_score <- readr::read_csv("/home/shimw/project/conservation/phastcons/splicing_score.csv")
splicing_score = left_join(splicing_score,df_data)
splicing_score_nacc = left_join(splicing_score,mouse_pairs_top20_cor)
splicing_score_nacc = splicing_score_nacc[,c("gene","score","abs_logF","cluster","NACC")]
splicing_score_nacc$type = "Splice Sites"
cor.test(splicing_score_nacc$score, splicing_score_nacc$abs_logF)#-0.0846116 p-value < 2.2e-16
splicing_score_nacc$cols <- densCols(splicing_score_nacc$score,log10(splicing_score_nacc$abs_logF), colramp=colorRampPalette(c("#809fff","green","yellow", "red")), nbin = 100)

all_score_logf = list(exon_scores_nacc,tss_scores_nacc,fiveUTRs_scores_nacc,
          threeUTRs_scores_nacc,intron_scores_nacc,splicing_score_nacc)%>%
  dplyr::bind_rows()
cor.test(exon_scores_nacc$score, exon_scores_nacc$score)#-0.1613364 p-value < 2.2e-16
cor.test(tss_scores_nacc$score, tss_scores_nacc$abs_logF)#-0.08863388 p-value < 2.2e-16 
cor.test(fiveUTRs_scores_nacc$score, fiveUTRs_scores_nacc$abs_logF)#-0.08115137 p-value < 2.2e-16
cor.test(threeUTRs_scores_nacc$score, threeUTRs_scores_nacc$abs_logF)#-0.1249493 p-value < 2.2e-16
cor.test(intron_scores_nacc$score, intron_scores_nacc$abs_logF)#-0.06667509 p-value = 7.797e-14
cor.test(splicing_score_nacc$score, splicing_score_nacc$abs_logF)#-0.0846116 p-value < 2.2e-16

all_score_logf$type = factor(all_score_logf$type,
                        levels = c("Exons","Promoter","5′ UTR",
                                   "3′ UTR","Introns","Splice Sites"))


f_labels <- data.frame(type =c("Exons","Promoter","5′ UTR",
                              "3′ UTR","Introns","Splice Sites"),
                       label = c("r= -0.16, p-value < 2.2e-16",
                                 "r= -0.09, p-value < 2.2e-16",
                                 "r= -0.08, p-value < 2.2e-16",
                                 "r= -0.12, p-value < 2.2e-16",
                                 "r= -0.07, p-value = 7.8e-14",
                                 "r= -0.08, p-value < 2.2e-16"))
f_labels$type = factor(f_labels$type,
                       levels = c("Exons","Promoter","5′ UTR",
                                  "3′ UTR","Introns","Splice Sites"))
p1 <- ggplot(data = all_score_logf,aes(x=score,y=log10(abs_logF))) +
  geom_point(size=0.1,color = all_score_logf$cols) +
  # geom_smooth(method=lm,se = F,color="red")+
  geom_text(x = 0.5, y = 1.8, aes(label = label), data = f_labels)+
  theme_classic(base_family = "sans",base_size = 18,base_line_size = 1)+
  scale_y_continuous(expand = expansion(mult=c(0,0.25)))+
  facet_wrap(~ type, ncol=3)+
  labs(y="Absolute deviation of fold change",x="Phastcons score",title="")+
  theme_classic(base_line_size = 1) +
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=12, color="black"),
        title = element_text(size = 18,
                             color = "black",
                             hjust = 0.5),
        axis.text.x = element_text(size = 12,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.text.y = element_text(size = 12,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.ticks.length.x = unit(1,"mm"),
        strip.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        strip.text = element_text(size = 12))

# ggsave("/home/shimw/project/conservation/NACC/all_logf_phastcons_cor.pdf",p1,width=8.2,height = 5.5)
####grDevices这个函数，可以保存这个撇
cairo_pdf("/home/shimw/project/conservation/NACC/all_logf_phastcons_cor.pdf",width=8,height = 5.5)
p1
dev.off()
########
##nacc和score的相关性
#####上面之前已经算过了，直接计算相关性和
cor.test(exon_scores_nacc$NACC, exon_scores_nacc$score)#-0.06552183 p-value = 1.882e-14 
exon_scores_nacc$cols <- densCols(exon_scores_nacc$NACC,exon_scores_nacc$score, colramp=colorRampPalette(c("#809fff","green","yellow", "red")), nbin = 100)

cor.test(tss_scores_nacc$NACC, tss_scores_nacc$score)#-0.05571977 p-value = 1.605e-07 
tss_scores_nacc$cols <- densCols(tss_scores_nacc$NACC,tss_scores_nacc$score, colramp=colorRampPalette(c("#809fff","green","yellow", "red")), nbin = 100)

cor.test(fiveUTRs_scores_nacc$NACC, fiveUTRs_scores_nacc$score)#-0.0582724 p-value = 1.437e-11
fiveUTRs_scores_nacc$cols <- densCols(fiveUTRs_scores_nacc$NACC,fiveUTRs_scores_nacc$score, colramp=colorRampPalette(c("#809fff","green","yellow", "red")), nbin = 100)

cor.test(threeUTRs_scores_nacc$NACC, threeUTRs_scores_nacc$score)#-0.0613889 p-value = 1.067e-11
threeUTRs_scores_nacc$cols <- densCols(threeUTRs_scores_nacc$NACC,threeUTRs_scores_nacc$score, colramp=colorRampPalette(c("#809fff","green","yellow", "red")), nbin = 100)

cor.test(intron_scores_nacc$NACC, intron_scores_nacc$score)#-0.01518848 p-value = 0.089
intron_scores_nacc$cols <- densCols(intron_scores_nacc$NACC,intron_scores_nacc$score, colramp=colorRampPalette(c("#809fff","green","yellow", "red")), nbin = 100)

cor.test(splicing_score_nacc$NACC, splicing_score_nacc$score)#-0.01112599 p-value = 0.1987
splicing_score_nacc$cols <- densCols(splicing_score_nacc$NACC,splicing_score_nacc$score, colramp=colorRampPalette(c("#809fff","green","yellow", "red")), nbin = 100)

all_score_nacc = list(exon_scores_nacc,tss_scores_nacc,fiveUTRs_scores_nacc,
                      threeUTRs_scores_nacc,intron_scores_nacc,splicing_score_nacc)%>%
  dplyr::bind_rows()
cor.test(exon_scores_nacc$NACC, exon_scores_nacc$score)#-0.06552183 p-value = 1.882e-14 
cor.test(tss_scores_nacc$NACC, tss_scores_nacc$score)#-0.05571977 p-value = 1.605e-07 
cor.test(fiveUTRs_scores_nacc$NACC, fiveUTRs_scores_nacc$score)#-0.0582724 p-value = 1.437e-11
cor.test(threeUTRs_scores_nacc$NACC, threeUTRs_scores_nacc$score)#-0.0613889 p-value = 1.067e-11
cor.test(intron_scores_nacc$NACC, intron_scores_nacc$score)#-0.01518848 p-value = 0.089
cor.test(splicing_score_nacc$NACC, splicing_score_nacc$score)#-0.01112599 p-value = 0.1987

all_score_nacc$type = factor(all_score_nacc$type,
                             levels = c("Exons","Promoter","5′ UTR",
                                        "3′ UTR","Introns","Splice Sites"))


f_labels <- data.frame(type =c("Exons","Promoter","5′ UTR",
                               "3′ UTR","Introns","Splice Sites"),
                       label = c("r= -0.07, p-value = 1.9e-14",
                                 "r= -0.06, p-value = 1.6e-7",
                                 "r= -0.06, p-value = 1.4e-11",
                                 "r= -0.06, p-value = 1.1e-11",
                                 "r= -0.02, p-value = 0.09",
                                 "r= -0.01, p-value = 0.20"))
f_labels$type = factor(f_labels$type,
                       levels = c("Exons","Promoter","5′ UTR",
                                  "3′ UTR","Introns","Splice Sites"))
p1 <- ggplot(data = all_score_nacc,aes(x=score,y=NACC)) +
  geom_point(size=0.1,color = all_score_nacc$cols) +
  # geom_smooth(method=lm,se = F,color="red")+
  geom_text(x = 0.5, y = 7, aes(label = label), data = f_labels)+
  theme_classic(base_family = "sans",base_size = 18,base_line_size = 1)+
  scale_y_continuous(expand = expansion(mult=c(0,0.2)))+
  facet_wrap(~ type, ncol=3)+
  labs(y="Phastcons score",x="NACC",title="")+
  theme_classic(base_line_size = 1) +
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=12, color="black"),
        title = element_text(size = 18,
                             color = "black",
                             hjust = 0.5),
        axis.text.x = element_text(size = 12,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.text.y = element_text(size = 12,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.ticks.length.x = unit(1,"mm"),
        strip.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        strip.text = element_text(size = 12))

cairo_pdf("/home/shimw/project/conservation/NACC/all_nacc_phastcons_cor.pdf",width=7.8,height = 5.5)
p1
dev.off()
##################################################################





df_intestine = read.csv("/home/shimw/project/enhancer_map/conserved/colon-specific.csv")
df_species = read.csv("/home/shimw/project/enhancer_map/conserved/species-specific.csv")#744
library("biomaRt")

ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "ensembl.org")
df_species <- getBM(attributes=c('external_gene_name','ensembl_gene_id'), 
      filters = 'external_gene_name', 
      values = df_species$gene_id, 
      mart = ensembl)


exon_scores_nacc
####物种
summary(exon_scores_nacc[exon_scores_nacc$Gene.stable.ID%in%df_species$ensembl_gene_id,]$NACC)
summary(exon_scores_nacc[!exon_scores_nacc$Gene.stable.ID%in%df_species$ensembl_gene_id,]$NACC)
####species特异基因无显著差异
t.test(exon_scores_nacc[exon_scores_nacc$Gene.stable.ID%in%df_species$ensembl_gene_id,]$NACC,
       exon_scores_nacc[!exon_scores_nacc$Gene.stable.ID%in%df_species$ensembl_gene_id,]$NACC)

####肠道
summary(exon_scores_nacc[exon_scores_nacc$gene%in%df_intestine$Gene,]$NACC)
summary(exon_scores_nacc[!exon_scores_nacc$gene%in%df_intestine$Gene,]$NACC)
####肠道特异基因无显著差异
t.test(exon_scores_nacc[exon_scores_nacc$gene%in%df_intestine$Gene,]$NACC,
       exon_scores_nacc[!exon_scores_nacc$gene%in%df_intestine$Gene,]$NACC)


#### 管家基因
library("biomaRt")
HK_genes <- readr::read_tsv("/home/shimw/project/conservation/HK_genes.txt",col_names = c("symbol","ncbi_id"))
ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "grch37.ensembl.org")
symboltoid = getBM(attributes=c('ensembl_gene_id', 'external_gene_name'), 
                   filters = 'external_gene_name', 
                   values = HK_genes$symbol, 
                   mart = ensembl)
exon_scores_nacc[exon_scores_nacc$Gene.stable.ID%in%symboltoid$ensembl_gene_id,]
summary(exon_scores_nacc[exon_scores_nacc$Gene.stable.ID%in%symboltoid$ensembl_gene_id,]$NACC)
summary(exon_scores_nacc[!exon_scores_nacc$Gene.stable.ID%in%symboltoid$ensembl_gene_id,]$NACC)
t.test(exon_scores_nacc[exon_scores_nacc$Gene.stable.ID%in%symboltoid$ensembl_gene_id,]$NACC,
                        exon_scores_nacc[!exon_scores_nacc$Gene.stable.ID%in%symboltoid$ensembl_gene_id,]$NACC)



exon_scores_nacc$is_house_keeping <- if_else(exon_scores_nacc$Gene.stable.ID%in%symboltoid$ensembl_gene_id,"House keeping","General")
my_comparisons <- list(c("House keeping", "General"))
p1 <- ggplot(data = na.omit(exon_scores_nacc),aes(x = is_house_keeping,y = NACC, fill = is_house_keeping)) +
  geom_violin(adjust = .5) +geom_boxplot(width=0.1, fill="white")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif", tip.length=0.02,method.args = list(alternative = "less"), method = "t.test") +
  theme_classic(base_family = "sans",base_size = 18,base_line_size = 1.1)+
  labs(title = NULL, y="NACC score", x=NULL, fill=NULL) +
  scale_fill_manual(values=c("House keeping" = "#CB9C7A", "General" = "#8696a7"))+
  scale_x_discrete(labels = c("House keeping" = "House keeping", "General" = "General"))+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5,size=16),
        axis.title=element_text(size=18),
        axis.text=element_text(size=15))
ggsave("/home/shimw/project/conservation/NACC/nacc_housekeeping.pdf",p1,width=3.7,height = 4)


exon_scores_nacc$is_specific <- if_else(exon_scores_nacc$Gene.stable.ID%in%df_species$ensembl_gene_id,"Species specific","General")
my_comparisons <- list(c("Species specific", "General"))
p2 <- ggplot(data = na.omit(exon_scores_nacc),aes(x = is_specific,y = NACC, fill = is_specific)) +
  geom_violin(adjust = .5) +geom_boxplot(width=0.1, fill="white")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif", tip.length=0.02, method = "t.test") +
  theme_classic(base_family = "sans",base_size = 18,base_line_size = 1.1)+
  labs(title = NULL, y=NULL, x=NULL, fill=NULL) +
  scale_fill_manual(values=c("Species specific" = "#CB9C7A", "General" = "#8696a7"))+
  scale_x_discrete(labels = c("Species specific" = "Species specific", "General" = "General"))+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5,size=16),
        axis.title=element_text(size=18),
        axis.text=element_text(size=15))


exon_scores_nacc$is_Intestine <- if_else(exon_scores_nacc$gene%in%df_intestine$Gene,"Intestine specific","General")
my_comparisons <- list(c("Intestine specific", "General"))
p3 <- ggplot(data = na.omit(exon_scores_nacc),aes(x = is_Intestine,y = NACC, fill = is_Intestine)) +
  geom_violin(adjust = .5) +geom_boxplot(width=0.1, fill="white")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif", tip.length=0.02, method = "t.test") +
  theme_classic(base_family = "sans",base_size = 18,base_line_size = 1.1)+
  labs(title = NULL, y=NULL, x=NULL, fill=NULL) +
  scale_fill_manual(values=c("Intestine specific" = "#CB9C7A", "General" = "#8696a7"))+
  scale_x_discrete(labels = c("Intestine specific" = "Intestine specific", "General" = "General"))+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5,size=16),
        axis.title=element_text(size=18),
        axis.text=element_text(size=15))
library("cowplot")

p1+p2+p3

ggdraw() +
  draw_plot(p1, x = 0, y = .5, width = .5, height = .5) +
  draw_plot(p2, x = .3, y = .5, width = .5, height = .5) +
  draw_plot(p3, x = 0.6, y = 0.5, width = 0.5, height = 0.5)



Cairo::CairoPDF("/home/shimw/project/conservation/NACC/nacc_ish_all.pdf",width=11,height = 3.7)
plot_grid(p1,p2,p3,ncol =3 )
dev.off()



