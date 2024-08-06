target_disease <- readr::read_tsv("/home/shimw/project/conservation/P1-06-Target_disease.txt",
                                  col_names = c("TARGETID","B","INDICATI"))
length(unique(target_disease$TARGETID))
target_disease = target_disease[grepl("colorectal cancer",target_disease$INDICATI, ignore.case = T),]
length(unique(target_disease$TARGETID))

success_target = readr::read_tsv("/home/shimw/project/conservation/P2-02-TTD_uniprot_successful.txt",
                                 col_names = c("A","B","C"))

success_target = tibble("TARGETID"=success_target[success_target$B=="TARGETID",]$C,
       "UNIPROID"=success_target[success_target$B=="UNIPROID",]$C,
       "TARGNAME"=success_target[success_target$B=="TARGNAME",]$C,
       "TARGTYPE"=success_target[success_target$B=="TARGTYPE",]$C)

success_target_re = success_target%>%
  tidyr::separate_rows(UNIPROID,sep = "; ")%>%
  tidyr::separate_rows(UNIPROID,sep = "-")


readr::write_csv(success_target_re[,"UNIPROID"],"/home/shimw/project/conservation/success_target_map.csv")


###clinical
clinical_target = readr::read_tsv("/home/shimw/project/conservation/P2-03-TTD_uniprot_clinical.txt",
                                 col_names = c("A","B","C"))

clinical_target = tibble("TARGETID"=clinical_target[clinical_target$B=="TARGETID",]$C,
                        "UNIPROID"=clinical_target[clinical_target$B=="UNIPROID",]$C,
                        "TARGNAME"=clinical_target[clinical_target$B=="TARGNAME",]$C,
                        "TARGTYPE"=clinical_target[clinical_target$B=="TARGTYPE",]$C)
clinical_target_re = clinical_target%>%
  tidyr::separate_rows(UNIPROID,sep = "; ")%>%
  tidyr::separate_rows(UNIPROID,sep = "-")
readr::write_csv(clinical_target[,"UNIPROID"],"/home/shimw/project/conservation/clinical_target_map.csv")

dim(ttt)
###https://www.uniprot.org/uploadlists/
##手动转一下

##success
success_target = readr::read_tsv("/home/shimw/project/conservation/success_map.txt")
success_target = success_target[success_target$From%in%success_target_re[success_target_re$TARGETID%in%target_disease$TARGETID,]$UNIPROID,]
df_data = read.csv("/home/panxl/CRC2/all_cluater")
table(df_data[df_data$Gene.stable.ID%in%success_target$To,]$cluster)
df_data$cluster = stringr::str_sub(df_data$cluster,0,9)

cgroups = c("Cluster 1","Cluster 2","Cluster 3","Cluster 4","Cluster 5")
odds_ratios <- purrr::map(cgroups,function(cgroup){
  genes_fisher_cluster = matrix(c(nrow(df_data[df_data$cluster==cgroup&df_data$Gene.stable.ID%in%success_target$To,]),
                                  nrow(df_data[df_data$cluster==cgroup&!(df_data$Gene.stable.ID%in%success_target$To),]),
                                  nrow(df_data[df_data$cluster!=cgroup&df_data$Gene.stable.ID%in%success_target$To,]),
                                  nrow(df_data[df_data$cluster!=cgroup&!(df_data$Gene.stable.ID%in%success_target$To),])),
                                nrow = 2,
                                dimnames =list(c("success_target", "non_success_target"),
                                               c("cluster1", "Not cluster1"))
                                
  )
  OR_1 = fisher.test(genes_fisher_cluster, conf.level = 0.95)
  return(tibble("boxLabels"=stringr::str_replace(cgroup," ",""),"boxOdds"=OR_1$estimate,"boxCILow"=OR_1$conf.int[1],"boxCIHigh"=OR_1$conf.int[2]))
})%>%dplyr::bind_rows()

p <- ggplot(odds_ratios, aes(x = boxOdds, y = boxLabels)) + 
  geom_vline(aes(xintercept = 1.0), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = 1, height = 
                   .15, color = c("#CB9C7A", "#8696a7", "#CDB97D", "#7b8b6f", "#A59B95")) +
  geom_point(size = 2, color = c("#CB9C7A", "#8696a7", "#CDB97D", "#7b8b6f", "#A59B95")) +
  scale_color_uchicago() +
  xlab("Odds ratio") +
  ylab("Gene's clusters") + 
  theme(axis.title = element_text(size=24),
        axis.title.y = element_blank(),
        #axis.text = element_text(size=24, color="black"),
        axis.text.x = element_text(size = 20,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.text.y = element_text(size = 20,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        legend.title = element_blank(),
        legend.text = element_text(size=12), legend.key.size = unit(0.6,"cm"),
        axis.ticks.length.x = unit(2,"mm"),
        axis.ticks.length.y = unit(2,"mm")
  )

##clinical
clinical_target = readr::read_tsv("/home/shimw/project/conservation/clinical_map.txt")
clinical_target = clinical_target[clinical_target$From%in%clinical_target_re[clinical_target_re$TARGETID%in%target_disease$TARGETID,]$UNIPROID,]


df_data = read.csv("/home/panxl/CRC2/all_cluater")
table(df_data[df_data$Gene.stable.ID%in%clinical_target$To,]$cluster)
df_data$cluster = stringr::str_sub(df_data$cluster,0,9)

cgroups = c("Cluster 1","Cluster 2","Cluster 3","Cluster 4","Cluster 5")
odds_ratios <- purrr::map(cgroups,function(cgroup){
  genes_fisher_cluster = matrix(c(nrow(df_data[df_data$cluster==cgroup&df_data$Gene.stable.ID%in%clinical_target$To,]),
                                  nrow(df_data[df_data$cluster==cgroup&!(df_data$Gene.stable.ID%in%clinical_target$To),]),
                                  nrow(df_data[df_data$cluster!=cgroup&df_data$Gene.stable.ID%in%clinical_target$To,]),
                                  nrow(df_data[df_data$cluster!=cgroup&!(df_data$Gene.stable.ID%in%clinical_target$To),])),
                                nrow = 2,
                                dimnames =list(c("clinical_target", "non_clinical_target"),
                                               c("cluster1", "Not cluster1"))
                                
  )
  OR_1 = fisher.test(genes_fisher_cluster, conf.level = 0.95)
  return(tibble("boxLabels"=stringr::str_replace(cgroup," ",""),"boxOdds"=OR_1$estimate,"boxCILow"=OR_1$conf.int[1],"boxCIHigh"=OR_1$conf.int[2]))
})%>%dplyr::bind_rows()

p <- ggplot(odds_ratios, aes(x = boxOdds, y = boxLabels)) + 
  geom_vline(aes(xintercept = 1.0), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = 1, height = 
                   .15, color = c("#CB9C7A", "#8696a7", "#CDB97D", "#7b8b6f", "#A59B95")) +
  geom_point(size = 2, color = c("#CB9C7A", "#8696a7", "#CDB97D", "#7b8b6f", "#A59B95")) +
  scale_color_uchicago() +
  xlab("Odds ratio") +
  ylab("Gene's clusters") + 
  theme(axis.title = element_text(size=24),
        axis.title.y = element_blank(),
        #axis.text = element_text(size=24, color="black"),
        axis.text.x = element_text(size = 20,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.text.y = element_text(size = 20,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        legend.title = element_blank(),
        legend.text = element_text(size=12), legend.key.size = unit(0.6,"cm"),
        axis.ticks.length.x = unit(2,"mm"),
        axis.ticks.length.y = unit(2,"mm")
  )


table(df_data[df_data$Gene.stable.ID%in%success_target$To,]$cluster)/(table(df_data[df_data$Gene.stable.ID%in%clinical_target$To,]$cluster)+table(df_data[df_data$Gene.stable.ID%in%success_target$To,]$cluster))
table(df_data[df_data$Gene.stable.ID%in%clinical_target$To,]$cluster)
table(df_data[df_data$Gene.stable.ID%in%success_target$To,]$cluster)

table(success_target$To%in%clinical_target$To)


############################
target_disease <- readr::read_tsv("/home/shimw/project/conservation/P1-06-Target_disease.txt",
                                  col_names = c("TARGETID","B","INDICATI"))
length(unique(target_disease$TARGETID))
target_disease = target_disease[grepl("colorectal cancer",target_disease$INDICATI, ignore.case = T),]
length(unique(target_disease$TARGETID))

success_target = readr::read_tsv("/home/shimw/project/conservation/P2-02-TTD_uniprot_successful.txt",
                                 col_names = c("A","B","C"))
success_target = tibble("TARGETID"=success_target[success_target$B=="TARGETID",]$C,
                        "UNIPROID"=success_target[success_target$B=="UNIPROID",]$C,
                        "TARGNAME"=success_target[success_target$B=="TARGNAME",]$C,
                        "TARGTYPE"=success_target[success_target$B=="TARGTYPE",]$C)
success_target_re = success_target%>%
  tidyr::separate_rows(UNIPROID,sep = "; ")%>%
  tidyr::separate_rows(UNIPROID,sep = "-")
readr::write_csv(success_target_re[,"UNIPROID"],"/home/shimw/project/conservation/success_target_map.csv")

clinical_target = readr::read_tsv("/home/shimw/project/conservation/P2-03-TTD_uniprot_clinical.txt",
                                  col_names = c("A","B","C"))
clinical_target = tibble("TARGETID"=clinical_target[clinical_target$B=="TARGETID",]$C,
                         "UNIPROID"=clinical_target[clinical_target$B=="UNIPROID",]$C,
                         "TARGNAME"=clinical_target[clinical_target$B=="TARGNAME",]$C,
                         "TARGTYPE"=clinical_target[clinical_target$B=="TARGTYPE",]$C)
clinical_target_re = clinical_target%>%
  tidyr::separate_rows(UNIPROID,sep = "; ")%>%
  tidyr::separate_rows(UNIPROID,sep = "-")
readr::write_csv(clinical_target[,"UNIPROID"],"/home/shimw/project/conservation/clinical_target_map.csv")


success_target = readr::read_tsv("/home/shimw/project/conservation/success_map.txt")
success_target = success_target[success_target$From%in%success_target_re[success_target_re$TARGETID%in%target_disease$TARGETID,]$UNIPROID,]
clinical_target = readr::read_tsv("/home/shimw/project/conservation/clinical_map.txt")
clinical_target = clinical_target[clinical_target$From%in%clinical_target_re[clinical_target_re$TARGETID%in%target_disease$TARGETID,]$UNIPROID,]

table(df_data[df_data$Gene.stable.ID%in%success_target$To,]$cluster)/(table(df_data[df_data$Gene.stable.ID%in%clinical_target$To,]$cluster)+table(df_data[df_data$Gene.stable.ID%in%success_target$To,]$cluster))
table(df_data[df_data$Gene.stable.ID%in%clinical_target$To,]$cluster)
table(df_data[df_data$Gene.stable.ID%in%success_target$To,]$cluster)


success_fisher = matrix(c(17,9,11,10),
                              nrow = 2,
                              dimnames =list(c("success", "non_success"),
                                             c("cluster13", "cluster24")))
OR_1 = fisher.test(success_fisher, conf.level = 0.95)


