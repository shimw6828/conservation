rm(list = ls()); gc()
library(biomaRt)
library(dplyr)
library(tibble)


library(GenomicFeatures)
library(GenomicRanges)


txdb = makeTxDbFromGFF("/NAS/luozh/CRC_conservation/reference/gencode.v32lift37.annotation.gtf", format="gtf")
get_genes_resion = function(txdb, upstream, downstream, output_file){
  gene <- genes(txdb)
  gene <- sort(gene)
  #gene = gene[seqnames(gene) != "chrM",]
  df <- data.frame(chr=seqnames(gene),
                   start=start(gene) - upstream,
                   end=end(gene) + downstream,
                   strands=strand(gene),
                   names=names(gene))
  df$names = gsub("\\..*","",df$names)
  return(df)
  # write.table(df, file=output_file, quote=F, sep="\t", row.names=F, col.names=F)
}

human_gene_position = get_genes_resion(txdb = txdb, upstream = 3000, downstream = 3000)

ensembl37 <- useMart(host="https://grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
## 用于获取绘制track图时基因的坐标，上下游扩张10kb
track_df <- readr::read_csv("~/project/conservation/track/state_change_gene.csv")

attributes <- c("ensembl_gene_id", "external_gene_name")
results <- getBM(attributes=attributes, filters="ensembl_gene_id", values=track_df$gene, mart=ensembl37)

track_df_3k <- left_join(track_df,results, by = c("gene" = "ensembl_gene_id"))%>%
  left_join(human_gene_position, by=c("gene"="names"))

readr::write_csv(track_df_3k, "~/project/conservation/track/state_change_gene_3k.csv")
## 获取KEGG注释，方便之后挑选基因
library(KEGGREST)


get_pathways <- function(gene_symbol) {
  result <- tryCatch({
    kegg_info <- keggGet(paste0("hsa:", gene_symbol))
    if (!is.null(kegg_info[[1]]$PATHWAY)) {
      tibble("external_gene_name"=gene_symbol, "pathway" = kegg_info[[1]]$PATHWAY)%>%
      return()
    } else {
      tibble("external_gene_name"=gene_symbol, "pathway" = "")%>%
      return()
    }
  }, 
  error = function(e) {
    message(paste("Error fetching data for", gene_symbol, ": ", e$message))
    tibble("external_gene_name"=gene_symbol, "pathway" = "")%>%
    return()
  })
  return(result)
}

pathways_list <- lapply(track_df_3k$external_gene_name, get_pathways)%>%
  bind_rows()

track_df_kegg <- left_join(track_df_3k,pathways_list)
readr::write_tsv(track_df_kegg, "~/project/conservation/track/state_change_gene_kegg.tsv")





## mouse##############
txdb = makeTxDbFromGFF("/home/shimw/project/conservation/gencode.vM20.basic.annotation.gtf", format="gtf")
get_genes_resion = function(txdb, upstream, downstream, output_file){
  gene <- genes(txdb)
  gene <- sort(gene)
  #gene = gene[seqnames(gene) != "chrM",]
  df <- data.frame(chr=seqnames(gene),
                   start=start(gene) - upstream,
                   end=end(gene) + downstream,
                   strands=strand(gene),
                   names=names(gene))
  df$names = gsub("\\..*","",df$names)
  return(df)
  # write.table(df, file=output_file, quote=F, sep="\t", row.names=F, col.names=F)
}

mouse_gene_position = get_genes_resion(txdb = txdb, upstream = 3000, downstream = 3000)

ensembl37 <- useMart(host="https://nov2020.archive.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
## 用于获取绘制track图时基因的坐标，上下游扩张10kb
track_df <- readr::read_csv("~/project/conservation/track/mouse_state_change_gene.csv")

attributes <- c("ensembl_gene_id", "external_gene_name")
results <- getBM(attributes=attributes, filters="ensembl_gene_id", values=track_df$gene, mart=ensembl37)

track_df_3k <- left_join(track_df,results, by = c("gene" = "ensembl_gene_id"))%>%
  left_join(gene_position, by=c("gene"="names"))

readr::write_csv(track_df_3k, "~/project/conservation/track/mouse_state_change_gene_3k.csv")
## 获取KEGG注释，方便之后挑选基因
library(KEGGREST)


get_pathways <- function(gene_symbol) {
  result <- tryCatch({
    kegg_info <- keggGet(paste0("mmu:", gene_symbol))
    if (!is.null(kegg_info[[1]]$PATHWAY)) {
      tibble("external_gene_name"=gene_symbol, "pathway" = kegg_info[[1]]$PATHWAY)%>%
        return()
    } else {
      tibble("external_gene_name"=gene_symbol, "pathway" = "")%>%
        return()
    }
  }, 
  error = function(e) {
    message(paste("Error fetching data for", gene_symbol, ": ", e$message))
    tibble("external_gene_name"=gene_symbol, "pathway" = "")%>%
      return()
  })
  return(result)
}

pathways_list <- lapply(track_df_3k$external_gene_name, get_pathways)%>%
  bind_rows()

track_df_kegg <- left_join(track_df_3k,pathways_list)
readr::write_tsv(track_df_kegg, "~/project/conservation/track/mouse_state_change_gene_kegg.tsv")


## expression track ##############
## human
## 从上面获取基因的位置信息
cluster_PeaktoGene <- readr::read_csv("/home/shimw/project/conservation/enhancer_count/human_enhancer_DE_annotated.csv")%>%
  dplyr::arrange(h_padj)

unique(cluster_PeaktoGene[1:200,]$h_gene)
ensembl37 <- useMart(host="https://asia.ensembl.org/", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
attributes <- c("ensembl_gene_id", "external_gene_name")
results <- getBM(attributes=attributes, filters="ensembl_gene_id", values=unique(cluster_PeaktoGene[1:200,]$h_gene), mart=ensembl37)
track_df_3k <- cluster_PeaktoGene[1:200,]%>%
  select("h_gene", "cluster")%>%
  rename("gene"=h_gene)%>%
  left_join(results, by = c("gene" = "ensembl_gene_id"))%>%
  left_join(gene_position, by=c("gene"="names"))%>%
  distinct()
readr::write_csv(track_df_3k, "~/project/conservation/track/expression_change_gene_3k.csv")

library(KEGGREST)


get_pathways <- function(gene_symbol) {
  result <- tryCatch({
    kegg_info <- keggGet(paste0("hsa:", gene_symbol))
    if (!is.null(kegg_info[[1]]$PATHWAY)) {
      tibble("external_gene_name"=gene_symbol, "pathway" = kegg_info[[1]]$PATHWAY)%>%
        return()
    } else {
      tibble("external_gene_name"=gene_symbol, "pathway" = "")%>%
        return()
    }
  }, 
  error = function(e) {
    message(paste("Error fetching data for", gene_symbol, ": ", e$message))
    tibble("external_gene_name"=gene_symbol, "pathway" = "")%>%
      return()
  })
  return(result)
}
pathways_list <- lapply(track_df_3k$external_gene_name, get_pathways)%>%
  bind_rows()

track_df_kegg <- left_join(track_df_3k,pathways_list)
readr::write_tsv(track_df_kegg, "~/project/conservation/track/expression_change_gene_kegg.tsv")

## mouse
cluster_PeaktoGene <- readr::read_csv("/home/shimw/project/conservation/enhancer_count/mouse_enhancer_DE_annotated.csv")%>%
  dplyr::arrange(m_padj)
unique(cluster_PeaktoGene[1:200,]$m_gene)
ensembl37 <- useMart(host="https://nov2020.archive.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")

attributes <- c("ensembl_gene_id", "external_gene_name")
results <- getBM(attributes=attributes, filters="ensembl_gene_id", values=unique(cluster_PeaktoGene[1:200,]$m_gene), mart=ensembl37)
track_df_3k <- cluster_PeaktoGene[1:200,]%>%
  select("m_gene", "cluster")%>%
  rename("gene"=m_gene)%>%
  left_join(results, by = c("gene" = "ensembl_gene_id"))%>%
  left_join(gene_position, by=c("gene"="names"))%>%
  distinct()
readr::write_csv(track_df_3k, "~/project/conservation/track/mouse_expression_change_gene_3k.csv")
get_pathways <- function(gene_symbol) {
  result <- tryCatch({
    kegg_info <- keggGet(paste0("mmu:", gene_symbol))
    if (!is.null(kegg_info[[1]]$PATHWAY)) {
      tibble("external_gene_name"=gene_symbol, "pathway" = kegg_info[[1]]$PATHWAY)%>%
        return()
    } else {
      tibble("external_gene_name"=gene_symbol, "pathway" = "")%>%
        return()
    }
  }, 
  error = function(e) {
    message(paste("Error fetching data for", gene_symbol, ": ", e$message))
    tibble("external_gene_name"=gene_symbol, "pathway" = "")%>%
      return()
  })
  return(result)
}
pathways_list <- lapply(track_df_3k$external_gene_name, get_pathways)%>%
  bind_rows()

track_df_kegg <- left_join(track_df_3k,pathways_list)
readr::write_tsv(track_df_kegg, "~/project/conservation/track/mouse_expression_change_gene_kegg.tsv")

############################ Chemokine signaling pathway in mouse ############################
load("/home/panxl/CRC2/hom_kegg.RData")
hom=read.csv("/home/panxl/CRC2/all_cluater")
hom_kegg=hom_kegg@result
l=hom_kegg[hom_kegg$ID=="hsa04062",]$geneID
ll=data.frame(ENTREZID=unlist(strsplit(l,split = "/")))
ENSEMBL_id =bitr(ll$ENTREZID,
                 fromType="ENTREZID",
                 toType= "ENSEMBL",
                 OrgDb = "org.Hs.eg.db", drop = TRUE)
gene_id=hom[hom$Gene.stable.ID%in%ENSEMBL_id$ENSEMBL,]
gene_id$Gene.stable.ID
gene_id$Mouse.gene.stable.ID

ensembl37 <- useMart(host="https://grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
ensembl10 <- useMart(host="https://nov2020.archive.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
attributes <- c("ensembl_gene_id", "external_gene_name")
#mouse_gene_position 从上面的代码中获得
#human_gene_position 
results <- getBM(attributes=attributes, filters="ensembl_gene_id", values=unique(gene_id$Gene.stable.ID), mart=ensembl37)%>%
  dplyr::rename("gene"=ensembl_gene_id)
human_Chemokine_gene <- left_join(results,human_gene_position, by = c("gene" = "names") )
readr::write_csv(human_Chemokine_gene, "~/project/conservation/track/human_Chemokine_gene_3k.csv")

results <- getBM(attributes=attributes, filters="ensembl_gene_id", values=unique(gene_id$Mouse.gene.stable.ID), mart=ensembl10)%>%
  dplyr::rename("gene"=ensembl_gene_id)
mouse_Chemokine_gene <- left_join(results,mouse_gene_position, by = c("gene" = "names") )
readr::write_csv(mouse_Chemokine_gene, "~/project/conservation/track/mouse_Chemokine_gene_3k.csv")







