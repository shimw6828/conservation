library(GenomicFeatures)

orth = readr::read_tsv("/home/shimw/project/conservation/mm10_hg19.orthologs.top.txt")
orth <- orth[,c("lnc","ortholog")]


mouse_txdb<-makeTxDbFromGFF(file="/home/shimw/project/conservation/gencode.vM27.long_noncoding_RNAs.gtf",format="gtf")
kk = transcriptsBy(mouse_txdb,by=c("gene"))
mouse_idmap = purrr::map(names(kk),function(nn){
  data.frame("gene_id"=nn,"txid"=kk[[nn]]$tx_name)
})%>%dplyr::bind_rows()
row.names(mouse_idmap) <- mouse_idmap$txid
orth$mouse_gene <- stringr::str_split_fixed(mouse_idmap[orth$lnc,"gene_id"],"\\.",2)[,1]

human_txdb<-makeTxDbFromGFF(file="/home/shimw/project/conservation/gencode.v19.long_noncoding_RNAs.gtf",format="gtf")
tt = transcriptsBy(human_txdb,by=c("gene"))
human_idmap = purrr::map(names(tt),function(nn){
  data.frame("gene_id"=nn,"txid"=tt[[nn]]$tx_name)
})%>%dplyr::bind_rows()

row.names(human_idmap) <- human_idmap$txid
orth$human_gene <- stringr::str_split_fixed(human_idmap[orth$ortholog,"gene_id"],"\\.",2)[,1]

stringr::str_split_fixed(human_idmap[orth$ortholog,"gene_id"],"\\.",2)[,1]
readr::write_csv(orth,"/home/shimw/project/conservation/lnc_orth.csv")

