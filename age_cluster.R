age <- read.csv("/home/shimw/project/conservation/conserve_age_result.csv")
age$cluster = stringr::str_sub(age$cluster,0,9)
age$Human.gene.age = as.numeric(age$Human.gene.age)
age$Human.gene.age[is.na(age$Human.gene.age)] = 4290
table(age[age$cluster=="Cluster 1",]$Human.gene.age)
table(age[age$cluster=="Cluster 2",]$Human.gene.age)
table(age[age$cluster=="Cluster 3",]$Human.gene.age)
table(age[age$cluster=="Cluster 4",]$Human.gene.age)
table(age[age$cluster=="Cluster 5",]$Human.gene.age)



kk = purrr::map(c("Cluster 1","Cluster 2","Cluster 3","Cluster 4","Cluster 5"),function(cc){
  dd = as.data.frame(table(age[age$cluster==cc,]$Human.gene.age))
  names(dd) = c("age","number")
  dd$cluster = cc
  return(dd)
})%>%bind_rows

kk$age = as.numeric(as.vector(kk$age))
ggplot(data = kk, mapping = aes(x = age, y = number, color = cluster)) + geom_line()


table(age[age$Human.gene.age<100,]$cluster)
table(age$cluster)




