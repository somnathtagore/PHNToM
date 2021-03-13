
# #egfr samples
# tcga.tumor<-read.table(file="~/Documents/viper-docs/benchmarking_1/signature_with_gtex/batch_corrected_data/data/batch-correct/tcga-tumor/luad-tcga-tumor.txt", header=TRUE)
# #colnames(tcga.tumor)<-gsub('.{10}$','',colnames(tcga.tumor))
# colnames(tcga.tumor)
# dim(tcga.tumor)
# colnames(tcga.tumor)<-gsub('.{16}$','',colnames(tcga.tumor))
# 
# egfr.samples<- 
#   tcga.tumor[,intersect(EGFR.Tumor_Sample_Barcode,colnames(tcga.tumor))]
# dim(egfr.samples) #78

#establish EGFR GOF protein activity signature
#high confidence EGFR-mutated LUAD vs GTEX colon
#L858R and Exon19 deletions
#bootstrap msviper 


brca.tcga.tumor<-read.table(file="~/Documents/viper-docs/benchmarking_1/signature_with_gtex/batch_corrected_data/data/batch-correct/tcga-tumor/brca-tcga-tumor.txt", header=TRUE)
#colnames(tcga.tumor)<-gsub('.{10}$','',colnames(tcga.tumor))
colnames(brca.tcga.tumor)
dim(brca.tcga.tumor)
colnames(brca.tcga.tumor)<-gsub('.{16}$','',colnames(brca.tcga.tumor))

####NEW NEW NEW
# #K111E
# K111E.erbb2.gene.exp.brca <- 
#   brca.tcga.tumor[,intersect(colnames(brca.tcga.tumor),K111E.erbb2.brca.samples)]
# # brca.tcga.tumor[,intersect(colnames(brca.tcga.tumor),G118D.samples)]
# dim(K111E.erbb2.gene.exp.brca) #1

#L755S
L755S.erbb2.gene.exp.brca <- 
  brca.tcga.tumor[,intersect(colnames(brca.tcga.tumor),colnames(L755S.erbb2.brca.samples.erbb2))]
# brca.tcga.tumor[,intersect(colnames(brca.tcga.tumor),G118D.samples)]
dim(L755S.erbb2.gene.exp.brca) #15

#V777L
V777L.erbb2.gene.exp.brca <- 
  brca.tcga.tumor[,intersect(colnames(brca.tcga.tumor),colnames(V777L.erbb2.brca.samples.erbb2))]
# brca.tcga.tumor[,intersect(colnames(brca.tcga.tumor),G118D.samples)]
dim(V777L.erbb2.gene.exp.brca) #10


#D769Y
D769Y.erbb2.gene.exp.brca <- 
  brca.tcga.tumor[,intersect(colnames(brca.tcga.tumor),colnames(D769Y.erbb2.brca.samples.erbb2))]
# brca.tcga.tumor[,intersect(colnames(brca.tcga.tumor),G118D.samples)]
dim(D769Y.erbb2.gene.exp.brca) #10

current.stouffer.only.L755S.erbb2.gene.exp.brca <- apply(L755S.erbb2.gene.exp.brca,1,function(x){
  y <- sum(x)/sqrt(length(x))
  #   #y <- sum(w*x)/sqrt(sum(w^2))
  return(y)
})
current.stouffer.only.V777L.erbb2.gene.exp.brca <- apply(V777L.erbb2.gene.exp.brca,1,function(x){
  y <- sum(x)/sqrt(length(x))
  #   #y <- sum(w*x)/sqrt(sum(w^2))
  return(y)
})
current.stouffer.only.D769Y.erbb2.gene.exp.brca <- apply(D769Y.erbb2.gene.exp.brca,1,function(x){
  y <- sum(x)/sqrt(length(x))
  #   #y <- sum(w*x)/sqrt(sum(w^2))
  return(y)
})


#all lof
all.high.conf.lof.erbb2.gene.exp.brca <- cbind(
  current.stouffer.only.L755S.erbb2.gene.exp.brca,
  current.stouffer.only.V777L.erbb2.gene.exp.brca,
  current.stouffer.only.D769Y.erbb2.gene.exp.brca
)

colnames(all.high.conf.lof.erbb2.gene.exp.brca)<-
  cbind(
    'L755S',
    'V777L',
    'D769Y'
    
  )

colnames(all.high.conf.lof.erbb2.gene.exp.brca)
dim(all.high.conf.lof.erbb2.gene.exp.brca)

# gtex data
# gtex data
colon.gtex.tissue <- read.table("~/Documents/viper-docs/benchmarking_1/signature_with_gtex/batch_corrected_data/data/batch-correct/gtex/colon-gtex.txt", header=TRUE)
dim(colon.gtex.tissue)

#check for identical genes
#all lof
identical(rownames(all.high.conf.lof.erbb2.gene.exp.brca),rownames(colon.gtex.tissue))
table(rownames(all.high.conf.lof.erbb2.gene.exp.brca)%in%rownames(colon.gtex.tissue))
common_genes <- intersect(rownames(all.high.conf.lof.erbb2.gene.exp.brca),rownames(colon.gtex.tissue))
length(common_genes)

all.high.conf.lof.erbb2.gene.exp.brca <- all.high.conf.lof.erbb2.gene.exp.brca[match(common_genes,rownames(all.high.conf.lof.erbb2.gene.exp.brca)),]
colon.gtex.tissue <- colon.gtex.tissue[match(common_genes,rownames(colon.gtex.tissue)),]
identical(rownames(colon.gtex.tissue),rownames(all.high.conf.lof.erbb2.gene.exp.brca))

# #G118D
# identical(rownames(only.G118D.no.low.amp.erbb2.gene.exp.brca),rownames(breast.gtex.tissue))
# table(rownames(only.G118D.no.low.amp.erbb2.gene.exp.brca)%in%rownames(breast.gtex.tissue))
# common_genes <- intersect(rownames(only.G118D.no.low.amp.erbb2.gene.exp.brca),rownames(breast.gtex.tissue))
# length(common_genes)
# 
# only.G118D.no.low.amp.erbb2.gene.exp.brca <- only.G118D.no.low.amp.erbb2.gene.exp.brca[match(common_genes,rownames(only.G118D.no.low.amp.erbb2.gene.exp.brca)),]
# breast.gtex.tissue <- breast.gtex.tissue[match(common_genes,rownames(breast.gtex.tissue)),]
# identical(rownames(breast.gtex.tissue),rownames(only.G118D.no.low.amp.erbb2.gene.exp.brca))
# 
# #N345K
# identical(rownames(only.N345K.no.low.amp.erbb2.gene.exp.brca),rownames(breast.gtex.tissue))
# table(rownames(only.N345K.no.low.amp.erbb2.gene.exp.brca)%in%rownames(breast.gtex.tissue))
# common_genes <- intersect(rownames(only.N345K.no.low.amp.erbb2.gene.exp.brca),rownames(breast.gtex.tissue))
# length(common_genes)
# 
# only.N345K.no.low.amp.erbb2.gene.exp.brca <- only.N345K.no.low.amp.erbb2.gene.exp.brca[match(common_genes,rownames(only.N345K.no.low.amp.erbb2.gene.exp.brca)),]
# breast.gtex.tissue <- breast.gtex.tissue[match(common_genes,rownames(breast.gtex.tissue)),]
# identical(rownames(breast.gtex.tissue),rownames(only.N345K.no.low.amp.erbb2.gene.exp.brca))
# 
# #H1047L
# identical(rownames(only.H1047L.no.low.amp.erbb2.gene.exp.brca),rownames(breast.gtex.tissue))
# table(rownames(only.H1047L.no.low.amp.erbb2.gene.exp.brca)%in%rownames(breast.gtex.tissue))
# common_genes <- intersect(rownames(only.H1047L.no.low.amp.erbb2.gene.exp.brca),rownames(breast.gtex.tissue))
# length(common_genes)
# 
# only.H1047L.no.low.amp.erbb2.gene.exp.brca <- only.H1047L.no.low.amp.erbb2.gene.exp.brca[match(common_genes,rownames(only.H1047L.no.low.amp.erbb2.gene.exp.brca)),]
# breast.gtex.tissue <- breast.gtex.tissue[match(common_genes,rownames(breast.gtex.tissue)),]
# identical(rownames(breast.gtex.tissue),rownames(only.H1047L.no.low.amp.erbb2.gene.exp.brca))
# 
# #H1047R
# identical(rownames(only.H1047R.no.low.amp.erbb2.gene.exp.brca),rownames(breast.gtex.tissue))
# table(rownames(only.H1047R.no.low.amp.erbb2.gene.exp.brca)%in%rownames(breast.gtex.tissue))
# common_genes <- intersect(rownames(only.H1047R.no.low.amp.erbb2.gene.exp.brca),rownames(breast.gtex.tissue))
# length(common_genes)
# 
# only.H1047R.no.low.amp.erbb2.gene.exp.brca <- only.H1047R.no.low.amp.erbb2.gene.exp.brca[match(common_genes,rownames(only.H1047R.no.low.amp.erbb2.gene.exp.brca)),]
# breast.gtex.tissue <- breast.gtex.tissue[match(common_genes,rownames(breast.gtex.tissue)),]
# identical(rownames(breast.gtex.tissue),rownames(only.H1047R.no.low.amp.erbb2.gene.exp.brca))
# 
# #C420R
# identical(rownames(only.C420R.no.low.amp.erbb2.gene.exp.brca),rownames(breast.gtex.tissue))
# table(rownames(only.C420R.no.low.amp.erbb2.gene.exp.brca)%in%rownames(breast.gtex.tissue))
# common_genes <- intersect(rownames(only.C420R.no.low.amp.erbb2.gene.exp.brca),rownames(breast.gtex.tissue))
# length(common_genes)
# 
# only.C420R.no.low.amp.erbb2.gene.exp.brca <- only.C420R.no.low.amp.erbb2.gene.exp.brca[match(common_genes,rownames(only.C420R.no.low.amp.erbb2.gene.exp.brca)),]
# breast.gtex.tissue <- breast.gtex.tissue[match(common_genes,rownames(breast.gtex.tissue)),]
# identical(rownames(breast.gtex.tissue),rownames(only.C420R.no.low.amp.erbb2.gene.exp.brca))
# 
# #E453K
# identical(rownames(only.E453K.no.low.amp.erbb2.gene.exp.brca),rownames(breast.gtex.tissue))
# table(rownames(only.E453K.no.low.amp.erbb2.gene.exp.brca)%in%rownames(breast.gtex.tissue))
# common_genes <- intersect(rownames(only.E453K.no.low.amp.erbb2.gene.exp.brca),rownames(breast.gtex.tissue))
# length(common_genes)
# 
# only.E453K.no.low.amp.erbb2.gene.exp.brca <- only.E453K.no.low.amp.erbb2.gene.exp.brca[match(common_genes,rownames(only.E453K.no.low.amp.erbb2.gene.exp.brca)),]
# breast.gtex.tissue <- breast.gtex.tissue[match(common_genes,rownames(breast.gtex.tissue)),]
# identical(rownames(breast.gtex.tissue),rownames(only.E453K.no.low.amp.erbb2.gene.exp.brca))
# 
# #E545K
# identical(rownames(only.E545K.no.low.amp.erbb2.gene.exp.brca),rownames(breast.gtex.tissue))
# table(rownames(only.E545K.no.low.amp.erbb2.gene.exp.brca)%in%rownames(breast.gtex.tissue))
# common_genes <- intersect(rownames(only.E545K.no.low.amp.erbb2.gene.exp.brca),rownames(breast.gtex.tissue))
# length(common_genes)
# 
# only.E545K.no.low.amp.erbb2.gene.exp.brca <- only.E545K.no.low.amp.erbb2.gene.exp.brca[match(common_genes,rownames(only.E545K.no.low.amp.erbb2.gene.exp.brca)),]
# breast.gtex.tissue <- breast.gtex.tissue[match(common_genes,rownames(breast.gtex.tissue)),]
# identical(rownames(breast.gtex.tissue),rownames(only.E545K.no.low.amp.erbb2.gene.exp.brca))
# 
# #E542K
# identical(rownames(only.E542K.no.low.amp.erbb2.gene.exp.brca),rownames(breast.gtex.tissue))
# table(rownames(only.E542K.no.low.amp.erbb2.gene.exp.brca)%in%rownames(breast.gtex.tissue))
# common_genes <- intersect(rownames(only.E542K.no.low.amp.erbb2.gene.exp.brca),rownames(breast.gtex.tissue))
# length(common_genes)
# 
# only.E542K.no.low.amp.erbb2.gene.exp.brca <- only.E542K.no.low.amp.erbb2.gene.exp.brca[match(common_genes,rownames(only.E542K.no.low.amp.erbb2.gene.exp.brca)),]
# breast.gtex.tissue <- breast.gtex.tissue[match(common_genes,rownames(breast.gtex.tissue)),]
# identical(rownames(breast.gtex.tissue),rownames(only.E542K.no.low.amp.erbb2.gene.exp.brca))
# 
# #Q546K
# identical(rownames(only.Q546K.no.low.amp.erbb2.gene.exp.brca),rownames(breast.gtex.tissue))
# table(rownames(only.Q546K.no.low.amp.erbb2.gene.exp.brca)%in%rownames(breast.gtex.tissue))
# common_genes <- intersect(rownames(only.Q546K.no.low.amp.erbb2.gene.exp.brca),rownames(breast.gtex.tissue))
# length(common_genes)
# 
# only.Q546K.no.low.amp.erbb2.gene.exp.brca <- only.Q546K.no.low.amp.erbb2.gene.exp.brca[match(common_genes,rownames(only.Q546K.no.low.amp.erbb2.gene.exp.brca)),]
# breast.gtex.tissue <- breast.gtex.tissue[match(common_genes,rownames(breast.gtex.tissue)),]
# identical(rownames(breast.gtex.tissue),rownames(only.Q546K.no.low.amp.erbb2.gene.exp.brca))

#Normalization starts
#log2cpm+1 normalization

# tpm normalization
for(i in 1:ncol(all.high.conf.lof.erbb2.gene.exp.brca)){
  all.high.conf.lof.erbb2.gene.exp.brca[,i] <- 1E6*all.high.conf.lof.erbb2.gene.exp.brca[,i]/sum(all.high.conf.lof.erbb2.gene.exp.brca[,i])
}

# for(i in 1:ncol(only.G118D.no.low.amp.erbb2.gene.exp.brca)){
#   only.G118D.no.low.amp.erbb2.gene.exp.brca[,i] <- 1E6*only.G118D.no.low.amp.erbb2.gene.exp.brca[,i]/sum(only.G118D.no.low.amp.erbb2.gene.exp.brca[,i])
# }
# for(i in 1:ncol(only.N345K.no.low.amp.erbb2.gene.exp.brca)){
#   only.N345K.no.low.amp.erbb2.gene.exp.brca[,i] <- 1E6*only.N345K.no.low.amp.erbb2.gene.exp.brca[,i]/sum(only.N345K.no.low.amp.erbb2.gene.exp.brca[,i])
# }
# for(i in 1:ncol(only.H1047L.no.low.amp.erbb2.gene.exp.brca)){
#   only.H1047L.no.low.amp.erbb2.gene.exp.brca[,i] <- 1E6*only.H1047L.no.low.amp.erbb2.gene.exp.brca[,i]/sum(only.H1047L.no.low.amp.erbb2.gene.exp.brca[,i])
# }
# for(i in 1:ncol(only.H1047R.no.low.amp.erbb2.gene.exp.brca)){
#   only.H1047R.no.low.amp.erbb2.gene.exp.brca[,i] <- 1E6*only.H1047R.no.low.amp.erbb2.gene.exp.brca[,i]/sum(only.H1047R.no.low.amp.erbb2.gene.exp.brca[,i])
# }
# for(i in 1:ncol(only.C420R.no.low.amp.erbb2.gene.exp.brca)){
#   only.C420R.no.low.amp.erbb2.gene.exp.brca[,i] <- 1E6*only.C420R.no.low.amp.erbb2.gene.exp.brca[,i]/sum(only.C420R.no.low.amp.erbb2.gene.exp.brca[,i])
# }
# for(i in 1:ncol(only.E453K.no.low.amp.erbb2.gene.exp.brca)){
#   only.E453K.no.low.amp.erbb2.gene.exp.brca[,i] <- 1E6*only.E453K.no.low.amp.erbb2.gene.exp.brca[,i]/sum(only.E453K.no.low.amp.erbb2.gene.exp.brca[,i])
# }
# for(i in 1:ncol(only.E545K.no.low.amp.erbb2.gene.exp.brca)){
#   only.E545K.no.low.amp.erbb2.gene.exp.brca[,i] <- 1E6*only.E545K.no.low.amp.erbb2.gene.exp.brca[,i]/sum(only.E545K.no.low.amp.erbb2.gene.exp.brca[,i])
# }
# for(i in 1:ncol(only.E542K.no.low.amp.erbb2.gene.exp.brca)){
#   only.E542K.no.low.amp.erbb2.gene.exp.brca[,i] <- 1E6*only.E542K.no.low.amp.erbb2.gene.exp.brca[,i]/sum(only.E542K.no.low.amp.erbb2.gene.exp.brca[,i])
# }
# for(i in 1:ncol(only.Q546K.no.low.amp.erbb2.gene.exp.brca)){
#   only.Q546K.no.low.amp.erbb2.gene.exp.brca[,i] <- 1E6*only.Q546K.no.low.amp.erbb2.gene.exp.brca[,i]/sum(only.Q546K.no.low.amp.erbb2.gene.exp.brca[,i])
# }

# ###log2cpm+1
# only.H1047L.no.low.amp.erbb2.gene.exp.brca.log2cpm <- apply(only.H1047L.no.low.amp.erbb2.gene.exp.brca,2,function(x){
#   y <- 1E6*x/sum(x) + 1
#   z <- log(y,2)
#   return(z)
# })
# 
# only.H1047L.no.low.amp.erbb2.gene.exp.brca.log2cpm.ges <- t(apply(only.H1047L.no.low.amp.erbb2.gene.exp.brca.log2cpm,1,function(x){
#   y <- (x - mean(x))/sd(x)
#   return(y)
# }))
# 
#for(i in 1:ncol(tcga.normal)){
#  tcga.normal[,i] <- 1E6*tcga.normal[,i]/sum(tcga.normal[,i])
#}
for(i in 1:ncol(colon.gtex.tissue)){
  colon.gtex.tissue[,i] <- 1E6*colon.gtex.tissue[,i]/sum(colon.gtex.tissue[,i])
}

# breast.gtex.tissue.log2cpm <- apply(breast.gtex.tissue,2,function(x){
#   y <- 1E6*x/sum(x) + 1
#   z <- log(y,2)
#   return(z)
# })
# breast.gtex.tissue.log2cpm.ges <- t(apply(breast.gtex.tissue.log2cpm,1,function(x){
#   y <- (x - mean(x))/sd(x)
#   return(y)
# }))
#check for identical genes
#all gof
identical(rownames(all.high.conf.lof.erbb2.gene.exp.brca),rownames(colon.gtex.tissue))
table(rownames(all.high.conf.lof.erbb2.gene.exp.brca)%in%rownames(colon.gtex.tissue))
common_genes <- intersect(rownames(all.high.conf.lof.erbb2.gene.exp.brca),rownames(colon.gtex.tissue))
length(common_genes)

all.high.conf.lof.erbb2.gene.exp.brca <- all.high.conf.lof.erbb2.gene.exp.brca[match(common_genes,rownames(all.high.conf.lof.erbb2.gene.exp.brca)),]
colon.gtex.tissue <- colon.gtex.tissue[match(common_genes,rownames(colon.gtex.tissue)),]
identical(rownames(colon.gtex.tissue),rownames(all.high.conf.lof.erbb2.gene.exp.brca))

# #G118D
# identical(rownames(only.G118D.no.low.amp.erbb2.gene.exp.brca),rownames(breast.gtex.tissue))
# table(rownames(only.G118D.no.low.amp.erbb2.gene.exp.brca)%in%rownames(breast.gtex.tissue))
# common_genes <- intersect(rownames(only.G118D.no.low.amp.erbb2.gene.exp.brca),rownames(breast.gtex.tissue))
# length(common_genes)
# 
# only.G118D.no.low.amp.erbb2.gene.exp.brca <- only.G118D.no.low.amp.erbb2.gene.exp.brca[match(common_genes,rownames(only.G118D.no.low.amp.erbb2.gene.exp.brca)),]
# breast.gtex.tissue <- breast.gtex.tissue[match(common_genes,rownames(breast.gtex.tissue)),]
# identical(rownames(breast.gtex.tissue),rownames(only.G118D.no.low.amp.erbb2.gene.exp.brca))
# 
# #N345K
# identical(rownames(only.N345K.no.low.amp.erbb2.gene.exp.brca),rownames(breast.gtex.tissue))
# table(rownames(only.N345K.no.low.amp.erbb2.gene.exp.brca)%in%rownames(breast.gtex.tissue))
# common_genes <- intersect(rownames(only.N345K.no.low.amp.erbb2.gene.exp.brca),rownames(breast.gtex.tissue))
# length(common_genes)
# 
# only.N345K.no.low.amp.erbb2.gene.exp.brca <- only.N345K.no.low.amp.erbb2.gene.exp.brca[match(common_genes,rownames(only.N345K.no.low.amp.erbb2.gene.exp.brca)),]
# breast.gtex.tissue <- breast.gtex.tissue[match(common_genes,rownames(breast.gtex.tissue)),]
# identical(rownames(breast.gtex.tissue),rownames(only.N345K.no.low.amp.erbb2.gene.exp.brca))
# 
# #H1047L
# # only.H1047L.no.low.amp.erbb2.gene.exp.brca<-only.H1047L.no.low.amp.erbb2.gene.exp.brca.log2cpm.ges
# # breast.gtex.tissue<-breast.gtex.tissue.log2cpm.ges
# identical(rownames(only.H1047L.no.low.amp.erbb2.gene.exp.brca),rownames(breast.gtex.tissue))
# table(rownames(only.H1047L.no.low.amp.erbb2.gene.exp.brca)%in%rownames(breast.gtex.tissue))
# common_genes <- intersect(rownames(only.H1047L.no.low.amp.erbb2.gene.exp.brca),rownames(breast.gtex.tissue))
# length(common_genes)
# 
# only.H1047L.no.low.amp.erbb2.gene.exp.brca <- only.H1047L.no.low.amp.erbb2.gene.exp.brca[match(common_genes,rownames(only.H1047L.no.low.amp.erbb2.gene.exp.brca)),]
# breast.gtex.tissue <- breast.gtex.tissue[match(common_genes,rownames(breast.gtex.tissue)),]
# identical(rownames(breast.gtex.tissue),rownames(only.H1047L.no.low.amp.erbb2.gene.exp.brca))
# 
# #H1047R
# identical(rownames(only.H1047R.no.low.amp.erbb2.gene.exp.brca),rownames(breast.gtex.tissue))
# table(rownames(only.H1047R.no.low.amp.erbb2.gene.exp.brca)%in%rownames(breast.gtex.tissue))
# common_genes <- intersect(rownames(only.H1047R.no.low.amp.erbb2.gene.exp.brca),rownames(breast.gtex.tissue))
# length(common_genes)
# 
# only.H1047R.no.low.amp.erbb2.gene.exp.brca <- only.H1047R.no.low.amp.erbb2.gene.exp.brca[match(common_genes,rownames(only.H1047R.no.low.amp.erbb2.gene.exp.brca)),]
# breast.gtex.tissue <- breast.gtex.tissue[match(common_genes,rownames(breast.gtex.tissue)),]
# identical(rownames(breast.gtex.tissue),rownames(only.H1047R.no.low.amp.erbb2.gene.exp.brca))
# 
# #C420R
# identical(rownames(only.C420R.no.low.amp.erbb2.gene.exp.brca),rownames(breast.gtex.tissue))
# table(rownames(only.C420R.no.low.amp.erbb2.gene.exp.brca)%in%rownames(breast.gtex.tissue))
# common_genes <- intersect(rownames(only.C420R.no.low.amp.erbb2.gene.exp.brca),rownames(breast.gtex.tissue))
# length(common_genes)
# 
# only.C420R.no.low.amp.erbb2.gene.exp.brca <- only.C420R.no.low.amp.erbb2.gene.exp.brca[match(common_genes,rownames(only.C420R.no.low.amp.erbb2.gene.exp.brca)),]
# breast.gtex.tissue <- breast.gtex.tissue[match(common_genes,rownames(breast.gtex.tissue)),]
# identical(rownames(breast.gtex.tissue),rownames(only.C420R.no.low.amp.erbb2.gene.exp.brca))
# 
# #E453K
# identical(rownames(only.E453K.no.low.amp.erbb2.gene.exp.brca),rownames(breast.gtex.tissue))
# table(rownames(only.E453K.no.low.amp.erbb2.gene.exp.brca)%in%rownames(breast.gtex.tissue))
# common_genes <- intersect(rownames(only.E453K.no.low.amp.erbb2.gene.exp.brca),rownames(breast.gtex.tissue))
# length(common_genes)
# 
# only.E453K.no.low.amp.erbb2.gene.exp.brca <- only.E453K.no.low.amp.erbb2.gene.exp.brca[match(common_genes,rownames(only.E453K.no.low.amp.erbb2.gene.exp.brca)),]
# breast.gtex.tissue <- breast.gtex.tissue[match(common_genes,rownames(breast.gtex.tissue)),]
# identical(rownames(breast.gtex.tissue),rownames(only.E453K.no.low.amp.erbb2.gene.exp.brca))
# 
# #E545K
# identical(rownames(only.E545K.no.low.amp.erbb2.gene.exp.brca),rownames(breast.gtex.tissue))
# table(rownames(only.E545K.no.low.amp.erbb2.gene.exp.brca)%in%rownames(breast.gtex.tissue))
# common_genes <- intersect(rownames(only.E545K.no.low.amp.erbb2.gene.exp.brca),rownames(breast.gtex.tissue))
# length(common_genes)
# 
# only.E545K.no.low.amp.erbb2.gene.exp.brca <- only.E545K.no.low.amp.erbb2.gene.exp.brca[match(common_genes,rownames(only.E545K.no.low.amp.erbb2.gene.exp.brca)),]
# breast.gtex.tissue <- breast.gtex.tissue[match(common_genes,rownames(breast.gtex.tissue)),]
# identical(rownames(breast.gtex.tissue),rownames(only.E545K.no.low.amp.erbb2.gene.exp.brca))
# 
# #E542K
# identical(rownames(only.E542K.no.low.amp.erbb2.gene.exp.brca),rownames(breast.gtex.tissue))
# table(rownames(only.E542K.no.low.amp.erbb2.gene.exp.brca)%in%rownames(breast.gtex.tissue))
# common_genes <- intersect(rownames(only.E542K.no.low.amp.erbb2.gene.exp.brca),rownames(breast.gtex.tissue))
# length(common_genes)
# 
# only.E542K.no.low.amp.erbb2.gene.exp.brca <- only.E542K.no.low.amp.erbb2.gene.exp.brca[match(common_genes,rownames(only.E542K.no.low.amp.erbb2.gene.exp.brca)),]
# breast.gtex.tissue <- breast.gtex.tissue[match(common_genes,rownames(breast.gtex.tissue)),]
# identical(rownames(breast.gtex.tissue),rownames(only.E542K.no.low.amp.erbb2.gene.exp.brca))
# 
# #Q546K
# identical(rownames(only.Q546K.no.low.amp.erbb2.gene.exp.brca),rownames(breast.gtex.tissue))
# table(rownames(only.Q546K.no.low.amp.erbb2.gene.exp.brca)%in%rownames(breast.gtex.tissue))
# common_genes <- intersect(rownames(only.Q546K.no.low.amp.erbb2.gene.exp.brca),rownames(breast.gtex.tissue))
# length(common_genes)
# 
# only.Q546K.no.low.amp.erbb2.gene.exp.brca <- only.Q546K.no.low.amp.erbb2.gene.exp.brca[match(common_genes,rownames(only.Q546K.no.low.amp.erbb2.gene.exp.brca)),]
# breast.gtex.tissue <- breast.gtex.tissue[match(common_genes,rownames(breast.gtex.tissue)),]
# identical(rownames(breast.gtex.tissue),rownames(only.Q546K.no.low.amp.erbb2.gene.exp.brca))

# bootstrap msVIPER analysis
library(ggplot2)

#test.mat <- as.matrix(only.G118D.no.low.amp.erbb2.gene.exp.brca)
#test.mat <- as.matrix(only.N345K.no.low.amp.erbb2.gene.exp.brca)
#test.mat <- as.matrix(only.H1047L.no.low.amp.erbb2.gene.exp.brca)
#test.mat <- as.matrix(only.H1047R.no.low.amp.erbb2.gene.exp.brca)
#test.mat <- as.matrix(only.C420R.no.low.amp.erbb2.gene.exp.brca)
#test.mat <- as.matrix(only.E453K.no.low.amp.erbb2.gene.exp.brca)
#test.mat <- as.matrix(only.E545K.no.low.amp.erbb2.gene.exp.brca)
#test.mat <- as.matrix(only.E542K.no.low.amp.erbb2.gene.exp.brca)
#test.mat <- as.matrix(only.Q546K.no.low.amp.erbb2.gene.exp.brca)
test.mat <- as.matrix(all.high.conf.lof.erbb2.gene.exp.brca)
test.mat <- test.mat[is.finite(rowSums(test.mat)),]
any(is.na(test.mat))

ref.mat <- as.matrix(colon.gtex.tissue)
ref.mat <- ref.mat[is.finite(rowSums(ref.mat)),]
any(is.na(ref.mat))

#ref.mat <- as.matrix(alk.wt.samples)

# name conversion
entrez.names <- rownames(test.mat)
hugo.names <- as.character(n1platformReports::entrez2gene(entrez.names))
hugo.names[is.na(hugo.names)] <- paste("None",seq(from = 1, to = sum(is.na(hugo.names)), by = 1), sep = "_")
names(hugo.names) <- entrez.names

#load regulon and prune
load('../../brca-gtex/pik3ca_brca/rda/brca_regulon.rda')
#pregul <- pruneRegulon(regul, cutoff = 50)

pregul.brca <- pruneRegulon(regul, cutoff = 50)
current.interactome <- pregul.brca

new.interactome<-current.interactome[which(names(current.interactome)%in%rownames(tf.cotf))]
current.interactome<-new.interactome

# compute signature
current.sig.erbb2.brca <- bootstrapTtest(x = test.mat, y = ref.mat, per = 100)

# compute nullmodel
current.nullmodel.erbb2.brca <- ttestNull(x = test.mat, y = ref.mat, per = 1000)

# compute differential protein activity
current.mrs.erbb2.brca <- msviper(ges = current.sig.erbb2.brca, regulon = current.interactome, nullmodel = current.nullmodel.erbb2.brca)
current.integration.method <- "mean"

current.mrs.erbb2.brca <- bootstrapmsviper(current.mrs.erbb2.brca, method = current.integration.method)
#saveRDS(current.mrs.erbb2.brca,"high.conf.mut.erbb2.G118D.brca.bootstrap.msviper.rds")
#saveRDS(current.mrs.erbb2.brca,"high.conf.mut.erbb2.N345K.brca.bootstrap.msviper.rds")
#saveRDS(current.mrs.erbb2.brca,"high.conf.mut.erbb2.H1047L.brca.bootstrap.msviper.rds")
#saveRDS(current.mrs.erbb2.brca,"high.conf.mut.erbb2.H1047R.brca.bootstrap.msviper.rds")
#saveRDS(current.mrs.erbb2.brca,"high.conf.mut.erbb2.C420R.brca.bootstrap.msviper.rds")
#saveRDS(current.mrs.erbb2.brca,"high.conf.mut.erbb2.E453K.brca.bootstrap.msviper.rds")
#saveRDS(current.mrs.erbb2.brca,"high.conf.mut.erbb2.E545K.brca.bootstrap.msviper.rds")
#saveRDS(current.mrs.erbb2.brca,"high.conf.mut.erbb2.E542K.brca.bootstrap.msviper.rds")
#saveRDS(current.mrs.erbb2.brca,"high.conf.mut.erbb2.Q546K.brca.bootstrap.msviper.rds")
saveRDS(current.mrs.erbb2.brca,"all.high.conf.mut.lof.erbb2.brca.bootstrap.msviper.rds")
#saveRDS(current.mrs.erbb2.brca,"all.high.conf.mut.gof.erbb2.brca.bootstrap.msviper.cluster.1.rds")
#saveRDS(current.mrs.erbb2.brca,"all.high.conf.mut.gof.erbb2.brca.bootstrap.msviper.cluster.2.rds")
#saveRDS(current.mrs.erbb2.brca,"all.high.conf.mut.gof.erbb2.brca.bootstrap.msviper.cluster.3.rds")
#saveRDS(current.mrs.erbb2.brca,"all.high.conf.mut.lof.yp5.brca.bootstrap.msviper.cluster.4.rds")
current.msa.erbb2.brca <- msviperAnnot(current.mrs.erbb2.brca, hugo.names)
summary(current.msa.erbb2.brca,20)
# current.mrs.erbb2.brca.ledge<-ledge(current.msa.erbb2.brca)
# summary(current.mrs.erbb2.brca.ledge,50)
# write.csv(summary(current.mrs.erbb2.brca.ledge,2313),file="erbb2.brca.ledge.csv")
#current.mrs.erbb2.ucec<-readRDS(file="all.high.conf.mut.gof.erbb2.ucec.bootstrap.msviper.rds")
#dim(current.mrs$signature)
c.s.all.high.conf.samples.erbb2.brca<-
  cbind(
    current.mrs.erbb2.brca$es$nes)

colnames(c.s.all.high.conf.samples.erbb2.brca)<- c(
  
  "High.Conf.GOF"
)
saveRDS(c.s.all.high.conf.samples.erbb2.brca,file="consensus_signature_erbb2_brca.rds")
# c.s.all.high.conf.samples.brca.cluster.2<-
#   cbind(
#     current.mrs.erbb2.brca$es$nes)
# 
# colnames(c.s.all.high.conf.samples.brca.cluster.2)<- c(
#   
#   "High.Conf.GOF"
# )
# saveRDS(c.s.all.high.conf.samples.brca.cluster.2,file="consensus_signature_cluster_2.rds")
# c.s.all.high.conf.samples.brca.cluster.3<-
#   cbind(
#     current.mrs.erbb2.brca$es$nes)
# 
# colnames(c.s.all.high.conf.samples.brca.cluster.3)<- c(
#   
#   "High.Conf.GOF"
# )
# saveRDS(c.s.all.high.conf.samples.brca.cluster.3,file="consensus_signature_cluster_3.rds")
# c.s.all.high.conf.samples.brca.cluster.4<-
#   cbind(
#     current.mrs.erbb2.brca$es$nes)
# 
# colnames(c.s.all.high.conf.samples.brca.cluster.4)<- c(
#   
#   "High.Conf.GOF"
# )
# saveRDS(c.s.all.high.conf.samples.brca.cluster.4,file="consensus_signature_cluster_4.rds")

# #restrict to tf/cotf
# common.tf.cotf <- intersect(rownames(current.mrs$signature),rownames(tf.cotf))
# length(common.tf.cotf)
# 
# current.mrs.tf.cotf <- current.mrs$signature[match(common.tf.cotf,rownames(current.mrs$signature)),]
# dim(current.mrs.tf.cotf)
# saveRDS(current.mrs.tf.cotf,"high.conf.mut.egfr.bootstrap.msviper.tf.cotf.rds")

# analysis of unknown mutations


test.1<-erbb2_brca_all_labels
for(i in 1:length(test.1)) { 
  nam <- noquote(paste0(test.1[i], '.erbb2.gene.exp.brca'))
  nam.1 <- noquote(paste0(test.1[i], '.erbb2.brca.samples'))
  assign(nam, (brca.tcga.tumor[,intersect(colnames(brca.tcga.tumor),get(nam.1))]))
  # assign(nam, (ucec.tcga.tumor[,intersect(colnames(ucec.tcga.tumor),get(nam))]))
  
}



known.unknown.erbb2.brca.mut.samples.1<-cbind(
  S310F.erbb2.gene.exp.brca,
  A1160V.erbb2.gene.exp.brca,
  A355Qfs.76.erbb2.gene.exp.brca,
  ADGRA3.ERBB2.erbb2.gene.exp.brca,
  AP2B1.ERBB2.erbb2.gene.exp.brca,
  CACNB1.ERBB2.erbb2.gene.exp.brca,
  CORO1B.ERBB2.erbb2.gene.exp.brca,
  D769H.erbb2.gene.exp.brca,
  E405D.erbb2.gene.exp.brca,
  E939G.erbb2.gene.exp.brca,
  EME1.ERBB2.erbb2.gene.exp.brca,
  ERBB2.ABI3.erbb2.gene.exp.brca,
  ERBB2.ARL5C.erbb2.gene.exp.brca,
  ERBB2.MED1.erbb2.gene.exp.brca,
  ERBB2.PPP1R1B.erbb2.gene.exp.brca,
  ERBB2.PSMB3.erbb2.gene.exp.brca,
  ERBB2.SDF4.erbb2.gene.exp.brca,
  ERBB2.SLC29A3.erbb2.gene.exp.brca,
  ERBB2.SRCIN1.erbb2.gene.exp.brca,
  ERBB2.ZAN.erbb2.gene.exp.brca,
  G309A.erbb2.gene.exp.brca,
  G621Afs.31.erbb2.gene.exp.brca,
  H470Q.erbb2.gene.exp.brca,
  L755M.erbb2.gene.exp.brca,
  L755W.erbb2.gene.exp.brca,
  MYO18A.ERBB2.erbb2.gene.exp.brca,
  PIP4K2B.ERBB2.erbb2.gene.exp.brca,
  R340G.erbb2.gene.exp.brca,
  R678Q.erbb2.gene.exp.brca,
  S305C.erbb2.gene.exp.brca,
  T306M.erbb2.gene.exp.brca,
  V797A.erbb2.gene.exp.brca,
  V842I.erbb2.gene.exp.brca,
  W482Gfs.74.erbb2.gene.exp.brca,
  Y1127_A1129del.erbb2.gene.exp.brca,
  L755S.erbb2.gene.exp.brca,
  V777L.erbb2.gene.exp.brca,
  D769Y.erbb2.gene.exp.brca
)


colnames(known.unknown.erbb2.brca.mut.samples.1)<-
  cbind(
     'S310F_TCGA.D8.A27G','S310F_TCGA.C8.A274',
         'A1160V_TCGA.BH.A0DZ',
         'A355Qfs*76_TCGA.EW.A2FV',
         'ADGRA3.ERBB2_TCGA.D8.A1X5',
         'AP2B1.ERBB2_TCGA.A8.A08X',
         'CACNB1.ERBB2_TCGA.D8.A1XJ',
         'CORO1B.ERBB2_TCGA.C8.A12Q',
         'D769H_TCGA.C8.A135',
         'E405D_TCGA.EW.A2FR',
         'E939G_TCGA.LL.A740',
         'EME1.ERBB2_TCGA.BH.A1F2',
         'ERBB2.ABI3_TCGA.C8.A132',
         'ERBB2.ARL5C_TCGA.AO.A0JM',
         'ERBB2.MED1_TCGA.D8.A1XJ',
         'ERBB2.PPP1R1B_TCGA.C8.A12Z',
         'ERBB2.PSMB3_TCGA.A2.A1G1',
         'ERBB2.SDF4_TCGA.C8.A8HP',
         'ERBB2.SLC29A3_TCGA.AR.A254',
         'ERBB2.SRCIN1_TCGA.BH.A1EV',
         'ERBB2.ZAN_TCGA.OL.A5RY',
         'G309A_TCGA.A8.A06Z',
         'G621Afs*31_TCGA.D8.A27V',
         'H470Q_TCGA.EW.A1PD',
         'L755M_TCGA.A2.A0T6',
         'L755W_TCGA.A2.A0T6',
         'MYO18A.ERBB2_TCGA.C8.A275',
         'PIP4K2B.ERBB2_TCGA.A8.A08X',
         'R340G_TCGA.BH.A1FU',
         'R678Q_TCGA.A2.A0T6',
         'S305C_TCGA.C8.A3M7',
         'T306M_TCGA.AN.A046',
         'V797A_TCGA.AO.A128',
         'V842I_TCGA.A8.A08Z',
         'W482Gfs*74_TCGA.EW.A2FV',
         'Y1127_A1129del_TCGA.D8.A1JG',
         #known
         'L755S_TCGA.A8.A0A6','L755S_TCGA.AC.A3YI','L755S_TCGA.A8.A0AB','L755S_TCGA.D8.A1XM','L755S_TCGA.BH.A18P',
     'V777L_TCGA.OL.A5D6','V777L_TCGA.4H.AAAK','V777L_TCGA.BH.A0C1',
     'D769Y_TCGA.E9.A1R5','D769Y_TCGA.BH.A1FE'
  )

curr.stouf<-c('S310F',
              'L755S',
              'V777L',
              'D769Y'
)
for(i in 1:length(curr.stouf)){
  nam <- noquote(paste0('current.stouffer.',curr.stouf[i], '.erbb2.gene.exp.brca'))
  nam1 <- noquote(paste0(curr.stouf[i], '.erbb2.gene.exp.brca'))
  assign(nam, apply(get(nam1),1,function(x){sum(x)/sqrt(length(x))}))
  
}

known.unknown.erbb2.brca.mut.samples.1<-cbind(
  current.stouffer.S310F.erbb2.gene.exp.brca,
  A1160V.erbb2.gene.exp.brca,
  A355Qfs.76.erbb2.gene.exp.brca,
  ADGRA3.ERBB2.erbb2.gene.exp.brca,
  AP2B1.ERBB2.erbb2.gene.exp.brca,
  CACNB1.ERBB2.erbb2.gene.exp.brca,
  CORO1B.ERBB2.erbb2.gene.exp.brca,
  D769H.erbb2.gene.exp.brca,
  E405D.erbb2.gene.exp.brca,
  E939G.erbb2.gene.exp.brca,
  EME1.ERBB2.erbb2.gene.exp.brca,
  ERBB2.ABI3.erbb2.gene.exp.brca,
  ERBB2.ARL5C.erbb2.gene.exp.brca,
  ERBB2.MED1.erbb2.gene.exp.brca,
  ERBB2.PPP1R1B.erbb2.gene.exp.brca,
  ERBB2.PSMB3.erbb2.gene.exp.brca,
  ERBB2.SDF4.erbb2.gene.exp.brca,
  ERBB2.SLC29A3.erbb2.gene.exp.brca,
  ERBB2.SRCIN1.erbb2.gene.exp.brca,
  ERBB2.ZAN.erbb2.gene.exp.brca,
  G309A.erbb2.gene.exp.brca,
  G621Afs.31.erbb2.gene.exp.brca,
  H470Q.erbb2.gene.exp.brca,
  L755M.erbb2.gene.exp.brca,
  L755W.erbb2.gene.exp.brca,
  MYO18A.ERBB2.erbb2.gene.exp.brca,
  PIP4K2B.ERBB2.erbb2.gene.exp.brca,
  R340G.erbb2.gene.exp.brca,
  R678Q.erbb2.gene.exp.brca,
  S305C.erbb2.gene.exp.brca,
  T306M.erbb2.gene.exp.brca,
  V797A.erbb2.gene.exp.brca,
  V842I.erbb2.gene.exp.brca,
  W482Gfs.74.erbb2.gene.exp.brca,
  Y1127_A1129del.erbb2.gene.exp.brca,
  current.stouffer.L755S.erbb2.gene.exp.brca,
  current.stouffer.V777L.erbb2.gene.exp.brca,
  current.stouffer.D769Y.erbb2.gene.exp.brca
)


colnames(known.unknown.erbb2.brca.mut.samples.1)<-
  cbind(
    'S310F',
    'A1160V',
    'A355Qfs*76',
    'ADGRA3.ERBB2',
    'AP2B1.ERBB2',
    'CACNB1.ERBB2',
    'CORO1B.ERBB2',
    'D769H',
    'E405D',
    'E939G',
    'EME1.ERBB2',
    'ERBB2.ABI3',
    'ERBB2.ARL5C',
    'ERBB2.MED1',
    'ERBB2.PPP1R1B',
    'ERBB2.PSMB3',
    'ERBB2.SDF4',
    'ERBB2.SLC29A3',
    'ERBB2.SRCIN1',
    'ERBB2.ZAN',
    'G309A',
    'G621Afs*31',
    'H470Q',
    'L755M',
    'L755W',
    'MYO18A.ERBB2',
    'PIP4K2B.ERBB2',
    'R340G',
    'R678Q',
    'S305C',
    'T306M',
    'V797A',
    'V842I',
    'W482Gfs*74',
    'Y1127_A1129del',
    #known
    'L755S',
    'V777L',
    'D769Y'
  )

unknown.erbb2.brca.mut.samples<-known.unknown.erbb2.brca.mut.samples.1
colnames(unknown.erbb2.brca.mut.samples)
dim(unknown.erbb2.brca.mut.samples)

#check for identical genes
colon.gtex.tissue <- read.table("~/Documents/viper-docs/benchmarking_1/signature_with_gtex/batch_corrected_data/data/batch-correct/gtex/breast-gtex.txt", header=TRUE)
dim(colon.gtex.tissue)


identical(rownames(unknown.erbb2.brca.mut.samples),rownames(colon.gtex.tissue))
table(rownames(unknown.erbb2.brca.mut.samples)%in%rownames(colon.gtex.tissue))
common_genes <- intersect(rownames(unknown.erbb2.brca.mut.samples),rownames(colon.gtex.tissue))
length(common_genes)

unknown.erbb2.brca.mut.samples <- unknown.erbb2.brca.mut.samples[match(common_genes,rownames(unknown.erbb2.brca.mut.samples)),]
colon.gtex.tissue <- colon.gtex.tissue[match(common_genes,rownames(colon.gtex.tissue)),]
identical(rownames(colon.gtex.tissue),rownames(unknown.erbb2.brca.mut.samples))

#Normalization starts
#log2cpm+1 normalization

# tpm normalization
for(i in 1:ncol(unknown.erbb2.brca.mut.samples)){
  unknown.erbb2.brca.mut.samples[,i] <- 1E6*unknown.erbb2.brca.mut.samples[,i]/sum(unknown.erbb2.brca.mut.samples[,i])
}
#for(i in 1:ncol(tcga.normal)){
#  tcga.normal[,i] <- 1E6*tcga.normal[,i]/sum(tcga.normal[,i])
#}
for(i in 1:ncol(colon.gtex.tissue)){
  colon.gtex.tissue[,i] <- 1E6*colon.gtex.tissue[,i]/sum(colon.gtex.tissue[,i])
}

#check for identical genes
identical(rownames(unknown.erbb2.brca.mut.samples),rownames(colon.gtex.tissue))
table(rownames(unknown.erbb2.brca.mut.samples)%in%rownames(colon.gtex.tissue))
common_genes <- intersect(rownames(unknown.erbb2.brca.mut.samples),rownames(colon.gtex.tissue))
length(common_genes)

unknown.erbb2.brca.mut.samples <- unknown.erbb2.brca.mut.samples[match(common_genes,rownames(unknown.erbb2.brca.mut.samples)),]
colon.gtex.tissue <- colon.gtex.tissue[match(common_genes,rownames(colon.gtex.tissue)),]
identical(rownames(colon.gtex.tissue),rownames(unknown.erbb2.brca.mut.samples))

#vipermat.unknown.mut.erbb2.brca.tf.cotf
unknown.erbb2.brca.mut.samples <- unknown.erbb2.brca.mut.samples[is.finite(rowSums(unknown.erbb2.brca.mut.samples)),]
any(is.na(unknown.erbb2.brca.mut.samples))
colon.gtex.tissue <- colon.gtex.tissue[is.finite(rowSums(colon.gtex.tissue)),]
any(is.na(colon.gtex.tissue))

vipsig.unknown.mut.erbb2.brca <- viperSignature(eset = as.matrix(unknown.erbb2.brca.mut.samples), 
                                                ref = as.matrix(colon.gtex.tissue), per = 1000, method = "ttest", 
                                                verbose = TRUE)

vipermat.unknown.mut.erbb2.brca<- viper(vipsig.unknown.mut.erbb2.brca, pregul.brca,minsize = 25)
colnames(vipermat.unknown.mut.erbb2.brca)
#restrict to tf/cotf
common.tf.cotf <- intersect(rownames(vipermat.unknown.mut.erbb2.brca),rownames(tf.cotf))
length(common.tf.cotf)

vipermat.unknown.mut.erbb2.brca.tf.cotf <- vipermat.unknown.mut.erbb2.brca[match(common.tf.cotf,rownames(vipermat.unknown.mut.erbb2.brca)),]
dim(vipermat.unknown.mut.erbb2.brca.tf.cotf)
#saveRDS(vipermat.unknown.mut.erbb2.brca.tf.cotf,file="vipermat.unknown.mut.erbb2.with.G118D.brca.tf.cotf.rds")
#saveRDS(vipermat.unknown.mut.erbb2.brca.tf.cotf,file="vipermat.unknown.mut.erbb2.with.N345K.brca.tf.cotf.rds")
#saveRDS(vipermat.unknown.mut.erbb2.brca.tf.cotf,file="vipermat.unknown.mut.erbb2.with.H1047L.brca.tf.cotf.rds")
#saveRDS(vipermat.unknown.mut.erbb2.brca.tf.cotf,file="vipermat.unknown.mut.erbb2.with.H1047R.brca.tf.cotf.rds")
#saveRDS(vipermat.unknown.mut.erbb2.brca.tf.cotf,file="vipermat.unknown.mut.erbb2.with.C420R.brca.tf.cotf.rds")
#saveRDS(vipermat.unknown.mut.erbb2.brca.tf.cotf,file="vipermat.unknown.mut.erbb2.with.E453K.brca.tf.cotf.rds")
#saveRDS(vipermat.unknown.mut.erbb2.brca.tf.cotf,file="vipermat.unknown.mut.erbb2.with.E545K.brca.tf.cotf.rds")
#saveRDS(vipermat.unknown.mut.erbb2.brca.tf.cotf,file="vipermat.unknown.mut.erbb2.with.E542K.brca.tf.cotf.rds")
#saveRDS(vipermat.unknown.mut.erbb2.brca.tf.cotf,file="vipermat.unknown.mut.erbb2.with.Q546K.brca.tf.cotf.rds")
#saveRDS(vipermat.unknown.mut.erbb2.brca.tf.cotf,file="vipermat.unknown.mut.erbb2.with.all.high.conf.gof.brca.tf.cotf.rds")
#saveRDS(vipermat.unknown.mut.erbb2.brca.tf.cotf,file="vipermat.unknown.mut.erbb2.with.all.high.conf.gof.brca.tf.cotf.smaples.separated.rds")
#saveRDS(vipermat.unknown.mut.erbb2.brca.tf.cotf,file="vipermat.unknown.mut.erbb2.with.all.high.conf.gof.brca.tf.cotf.samaples.separated.cluster.1.rds")
#saveRDS(vipermat.unknown.mut.erbb2.brca.tf.cotf,file="vipermat.unknown.mut.erbb2.with.all.high.conf.gof.brca.tf.cotf.samaples.separated.cluster.2.rds")
#saveRDS(vipermat.unknown.mut.erbb2.brca.tf.cotf,file="vipermat.unknown.mut.erbb2.with.all.high.conf.gof.brca.tf.cotf.samaples.separated.cluster.3.rds")
#saveRDS(vipermat.unknown.mut.erbb2.brca.tf.cotf,file="vipermat.unknown.mut.erbb2.with.all.high.conf.gof.brca.tf.cotf.samaples.separated.cluster.4.rds")
#saveRDS(vipermat.unknown.mut.erbb2.brca.tf.cotf,file="vipermat.unknown.mut.erbb2.with.all.high.conf.lof.brca.tf.cotf.rds")
#saveRDS(vipermat.unknown.mut.erbb2.brca.tf.cotf,file="vipermat.known.unknown.mut.erbb2.with.all.high.conf.lof.brca.tf.cotf.rds")
saveRDS(vipermat.unknown.mut.erbb2.brca.tf.cotf,file="vipermat.known.unknown.mut.erbb2.with.all.high.conf.gof.brca.tf.cotf.merge.mut.rds")

#saveRDS(c.s.all.high.conf.samples.brca.cluster.1,file="consensus_signature_cluster_1.rds")

# any(is.na(vipermat.unknown.mut.erbb2.brca.tf.cotf))
# 
# vipermat.unknown.mut.erbb2.brca.tf.cotf<-readRDS(file="vipermat.unknown.mut.erbb2.with.all.high.conf.gof.brca.tf.cotf.smaples.separated.rds")

#all.neomorphs
plot.data<-readRDS(file="erbb2.brca.gof.lof.neo.plot.data.rds")
'S310F_TCGA.C8.A274'
'S310F_TCGA.D8.A27G'
'S305C_TCGA.C8.A3M7'
'D769Y_TCGA.E9.A1R5'

neomorphs.erbb2.brca<-cbind(
  plot.data[,colnames(plot.data)=='S310F_TCGA.C8.A274'],
  plot.data[,colnames(plot.data)=='S310F_TCGA.D8.A27G'],
  plot.data[,colnames(plot.data)=='S305C_TCGA.C8.A3M7'],
  plot.data[,colnames(plot.data)=='D769Y_TCGA.E9.A1R5']
)
colnames(neomorphs.erbb2.brca)<-cbind(
  'S310F_TCGA.C8.A274',
  'S310F_TCGA.D8.A27G',
  'S305C_TCGA.C8.A3M7',
  'D769Y_TCGA.E9.A1R5'
) 
saveRDS(neomorphs.erbb2.brca,file="neomorphs.erbb2.brca.rds")

#compare neomorphs

neomorphs.erbb2.brca<-readRDS(file="neomorphs.erbb2.brca.rds")
dim(neomorphs.erbb2.brca)

library(fpc)
library(cluster)

#tumor.clust <- pamk(tumor.dist)$pamobject

#plot(silhouette(tumor.clust))

#any(is.na(vipermat))
tumor.vpmat<-neomorphs.erbb2.brca
dim(tumor.vpmat)
any(is.na(tumor.vpmat))

# #viper similarity
# vipsim <- viperSimilarity(tumor.vpmat)
# dist.matr <- as.dist(vipsim)
# 
# #pamk clustering
# library(fpc)
# pam.clust<-pamk(dist.matr,krange=1:3,criterion="asw", usepam=TRUE,
#                 scaling=FALSE, alpha=0.001, diss=inherits(data, "dist"),
#                 critout=FALSE, ns=10, seed=NULL)
# pam.clust$nc
# 
# current.silinfo <- as.data.frame(pam.clust$pamobject$silinfo$widths)
# 
# library(cluster)
# current.bootsil <- lapply(1:1000,function(x,ref.clust,ref.dist){
#   y <- sample(x = ref.clust, size = length(ref.clust), replace = FALSE)
#   z <- silhouette(x = y, dist = dist.matr)[,3]
#   #z <- silhouette(x = y, dist = dist.matr)[,1]
#   return(z)
# }, ref.clust = as.numeric(pam.clust$pamobject$clustering), ref.dist = dist.matr)
# current.bootsil <- unlist(current.bootsil)
# 
# # compute true silhouette score right-tailed p-values
# current.silinfo$pvalue <- pnorm(q = current.silinfo$sil_width, 
#                                 mean = mean(current.bootsil), sd = sd(current.bootsil), 
#                                 lower.tail = FALSE) - pnorm(q = 1, 
#                                                             mean = mean(current.bootsil), 
#                                                             sd = sd(current.bootsil), 
#                                                             lower.tail = FALSE)
# 
# # adjust p-values for FDR 
# current.silinfo$padj <- p.adjust(current.silinfo$pvalue, method = "BH")
# 
# # assign samples with FDR-adjusted p-values > 0.05 to an outlier cluster 0
# # current.silinfo[which(current.silinfo$padj > 0.05),"cluster"] <- 0
# 
# # order the silhouette score info by cluster and silhoeutte score
# current.silinfo <- current.silinfo[order(current.silinfo$cluster, current.silinfo$sil_width, decreasing = TRUE),]
# 
# # plot distance matrix after rearranging the results by silhouette score and cluster and annotation responders and non-responders
# 
# annot.col <- current.silinfo[,c("cluster","sil_width")]
# colnames(annot.col) <- c("Cluster","SilScore")
# annot.col$Cluster <- factor(annot.col$Cluster, 
#                             levels = sort(unique(annot.col$Cluster), decreasing = TRUE))
# annot.col <- annot.col[,c("SilScore","Cluster")]
# 
# #annot.col.brca<-annot.col
# #saveRDS(annot.col.brca,file="annot.col.brca.rds")
# #vipermat <- vipermat[is.finite(rowSums(vipermat)),]
# #plot.data<-vipermat
# #rownames(annot.col)<-gsub('^.{12}','',rownames(annot.col))
# rownames(annot.col)
# plot.data <- tumor.vpmat[,match(rownames(annot.col),colnames(tumor.vpmat))] 
# 
# any(is.na(plot.data))

#plot.data <- plot.data[,is.finite(colSums(plot.data))]
dim(plot.data)
#plot.data[is.na(plot.data)] <- 0

mat_breaks <- quantile_breaks(tumor.vpmat,n = 100)

hm.3<-pheatmap(tumor.vpmat, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))
               (length(mat_breaks)), 
               fontsize_row = 5, fontsize_col = 8, 
               #main = "BRAF mutation cluster brca Test 1",
               main = "Predicted erbb2 Neomorphs: brca",
               show_rownames = FALSE, show_colnames = TRUE, breaks = mat_breaks, 
               cluster_cols = TRUE, cluster_rows = TRUE)
save_pheatmap_pdf(hm.3,file="neomorph_comparison_erbb2_brca_all_tf_cotf.pdf")
#, annotation_col = annot.col)
#save_pheatmap_pdf(test1,"BRAF_mutation_cluster_brca_test1.pdf")
# save_pheatmap_pdf(test1,"BRCA_Unsupervised_Internal.pdf")
# saveRDS(plot.data,"BRCA.Unsup.Internal.rds")

# keep.index <- apply(plot.data,1,function(x,ref){
#   y <- pairwise.t.test(x = x,g = ref)
#   z <- min(y$p.value, na.rm = TRUE)
#   return(z)
# }, ref = annot.col$Cluster)
# 
# 
# keep.num <- length(rownames(plot.data))
# 
# final.keep.index <- sort(keep.index, decreasing = FALSE)[1:keep.num]
# 
# #rownames(plot.data)<-getSYMBOL(rownames(plot.data),data='org.Hs.eg')
# final.plot.data <- plot.data[match(names(final.keep.index),rownames(plot.data)),]
# dim(final.plot.data)
# rownames(final.plot.data)

# #mat_breaks <- quantile_breaks(final.plot.data,n = 100)
# test2<-pheatmap(final.plot.data, color = colorRampPalette(rev(brewer.pal(n = 7, 
#                                                                          name = "RdBu")))
#                 (length(mat_breaks)), fontsize_row = 5, fontsize_col = 5,
#                 #main = "BRAF mutation cluster brca Test 2",
#                 main = "brca Unsupervised (+GTEx) test2",
#                 show_rownames = FALSE, show_colnames = FALSE, breaks = mat_breaks, 
#                 cluster_cols = FALSE, cluster_rows = TRUE, annotation_col = annot.col)
# #save_pheatmap_pdf(test2,"BRAF_mutation_cluster_brca_test2.pdf")
# save_pheatmap_pdf(test2,"brca_Unsupervised_GTEx_test2.pdf")
# dim(final.plot.data)
# any(is.na(final.plot.data))

tf.cotf <- read.table(file="~/Desktop/tf-cotf-homo-current.txt")
row.names(tf.cotf) <- tf.cotf[,1]
tf.cotf<-tf.cotf[,-1]
rownames(tf.cotf)

#rownames(tf.cotf)<-getSYMBOL(rownames(tf.cotf),data='org.Hs.eg')
final.plot.data<-tumor.vpmat
common.tf.cotf <- intersect(rownames(final.plot.data),rownames(tf.cotf))
length(common.tf.cotf)

plot.data.tf.cotf <- final.plot.data[match(common.tf.cotf,rownames(final.plot.data)),]
dim(plot.data.tf.cotf)
#saveRDS(plot.data.tf.cotf,file="braf.brca.plot.data.rds")
#mat_breaks <- quantile_breaks(plot.data.tf.cotf,n = 100)

plot.data.tf.cotf <- plot.data.tf.cotf[is.finite(rowSums(plot.data.tf.cotf)),]
# pheatmap(plot.data.tf.cotf, color = colorRampPalette(rev(brewer.pal(n = 7, 
#                                                                            name = "RdBu")))
#                 (length(mat_breaks)), fontsize_row = 5, fontsize_col = 5,
#                 #main = "BRAF mutation cluster brca TF/Co-TF",
#                 main = "BRCA Unsupervised\n TF/Co-TF (+Internal)",
#                 show_rownames = FALSE, show_colnames = FALSE, breaks = mat_breaks, 
#                 cluster_cols = FALSE, cluster_rows = TRUE, annotation_col = annot.col)
# #save_pheatmap_pdf(test3,"BRAF_mutation_cluster_brca_tf_cotfs.pdf")
# save_pheatmap_pdf(test3,"BRCA_tf_cotfs_Internal.pdf")
# saveRDS(plot.data.tf.cotf,"BRCA.Unsup.all.tf.cotf.Internal.rds")
# dim(plot.data.tf.cotf)

# zscore_data_egfr <- t(apply(plot.data.tf.cotf,1,function(x){
#   y <- (x - mean(x))/sd(x)
#   return(y)
# }))
# 
# current.stouffer <- t(apply(zscore_data_egfr,1,function(x,clust,weight){
#   temp.data <- data.frame(nes = x, cluster = clust, weight = weight)
#   temp.data$w_z <- temp.data$nes*temp.data$weight
#   temp.data$w_w <- temp.data$weight*temp.data$weight
#   temp.agg <- aggregate(.~cluster, temp.data,sum)
#   temp.agg$res <- temp.agg$w_z/sqrt(temp.agg$w_w)
#   temp.res <- temp.agg$res
#   return(temp.res)
# }, clust = annot.col$Cluster, weight = annot.col$SilScore))
# colnames(current.stouffer) <- levels(annot.col$Cluster)
# 
# rank.stouffer <- apply(current.stouffer,2,rank)
# 
# keep.num <- 100
# keep.index <- apply(rank.stouffer,1,function(x){
#   y <- ((max(x) > (nrow(rank.stouffer) + 1 - keep.num)) | min(x) <= keep.num)
#   return(y)
# })
# 
# temp.plot.data.2 <- plot.data.tf.cotf[keep.index,]
# 
# mat_breaks <- quantile_breaks(temp.plot.data.2,n = 100)
# library(org.Hs.eg.db)
# library(annotate)
# rownames(temp.plot.data.2)<-getSYMBOL(rownames(temp.plot.data.2),data='org.Hs.eg')
# 
# pheatmap(temp.plot.data.2, color = colorRampPalette(rev(brewer.pal(n = 7, 
#                                                                    name = "RdBu")))
#          (length(mat_breaks)), fontsize_row = 4, fontsize_col = 5,
#          main = "Predicted erbb2 Neomorphs: BRCA",
#          show_rownames = TRUE, show_colnames = FALSE, breaks = mat_breaks, 
#          cluster_cols = FALSE, cluster_rows = TRUE, annotation_col = annot.col)
# 
# var.keep.index <- apply(plot.data.tf.cotf,1,sd)
# var.keep.index <- (rank(var.keep.index) > length(var.keep.index) - 100)
# 
# temp.plot.data.3 <- plot.data.tf.cotf[var.keep.index,]
# 
# mat_breaks <- quantile_breaks(temp.plot.data.3,n = 100)
# #rownames(temp.plot.data.3)<-getSYMBOL(rownames(temp.plot.data.3),data='org.Hs.eg')
# 
# # library(org.Hs.eg.db)
# # library(annotate)
# rownames()<-getSYMBOL(rownames(temp.plot.data.3),data='org.Hs.eg')
# #rev_names<-rev(rownames(temp.plot.data.3))
# #temp.plot.data.3.1 <- temp.plot.data.3[match(rev_names,rownames(temp.plot.data.3)),]

current.stouffer.erbb2 <- apply(plot.data.tf.cotf,1,function(x){
  y <- sum(x)/sqrt(length(x))
  #y <- sum(w*x)/sqrt(sum(w^2))
  return(y)
})

# create a synthetic regulon with the top 50 and bottom 50 proteins from the consensus protein activity signature
reference <- c(sort(current.stouffer.erbb2,decreasing = TRUE)[1:50],rev(sort(current.stouffer.erbb2,decreasing = FALSE)[1:50]))

# order the samples in the protein activity matrix by their enrichment and the rows by the protein activity in the consensus signature
protein.sort <- names(sort(reference,decreasing = TRUE))

# view results with a heatmap
temp.plot.data.3.1 <- plot.data.tf.cotf[match(protein.sort,rownames(plot.data.tf.cotf)),]

mat_breaks <- quantile_breaks(temp.plot.data.3.1,n = 100)
rownames(temp.plot.data.3.1)<-getSYMBOL(rownames(temp.plot.data.3.1),data='org.Hs.eg')

hm.4<-pheatmap(temp.plot.data.3.1,kmeans_k = NA,  color = colorRampPalette(rev(brewer.pal(n = 7, 
                                                                                          name = "RdBu")))
               (length(mat_breaks)), fontsize_row = 6, fontsize_col = 8,
               main = "Predicted erbb2 Neomorphs: brca",
               show_rownames = TRUE, show_colnames = TRUE, breaks = mat_breaks, 
               cluster_cols = TRUE, cluster_rows = FALSE)
save_pheatmap_pdf(hm.4,file="neomorph_comparison_erbb2_brca_topbottom_tf_cotf.pdf")

#pathway enrichment
#1. do for the top50 bottom50 neomorphic tf/cotf
mapid<-mapIds(org.Hs.eg.db,rownames(temp.plot.data.3.1),'ENTREZID','SYMBOL')
#2. do for the consensus signature
consensus.sig<-readRDS(file="consensus_signature_erbb2_brca.rds")
consensus.sig <- consensus.sig[,1]
reference <- c(sort(consensus.sig,decreasing = TRUE)[1:50],rev(sort(consensus.sig,decreasing = FALSE)[1:50]))

# order the samples in the protein activity matrix by their enrichment and the rows by the protein activity in the consensus signature
protein.sort.1 <- names(sort(reference,decreasing = TRUE))

mapid<-mapIds(org.Hs.eg.db,rownames(temp.plot.data.3.1),'ENTREZID','SYMBOL')
#mapid<-mapIds(org.Hs.eg.db,protein.sort.1,'ENTREZID','SYMBOL')

write.csv(mapid,file="mapid_neomorphs_erbb2_brca.csv")
mapid.1<-read.csv("mapid_neomorphs_erbb2_brca.csv",header=TRUE, sep=",")
mapid.1[,2]
mapid.1[,2] <- protein.sort.1

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ReactomePA")

x <- enrichPathway(gene=protein.sort.1,pvalueCutoff=0.05, readable=T)
y <- as.data.frame(x)
head(y)
gp<-ggplot(y, # you can replace the numbers to the row number of pathway of your interest
           aes(x = Count, y = Description, label=y$geneID)) + 
  geom_point(aes(size = Count, color = p.adjust)) +
  # geom_text(aes(label=y$geneID),hjust=0, vjust=0)+
  #  theme_bw(base_size = 10) +
  scale_colour_gradient(limits=c(0, 0.10), low="blue", high="red") +
  ylab(NULL) +
  #ggtitle("GO pathway enrichment: erbb2 (brca) Neomorphs Top50, Bottom50 TFs/CoTFs")
  ggtitle("GO pathway enrichment: erbb2 (brca) Consensus Signature Top50, Bottom50 TFs/CoTFs")
#  gp + 
#    geom_text_repel(
#                   size          = 3,
#                   direction     = "x")
# 
# +
#   

#plot1
dp.1<-gp+  geom_label_repel(aes(label=y$geneID),
                            box.padding   = 0.3, 
                            point.padding = 0.6,
                            size=3,
                            direction = "x")#,
#pdf(paste('GO_Enrich_dotplot_erbb2_brca_neomorph.pdf', sep=""),width = 15, height = 10, family = "Times", pointsize = 8)
pdf(paste('GO_Enrich_dotplot_erbb2_brca_consensus.pdf', sep=""),width = 15, height = 10, family = "Times", pointsize = 8)
print(dp.1)
dev.off()
gp+geom_label_repel(mapping = NULL, data = NULL, stat = "identity",
                    position = "identity", parse = FALSE,  box.padding = 0.25,
                    label.padding = 0.25, point.padding = 1e-06, label.r = 0.15,
                    label.size = 0.25, segment.colour = NULL, segment.color = NULL,
                    segment.size = 0.5, segment.alpha = NULL, min.segment.length = 0.5,
                    arrow = NULL, force = 1, max.iter = 2000, #nudge_x = 0, nudge_y = 0,
                    xlim = c(NA, NA), ylim = c(NA, NA), na.rm = FALSE, show.legend = NA,
                    direction = c("both", "y", "x"), seed = NA, inherit.aes = TRUE)

#  box.padding   = 1.5,
# point.padding = 0.5,
#  force         = 50,
#  segment.size  = 0.2)
#  segment.color = "grey50",
#  direction     = "x") 

gp<-ggplot(y, # you can replace the numbers to the row number of pathway of your interest
           aes(x = Count, y = Description, label=y$geneID)) + 
  geom_point(aes(size = Count, color = p.adjust)) 

gp + 
  geom_label_repel(aes(label=y$geneID),
                   box.padding   = 0.3, 
                   point.padding = 0.6,
                   segment.color = "red")

library("ReactomePA")
rpa <- enrichPathway(gene=protein.sort.1, organism = "human",
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "BH",
                     readable=T)
# rpa <- enrichPathway(gene=mapid.1[,2], organism = "human", pvalueCutoff = 0.05,
#               pAdjustMethod = "BH", qvalueCutoff = 0.2, universe, minGSSize = 10,
#               maxGSSize = 500, readable = T)
head(as.data.frame(rpa))
write.csv(rpa,file="pathway_enrich_erbb2_brca.csv")
barplot(rpa, showCategory=100)
dp<-dotplot(rpa, x = "Count",
            color = "p.adjust", showCategory = 100, split = NULL, font.size = 12,
            title = "rpa")
emapplot(rpa,showCategory = 10, color = "p.adjust", layout = "kk")
#plot2
test.1<-sort(unique(dp$data[,9]))
for(i in 1:length(test.1)){
  cnp<-cnetplot(rpa, showCategory = dp$data[,9]==test.1[i],categorySize="pvalue", foldChange=protein.sort.1,
                layout = "kk")
  pdf(paste('GO_Enrich_cnetplot_erbb2_brca_consensus_',test.1[i],'.pdf', sep=""),width = 15, height = 10, family = "Times", pointsize = 8)
  #pdf(paste('GO_Enrich_cnetplot_erbb2_brca_neomorph_',test.1[i],'.pdf', sep=""),width = 15, height = 10, family = "Times", pointsize = 8)
  print(cnp)
  dev.off()
}
cnetplot(rpa, categorySize="pvalue", foldChange=mapid.1[,2],
         layout = "kk")
# cnetplot.enrichResult(rpa, showCategory = 50,categorySize="pvalue", foldChange=mapid.1[,2],
#                       layout = "kk",node_label = "all")

library(DOSE)

edo <- enrichDGN(mapid.1[,2])
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')

p1 <- heatplot(edox)
p2 <- heatplot(edox, foldChange=mapid.1[,2])
cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])

p1 <- cnetplot(edox, foldChange=mapid.1[,2])
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=mapid.1[,2])
p3 <- cnetplot(edox, foldChange=mapid.1[,2], circular = TRUE, colorEdge = TRUE)
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))


p1 <- cnetplot(edox, node_label="category") 
p2 <- cnetplot(edox, node_label="gene") 
p3 <- cnetplot(edox, node_label="all") 
p4 <- cnetplot(edox, node_label="none") 
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])

library(clusterProfiler)
ridgeplot(edo)
####ggplot dotplot

library(UpSetR)
upset(edo)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("clusterProfiler")
# library("clusterProfiler")
# res <- compareCluster(mapid.1[,2], fun="enrichKEGG",organism="hsa", pvalueCutoff=0.05)
# dotplot(res)

