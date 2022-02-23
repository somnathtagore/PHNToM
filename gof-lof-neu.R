#Define GOF/LOF signature

#ERBB2 in BRCA as a case study
save_pheatmap_pdf <- function(x, filename, width=10, height=10) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
library(pheatmap)
library(RColorBrewer)

quantile_breaks <- function(xs, n = 10){
  breaks <- quantile(xs, probs = seq(from = 0, to = 1, length.out = n))
  unique(breaks)
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
library(viper)

# tcga data
brca.tumor<-read.table(file="../../batch_corrected_data/data/batch-correct/tcga-tumor/brca-tcga-tumor.txt", header=TRUE)
dim(brca.tumor)

# tcga normal data
brca.normal<-read.table(file="../../batch_corrected_data/data/batch-correct/tcga-normal/brca-tcga-normal.txt", header=TRUE)
dim(brca.normal)

# # n1data
# load("brca-tissue_rawcounts.rda")
# brca.samples<-rawcounts
# dim(brca.samples)

#normal and tumor samples
# brca.normal <- brca.samples[is.finite(rowSums(brca.samples)),grep('TCGA-\\w+-\\w+-11',colnames(brca.samples))] #wt
# brca.tumor <- brca.samples[is.finite(rowSums(brca.samples)),grep('TCGA-\\w+-\\w+-01',colnames(brca.samples))] #wt

# head(brca.normal)
# head(brca.tumor)

tcga.normal<-brca.normal
tcga.tumor<-brca.tumor
#removing #N/A from samples
tcga.normal <- tcga.normal[is.finite(rowSums(tcga.normal)),]
tcga.tumor <- tcga.tumor[is.finite(rowSums(tcga.tumor)),]

any(is.na(tcga.normal))
any(is.na(tcga.tumor))

#check for identical genes
identical(rownames(tcga.tumor),rownames(tcga.normal))
table(rownames(tcga.tumor)%in%rownames(tcga.normal))
common_genes <- intersect(rownames(tcga.tumor),rownames(tcga.normal))
length(common_genes)

tcga.tumor <- tcga.tumor[match(common_genes,rownames(tcga.tumor)),]
tcga.normal <- tcga.normal[match(common_genes,rownames(tcga.normal)),]
identical(rownames(tcga.normal),rownames(tcga.tumor))

mixed.tcga.tumor.normal <- cbind(tcga.tumor, tcga.normal)
dim(mixed.tcga.tumor.normal)

#brca.samples <- brca.samples[is.finite(rowSums(brca.samples)),]
any(is.na(mixed.tcga.tumor.normal))

#log2cpm+1 normalization
mixed.tcga.tumor.normal.log2cpm <- apply(mixed.tcga.tumor.normal,2,function(x){
  y <- 1E6*x/sum(x) + 1
  z <- log(y,2)
  return(z)
})
# 
mixed.tcga.tumor.normal.log2cpm.ges <- t(apply(mixed.tcga.tumor.normal.log2cpm,1,function(x){
  y <- (x - mean(x))/sd(x)
  return(y)
}))
#tpm normalization
# for(i in 1:ncol(mixed.tcga.tumor.normal)){
#   mixed.tcga.tumor.normal[,i] <- 1E6*mixed.tcga.tumor.normal[,i]/sum(mixed.tcga.tumor.normal[,i])
#  }
#mixed.tcga.tumor.normal<-mixed.tcga.tumor.normal.log2cpm.ges
#brca.samples <- brca.samples[is.finite(rowSums(brca.samples)),]

#viper without null
load('../../brca-gtex/pik3ca_brca/brca_regulon.rda')

pregul.brca <- pruneRegulon(regul, cutoff = 50)
mixed.tcga.tumor.normal.ges.vpmat <- viper(eset = mixed.tcga.tumor.normal.log2cpm.ges, regulon = pregul.brca, method = "none")
# vpsig_internal <- viperSignature(eset = brca.samples, ref = brca.samples, per = 1000, method = "ttest", verbose = TRUE)
# mixed.tcga.tumor.normal.ges.vpmat <- viper(vpsig_internal, regulon = pregul.brca)

brca.mixed.tcga.tumor.normal.int.withoutNull <- mixed.tcga.tumor.normal.ges.vpmat
saveRDS(brca.mixed.tcga.tumor.normal.int.withoutNull,file="brca.mixed.tcga.tumor.normal.int.withoutNull.rds")
plot(density(brca.mixed.tcga.tumor.normal.int.withoutNull))

brca.mixed.tcga.tumor.normal.int.withoutNull<-readRDS(file="./pik3ca_brca/brca.mixed.tcga.tumor.normal.int.withoutNull.rds")
#tumor.dist <- as.dist(viperSimilarity(tumor.vpmat,method = "two.sided"))

library(fpc)
library(cluster)

#tumor.clust <- pamk(tumor.dist)$pamobject

#plot(silhouette(tumor.clust))

#any(is.na(vipermat))
tumor.vpmat<-brca.mixed.tcga.tumor.normal.int.withoutNull
dim(tumor.vpmat)
any(is.na(tumor.vpmat))

#viper similarity
vipsim <- viperSimilarity(tumor.vpmat)
dist.matr <- as.dist(vipsim)

#pamk clustering
library(fpc)
pam.clust<-pamk(dist.matr,krange=2:10,criterion="asw", usepam=TRUE,
                scaling=FALSE, alpha=0.001, diss=inherits(data, "dist"),
                critout=FALSE, ns=10, seed=NULL)
pam.clust$nc

current.silinfo <- as.data.frame(pam.clust$pamobject$silinfo$widths)

library(cluster)
current.bootsil <- lapply(1:1000,function(x,ref.clust,ref.dist){
  y <- sample(x = ref.clust, size = length(ref.clust), replace = FALSE)
  z <- silhouette(x = y, dist = dist.matr)[,3]
  return(z)
}, ref.clust = as.numeric(pam.clust$pamobject$clustering), ref.dist = dist.matr)
current.bootsil <- unlist(current.bootsil)

# compute true silhouette score right-tailed p-values
current.silinfo$pvalue <- pnorm(q = current.silinfo$sil_width, 
                                mean = mean(current.bootsil), sd = sd(current.bootsil), 
                                lower.tail = FALSE) - pnorm(q = 1, 
                                                            mean = mean(current.bootsil), 
                                                            sd = sd(current.bootsil), 
                                                            lower.tail = FALSE)

# adjust p-values for FDR 
current.silinfo$padj <- p.adjust(current.silinfo$pvalue, method = "BH")

# assign samples with FDR-adjusted p-values > 0.05 to an outlier cluster 0
# current.silinfo[which(current.silinfo$padj > 0.05),"cluster"] <- 0

# order the silhouette score info by cluster and silhoeutte score
current.silinfo <- current.silinfo[order(current.silinfo$cluster, current.silinfo$sil_width, decreasing = TRUE),]

# plot distance matrix after rearranging the results by silhouette score and cluster and annotation responders and non-responders

annot.col <- current.silinfo[,c("cluster","sil_width")]
colnames(annot.col) <- c("Cluster","SilScore")
annot.col$Cluster <- factor(annot.col$Cluster, 
                            levels = sort(unique(annot.col$Cluster), decreasing = TRUE))
annot.col <- annot.col[,c("SilScore","Cluster")]

#annot.col.coad<-annot.col
#saveRDS(annot.col.coad,file="annot.col.coad.rds")
#vipermat <- vipermat[is.finite(rowSums(vipermat)),]
#plot.data<-vipermat
#rownames(annot.col)<-gsub('^.{12}','',rownames(annot.col))
rownames(annot.col)
plot.data <- tumor.vpmat[,match(rownames(annot.col),colnames(tumor.vpmat))] 

any(is.na(plot.data))

plot.data <- plot.data[,is.finite(colSums(plot.data))]
dim(plot.data)
plot.data[is.na(plot.data)] <- 0

mat_breaks <- quantile_breaks(plot.data,n = 100)

test1<-pheatmap(plot.data, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))
                (length(mat_breaks)), 
                fontsize_row = 5, fontsize_col = 5, 
                #main = "BRAF mutation cluster COAD Test 1",
                main = "BRCA Unsupervised\n Internal",
                show_rownames = FALSE, show_colnames = FALSE, breaks = mat_breaks, 
                cluster_cols = FALSE, cluster_rows = TRUE, annotation_col = annot.col)
#save_pheatmap_pdf(test1,"BRAF_mutation_cluster_COAD_test1.pdf")
save_pheatmap_pdf(test1,"BRCA_Unsupervised_Internal.pdf")
saveRDS(plot.data,"BRCA.Unsup.Internal.rds")

keep.index <- apply(plot.data,1,function(x,ref){
  y <- pairwise.t.test(x = x,g = ref)
  z <- min(y$p.value, na.rm = TRUE)
  return(z)
}, ref = annot.col$Cluster)


keep.num <- length(rownames(plot.data))

final.keep.index <- sort(keep.index, decreasing = FALSE)[1:keep.num]

#rownames(plot.data)<-getSYMBOL(rownames(plot.data),data='org.Hs.eg')
final.plot.data <- plot.data[match(names(final.keep.index),rownames(plot.data)),]
dim(final.plot.data)
rownames(final.plot.data)

# #mat_breaks <- quantile_breaks(final.plot.data,n = 100)
# test2<-pheatmap(final.plot.data, color = colorRampPalette(rev(brewer.pal(n = 7, 
#                                                                          name = "RdBu")))
#                 (length(mat_breaks)), fontsize_row = 5, fontsize_col = 5,
#                 #main = "BRAF mutation cluster COAD Test 2",
#                 main = "LUAD Unsupervised (+GTEx) test2",
#                 show_rownames = FALSE, show_colnames = FALSE, breaks = mat_breaks, 
#                 cluster_cols = FALSE, cluster_rows = TRUE, annotation_col = annot.col)
# #save_pheatmap_pdf(test2,"BRAF_mutation_cluster_COAD_test2.pdf")
# save_pheatmap_pdf(test2,"LUAD_Unsupervised_GTEx_test2.pdf")
# dim(final.plot.data)
# any(is.na(final.plot.data))

tf.cotf <- read.table(file="~/Desktop/tf-cotf-homo-current.txt")
row.names(tf.cotf) <- tf.cotf[,1]
tf.cotf<-tf.cotf[,-1]
rownames(tf.cotf)

#rownames(tf.cotf)<-getSYMBOL(rownames(tf.cotf),data='org.Hs.eg')
common.tf.cotf <- intersect(rownames(final.plot.data),rownames(tf.cotf))
length(common.tf.cotf)

plot.data.tf.cotf <- final.plot.data[match(common.tf.cotf,rownames(final.plot.data)),]
dim(plot.data.tf.cotf)
#saveRDS(plot.data.tf.cotf,file="braf.coad.plot.data.rds")
#mat_breaks <- quantile_breaks(plot.data.tf.cotf,n = 100)

plot.data.tf.cotf <- plot.data.tf.cotf[is.finite(rowSums(plot.data.tf.cotf)),]
test3<-pheatmap(plot.data.tf.cotf, color = colorRampPalette(rev(brewer.pal(n = 7, 
                                                                           name = "RdBu")))
                (length(mat_breaks)), fontsize_row = 5, fontsize_col = 5,
                #main = "BRAF mutation cluster COAD TF/Co-TF",
                main = "BRCA Unsupervised\n TF/Co-TF (+Internal)",
                show_rownames = FALSE, show_colnames = FALSE, breaks = mat_breaks, 
                cluster_cols = FALSE, cluster_rows = TRUE, annotation_col = annot.col)
#save_pheatmap_pdf(test3,"BRAF_mutation_cluster_COAD_tf_cotfs.pdf")
save_pheatmap_pdf(test3,"BRCA_tf_cotfs_Internal.pdf")
saveRDS(plot.data.tf.cotf,"BRCA.Unsup.all.tf.cotf.Internal.rds")
dim(plot.data.tf.cotf)

zscore_data_egfr <- t(apply(plot.data.tf.cotf,1,function(x){
  y <- (x - mean(x))/sd(x)
  return(y)
}))

current.stouffer <- t(apply(zscore_data_egfr,1,function(x,clust,weight){
  temp.data <- data.frame(nes = x, cluster = clust, weight = weight)
  temp.data$w_z <- temp.data$nes*temp.data$weight
  temp.data$w_w <- temp.data$weight*temp.data$weight
  temp.agg <- aggregate(.~cluster, temp.data,sum)
  temp.agg$res <- temp.agg$w_z/sqrt(temp.agg$w_w)
  temp.res <- temp.agg$res
  return(temp.res)
}, clust = annot.col$Cluster, weight = annot.col$SilScore))
colnames(current.stouffer) <- levels(annot.col$Cluster)

rank.stouffer <- apply(current.stouffer,2,rank)

keep.num <- 2000
keep.index <- apply(rank.stouffer,1,function(x){
  y <- ((max(x) > (nrow(rank.stouffer) + 1 - keep.num)) | min(x) <= keep.num)
  return(y)
})

temp.plot.data.2 <- plot.data.tf.cotf[keep.index,]

mat_breaks <- quantile_breaks(temp.plot.data.2,n = 100)
library(org.Hs.eg.db)
library(annotate)
rownames(temp.plot.data.2)<-getSYMBOL(rownames(temp.plot.data.2),data='org.Hs.eg')

test3<-pheatmap(temp.plot.data.2, color = colorRampPalette(rev(brewer.pal(n = 7, 
                                                                          name = "RdBu")))
                (length(mat_breaks)), fontsize_row = 4, fontsize_col = 5,
                main = "BRCA Unsupervised TF/Co-TF \n(+Internal).All",
                show_rownames = TRUE, show_colnames = FALSE, breaks = mat_breaks, 
                cluster_cols = FALSE, cluster_rows = TRUE, annotation_col = annot.col)
save_pheatmap_pdf(test3,"BRCA_Unsupervised_TF_Co-TF_Int.All.pdf")
saveRDS(temp.plot.data.2,"BRCA.Unsup.all.tf.cotf.Internal.1.rds")


var.keep.index <- apply(plot.data.tf.cotf,1,sd)
var.keep.index <- (rank(var.keep.index) > length(var.keep.index) - 100)

temp.plot.data.3 <- plot.data.tf.cotf[var.keep.index,]

mat_breaks <- quantile_breaks(temp.plot.data.3,n = 100)
#rownames(temp.plot.data.3)<-getSYMBOL(rownames(temp.plot.data.3),data='org.Hs.eg')

# library(org.Hs.eg.db)
# library(annotate)
rownames(temp.plot.data.3)<-getSYMBOL(rownames(temp.plot.data.3),data='org.Hs.eg')

test3<-pheatmap(temp.plot.data.3, color = colorRampPalette(rev(brewer.pal(n = 7, 
                                                                          name = "RdBu")))
                (length(mat_breaks)), fontsize_row = 4, fontsize_col = 5,
                main = "BRCA Unsupervised TF/Co-TF \n(+Internal) Top100",
                show_rownames = TRUE, show_colnames = FALSE, breaks = mat_breaks, 
                cluster_cols = FALSE, cluster_rows = TRUE, annotation_col = annot.col)
save_pheatmap_pdf(test3,"BRCA_Unsupervised_TF_Co-TF_Int.top100.pdf")
saveRDS(temp.plot.data.3,"BRCA.Unsup.top100.tf.cotf.Int.rds")

#Cluster-specific semi-supervised analysis

#samples in clusters
colnames(brca.tumor)<-gsub('.{16}$','',colnames(brca.tumor))


brca.samples.cl.1<-rownames(annot.col[which(annot.col$Cluster=='1'),])#273
#brca.samples.cl.1<-gsub('.{16}$','',brca.samples.cl.1)
brca.samples.cl.2<-rownames(annot.col[which(annot.col$Cluster=='2'),])#242
#brca.samples.cl.2<-gsub('.{16}$','',brca.samples.cl.2)
brca.samples.cl.3<-rownames(annot.col[which(annot.col$Cluster=='3'),])#244
#brca.samples.cl.3<-gsub('.{16}$','',brca.samples.cl.3)
brca.samples.cl.4<-rownames(annot.col[which(annot.col$Cluster=='4'),])#223
#brca.samples.cl.4<-gsub('.{16}$','',brca.samples.cl.4)
#colnames(mixed.tcga.tumor.normal)<-gsub('.{16}$','',colnames(mixed.tcga.tumor.normal))

brca.mixed.cl.1 <- 
  mixed.tcga.tumor.normal[,intersect(brca.samples.cl.1,colnames(mixed.tcga.tumor.normal))]
dim(brca.mixed.cl.1)
saveRDS(brca.mixed.cl.1,file = "brca.mixed.cl.1.rds")

brca.mixed.cl.2 <- 
  mixed.tcga.tumor.normal[,intersect(brca.samples.cl.2,colnames(mixed.tcga.tumor.normal))]
dim(brca.mixed.cl.2)
saveRDS(brca.mixed.cl.2,file = "brca.mixed.cl.2.rds")

brca.mixed.cl.3 <- 
  mixed.tcga.tumor.normal[,intersect(brca.samples.cl.3,colnames(mixed.tcga.tumor.normal))]
dim(brca.mixed.cl.3)
saveRDS(brca.mixed.cl.3,file = "brca.mixed.cl.3.rds")

brca.mixed.cl.4 <- 
  mixed.tcga.tumor.normal[,intersect(brca.samples.cl.4,colnames(mixed.tcga.tumor.normal))]
dim(brca.mixed.cl.4)
saveRDS(brca.mixed.cl.4,file = "brca.mixed.cl.4.rds")

#match samples in clusters with samples with mutations

outersect <- function(x,y){
  sort(c(setdiff(x,y),setdiff(y,x)))
}
annot.1<-rownames(annot.col[which(annot.col$Cluster=='1'),])
#annot.1<-gsub('.{16}$','',annot.1)
annot.2<-rownames(annot.col[which(annot.col$Cluster=='2'),])
#annot.2<-gsub('.{16}$','',annot.2)
annot.3<-rownames(annot.col[which(annot.col$Cluster=='3'),])
#annot.3<-gsub('.{16}$','',annot.3)
annot.4<-rownames(annot.col[which(annot.col$Cluster=='4'),])
#annot.4<-gsub('.{16}$','',annot.4)

brca.mixed.cl.1<-readRDS(file = "../pik3ca_brca/rds/brca.mixed.cl.1.rds")
brca.mixed.cl.2<-readRDS(file = "../pik3ca_brca/rds/brca.mixed.cl.2.rds")
brca.mixed.cl.3<-readRDS(file = "../pik3ca_brca/rds/brca.mixed.cl.3.rds")
brca.mixed.cl.4<-readRDS(file = "../pik3ca_brca/rds/brca.mixed.cl.4.rds")

erbb2.brca.Tumor_Sample_Barcode<-read.csv("~/Documents/Mutation_mapper/brca_mutations/erbb2_brca_mutations.csv",header=TRUE, sep=",")
rownames(erbb2.brca.Tumor_Sample_Barcode) <- erbb2.brca.Tumor_Sample_Barcode[,1]
erbb2.brca.Tumor_Sample_Barcode<-erbb2.brca.Tumor_Sample_Barcode[,-1]
rownames(erbb2.brca.Tumor_Sample_Barcode)

rownames(erbb2.brca.Tumor_Sample_Barcode)<-gsub("\\-","\\.",rownames(erbb2.brca.Tumor_Sample_Barcode))
rownames(erbb2.brca.Tumor_Sample_Barcode)<-gsub('.{3}$','',rownames(erbb2.brca.Tumor_Sample_Barcode))
colnames(brca.tumor)<-gsub('.{16}$','',colnames(brca.tumor))
erbb2.brca.samples<-brca.tumor[,intersect(rownames(erbb2.brca.Tumor_Sample_Barcode),colnames(brca.tumor))]
dim(erbb2.brca.samples)
saveRDS(erbb2.brca.samples,file="erbb2.brca.samples.rds")

# erbb2.brca.cl.1<- 
#   vipermat.brca.cl.1.wt.1[,intersect(rownames(erbb2.brca.Tumor_Sample_Barcode),colnames(vipermat.brca.cl.1.wt.1))]
# dim(erbb2.brca.cl.1)
brca.samples.cl.1<-colnames(brca.mixed.cl.1)
brca.samples.cl.1<-gsub('.{16}$','',brca.samples.cl.1)
brca.tum.cl.1 <- 
  brca.tumor[,intersect(brca.samples.cl.1,colnames(brca.tumor))]
dim(brca.tum.cl.1)
colnames(brca.tum.cl.1)<-gsub('.{16}$','',colnames(brca.tum.cl.1))
erbb2.brca.tum.cl.1 <- 
  brca.tum.cl.1[,intersect(colnames(brca.tum.cl.1),rownames(erbb2.brca.Tumor_Sample_Barcode))]
dim(erbb2.brca.tum.cl.1)
erbb2.brca.wt.cl.1 <- 
  brca.tum.cl.1[,outersect(colnames(brca.tum.cl.1),colnames(erbb2.brca.tum.cl.1))]
dim(erbb2.brca.wt.cl.1)
saveRDS(erbb2.brca.wt.cl.1,file = "erbb2.brca.wt.cl.1.rds")
saveRDS(erbb2.brca.tum.cl.1,file = "erbb2.brca.tumor.cl.1.rds")

brca.samples.cl.2<-colnames(brca.mixed.cl.2)
brca.samples.cl.2<-gsub('.{16}$','',brca.samples.cl.2)
brca.tum.cl.2 <- 
  brca.tumor[,intersect(brca.samples.cl.2,colnames(brca.tumor))]
dim(brca.tum.cl.2)
colnames(brca.tum.cl.2)<-gsub('.{16}$','',colnames(brca.tum.cl.2))
erbb2.brca.tum.cl.2 <- 
  brca.tum.cl.2[,intersect(colnames(brca.tum.cl.2),rownames(erbb2.brca.Tumor_Sample_Barcode))]
dim(erbb2.brca.tum.cl.2)
erbb2.brca.wt.cl.2 <- 
  brca.tum.cl.2[,outersect(colnames(brca.tum.cl.2),colnames(erbb2.brca.tum.cl.2))]
dim(erbb2.brca.wt.cl.2)
saveRDS(erbb2.brca.wt.cl.2,file = "erbb2.brca.wt.cl.2.rds")
saveRDS(erbb2.brca.tum.cl.2,file = "erbb2.brca.tumor.cl.2.rds")

brca.samples.cl.3<-colnames(brca.mixed.cl.3)
brca.samples.cl.3<-gsub('.{16}$','',brca.samples.cl.3)
brca.tum.cl.3 <- 
  brca.tumor[,intersect(brca.samples.cl.3,colnames(brca.tumor))]
dim(brca.tum.cl.3)
colnames(brca.tum.cl.3)<-gsub('.{16}$','',colnames(brca.tum.cl.3))
erbb2.brca.tum.cl.3 <- 
  brca.tum.cl.3[,intersect(colnames(brca.tum.cl.3),rownames(erbb2.brca.Tumor_Sample_Barcode))]
dim(erbb2.brca.tum.cl.3)
erbb2.brca.wt.cl.3 <- 
  brca.tum.cl.3[,outersect(colnames(brca.tum.cl.3),colnames(erbb2.brca.tum.cl.3))]
dim(erbb2.brca.wt.cl.3)
saveRDS(erbb2.brca.wt.cl.3,file = "erbb2.brca.wt.cl.3.rds")
saveRDS(erbb2.brca.tum.cl.3,file = "erbb2.brca.tumor.cl.3.rds")

brca.samples.cl.4<-colnames(brca.mixed.cl.4)
brca.samples.cl.4<-gsub('.{16}$','',brca.samples.cl.4)
brca.tum.cl.4 <- 
  brca.tumor[,intersect(brca.samples.cl.4,colnames(brca.tumor))]
dim(brca.tum.cl.4)
colnames(brca.tum.cl.4)<-gsub('.{16}$','',colnames(brca.tum.cl.4))
erbb2.brca.tum.cl.4 <- 
  brca.tum.cl.4[,intersect(colnames(brca.tum.cl.4),rownames(erbb2.brca.Tumor_Sample_Barcode))]
dim(erbb2.brca.tum.cl.4)
erbb2.brca.wt.cl.4 <- 
  brca.tum.cl.4[,outersect(colnames(brca.tum.cl.4),colnames(erbb2.brca.tum.cl.4))]
dim(erbb2.brca.wt.cl.4)
saveRDS(erbb2.brca.wt.cl.4,file = "erbb2.brca.wt.cl.4.rds")
saveRDS(erbb2.brca.tum.cl.4,file = "erbb2.brca.tumor.cl.4.rds")

#cluster specific protein activity to identify erbb2.wt closest to gtex breast
#cluster 1

dim(erbb2.brca.wt.cl.1)
dim(erbb2.brca.tum.cl.1)
breast_gtex<-read.table(file="~/Documents/viper-docs/benchmarking_1/signature_with_gtex/batch_corrected_data/data/batch-correct/gtex/breast-gtex.txt", header=TRUE)
brca.normal<-breast_gtex
dim(brca.normal)

identical(rownames(erbb2.brca.tum.cl.1),rownames(erbb2.brca.wt.cl.1))
table(rownames(erbb2.brca.tum.cl.1)%in%rownames(erbb2.brca.wt.cl.1))
common_genes <- intersect(rownames(erbb2.brca.tum.cl.1),rownames(erbb2.brca.wt.cl.1))
length(common_genes)

erbb2.brca.tum.cl.1 <- erbb2.brca.tum.cl.1[match(common_genes,rownames(erbb2.brca.tum.cl.1)),]
erbb2.brca.wt.cl.1 <- erbb2.brca.wt.cl.1[match(common_genes,rownames(erbb2.brca.wt.cl.1)),]
identical(rownames(erbb2.brca.wt.cl.1),rownames(erbb2.brca.tum.cl.1))

identical(rownames(brca.normal),rownames(erbb2.brca.wt.cl.1))
table(rownames(brca.normal)%in%rownames(erbb2.brca.wt.cl.1))
common_genes <- intersect(rownames(brca.normal),rownames(erbb2.brca.wt.cl.1))
length(common_genes)

brca.normal <- brca.normal[match(common_genes,rownames(brca.normal)),]
erbb2.brca.wt.cl.1 <- erbb2.brca.wt.cl.1[match(common_genes,rownames(erbb2.brca.wt.cl.1)),]
identical(rownames(erbb2.brca.wt.cl.1),rownames(brca.normal))

mixed.erbb2.brca.cl.1 <- cbind(brca.normal, erbb2.brca.wt.cl.1,erbb2.brca.tum.cl.1)
dim(mixed.erbb2.brca.cl.1)

#log2cpm+1 normalization
mixed.erbb2.brca.cl.1.log2cpm <- apply(mixed.erbb2.brca.cl.1,2,function(x){
  y <- 1E6*x/sum(x) + 1
  z <- log(y,2)
  return(z)
})

mixed.erbb2.brca.cl.1.log2cpm.ges <- t(apply(mixed.erbb2.brca.cl.1.log2cpm,1,function(x){
  y <- (x - mean(x))/sd(x)
  return(y)
}))
brca.normal.log2cpm <- apply(brca.normal,2,function(x){
  y <- 1E6*x/sum(x) + 1
  z <- log(y,2)
  return(z)
})
brca.normal.log2cpm.ges <- t(apply(brca.normal.log2cpm,1,function(x){
  y <- (x - mean(x))/sd(x)
  return(y)
}))
# for(i in 1:ncol(mixed.brca.cl.2)){
#   mixed.brca.cl.2[,i] <- 1E6*mixed.brca.cl.2[,i]/sum(mixed.brca.cl.2[,i])
# }
# 
# for(i in 1:ncol(brca.normal)){
#   brca.normal[,i] <- 1E6*brca.normal[,i]/sum(brca.normal[,i])
# }

#mixed.erbb2.brca.cl.1.log2cpm.ges <- mixed.erbb2.brca.cl.1.log2cpm.ges[is.finite(rowSums(mixed.erbb2.brca.cl.1.log2cpm.ges)),]
#brca.normal.log2cpm.ges <- brca.normal.log2cpm.ges[is.finite(rowSums(brca.normal.log2cpm.ges)),]

#viper using gtex
load('../pik3ca_brca/brca_regulon.rda')
pregul.brca<-regul
vipsig.mixed.erbb2.brca.cl.1.gtex <- viperSignature(eset = as.matrix(mixed.erbb2.brca.cl.1.log2cpm.ges), 
                                                    ref = as.matrix(brca.normal.log2cpm.ges), per = 1000, method = "ttest", 
                                                    verbose = TRUE)

vipermat.mixed.erbb2.brca.cl.1.gtex <- viper(vipsig.mixed.erbb2.brca.cl.1.gtex, pregul.brca,minsize = 25)
saveRDS(vipermat.mixed.erbb2.brca.cl.1.gtex,file="vipermat.mixed.erbb2.brca.cl.1.gtex.rds")
plot(density(vipermat.mixed.erbb2.brca.cl.1.gtex))

#cluster 2
dim(erbb2.brca.wt.cl.2)
dim(erbb2.brca.tum.cl.2)
breast_gtex<-read.table(file="~/Documents/viper-docs/benchmarking_1/signature_with_gtex/batch_corrected_data/data/batch-correct/gtex/breast-gtex.txt", header=TRUE)
brca.normal<-breast_gtex

identical(rownames(erbb2.brca.tum.cl.2),rownames(erbb2.brca.wt.cl.2))
table(rownames(erbb2.brca.tum.cl.2)%in%rownames(erbb2.brca.wt.cl.2))
common_genes <- intersect(rownames(erbb2.brca.tum.cl.2),rownames(erbb2.brca.wt.cl.2))
length(common_genes)

erbb2.brca.tum.cl.2 <- erbb2.brca.tum.cl.2[match(common_genes,rownames(erbb2.brca.tum.cl.2)),]
erbb2.brca.wt.cl.2 <- erbb2.brca.wt.cl.2[match(common_genes,rownames(erbb2.brca.wt.cl.2)),]
identical(rownames(erbb2.brca.wt.cl.2),rownames(erbb2.brca.tum.cl.2))

identical(rownames(brca.normal),rownames(erbb2.brca.wt.cl.2))
table(rownames(brca.normal)%in%rownames(erbb2.brca.wt.cl.2))
common_genes <- intersect(rownames(brca.normal),rownames(erbb2.brca.wt.cl.2))
length(common_genes)

brca.normal <- brca.normal[match(common_genes,rownames(brca.normal)),]
erbb2.brca.wt.cl.2 <- erbb2.brca.wt.cl.2[match(common_genes,rownames(erbb2.brca.wt.cl.2)),]
identical(rownames(erbb2.brca.wt.cl.2),rownames(brca.normal))

mixed.erbb2.brca.cl.2 <- cbind(brca.normal, erbb2.brca.wt.cl.2,erbb2.brca.tum.cl.2)
dim(mixed.erbb2.brca.cl.2)

#log2cpm+1 normalization
mixed.erbb2.brca.cl.2.log2cpm <- apply(mixed.erbb2.brca.cl.2,2,function(x){
  y <- 1E6*x/sum(x) + 1
  z <- log(y,2)
  return(z)
})

mixed.erbb2.brca.cl.2.log2cpm.ges <- t(apply(mixed.erbb2.brca.cl.2.log2cpm,1,function(x){
  y <- (x - mean(x))/sd(x)
  return(y)
}))
brca.normal.log2cpm <- apply(brca.normal,2,function(x){
  y <- 1E6*x/sum(x) + 1
  z <- log(y,2)
  return(z)
})
brca.normal.log2cpm.ges <- t(apply(brca.normal.log2cpm,1,function(x){
  y <- (x - mean(x))/sd(x)
  return(y)
}))
# for(i in 1:ncol(mixed.brca.cl.2)){
#   mixed.brca.cl.2[,i] <- 1E6*mixed.brca.cl.2[,i]/sum(mixed.brca.cl.2[,i])
# }
# 
# for(i in 1:ncol(brca.normal)){
#   brca.normal[,i] <- 1E6*brca.normal[,i]/sum(brca.normal[,i])
# }

# mixed.erbb2.brca.cl.1.log2cpm.ges <- mixed.erbb2.brca.cl.1.log2cpm.ges[is.finite(rowSums(mixed.erbb2.brca.cl.1.log2cpm.ges)),]
# brca.normal.log2cpm.ges <- brca.normal.log2cpm.ges[is.finite(rowSums(brca.normal.log2cpm.ges)),]

#viper using gtex
vipsig.mixed.erbb2.brca.cl.2.gtex <- viperSignature(eset = as.matrix(mixed.erbb2.brca.cl.2.log2cpm.ges), 
                                                    ref = as.matrix(brca.normal.log2cpm.ges), per = 1000, method = "ttest", 
                                                    verbose = TRUE)

vipermat.mixed.erbb2.brca.cl.2.gtex <- viper(vipsig.mixed.erbb2.brca.cl.2.gtex, pregul.brca,minsize = 25)
saveRDS(vipermat.mixed.erbb2.brca.cl.2.gtex,file="vipermat.mixed.erbb2.brca.cl.2.gtex.rds")
plot(density(vipermat.mixed.erbb2.brca.cl.2.gtex))


#cluster 3
dim(erbb2.brca.wt.cl.3)
dim(erbb2.brca.tum.cl.3)
breast_gtex<-read.table(file="~/Documents/viper-docs/benchmarking_1/signature_with_gtex/batch_corrected_data/data/batch-correct/gtex/breast-gtex.txt", header=TRUE)
brca.normal<-breast_gtex

identical(rownames(erbb2.brca.tum.cl.3),rownames(erbb2.brca.wt.cl.3))
table(rownames(erbb2.brca.tum.cl.3)%in%rownames(erbb2.brca.wt.cl.3))
common_genes <- intersect(rownames(erbb2.brca.tum.cl.3),rownames(erbb2.brca.wt.cl.3))
length(common_genes)

erbb2.brca.tum.cl.3 <- erbb2.brca.tum.cl.3[match(common_genes,rownames(erbb2.brca.tum.cl.3)),]
erbb2.brca.wt.cl.3 <- erbb2.brca.wt.cl.3[match(common_genes,rownames(erbb2.brca.wt.cl.3)),]
identical(rownames(erbb2.brca.wt.cl.3),rownames(erbb2.brca.tum.cl.3))

identical(rownames(brca.normal),rownames(erbb2.brca.wt.cl.3))
table(rownames(brca.normal)%in%rownames(erbb2.brca.wt.cl.3))
common_genes <- intersect(rownames(brca.normal),rownames(erbb2.brca.wt.cl.3))
length(common_genes)

brca.normal <- brca.normal[match(common_genes,rownames(brca.normal)),]
erbb2.brca.wt.cl.3 <- erbb2.brca.wt.cl.3[match(common_genes,rownames(erbb2.brca.wt.cl.3)),]
identical(rownames(erbb2.brca.wt.cl.3),rownames(brca.normal))

mixed.erbb2.brca.cl.3 <- cbind(brca.normal, erbb2.brca.wt.cl.3,erbb2.brca.tum.cl.3)
dim(mixed.erbb2.brca.cl.3)

#log2cpm+1 normalization
mixed.erbb2.brca.cl.3.log2cpm <- apply(mixed.erbb2.brca.cl.3,2,function(x){
  y <- 1E6*x/sum(x) + 1
  z <- log(y,2)
  return(z)
})

mixed.erbb2.brca.cl.3.log2cpm.ges <- t(apply(mixed.erbb2.brca.cl.3.log2cpm,1,function(x){
  y <- (x - mean(x))/sd(x)
  return(y)
}))
brca.normal.log2cpm <- apply(brca.normal,2,function(x){
  y <- 1E6*x/sum(x) + 1
  z <- log(y,2)
  return(z)
})
brca.normal.log2cpm.ges <- t(apply(brca.normal.log2cpm,1,function(x){
  y <- (x - mean(x))/sd(x)
  return(y)
}))
# for(i in 1:ncol(mixed.brca.cl.2)){
#   mixed.brca.cl.2[,i] <- 1E6*mixed.brca.cl.2[,i]/sum(mixed.brca.cl.2[,i])
# }
# 
# for(i in 1:ncol(brca.normal)){
#   brca.normal[,i] <- 1E6*brca.normal[,i]/sum(brca.normal[,i])
# }

# mixed.erbb2.brca.cl.1.log2cpm.ges <- mixed.erbb2.brca.cl.1.log2cpm.ges[is.finite(rowSums(mixed.erbb2.brca.cl.1.log2cpm.ges)),]
# brca.normal.log2cpm.ges <- brca.normal.log2cpm.ges[is.finite(rowSums(brca.normal.log2cpm.ges)),]

#viper using gtex
vipsig.mixed.erbb2.brca.cl.3.gtex <- viperSignature(eset = as.matrix(mixed.erbb2.brca.cl.3.log2cpm.ges), 
                                                    ref = as.matrix(brca.normal.log2cpm.ges), per = 1000, method = "ttest", 
                                                    verbose = TRUE)

vipermat.mixed.erbb2.brca.cl.3.gtex <- viper(vipsig.mixed.erbb2.brca.cl.3.gtex, pregul.brca,minsize = 25)
saveRDS(vipermat.mixed.erbb2.brca.cl.3.gtex,file="vipermat.mixed.erbb2.brca.cl.3.gtex.rds")
plot(density(vipermat.mixed.erbb2.brca.cl.3.gtex))

#cluster 4
dim(erbb2.brca.wt.cl.4)
dim(erbb2.brca.tum.cl.4)
breast_gtex<-read.table(file="~/Documents/viper-docs/benchmarking_1/signature_with_gtex/batch_corrected_data/data/batch-correct/gtex/breast-gtex.txt", header=TRUE)
brca.normal<-breast_gtex

identical(rownames(erbb2.brca.tum.cl.4),rownames(erbb2.brca.wt.cl.4))
table(rownames(erbb2.brca.tum.cl.4)%in%rownames(erbb2.brca.wt.cl.4))
common_genes <- intersect(rownames(erbb2.brca.tum.cl.4),rownames(erbb2.brca.wt.cl.4))
length(common_genes)

erbb2.brca.tum.cl.4 <- erbb2.brca.tum.cl.4[match(common_genes,rownames(erbb2.brca.tum.cl.4)),]
erbb2.brca.wt.cl.4 <- erbb2.brca.wt.cl.4[match(common_genes,rownames(erbb2.brca.wt.cl.4)),]
identical(rownames(erbb2.brca.wt.cl.4),rownames(erbb2.brca.tum.cl.4))

identical(rownames(brca.normal),rownames(erbb2.brca.wt.cl.4))
table(rownames(brca.normal)%in%rownames(erbb2.brca.wt.cl.4))
common_genes <- intersect(rownames(brca.normal),rownames(erbb2.brca.wt.cl.4))
length(common_genes)

brca.normal <- brca.normal[match(common_genes,rownames(brca.normal)),]
erbb2.brca.wt.cl.4 <- erbb2.brca.wt.cl.4[match(common_genes,rownames(erbb2.brca.wt.cl.4)),]
identical(rownames(erbb2.brca.wt.cl.4),rownames(brca.normal))

mixed.erbb2.brca.cl.4 <- cbind(brca.normal, erbb2.brca.wt.cl.4,erbb2.brca.tum.cl.4)
dim(mixed.erbb2.brca.cl.4)

#log2cpm+1 normalization
mixed.erbb2.brca.cl.4.log2cpm <- apply(mixed.erbb2.brca.cl.4,2,function(x){
  y <- 1E6*x/sum(x) + 1
  z <- log(y,2)
  return(z)
})

mixed.erbb2.brca.cl.4.log2cpm.ges <- t(apply(mixed.erbb2.brca.cl.4.log2cpm,1,function(x){
  y <- (x - mean(x))/sd(x)
  return(y)
}))
brca.normal.log2cpm <- apply(brca.normal,2,function(x){
  y <- 1E6*x/sum(x) + 1
  z <- log(y,2)
  return(z)
})
brca.normal.log2cpm.ges <- t(apply(brca.normal.log2cpm,1,function(x){
  y <- (x - mean(x))/sd(x)
  return(y)
}))
# for(i in 1:ncol(mixed.brca.cl.2)){
#   mixed.brca.cl.2[,i] <- 1E6*mixed.brca.cl.2[,i]/sum(mixed.brca.cl.2[,i])
# }
# 
# for(i in 1:ncol(brca.normal)){
#   brca.normal[,i] <- 1E6*brca.normal[,i]/sum(brca.normal[,i])
# }

# mixed.erbb2.brca.cl.1.log2cpm.ges <- mixed.erbb2.brca.cl.1.log2cpm.ges[is.finite(rowSums(mixed.erbb2.brca.cl.1.log2cpm.ges)),]
# brca.normal.log2cpm.ges <- brca.normal.log2cpm.ges[is.finite(rowSums(brca.normal.log2cpm.ges)),]

#viper using gtex
vipsig.mixed.erbb2.brca.cl.4.gtex <- viperSignature(eset = as.matrix(mixed.erbb2.brca.cl.4.log2cpm.ges), 
                                                    ref = as.matrix(brca.normal.log2cpm.ges), per = 1000, method = "ttest", 
                                                    verbose = TRUE)

vipermat.mixed.erbb2.brca.cl.4.gtex <- viper(vipsig.mixed.erbb2.brca.cl.4.gtex, pregul.brca,minsize = 25)
saveRDS(vipermat.mixed.erbb2.brca.cl.4.gtex,file="vipermat.mixed.erbb2.brca.cl.4.gtex.rds")
plot(density(vipermat.mixed.erbb2.brca.cl.4.gtex))

#GMM on PIK3CA data
erbb2.gmm.cl.1<-vipermat.mixed.erbb2.brca.cl.1.gtex[match("2064",rownames(vipermat.mixed.erbb2.brca.cl.1.gtex)),]
saveRDS(erbb2.gmm.cl.1,file="erbb2.gmm.cl.1.rds")
erbb2.gmm.cl.2<-vipermat.mixed.erbb2.brca.cl.2.gtex[match("2064",rownames(vipermat.mixed.erbb2.brca.cl.2.gtex)),]
saveRDS(erbb2.gmm.cl.2,file="erbb2.gmm.cl.2.rds")
erbb2.gmm.cl.3<-vipermat.mixed.erbb2.brca.cl.3.gtex[match("2064",rownames(vipermat.mixed.erbb2.brca.cl.3.gtex)),]
saveRDS(erbb2.gmm.cl.3,file="erbb2.gmm.cl.3.rds")
erbb2.gmm.cl.4<-vipermat.mixed.erbb2.brca.cl.4.gtex[match("2064",rownames(vipermat.mixed.erbb2.brca.cl.4.gtex)),]
saveRDS(erbb2.gmm.cl.4,file="erbb2.gmm.cl.4.rds")

#compute the protein activity of all mutant samples vs. the WT
#with the protein activity most similar to GTEx
#cluster 1
dim(erbb2.brca.wt.cl.1)
dim(erbb2.brca.tum.cl.1)

identical(rownames(erbb2.brca.tum.cl.1),rownames(erbb2.brca.wt.cl.1))
table(rownames(erbb2.brca.tum.cl.1)%in%rownames(erbb2.brca.wt.cl.1))
common_genes <- intersect(rownames(erbb2.brca.tum.cl.1),rownames(erbb2.brca.wt.cl.1))
length(common_genes)

erbb2.brca.tum.cl.1 <- erbb2.brca.tum.cl.1[match(common_genes,rownames(erbb2.brca.tum.cl.1)),]
erbb2.brca.wt.cl.1 <- erbb2.brca.wt.cl.1[match(common_genes,rownames(erbb2.brca.wt.cl.1)),]
identical(rownames(erbb2.brca.wt.cl.1),rownames(erbb2.brca.tum.cl.1))

erbb2.brca.wt.mut.cl.1 <- cbind(erbb2.brca.wt.cl.1,erbb2.brca.tum.cl.1)
dim(erbb2.brca.wt.mut.cl.1)

#log2cpm+1 normalization
erbb2.brca.wt.mut.cl.1.log2cpm <- apply(erbb2.brca.wt.mut.cl.1,2,function(x){
  y <- 1E6*x/sum(x) + 1
  z <- log(y,2)
  return(z)
})

erbb2.brca.wt.mut.cl.1.log2cpm.ges <- t(apply(erbb2.brca.wt.mut.cl.1.log2cpm,1,function(x){
  y <- (x - mean(x))/sd(x)
  return(y)
}))
erbb2.brca.wt.cl.1.log2cpm <- apply(erbb2.brca.wt.cl.1,2,function(x){
  y <- 1E6*x/sum(x) + 1
  z <- log(y,2)
  return(z)
})

erbb2.brca.wt.cl.1.log2cpm.ges <- t(apply(erbb2.brca.wt.cl.1.log2cpm,1,function(x){
  y <- (x - mean(x))/sd(x)
  return(y)
}))
# for(i in 1:ncol(wt.mut.brca.cl.2)){
#   wt.mut.brca.cl.2[,i] <- 1E6*wt.mut.brca.cl.2[,i]/sum(wt.mut.brca.cl.2[,i])
# }
# 
# for(i in 1:ncol(brca.wt.cl.2)){
#   brca.wt.cl.2[,i] <- 1E6*brca.wt.cl.2[,i]/sum(brca.wt.cl.2[,i])
# }
# wt.mut.brca.cl.2<-wt.mut.brca.cl.2.log2cpm.ges
# brca.wt.cl.2<-brca.wt.cl.2.log2cpm.ges
# wt.mut.brca.cl.2 <- wt.mut.brca.cl.2[is.finite(rowSums(wt.mut.brca.cl.2)),]
# brca.wt.cl.2 <- brca.wt.cl.2[is.finite(rowSums(brca.wt.cl.2)),]

#viper using wt
vipsig.erbb2.brca.wt.tum.cl.1 <- viperSignature(eset = as.matrix(erbb2.brca.wt.mut.cl.1.log2cpm.ges), 
                                                ref = as.matrix(erbb2.brca.wt.cl.1.log2cpm.ges), per = 1000, method = "ttest", 
                                                verbose = TRUE)

vipermat.erbb2.brca.wt.tum.cl.1 <- viper(vipsig.erbb2.brca.wt.tum.cl.1, pregul.brca,minsize = 25)
saveRDS(vipermat.erbb2.brca.wt.tum.cl.1,file="vipermat.erbb2.brca.wt.tum.cl.1.rds")
plot(density(vipermat.erbb2.brca.wt.tum.cl.1))

#cluster 2
dim(erbb2.brca.wt.cl.2)
dim(erbb2.brca.tum.cl.2)

identical(rownames(erbb2.brca.tum.cl.2),rownames(erbb2.brca.wt.cl.2))
table(rownames(erbb2.brca.tum.cl.2)%in%rownames(erbb2.brca.wt.cl.2))
common_genes <- intersect(rownames(erbb2.brca.tum.cl.2),rownames(erbb2.brca.wt.cl.2))
length(common_genes)

erbb2.brca.tum.cl.2 <- erbb2.brca.tum.cl.2[match(common_genes,rownames(erbb2.brca.tum.cl.2)),]
erbb2.brca.wt.cl.2 <- erbb2.brca.wt.cl.2[match(common_genes,rownames(erbb2.brca.wt.cl.2)),]
identical(rownames(erbb2.brca.wt.cl.2),rownames(erbb2.brca.tum.cl.2))

erbb2.brca.wt.mut.cl.2 <- cbind(erbb2.brca.wt.cl.2,erbb2.brca.tum.cl.2)
dim(erbb2.brca.wt.mut.cl.2)

#log2cpm+1 normalization
erbb2.brca.wt.mut.cl.2.log2cpm <- apply(erbb2.brca.wt.mut.cl.2,2,function(x){
  y <- 1E6*x/sum(x) + 1
  z <- log(y,2)
  return(z)
})

erbb2.brca.wt.mut.cl.2.log2cpm.ges <- t(apply(erbb2.brca.wt.mut.cl.2.log2cpm,1,function(x){
  y <- (x - mean(x))/sd(x)
  return(y)
}))
erbb2.brca.wt.cl.2.log2cpm <- apply(erbb2.brca.wt.cl.2,2,function(x){
  y <- 1E6*x/sum(x) + 1
  z <- log(y,2)
  return(z)
})

erbb2.brca.wt.cl.2.log2cpm.ges <- t(apply(erbb2.brca.wt.cl.2.log2cpm,1,function(x){
  y <- (x - mean(x))/sd(x)
  return(y)
}))
# for(i in 1:ncol(wt.mut.brca.cl.2)){
#   wt.mut.brca.cl.2[,i] <- 1E6*wt.mut.brca.cl.2[,i]/sum(wt.mut.brca.cl.2[,i])
# }
# 
# for(i in 1:ncol(brca.wt.cl.2)){
#   brca.wt.cl.2[,i] <- 1E6*brca.wt.cl.2[,i]/sum(brca.wt.cl.2[,i])
# }
# wt.mut.brca.cl.2<-wt.mut.brca.cl.2.log2cpm.ges
# brca.wt.cl.2<-brca.wt.cl.2.log2cpm.ges
# wt.mut.brca.cl.2 <- wt.mut.brca.cl.2[is.finite(rowSums(wt.mut.brca.cl.2)),]
# brca.wt.cl.2 <- brca.wt.cl.2[is.finite(rowSums(brca.wt.cl.2)),]

#viper using wt
vipsig.erbb2.brca.wt.tum.cl.2 <- viperSignature(eset = as.matrix(erbb2.brca.wt.mut.cl.2.log2cpm.ges), 
                                                ref = as.matrix(erbb2.brca.wt.cl.2.log2cpm.ges), per = 1000, method = "ttest", 
                                                verbose = TRUE)

vipermat.erbb2.brca.wt.tum.cl.2 <- viper(vipsig.erbb2.brca.wt.tum.cl.2, pregul.brca,minsize = 25)
saveRDS(vipermat.erbb2.brca.wt.tum.cl.2,file="vipermat.erbb2.brca.wt.tum.cl.2.rds")
plot(density(vipermat.erbb2.brca.wt.tum.cl.2))


#cluster 3
dim(erbb2.brca.wt.cl.3)
dim(erbb2.brca.tum.cl.3)

identical(rownames(erbb2.brca.tum.cl.3),rownames(erbb2.brca.wt.cl.3))
table(rownames(erbb2.brca.tum.cl.3)%in%rownames(erbb2.brca.wt.cl.3))
common_genes <- intersect(rownames(erbb2.brca.tum.cl.3),rownames(erbb2.brca.wt.cl.3))
length(common_genes)

erbb2.brca.tum.cl.3 <- erbb2.brca.tum.cl.3[match(common_genes,rownames(erbb2.brca.tum.cl.3)),]
erbb2.brca.wt.cl.3 <- erbb2.brca.wt.cl.3[match(common_genes,rownames(erbb2.brca.wt.cl.3)),]
identical(rownames(erbb2.brca.wt.cl.3),rownames(erbb2.brca.tum.cl.3))

erbb2.brca.wt.mut.cl.3 <- cbind(erbb2.brca.wt.cl.3,erbb2.brca.tum.cl.3)
dim(erbb2.brca.wt.mut.cl.3)

#log2cpm+1 normalization
erbb2.brca.wt.mut.cl.3.log2cpm <- apply(erbb2.brca.wt.mut.cl.3,2,function(x){
  y <- 1E6*x/sum(x) + 1
  z <- log(y,2)
  return(z)
})

erbb2.brca.wt.mut.cl.3.log2cpm.ges <- t(apply(erbb2.brca.wt.mut.cl.3.log2cpm,1,function(x){
  y <- (x - mean(x))/sd(x)
  return(y)
}))
erbb2.brca.wt.cl.3.log2cpm <- apply(erbb2.brca.wt.cl.3,2,function(x){
  y <- 1E6*x/sum(x) + 1
  z <- log(y,2)
  return(z)
})

erbb2.brca.wt.cl.3.log2cpm.ges <- t(apply(erbb2.brca.wt.cl.3.log2cpm,1,function(x){
  y <- (x - mean(x))/sd(x)
  return(y)
}))
# for(i in 1:ncol(wt.mut.brca.cl.2)){
#   wt.mut.brca.cl.2[,i] <- 1E6*wt.mut.brca.cl.2[,i]/sum(wt.mut.brca.cl.2[,i])
# }
# 
# for(i in 1:ncol(brca.wt.cl.2)){
#   brca.wt.cl.2[,i] <- 1E6*brca.wt.cl.2[,i]/sum(brca.wt.cl.2[,i])
# }
# wt.mut.brca.cl.2<-wt.mut.brca.cl.2.log2cpm.ges
# brca.wt.cl.2<-brca.wt.cl.2.log2cpm.ges
# wt.mut.brca.cl.2 <- wt.mut.brca.cl.2[is.finite(rowSums(wt.mut.brca.cl.2)),]
# brca.wt.cl.2 <- brca.wt.cl.2[is.finite(rowSums(brca.wt.cl.2)),]

#viper using wt
vipsig.erbb2.brca.wt.tum.cl.3 <- viperSignature(eset = as.matrix(erbb2.brca.wt.mut.cl.3.log2cpm.ges), 
                                                ref = as.matrix(erbb2.brca.wt.cl.3.log2cpm.ges), per = 1000, method = "ttest", 
                                                verbose = TRUE)

vipermat.erbb2.brca.wt.tum.cl.3 <- viper(vipsig.erbb2.brca.wt.tum.cl.3, pregul.brca,minsize = 25)
saveRDS(vipermat.erbb2.brca.wt.tum.cl.3,file="vipermat.erbb2.brca.wt.tum.cl.3.rds")
plot(density(vipermat.erbb2.brca.wt.tum.cl.3))

#cluster 4
dim(erbb2.brca.wt.cl.4)
dim(erbb2.brca.tum.cl.4)

identical(rownames(erbb2.brca.tum.cl.4),rownames(erbb2.brca.wt.cl.4))
table(rownames(erbb2.brca.tum.cl.4)%in%rownames(erbb2.brca.wt.cl.4))
common_genes <- intersect(rownames(erbb2.brca.tum.cl.4),rownames(erbb2.brca.wt.cl.4))
length(common_genes)

erbb2.brca.tum.cl.4 <- erbb2.brca.tum.cl.4[match(common_genes,rownames(erbb2.brca.tum.cl.4)),]
erbb2.brca.wt.cl.4 <- erbb2.brca.wt.cl.4[match(common_genes,rownames(erbb2.brca.wt.cl.4)),]
identical(rownames(erbb2.brca.wt.cl.4),rownames(erbb2.brca.tum.cl.4))

erbb2.brca.wt.mut.cl.4 <- cbind(erbb2.brca.wt.cl.4,erbb2.brca.tum.cl.4)
dim(erbb2.brca.wt.mut.cl.4)

#log2cpm+1 normalization
erbb2.brca.wt.mut.cl.4.log2cpm <- apply(erbb2.brca.wt.mut.cl.4,2,function(x){
  y <- 1E6*x/sum(x) + 1
  z <- log(y,2)
  return(z)
})

erbb2.brca.wt.mut.cl.4.log2cpm.ges <- t(apply(erbb2.brca.wt.mut.cl.4.log2cpm,1,function(x){
  y <- (x - mean(x))/sd(x)
  return(y)
}))
erbb2.brca.wt.cl.4.log2cpm <- apply(erbb2.brca.wt.cl.4,2,function(x){
  y <- 1E6*x/sum(x) + 1
  z <- log(y,2)
  return(z)
})

erbb2.brca.wt.cl.4.log2cpm.ges <- t(apply(erbb2.brca.wt.cl.4.log2cpm,1,function(x){
  y <- (x - mean(x))/sd(x)
  return(y)
}))
# for(i in 1:ncol(wt.mut.brca.cl.2)){
#   wt.mut.brca.cl.2[,i] <- 1E6*wt.mut.brca.cl.2[,i]/sum(wt.mut.brca.cl.2[,i])
# }
# 
# for(i in 1:ncol(brca.wt.cl.2)){
#   brca.wt.cl.2[,i] <- 1E6*brca.wt.cl.2[,i]/sum(brca.wt.cl.2[,i])
# }
# wt.mut.brca.cl.2<-wt.mut.brca.cl.2.log2cpm.ges
# brca.wt.cl.2<-brca.wt.cl.2.log2cpm.ges
# wt.mut.brca.cl.2 <- wt.mut.brca.cl.2[is.finite(rowSums(wt.mut.brca.cl.2)),]
# brca.wt.cl.2 <- brca.wt.cl.2[is.finite(rowSums(brca.wt.cl.2)),]

#viper using wt
vipsig.erbb2.brca.wt.tum.cl.4 <- viperSignature(eset = as.matrix(erbb2.brca.wt.mut.cl.4.log2cpm.ges), 
                                                ref = as.matrix(erbb2.brca.wt.cl.4.log2cpm.ges), per = 1000, method = "ttest", 
                                                verbose = TRUE)

vipermat.erbb2.brca.wt.tum.cl.4 <- viper(vipsig.erbb2.brca.wt.tum.cl.4, pregul.brca,minsize = 25)
saveRDS(vipermat.erbb2.brca.wt.tum.cl.4,file="vipermat.erbb2.brca.wt.tum.cl.4.rds")
plot(density(vipermat.erbb2.brca.wt.tum.cl.4))

#extracting erbb2.mut samples from all the 4 clusters
# erbb2.brca.Tumor_Sample_Barcode<-read.csv("~/Documents/Mutation_mapper/erbb2_mutations_brca_csv.csv",header=TRUE, sep=",")
# rownames(erbb2.brca.Tumor_Sample_Barcode) <- erbb2.brca.Tumor_Sample_Barcode[,1]
# erbb2.brca.Tumor_Sample_Barcode<-erbb2.brca.Tumor_Sample_Barcode[,-1]
# rownames(erbb2.brca.Tumor_Sample_Barcode)
# 
# rownames(erbb2.brca.Tumor_Sample_Barcode)<-gsub("\\-","\\.",rownames(erbb2.brca.Tumor_Sample_Barcode))
vipermat.erbb2.brca.wt.tum.cl.1.1<-vipermat.erbb2.brca.wt.tum.cl.1
# colnames(vipermat.brca.cl.1.wt.1)<-gsub('.{16}$','',colnames(vipermat.brca.cl.1.wt.1))
vipermat.erbb2.brca.tum.cl.1.1<- 
  vipermat.erbb2.brca.wt.tum.cl.1.1[,intersect(rownames(erbb2.brca.Tumor_Sample_Barcode),colnames(vipermat.erbb2.brca.wt.tum.cl.1.1))]
dim(vipermat.erbb2.brca.tum.cl.1.1)
vipermat.erbb2.brca.wt.cl.1.1<- 
  vipermat.erbb2.brca.wt.tum.cl.1.1[,outersect(colnames(vipermat.erbb2.brca.wt.tum.cl.1.1),colnames(vipermat.erbb2.brca.tum.cl.1.1))]
dim(vipermat.erbb2.brca.wt.cl.1.1)

vipermat.erbb2.brca.wt.tum.cl.2.1<-vipermat.erbb2.brca.wt.tum.cl.2
# colnames(vipermat.brca.cl.1.wt.1)<-gsub('.{16}$','',colnames(vipermat.brca.cl.1.wt.1))
vipermat.erbb2.brca.tum.cl.2.1<- 
  vipermat.erbb2.brca.wt.tum.cl.2.1[,intersect(rownames(erbb2.brca.Tumor_Sample_Barcode),colnames(vipermat.erbb2.brca.wt.tum.cl.2.1))]
dim(vipermat.erbb2.brca.tum.cl.2.1)
vipermat.erbb2.brca.wt.cl.2.1<- 
  vipermat.erbb2.brca.wt.tum.cl.2.1[,outersect(colnames(vipermat.erbb2.brca.wt.tum.cl.2.1),colnames(vipermat.erbb2.brca.tum.cl.2.1))]
dim(vipermat.erbb2.brca.wt.cl.2.1)

vipermat.erbb2.brca.wt.tum.cl.3.1<-vipermat.erbb2.brca.wt.tum.cl.3
# colnames(vipermat.brca.cl.1.wt.1)<-gsub('.{16}$','',colnames(vipermat.brca.cl.1.wt.1))
vipermat.erbb2.brca.tum.cl.3.1<- 
  vipermat.erbb2.brca.wt.tum.cl.3.1[,intersect(rownames(erbb2.brca.Tumor_Sample_Barcode),colnames(vipermat.erbb2.brca.wt.tum.cl.3.1))]
dim(vipermat.erbb2.brca.tum.cl.3.1)
vipermat.erbb2.brca.wt.cl.3.1<- 
  vipermat.erbb2.brca.wt.tum.cl.3.1[,outersect(colnames(vipermat.erbb2.brca.wt.tum.cl.3.1),colnames(vipermat.erbb2.brca.tum.cl.3.1))]
dim(vipermat.erbb2.brca.wt.cl.3.1)

vipermat.erbb2.brca.wt.tum.cl.4.1<-vipermat.erbb2.brca.wt.tum.cl.4
# colnames(vipermat.brca.cl.1.wt.1)<-gsub('.{16}$','',colnames(vipermat.brca.cl.1.wt.1))
vipermat.erbb2.brca.tum.cl.4.1<- 
  vipermat.erbb2.brca.wt.tum.cl.4.1[,intersect(rownames(erbb2.brca.Tumor_Sample_Barcode),colnames(vipermat.erbb2.brca.wt.tum.cl.4.1))]
dim(vipermat.erbb2.brca.tum.cl.4.1)
vipermat.erbb2.brca.wt.cl.4.1<- 
  vipermat.erbb2.brca.wt.tum.cl.4.1[,outersect(colnames(vipermat.erbb2.brca.wt.tum.cl.4.1),colnames(vipermat.erbb2.brca.tum.cl.4.1))]
dim(vipermat.erbb2.brca.wt.cl.4.1)

#combine all erbb2.mut samples into a single matrix
mix.erbb2.brca.cl.1.2.3.4<-cbind(vipermat.erbb2.brca.tum.cl.1.1,
                                 vipermat.erbb2.brca.tum.cl.2.1,
                                 vipermat.erbb2.brca.tum.cl.3.1,
                                 vipermat.erbb2.brca.tum.cl.4.1)
dim(mix.erbb2.brca.cl.1.2.3.4)
saveRDS(mix.erbb2.brca.cl.1.2.3.4,file="mix.erbb2.brca.cl.1.2.3.4.rds")
mix.erbb2.brca.cl.1.2.3.4<-readRDS(file="mix.erbb2.brca.cl.1.2.3.4.rds")

mix.erbb2.brca.cl.1.2.3.4.wt<-cbind(vipermat.erbb2.brca.wt.cl.1.1,
                                    vipermat.erbb2.brca.wt.cl.2.1,
                                    vipermat.erbb2.brca.wt.cl.3.1,
                                    vipermat.erbb2.brca.wt.cl.4.1)
dim(mix.erbb2.brca.cl.1.2.3.4.wt)
saveRDS(mix.erbb2.brca.cl.1.2.3.4.wt,file="mix.erbb2.brca.cl.1.2.3.4.wt.rds")
mix.erbb2.brca.cl.1.2.3.4.wt<-readRDS(file="mix.erbb2.brca.cl.1.2.3.4.wt.rds")

# mix.pik3ca.brca.cl.1.2<-cbind(pik3ca.brca.cl.1,pik3ca.brca.cl.2)
# dim(mix.pik3ca.brca.cl.1.2)
# saveRDS(mix.pik3ca.brca.cl.1.2,file="mix.pik3ca.brca.cl.1.2.rds")

#downstream analysis
#mutations
mutation_file_erbb2_brca<-read.csv("~/Documents/Mutation_mapper/brca_mutations/erbb2_brca_mutations.csv",header=TRUE, sep=",")
rownames(mutation_file_erbb2_brca) <- mutation_file_erbb2_brca[,1]
mutation_file_erbb2_brca<-mutation_file_erbb2_brca[,-1]
rownames(mutation_file_erbb2_brca)

rownames(mutation_file_erbb2_brca)
colnames(mutation_file_erbb2_brca)

rownames(mutation_file_erbb2_brca)<-gsub("\\-","\\.",rownames(mutation_file_erbb2_brca))
rownames(mutation_file_erbb2_brca)<-gsub('.{3}$','',rownames(mutation_file_erbb2_brca))
# annot.col.1 <- annot.col[,match(rownames(mutation_file),
#                                 colnames(vipermat.egfr))] 
#cbind(annot.col.1,)
vipermat.erbb2.brca<-mix.erbb2.brca.cl.1.2.3.4
common.rows <- intersect(colnames(vipermat.erbb2.brca),rownames(mutation_file_erbb2_brca))
length(common.rows)
# annot.col.1<-annot.col
mutation.file.1<-mutation_file_erbb2_brca
# #mutation.file.1<-mutation_file[,c(7:8)]
# #mutation.file.1<-mutation_file[,c(7:10,18,27:28,30,31)]
# mutation.file.1<-mutation_file[,c(7:8,18)]
# annot.col.1 <- annot.col.1[match(common.rows,rownames(annot.col.1)),]
mutation.file.1 <- mutation.file.1[match(common.rows,rownames(mutation.file.1)),]
#mutation.file.1 <- common.rows%in%mutation.file.1
dim(mutation.file.1)

annot.col.1.mutation.file.1<-mutation.file.1
dim(annot.col.1.mutation.file.1)
head(annot.col.1.mutation.file.1)
saveRDS(annot.col.1.mutation.file.1,file="annot.col.1.mutation.file.1.erbb2.brca.gof.new.rds")
write.csv(annot.col.1.mutation.file.1,file="annot.col.1.mutation.file.1.erbb2.brca.gof.new.csv")


test.1<-colnames(annot.col.1.mutation.file.1)
#test.1<-c('R273H')
for(i in 1:length(test.1)) { 
  
  nam <- noquote(paste0(test.1[i], '.erbb2.brca.samples'))
  # test.1[i]<-noquote(test.1[i])
  #t[i]<-noquote(test.1[i])
  nam.1 <- noquote(paste0(test.1[i]))
  #get(nam.1)
  assign(nam, (rownames(annot.col.1.mutation.file.1[which(annot.col.1.mutation.file.1[,nam.1]=='1'),])))
  # assign(nam, (brca.tcga.tumor[,intersect(colnames(brca.tcga.tumor),get(nam))]))
  
}
A159P.erbb2.brca.samples
R273C.erbb2.brca.samples

#E746_A750del.samples<-rownames(annot.col.1.mutation.file.1[which(annot.col.1.mutation.file.1$E746_A750del=='1'),])

vipermat.plot.data.Int.gtex.c1.samples.2.copy<-vipermat.erbb2.brca
test.1<-colnames(annot.col.1.mutation.file.1)
#test.1<-c('R273H')
for(i in 1:length(test.1)) { 
  
  nam <- noquote(paste0(test.1[i], '.erbb2.brca.samples.erbb2'))
  # test.1[i]<-noquote(test.1[i])
  #t[i]<-noquote(test.1[i])
  nam.1 <- noquote(paste0(test.1[i],'.erbb2.brca.samples'))
  #get(nam.1)
  # assign(nam, (rownames(annot.col.1.mutation.file.1[which(annot.col.1.mutation.file.1[,nam.1]=='1'),])))
  assign(nam, (vipermat.plot.data.Int.gtex.c1.samples.2.copy[,intersect(get(nam.1),colnames(vipermat.plot.data.Int.gtex.c1.samples.2.copy))]))
  
}


#domain-wise mutations
#domain-wise mutations
#domain-wise mutations
erbb2_brca_all_labels<-cbind('S310F','A1160V','A355Qfs.76',
                             'ADGRA3.ERBB2','AP2B1.ERBB2','CACNB1.ERBB2',
                             'CORO1B.ERBB2','D769H','E405D',
                             'E939G','EME1.ERBB2','ERBB2.ABI3',
                             'ERBB2.ARL5C','ERBB2.MED1','ERBB2.PPP1R1B'
                             ,'ERBB2.PSMB3','ERBB2.SDF4','ERBB2.SLC29A3',
                             'ERBB2.SRCIN1','ERBB2.ZAN','G309A',
                             'G621Afs.31','H470Q','L755M','L755W',
                             'MYO18A.ERBB2','PIP4K2B.ERBB2','R340G',
                             'R678Q','S305C','T306M','V797A','V842I'
                             ,'W482Gfs.74','Y1127_A1129del',
                             'L755S','V777L','D769Y'
)
Recep_L.52.172.ERBB2.brca.domain.samples<-cbind()
Furin_like.190.343.ERBB2.brca.domain.samples<-cbind('S310F','G309A','R340G','S305C',
                                                    'T306M')
Recep_L.366.485.ERBB2.brca.domain.samples<-cbind('E405D','H470Q','W482Gfs.74')
GF_recep_IV.511.642.ERBB2.brca.domain.samples<-cbind()
Pkinase_Tyr.721.975.ERBB2.brca.domain.samples<-cbind('D769H','L755M','L755W','V797A',
                                                     'V842I','L755S','V777L','D769Y',
                                                     'E939G')

cbind_all<-cbind(Recep_L.52.172.ERBB2.brca.domain.samples,
                 Furin_like.190.343.ERBB2.brca.domain.samples,
                 Recep_L.366.485.ERBB2.brca.domain.samples,
                 GF_recep_IV.511.642.ERBB2.brca.domain.samples,
                 Pkinase_Tyr.721.975.ERBB2.brca.domain.samples)
som3<-outersect(erbb2_brca_all_labels,cbind_all)
som3
inter.domain.samples<-cbind( "A1160V",         "A355Qfs.76",     "ADGRA3.ERBB2",   "AP2B1.ERBB2",    "CACNB1.ERBB2",  
                             "CORO1B.ERBB2",   "EME1.ERBB2",     "ERBB2.ABI3",     "ERBB2.ARL5C",    "ERBB2.MED1",    
                             "ERBB2.PPP1R1B",  "ERBB2.PSMB3",    "ERBB2.SDF4",     "ERBB2.SLC29A3",  "ERBB2.SRCIN1",  
                            "ERBB2.ZAN",      "G621Afs.31",     "MYO18A.ERBB2",   "PIP4K2B.ERBB2",  "R678Q",         
                             "Y1127_A1129del")
########################
Furin_like.190.343.ERBB2.brca.domain.samples<-cbind(S310F.erbb2.brca.samples.erbb2,G309A.erbb2.brca.samples.erbb2,R340G.erbb2.brca.samples.erbb2,S305C.erbb2.brca.samples.erbb2,
                                                    T306M.erbb2.brca.samples.erbb2)
Recep_L.366.485.ERBB2.brca.domain.samples<-cbind(E405D.erbb2.brca.samples.erbb2,
                                                 H470Q.erbb2.brca.samples.erbb2,
                                                 W482Gfs.74.erbb2.brca.samples.erbb2)
Pkinase_Tyr.721.975.ERBB2.brca.domain.samples<-cbind(D769H.erbb2.brca.samples.erbb2,L755M.erbb2.brca.samples.erbb2,L755W.erbb2.brca.samples.erbb2,V797A.erbb2.brca.samples.erbb2,
                                                     V842I.erbb2.brca.samples.erbb2,L755S.erbb2.brca.samples.erbb2,V777L.erbb2.brca.samples.erbb2,D769Y.erbb2.brca.samples.erbb2,
                                                     E939G.erbb2.brca.samples.erbb2)
inter.domain.samples<-cbind( A1160V.erbb2.brca.samples.erbb2,         A355Qfs.76.erbb2.brca.samples.erbb2,     ADGRA3.ERBB2.erbb2.brca.samples.erbb2,   AP2B1.ERBB2.erbb2.brca.samples.erbb2,    CACNB1.ERBB2.erbb2.brca.samples.erbb2,  
                             CORO1B.ERBB2.erbb2.brca.samples.erbb2,   EME1.ERBB2.erbb2.brca.samples.erbb2,     ERBB2.ABI3.erbb2.brca.samples.erbb2,     ERBB2.ARL5C.erbb2.brca.samples.erbb2,    ERBB2.MED1.erbb2.brca.samples.erbb2,    
                             ERBB2.PPP1R1B.erbb2.brca.samples.erbb2,  ERBB2.PSMB3.erbb2.brca.samples.erbb2,    ERBB2.SDF4.erbb2.brca.samples.erbb2,     ERBB2.SLC29A3.erbb2.brca.samples.erbb2,  ERBB2.SRCIN1.erbb2.brca.samples.erbb2,  
                             ERBB2.ZAN.erbb2.brca.samples.erbb2,      G621Afs.31.erbb2.brca.samples.erbb2,     MYO18A.ERBB2.erbb2.brca.samples.erbb2,   PIP4K2B.ERBB2.erbb2.brca.samples.erbb2,  R678Q.erbb2.brca.samples.erbb2,         
                             Y1127_A1129del.erbb2.brca.samples.erbb2)

current.stouffer.Furin_like.190.343.ERBB2.brca.domain.samples<- apply(Furin_like.190.343.ERBB2.brca.domain.samples,1,function(x){
  y <- sum(x)/sqrt(length(x))
  #y <- sum(w*x)/sqrt(sum(w^2))
  return(y)
})
current.stouffer.Recep_L.366.485.ERBB2.brca.domain.samples<- apply(Recep_L.366.485.ERBB2.brca.domain.samples,1,function(x){
  y <- sum(x)/sqrt(length(x))
  #y <- sum(w*x)/sqrt(sum(w^2))
  return(y)
})
current.stouffer.Pkinase_Tyr.721.975.ERBB2.brca.domain.samples <- apply(Pkinase_Tyr.721.975.ERBB2.brca.domain.samples,1,function(x){
  y <- sum(x)/sqrt(length(x))
  #y <- sum(w*x)/sqrt(sum(w^2))
  return(y)
})
current.stouffer.inter.domain.samples <- apply(inter.domain.samples,1,function(x){
  y <- sum(x)/sqrt(length(x))
  #y <- sum(w*x)/sqrt(sum(w^2))
  return(y)
})

current.stouffer.domains.erbb2.brca<-
  cbind(current.stouffer.Furin_like.190.343.ERBB2.brca.domain.samples,
        current.stouffer.Recep_L.366.485.ERBB2.brca.domain.samples,
        current.stouffer.Pkinase_Tyr.721.975.ERBB2.brca.domain.samples,
        current.stouffer.inter.domain.samples)

colnames(current.stouffer.domains.erbb2.brca)<- c('Furin_like.190.343.ERBB2.brca',
                                                  'Recep_L.366.485.ERBB2.brca',
                                                  'Pkinase_Tyr.721.975.ERBB2',
                                                  'Inter.domain.samples')

# colnames.current.stouffer.domains.egfr<- c("ECD.I \n (L62R \n I91V)",
#                                            "ECD.II \n (L210L \n R222L \n V300M)",
#                                            "ECD.III\n(L387M\nQ432H\nG465R)",
#                                            "TM\n(G652G)",
#                                            "Exon18 \n (E709_T710delinsD \n G719A\nG719C)",
#                                            "Exon19\n(E746_A750del\nL747_E749del\nL747_A750del\nL747_T751del\nT751P\nT751Nfs.15\nS752Pfs.3\nP753P\nK754_I759del\nK754I",
#                                            "Exon20\n(A767_V769du\nS768I\nP772_H773du\nH773du\nT790M)",
#                                            "Exon21 \n (L833V \n I853I \n L858R \n L861Q \n E866K \n G901V \n S921R)",
#                                            "CTT\n(R1052I\nD1083Efs.11)")

common.tf.cotf.1 <- intersect(rownames(current.stouffer.domains.erbb2.brca),rownames(tf.cotf))
length(common.tf.cotf.1)

all.domains.mutations.erbb2.brca <- current.stouffer.domains.erbb2.brca[match(common.tf.cotf.1,rownames(current.stouffer.domains.erbb2.brca)),]
dim(all.domains.mutations.erbb2.brca)

#rownames(all.domains.mutations.erbb2.brca)<-getSYMBOL(rownames(all.domains.mutations.erbb2.brca),data=org.Hs.eg)
zscore_data.1 <- t(apply(all.domains.mutations.erbb2.brca,1,function(x){
  y <- (x - mean(x))/sd(x)
  return(y)
}))

rank.zscore_data.1 <- apply(zscore_data.1,2,rank)

keep.num.1 <- 2352
keep.index.1 <- apply(rank.zscore_data.1,1,function(x){
  y <- ((max(x) > (nrow(rank.zscore_data.1) + 1 - keep.num.1)) | min(x) <= keep.num.1)
  return(y)
})

temp.plot.data.2.1 <- all.domains.mutations.erbb2.brca[keep.index.1,]

mat_breaks <- quantile_breaks(temp.plot.data.2.1,n = 100)

test3<-pheatmap(temp.plot.data.2.1, color = colorRampPalette(rev(brewer.pal(n = 7, 
                                                                            name = "RdBu")))
                (length(mat_breaks)), fontsize_row = 5, fontsize_col = 10,
                main = "erbb2.brca.TF.CoTF\nMutations.domain.wise.all",
                show_rownames = FALSE, 
                show_colnames = TRUE, 
                breaks = mat_breaks, 
                #labels_col = as.character(colnames.current.stouffer.domains.egfr),
                cellheight = 0.3,
                #  margins=c(100,100),
                cluster_cols = TRUE, cluster_rows = TRUE,angle_col ="0")
save_pheatmap_pdf(test3,"erbb2.brca.Mutations.domain.wise.all.TF.CoTF.v3.pdf")
saveRDS(temp.plot.data.2.1,file="erbb2.tf.cotf.all.domain.wise.v3.rds")

var.keep.index.1 <- apply(all.domains.mutations.erbb2.brca,1,sd)
var.keep.index.1 <- (rank(var.keep.index.1) > length(var.keep.index.1) - 100)

temp.plot.data.3.1 <- all.domains.mutations.erbb2.brca[var.keep.index.1,]

mat_breaks <- quantile_breaks(temp.plot.data.3.1,n = 100)

rownames(temp.plot.data.3.1)<-getSYMBOL(rownames(temp.plot.data.3.1),data='org.Hs.eg')

test3<-pheatmap(temp.plot.data.3.1, color = colorRampPalette(rev(brewer.pal(n = 7, 
                                                                            name = "RdBu")))
                (length(mat_breaks)), fontsize_row = 5, fontsize_col = 10,
                main = "erbb2.brca.Top.100\nMutations.domain.wise",
                show_rownames = TRUE, show_colnames = TRUE, breaks = mat_breaks, 
                #labels_col = as.character(colnames.current.stouffer.domains.egfr),
                cellheight = 7.7,
                cluster_cols = TRUE, cluster_rows = TRUE,angle_col ="0")
save_pheatmap_pdf(test3,"erbb2.brca.Mutations.domain.wise.top100.v3.pdf")
saveRDS(temp.plot.data.3.1,file="erbb2.tf.cotf.top100.domain.wise.v3.rds")

#mutation comparison - high confidence known gof vs unknown
#P53_TAD_6_29_erbb2_luad<-cbind()

high.conf.gof.erbb2_brca_labels<-cbind(
  'L755S',
  'V777L',
  'D769Y')
high.conf.gof.erbb2_brca<-cbind(L755S.erbb2.brca.samples.erbb2,
                                V777L.erbb2.brca.samples.erbb2,
                                D769Y.erbb2.brca.samples.erbb2)

Pkinase_Tyr.721.975.ERBB2.brca.domain.samples.others<-cbind(L755M.erbb2.brca.samples.erbb2,L755W.erbb2.brca.samples.erbb2,V797A.erbb2.brca.samples.erbb2,
                                                     V842I.erbb2.brca.samples.erbb2,V777L.erbb2.brca.samples.erbb2,
                                                     E939G.erbb2.brca.samples.erbb2)


current.stouffer.high.conf.gof.erbb2_brca <- apply(high.conf.gof.erbb2_brca,1,function(x){
  y <- sum(x)/sqrt(length(x))
  #y <- sum(w*x)/sqrt(sum(w^2))
  return(y)
})

current.stouffer.Pkinase_Tyr.721.975.ERBB2.brca.domain.samples.others <- apply(Pkinase_Tyr.721.975.ERBB2.brca.domain.samples.others,1,function(x){
  y <- sum(x)/sqrt(length(x))
  #y <- sum(w*x)/sqrt(sum(w^2))
  return(y)
})

current.stouffer.high.conf.domains.erbb2.brca<-
  cbind(current.stouffer.high.conf.gof.erbb2_brca,
        current.stouffer.Furin_like.190.343.ERBB2.brca.domain.samples,
        current.stouffer.Recep_L.366.485.ERBB2.brca.domain.samples,
        current.stouffer.Pkinase_Tyr.721.975.ERBB2.brca.domain.samples.others,
        current.stouffer.inter.domain.samples)

colnames(current.stouffer.high.conf.domains.erbb2.brca)<- c("High.Conf.gof \n(L755S,V777L,D769Y)",
                                                            'Furin_like.190.343.ERBB2.brca',
                                                            'Recep_L.366.485.ERBB2.brca',
                                                            'Pkinase_Tyr.721.975.ERBB2.others',
                                                            'Inter.domain.samples')

common.tf.cotf.1 <- intersect(rownames(current.stouffer.high.conf.domains.erbb2.brca),rownames(tf.cotf))
length(common.tf.cotf.1)

all.domains.high.conf.mutations.erbb2.brca <- current.stouffer.high.conf.domains.erbb2.brca[match(common.tf.cotf.1,rownames(current.stouffer.high.conf.domains.erbb2.brca)),]
dim(all.domains.high.conf.mutations.erbb2.brca)

rownames(all.domains.high.conf.mutations.erbb2.brca)<-getSYMBOL(rownames(all.domains.high.conf.mutations.erbb2.brca),data='org.Hs.eg')
zscore_data.1 <- t(apply(all.domains.high.conf.mutations.erbb2.brca,1,function(x){
  y <- (x - mean(x))/sd(x)
  return(y)
}))

rank.zscore_data.1 <- apply(zscore_data.1,2,rank)

keep.num.1 <- 2352
keep.index.1 <- apply(rank.zscore_data.1,1,function(x){
  y <- ((max(x) > (nrow(rank.zscore_data.1) + 1 - keep.num.1)) | min(x) <= keep.num.1)
  return(y)
})

temp.plot.data.2.1 <- all.domains.high.conf.mutations.erbb2.brca[keep.index.1,]

mat_breaks <- quantile_breaks(temp.plot.data.2.1,n = 100)

test3<-pheatmap(temp.plot.data.2.1, color = colorRampPalette(rev(brewer.pal(n = 7, 
                                                                            name = "RdBu")))
                (length(mat_breaks)), fontsize_row = 5, fontsize_col = 8,
                main = "erbb2.brca.TF.CoTF\nHigh.Confidence.gof.Mutations.domain.wise.all\nKnown gof VS rest",
                show_rownames = FALSE, 
                show_colnames = TRUE, 
                breaks = mat_breaks, 
                #labels_col = as.character(colnames.current.stouffer.domains.egfr),
                cellheight = 0.30,
                #  margins=c(100,100),
                cluster_cols = TRUE, cluster_rows = TRUE)
save_pheatmap_pdf(test3,"erbb2.brca.Mutations.domain.wise.all.TF.CoTF.high.conf.v3.pdf")
saveRDS(temp.plot.data.2.1,file="erbb2.brca.Mutations.domain.wise.all.TF.CoTF.high.conf.v3.rds")

var.keep.index.1 <- apply(all.domains.high.conf.mutations.erbb2.brca,1,sd)
var.keep.index.1 <- (rank(var.keep.index.1) > length(var.keep.index.1) - 2352)
#var.keep.index.1 <- (rank(var.keep.index.1) < length(var.keep.index.1) - 100)

temp.plot.data.3.1 <- all.domains.high.conf.mutations.erbb2.brca[var.keep.index.1,]

mat_breaks <- quantile_breaks(temp.plot.data.3.1,n = 100)

#tpp 50 bottom50
common.tf.cotf.1 <- intersect(names(current.stouffer.high.conf.gof.erbb2_brca),rownames(tf.cotf))
cc.new<-current.stouffer.high.conf.gof.erbb2_brca[common.tf.cotf.1]
reference <- c(sort(cc.new,decreasing = TRUE)[1:50],rev(sort(cc.new,decreasing = FALSE)[1:50]))

# order the samples in the protein activity matrix by their enrichment and the rows by the protein activity in the consensus signature
protein.sort <- names(sort(reference,decreasing = TRUE))

# view results with a heatmap
protein.sort<-getSYMBOL(protein.sort,data='org.Hs.eg')

temp.plot.data.3.1 <- temp.plot.data.3.1[match(protein.sort,rownames(temp.plot.data.3.1)),]

mat_breaks <- quantile_breaks(temp.plot.data.3.1,n = 100)

#rownames(temp.plot.data.3.1)<-getSYMBOL(rownames(temp.plot.data.3.1),data='org.Hs.eg')

test3<-pheatmap(temp.plot.data.3.1,kmeans_k = NA,  color = colorRampPalette(rev(brewer.pal(n = 7, 
                                                                                           name = "RdBu")))
                (length(mat_breaks)), fontsize_row = 5, fontsize_col = 8,
                main = "erbb2.brca.Top.100\nHigh.Confidence.gof.Mutations.domain.wise.all\nKnown gof VS rest",
                show_rownames = TRUE, show_colnames = TRUE, breaks = mat_breaks, 
                cluster_cols = TRUE, cluster_rows = FALSE,cellheight = 7)

# test3<-pheatmap(temp.plot.data.3.1, color = colorRampPalette(rev(brewer.pal(n = 7, 
#                                                                             name = "RdBu")))
#                 (length(mat_breaks)), fontsize_row = 4, fontsize_col = 5,
#                 main = "erbb2.UCEC.Top.100\nHigh.Confidence.gof.Mutations.domain.wise.all\nKnown gof VS rest",
#                 show_rownames = TRUE, show_colnames = TRUE, breaks = mat_breaks, 
#                 #labels_col = as.character(colnames.current.stouffer.domains.egfr),
#                 cellheight = 3.1,
#                 cluster_cols = TRUE, cluster_rows = TRUE,angle_col ="0")
save_pheatmap_pdf(test3,"erbb2.brca.Mutations.domain.wise.top100.high.conf.v3.pdf")
saveRDS(temp.plot.data.3.1,file="erbb2.brca.Mutations.domain.wise.top100.high.conf.v3.rds")


