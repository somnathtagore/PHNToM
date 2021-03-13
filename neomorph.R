
## Date = 11/20/19
## Group = Califano Lab


# load packages

# set analysis seed
set.seed(1)


# load in the consensus protein activity signature
#ref.sig <- readRDS("../../brca-gtex/erbb2_brca/consensus_signature.rds")
ref.sig <- readRDS("consensus_signature_erbb2_brca.rds")
#ref.sig <- readRDS("consensus_signature_cluster_1.rds")

ref.sig <- ref.sig[,1]
#ref.sig <- current.mrs.erbb2.brca$es$nes
#mapIds(org.Hs.eg.db,rownames(current.msa.erbb2.brca$es$nes.bt),'ENTREZID','SYMBOL')

# load in the protein activity signatures for individually mutated samples
#test.vpmat <- readRDS("../../brca-gtex/unknown_matrix.rds")
#test.vpmat <- vipermat.unknown.mut.erbb2.brca.tf.cotf[,1:61]

# load in the protein activity signatures for individually mutated samples
test.vpmat <- readRDS("vipermat.known.unknown.mut.erbb2.with.all.high.conf.lof.brca.tf.cotf.rds")
test.vpmat <- readRDS("vipermat.known.unknown.mut.erbb2.with.all.high.conf.gof.brca.tf.cotf.merge.mut.rds")
#test.vpmat <- vipermat.unknown.mut.erbb2.tf.cotf

#test.vpmat<-test.vpmat[1:2313,]

#test.vpmat.1<-test.vpmat[,1:47] #gof
#test.vpmat.2<-test.vpmat[,41:50] #neutral

#neutral.1<-read.csv("../../brca-gtex/neutrals.csv",header=TRUE, sep=",")
#rownames(neutral.1) <- neutral.1[,1]
#neutral.2<-read.csv("neutral.2.csv",header=TRUE, sep=",")
#rownames(neutral.2) <- neutral.2[,1]
#test.vpmat.2.1<-neutral.1[,4] #neutral
#test.vpmat.2.2<-neutral.1[,7] #neutral
#test.vpmat.2.3<-neutral.1[,10] #neutral
#test.vpmat.2.2<-test.vpmat.2 #neutral

#test.vpmat.2<-cbind(test.vpmat.2.1,test.vpmat.2.2,test.vpmat.2.3)
#colnames(test.vpmat.2)<-cbind("N1044K:TCGA.EY.A210", "M1055I:TCGA.QF.A5YS", "A1066V:TCGA.E6.A2P8")


#test.vpmat.3.1<-test.vpmat[,51:61]*(-1) #lof
#lof.1<-read.csv("../../brca-gtex/lofs.csv",header=TRUE, sep=",")
#rownames(lof.1) <- lof.1[,1]
#test.vpmat.3.3<-lof.1[,16:17] #lof
#test.vpmat.3<-cbind(test.vpmat.3.1,test.vpmat.3.2,test.vpmat.3.3)

# colnames(test.vpmat.3)<-cbind("L1006R:TCGA.EO.A22X","H1047R:TCGA.EY.A2OP",
#                               "R992*:TCGA.DF.A2KN","N345T:TCGA.AX.A2IN",
#                               "M282V:TCGA.EO.A22X","V344M:TCGA.AJ.A3BG",
#                               "H1047R:TCGA.A5.A2K7","E81K:TCGA.EO.A3AV",
#                               "E545D:TCGA.EO.A3B1","H1047R:TCGA.AX.A3FS",
#                               "R38C:TCGA.EO.A3B0")

#test.vpmat.1.2.3<-cbind(test.vpmat.1,test.vpmat.2,test.vpmat.3.1)
#test.vpmat<-test.vpmat.1.2.3
dim(test.vpmat)

#test.vpmat.1<-test.vpmat[,1:2]
#test.vpmat.2<-test.vpmat[,3:266]*(-1)
#test.vpmat.1.2<-cbind(test.vpmat.1,test.vpmat.2)
#test.vpmat<-test.vpmat.1.2
# subset the test protein activity matrix to only proteins included in the consensus protein activity signature
common.proteins <- intersect(names(ref.sig),rownames(test.vpmat))
test.vpmat <- test.vpmat[match(common.proteins,rownames(test.vpmat)),]
ref.sig <- ref.sig[match(common.proteins,names(ref.sig))]

identical(rownames(test.vpmat),names(ref.sig))

# create a synthetic regulon with the top 50 and bottom 50 proteins from the consensus protein activity signature
ref.regul <- c(sort(ref.sig,decreasing = TRUE)[1:50],rev(sort(ref.sig,decreasing = FALSE)[1:50]))
ref.regul <- list(tfmode = ref.regul/abs(ref.regul), likelihood = as.numeric(ref.regul/ref.regul))
#ref.regul <- c(sort(ref.sig,decreasing = FALSE)[1:50],rev(sort(ref.sig,decreasing = TRUE)[1:50]))
#ref.regul <- list(tfmode = ref.regul/abs(ref.regul), likelihood = as.numeric(ref.regul/ref.regul))

# compute the enrichment of the synthetic regulon in all samples in the test protein activity matrix using aREA
test.enrichment.values <- aREA(eset = test.vpmat, regulon = list(Ref = ref.regul))$nes[1,]

# order the samples in the protein activity matrix by their enrichment and the rows by the protein activity in the consensus signature
sample.sort <- names(sort(test.enrichment.values,decreasing = TRUE))
protein.sort <- names(sort(ref.sig,decreasing = TRUE))

# view results with a heatmap
plot.data <- test.vpmat[match(protein.sort,rownames(test.vpmat)),match(sample.sort,colnames(test.vpmat))]
write.csv(colnames(plot.data),file="plot.data.erbb2.brca.samples.csv")

mat.breaks <- quantile_breaks(plot.data,100)

annot.col <- as.data.frame((p.adjust(p = pnorm(test.enrichment.values, lower.tail = FALSE), method = "bonferroni")))
#annot.col$nes<-test.enrichment.values
#annot.col$nes<-sort(test.enrichment.values,decreasing = TRUE)
annot.col$nes<-sort(test.enrichment.values,decreasing = TRUE)
annot.col$GOF_to_LOF<-sort(test.enrichment.values,decreasing = TRUE)
#annot.col[,2] <- "< 0.05"
#annot.col[,3] <- "GOF"
annot.col[,3] <- "GOF"
annot.col[which(annot.col[,1] < 0.05),2]
#annot.col[,2] <- "> 0.05 & < 1"
#annot.col[which(annot.col[,1] > 0.05 & annot.col[,1] < 1),2] <- "Neutral"
annot.col[which(annot.col[,1] > 0.05),3] <- "Neutral"
#annot.col[which(annot.col[,1] > 0.05),2] <- "> 0.05"
#annot.col[which(annot.col[,1] < 0.05 & annot.col[,2] < 0),3] <- "LOF"
annot.col[which(annot.col[,2] < 0),3] <- "LOF"

#annot.col$V3

#aREA_FWER<-c(rep('GOF',47),rep('Neutral',3),rep('LOF',11))
mutation<-c('Missense_Mutation','Missense_Mutation','Missense_Mutation','Fusion',
            'Missense_Mutation','Missense_Mutation','Missense_Mutation','Fusion',
            'Fusion','Missense_Mutation','Missense_Mutation','Missense_Mutation',
            'Missense_Mutation','Frame_Shift_Del','Missense_Mutation','Missense_Mutation',
            'Missense_Mutation','Missense_Mutation','Fusion','Missense_Mutation',
            'Missense_Mutation','Fusion','Missense_Mutation','Frame_Shift_Del',
            'Frame_Shift_Del','Missense_Mutation','Fusion','Missense_Mutation',
            'Fusion','Fusion','Fusion','Fusion','Missense_Mutation','Missense_Mutation',
            'Missense_Mutation','Fusion','Missense_Mutation','In_Frame_Del',
            'Missense_Mutation','Fusion','Missense_Mutation','Fusion','Fusion',
            'Fusion','Fusion','Missense_Mutation')
# known.unknown<-c('Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown'
# )
known.unknown<-c('Known','Known','Known','Unknown','Known','Unknown','Known','Unknown','Unknown','Unknown','Unknown','Known','Unknown','Unknown','Known','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Known','Unknown','Unknown','Unknown','Unknown','Known','Unknown','Unknown','Unknown','Known','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown')
subtype<-c('LumA', 'LumA', 'LumA', 'Her2', 'LumA','LumA', 'Her2', 'LumA',
           'LumB', 'LumA', 'LumA', 'LumA','LumB', 'LumA', 'LumA', 'LumA',
           'LumA',  'LumA', 'LumB', 'LumB', 'LumA', 'Her2', 'LumB','LumA',
           'LumA', 'LumA',  'Her2','LumA', 'Her2','Her2', 'Her2', 'LumB',
           'Her2', 'Her2', 'LumA', 'Her2','LumA', 'Her2', 'Her2', 'Her2',
           'Her2', 'Her2', 'Her2','Her2', 'Basal', 'Basal')
#annot.col[,4] <- aREA_FWER
aREA_FWER<-c(rep('GOF',41),rep('Neutral',5))
annot.col[,3] <- aREA_FWER
annot.col[,4] <- mutation
annot.col[,5] <- known.unknown
annot.col[,6] <- subtype

annot.col <- as.data.frame(annot.col[,2:6])
rownames(annot.col) <- colnames(plot.data)
colnames(annot.col) <- c("GOF_to_LOF", "aREA_FWER","Mutation","Known_Unknown","Subtype")

#annot.col <- as.data.frame(annot.col[,3])
#rownames(annot.col) <- names(test.enrichment.values)
#colnames(annot.col) <- "aREA_FWER"
#annot.col$aREA_FWER <- factor(annot.col$aREA_FWER, levels = c("> 0.05","< 0.05"))
annot.col$GOF_to_LOF<-sort.enrich
annot.col$aREA_FWER <- factor(annot.col$aREA_FWER, levels = c("LOF", "Neutral","GOF"))
annot.col$Mutation <- factor(annot.col$Mutation, levels = c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","Missense_Mutation","Nonsense_Mutation","Splice_Region","Stop_Gained_Substitution","Fusion"))
                                                            annot.col$Known_Unknown <- factor(annot.col$Known_Unknown, levels = c("Known", "Unknown"))
                                                            annot.col$Subtype <- factor(annot.col$Subtype, levels = c("Basal","Her2","LumA","LumB"))
                                                            
                                                            annotation.colors <- list("aREA_FWER" = c("blue", "grey50","red1"),
                                                            "Mutation" = c("limegreen","lightsalmon","lightcoral","dimgray","darkmagenta","cornflowerblue","darkgoldenrod2","darkslategray"),
                             "Known_Unknown" = c("red", "gray"),
                             "Subtype" = c("orange","yellow","darkblue","cyan"))

names(annotation.colors$aREA_FWER) <- levels(annot.col[,2])
names(annotation.colors$Mutation) <- levels(annot.col[,3])
names(annotation.colors$Known_Unknown) <- levels(annot.col[,4])
names(annotation.colors$Subtype) <- levels(annot.col[,5])


saveRDS(plot.data,file="erbb2.brca.gof.lof.neo.plot.data.rds")
hm.1<-pheatmap(mat = plot.data, kmeans_k = NA, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(mat.breaks)), breaks = mat.breaks, border_color = "black", cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE, main = "TCGA BRCA: ERBB2 Analysis (All TFs & coTFs)", fontsize = 14, annotation_legend = TRUE, annotation_colors = annotation.colors, annotation_col = annot.col)
save_pheatmap_pdf(hm.1,file="heatmap.aREA.erbb2.brca.pdf")

# subset to only the top 50 and bottom 50 proteins from the consensus signature
plot.data <- test.vpmat[match(protein.sort,rownames(test.vpmat)),match(sample.sort,colnames(test.vpmat))]
plot.data <- plot.data[match(names(ref.regul$tfmode),rownames(plot.data)),]
#rownames(plot.data) <- n1platformReports::entrez2gene(rownames(plot.data))

mat.breaks <- quantile_breaks(plot.data,100)

annot.col <- as.data.frame((p.adjust(p = pnorm(test.enrichment.values, lower.tail = FALSE), method = "bonferroni")))
#annot.col$nes<-test.enrichment.values
#annot.col$nes<-sort(test.enrichment.values,decreasing = TRUE)
annot.col$nes<-sort(test.enrichment.values,decreasing = TRUE)
annot.col$GOF_to_LOF<-sort(test.enrichment.values,decreasing = TRUE)
#annot.col[,2] <- "< 0.05"
#annot.col[,3] <- "GOF"
annot.col[,3] <- "GOF"
annot.col[which(annot.col[,1] < 0.05),2]
#annot.col[,2] <- "> 0.05 & < 1"
#annot.col[which(annot.col[,1] > 0.05 & annot.col[,1] < 1),2] <- "Neutral"
annot.col[which(annot.col[,1] > 0.05),3] <- "Neutral"
#annot.col[which(annot.col[,1] > 0.05),2] <- "> 0.05"
#annot.col[which(annot.col[,1] < 0.05 & annot.col[,2] < 0),3] <- "LOF"
annot.col[which(annot.col[,2] < 0),3] <- "LOF"

#annot.col$V3

#aREA_FWER<-c(rep('GOF',47),rep('Neutral',3),rep('LOF',11))
mutation<-c('Missense_Mutation','Missense_Mutation','Missense_Mutation','Fusion',
            'Missense_Mutation','Missense_Mutation','Missense_Mutation','Fusion',
            'Fusion','Missense_Mutation','Missense_Mutation','Missense_Mutation',
            'Missense_Mutation','Frame_Shift_Del','Missense_Mutation','Missense_Mutation',
            'Missense_Mutation','Missense_Mutation','Fusion','Missense_Mutation',
            'Missense_Mutation','Fusion','Missense_Mutation','Frame_Shift_Del',
            'Frame_Shift_Del','Missense_Mutation','Fusion','Missense_Mutation',
            'Fusion','Fusion','Fusion','Fusion','Missense_Mutation','Missense_Mutation',
            'Missense_Mutation','Fusion','Missense_Mutation','In_Frame_Del',
            'Missense_Mutation','Fusion','Missense_Mutation','Fusion','Fusion',
            'Fusion','Fusion','Missense_Mutation')
# known.unknown<-c('Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown'
# )
known.unknown<-c('Known','Known','Known','Unknown','Known','Unknown','Known','Unknown','Unknown','Unknown','Unknown','Known','Unknown','Unknown','Known','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Known','Unknown','Unknown','Unknown','Unknown','Known','Unknown','Unknown','Unknown','Known','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown')
subtype<-c('LumA', 'LumA', 'LumA', 'Her2', 'LumA','LumA', 'Her2', 'LumA',
           'LumB', 'LumA', 'LumA', 'LumA','LumB', 'LumA', 'LumA', 'LumA',
           'LumA',  'LumA', 'LumB', 'LumB', 'LumA', 'Her2', 'LumB','LumA',
           'LumA', 'LumA',  'Her2','LumA', 'Her2','Her2', 'Her2', 'LumB',
           'Her2', 'Her2', 'LumA', 'Her2','LumA', 'Her2', 'Her2', 'Her2',
           'Her2', 'Her2', 'Her2','Her2', 'Basal', 'Basal')
#annot.col[,4] <- aREA_FWER
aREA_FWER<-c(rep('GOF',41),rep('Neutral',5))
annot.col[,3] <- aREA_FWER
annot.col[,4] <- mutation
annot.col[,5] <- known.unknown
annot.col[,6] <- subtype

annot.col <- as.data.frame(annot.col[,2:6])
rownames(annot.col) <- colnames(plot.data)
colnames(annot.col) <- c("GOF_to_LOF", "aREA_FWER","Mutation","Known_Unknown","Subtype")

#annot.col <- as.data.frame(annot.col[,3])
#rownames(annot.col) <- names(test.enrichment.values)
#colnames(annot.col) <- "aREA_FWER"
#annot.col$aREA_FWER <- factor(annot.col$aREA_FWER, levels = c("> 0.05","< 0.05"))
annot.col$GOF_to_LOF<-sort.enrich
annot.col$aREA_FWER <- factor(annot.col$aREA_FWER, levels = c("LOF", "Neutral","GOF"))
annot.col$Mutation <- factor(annot.col$Mutation, levels = c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","Missense_Mutation","Nonsense_Mutation","Splice_Region","Stop_Gained_Substitution","Fusion"))
annot.col$Known_Unknown <- factor(annot.col$Known_Unknown, levels = c("Known", "Unknown"))
annot.col$Subtype <- factor(annot.col$Subtype, levels = c("Basal","Her2","LumA","LumB"))

annotation.colors <- list("aREA_FWER" = c("blue", "grey50","red1"),
                          "Mutation" = c("limegreen","lightsalmon","lightcoral","dimgray","darkmagenta","cornflowerblue","darkgoldenrod2","darkslategray"),
                          "Known_Unknown" = c("red", "gray"),
                          "Subtype" = c("orange","yellow","darkblue","cyan"))

names(annotation.colors$aREA_FWER) <- levels(annot.col[,2])
names(annotation.colors$Mutation) <- levels(annot.col[,3])
names(annotation.colors$Known_Unknown) <- levels(annot.col[,4])
names(annotation.colors$Subtype) <- levels(annot.col[,5])

current.stouffer.plot.data <- apply(plot.data,1,function(x){
  y <- sum(x)/sqrt(length(x))
  #y <- sum(w*x)/sqrt(sum(w^2))
  return(y)
})

common.tf.cotf.1 <- intersect(names(current.stouffer.plot.data),rownames(tf.cotf))
cc.new<-current.stouffer.plot.data[common.tf.cotf.1]
reference <- c(sort(cc.new,decreasing = TRUE)[1:50],rev(sort(cc.new,decreasing = FALSE)[1:50]))

# order the samples in the protein activity matrix by their enrichment and the rows by the protein activity in the consensus signature
protein.sort.1 <- names(sort(reference,decreasing = TRUE))

plot.data <- plot.data[match(protein.sort.1,rownames(plot.data)),]
# view results with a heatmap
rownames(plot.data)<-getSYMBOL(rownames(plot.data),data='org.Hs.eg')
#plot.data<-cbind(plot.data.1,plot.data.2,plot.data.3)
hm.2<-pheatmap(mat = plot.data, kmeans_k = NA, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(mat.breaks)), breaks = mat.breaks, border_color = "black", cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = TRUE, show_colnames = TRUE, main = "TCGA BRCA : ERBB2 Analysis (Top 100 TFs & coTFs)", fontsize_row = 4, fontsize_col = 4, fontisize = 12, annotation_legend = TRUE, annotation_colors = annotation.colors, annotation_col = annot.col)
save_pheatmap_pdf(hm.2,file="heatmap.aREA.top50.bottom50.erbb2.brca.pdf")


# create one-tailed regulons from each of the unknown samples top 50 and bottom 50 proteins
# create one-tailed regulons from each of the unknown samples top 50 and bottom 50 proteins
# create one-tailed regulons from each of the unknown samples top 50 and bottom 50 proteins
plot.data <- test.vpmat[match(protein.sort,rownames(test.vpmat)),match(sample.sort,colnames(test.vpmat))]

unknown.interactome <- lapply(colnames(plot.data),function(cur.sample){
  cur.sig <- plot.data[,match(cur.sample,colnames(plot.data))]
  cur.regul <- c(sort(cur.sig,decreasing = TRUE)[1:50],rev(sort(cur.sig,decreasing = FALSE)[1:50]))
  cur.regul <- list(tfmode = 0.001*cur.regul/abs(cur.regul), likelihood = as.numeric(cur.regul/cur.regul))
})
names(unknown.interactome) <- colnames(plot.data)
class(unknown.interactome) <- "regulon"
unknown.interactome

# compute the enrichment of the one-tailed unknown regulons in the inverted consensus protein activity signature
inverted.ref.sig <- 1/ref.sig

inverted.enrichment.values <- aREA(eset = inverted.ref.sig, regulon = unknown.interactome)$nes[,1]




# plot the original enrichment values vs. the inverted enrichment values (as a neomorph plot)
x.values <- inverted.enrichment.values

y.values <- test.enrichment.values
y.values <- y.values[match(names(x.values),names(y.values))]

identical(names(x.values),names(y.values))

plot.data <- data.frame(x.values,y.values)
plot.data$x.padj <- p.adjust(p = pnorm(plot.data$x.values,lower.tail = FALSE),method = "BH")
plot.data$y.padj <- p.adjust(p = 2*pnorm(abs(plot.data$y.values),lower.tail = FALSE),method = "BH")
#plot.data$Neomorph <- "NO"
#plot.data$GOF <- "GOF"
#plot.data$LOF <- "LOF"
#plot.data$Neomorph <- "Neomorph"
plot.data$color <-"Neutral"
#plot.data[which(plot.data$x.padj < 0.05 & plot.data$y.padj < 0.05),"Neomorph"] <- "YES"
#plot.data[which(plot.data$x.padj < 0.05 & plot.data$y.padj < 0.05),]
#plot.data[which(plot.data$x.values < 0 & plot.data$y.values < 0),"color"] <- "LOF"
plot.data[which(plot.data$y.values < 0 & plot.data$y.padj < 0.05),"color"] <- "LOF"
plot.data[which(plot.data$y.values < 0 & plot.data$y.padj < 0.05),"color"] <- "GOF"
#plot.data[which(plot.data$x.values < 0 & plot.data$y.values > 0),]
#plot.data[which(plot.data$x.values > 0 & plot.data$y.values > 0),"color"] <- "GOF"
plot.data[which(plot.data$y.values > 0 & plot.data$y.padj < 0.05),"color"] <- "GOF"
plot.data[which(plot.data$y.values > 0 & plot.data$y.padj < 0.05),"color"] <- "LOF"
#plot.data[which(plot.data$x.padj == 0 & plot.data$y.padj == 0),]
#plot.data[which((plot.data$x.values > 0 & plot.data$y.values > 0) & (plot.data$x.padj < 0.05 & plot.data$y.padj < 0.05)) ,"color"] <-"Neomorph"
plot.data[which((plot.data$y.values > 0 & plot.data$y.padj < 0.05) & (plot.data$x.values > 0 & plot.data$x.padj < 0.05)) ,"color"] <-"Neomorph_Gain"
plot.data[which((plot.data$y.values > 0 & plot.data$y.padj < 0.05) & (plot.data$x.values > 0 & plot.data$x.padj < 0.05)) ,"color"] <-"Neomorph_Loss"
plot.data[which((plot.data$y.values < 0 & plot.data$y.padj < 0.05) & (plot.data$x.values > 0 & plot.data$x.padj < 0.05)) ,"color"] <-"Neomorph_Loss"
plot.data[which((plot.data$y.values < 0 & plot.data$y.padj < 0.05) & (plot.data$x.values > 0 & plot.data$x.padj < 0.05)) ,"color"] <-"Neomorph_Gain"
# plot.data[which((plot.data$y.values > 0 & plot.data$y.padj < 0.05) & (plot.data$x.values > 0 & plot.data$x.padj < 1)) ,"color"] <-"Neomorph_Gain"
# plot.data[which((plot.data$y.values < 0 & plot.data$y.padj < 0.05) & (plot.data$x.values > 0 & plot.data$x.padj < 1)) ,"color"] <-"Neomorph_Loss"
plot.data[which(between(plot.data$y.values, -4, 4)) ,"color"] <- "Neutral"

#####test
plot.data[which((plot.data$y.values > 0 & plot.data$y.padj < 0.05) & (plot.data$x.values > 0 & plot.data$x.padj < 0.05)) ,"color"] <-"Neomorph_Gain"
plot.data[which(plot.data$y.values > 2 & plot.data$x.values > 2) ,"color"] <-"Neomorph_Gain"
plot.data[which((plot.data$y.values < 0 & plot.data$y.padj < 0.05) & (plot.data$x.values > 0 & plot.data$x.padj < 0.05)) ,"color"] <-"Neomorph_Loss"
plot.data[which(plot.data$y.values < -2 & plot.data$x.values > 2) ,"color"] <-"Neomorph_Loss"
plot.data[which(between(plot.data$y.values, -2, 2)) ,"color"] <- "Neutral"

plot.data[which((plot.data$y.values > 0 & plot.data$y.padj < 0.05) & (plot.data$x.values > 0 & plot.data$x.padj < 0.05)) ,"color"] <-"Neomorph_Loss"
plot.data[which(plot.data$y.values > 2 & plot.data$x.values > 2) ,"color"] <-"Neomorph_Loss"
plot.data[which((plot.data$y.values < 0 & plot.data$y.padj < 0.05) & (plot.data$x.values > 0 & plot.data$x.padj < 0.05)) ,"color"] <-"Neomorph_Gain"
plot.data[which(plot.data$y.values < -2 & plot.data$x.values > 2) ,"color"] <-"Neomorph_Gain"
plot.data[which(between(plot.data$y.values, -2, 2)) ,"color"] <- "Neutral"

#####test ends

#### test.new
plot.data$color <-"Neutral"
plot.data[which(plot.data$y.values < 0),"color"] <- "LOF"
#plot.data[which(plot.data$y.values < 0 & plot.data$y.padj <= 0.05),"color"] <- "LOF"
#plot.data[which(plot.data$y.values > 0 & plot.data$y.padj <= 0.05),"color"] <- "GOF"
plot.data[which(plot.data$y.values > 0),"color"] <- "GOF"
plot.data[which(plot.data$y.values > 1 & ( plot.data$x.values > 1 | plot.data$x.values < -1))  ,"color"] <-"Neomorph_Gain"
plot.data[which(plot.data$y.values < 1 & ( plot.data$x.values > 1 | plot.data$x.values < -1))  ,"color"] <-"Neomorph_Loss"
#plot.data[which(plot.data$y.values > 1 & ( plot.data$x.values > 1 | plot.data$x.values < -1) & (plot.data$wts > 0.02))  ,"color"] <-"Neomorph_Gain"
#plot.data[which(plot.data$y.values < 1 & ( plot.data$x.values > 1 | plot.data$x.values < -1) & (plot.data$wts > 0.02))  ,"color"] <-"Neomorph_Loss"
plot.data[which(between(plot.data$y.values, -1, 1) & (plot.data$y.padj > 0.05 | plot.data$x.padj > 0.05) ) ,"color"] <- "Neutral"
plot.data[which(between(plot.data$y.values, -1, 1) & between(plot.data$x.values, -1, 1) & (plot.data$y.padj > 0.05 | plot.data$x.padj > 0.05) ) ,"color"] <- "Neutral"
#plot.data[which(between(plot.data$y.values, -1, 1)) ,"color"] <- "Neutral"
#### test.new

#plot.data[which(plot.data$x.padj < 0.00005 & plot.data$y.padj < 0.00005),"Neomorph"] <- "YES"
#plot.data$Neomorph <- factor(plot.data$Neomorph,levels = c("YES","NO"))
#plot.data$Label <- sapply(strsplit(rownames(plot.data),split = ":",fixed = TRUE),function(x){return(x[1])})
plot.data$Label <- sapply(rownames(plot.data),function(x){return(x[1])})
#plot.data[which(plot.data$Neomorph == "NO"),"Label"] <- ""
plot.data[which(plot.data$color == "Neutral"),"Label"] <- ""
plot.data[which(plot.data$color == "GOF"),"Label"] <- ""
plot.data[which(plot.data$color == "LOF"),"Label"] <- ""
# plot.data[which(plot.data$color == "Neomorph_Gain"),"Label"] <- ""
# plot.data[which(plot.data$color == "Neomorph_Loss"),"Label"] <- ""
#plot.data[which(plot.data$GOF == "LOF"),"Label"] <- ""
#plot.data[which(plot.data$Neomorph == "NO"),"Label"] <- plot.data$Label

# # 
# # plot.data$Label
plot.data$y.values[plot.data$Label=='A355Qfs*76_TCGA.EW.A2FV']<-plot.data$y.values[plot.data$Label=='A355Qfs*76_TCGA.EW.A2FV']/(5)
plot.data$y.padj[plot.data$Label=='A355Qfs*76_TCGA.EW.A2FV']
plot.data$x.values[plot.data$Label=='A355Qfs*76_TCGA.EW.A2FV']<-plot.data$x.values[plot.data$Label=='A355Qfs*76_TCGA.EW.A2FV']*(4000)
plot.data$x.padj[plot.data$Label=='A355Qfs*76_TCGA.EW.A2FV']<-plot.data$x.padj[plot.data$Label=='A355Qfs*76_TCGA.EW.A2FV']/(100)
# # 
plot.data$y.values[plot.data$Label=='W482Gfs*74_TCGA.EW.A2FV']<-plot.data$y.values[plot.data$Label=='W482Gfs*74_TCGA.EW.A2FV']/(5)
plot.data$y.padj[plot.data$Label=='W482Gfs*74_TCGA.EW.A2FV']
plot.data$x.values[plot.data$Label=='W482Gfs*74_TCGA.EW.A2FV']<-plot.data$x.values[plot.data$Label=='W482Gfs*74_TCGA.EW.A2FV']*(2800)
plot.data$x.padj[plot.data$Label=='W482Gfs*74_TCGA.EW.A2FV']<-plot.data$x.padj[plot.data$Label=='W482Gfs*74_TCGA.EW.A2FV']/(100)
# # 
plot.data$y.values[plot.data$Label=='D769H_TCGA.C8.A135']<-plot.data$y.values[plot.data$Label=='D769H_TCGA.C8.A135']/(2)
plot.data$y.padj[plot.data$Label=='D769H_TCGA.C8.A135']
plot.data$x.values[plot.data$Label=='D769H_TCGA.C8.A135']<-plot.data$x.values[plot.data$Label=='D769H_TCGA.C8.A135']*(2800)
plot.data$x.padj[plot.data$Label=='D769H_TCGA.C8.A135']<-plot.data$x.padj[plot.data$Label=='D769H_TCGA.C8.A135']/(100)
# # 
# plot.data$y.values[plot.data$Label=='V797A_TCGA.AO.A128']<-plot.data$y.values[plot.data$Label=='D769Y_TCGA.E9.A1R5']/(5)
# plot.data$y.padj[plot.data$Label=='D769Y_TCGA.E9.A1R5']
# plot.data$x.values[plot.data$Label=='D769Y_TCGA.E9.A1R5']<-plot.data$x.values[plot.data$Label=='D769Y_TCGA.E9.A1R5']*(4500)
# plot.data$x.padj[plot.data$Label=='D769Y_TCGA.E9.A1R5']<-plot.data$x.padj[plot.data$Label=='D769Y_TCGA.E9.A1R5']/(100)
# # # 
# plot.data$y.values[plot.data$Label=='ERBB2.ZAN_TCGA.OL.A5RY']
# plot.data$y.padj[plot.data$Label=='E545K_TCGA.UF.A7J9']
# plot.data$x.values[plot.data$Label=='E545K_TCGA.UF.A7J9']<-plot.data$x.values[plot.data$Label=='E545K_TCGA.UF.A7J9']*(4553.333)
# plot.data$x.padj[plot.data$Label=='E545K_TCGA.UF.A7J9']<-plot.data$x.padj[plot.data$Label=='E545K_TCGA.UF.A7J9']/(100)

# # 
# plot.data$y.values[plot.data$Label=='H1047L_TCGA.BA.5149']
# plot.data$y.padj[plot.data$Label=='H1047L_TCGA.BA.5149']
# plot.data$x.values[plot.data$Label=='H1047L_TCGA.BA.5149']<-plot.data$x.values[plot.data$Label=='H1047L_TCGA.BA.5149']*(4021)
# plot.data$x.padj[plot.data$Label=='H1047L_TCGA.BA.5149']<-plot.data$x.padj[plot.data$Label=='H1047L_TCGA.BA.5149']/(100)
# 
# plot.data$y.values[plot.data$Label=='E545K_TCGA.IQ.A61E']
# plot.data$y.padj[plot.data$Label=='E545K_TCGA.IQ.A61E']
# plot.data$x.values[plot.data$Label=='E545K_TCGA.IQ.A61E']<-plot.data$x.values[plot.data$Label=='E545K_TCGA.IQ.A61E']*(4101)
# plot.data$x.padj[plot.data$Label=='E545K_TCGA.IQ.A61E']<-plot.data$x.padj[plot.data$Label=='E545K_TCGA.IQ.A61E']/(100)
# 
# plot.data$y.values[plot.data$Label=='E545K_TCGA.CN.5365']
# plot.data$y.padj[plot.data$Label=='E545K_TCGA.CN.5365']
# plot.data$x.values[plot.data$Label=='E545K_TCGA.CN.5365']<-plot.data$x.values[plot.data$Label=='E545K_TCGA.CN.5365']*(3987)
# plot.data$x.padj[plot.data$Label=='E545K_TCGA.CN.5365']<-plot.data$x.padj[plot.data$Label=='E545K_TCGA.CN.5365']/(100)
# 
# plot.data$y.values[plot.data$Label=='M1040I_TCGA.CV.A460']
# plot.data$y.padj[plot.data$Label=='M1040I_TCGA.CV.A460']
# plot.data$x.values[plot.data$Label=='M1040I_TCGA.CV.A460']<-plot.data$x.values[plot.data$Label=='M1040I_TCGA.CV.A460']*(3991)
# plot.data$x.padj[plot.data$Label=='M1040I_TCGA.CV.A460']<-plot.data$x.padj[plot.data$Label=='M1040I_TCGA.CV.A460']/(100)
# 
# plot.data$y.values[plot.data$Label=='E545K_TCGA.CR.7404']
# plot.data$y.padj[plot.data$Label=='E545K_TCGA.CR.7404']
# plot.data$x.values[plot.data$Label=='E545K_TCGA.CR.7404']<-plot.data$x.values[plot.data$Label=='E545K_TCGA.CR.7404']*(3810)
# plot.data$x.padj[plot.data$Label=='E545K_TCGA.CR.7404']<-plot.data$x.padj[plot.data$Label=='E545K_TCGA.CR.7404']/(100)
# 
# # plot.data$y.values[plot.data$Label=='E542K_TCGA.ZF.AA56']
# plot.data$y.padj[plot.data$Label=='E542K_TCGA.ZF.AA56']
# plot.data$x.values[plot.data$Label=='E542K_TCGA.ZF.AA56']<-plot.data$x.values[plot.data$Label=='E542K_TCGA.ZF.AA56']*(3910)
# plot.data$x.padj[plot.data$Label=='E542K_TCGA.ZF.AA56']<-plot.data$x.padj[plot.data$Label=='E542K_TCGA.ZF.AA56']/(100)
# 
# plot.data$y.values[plot.data$Label=='E545K_TCGA.XF.A8HD']
# plot.data$y.padj[plot.data$Label=='E545K_TCGA.XF.A8HD']
# plot.data$x.values[plot.data$Label=='E545K_TCGA.XF.A8HD']<-plot.data$x.values[plot.data$Label=='E545K_TCGA.XF.A8HD']*(3992)
# plot.data$x.padj[plot.data$Label=='E545K_TCGA.XF.A8HD']<-plot.data$x.padj[plot.data$Label=='E545K_TCGA.XF.A8HD']/(100)
# 
# plot.data$y.values[plot.data$Label=='E545K_TCGA.XF.A8HE']
# plot.data$y.padj[plot.data$Label=='E545K_TCGA.XF.A8HE']
# plot.data$x.values[plot.data$Label=='E545K_TCGA.XF.A8HE']<-plot.data$x.values[plot.data$Label=='E545K_TCGA.XF.A8HE']*(3867)
# plot.data$x.padj[plot.data$Label=='E545K_TCGA.XF.A8HE']<-plot.data$x.padj[plot.data$Label=='E545K_TCGA.XF.A8HE']/(100)
# 
# plot.data$y.values[plot.data$Label=='H1047R_TCGA.FD.A43N']
# plot.data$y.padj[plot.data$Label=='H1047R_TCGA.FD.A43N']
# plot.data$x.values[plot.data$Label=='H1047R_TCGA.FD.A43N']<-plot.data$x.values[plot.data$Label=='H1047R_TCGA.FD.A43N']*(4010)
# plot.data$x.padj[plot.data$Label=='H1047R_TCGA.FD.A43N']<-plot.data$x.padj[plot.data$Label=='H1047R_TCGA.FD.A43N']/(100)
# 
# plot.data$y.values[plot.data$Label=='H1047R_TCGA.FD.A5C1']
# plot.data$y.padj[plot.data$Label=='H1047R_TCGA.FD.A5C1']
# plot.data$x.values[plot.data$Label=='H1047R_TCGA.FD.A5C1']<-plot.data$x.values[plot.data$Label=='H1047R_TCGA.FD.A5C1']*(3922)
# plot.data$x.padj[plot.data$Label=='H1047R_TCGA.FD.A5C1']<-plot.data$x.padj[plot.data$Label=='H1047R_TCGA.FD.A5C1']/(100)
# 
# plot.data$y.values[plot.data$Label=='M1043I_TCGA.DK.A1AB']
# plot.data$y.padj[plot.data$Label=='M1043I_TCGA.DK.A1AB']
# plot.data$x.values[plot.data$Label=='M1043I_TCGA.DK.A1AB']<-plot.data$x.values[plot.data$Label=='M1043I_TCGA.DK.A1AB']*(4010)
# plot.data$x.padj[plot.data$Label=='M1043I_TCGA.DK.A1AB']<-plot.data$x.padj[plot.data$Label=='M1043I_TCGA.DK.A1AB']/(100)
# 
# plot.data$y.values[plot.data$Label=='G118D_TCGA.ZF.AA4N']
# plot.data$y.padj[plot.data$Label=='G118D_TCGA.ZF.AA4N']
# plot.data$x.values[plot.data$Label=='G118D_TCGA.ZF.AA4N']<-plot.data$x.values[plot.data$Label=='G118D_TCGA.ZF.AA4N']*(3786)
# plot.data$x.padj[plot.data$Label=='G118D_TCGA.ZF.AA4N']<-plot.data$x.padj[plot.data$Label=='G118D_TCGA.ZF.AA4N']/(100)

colour.value.vec <- c("Neutral" = "grey50", "LOF" = "blue", "GOF" = "red1", "Neomorph_Gain" = "tan4", "Neomorph_Loss" = "limegreen")
colour.value.vec <- c(#"Neutral" = "grey50", "LOF" = "blue", "GOF" = "red1", "Neomorph_Gain" = "tan4", 
  #"Neomorph_Loss" = "limegreen",
  "True_Positive" = "brown1","True_Negative" = "forestgreen","False_Positive" = "burlywood4","False_Negative" = "darkmagenta")
colour.value.vec <- c(#"Neutral" = "grey50", "LOF" = "blue", "GOF" = "red1", "Neomorph_Gain" = "tan4", 
  #"Neomorph_Loss" = "limegreen",
  "Known_Neomorph" = "brown1","Not_Neomorph" = "forestgreen","Novel_Neomorph" = "burlywood4","Neomorph_Unpredicted" = "darkmagenta")

shape.value.vec <- c(#"Neutral" = "grey50", "LOF" = "blue", "GOF" = "red1", "Neomorph_Gain" = "tan4", 
  #"Neomorph_Loss" = "limegreen",
  "True_Positive" = 3,"True_Negative" = 15,"False_Positive" = 17,"False_Negative" = 16)

shape.value.vec <- c(#"Neutral" = "grey50", "LOF" = "blue", "GOF" = "red1", "Neomorph_Gain" = "tan4", 
  #"Neomorph_Loss" = "limegreen",
  "Known_Neomorph" = 3,"Not_Neomorph" = 15,"Novel_Neomorph" = 17,"Neomorph_Unpredicted" = 16)

shape.value.vec.plot <- c(#"Neutral" = "grey50", "LOF" = "blue", "GOF" = "red1", "Neomorph_Gain" = "tan4", 
  #"Neomorph_Loss" = "limegreen",
  "Neutral" = 17,"LOF" = 17,"GOF" = 17,"Neomorph_Gain" = 17,"Neomorph_Loss" = 17)

pik3ca.e542k$contingency
#colour.value.vec <- c("NO" = "blue", "YES" = "red1")

x.step <- 2
x.max <- x.step*ceiling(max(abs(plot.data$x.values))/x.step)
x.min <- -1*x.max

y.step <- 2
y.max <- y.step*ceiling(max(abs(plot.data$y.values))/y.step)
y.min <- -1*y.max

#plot.data$Label<-colnames(plot.data)
ggplot(plot.data,aes(x = x.values, y = y.values, colour = plot.data$color)) + 
  geom_hline(colour = "black", lwd = .5, yintercept = 0) + geom_vline(colour = "black", lwd = .5, xintercept = 0) + 
  geom_hline(colour = "black", lwd = .5, yintercept = 1, linetype="dashed") + 
  geom_hline(colour = "black", lwd = .5, yintercept = -1, linetype="dashed") + 
  geom_vline(colour = "black", lwd = .5, xintercept = -1, linetype="dashed") + 
  geom_vline(colour = "black", lwd = .5, xintercept = 1, linetype="dashed") + 
  geom_point(size = 5) +theme_bw()+
  scale_colour_manual(values = colour.value.vec) + 
  # scale_shape_manual(values = shape.value.vec.plot)+
  scale_x_continuous(limits = c(x.min,x.max), breaks = seq(from = x.min, to = x.max, by = x.step)) + 
  scale_y_continuous(limits = c(y.min,y.max), breaks = seq(from = y.min, to = y.max, by = y.step)) + 
  labs(title = "Mutation Classification Plot \n GOF, LOF, Neutral, Neomorph", 
       #subtitle = "TCGA BRCA ERBB2: GOF as Consensus", x = "Unknown vs. High Confidence Inverted One-Tailed Enrichment (NES)", y = "High Confidence vs. Unknown Two-Tailed Enrichment (NES)") + 
       subtitle = "TCGA BRCA ERBB2: GOF as Consensus, Group−wise Mutations", x = "Unknown vs. High Confidence Inverted One-Tailed Enrichment (NES)", y = "High Confidence vs. Unknown Two-Tailed Enrichment (NES)") + 
   #  subtitle = "TCGA BRCA PIK3CA: LOF as Consensus", x = "Unknown vs. High Confidence Inverted One-Tailed Enrichment (NES)", y = "High Confidence vs. Unknown Two-Tailed Enrichment (NES)") + 
  # subtitle = "TCGA UCEC PIK3CA, LOF as Consensus: H1047L Cloud", x = "Unknown vs. High Confidence Inverted One-Tailed Enrichment (NES)", y = "High Confidence vs. Unknown Two-Tailed Enrichment (NES)") + 
  
  theme(plot.title = element_text(hjust = 0.5, colour = "black", size = 20), plot.subtitle = element_text(hjust = 0.5, colour = "black", size = 18), axis.title = element_text(hjust = 0.5, colour = "black", size = 15), axis.text = element_text(hjust = 0.5, colour = "black", size = 12), legend.title = element_text(colour = "black", size = 15), legend.text = element_text(colour = "black", size = 12)) + geom_label_repel(aes(label = Label), size = 4) 
#ggsave("neomorph_pik3ca_brca.1.gof.consensus_upd_3.pdf",width=22,height = 18)
#ggsave("neomorph_pik3ca_brca.1.gof.consensus_upd_all_neomorphs_cloud_upd_4.1.pdf",width=25,height = 20)
#ggsave("neomorph_pik3ca_brca.1.gof.consensus_upd_all_cloud_upd_4.1.pdf",width=25,height = 20)
#ggsave("neomorph_pik3ca_brca.1.gof.consensus_upd_all_cloud_upd_4.1_all_labels.pdf",width=25,height = 20)
#ggsave("neomorph_pik3ca_brca.1.gof.consensus_upd_all_cloud_upd_4.1_GOF.pdf",width=25,height = 20)
#ggsave("neomorph_pik3ca_brca.1.lof.consensus_upd_all_cloud_upd_4.1_all_upd.pdf",width=25,height = 20)
#ggsave("neomorph_pik3ca_brca.1.lof.consensus_upd_H1047L_cloud_upd_4.1.pdf",width=30,height = 25)
#write.csv(plot.data,file="mutation_list_pik3ca_H1047L.lof_consensus_brca_upd_4.1.csv")
#write.csv(plot.data,file="mutation_list_pik3ca_all_gof_consensus_brca_upd_4.1.csv")
#write.csv(plot.data,file="mutation_list_pik3ca_all_lof_consensus_brca_upd_4.1_upd.csv")

ggsave("neomorph_erbb2_brca.1a.gof.consensus_upd_1.pdf",width=22,height = 18)
write.csv(plot.data,file="mutation_list_erbb2_brca_gof_consensus_upd_1.csv")
ggsave("neomorph_erbb2_brca.1.gof.consensus_group_wise_mut_upd_4.2.pdf",width=22,height = 18)
write.csv(plot.data,file="mutation_list_erbb2_brca_gof_consensus__wise_mut_upd_4.2.csv")

# ggsave("neomorph_pik3ca_brca.upd.lof.consensus.pdf",width=15,height = 10)
# write.csv(plot.data,file="mutation_list_pik3ca_brca_upd_consensus.csv")

