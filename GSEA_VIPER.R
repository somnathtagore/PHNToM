# 2022/12/13 GSEA-VIPER script ftom Somnath modified by AN

# Previous documentation
# Name: aREA_PHaNToM Analysis Script
# Author: Somnath Tagore
# Date Created: 03/03/2020
# Date Last Modified: 01/24/2021
# Description: aREA enrichment analysis plot (Signature vs Geneset or vice-versa)


# a function to do a ggplot version of the aREA GSEA plot
GSEA_viper_plot <- function(ref.sig,test.sig,uniform.wts = FALSE,gsea.score = 1,tcm.size,color.positive = "salmon", color.negative = "cornflowerblue",main.title = "Analytical Rank-Based Enrichment Analysis", y.label = "KS Enrichment Score", x.label = "Signature", x.step = 500, y.step = 0.1,gsea.score = 1){
  
  # load packages
  library(ggplot2)
  library(viper)
  
  # specify the gsea score to use (1 is a KS test)
  gsea.score.value <- gsea.score
  
  # compute the aREA enrichment
  test.regul <- list(tfmode = sign(test.sig), likelihood = as.numeric(abs(test.sig)))
  
  if(!is.null(tcm.size)){
    subset.targets <- test.regul$tfmode*test.regul$likelihood
    subset.targets <- c(sort(subset.targets,decreasing = TRUE)[seq(from = 1, to = round(tcm.size/2), by = 1)],rev(sort(subset.targets,decreasing = FALSE)[seq(from = 1, to = round(tcm.size/2), by = 1)])) 
    # modify the geneset if desired
    if(uniform.wts){
      subset.targets <- sign(subset.targets)
    }
    new.regul <- test.regul
    new.regul$tfmode <- sign(subset.targets)
    new.regul$likelihood <- as.numeric(abs(subset.targets))
    test.geneset.regul <- list(new.regul)
    class(test.geneset.regul) <- "regulon"
    test.sig.req <- new.regul$tfmode
  } else {
    subset.targets <- test.regul$tfmode*test.regul$likelihood
    # modify the geneset if desired
    if(uniform.wts){
      subset.targets <- sign(subset.targets)
    }
    new.regul <- test.regul
    new.regul$tfmode <- sign(subset.targets)
    new.regul$likelihood <- as.numeric(abs(subset.targets))
    test.geneset.regul <- list(new.regul)
    class(test.geneset.regul) <- "regulon"
    test.sig.req <- new.regul$tfmode
  }
  
  current.enrich.area <- aREA(eset = ref.sig, regulon = test.geneset.regul)
  
  # define the classical gsea scoring function
  gsea.es <- function(rlist, x, score) {
    nr <- sum(abs(rlist[x])^score)
    nh <- length(rlist)-length(x)
    es <- rep(-(1/nh),length(rlist))
    es[x] <- abs(rlist[x])^score/nr
    return(cumsum(es))
  }
  
  # compute positive target enrichment scores and ranks
  ref.sig <- sort(ref.sig)
  up.es <- gsea.es(rlist = (ref.sig), x = which(names(ref.sig) %in% names(test.sig.req[test.sig.req >= 0])), score = gsea.score.value)
  up.es <- -1*up.es
  up.targets <- test.sig.req[which(test.sig.req >= 0)]
  up.target.ranks <- match(names(up.targets),names(ref.sig))
  
  # compute negative target enrichment scores and ranks
  down.es <- gsea.es(rlist = (ref.sig), x = which(names(ref.sig) %in% names(test.sig.req[test.sig.req < 0])), score = gsea.score.value)
  down.es <- -1*down.es
  down.targets <- test.sig.req[which(test.sig.req < 0)]
  down.target.ranks <- match(names(down.targets),names(ref.sig))
  
  # identify the absolute maximal positive target enrichment score and negative target enrichment score
  max.up.es.value <- max(abs(up.es))
  max.up.es.rank <- which(abs(up.es) == max.up.es.value)
  max.up.es.value <- up.es[max.up.es.rank]
  
  max.down.es.value <- max(abs(down.es))
  max.down.es.rank <- which(abs(down.es) == max.down.es.value)
  max.down.es.value <- down.es[max.down.es.rank]
  
  # plot the gsea enrichment plot
  plot.data <- data.frame(x.values = rep(seq(from = 0, to = (length(up.es) + 1), by = 1),times = 2), y.values = c(0,up.es,0,0,down.es,0), type = c(rep("UP", times = length(up.es)+2),rep("DOWN",times = length(down.es)+2)))
  plot.data$type <- factor(plot.data$type, levels = c("UP","DOWN"))
  #plot.data$up<-up.target.ranks
  #plot.data$down<-down.target.ranks
  plot.data.geneset <- data.frame(x.values = c(up.target.ranks,down.target.ranks), y.start = c(rep(1.2,times = length(up.target.ranks)),rep(1.2,times = length(down.target.ranks))), y.end = c(rep(1.3,times = length(up.target.ranks)),rep(1.1,times = length(down.target.ranks))), type = c(rep("UP",times = length(up.target.ranks)),rep("DOWN",times = length(down.target.ranks))))
  plot.data.geneset$type <- factor(plot.data.geneset$type, levels = c("UP","DOWN"))
  
  if(length(up.targets) == 0){
    plot.data.geneset$opacity <- c((as.numeric(abs(down.targets))/max(abs(down.targets))))
  } else if(length(down.targets) == 0){
    plot.data.geneset$opacity <- c((as.numeric(up.targets)/max(up.targets)))
  } else {
    plot.data.geneset$opacity <- c((as.numeric(up.targets)/max(up.targets)),(as.numeric(abs(down.targets))/max(abs(down.targets))))
  }
  
  area.es.value <- signif(current.enrich.area$es,3)
  area.nes.value <- signif(current.enrich.area$nes,3)
  area.log10.pvalue <- (log(2) + pnorm(q = abs(current.enrich.area$nes), lower.tail = FALSE, log.p = TRUE))/log(10)
  area.pvalue.power <- floor(area.log10.pvalue)
  area.pvalue.num <- signif((10^(area.log10.pvalue - floor(area.log10.pvalue))),3)
  
  
  # current.subtitle <- paste("aREA ES = ",area.es.value," : aREA NES = ",area.nes.value," : aREA p-value = ",area.pvalue.num,"e",area.pvalue.power,sep = "", " : DETOR score = ",round((tcm.size-(tcm.size-(tcm.size-(length(up.target.ranks[up.target.ranks<max.up.es.rank])+ length(down.target.ranks[down.target.ranks>max.down.es.rank])))))/(tcm.size-(tcm.size-(length(up.target.ranks[up.target.ranks<max.up.es.rank])+ length(down.target.ranks[down.target.ranks>max.down.es.rank])))),digits=2)," up.rank: ",max.up.es.rank," down.rank: ",max.down.es.rank)
  # ledge_1 <- (length(test.regul$tfmode)-max.up.es.rank)+max.down.es.rank
  # ledge_c <- length(test.regul$tfmode)-((length(test.regul$tfmode)-max.up.es.rank)+max.down.es.rank)
  # b_detor <- tcm.size-(tcm.size-(tcm.size-(tcm.size-(length(up.target.ranks[up.target.ranks<max.up.es.rank])+ length(down.target.ranks[down.target.ranks>max.down.es.rank])))))
  # a_plus_c_detor <- tcm.size-(tcm.size-(tcm.size-(length(up.target.ranks[up.target.ranks<max.up.es.rank])+ length(down.target.ranks[down.target.ranks>max.down.es.rank]))))
  # up_test<-up.target.ranks[up.target.ranks<max.up.es.rank]
  # down_test<-down.target.ranks[down.target.ranks>max.down.es.rank]
  # not_up_test<-up.target.ranks[up.target.ranks>max.up.es.rank]
  # not_down_test<-down.target.ranks[down.target.ranks<max.down.es.rank]
  # 
  # write.csv(up_test,file=paste0('up_test.csv'))
  # write.csv(down_test,file=paste0("down_test.csv"))
  # write.csv(not_up_test,file=paste0("not_up_test.csv"))
  # write.csv(not_down_test,file=paste0("not_down_test.csv"))
  # write.csv(plot.data.geneset,file=paste0("plot.data.geneset.csv"))
  # write.csv(plot.data,file=paste0("plot.data.csv"))
  # write.csv(up.target.ranks,file=paste0("up.target.ranks.csv"))
  # write.csv(down.target.ranks,file=paste0("down.target.ranks.csv"))
  # 
  
  # current.subtitle <- paste("aREA ES = ",area.es.value," : aREA NES = ",area.nes.value," : aREA p-value = ",area.pvalue.num,"e",
  #                           area.pvalue.power,sep = "", " : DETOR score = ",
  #                           round((tcm.size-(tcm.size-(tcm.size-(length(up.target.ranks[up.target.ranks<max.up.es.rank])+ 
  #                                                                  length(down.target.ranks[down.target.ranks>max.down.es.rank])))))/(tcm.size-(tcm.size-(length(up.target.ranks[up.target.ranks<max.up.es.rank])+ 
  #                                                                                                                                                           length(down.target.ranks[down.target.ranks>max.down.es.rank])))),digits=2))
  # 
  
  
  
  x.max <- ceiling(max(plot.data$x.values)/x.step)*x.step
  
  # activeES <- max(plot.data$)
  
  current.subtitle <- paste("ES up:", round(max.up.es.value,1)," down:",round(max.down.es.value,1),  " aREA NES:",area.nes.value, " p-value:",area.pvalue.num,"e",area.pvalue.power)
  
  final.plot <- ggplot(plot.data,aes(x = x.values, y = y.values, colour = type)) + geom_rect(xmin = -10, xmax = x.max + 10, ymin = 1.005, ymax = 1.5, fill = "white", inherit.aes = FALSE) + geom_rect(xmin = -10, xmax = x.max + 10, ymin = -1.005, ymax = -1.5, fill = "white", inherit.aes = FALSE) + geom_segment(x = max.up.es.rank, xend = max.up.es.rank, y = 0, yend = max.up.es.value, colour = color.positive, lwd = .5) + geom_segment(x = max.down.es.rank, xend = max.down.es.rank, y = 0, yend = max.down.es.value, colour = color.negative, lwd = .5) + geom_segment(x = 0, y = 0, xend = max(plot.data$x.values), yend = 0, colour = "black", lwd = 1) + geom_line(lwd = 1) + scale_x_continuous(limits = c(0,x.max), breaks = seq(from = 0, to = x.max, by = x.step)) + theme_bw() + guides(colour = FALSE, alpha = FALSE, fill = FALSE) + labs(title = main.title, subtitle = current.subtitle, x = x.label, y = y.label) + theme(plot.title = element_text(hjust = 0.5, colour = "black", size = 20), plot.subtitle = element_text(hjust = 0.5, colour = "black", size = 18), axis.title.x = element_text(colour = "black", size = 18), axis.title.y = element_text(colour = "black", size = 18), axis.text.x = element_text(colour = "black", size = 15), axis.text.y = element_text(colour = "black", size = 15), panel.grid.major.x = element_line(colour = "grey", size = .5, linetype = "dotted"), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_line(colour = "grey", size = .5, linetype = "dotted"), panel.grid.minor.y = element_blank()) + geom_segment(data = plot.data.geneset, aes(x = x.values, xend = x.values, y = y.start, yend = y.end, colour = type), lwd = .75, inherit.aes = FALSE) + geom_segment(x = 0, y = 1.1, xend = max(plot.data$x.values), yend = 1.1, colour = "black", lwd = 1) + geom_segment(x = 0, y = 1.2, xend = max(plot.data$x.values), yend = 1.2, colour = "black", lwd = 1) + geom_segment(x = 0, y = 1.3, xend = max(plot.data$x.values), yend = 1.3, colour = "black", lwd = 1) + geom_segment(x = 0, y = 1.098, xend = 0, yend = 1.302, colour = "black", lwd = 1) + geom_segment(x = max(plot.data$x.values), y = 1.098, xend = max(plot.data$x.values), yend = 1.302, colour = "black", lwd = 1) + scale_colour_manual(values = c(color.positive,color.negative)) + scale_y_continuous(limits = c(-1.05,1.3), breaks = seq(from = -1, to = 1, by = y.step)) 
  
  return(final.plot)
  
}


# *************************************USAGE
# ref.signature <- 'CCl4'
# query.geneset <- 'BDL'
# a.p<-aREA_GSEA_Plot(ref.sig,test.sig,uniform.wts = FALSE,
#                     tcm.size,color.positive = "red", color.negative = "blue",
#                     main.title = paste0('GSEA Enrichment Analysis\n',ref.signature,' as reference signature and ', query.geneset,' as query gene set'),
#                     y.label = "Enrichment Score", x.label = paste0("Rank ordered regulators by activity in reference signature"),
#                     x.step = 500, y.step = 0.5,gsea.score = 1)

