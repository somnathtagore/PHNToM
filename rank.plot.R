library(ggplot2)
theme_set(theme_bw())

cty_mpg<-read.csv("~/Documents/cty_mpg.csv",header=TRUE, sep=",")
rownames(cty_mpg) <- cty_mpg[,1]
cty_mpg<-cty_mpg[,-1]

cty_mpg<-as.data.frame(cty_mpg)
cty_mpg
cty_mpg <- cty_mpg[order(cty_mpg$Samples), ]  # sort
cty_mpg$Proteins <- factor(cty_mpg$Proteins, levels = cty_mpg$Proteins)  # to retain the order in plot.
head(cty_mpg, 4)
# Plot
ggplot(cty_mpg, aes(x=Proteins, y=Samples)) + 
  geom_point(size=3) + 
  geom_segment(aes(x=Proteins, 
                   xend=Proteins, 
                   y=0, 
                   yend=Samples)) + 
  geom_label_repel(aes(label = Samples)) +
  labs(title="Mutation Rank Plot: PIK3CA-WT (BRCA)", 
       subtitle="Mutation Mimicry: Like PIK3CA-MUT LOF") + 
  theme(axis.text.x = element_text(angle=90, vjust=0.6))
ggsave("mimic_PIK3CA_brca_wt_LOF.pdf",width=8,height = 8)


# ########
# 
# library(ggplot2) 
# install.packages("ggplotify")
# library(treemapify)
# library(ggplotify)
# proglangs <- read.csv("https://raw.githubusercontent.com/selva86/datasets/master/proglanguages.csv")
# 
# # plot
# treeMapCoordinates <- treemapify(proglangs,
#                                  area = "value",
#                                  fill = "parent",
#                                  label = "id",
#                                  group = "parent")
# 
# treeMapPlot <- ggplotify(treeMapCoordinates) + 
#   scale_x_continuous(expand = c(0, 0)) +
#   scale_y_continuous(expand = c(0, 0)) +
#   scale_fill_brewer(palette = "Dark2")
# 
# print(treeMapPlot)
