library(tidyverse)
library(igraph)
library(reshape2)

tmp <- data.frame(x=gl(2,3, labels=letters[24:25]),
                  y=gl(3,1,6, labels=letters[1:3]), 
                  z=c(1,2,3,3,3,2))

spread(tmp, y, z)
KingIn <- read.table('/home/bjf79/Chimp_eQTL_repo/PopulationSubstructure/plink/king.kin', sep='\t', header=TRUE)
ToPlot <- select(KingIn, ID1, ID2, Kinship) %>%
  graph.data.frame(directed=FALSE) %>%
  as_adjacency_matrix(names=TRUE,sparse=FALSE,attr="Kinship",type='both') %>%
  melt()
ToPlot$value[which(ToPlot$Var1 == ToPlot$Var2)] <- 0.5

ggplot(ToPlot, aes(x = Var1, y = Var2, fill = value)) +
  scale_fill_gradient2(low="navy", mid="white", high="red") +
  geom_raster() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))