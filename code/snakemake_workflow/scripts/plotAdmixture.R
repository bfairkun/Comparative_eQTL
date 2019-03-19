#!/usr/bin/env Rscript 

# Usage:
# plotAdmixture.R [data.K.Q] [data.fam] [output.pdf]

library(reshape2)
library(ggplot2)
library(plyr)

args = commandArgs(trailingOnly=TRUE)

Q.df <- read.table(args[1])
fam.df <- read.table(args[2], col.names = c('fam','ind','father','mother','sex','pheno'))


combined <- cbind(id=paste(fam.df$fam, fam.df$ind), Q.df)
toplot <- melt(combined)
ggplot(toplot, aes(x = id, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
    theme(text = element_text(size=6),
                  axis.text.x = element_text(angle=90, hjust=1))

combined <- cbind(fam.df[c(1,2)], Q.df)
combined$ind <- reorder(combined$ind, combined$V1)

# a hack for now to change fam to common subspecies names
combined$fam <- mapvalues(combined$fam, from=c("Pan_troglodytes_schweinfurthii", "Pan_troglodytes_ellioti", "Pan_troglodytes_ThisStudy", "Pan_troglodytes", "Pan_troglodytes_troglodytes", "Pan_troglodytes_verus"), to=c("Eastern", "Nigeria\nCameroon", "This Study", "Eastern", "Central", "Western"))

toplot <- melt(combined)
ggplot(toplot, aes(x = ind, y = value, fill = variable)) +
    geom_bar(stat = "identity") +
    scale_y_continuous(limits = c(-0.001,1.001), expand = c(0, 0)) +
    facet_grid(~fam, scales="free_x", space="free_x") +
    theme(legend.position="none") +
    theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_text(size=rel(0.5), angle=70, hjust=1))
ggsave(filename = args[3])

