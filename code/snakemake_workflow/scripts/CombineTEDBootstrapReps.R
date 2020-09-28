library(tidyverse)

df <- list.files("AddressReviews/TED_Bootstraps", pattern = ".rds", full.names=T) %>%
  map(readRDS) %>% 
  bind_rows(.id="file_source")

df <- df %>% 
    mutate(
        replicate=as.numeric( replicate),
        file_source=as.numeric(file_source) -1) %>%
    mutate(FullReplicate = file_source*20 + replicate) %>%
    dplyr::select(-file_source, -replicate)

df.summary <- df %>%
    group_by(CellType, Species, gene) %>%
    summarize(mu.SE=sd(mu, na.rm=T),
        resid.SE=sd(resid, na.rm=T))

write_delim(df.summary, "../../output/CellTypeDispersion.SE.tsv.gz", delim='\t')
