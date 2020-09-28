library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
initial.seed <- as.numeric(args[1])
chunk.size <- as.numeric(args[2])

# initial.seed <- 0
# chunk.size <- 2

Normalized.Expression.Per.CellType <- readRDS(file = "MiscOutput/NormalizedExpressionPerCellType.rds")

Normalized.Expression.Per.CellType.chimp <- Normalized.Expression.Per.CellType %>%
  dplyr::select(-TotalBulkCount) %>% filter(Species=="Chimp")
N.chimp <- Normalized.Expression.Per.CellType.chimp$Ind %>% unique() %>% length()

Normalized.Expression.Per.CellType.human <- Normalized.Expression.Per.CellType %>%
  dplyr::select(-TotalBulkCount) %>% filter(Species=="Human")
N.human <- Normalized.Expression.Per.CellType.human$Ind %>% unique() %>% length()


sample_n_groups = function(tbl, size, replace = FALSE, weight = NULL) {
  # regroup when done
  grps = tbl %>% groups %>% lapply(as.character) %>% unlist
  # check length of groups non-zero
  keep = tbl %>% summarise() %>% ungroup() %>% sample_n(size, replace, weight)
  # keep only selected groups, regroup because joins change count.
  # regrouping may be unnecessary but joins do something funky to grouping variable
  tbl %>% right_join(keep, by=grps) %>% group_by_(.dots = grps)
}

N.Bootstrap.Reps <- chunk.size

datalist.chimp = list()
datalist.human = list()

set.seed(initial.seed)
ptm <- proc.time()
for (i in 1:N.Bootstrap.Reps) {
  print(i)
  datalist.chimp[[i+initial.seed]] <- Normalized.Expression.Per.CellType.chimp %>%
  group_by(Ind) %>%
  sample_n_groups(N.chimp, replace = T) %>%
  ungroup() %>%
  group_by(gene, CellType) %>%
  summarise(mu=mean(Log.CPM.Expression), log.var=log(var(Log.CPM.Expression))) %>%
  group_by(CellType) %>%
  do(data.frame(., resid = residuals(loess(log.var ~ mu, data=., degree=1, na.action="na.exclude")))) %>%
  dplyr::select(-log.var)
  datalist.human[[i+initial.seed]] <- Normalized.Expression.Per.CellType.chimp %>%
  group_by(Ind) %>%
  sample_n_groups(N.human, replace = T) %>%
  ungroup() %>%
  group_by(gene, CellType) %>%
  summarise(mu=mean(Log.CPM.Expression), log.var=log(var(Log.CPM.Expression))) %>%
  group_by(CellType) %>%
  do(data.frame(., resid = residuals(loess(log.var ~ mu, data=., degree=1, na.action="na.exclude")))) %>%
  dplyr::select(-log.var)
}
proc.time() - ptm

dplyr::bind_rows(
    dplyr::bind_rows(datalist.chimp, .id='replicate') %>%
    mutate(Species="Chimp"),
    dplyr::bind_rows(datalist.human, .id='replicate') %>%
    mutate(Species="Human")) %>%
    saveRDS(paste0("AddressReviews/TED_Bootstraps/",initial.seed,".rds" ))
