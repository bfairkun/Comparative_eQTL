library(tidyverse)

Fam.filepath <- "AddressReviews/plink/GTEX.v8.fam"
GTEXCovs.filepath <- "GTEX_renalysis/data/GTEx_Analysis_v8_eQTL_covariates/Heart_Left_Ventricle.v8.covariates.txt"
Pheno.filepath <- "../../output/CellProportionPhenotypesNormalizedForGWAS.tab"

Fam.out.filepath <- "AddressReviews/GWAS/GTEX.v8.fam"
Cov.out.filepath <- "AddressReviews/GWAS/GTEX.v8.cov"

Fam <- read.delim(Fam.filepath, header=F, sep=" ")
GTEXCovs <- read.table(GTEXCovs.filepath, sep="\t", header=T) %>%
    filter(!str_detect(id, "InferredCov") ) %>%
    drop_na() %>%
    rownames_to_column() %>%
    dplyr::select(-rowname) %>%
    column_to_rownames("id") %>% as.data.frame() %>% t() %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    separate("rowname", into=c("Fam", "IndID"), sep="\\.")
Pheno <- read.delim(Pheno.filepath) %>%
    separate(Ind, into=c("Fam", "IndID"), remove=F, sep="-")

Fam %>%
    left_join(Pheno, by=c("V2"="IndID")) %>%
    dplyr::select(V1,V2,V3,V4,V5,V6=PC1) %>%
    write_delim(Fam.out.filepath, col_names=F)

Fam %>%
    left_join(GTEXCovs, by=c("V2"="IndID")) %>%
    mutate(intercept=1) %>%
    dplyr::select(intercept, everything(), -Fam, -V1, -V2, -V3, -V4, -V5, -V6) %>%
    write_delim(Cov.out.filepath, col_names=F, delim='\t')
