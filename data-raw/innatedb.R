# Load packages
library(tidyverse)


# Load and filter InnateDB data
innatedb_init <- read_tsv(
  file = "https://innatedb.com/download/interactions/all.mitab.gz",
  col_types = cols()
)

innatedb_filtered <- innatedb_init %>%
  filter(
    # We just want human interactions
    ncbi_taxid_A == "taxid:9606(Human)" & ncbi_taxid_B == "taxid:9606(Human)",

    # Only interested in protein-protein interactions
    interactor_type_A == 'psi-mi:"MI:0326"(protein)' &
      interactor_type_B == 'psi-mi:"MI:0326"(protein)'
  )

innatedb_trimmed <- innatedb_filtered %>%
  dplyr::select(
    "ensembl_gene_A" = alt_identifier_A,
    "ensembl_gene_B" = alt_identifier_B
  ) %>%
  mutate(across(everything(), str_remove, pattern = "ensembl\\:"))


# Get gene mapping from biomaRt
biomart_mapping <- biomaRt::getBM(
  mart       = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl"),
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  verbose    = TRUE
)

innatedb <- innatedb_trimmed %>%
  left_join(biomart_mapping, by = c("ensembl_gene_A" = "ensembl_gene_id")) %>%
  rename("hgnc_symbol_A" = hgnc_symbol) %>%
  left_join(biomart_mapping, by = c("ensembl_gene_B" = "ensembl_gene_id")) %>%
  rename("hgnc_symbol_B" = hgnc_symbol) %>%
  relocate(ends_with("A"))

usethis::use_data(innatedb, overwrite = TRUE)
