# Load the required packages
library(biomaRt)
library(dplyr)

# Use `biomaRt::getBM()` to create the conversion table, with the three most
# common human ID types.
biomart_table_human_1 <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"),
  mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
)

# Replace empty values with NA
biomart_table_human_2 <- biomart_table_human_1 %>% replace(. == "", NA)

# Keep only one row for each Ensembl gene
biomart_id_mapping_human <- biomart_table_human_2 %>%
  rename("entrez_gene_id" = entrezgene_id) %>%
  arrange(ensembl_gene_id, hgnc_symbol, entrez_gene_id) %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

usethis::use_data(biomart_id_mapping_human, overwrite = TRUE)
