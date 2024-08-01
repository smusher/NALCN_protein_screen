library(tidyverse)
library(dbplyr)

orthodb_con = DBI::dbConnect(RSQLite::SQLite(), dbname = "coevolution_datasets/orthodb_data")

species = tbl(orthodb_con, "species")
level2species = tbl(orthodb_con, "level2species")
OGs = tbl(orthodb_con, "OGs")
OG2genes = tbl(orthodb_con, "OG2genes")
OG_xrefs = tbl(orthodb_con, "OG_xrefs")
genes = tbl(orthodb_con, "genes")
gene_xrefs = tbl(orthodb_con, "gene_xrefs")

eukaryotic_ogs = OGs %>%
  filter(V2 == 2759) %>%
  select(orthodb_og_id = V1) %>%
  left_join(OG2genes %>% rename(orthodb_og_id = V1, orthodb_gene_id = V2)) %>%
  left_join(genes %>% select(orthodb_gene_id = V1, orthodb_species_id = V2, gene_name = V7)) %>%
  left_join(species %>% select(ncbi_species_id = V1, orthodb_species_id = V2, scientific_name = V3)) %>%
  collect()

human_genes = eukaryotic_ogs %>%
  filter(scientific_name == "Homo sapiens") %>%
  pull(orthodb_gene_id) %>%
  unique()

orthodb_human_crossref = gene_xrefs %>%
  filter(V1 %in% human_genes, V3 == "UniProt") %>%
  collect()

write_csv(orthodb_human_crossref %>% select(orthodb_gene_id = V1, uniprot_gene_id = V2), "coevolution_datasets/orthodb_human_crossref.csv.gz")
write_csv(eukaryotic_ogs, "coevolution_datasets/orthodb_eukaryotic_ogs.csv.gz")

orthodb_eukaryotic_species = level2species %>%
  filter(V1 == 2759) %>%
  select(orthodb_species_id = V2) %>%
  left_join(species %>% select(ncbi_species_id = V1, orthodb_species_id = V2, scientific_name = V3)) %>%
  collect() %>%
  write_csv("coevolution_datasets/orthodb_eukaryotic_species.csv.gz")

eukaryotic_ogs = read_csv("coevolution_datasets/orthodb_eukaryotic_ogs.csv.gz")
human_crossref = read_csv("coevolution_datasets/orthodb_human_crossref.csv.gz")
orthodb_eukaryotic_species = read_csv("coevolution_datasets/orthodb_eukaryotic_species.csv.gz")

human_orthologs = eukaryotic_ogs %>%
  filter(scientific_name == "Homo sapiens") %>%
  select(orthodb_og_id) %>%
  unique() %>%
  left_join(
    eukaryotic_ogs %>%
    select(orthodb_og_id, ncbi_species_id) %>%
    unique()
  ) %>%
  mutate(ortholog_present = TRUE)

orthodb_human_orthologs_euk = expand.grid(unique(orthodb_eukaryotic_species$ncbi_species_id), unique(human_orthologs$orthodb_og_id)) %>%
  as_tibble() %>%
  rename(ncbi_species_id = Var1, orthodb_og_id = Var2) %>%
  left_join(human_orthologs) %>%
  mutate(ortholog_present = replace_na(ortholog_present, FALSE)) %>%
  left_join(human_orthologs %>% select(orthodb_og_id, human_gene_name) %>% unique())

write_csv(orthodb_human_orthologs_euk, "coevolution_datasets/orthodb_human_orthologs_euk.csv.gz")

orthodb_human_orthologs_euk %>%
  filter(human_gene_name %in% c("NALCN", "UNC79", "UNC80", "NALF1")) %>%
  select(ncbi_species_id, human_gene_name, ortholog_present) %>%
  pivot_wider(names_from = human_gene_name, values_from = ortholog_present) %>%
  left_join(orthodb_human_orthologs_euk %>% select(ncbi_species_id, human_gene_name, ortholog_present)) %>%
  mutate(across(c(NALCN,UNC79,UNC80,NALF1), ~ !xor(.x, ortholog_present))) %>%
  group_by(human_gene_name) %>%
  summarise(across(c(NALCN,UNC79,UNC80,NALF1), ~ sum(.x)/length(.x))) %>%
  pivot_longer(-human_gene_name, names_to = "complex_gene", values_to = "orthodb_similarity_score") %>%
  write_csv("coevolution_datasets/orthodb_similarity_scores_euk.csv.gz")
