library(tidyverse)

ncbi_ids = read_csv("coevolution_datasets/ncbi_ranked_lineage.csv.gz")

euk_species = read_tsv("coevolution_datasets/oma_species.tsv", skip = 3, col_names = c("oma_species_id", NA, "ncbi_species_id")) %>%
  select(oma_species_id, ncbi_species_id) %>%
  filter(ncbi_species_id %in% ncbi_ids$ncbi_species_id)

groups = read_tsv("coevolution_datasets/oma_groups.tsv", skip = 3, col_names = c("og_id", "oma_fingerprint")) %>%
  select(-oma_fingerprint) %>%
  pivot_longer(-og_id, names_to = NULL, values_to = "oma_gene_id") %>%
  drop_na() %>%
  mutate(oma_gene_id = str_replace_all(oma_gene_id, "\\t", "_")) %>%
  separate_longer_delim(oma_gene_id, delim = "_") %>%
  mutate(oma_species_id = str_extract(oma_gene_id, "[A-Z0-9]{5}")) %>%
  drop_na() %>%
  filter(oma_species_id %in% euk_species$oma_species_id)

human_crossref = read_csv("coevolution_datasets/omadb_human_crossref.csv.gz")

omadb_human_orthologs = human_crossref %>%
  select(-oma_gene_id) %>%
  left_join(groups) %>%
  unique() %>%
  left_join(euk_species) %>%
  select(-oma_species_id)

omadb_human_orthologs_euk = expand.grid(unique(euk_species$ncbi_species_id), unique(omadb_human_orthologs$og_id)) %>%
  as_tibble() %>%
  rename(ncbi_species_id = Var1, og_id = Var2) %>%
  left_join(omadb_human_orthologs %>% select(ncbi_species_id, og_id) %>% mutate(ortholog_present = TRUE) %>% unique()) %>%
  mutate(ortholog_present = replace_na(ortholog_present, FALSE)) %>%
  left_join(human_crossref %>% select(-oma_gene_id) %>% unique()) %>%
  select(-og_id)

write_csv(omadb_human_orthologs_euk, "coevolution_datasets/omadb_human_orthologs_euk.csv.gz")  

omadb_human_orthologs_euk %>%
  filter(human_gene_name %in% c("NALCN", "UNC79", "UNC80", "NALF1")) %>%
  select(ncbi_species_id, human_gene_name, ortholog_present) %>%
  pivot_wider(names_from = human_gene_name, values_from = ortholog_present) %>%
  left_join(omadb_human_orthologs_euk %>% select(ncbi_species_id, human_gene_name, ortholog_present)) %>%
  mutate(across(c(NALCN,UNC79,UNC80,NALF1), ~ !xor(.x, ortholog_present))) %>%
  group_by(human_gene_name) %>%
  summarise(across(c(NALCN,UNC79,UNC80,NALF1), ~ sum(.x)/length(.x))) %>%
  pivot_longer(-human_gene_name, names_to = "complex_gene", values_to = "omadb_similarity_score") %>%
  write_csv("coevolution_datasets/omadb_similarity_scores_euk.csv.gz")
