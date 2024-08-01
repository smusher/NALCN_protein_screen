library(tidyverse)

# The eggnog_eukaryotic_ogs.tsv file can be downloaded from http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/2759/2759_members.tsv.gz

members = read_tsv("coevolution_datasets/eggnog_eukaryotic_ogs.tsv", col_names = c(NA, "eggnog_og_id", NA, NA, "contents"))
list = read_tsv("coevolution_datasets/eggnog_uniprot_mapping.tsv", col_names = c("uniprot_gene_id", "eggnog_og_id"))

human_crossref = read_csv("coevolution_datasets/eggnog_human_crossref.csv.gz")

mapping = human_crossref %>%
  select(-uniprot_gene_id) %>%
  left_join(list) %>%
  separate_longer_delim(eggnog_og_id, delim = ",") %>%
  mutate(eggnog_og_id = str_remove(eggnog_og_id, "ENOG50"))

eggnog_human_orthologs = left_join(
  mapping %>%
    select(human_gene_name, eggnog_og_id) %>%
    unique(),
  members %>%
    select(eggnog_og_id, contents)
  ) %>%
  separate_longer_delim(contents, delim = ",") %>%
  mutate(contents = str_remove(contents, "\\..*") %>% as.integer()) %>%
  rename(ncbi_species_id = contents) %>%
  unique()
  
eggnog_human_orthologs_euk = expand.grid(unique(eggnog_human_orthologs$ncbi_species_id), unique(eggnog_human_orthologs$eggnog_og_id)) %>%
  as_tibble() %>%
  rename(ncbi_species_id = Var1, eggnog_og_id = Var2) %>%
  left_join(eggnog_human_orthologs %>% select(-human_gene_name) %>% mutate(ortholog_present = TRUE) %>% unique()) %>%
  mutate(ortholog_present = replace_na(ortholog_present, FALSE)) %>%
  left_join(eggnog_human_orthologs %>% select(eggnog_og_id, human_gene_name) %>% unique()) %>%
  select(-eggnog_og_id)

write_csv(eggnog_human_orthologs_euk, "coevolution_datasets/eggnog_human_orthologs_euk.csv.gz")  

eggnog_human_orthologs_euk %>%
  filter(human_gene_name %in% c("NALCN", "UNC79", "UNC80", "NALF1")) %>%
  select(ncbi_species_id, human_gene_name, ortholog_present) %>%
  pivot_wider(names_from = human_gene_name, values_from = ortholog_present) %>%
  left_join(eggnog_human_orthologs_euk %>% select(ncbi_species_id, human_gene_name, ortholog_present)) %>%
  mutate(across(c(NALCN,UNC79,UNC80,NALF1), ~ !xor(.x, ortholog_present))) %>%
  group_by(human_gene_name) %>%
  summarise(across(c(NALCN,UNC79,UNC80,NALF1), ~ sum(.x)/length(.x))) %>%
  pivot_longer(-human_gene_name, names_to = "complex_gene", values_to = "eggnog_similarity_score") %>%
  write_csv("coevolution_datasets/eggnog_similarity_scores_euk.csv.gz")
