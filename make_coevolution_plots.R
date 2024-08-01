library(tidyverse)
set.seed(2024)

coevo_scores = bind_rows(
  read_csv("coevolution_datasets/omadb_similarity_scores_euk.csv.gz") %>%
    rename(similarity_score = omadb_similarity_score) %>%
    mutate(db = "omadb"),
  read_csv("coevolution_datasets/orthodb_similarity_scores_euk.csv.gz") %>%
    rename(similarity_score = orthodb_similarity_score) %>%
    mutate(db = "orthodb"),
  read_csv("coevolution_datasets/eggnog_similarity_scores_euk.csv.gz") %>%
    rename(similarity_score = eggnog_similarity_score) %>%
    mutate(db = "eggnog")
)

ncbi_ranked_lin = read_csv("coevolution_datasets/ncbi_ranked_lineage.csv.gz") %>%
  filter(superkingdom == "Eukaryota") %>%
  select(phylum, class, order, ncbi_species_id, scientific_name)

core_list = c(
  "NALCN",
  "UNC79",
  "UNC80",
  "FAM155A",
  "NALF1"
)

hit_list = c(
  "RIMS2",
  "CNTNAP2",
  "JPH2",
  "JPH4",
  "RIMBP2",
  "RPH3A",
  "ERC1",
  "RAB3A",
  "SNAP25",
  "KCTD16",
  "RIMS4",
  "TBC1D24",
  "VAMP2",
  "SYN2",
  "SPOCK2",
  "VAMP1",
  "CPLX1",
  "CPLX2",
  "STOML3",
  "STXBP1",
  "STX1A",
  "STX1B",
  "STOM",
  "SYT1"
)

snare_list = c(
  "STX1A",
  "STX1B",
  "STXBP1",
  "SNAP25",
  "VAMP2",
  "SYT1"
)

essential_list = read_csv("hart_essential_genes.csv") %>% pull(human_gene_name)
neuronal_list = read_tsv("goterm_annotation_nervous_system.tsv") %>% pull(SYMBOL) %>% unique()

coevo_barcodes = bind_rows(
  read_csv("coevolution_datasets/omadb_human_orthologs_euk.csv.gz") %>%
    filter(human_gene_name %in% c(core_list, hit_list, essential_list)) %>%
    mutate(db = "omadb"),
  read_csv("coevolution_datasets/orthodb_human_orthologs_euk.csv.gz") %>%
    filter(human_gene_name %in% c(core_list, hit_list, essential_list)) %>%
    mutate(db = "orthodb"),
  read_csv("coevolution_datasets/eggnog_human_orthologs_euk.csv.gz") %>%
    filter(human_gene_name %in% c(core_list, hit_list, essential_list)) %>%
    mutate(db = "eggnog")
) %>%
  left_join(coevo_scores %>% group_by(human_gene_name, db) %>% summarise(similarity_score = median(similarity_score))) %>%
  left_join(ncbi_ranked_lin)

coevo_barcodes %>%
  select(ncbi_species_id, db) %>%
  unique() %>%
  mutate(contains = TRUE) %>%
  pivot_wider(names_from = db, values_from = contains) %>%
  unite(dbs, c(omadb, orthodb, eggnog)) %>%
  count(dbs)

shared_species = coevo_barcodes %>%
  select(ncbi_species_id, db) %>%
  unique() %>%
  mutate(contains = TRUE) %>%
  pivot_wider(names_from = db, values_from = contains) %>%
  filter(omadb + orthodb + eggnog == 3) %>%
  pull(ncbi_species_id)

stat_mode <- function(x) {
  ux = unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

ggplot() +
  geom_tile(data = coevo_barcodes %>% filter(human_gene_name %in% core_list), aes(x = fct_reorder(factor(ncbi_species_id), ortholog_present), y = fct_reorder(human_gene_name, similarity_score), fill = ortholog_present)) +
  geom_tile(data = coevo_barcodes %>%
               filter(human_gene_name %in% essential_list) %>%
               group_by(ncbi_species_id, phylum, db) %>%
               summarise(mode_present = stat_mode(ortholog_present)),
            aes(x = factor(ncbi_species_id), fill = mode_present, y = "Average")) +
  facet_wrap(vars(db), scales = "free_x") +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggplot(coevo_scores %>% filter(human_gene_name %in% essential_list), aes(x = similarity_score, after_stat(density))) +
  geom_histogram(data = coevo_scores %>% filter(human_gene_name %in% neuronal_list), bins=21, fill = "white", colour = "black") +
  geom_freqpoly(bins = 31) +
  theme_classic() +
  scale_x_continuous(limits = c(0, 1)) +
  facet_wrap(vars(db), scales = "free_y")

ggplot(coevo_scores %>% filter(human_gene_name %in% essential_list), aes(x = similarity_score, after_stat(density))) +
  geom_histogram(data = coevo_scores %>% filter(human_gene_name %in% hit_list), bins=21, fill = "white", colour = "black") +
  geom_freqpoly(bins = 31) +
  theme_classic() +
  scale_x_continuous(limits = c(0, 1)) +
  facet_wrap(vars(db), scales = "free_y")

ggplot() +
  geom_tile(data = coevo_barcodes %>% filter(human_gene_name %in% core_list, ncbi_species_id %in% shared_species), aes(x = fct_reorder(scientific_name, ortholog_present), y = fct_reorder(human_gene_name, similarity_score), fill = ortholog_present)) +
  geom_tile(data = coevo_barcodes %>%
              filter(human_gene_name %in% essential_list, ncbi_species_id %in% shared_species) %>%
              group_by(scientific_name, phylum, db) %>%
              summarise(mode_present = stat_mode(ortholog_present)),
            aes(x = factor(scientific_name), fill = mode_present, y = "Average")) +
  facet_grid(rows = vars(db)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 5, hjust = 1, vjust = 0.3))
