library(tidyverse)

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
  "SYT1",
  "NALCN",
  "UNC79",
  "UNC80",
  "FAM155A",
  "FAM155B",
  "NALF1",
  "NALF2"
)

associations = read_csv("open_targets.csv")

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
) %>%
  mutate(complex_gene = str_replace_all(complex_gene, c("NALF1" = "FAM155A", "NALF2" = "FAM155B"))) %>%
  filter(human_gene_name %in% associations$Gene_name)

expression_correlations = read_csv("expression_datasets/combined_expression_correlations.csv.gz") %>%
  rename(human_gene_name = `Gene name`) %>%
  filter(source != "allen_smartseq", human_gene_name %in% associations$Gene_name)

standardised_expression = expression_correlations %>%
  select(-Complex_nTPM_scc) %>%
  pivot_longer(ends_with("scc"), names_to = "complex_gene", values_to = "scc") %>%
  drop_na() %>%
  mutate(complex_gene = str_replace_all(complex_gene, c("NALF1" = "FAM155A", "NALF2" = "FAM155B"))) %>%
  mutate(complex_gene = str_remove(complex_gene, "_scc")) %>%
  pivot_wider(names_from = source, values_from = scc) %>%
  group_by(complex_gene) %>%
  mutate(across(is.numeric, ~ (.x - mean(.x, na.rm = TRUE)) / sd(.x, na.rm = TRUE)))

standardised_coevolution = coevo_scores %>%
  drop_na() %>%
  pivot_wider(names_from = db, values_from = similarity_score) %>%
  group_by(complex_gene) %>%
  mutate(
    across(is.numeric, ~ .x - mean(.x, na.rm = TRUE)),
    across(is.numeric, ~ .x / sd(.x, na.rm = TRUE))
    )

combined_scores = left_join(standardised_coevolution, standardised_expression) %>%
  rowwise() %>%
  mutate(na_sum = sum(is.na(c(omadb, orthodb, eggnog, HPA_tissue, HPA_cells, cao_2020, allen_10x)))) %>%
  filter(na_sum <= 2) %>%
  select(-na_sum) %>%
  mutate(across(is.numeric, ~ replace_na(.x, 0)))

ggplot(combined_scores %>% pivot_longer(is.numeric, names_to = "source", values_to = "std_cor"), aes(std_cor)) +
  geom_histogram() +
  facet_wrap(vars(complex_gene, source), scales = "free") +
  theme_classic()

combined_scores_wide = combined_scores %>%
  pivot_longer(is.numeric, names_to = "source", values_to = "score") %>%
  unite("gene_source", c(complex_gene, source)) %>%
  pivot_wider(names_from = gene_source, values_from = score)

data_matrix = combined_scores_wide %>%
  ungroup() %>%
  select(where(is.numeric)) %>%
  as.data.frame()

rownames(data_matrix) = combined_scores_wide %>% pull(human_gene_name)

umap_settings = umap::umap.defaults
umap_settings$min_dist = 0.1
umap_settings$n_neighbors = 50
umap_settings$n_epochs = 500

combined_umap = data_matrix %>%
  umap::umap(config = umap_settings)

combined_umap_result = bind_cols(combined_umap$layout %>% as_tibble(), combined_scores_wide)

combined_umap_result = combined_umap_result %>%
  filter(human_gene_name %in% c("NALCN", "UNC79")) %>%
  select(V1, V2) %>%
  summarise(V3 = mean(V1), V4 = mean(V2)) %>%
  bind_cols(combined_umap_result) %>%
  rowwise() %>%
  mutate(
    euclid_distance = sqrt(sum((c(V1, V2) - c(V3, V4))^2))
  )

hit_list_points = combined_umap_result %>%
  filter(human_gene_name %in% hit_list) %>%
  select(V1, V2, human_gene_name) %>%
  unique()

ggplot(combined_umap_result, aes(x=V1, y = V2)) +
  geom_point(aes(fill = euclid_distance), shape = 21, colour = "white") +
  ggrepel::geom_text_repel(data = hit_list_points, aes(label = human_gene_name), max.overlaps = Inf, min.segment.length = 0) +
  theme_classic() +
  theme(legend.position = "left") +
  coord_fixed() +
  scale_fill_fermenter(palette = "RdYlBu", breaks = c(0:5)) +
  labs(
    fill = "Euclidean distance\nfrom NALCN/UNC79/UNC80",
    title = "UMAP of gene expression correlations"
  )
