library(tidyverse)
library(scales)
library(patchwork)

id_mapping = read_csv("expression_datasets/id_mapping.csv")

hpa_tissue = read_tsv("expression_datasets/hpa_tissue_consensus.tsv")

hpa_cells = read_tsv("expression_datasets/hpa_single_cell_type.tsv")

cao_2020 = read_csv("expression_datasets/cao_2020_single_cell_type.csv") %>%
  separate(RowID, into = c("Ensembl_id", NA), sep = "\\.") %>%
  left_join(id_mapping)

allen_smartseq = read_csv("expression_datasets/allen_brain_smartseq_medians.csv") %>%
  rename(`Mouse gene` = "feature") %>%
  left_join(id_mapping)

allen_10x = read_csv("expression_datasets/allen_brain_10x_medians.csv") %>%
  rename(`Mouse gene` = "feature") %>%
  left_join(id_mapping)

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

essential_list = read_csv("hart_essential_genes.csv") %>% pull(human_gene_name)
neuronal_list = read_tsv("goterm_annotation_nervous_system.tsv") %>% pull(SYMBOL) %>% unique()

hpa_tissue %>%
  select(-Gene) %>%
  filter(`Gene name` %in% c("NALCN", "UNC79", "UNC80", "FAM155A", "FAM155B")) %>%
  pivot_wider(names_from = `Gene name`, values_from = nTPM) %>%
  group_by(Tissue) %>%
  mutate(Complex_nTPM = sum(NALCN, UNC79, UNC80, FAM155A, FAM155B)) %>%
  left_join(hpa_tissue %>% select(-Gene)) -> tissue_df

tissue_df %>%
  group_by(`Gene name`) %>%
  summarise(across(c(NALCN,UNC79,UNC80,FAM155A,FAM155B,Complex_nTPM), ~ cor(.x, nTPM, method = "spearman", use = "complete.obs"), .names = "{.col}_scc")) %>%
  left_join(tissue_df) -> tissue_df

tissue_df %>%
  filter(!(`Gene name` %in% c("NALCN", "UNC79", "UNC80", "FAM155A", "FAM155B")), `Gene name` %in% hit_list) %>%
  ggplot(aes(x = Complex_nTPM, y = nTPM, fill = UNC79_scc)) +
  geom_point(shape = 21, size = 3) +
  facet_wrap(vars(`Gene name`), scales = "free") +
  theme_classic() +
  scale_x_log10(labels = label_number()) +
  scale_y_log10(labels = label_number()) +
  scale_fill_distiller(palette = "RdYlBu")

hpa_cells %>%
  select(-Gene) %>%
  filter(`Gene name` %in% c("NALCN", "UNC79", "UNC80", "FAM155A", "FAM155B")) %>%
  pivot_wider(names_from = `Gene name`, values_from = nTPM) %>%
  group_by(`Cell type`) %>%
  mutate(Complex_nTPM = sum(NALCN, UNC79, UNC80, FAM155A, FAM155B)) %>%
  left_join(hpa_cells %>% select(-Gene)) -> cells_df

cells_df %>%
  group_by(`Gene name`) %>%
  summarise(across(c(NALCN,UNC79,UNC80,FAM155A,FAM155B,Complex_nTPM), ~ cor(.x, nTPM, method = "spearman", use = "complete.obs"), .names = "{.col}_scc")) %>%
  left_join(cells_df) -> cells_df

cells_df %>%
  filter(!(`Gene name`%in% c("NALCN", "UNC79", "UNC80", "FAM155A", "FAM155B")), `Gene name` %in% hit_list) %>%
  ggplot(aes(x = Complex_nTPM, y = nTPM, fill = UNC79_scc)) +
  geom_point(shape = 21, size = 3) +
  facet_wrap(vars(`Gene name`), scales = "free") +
  theme_classic() +
  scale_x_log10(labels = label_number()) +
  scale_y_log10(labels = label_number()) +
  scale_fill_distiller(palette = "RdYlBu")

cao_2020 %>%
  select(-Ensembl_id, -`Mouse gene`) %>%
  filter(`Gene name` %in% c("NALCN", "UNC79", "UNC80", "NALF1", "NALF2")) %>%
  pivot_longer(-`Gene name`, names_to = "Cell type", values_to = "nTPM") %>%
  pivot_wider(names_from = `Gene name`, values_from = nTPM) %>%
  group_by(`Cell type`) %>%
  mutate(Complex_nTPM = sum(NALCN, UNC79, UNC80, NALF1, NALF2)) %>%
  left_join(cao_2020 %>%
              select(-Ensembl_id,-`Mouse gene`) %>%
              pivot_longer(-`Gene name`, names_to = "Cell type", values_to = "nTPM")
            ) %>%
  rename(FAM155A = NALF1, FAM155B = NALF2) -> cao_df

cao_df %>%
  group_by(`Gene name`) %>%
  summarise(across(c(NALCN,UNC79,UNC80,FAM155A,FAM155B,Complex_nTPM), ~ cor(.x, nTPM, method = "spearman", use = "complete.obs"), .names = "{.col}_scc")) %>%
  left_join(cao_df) -> cao_df

cao_df %>%
  filter(!(`Gene name`%in% c("NALCN", "UNC79", "UNC80", "FAM155A", "FAM155B")), `Gene name` %in% hit_list) %>%
  ggplot(aes(x = UNC80, y = nTPM, fill = UNC80_scc)) +
  geom_point(shape = 21, size = 3) +
  facet_wrap(vars(`Gene name`), scales = "free") +
  theme_classic() +
  scale_x_log10(labels = label_number()) +
  scale_y_log10(labels = label_number()) +
  scale_fill_distiller(palette = "RdYlBu")

allen_smartseq %>%
  select(-Ensembl_id, -`Mouse gene`) %>%
  filter(`Gene name` %in% c("NALCN", "UNC79", "UNC80", "FAM155A")) %>%
  pivot_longer(-`Gene name`, names_to = "Cell type", values_to = "nTPM") %>%
  pivot_wider(names_from = `Gene name`, values_from = nTPM) %>%
  group_by(`Cell type`) %>%
  mutate(Complex_nTPM = sum((2^NALCN)-1, (2^UNC79)-1, (2^UNC80)-1, (2^FAM155A)-1)) %>%
  left_join(allen_smartseq %>%
              select(-Ensembl_id,-`Mouse gene`) %>%
              pivot_longer(-`Gene name`, names_to = "Cell type", values_to = "nTPM") %>%
              mutate(nTPM = 2^nTPM)
  ) -> allen_smartseq_df

allen_smartseq_df %>%
  group_by(`Gene name`) %>%
  summarise(across(c(NALCN,UNC79,UNC80,FAM155A,Complex_nTPM), ~ cor(.x, nTPM, method = "spearman", use = "complete.obs"), .names = "{.col}_scc")) %>%
  left_join(allen_smartseq_df) -> allen_smartseq_df

allen_smartseq_df %>%
  filter(!(`Gene name`%in% c("NALCN", "UNC79", "UNC80", "FAM155A")), `Gene name` %in% hit_list) %>%
  ggplot(aes(x = UNC80, y = nTPM, fill = UNC80_scc)) +
  geom_point(shape = 21, size = 3) +
  facet_wrap(vars(`Gene name`), scales = "free") +
  theme_classic() +
  scale_x_log10(labels = label_number()) +
  scale_y_log10(labels = label_number()) +
  scale_fill_distiller(palette = "RdYlBu")

allen_10x %>%
  select(-Ensembl_id, -`Mouse gene`) %>%
  filter(`Gene name` %in% c("NALCN", "UNC79", "UNC80", "FAM155A")) %>%
  pivot_longer(-`Gene name`, names_to = "Cell type", values_to = "nTPM") %>%
  pivot_wider(names_from = `Gene name`, values_from = nTPM) %>%
  group_by(`Cell type`) %>%
  mutate(Complex_nTPM = sum((2^NALCN)-1, (2^UNC79)-1, (2^UNC80)-1, (2^FAM155A)-1)) %>%
  left_join(allen_10x %>%
              select(-Ensembl_id,-`Mouse gene`) %>%
              pivot_longer(-`Gene name`, names_to = "Cell type", values_to = "nTPM")%>%
              mutate(nTPM = (2^nTPM)-1)
  ) -> allen_10x_df

allen_10x_df %>%
  group_by(`Gene name`) %>%
  summarise(across(c(NALCN,UNC79,UNC80,FAM155A,Complex_nTPM), ~ cor(.x, nTPM, method = "spearman", use = "complete.obs"), .names = "{.col}_scc")) %>%
  left_join(allen_10x_df) -> allen_10x_df

allen_10x_df %>%
  filter(!(`Gene name`%in% c("NALCN", "UNC79", "UNC80", "FAM155A")), `Gene name` %in% hit_list) %>%
  ggplot(aes(x = UNC80, y = nTPM, fill = UNC80_scc)) +
  geom_point(shape = 21, size = 3) +
  facet_wrap(vars(`Gene name`), scales = "free") +
  theme_classic() +
  scale_x_log10(labels = label_number()) +
  scale_y_log10(labels = label_number()) +
  scale_fill_distiller(palette = "RdYlBu")

combined_expression = bind_rows(
  tissue_df %>% mutate(source = "HPA_tissue"),
  cells_df %>% mutate(source = "HPA_cells"),
  cao_df %>% mutate(source = "cao_2020"),
  allen_smartseq_df %>% mutate(source = "allen_smartseq"),
  allen_10x_df %>% mutate(source = "allen_10x")
) %>%
  select(`Gene name`, ends_with("scc"), source) %>%
  unique()

#write_csv(combined_expression, "expression_datasets/combined_expression_correlations.csv.gz")
combined_expression = read_csv("expression_datasets/combined_expression_correlations.csv.gz")

combined_expression %>%
  group_by(`Gene name`) %>%
  summarise(correlation = mean(Complex_nTPM_scc, na.rm = TRUE)) %>%
  drop_na() %>%
  ggplot(aes(y = fct_reorder(`Gene name`, correlation), x = correlation)) +
  geom_point() +
  geom_point(data=combined_expression %>%
               filter(`Gene name` %in% hit_list) %>%
               group_by(`Gene name`) %>%
               summarise(correlation = mean(Complex_nTPM_scc, na.rm = TRUE)), colour = "red") +
  ggrepel::geom_label_repel(data=combined_expression %>%
               filter(`Gene name` %in% hit_list) %>%
               group_by(`Gene name`) %>%
               summarise(correlation = mean(Complex_nTPM_scc, na.rm = TRUE)), aes(label = `Gene name`), max.overlaps = Inf) +
  theme_classic() +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
  coord_cartesian(clip = "off")

bind_rows(
  tissue_df %>% mutate(source = "HPA_tissue"),
  cells_df %>% mutate(source = "HPA_cells"),
  cao_df %>% mutate(source = "cao_2020"),
  allen_10x_df %>% mutate(source = "allen_10x") %>% filter(UNC80 >= 5)
) %>%
  filter(`Gene name` %in% c("SNAP25", "UBA1"))  %>%
  group_by(source, `Gene name`) %>%
  mutate(across(c(UNC80, nTPM), ~ rank(.x))) -> temp_df

ggplot(temp_df, aes(x = UNC80, y = nTPM)) +
  geom_smooth(method = "lm") +
  geom_point(shape = 21, size = 3, fill = "white") +
  geom_point(data = temp_df %>% filter(Tissue == "skin"), colour = "red") +
  geom_point(data = temp_df %>% filter(Tissue == "hippocampal formation"), colour = "blue") +
  geom_point(data = temp_df %>% filter(`Cell type` == "Adipocytes"), colour = "red") +
  geom_point(data = temp_df %>% filter(`Cell type` == "Inhibitory neurons"), colour = "blue") +
  geom_point(data = temp_df %>% filter(`Cell type` == "Lung-Lymphoid cells"), colour = "red") +
  geom_point(data = temp_df %>% filter(`Cell type` == "Cerebellum-Purkinje neurons"), colour = "blue") +
  geom_point(data = temp_df %>% filter(`Cell type` == "14_Lamp5"), colour = "red") +
  geom_point(data = temp_df %>% filter(`Cell type` == "293_L6 CT CTX"), colour = "blue") +
  facet_wrap(vars(source,`Gene name`), scales = "free") +
  theme_classic() +
  coord_cartesian(clip = "off")

temp_df %>%
  split(~temp_df$`Gene name` + temp_df$source) %>%
  map(\(df) lm(nTPM ~ UNC80, data = df)) %>%
  map(summary) %>%
  map_dbl("r.squared")

bind_rows(
  tissue_df %>% mutate(source = "HPA_tissue"),
  cells_df %>% mutate(source = "HPA_cells"),
  cao_df %>% mutate(source = "cao_2020"),
  allen_10x_df %>% mutate(source = "allen_10x")
) %>%
  filter(`Gene name` %in% c("NALCN", "UNC79", "UNC80", "FAM155A"))  %>%
  select(gene_a = `Gene name`, Tissue, NALCN, UNC79, UNC80, FAM155A, ntpm_a = nTPM, source) %>%
  pivot_longer(NALCN:FAM155A, names_to = "gene_b", values_to = "ntpm_b") %>%
  group_by(source, gene_a, gene_b) %>%
  mutate(
    scc = cor(ntpm_a, ntpm_b, method = "spearman", use = "complete.obs"),
    across(c(ntpm_a, ntpm_b), ~ rank(.x))
    ) -> core_complex_df

ggplot(core_complex_df %>% filter(source == "HPA_tissue"), aes(x = ntpm_a, y = ntpm_b)) +
  geom_smooth(method = "lm") +
  geom_point(shape = 21, size = 3, aes(fill = scc)) +
  facet_grid(rows = vars(gene_a), cols = vars(gene_b)) +
  theme_classic() +
  coord_cartesian(clip = "off") +
  scale_fill_fermenter(palette = "RdYlBu", limits = c(-1, 1), breaks = seq(from=-1, to=1, length.out=9)) -> plot_tissue

ggplot(core_complex_df %>% filter(source == "HPA_cells"), aes(x = ntpm_a, y = ntpm_b)) +
  geom_smooth(method = "lm") +
  geom_point(shape = 21, size = 3, aes(fill = scc)) +
  facet_grid(rows = vars(gene_a), cols = vars(gene_b)) +
  theme_classic() +
  coord_cartesian(clip = "off") +
  scale_fill_fermenter(palette = "RdYlBu", limits = c(-1, 1), breaks = seq(from=-1, to=1, length.out=9)) -> plot_cells

ggplot(core_complex_df %>% filter(source == "cao_2020"), aes(x = ntpm_a, y = ntpm_b)) +
  geom_smooth(method = "lm") +
  geom_point(shape = 21, size = 3, aes(fill = scc)) +
  facet_grid(rows = vars(gene_a), cols = vars(gene_b)) +
  theme_classic() +
  coord_cartesian(clip = "off") +
  scale_fill_fermenter(palette = "RdYlBu", limits = c(-1, 1), breaks = seq(from=-1, to=1, length.out=9)) -> plot_cao

ggplot(core_complex_df %>% filter(source == "allen_10x"), aes(x = ntpm_a, y = ntpm_b)) +
  geom_smooth(method = "lm") +
  geom_point(shape = 21, size = 3, aes(fill = scc)) +
  facet_grid(rows = vars(gene_a), cols = vars(gene_b)) +
  theme_classic() +
  coord_cartesian(clip = "off") +
  scale_fill_fermenter(palette = "RdYlBu", limits = c(-1, 1), breaks = seq(from=-1, to=1, length.out=9)) -> plot_allen

plot_tissue + plot_cells + plot_cao + plot_allen

control_compare = combined_expression %>%
  select(-Complex_nTPM_scc, -FAM155B_scc) %>%
  pivot_longer(ends_with("scc"), names_to = "core_complex_gene", values_to = "scc")

ggplot(control_compare, aes(x = scc, after_stat(density))) +
  geom_freqpoly(bins = 31) +
  geom_histogram(data = control_compare %>% filter(`Gene name` %in% hit_list), bins=21) +
  theme_classic() +
  facet_wrap(vars(source))

neuronal_compare = combined_expression %>%
  select(-Complex_nTPM_scc, -FAM155B_scc) %>%
  pivot_longer(ends_with("scc"), names_to = "core_complex_gene", values_to = "scc") %>%
  filter(`Gene name` %in% c(essential_list, neuronal_list), source != "allen_smartseq")

ggplot(neuronal_compare %>% filter(`Gene name` %in% essential_list), aes(x = scc, after_stat(density))) +
  geom_freqpoly(bins = 31) +
  geom_histogram(data = neuronal_compare %>% filter(`Gene name` %in% neuronal_list), bins=21) +
  theme_classic() +
  facet_wrap(vars(source), nrow = 1) +
  scale_x_continuous(limits = c(-1, 1))

neuronal_compare %>%
  drop_na() %>%
  filter(source != "allen_10x") %>%
  mutate(list = case_when(`Gene name` %in% essential_list ~ "essential", `Gene name` %in% neuronal_list ~ "neuronal")) %>%
  group_by(list) %>%
  summarise(median_scc = median(scc))

