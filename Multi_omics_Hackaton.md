---
title: "Multi-Omics Analysis of Adenomas"
author: "Deborah Mina Ikann"
date: "`r Sys.Date()`"

---

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE,
  warning = FALSE
)
```

```{r load-libraries}
# Data manipulation
library(readr)
library(data.table)
library(DT)

# Differential analysis
library(limma)
library(globaltest)

# Visualization
library(ggplot2)
library(qqman)
library(qqplotr)

# Multi-omics
library(mixOmics)
library(glmnet)

# Microbiome analysis
library(dplyr)
library(tidyr)
library(knitr)
library(kableExtra)
library(glmnet)
library(caret)
```

# 1. Data Import
```{r import-data}
# Read metabolomics
mtb <- read_tsv("mtb.tsv", show_col_types = FALSE)

# Save the first column as sample ID vector and set as rownames
first_col <- mtb[[1]]
mtb <- as.data.frame(mtb)
rownames(mtb) <- first_col
mtb <- mtb[ , -1, drop = FALSE]

# Read metadata
metadata <- read_tsv(
  "metadata.tsv",
  show_col_types = FALSE
)

# Read metagenomics (genera counts)
Meg <- read_tsv("genera.counts.tsv", show_col_types = FALSE)
Meg <- as.data.frame(Meg)
rownames(Meg) <- Meg[[1]]
Meg <- Meg[ , -1, drop = FALSE]

## === New: rename metagenomic columns to Family;Genus ===

tax_full <- colnames(Meg)

get_family_genus <- function(x) {
  parts <- strsplit(x, ";")[[1]]
  # take last two levels (family and genus)
  last_two <- tail(parts, 2)
  # remove rank prefixes like "f__" and "g__"
  last_two <- sub("^[a-z]__", "", last_two)
  # combine as "Family;Genus" (or use "_" instead of ";" if you prefer)
  paste(last_two, collapse = ";")
}

tax_fg <- vapply(tax_full, get_family_genus, character(1))

# Mapping table for later reference
tax_map_meg <- data.frame(
  taxon_short = tax_fg,
  taxon_full  = tax_full,
  stringsAsFactors = FALSE
)

# Apply short names to Meg
colnames(Meg) <- tax_fg
```

```{r verify-sample-ids}
cat("Metabolomics dimensions (samples x features):", dim(mtb), "\n")
cat("Metagenomics dimensions (samples x taxa):", dim(Meg), "\n")
cat("Metadata dimensions (samples x variables):", dim(metadata), "\n")

cat("Total missing values in metabolomics:", sum(is.na(mtb)), "\n")
cat("Total missing values in metagenomics:", sum(is.na(Meg)), "\n")

summary(metadata)
```

# 2. Data Preprocessing

```{r align-samples}
# Create aligned datasets and metadata (used throughout analysis)
align_data <- function(omics_data, metadata) {
  # Match samples
  aligned_data <- omics_data[rownames(omics_data) %in% metadata$Sample, , drop = FALSE]
  
  # Reorder metadata to match
  idx <- match(rownames(aligned_data), metadata$Sample)
  metadata_aligned <- metadata[idx, ]
  
  # Verify alignment
  stopifnot(all(metadata_aligned$Sample == rownames(aligned_data)))
  
  return(list(data = aligned_data, metadata = metadata_aligned))
}

# Align metabolomics data
mtb_aligned <- align_data(mtb, metadata)
mtb_DA <- mtb_aligned$data
metadata_aligned <- mtb_aligned$metadata

# Extract group factors
Group <- factor(metadata_aligned$Study.Group)
```

```{r}
# Extract sample IDs
ids      <- metadata$Sample
ids_mtb  <- rownames(mtb)
ids_meg  <- rownames(Meg)

# Check for mismatches
cat("Metadata IDs missing in mtb:", setdiff(ids, ids_mtb), "\n")
cat("Metadata IDs missing in Meg:", setdiff(ids, ids_meg), "\n")
cat("All IDs from metadata present in mtb:", all(ids %in% ids_mtb), "\n")
cat("All IDs from metadata present in Meg:", all(ids %in% ids_meg), "\n")

cat("Total number of genera in Meg:", ncol(Meg), "\n")
cat("Total number of metabolites:", ncol(mtb), "\n")
```
```{r}
align_data <- function(omics_data, metadata) {
  # Keep only samples present in metadata
  aligned_data <- omics_data[rownames(omics_data) %in% metadata$Sample, , drop = FALSE]
  
  # Reorder metadata to match omics_data row order
  idx <- match(rownames(aligned_data), metadata$Sample)
  metadata_aligned <- metadata[idx, ]
  
  # Verify alignment
  stopifnot(all(metadata_aligned$Sample == rownames(aligned_data)))
  
  return(list(data = aligned_data, metadata = metadata_aligned))
}
```


# 3. Metabolomics Analysis

## 3.1 Quality Control

```{r mtb-qc}
# Align metabolomics with metadata
mtb_aligned <- align_data(mtb, metadata)
mtb_DA <- mtb_aligned$data
metadata_aligned <- mtb_aligned$metadata

# Group factor
Group <- factor(metadata_aligned$Study.Group,
                levels = c("Control", "Adenoma", "Carcinoma"))

```

```{r}
# Boxplots
boxplot(t(mtb_DA), las = 2, col = "lightgreen", 
        main = "Metabolomics Boxplots (log2 scale)", cex.axis = 0.6)

# Variability distribution
sds <- sapply(mtb_DA, sd, na.rm = TRUE)
hist(sds, main = "Standard Deviation of Metabolites", col = "skyblue")

# Run order
metadata_aligned$RunOrder <- seq_len(nrow(metadata_aligned))

# Total and median intensity per sample
metadata_aligned$total_intensity  <- rowSums(mtb_DA, na.rm = TRUE)
metadata_aligned$median_intensity <- apply(mtb_DA, 1, median, na.rm = TRUE)

ggplot(metadata_aligned,
       aes(x = median_intensity)) +
  geom_histogram(bins = 30, colour = "black", fill = "grey80") +
  theme_bw() +
  labs(x = "Median log2 intensity per sample",
       y = "Count of samples")

```


## 3.2 Principal Component Analysis

```{r mtb-pca}
# PCA
pca.mtb <- pca(mtb_DA, ncomp = 5, center = TRUE, scale = TRUE)
pca.mtb

my_cols <- c("#1b9e77", "#d95f02", "#7570b3")

plotIndiv(
  pca.mtb,
  group = Group,
  ind.names = FALSE,
  ellipse = TRUE,
  legend = TRUE,
  pch = 2,
  cex = 1.4,
  col = my_cols,
  title = "PCA of Metabolomics Data"
)

plotVar(pca.mtb, var.names = FALSE, title = "Metabolite Loadings")


```

```{r}
# Loadings: one column per component, one row per metabolite
loadings_mtb <- pca.mtb$loadings$X   # matrix: metabolites × components

# Component 1
load_pc1 <- loadings_mtb[, 1]
# Component 2
load_pc2 <- loadings_mtb[, 2]

get_top_loadings <- function(load_vec, n = 10) {
  ord <- order(abs(load_vec), decreasing = TRUE)
  data.frame(
    Metabolite = names(load_vec)[ord][1:n],
    Loading    = load_vec[ord][1:n]
  )
}

top_pc1 <- get_top_loadings(load_pc1, n = 10)
top_pc2 <- get_top_loadings(load_pc2, n = 10)

top_pc1
top_pc2
```

## 3.3 Demographic characteristics

```{r covariate-summary-functions}

make_covariate_block <- function(covariate, cov_name, group) {
  tab <- table(group, covariate)
  test <- fisher.test(tab)

  df <- as.data.frame(tab)
  colnames(df) <- c("Group", "Level", "Count")

  df_wide <- df |>
    pivot_wider(
      names_from  = Group,
      values_from = Count,
      values_fill = 0
    ) |>
    mutate(
      Control   = as.character(Control),
      Adenoma   = as.character(Adenoma),
      Carcinoma = as.character(Carcinoma),
      Level     = paste0("  ", Level),
      `P value` = ""
    ) |>
    arrange(Level) |>
    select(Level, Control, Adenoma, Carcinoma, `P value`)

  header <- tibble(
    Level     = cov_name,
    Control   = "",
    Adenoma   = "",
    Carcinoma = "",
    `P value` = as.character(signif(test$p.value, 3))
  )

  bind_rows(header, df_wide)
}

Group   <- factor(metadata_aligned$Study.Group,
                  levels = c("Control", "Adenoma", "Carcinoma"))
AgeCat  <- factor(metadata_aligned$Age)
Sex     <- factor(metadata_aligned$Gender)
Race    <- factor(metadata_aligned$Race)
Smoking <- factor(metadata_aligned$`Chem ID / Smoking history`)

group_counts <- table(Group)

total_row <- tibble(
  Level     = "Total",
  Control   = paste0("n = ", group_counts["Control"]),
  Adenoma   = paste0("n = ", group_counts["Adenoma"]),
  Carcinoma = paste0("n = ", group_counts["Carcinoma"]),
  `P value` = ""
)

age_block   <- make_covariate_block(AgeCat,  "Age (yr)",        Group)
sex_block   <- make_covariate_block(Sex,     "Sex",             Group)
race_block  <- make_covariate_block(Race,    "Race",            Group)
smoke_block <- make_covariate_block(Smoking, "Smoking history", Group)

demog_table <- bind_rows(
  total_row,
  age_block,
  sex_block,
  race_block,
  smoke_block
)

kable_out <- knitr::kable(
  demog_table,
  col.names = c("Demographic characteristic", "Control", "Adenoma", "Carcinoma", "P value"),
  align = c("l", "c", "c", "c", "c"),
  caption = "No. of patients in group by demographic characteristic.",
  escape = FALSE,
  na = ""
)


kable_out |>
  kableExtra::kable_styling(full_width = FALSE) |>
  kableExtra::column_spec(2:4, width = "1.4cm")

```




## 3.4 Differential Expression Analysis (metabolomics with Limma)

```{r limma-dea-function}
# Function to perform limma differential analysis
perform_limma_analysis <- function(data, group, levels_order) {
  X <- t(as.matrix(data))  # features x samples
  group <- factor(group, levels = levels_order)
  stopifnot(ncol(X) == length(group))

  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)

  fit <- lmFit(X, design)
  cont.matrix <- makeContrasts(
    Adenoma_vs_Control = Adenoma - Control,
    Carcinoma_vs_Control = Carcinoma - Control,
    Carcinoma_vs_Adenoma = Carcinoma - Adenoma,
    levels = design
  )

  fit2 <- contrasts.fit(fit, cont.matrix)
  eBayes(fit2)
}

extract_toptable <- function(fit, contrast_name) {
  tt <- topTable(fit, coef = contrast_name, number = Inf, adjust.method = "BH")
  tt$FC <- 2^(tt$logFC)
  tt$negLog10FDR <- -log10(pmax(tt$adj.P.Val, .Machine$double.xmin))
  tt
}
```

```{r run-limma-dea}
# Perform differential analysis
fit2 <- perform_limma_analysis(
  mtb_DA, 
  metadata_aligned$Study.Group,
  c("Control", "Adenoma", "Carcinoma")
)

# Extract results for all contrasts
tt_adenoma <- extract_toptable(fit2, "Adenoma_vs_Control")
tt_carcinoma_control <- extract_toptable(fit2, "Carcinoma_vs_Control")
tt_carcinoma_adenoma <- extract_toptable(fit2, "Carcinoma_vs_Adenoma")

# Display top results
#cat("\nTop differential metabolites (Adenoma vs Control):\n")
datatable(tt_adenoma[tt_adenoma$adj.P.Val < 0.23, c("logFC", "adj.P.Val", "FC")])
```
```{r}
# Adenoma vs Control
ggplot(tt_adenoma, aes(x = logFC, y = negLog10FDR)) +
  geom_point(aes(color = ifelse(adj.P.Val < 0.23, 
                                 ifelse(logFC < 0, "Decreased", "Increased"), 
                                 "Not significant")), 
             alpha = 0.7, size = 2.5) +
  scale_color_manual(values = c("Decreased" = "red", 
                                 "Increased" = "blue", 
                                 "Not significant" = "grey")) +
  geom_hline(yintercept = -log10(0.23), linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.3) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(size = 12),
    legend.position = "top"
  ) +
  labs(
    x = "Log2 Fold Change",
    y = "-log10(FDR)",
    color = "",
    title = "Volcano Plot: Adenoma vs Control"
  )

# Carcinoma vs Control
ggplot(tt_carcinoma_control, aes(x = logFC, y = negLog10FDR)) +
  geom_point(aes(color = ifelse(adj.P.Val < 0.23, 
                                 ifelse(logFC < 0, "Decreased", "Increased"), 
                                 "Not significant")), 
             alpha = 0.7, size = 2.5) +
  scale_color_manual(values = c("Decreased" = "red", 
                                 "Increased" = "blue", 
                                 "Not significant" = "grey")) +
  geom_hline(yintercept = -log10(0.23), linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.3) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(size = 12),
    legend.position = "top"
  ) +
  labs(
    x = "Log2 Fold Change",
    y = "-log10(FDR)",
    color = "",
    title = "Volcano Plot: Carcinoma vs Control"
  )

# Carcinoma vs Adenoma
ggplot(tt_carcinoma_adenoma, aes(x = logFC, y = negLog10FDR)) +
  geom_point(aes(color = ifelse(adj.P.Val < 0.23, 
                                 ifelse(logFC < 0, "Decreased", "Increased"), 
                                 "Not significant")), 
             alpha = 0.7, size = 2.5) +
  scale_color_manual(values = c("Decreased" = "red", 
                                 "Increased" = "blue", 
                                 "Not significant" = "grey")) +
  geom_hline(yintercept = -log10(0.23), linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.3) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(size = 12),
    legend.position = "top"
  ) +
  labs(
    x = "Log2 Fold Change",
    y = "-log10(FDR)",
    color = "",
    title = "Volcano Plot: Carcinoma vs Adenoma"
  )

```


## 3.5 Global Testing

```{r globaltest-function}
# Match samples to metadata
mtb_DA <- mtb[rownames(mtb) %in% metadata$Sample, , drop = FALSE]

# Reorder metadata to match sample order in mtb_DA
idx <- match(rownames(mtb_DA), metadata$Sample)
metadata_aligned <- metadata[idx, ]

# Verify alignment
stopifnot(all(metadata_aligned$Sample == rownames(mtb_DA)))

# Prepare data (samples x metabolites)
X <- mtb_DA  # already samples x metabolites, no need to transpose
stopifnot(all(rownames(X) == metadata_aligned$Sample))

# Use metadata_aligned directly
meta <- metadata_aligned
meta$Study.Group <- factor(meta$Study.Group)

# Subset to Control and Adenoma
sub_AC <- meta$Study.Group %in% c("Control", "Adenoma","Carcinoma" )
X_AC <- X[sub_AC, , drop = FALSE]
meta_AC <- droplevels(meta[sub_AC, , drop = FALSE])

# Check levels
levels(meta_AC$Study.Group)

# Globaltest without covariates
gt_AC <- gt(meta_AC$Study.Group, X_AC)
summary(gt_AC)

# Globaltest with covariates (Age, Sex)
gt_AC_adj <- gt(Study.Group ~ Age + Gender, X_AC, data = meta_AC)
summary(gt_AC_adj)
```

# 4. Metagenomics Analysis

```{r prepare-metagenomics}
meg_aligned <- align_data(Meg, metadata)
meg <- meg_aligned$data
metadata_aligned_meg <- meg_aligned$metadata

# Filter low-abundance features
meg <- meg[, colSums(meg) > 10]
cat("Metagenomics data dimensions after filtering:", dim(meg), "\n")

# Remove samples with zero total counts
meg <- meg[rowSums(meg) > 10, ]

# Relative abundance
meg_rel <- sweep(meg, 1, rowSums(meg), FUN = "/")
row_sums_check <- rowSums(meg_rel)
summary(row_sums_check)
```
```{r}
sample_totals <- rowSums(meg)
summary(sample_totals)
hist(sample_totals, main = "Total metagenomic counts per sample",
     xlab = "Total counts", col = "grey")
```

# 4.1 Alpha Diversity
```{r}
library(vegan)
# Shannon and richness
shannon  <- diversity(meg, index = "shannon")
richness <- specnumber(meg)

alpha_df <- data.frame(
  Sample   = rownames(meg),
  Shannon  = shannon,
  Richness = richness
) %>%
  left_join(metadata_aligned_meg, by = "Sample")

ggplot(alpha_df, aes(x = Study.Group, y = Shannon, fill = Study.Group)) +
  geom_boxplot(alpha = 0.7) +
  theme_bw() +
  labs(x = "Group", y = "Shannon diversity", title = "Alpha diversity (Shannon)")

ggplot(alpha_df, aes(x = Study.Group, y = Richness, fill = Study.Group)) +
  geom_boxplot(alpha = 0.7) +
  theme_bw() +
  labs(x = "Group", y = "Richness", title = "Alpha diversity (Richness)")

# Kruskal-Wallis tests
kruskal_shannon  <- kruskal.test(Shannon ~ Study.Group, data = alpha_df)
kruskal_richness <- kruskal.test(Richness ~ Study.Group, data = alpha_df)

kruskal_shannon
kruskal_richness
```


## 4.1 Beta Diversity (PCoA)

```{r beta-diversity}

# Bray-Curtis distance
dist_bc <- vegdist(meg_rel, method = "bray")

pcoa_res <- cmdscale(dist_bc, eig = TRUE, k = 2, add = TRUE)

eig_vals <- pcoa_res$eig
pos_eig  <- eig_vals[eig_vals > 0]
var_expl <- 100 * eig_vals[1:2] / sum(pos_eig)

pcoa_df <- data.frame(
  Sample = rownames(meg_rel),
  PC1    = pcoa_res$points[, 1],
  PC2    = pcoa_res$points[, 2]
) %>%
  left_join(metadata_aligned_meg, by = "Sample")

p <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Study.Group)) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_color_manual(values = my_cols) +
  theme_bw() +
  labs(
    x = paste0("PCoA1 (", round(var_expl[1], 1), "%)"),
    y = paste0("PCoA2 (", round(var_expl[2], 1), "%)"),
    color = "Group",
    title = "PCoA of Metagenomic Profiles"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank()
  )

print(p)

# PERMANOVA
adonis_res <- adonis2(dist_bc ~ Study.Group + Age + Gender,
                      data = metadata_aligned_meg)
adonis_res

```
## 4.2 DEA 

```{r}
library(ANCOMBC)

# Counts: taxa x samples
meg_counts     <- t(meg)
abundance_data <- as.matrix(meg_counts)
storage.mode(abundance_data) <- "double"

metadata$Sample <- as.character(metadata$Sample)

common_ids <- intersect(colnames(abundance_data), metadata$Sample)
abundance_data <- abundance_data[, common_ids, drop = FALSE]

idx <- match(common_ids, metadata$Sample)

Sample_vec      <- metadata$Sample[idx]
StudyGroup_vec  <- metadata$Study.Group[idx]
Age_vec         <- metadata$Age[idx]
Gender_vec      <- metadata$Gender[idx]
Race_vec        <- metadata$Race[idx]
Smoking_vec     <- metadata$`Chem ID / Smoking history`[idx]

meta_data <- data.frame(
  Sample      = Sample_vec,
  Study.Group = factor(StudyGroup_vec, levels = c("Control", "Adenoma", "Carcinoma")),
  Age         = Age_vec,
  Gender      = factor(Gender_vec),
  Race        = factor(Race_vec),
  Smoking     = factor(Smoking_vec),
  stringsAsFactors = FALSE
)

rownames(meta_data) <- meta_data$Sample

identical(colnames(abundance_data), rownames(meta_data))

```

```{r}
out_meg <- ancombc2(
  data        = abundance_data,
  meta_data   = meta_data,
  fix_formula = "Study.Group + Age + Gender",     # you can also try: "Study.Group + Age + Gender"
  group       = "Study.Group",
  p_adj_method = "BH",
  prv_cut      = 0.10,
  lib_cut      = 1000,
  struc_zero   = TRUE,
  neg_lb       = TRUE,
  alpha        = 0.05,
  global       = TRUE
)

res <- out_meg$res
```


## Adenoma volcano
```{r}
volcano_adenoma <- res %>%
  transmute(
    taxon,
    logFC     = lfc_Study.GroupAdenoma,
    adj.P.val = q_Study.GroupAdenoma
  ) %>%
  mutate(
    negLog10FDR = -log10(adj.P.val),
    FC          = exp(logFC),
    sig         = adj.P.val < 0.05,
    direction   = dplyr::case_when(
      sig & logFC > 0  ~ "Up in Adenoma",
      sig & logFC < 0  ~ "Down in Adenoma",
      TRUE             ~ "Not significant"
    )
  )

ggplot(volcano_adenoma,
       aes(x = logFC, y = negLog10FDR, color = direction)) +
  geom_point(alpha = 0.7, size = 1.8) +
  scale_color_manual(
    values = c(
      "Up in Adenoma"   = "red",
      "Down in Adenoma" = "blue",
      "Not significant" = "grey70"
    )
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.3) +
  theme_bw() +
  labs(
    x = "Log fold-change (Adenoma vs reference)",
    y = "-log10(FDR)",
    color = "Regulation",
    title = "ANCOM-BC2 volcano plot: Adenoma contrast"
  )

sig_adenoma <- volcano_adenoma %>%
  filter(adj.P.val < 0.05) %>%
  arrange(adj.P.val) %>%
  select(taxon, logFC, FC, adj.P.val)

datatable(sig_adenoma)


```

```{r}
volcano_carcinoma <- res %>%
  transmute(
    taxon,
    logFC     = lfc_Study.GroupCarcinoma,
    adj.P.val = q_Study.GroupCarcinoma
  ) %>%
  mutate(
    negLog10FDR = -log10(adj.P.val),
    FC          = exp(logFC),
    sig         = adj.P.val < 0.05,
    direction   = dplyr::case_when(
      sig & logFC > 0  ~ "Up in Carcinoma",
      sig & logFC < 0  ~ "Down in Carcinoma",
      TRUE             ~ "Not significant"
    )
  )

ggplot(volcano_carcinoma,
       aes(x = logFC, y = negLog10FDR, color = direction)) +
  geom_point(alpha = 0.7, size = 1.8) +
  scale_color_manual(
    values = c(
      "Up in Carcinoma"   = "red",
      "Down in Carcinoma" = "blue",
      "Not significant"   = "grey70"
    )
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.3) +
  theme_bw() +
  labs(
    x = "Log fold-change (Carcinoma vs reference)",
    y = "-log10(FDR)",
    color = "Regulation",
    title = "ANCOM-BC2 volcano plot: Carcinoma contrast"
  )

sig_carcinoma <- volcano_carcinoma %>%
  filter(adj.P.val < 0.05) %>%
  arrange(adj.P.val) %>%
  select(taxon, logFC, FC, adj.P.val)

datatable(sig_carcinoma)
```


```{r}
df_fig_adenoma <- res %>%
  transmute(
    taxon,
    lfc = lfc_Study.GroupAdenoma,
    se  = se_Study.GroupAdenoma,
    q   = q_Study.GroupAdenoma
  ) %>%
  mutate(
    sig = q < 0.05,
    direction = dplyr::case_when(
      sig & lfc > 0 ~ "Up in Adenoma",
      sig & lfc < 0 ~ "Down in Adenoma",
      TRUE          ~ "Not significant"
    )
  ) %>%
  filter(sig) %>%
  arrange(lfc)

df_fig_adenoma$taxon <- factor(df_fig_adenoma$taxon,
                               levels = df_fig_adenoma$taxon)

ggplot(df_fig_adenoma,
       aes(x = taxon, y = lfc, fill = direction)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se),
                width = 0.2, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  theme_bw() +
  labs(
    x = "Taxon",
    y = "Log fold-change (Adenoma vs reference)",
    fill = "Regulation",
    title = "Significant taxa (FDR < 0.05) – Adenoma contrast"
  )

df_fig_carcinoma <- res %>%
  transmute(
    taxon,
    lfc = lfc_Study.GroupCarcinoma,
    se  = se_Study.GroupCarcinoma,
    q   = q_Study.GroupCarcinoma
  ) %>%
  mutate(
    sig = q < 0.05,
    direction = dplyr::case_when(
      sig & lfc > 0 ~ "Up in Carcinoma",
      sig & lfc < 0 ~ "Down in Carcinoma",
      TRUE          ~ "Not significant"
    )
  ) %>%
  filter(sig) %>%
  arrange(lfc)

df_fig_carcinoma$taxon <- factor(df_fig_carcinoma$taxon,
                                 levels = df_fig_carcinoma$taxon)

ggplot(df_fig_carcinoma,
       aes(x = taxon, y = lfc, fill = direction)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se),
                width = 0.2, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  theme_bw() +
  labs(
    x = "Taxon",
    y = "Log fold-change (Carcinoma vs reference)",
    fill = "Regulation",
    title = "Significant taxa (FDR < 0.05) – Carcinoma contrast"
  )

```


# 5. Multi-Omics Integration with DIABLO
```{r diablo-full, fig.width=8, fig.height=7}
library(mixOmics)

# 1. Common samples
common_samples <- Reduce(intersect, list(
  rownames(mtb_DA),
  rownames(meg_rel),
  metadata_aligned$Sample
))
common_samples <- sort(common_samples)

metabolome  <- mtb_DA[common_samples, , drop = FALSE]
metagenome  <- meg_rel[common_samples, , drop = FALSE]
meta_diablo <- metadata_aligned[match(common_samples, metadata_aligned$Sample), , drop = FALSE]

stopifnot(identical(rownames(metabolome), rownames(metagenome)))
stopifnot(identical(rownames(metabolome), meta_diablo$Sample))

Y <- factor(meta_diablo$Study.Group,
            levels = c("Control", "Adenoma", "Carcinoma"))

# 2. Optional near-zero variance filtering (if caret loaded)
if ("caret" %in% (.packages())) {
  nzv_mtb <- caret::nearZeroVar(metabolome)
  if (length(nzv_mtb) > 0) {
    metabolome <- metabolome[, -nzv_mtb, drop = FALSE]
  }
  
  nzv_meg <- caret::nearZeroVar(metagenome)
  if (length(nzv_meg) > 0) {
    metagenome <- metagenome[, -nzv_meg, drop = FALSE]
  }
}

X <- list(
  metabolomics = as.matrix(metabolome),
  metagenomics = as.matrix(metagenome)
)

stopifnot(identical(rownames(X$metabolomics), rownames(X$metagenomics)))
stopifnot(length(Y) == nrow(X$metabolomics))

# 3. Design matrix
design <- matrix(c(
  0,   0.1,
  0.1, 0
), nrow = 2, byrow = TRUE)
colnames(design) <- rownames(design) <- names(X)

ncomp_diablo <- 2

test.keepX <- list(
  metabolomics = c(3, 7),
  metagenomics = c(3, 7)
)

set.seed(123)

tune.diablo <- tune.block.splsda(
  X           = X,
  Y           = Y,
  ncomp       = ncomp_diablo,
  test.keepX  = test.keepX,
  design      = design,
  validation  = "Mfold",
  folds       = 10,
  nrepeat     = 10,
  dist        = "max.dist",
  measure     = "BER",
  progressBar = TRUE
)

choice.keepX <- tune.diablo$choice.keepX
choice.keepX

# 4. Final DIABLO model
final.diablo <- block.splsda(
  X             = X,
  Y             = Y,
  ncomp         = ncomp_diablo,
  keepX         = choice.keepX,
  design        = design,
  near.zero.var = TRUE
)

# 5. Performance
set.seed(123)

perf.diablo <- perf(
  final.diablo,
  validation = "Mfold",
  folds      = 10,
  nrepeat    = 10,
  dist       = "max.dist",
  progressBar = TRUE
)

perf.diablo$error.rate

# 6. AUC / ROC
auc.metabolomics.comp1 <- auroc(
  final.diablo,
  roc.block = "metabolomics",
  roc.comp  = 1
)
auc.metagenomics.comp1 <- auroc(
  final.diablo,
  roc.block = "metagenomics",
  roc.comp  = 1
)

auc.metabolomics.comp1$auc
auc.metagenomics.comp1$auc

auc.metabolomics.comp2 <- auroc(
  final.diablo,
  roc.block = "metabolomics",
  roc.comp  = 2
)
auc.metagenomics.comp2 <- auroc(
  final.diablo,
  roc.block = "metagenomics",
  roc.comp  = 2
)

auc.metabolomics.comp2$auc
auc.metagenomics.comp2$auc

# 7. Plots
plotIndiv(
  final.diablo,
  ind.names = FALSE,
  legend    = TRUE,
  group     = Y,
  ellipse   = TRUE,
  title     = "DIABLO sample plot"
)


par(mar = c(4, 10, 4, 4))
p <-circosPlot(
  final.diablo,
  cutoff       = 0.6,
  color.blocks = c("#1b9e77", "#d95f02"),
  size.variables = 0.6   # default is ~0.25; increase for bigger tex     # block labels (optional)
)


```
```{r}
## 8. Predictions for sensitivity / specificity and BER (per component)

pred <- predict(
  final.diablo,
  newdata = X,
  dist    = "max.dist"
)

pred_class <- pred$class$max.dist   # list with $metabolomics, $metagenomics [web:230]

## Helper: confusion matrix -> BER
conf_to_ber <- function(truth, predicted) {
  conf <- get.confusion_matrix(
    truth      = truth,
    all.levels = levels(truth),
    predicted  = predicted
  )
  ber <- get.BER(conf)
  list(confusion = conf, BER = ber)
}

## Helper: sensitivity & specificity (one-vs-rest)
calc_sens_spec <- function(truth, predicted, positive_level) {
  pos     <- truth == positive_level
  predpos <- predicted == positive_level
  
  TP <- sum(pos & predpos, na.rm = TRUE)
  FN <- sum(pos & !predpos, na.rm = TRUE)
  FP <- sum(!pos & predpos, na.rm = TRUE)
  TN <- sum(!pos & !predpos, na.rm = TRUE)
  
  sensitivity <- ifelse((TP + FN) > 0, TP / (TP + FN), NA)
  specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), NA)
  
  data.frame(
    class       = positive_level,
    sensitivity = sensitivity,
    specificity = specificity
  )
}

## Component 1 (using metabolomics predictions; you can change the rule)

comp1_metab <- pred_class$metabolomics[, "comp1"]

yhat_comp1 <- comp1_metab

res_comp1  <- conf_to_ber(Y, yhat_comp1)
conf_comp1 <- res_comp1$confusion
BER_comp1  <- res_comp1$BER

conf_comp1
BER_comp1

sens_spec_comp1 <- bind_rows(
  calc_sens_spec(Y, yhat_comp1, "Control"),
  calc_sens_spec(Y, yhat_comp1, "Adenoma"),
  calc_sens_spec(Y, yhat_comp1, "Carcinoma")
)

sens_spec_comp1

## Component 2 (same approach)

comp2_metab <- pred_class$metabolomics[, "comp2"]

yhat_comp2 <- comp2_metab

res_comp2  <- conf_to_ber(Y, yhat_comp2)
conf_comp2 <- res_comp2$confusion
BER_comp2  <- res_comp2$BER

conf_comp2
BER_comp2

sens_spec_comp2 <- bind_rows(
  calc_sens_spec(Y, yhat_comp2, "Control"),
  calc_sens_spec(Y, yhat_comp2, "Adenoma"),
  calc_sens_spec(Y, yhat_comp2, "Carcinoma")
)

sens_spec_comp2

## 9. Joint table for both components (only sensitivity and specificity)

sens_spec_comp1$component <- "comp1"
sens_spec_comp2$component <- "comp2"

sens_spec_all <- bind_rows(sens_spec_comp1, sens_spec_comp2) %>%
  select(class, component, sensitivity, specificity)

sens_spec_all

## Optional: wide format (one row per class, columns per component)
sens_spec_wide <- pivot_wider(
  sens_spec_all,
  id_cols    = class,
  names_from = component,
  values_from = c(sensitivity, specificity),
  names_sep  = "_"
)

sens_spec_wide

```




```{r}
# 8. Main contributing variables (component 1)
plotLoadings(final.diablo, block = "metabolomics", comp = 1, contrib = "max")
plotLoadings(final.diablo, block = "metagenomics", comp = 1, contrib = "max")

plotLoadings(final.diablo, block = "metabolomics", comp = 2, contrib = "max")
plotLoadings(final.diablo, block = "metagenomics", comp = 2, contrib = "max")
```
# 6 LASSO Classidication of metabolites 

```{r lasso_nested_cv_CA, message=FALSE, warning=FALSE}

library(glmnet)
library(dplyr)
library(pROC)

set.seed(123)

# 1. Prepare binary dataset (Control vs Adenoma only)
keep_idx <- metadata_aligned$Study.Group %in% c("Control", "Adenoma")

X_full <- as.matrix(mtb_DA[keep_idx, , drop = FALSE])
y_full <- factor(metadata_aligned$Study.Group[keep_idx],
                 levels = c("Control", "Adenoma"))

stopifnot(nrow(X_full) == length(y_full))

# 2. Fit full LASSO path (for coefficient paths)
fit_full <- glmnet(
  x           = X_full,
  y           = y_full,
  family      = "binomial",
  alpha       = 1,
  standardize = TRUE
)

# 3. Cross-validation with AUC on all data
set.seed(123)
cv_full_auc <- cv.glmnet(
  x            = X_full,
  y            = y_full,
  family       = "binomial",
  alpha        = 1,
  nfolds       = 10,
  type.measure = "auc",
  standardize  = TRUE
)

lambda_min <- cv_full_auc$lambda.min
lambda_1se <- cv_full_auc$lambda.1se

auc_min <- cv_full_auc$cvm[cv_full_auc$lambda == lambda_min]
auc_1se <- cv_full_auc$cvm[cv_full_auc$lambda == lambda_1se]

# 3b. Cross-validation with misclassification error on all data
set.seed(123)
cv_full_err <- cv.glmnet(
  x            = X_full,
  y            = y_full,
  family       = "binomial",
  alpha        = 1,
  nfolds       = 10,
  type.measure = "class",
  standardize  = TRUE
)

err_min <- cv_full_err$cvm[cv_full_err$lambda == lambda_min]
err_1se <- cv_full_err$cvm[cv_full_err$lambda == lambda_1se]

# 4. Simple train/test split for test metrics
set.seed(456)
n <- nrow(X_full)
train_idx <- sample(seq_len(n), size = floor(0.7 * n))
test_idx  <- setdiff(seq_len(n), train_idx)

X_train <- X_full[train_idx, , drop = FALSE]
y_train <- y_full[train_idx]

X_test  <- X_full[test_idx, , drop = FALSE]
y_test  <- y_full[test_idx]

# CV on training set (AUC)
cv_train <- cv.glmnet(
  x            = X_train,
  y            = y_train,
  family       = "binomial",
  alpha        = 1,
  nfolds       = 10,
  type.measure = "auc",
  standardize  = TRUE
)

lambda_1se_train <- cv_train$lambda.1se

# Predict probabilities on test set with lambda.1se
prob_test <- predict(
  cv_train,
  newx = X_test,
  s    = "lambda.1se",
  type = "response"
)
prob_test <- as.numeric(prob_test)

# Test AUC
roc_obj  <- roc(response = y_test,
                predictor = prob_test,
                levels = c("Control", "Adenoma"),
                direction = "<")
auc_test <- as.numeric(auc(roc_obj))

# Test sensitivity, specificity, misclassification error at threshold 0.5
thr <- 0.5
y_pred <- ifelse(prob_test >= thr, "Adenoma", "Control")
y_pred <- factor(y_pred, levels = levels(y_test))

tab <- table(Truth = y_test, Pred = y_pred)

TP <- tab["Adenoma", "Adenoma"]
FN <- tab["Adenoma", "Control"]
FP <- tab["Control", "Adenoma"]
TN <- tab["Control", "Control"]

sensitivity    <- TP / (TP + FN)
specificity    <- TN / (TN + FP)
misclass_test  <- 1 - sum(diag(tab)) / sum(tab)

# 5. Selected variables at lambda.min (or lambda.1se if you prefer)
coef_min <- coef(cv_full_auc, s = "lambda.min")
coef_min_mat <- as.matrix(coef_min)[-1, , drop = FALSE]  # drop intercept

selected_df <- data.frame(
  feature = rownames(coef_min_mat),
  coef    = as.numeric(coef_min_mat[, 1])
) |>
  dplyr::filter(coef != 0) |>
  dplyr::arrange(dplyr::desc(abs(coef)))

n_selected <- nrow(selected_df)



# Plot 1: coefficient paths
plot(fit_full,
     xvar  = "lambda",
     label = FALSE,
     main  = "LASSO coefficient paths")

# Plot 2: cross-validation curve
plot(cv_full_auc,
     main = "Cross-validation (AUC) for LASSO")

# Optional Plot 3: ROC curve on test set
plot(roc_obj,
     col = "red",
     lwd = 2,
     main = "ROC curve (test set)")


# 6. Simple summary table for slides
perf_table <- data.frame(
  Metric = c("Selected features (lambda.min)",
             "CV AUC (lambda.min)",
             "CV misclassification error (lambda.min)",
             "Test sensitivity (thr = 0.5)",
             "Test specificity (thr = 0.5)",
             "Test misclassification error (thr = 0.5)"),
  Value  = c(
    n_selected,
    round(auc_min, 3),
    round(err_min, 3),
    round(sensitivity, 3),
    round(specificity, 3),
    round(misclass_test, 3)
  )
)

perf_table

```



```{r}
# -------------------------
# 5. Variables selected at lambda.1se
# -------------------------

# Coefficient vector at lambda.1se (includes intercept in first row)
coef_1se <- coef(cv_full_auc, s = "lambda.min")

# Convert to a tidy data frame and drop intercept + zero coefficients
sel_1se <- as.matrix(coef_1se)
sel_1se <- sel_1se[-1, , drop = FALSE]  # remove intercept row

selected_df <- data.frame(
  feature = rownames(sel_1se),
  coef    = as.numeric(sel_1se[, 1])
)

selected_df <- selected_df %>%
  dplyr::filter(coef != 0) %>%         # keep only non-zero coefficients
  dplyr::arrange(dplyr::desc(abs(coef)))

# Print selected features at lambda.1se
print(selected_df)

cat("\nNumber of selected features at lambda.min:", nrow(selected_df), "\n")

```








# Session Info
```{r session-info}
#sessionInfo()
```
