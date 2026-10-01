library(GEOquery)
library(limma)
library(Biobase)

# extract clinical group (NGT / IGT / DM) from the title field
titles <- pData(gse27951)$title
group <- ifelse(grepl("_NGT", titles), "NGT",
                ifelse(grepl("_IGT", titles), "IGT",
                       ifelse(grepl("_DM", titles), "DM", NA)))

table(group)  # sanity check counts per group

pData(gse27951)$group <- group

# keep only NGT and DM samples for the primary comparison (drop IGT for now)
keep_samples <- group %in% c("NGT", "DM")
gse_sub <- gse27951[, keep_samples]

dim(gse_sub)
table(pData(gse_sub)$group)


exprs_data <- exprs(gse_sub)
print(range(exprs_data))

if (max(exprs_data) > 50) {
  exprs_data <- log2(exprs_data + 1)
  exprs(gse_sub) <- exprs_data
}

group_factor <- factor(pData(gse_sub)$group, levels = c("NGT", "DM"))
design <- model.matrix(~ group_factor)
colnames(design) <- c("Intercept", "DM_vs_NGT")

print(design)

fit <- lmFit(gse_sub, design)
fit <- eBayes(fit)

results <- topTable(fit, coef = "DM_vs_NGT", number = Inf, adjust.method = "BH")

print(dim(results))
print(head(results, 20))

write.csv(results, file = "C:/tri_disease_project/results/tables/GSE27951_DE_results.csv", row.names = TRUE)


library(dplyr)

# collapse multiple probes per gene -> keep the probe with lowest p-value per gene (standard practice)
results$Gene.symbol <- as.character(results$Gene.symbol)
results_clean <- results[results$Gene.symbol != "" & !is.na(results$Gene.symbol), ]

gene_level_results <- results_clean %>%
  group_by(Gene.symbol) %>%
  slice_min(P.Value, n = 1) %>%
  ungroup()

dim(gene_level_results)

write.csv(gene_level_results, file = "C:/tri_disease_project/results/tables/GSE27951_DE_gene_level.csv", row.names = FALSE)

# now check against your shared gene list and hub genes
shared_genes <- read.csv("C:/tri_disease_project/data/processed/shared_genes_protein_coding.csv")$hgnc_symbol
hub_genes <- read.csv("C:/tri_disease_project/results/tables/hub_genes_final.csv")$gene_symbol

shared_hits <- gene_level_results[gene_level_results$Gene.symbol %in% shared_genes, ]
hub_hits <- gene_level_results[gene_level_results$Gene.symbol %in% hub_genes, ]

cat("Shared genes found in this dataset:", nrow(shared_hits), "out of", length(shared_genes), "\n")
cat("Hub genes found in this dataset:", nrow(hub_hits), "out of", length(hub_genes), "\n\n")

cat("Hub genes - differential expression in DM vs NGT adipose tissue:\n")
print(hub_hits[order(hub_hits$P.Value), c("Gene.symbol", "logFC", "P.Value", "adj.P.Val")])

cat("\nShared genes significant at raw p<0.05:\n")
print(shared_hits[shared_hits$P.Value < 0.05, c("Gene.symbol", "logFC", "P.Value", "adj.P.Val")])

write.csv(shared_hits, file = "C:/tri_disease_project/results/tables/GSE27951_shared_genes_DE.csv", row.names = FALSE)
write.csv(hub_hits, file = "C:/tri_disease_project/results/tables/GSE27951_hub_genes_DE.csv", row.names = FALSE)

n_shared_tested <- nrow(shared_hits)
n_shared_sig <- sum(shared_hits$P.Value < 0.05)

n_all_tested <- nrow(gene_level_results)
n_all_sig <- sum(gene_level_results$P.Value < 0.05)

cat("Shared genes: ", n_shared_sig, "/", n_shared_tested, " = ", round(100*n_shared_sig/n_shared_tested, 1), "%\n", sep="")
cat("All genes: ", n_all_sig, "/", n_all_tested, " = ", round(100*n_all_sig/n_all_tested, 1), "%\n", sep="")

# binomial test: is the shared-gene hit rate higher than the background (all-gene) hit rate?
background_rate <- n_all_sig / n_all_tested
binom_test <- binom.test(n_shared_sig, n_shared_tested, p = background_rate, alternative = "greater")
print(binom_test)

gse24752_list <- getGEO("GSE24752", GSEMatrix = TRUE, AnnotGPL = TRUE)
names(gse24752_list)

gse24752 <- gse24752_list[[1]]

dim(gse24752)
pData(gse24752)[, c("title", "geo_accession")]

saveRDS(gse24752, file = "C:/tri_disease_project/data/raw/GSE24752_GPL570.rds")

titles_htn <- pData(gse24752)$title
group_htn <- ifelse(grepl("hypertensive", titles_htn), "HTN",
                    ifelse(grepl("normotensive", titles_htn), "Normal", NA))

table(group_htn)
pData(gse24752)$group <- group_htn

exprs_htn <- exprs(gse24752)
print(range(exprs_htn))

if (max(exprs_htn) > 50) {
  exprs_htn <- log2(exprs_htn + 1)
  exprs(gse24752) <- exprs_htn
}

group_factor_htn <- factor(pData(gse24752)$group, levels = c("Normal", "HTN"))
design_htn <- model.matrix(~ group_factor_htn)
colnames(design_htn) <- c("Intercept", "HTN_vs_Normal")

print(design_htn)

fit_htn <- lmFit(gse24752, design_htn)
fit_htn <- eBayes(fit_htn)

results_htn <- topTable(fit_htn, coef = "HTN_vs_Normal", number = Inf, adjust.method = "BH")

dim(results_htn)
head(results_htn, 10)

write.csv(results_htn, file = "C:/tri_disease_project/results/tables/GSE24752_DE_results.csv", row.names = TRUE)

results_htn$Gene.symbol <- as.character(results_htn$Gene.symbol)
results_htn_clean <- results_htn[results_htn$Gene.symbol != "" & !is.na(results_htn$Gene.symbol), ]

gene_level_htn <- results_htn_clean %>%
  group_by(Gene.symbol) %>%
  slice_min(P.Value, n = 1) %>%
  ungroup()

dim(gene_level_htn)

write.csv(gene_level_htn, file = "C:/tri_disease_project/results/tables/GSE24752_DE_gene_level.csv", row.names = FALSE)

shared_hits_htn <- gene_level_htn[gene_level_htn$Gene.symbol %in% shared_genes, ]
hub_hits_htn <- gene_level_htn[gene_level_htn$Gene.symbol %in% hub_genes, ]

cat("Shared genes found in this dataset:", nrow(shared_hits_htn), "out of", length(shared_genes), "\n")
cat("Hub genes found in this dataset:", nrow(hub_hits_htn), "out of", length(hub_genes), "\n\n")

cat("Hub genes - differential expression in HTN vs Normal PBC:\n")
print(hub_hits_htn[order(hub_hits_htn$P.Value), c("Gene.symbol", "logFC", "P.Value", "adj.P.Val")])

n_shared_tested_htn <- nrow(shared_hits_htn)
n_shared_sig_htn <- sum(shared_hits_htn$P.Value < 0.05)
n_all_tested_htn <- nrow(gene_level_htn)
n_all_sig_htn <- sum(gene_level_htn$P.Value < 0.05)

cat("\nShared genes: ", n_shared_sig_htn, "/", n_shared_tested_htn, " = ", round(100*n_shared_sig_htn/n_shared_tested_htn, 1), "%\n", sep="")
cat("All genes: ", n_all_sig_htn, "/", n_all_tested_htn, " = ", round(100*n_all_sig_htn/n_all_tested_htn, 1), "%\n", sep="")

background_rate_htn <- n_all_sig_htn / n_all_tested_htn
binom_test_htn <- binom.test(n_shared_sig_htn, n_shared_tested_htn, p = background_rate_htn, alternative = "greater")
print(binom_test_htn)

write.csv(shared_hits_htn, file = "C:/tri_disease_project/results/tables/GSE24752_shared_genes_DE.csv", row.names = FALSE)
write.csv(hub_hits_htn, file = "C:/tri_disease_project/results/tables/GSE24752_hub_genes_DE.csv", row.names = FALSE)

library(GEOquery)
library(limma)
library(Biobase)

gse75360_list <- getGEO("GSE75360", GSEMatrix = TRUE, AnnotGPL = TRUE)
names(gse75360_list)

gse75360 <- gse75360_list[[1]]

dim(gse75360)
pData(gse75360)[, c("title", "geo_accession")]

saveRDS(gse75360, file = "C:/tri_disease_project/data/raw/GSE75360.rds")

titles_75360 <- pData(gse75360)$title
group_75360 <- ifelse(grepl("hypertensive", titles_75360) & !grepl("normotensive", titles_75360), "HTN", "Normal")
ancestry_75360 <- ifelse(grepl("_AA_", titles_75360), "AA", "CAU")

table(group_75360)
table(group_75360, ancestry_75360)

pData(gse75360)$group <- group_75360
pData(gse75360)$ancestry <- ancestry_75360

exprs_75360 <- exprs(gse75360)
print(range(exprs_75360))

if (max(exprs_75360) > 50) {
  exprs_75360 <- log2(exprs_75360 + 1)
  exprs(gse75360) <- exprs_75360
}

group_factor_75360 <- factor(pData(gse75360)$group, levels = c("Normal", "HTN"))
design_75360 <- model.matrix(~ group_factor_75360)
colnames(design_75360) <- c("Intercept", "HTN_vs_Normal")

print(design_75360)

fit_75360 <- lmFit(gse75360, design_75360)
fit_75360 <- eBayes(fit_75360)

results_75360 <- topTable(fit_75360, coef = "HTN_vs_Normal", number = Inf, adjust.method = "BH")

dim(results_75360)
head(results_75360, 10)

write.csv(results_75360, file = "C:/tri_disease_project/results/tables/GSE75360_DE_results.csv", row.names = TRUE)

library(dplyr)

results_75360$Gene.symbol <- as.character(results_75360$Gene.symbol)
results_75360_clean <- results_75360[results_75360$Gene.symbol != "" & !is.na(results_75360$Gene.symbol), ]

gene_level_75360 <- results_75360_clean %>%
  group_by(Gene.symbol) %>%
  slice_min(P.Value, n = 1) %>%
  ungroup()

dim(gene_level_75360)

write.csv(gene_level_75360, file = "C:/tri_disease_project/results/tables/GSE75360_DE_gene_level.csv", row.names = FALSE)

shared_genes <- read.csv("C:/tri_disease_project/data/processed/shared_genes_protein_coding.csv")$hgnc_symbol
hub_genes <- read.csv("C:/tri_disease_project/results/tables/hub_genes_final.csv")$gene_symbol

shared_hits_75360 <- gene_level_75360[gene_level_75360$Gene.symbol %in% shared_genes, ]
hub_hits_75360 <- gene_level_75360[gene_level_75360$Gene.symbol %in% hub_genes, ]

cat("Shared genes found:", nrow(shared_hits_75360), "/", length(shared_genes), "\n")
cat("Hub genes found:", nrow(hub_hits_75360), "/", length(hub_genes), "\n\n")

cat("Hub genes - differential expression in HTN vs Normal PBMC (GSE75360):\n")
print(hub_hits_75360[order(hub_hits_75360$P.Value), c("Gene.symbol", "logFC", "P.Value", "adj.P.Val")])

n_shared_tested_75360 <- nrow(shared_hits_75360)
n_shared_sig_75360 <- sum(shared_hits_75360$P.Value < 0.05)
n_all_tested_75360 <- nrow(gene_level_75360)
n_all_sig_75360 <- sum(gene_level_75360$P.Value < 0.05)

cat("\nShared genes: ", n_shared_sig_75360, "/", n_shared_tested_75360, " = ", round(100*n_shared_sig_75360/n_shared_tested_75360, 1), "%\n", sep="")
cat("All genes: ", n_all_sig_75360, "/", n_all_tested_75360, " = ", round(100*n_all_sig_75360/n_all_tested_75360, 1), "%\n", sep="")

background_rate_75360 <- n_all_sig_75360 / n_all_tested_75360
binom_test_75360 <- binom.test(n_shared_sig_75360, n_shared_tested_75360, p = background_rate_75360, alternative = "greater")
print(binom_test_75360)

write.csv(shared_hits_75360, file = "C:/tri_disease_project/results/tables/GSE75360_shared_genes_DE.csv", row.names = FALSE)
write.csv(hub_hits_75360, file = "C:/tri_disease_project/results/tables/GSE75360_hub_genes_DE.csv", row.names = FALSE)

gene_level_24752 <- read.csv("C:/tri_disease_project/results/tables/GSE24752_DE_gene_level.csv")

shared_hits_24752 <- gene_level_24752[gene_level_24752$Gene.symbol %in% shared_genes, ]

n_shared_sig_24752 <- sum(shared_hits_24752$P.Value < 0.05)
n_shared_tested_24752 <- nrow(shared_hits_24752)

cat("GSE24752: ", n_shared_sig_24752, "/", n_shared_tested_24752, "\n")
cat("GSE75360: ", n_shared_sig_75360, "/", n_shared_tested_75360, "\n")

# pooled counts across both independent hypertension datasets
pooled_sig <- n_shared_sig_24752 + n_shared_sig_75360
pooled_tested <- n_shared_tested_24752 + n_shared_tested_75360

# pooled background rate (weighted average of both arrays' background rates)
n_all_sig_24752 <- sum(gene_level_24752$P.Value < 0.05)
n_all_tested_24752 <- nrow(gene_level_24752)

pooled_all_sig <- n_all_sig_24752 + n_all_sig_75360
pooled_all_tested <- n_all_tested_24752 + n_all_tested_75360
pooled_background_rate <- pooled_all_sig / pooled_all_tested

cat("\nPooled shared genes: ", pooled_sig, "/", pooled_tested, " = ", round(100*pooled_sig/pooled_tested, 1), "%\n", sep="")
cat("Pooled background rate: ", round(100*pooled_background_rate, 1), "%\n", sep="")

pooled_binom <- binom.test(pooled_sig, pooled_tested, p = pooled_background_rate, alternative = "greater")
print(pooled_binom)

# also run Fisher's method to combine the two independent p-values formally
p1 <- binom_test_75360$p.value
htn_24752_binom <- binom.test(n_shared_sig_24752, n_shared_tested_24752, 
                              p = n_all_sig_24752/n_all_tested_24752, alternative = "greater")
p2 <- htn_24752_binom$p.value

fisher_stat <- -2 * (log(p1) + log(p2))
fisher_p <- pchisq(fisher_stat, df = 4, lower.tail = FALSE)

cat("\nFisher's combined test across both independent hypertension datasets:\n")
cat("GSE24752 p =", p2, "\n")
cat("GSE75360 p =", p1, "\n")
cat("Fisher's combined statistic:", fisher_stat, "\n")
cat("Fisher's combined p-value:", fisher_p, "\n")

combined_htn_summary <- data.frame(
  Dataset = c("GSE24752 (n=3 vs 3)", "GSE75360 (n=10 vs 11)", "Pooled (both datasets)", "Fisher's combined method"),
  Shared_genes_significant = c(20, 74, 94, NA),
  Shared_genes_tested = c(648, 654, 1302, NA),
  Percent = c(round(100*20/648,1), round(100*74/654,1), round(100*94/1302,1), NA),
  P_value = c(0.1102, 0.1188, 0.0320, 0.0698)
)
print(combined_htn_summary)
write.csv(combined_htn_summary, file = "C:/tri_disease_project/results/tables/hypertension_combined_summary.csv", row.names = FALSE)

# subset to AA samples only, rerun the full DE pipeline
gse75360_AA <- gse75360[, pData(gse75360)$ancestry == "AA"]
table(pData(gse75360_AA)$group)

group_factor_AA <- factor(pData(gse75360_AA)$group, levels = c("Normal", "HTN"))
design_AA <- model.matrix(~ group_factor_AA)
colnames(design_AA) <- c("Intercept", "HTN_vs_Normal")

fit_AA <- lmFit(gse75360_AA, design_AA)
fit_AA <- eBayes(fit_AA)
results_AA <- topTable(fit_AA, coef = "HTN_vs_Normal", number = Inf, adjust.method = "BH")

results_AA$Gene.symbol <- as.character(results_AA$Gene.symbol)
results_AA_clean <- results_AA[results_AA$Gene.symbol != "" & !is.na(results_AA$Gene.symbol), ]
gene_level_AA <- results_AA_clean %>% group_by(Gene.symbol) %>% slice_min(P.Value, n = 1) %>% ungroup()

shared_hits_AA <- gene_level_AA[gene_level_AA$Gene.symbol %in% shared_genes, ]
n_shared_sig_AA <- sum(shared_hits_AA$P.Value < 0.05)
n_shared_tested_AA <- nrow(shared_hits_AA)
n_all_sig_AA <- sum(gene_level_AA$P.Value < 0.05)
n_all_tested_AA <- nrow(gene_level_AA)

cat("African American subgroup (n=5 vs 6):\n")
cat("Shared genes: ", n_shared_sig_AA, "/", n_shared_tested_AA, " = ", round(100*n_shared_sig_AA/n_shared_tested_AA,1), "%\n", sep="")
cat("Background: ", round(100*n_all_sig_AA/n_all_tested_AA,1), "%\n", sep="")

binom_AA <- binom.test(n_shared_sig_AA, n_shared_tested_AA, p = n_all_sig_AA/n_all_tested_AA, alternative = "greater")
print(binom_AA)

write.csv(gene_level_AA, file = "C:/tri_disease_project/results/tables/GSE75360_AA_DE_gene_level.csv", row.names = FALSE)

gse75360_CAU <- gse75360[, pData(gse75360)$ancestry == "CAU"]
table(pData(gse75360_CAU)$group)

group_factor_CAU <- factor(pData(gse75360_CAU)$group, levels = c("Normal", "HTN"))
design_CAU <- model.matrix(~ group_factor_CAU)
colnames(design_CAU) <- c("Intercept", "HTN_vs_Normal")

fit_CAU <- lmFit(gse75360_CAU, design_CAU)
fit_CAU <- eBayes(fit_CAU)
results_CAU <- topTable(fit_CAU, coef = "HTN_vs_Normal", number = Inf, adjust.method = "BH")

results_CAU$Gene.symbol <- as.character(results_CAU$Gene.symbol)
results_CAU_clean <- results_CAU[results_CAU$Gene.symbol != "" & !is.na(results_CAU$Gene.symbol), ]
gene_level_CAU <- results_CAU_clean %>% group_by(Gene.symbol) %>% slice_min(P.Value, n = 1) %>% ungroup()

shared_hits_CAU <- gene_level_CAU[gene_level_CAU$Gene.symbol %in% shared_genes, ]
n_shared_sig_CAU <- sum(shared_hits_CAU$P.Value < 0.05)
n_shared_tested_CAU <- nrow(shared_hits_CAU)
n_all_sig_CAU <- sum(gene_level_CAU$P.Value < 0.05)
n_all_tested_CAU <- nrow(gene_level_CAU)

cat("Caucasian subgroup (n=5 vs 5):\n")
cat("Shared genes: ", n_shared_sig_CAU, "/", n_shared_tested_CAU, " = ", round(100*n_shared_sig_CAU/n_shared_tested_CAU,1), "%\n", sep="")
cat("Background: ", round(100*n_all_sig_CAU/n_all_tested_CAU,1), "%\n", sep="")

binom_CAU <- binom.test(n_shared_sig_CAU, n_shared_tested_CAU, p = n_all_sig_CAU/n_all_tested_CAU, alternative = "greater")
print(binom_CAU)

# final ancestry comparison table
ancestry_summary <- data.frame(
  Subgroup = c("African American (n=5 vs 6)", "Caucasian (n=5 vs 5)"),
  Shared_significant = c(n_shared_sig_AA, n_shared_sig_CAU),
  Shared_tested = c(n_shared_tested_AA, n_shared_tested_CAU),
  Percent = c(round(100*n_shared_sig_AA/n_shared_tested_AA,1), round(100*n_shared_sig_CAU/n_shared_tested_CAU,1)),
  Background_pct = c(round(100*n_all_sig_AA/n_all_tested_AA,1), round(100*n_all_sig_CAU/n_all_tested_CAU,1)),
  P_value = c(0.00219, NA)
)
write.csv(ancestry_summary, file = "C:/tri_disease_project/results/tables/ancestry_stratified_summary.csv", row.names = FALSE)
print(ancestry_summary)

ancestry_summary$P_value[2] <- 0.03858
print(ancestry_summary)
write.csv(ancestry_summary, file = "C:/tri_disease_project/results/tables/ancestry_stratified_summary.csv", row.names = FALSE)