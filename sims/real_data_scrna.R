# Real-data scRNA-seq application (the CANONICAL motivation for clustering
# inference; Chen-Witten used 10x PBMC). Loaded programmatically via the scRNAseq
# package -> Pollen et al. (2014) developing-cortex cells with 4 inferred cell
# types (Interneuron / IPC / Neuron / RG). Manual preprocessing (lib-size
# normalize -> log1p -> top HVGs -> PCA) avoids the heavy DropletUtils/scater stack
# the inherited vignette needs. We k-means cluster the cells, then test each
# cluster pair: naive (invalid) / known-sigma (sigma_MED) PATH & UNION / studentized.
# Fixes the inherited vignette gap (it used only the path test + a dead data path).
source("sims/house_style.R")
suppressMessages({ library(KmeansInference); library(scRNAseq); library(SingleCellExperiment)
                   library(ggplot2); library(gridExtra) })
source("sims/kmeans_union_unknownvar.R")
seed <- 2021; per <- 55; QPC <- 10
sce <- PollenGliaData()
ct0 <- as.character(colData(sce)[["Inferred Cell Type"]])
counts <- as.matrix(assay(sce, "counts"))
## balanced subsample per cell type (keeps n tractable for union exploration)
set.seed(seed)
keep <- unlist(lapply(sort(unique(ct0)), function(t) { idx <- which(ct0 == t); sample(idx, min(per, length(idx))) }))
counts <- counts[, keep]; ct <- ct0[keep]
## preprocess: CP10k library-size normalize -> log1p -> top-500 HVG -> PCA -> top QPC PCs (standardized)
ls <- pmax(colSums(counts), 1); cn <- log1p(t(t(counts) / ls * median(ls)))
hvg <- order(apply(cn, 1, var), decreasing = TRUE)[1:500]
pca <- prcomp(t(cn[hvg, ]), center = TRUE, scale. = TRUE)
X <- scale(pca$x[, 1:QPC])                                 # n x QPC feature matrix
K <- as.integer(Sys.getenv("K", as.character(length(unique(ct)))))
cl <- kmeans_estimation(X, K, 10, seed, verbose = FALSE)$final_cluster
s_med <- KmeansInference:::.kmeans_estimate_MED(X)
cat(sprintf("Pollen scRNA: n=%d cells, %d genes -> q=%d PCs, K=%d, sigma_MED=%.3f\n",
            nrow(X), nrow(counts), QPC, K, s_med))
cat("cell types:", paste(names(table(ct)), table(ct), sep="="), "| k-means sizes:", paste(table(cl), collapse="/"), "\n")

prs <- combn(K, 2, simplify = FALSE)
fmt <- function(p) ifelse(p < 1e-4, sprintf("%.1e", p), sprintf("%.4f", p))
tab <- do.call(rbind, lapply(prs, function(pr) {
  c1 <- pr[1]; c2 <- pr[2]
  ft <- tryCatch(kmeans_inference_union(X, K, c1, c2, sig = s_med, seed = seed, max_paths = 3000), error=function(e) NULL)
  rf <- tryCatch(kmeans_union_unknownvar(X, K, c1, c2, seed = seed), error=function(e) NULL)
  if (is.null(ft) || is.null(rf)) return(NULL)
  data.frame(pair = sprintf("%d vs %d", c1, c2), n1 = sum(cl==c1), n2 = sum(cl==c2),
             naive = ft$p_naive, path_known = ft$pval_path, union_known = ft$pval_union,
             path_studF = rf$p_path, union_studF = rf$p_union)
}))
out <- tab; for (j in 4:8) out[[j]] <- fmt(tab[[j]])
cat("\n=== Pollen scRNA cluster-pair selective p-values ===\n"); print(out, row.names = FALSE)
write.csv(tab, sprintf("sims/results/scrna_pvalues_k%d.csv", K), row.names = FALSE)

## ---- figure: PCA scatter + p-value comparison (same grammar as penguins) -----
clpal <- c("#882255","#44AA99","#DDCC77","#332288","#AA4499")[seq_len(K)]; names(clpal) <- as.character(seq_len(K))
fdf <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], cluster = factor(cl), celltype = ct)
pA <- ggplot(fdf, aes(PC1, PC2, colour = cluster, shape = celltype)) +
  geom_point(size = 1.9, alpha = 0.85) +
  scale_colour_manual(values = clpal, name = "k-means cluster") +
  scale_shape_discrete(name = "inferred cell type") +
  labs(title = sprintf("Pollen scRNA-seq: k-means (K=%d)", K)) +
  theme_km() + theme(legend.position = "right",
                     legend.title = element_text(face = "bold", size = 9), legend.text = element_text(size = 8))
nlp <- function(p) -log10(pmax(p, 1e-300))
mlt <- do.call(rbind, lapply(seq_len(nrow(tab)), function(i) data.frame(
  pair = tab$pair[i], treat = c("med","med","studentized","studentized"), cond = c("path","union","path","union"),
  nlp = nlp(c(tab$path_known[i], tab$union_known[i], tab$path_studF[i], tab$union_studF[i])))))
mlt$treat <- factor(mlt$treat, levels = c("med","studentized")); mlt$cond <- factor(mlt$cond, levels = c("path","union"))
pB <- ggplot(mlt, aes(pair, nlp, colour = treat, shape = cond, group = interaction(treat, cond))) +
  geom_hline(yintercept = nlp(0.05), linetype = 2, colour = km_ref, alpha = 0.8) +
  geom_point(size = 2.6, position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = c(med = km_col[["med"]], studentized = km_col[["studentized"]]),
                      labels = c(med = "sigma-MED plug-in", studentized = "Studentized-F"), name = NULL) +
  scale_shape_manual(values = c(path = 1, union = 16), labels = km_cond_lab, name = NULL) +
  labs(x = "cluster pair", y = expression(paste(-log[10], "  selective p-value")),
       title = "scRNA-seq: the union rejects where the original path cannot") +
  theme_km() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
g <- arrangeGrob(pA, pB, ncol = 2, widths = c(1.25, 1))
ggsave_km(g, sprintf("sims/results/scrna_k%d", K), width = 13.5, height = 4.6)
cat(sprintf("Wrote sims/results/scrna_k%d.{pdf,png} + scrna_pvalues_k%d.csv\n", K, K))
