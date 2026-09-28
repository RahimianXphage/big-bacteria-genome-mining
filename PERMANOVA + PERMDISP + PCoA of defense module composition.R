# =============================================================================
# PERMANOVA + PERMDISP + PCoA of defense module composition
# Reference implementation: vegan::adonis2 and vegan::betadisper + permutest
#
# Input : Matrix.csv  (first column = genome names, remaining columns =
#                      binary presence/absence of defense modules)
# Output: Defense_summary.txt
#         Defense_stats.csv
#         Defense_PERMDISP_distances.csv
#         Defense_PERMDISP_pergenus.csv
#         Defense_PERMDISP_pairwise.csv   (only if overall p < 0.05)
# =============================================================================

suppressPackageStartupMessages(library(vegan))

# ---------- User settings ----------
FILE_PATH      <- "Matrix.csv"
PERMUTATIONS   <- 9999
SEED           <- 42
OUT_PREFIX     <- "Defense"
MIN_GROUP_N    <- 2
DEGENERATE_TOL <- 1e-8

set.seed(SEED)

# =============================================================================
# 1. Read and validate the input matrix
# =============================================================================
if (!file.exists(FILE_PATH)) {
  msg <- sprintf("Cannot find '%s' in the current working directory (%s).", FILE_PATH, getwd())
  stop(msg)
}

df <- read.csv(FILE_PATH, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(df) < 2) {
  stop("Matrix.csv must contain at least one genome column and one module column.")
}

colnames(df)[1] <- "Bac_name"
module_cols <- setdiff(colnames(df), "Bac_name")

for (c in module_cols) {
  df[[c]] <- as.integer(df[[c]])
}

cat(sprintf("Loaded %d genomes and %d defense modules.\n", nrow(df), length(module_cols)))

# =============================================================================
# 2. Drop genomes with no defense modules detected
# =============================================================================
X_numeric <- as.matrix(df[, module_cols])
all_zero  <- rowSums(X_numeric) == 0

if (any(all_zero)) {
  n_drop <- sum(all_zero)
  cat(sprintf("\nRemoving %d genome(s) with no defense modules detected:\n", n_drop))
  for (nm in df$Bac_name[all_zero]) {
    cat("  - ", nm, "\n", sep = "")
  }
  df <- df[!all_zero, ]
}
cat(sprintf("\nRetained %d genomes for downstream analysis.\n", nrow(df)))

# =============================================================================
# 3. Genus assignment
# =============================================================================
GENUS_PREFIXES <- c(
  "Candidatus Epulonipiscioides",
  "Candidatus Epulonipiscium",
  "Candidatus Parepulonipiscium",
  "Candidatus Thiomargarita",
  "Epulopiscium",
  "Beggiatoa",
  "Thioploca",
  "Thiomargarita"
)

assign_genus <- function(name) {
  name <- trimws(as.character(name))
  name <- sub("^[Uu]ncultured\\s+", "", name)
  for (p in GENUS_PREFIXES) {
    if (startsWith(name, p)) return(p)
  }
  m <- regmatches(name, regexpr("Candidatus\\s+\\w+", name))
  if (length(m) > 0) return(m)
  "Other"
}

df$Genus <- vapply(df$Bac_name, assign_genus, character(1))
df$Genus <- factor(df$Genus)

genus_tab <- table(df$Genus)
small_genera <- names(genus_tab)[genus_tab < MIN_GROUP_N]

if (length(small_genera) > 0) {
  cat(sprintf("\nDropping %d genus/genera with fewer than %d genome(s):\n",
              length(small_genera), MIN_GROUP_N))
  for (g in small_genera) {
    cat(sprintf("  - %s (n = %d)\n", g, genus_tab[g]))
  }
  df <- df[!(df$Genus %in% small_genera), ]
  df$Genus <- droplevels(df$Genus)
}

cat("\nGenus sample sizes:\n")
print(table(df$Genus))
cat(sprintf("\nFinal dataset: %d genomes in %d genera.\n", nrow(df), nlevels(df$Genus)))

# =============================================================================
# 4. Jaccard distance matrix
# =============================================================================
X <- as.matrix(df[, module_cols])
X[X > 0] <- 1L
D <- vegdist(X, method = "jaccard", binary = TRUE)

if (any(is.na(D))) {
  stop("Distance matrix still contains NA. Check the input matrix.")
}
if (any(D < 0)) {
  stop("Distance matrix contains negative values.")
}

# =============================================================================
# 5. PERMANOVA
# =============================================================================
cat("\n", strrep("=", 72), "\n", sep = "")
cat("PERMANOVA (vegan::adonis2)\n")
cat(strrep("=", 72), "\n", sep = "")

permanova_res <- adonis2(D ~ Genus, data = df,
                         permutations = PERMUTATIONS,
                         method = "jaccard",
                         by = "terms")
print(permanova_res)

rn <- rownames(permanova_res)
effect_row <- setdiff(rn, c("Residual", "Total"))
if (length(effect_row) == 0) {
  stop("Could not identify the effect row in the adonis2 output.")
}
effect_row <- effect_row[1]

perm_F   <- as.numeric(permanova_res[effect_row, "F"])
perm_R2  <- as.numeric(permanova_res[effect_row, "R2"])
perm_p   <- as.numeric(permanova_res[effect_row, "Pr(>F)"])
perm_df1 <- as.numeric(permanova_res[effect_row, "Df"])
perm_df2 <- as.numeric(permanova_res["Residual", "Df"])

cat(sprintf("\nExtracted PERMANOVA values (row '%s'):\n", effect_row))
cat(sprintf("  R2      = %.4f\n", perm_R2))
cat(sprintf("  F       = %.4f\n", perm_F))
cat(sprintf("  p-value = %.4f\n", perm_p))
cat(sprintf("  df      = %d (between), %d (within)\n", perm_df1, perm_df2))

# =============================================================================
# 6. PERMDISP
# =============================================================================
cat("\n", strrep("=", 72), "\n", sep = "")
cat("PERMDISP (vegan::betadisper, type = \"median\")\n")
cat(strrep("=", 72), "\n", sep = "")

bd <- betadisper(D, df$Genus, type = "median")
pt <- permutest(bd, permutations = PERMUTATIONS)
print(pt)

disp_F <- as.numeric(pt$tab$F[1])
disp_p <- as.numeric(pt$tab$`Pr(>F)`[1])

cat(sprintf("\nExtracted PERMDISP values:\n"))
cat(sprintf("  F       = %.4f\n", disp_F))
cat(sprintf("  p-value = %.4f\n", disp_p))

disp_df <- data.frame(
  Genus  = names(tapply(bd$distances, bd$group, mean)),
  n      = as.integer(table(bd$group)),
  mean   = as.numeric(tapply(bd$distances, bd$group, mean)),
  median = as.numeric(tapply(bd$distances, bd$group, median)),
  sd     = as.numeric(tapply(bd$distances, bd$group, sd))
)
disp_df <- disp_df[order(-disp_df$mean), ]
rownames(disp_df) <- NULL

cat("\nMean distance to spatial median per genus:\n")
print(disp_df, row.names = FALSE)

# =============================================================================
# 7. Pairwise PERMDISP
# =============================================================================
pairwise_res <- NULL
cat("\n", strrep("=", 72), "\n", sep = "")
cat("Pairwise PERMDISP (Holm-Bonferroni corrected)\n")
cat(strrep("=", 72), "\n", sep = "")

if (disp_p < 0.05) {

  genera <- levels(df$Genus)
  pairs  <- combn(genera, 2, simplify = FALSE)

  pairwise_list <- lapply(pairs, function(p) {

    g1  <- p[1]
    g2  <- p[2]
    idx <- which(df$Genus %in% c(g1, g2))

    if (length(idx) < 2 * MIN_GROUP_N) return(NULL)

    D_sub <- as.dist(as.matrix(D)[idx, idx])
    lab   <- droplevels(df$Genus[idx])

    if (nlevels(lab) < 2) return(NULL)
    if (any(is.na(D_sub))) return(NULL)

    m_sub <- as.matrix(D_sub)
    d_g1  <- m_sub[lab == g1, lab == g1]
    d_g2  <- m_sub[lab == g2, lab == g2]
    ss1   <- sum((d_g1[upper.tri(d_g1)])^2)
    ss2   <- sum((d_g2[upper.tri(d_g2)])^2)

    if (ss1 < DEGENERATE_TOL && ss2 < DEGENERATE_TOL) {
      cat(sprintf("  [skipping] %s vs %s: both groups have zero within-group dispersion\n", g1, g2))
      return(NULL)
    }

    bd_sub <- tryCatch(betadisper(D_sub, lab, type = "median"),
                       error = function(e) NULL)
    if (is.null(bd_sub)) return(NULL)

    pt_sub <- tryCatch(permutest(bd_sub, permutations = PERMUTATIONS,
                                 pairwise = FALSE),
                       error = function(e) NULL)
    if (is.null(pt_sub)) return(NULL)

    Fval <- as.numeric(pt_sub$tab$F[1])

    if (!is.finite(Fval) || Fval > 1e6) {
      cat(sprintf("  [skipping] %s vs %s: F-statistic is not finite\n", g1, g2))
      return(NULL)
    }

    data.frame(
      Genus_1     = g1,
      Genus_2     = g2,
      n_1         = sum(lab == g1),
      n_2         = sum(lab == g2),
      mean_dist_1 = mean(bd_sub$distances[lab == g1]),
      mean_dist_2 = mean(bd_sub$distances[lab == g2]),
      F           = Fval,
      p_raw       = as.numeric(pt_sub$tab$`Pr(>F)`[1])
    )
  })

  pairwise_res <- do.call(rbind, Filter(Negate(is.null), pairwise_list))

  if (!is.null(pairwise_res) && nrow(pairwise_res) > 0) {
    pairwise_res <- pairwise_res[order(pairwise_res$p_raw), ]
    m <- nrow(pairwise_res)
    p_holm <- sapply(seq_len(m), function(i) min(1, pairwise_res$p_raw[i] * (m - i + 1)))
    if (m > 1) {
      for (i in 2:m) p_holm[i] <- max(p_holm[i], p_holm[i - 1])
    }
    pairwise_res$p_holm <- p_holm
    rownames(pairwise_res) <- NULL
    print(pairwise_res, row.names = FALSE)

    write.csv(pairwise_res,
              sprintf("%s_PERMDISP_pairwise.csv", OUT_PREFIX),
              row.names = FALSE)
  } else {
    cat("No valid pairwise comparisons could be computed.\n")
  }

} else {
  cat("Overall PERMDISP is not significant; pairwise tests not run.\n")
}

# =============================================================================
# 8. Export per-genome distances
# =============================================================================
disp_out <- data.frame(
  Bac_name           = df$Bac_name,
  Genus              = df$Genus,
  distance_to_median = as.numeric(bd$distances)
)
write.csv(disp_out,
          sprintf("%s_PERMDISP_distances.csv", OUT_PREFIX),
          row.names = FALSE)

# =============================================================================
# 9. Export per-genus dispersion summary
# =============================================================================
write.csv(disp_df,
          sprintf("%s_PERMDISP_pergenus.csv", OUT_PREFIX),
          row.names = FALSE)

# =============================================================================
# 10. Export stats file for the Python figure script
# =============================================================================
stats_df <- data.frame(
  term  = c("permanova_R2", "permanova_F", "permanova_p",
            "permdisp_F",  "permdisp_p",  "n_permutations",
            "df_between",  "df_within",
            "n_genomes",   "n_genera"),
  value = c(perm_R2, perm_F, perm_p,
            disp_F,  disp_p,  PERMUTATIONS,
            perm_df1, perm_df2,
            nrow(df), nlevels(df$Genus))
)
write.csv(stats_df,
          sprintf("%s_stats.csv", OUT_PREFIX),
          row.names = FALSE)

# =============================================================================
# 11. Plain-text summary
# =============================================================================
sink(sprintf("%s_summary.txt", OUT_PREFIX))

cat(strrep("=", 72), "\n", sep = "")
cat("PERMANOVA (vegan::adonis2)\n")
cat(strrep("=", 72), "\n", sep = "")
cat(sprintf("R2      = %.4f\n", perm_R2))
cat(sprintf("F       = %.4f\n", perm_F))
cat(sprintf("p-value = %.4f  (%d permutations)\n", perm_p, PERMUTATIONS))
cat(sprintf("df      = %d (between), %d (within)\n", perm_df1, perm_df2))
cat(sprintf("Number of genomes = %d, number of genera = %d\n", nrow(df), nlevels(df$Genus)))
cat("\n")

cat(strrep("=", 72), "\n", sep = "")
cat("PERMDISP (vegan::betadisper, type = \"median\")\n")
cat(strrep("=", 72), "\n", sep = "")
cat(sprintf("F       = %.4f\n", disp_F))
cat(sprintf("p-value = %.4f  (%d permutations)\n", disp_p, PERMUTATIONS))
cat("\nMean distance to spatial median per genus:\n")

for (i in seq_len(nrow(disp_df))) {
  cat(sprintf("  %-38s n=%d  mean=%.4f  median=%.4f  sd=%.4f\n",
              disp_df$Genus[i], disp_df$n[i],
              disp_df$mean[i], disp_df$median[i], disp_df$sd[i]))
}

cat("\n")

if (disp_p < 0.05) {
  cat("Overall PERMDISP is significant.\n")
  cat("The PERMANOVA result should be interpreted with the caveat that within-genus dispersion is heterogeneous.\n")
  cat("See pairwise table for which genera differ.\n")
} else {
  cat("Overall PERMDISP is not significant.\n")
  cat("This supports interpretation of the PERMANOVA result as reflecting compositional differences.\n")
}

if (!is.null(pairwise_res) && nrow(pairwise_res) > 0) {
  cat("\n")
  cat(strrep("=", 72), "\n", sep = "")
  cat("Pairwise PERMDISP (Holm-Bonferroni corrected)\n")
  cat(strrep("=", 72), "\n", sep = "")
  for (i in seq_len(nrow(pairwise_res))) {
    cat(sprintf("  %-32s vs %-32s F=%7.3f  p_raw=%.4f  p_holm=%.4f\n",
                pairwise_res$Genus_1[i], pairwise_res$Genus_2[i],
                pairwise_res$F[i], pairwise_res$p_raw[i], pairwise_res$p_holm[i]))
  }
}

sink()

# =============================================================================
# 12. Console summary
# =============================================================================
cat("\nFiles written:\n")
cat(sprintf("  %s_summary.txt\n", OUT_PREFIX))
cat(sprintf("  %s_stats.csv\n", OUT_PREFIX))
cat(sprintf("  %s_PERMDISP_distances.csv\n", OUT_PREFIX))
cat(sprintf("  %s_PERMDISP_pergenus.csv\n", OUT_PREFIX))
if (!is.null(pairwise_res) && nrow(pairwise_res) > 0) {
  cat(sprintf("  %s_PERMDISP_pairwise.csv\n", OUT_PREFIX))
}
