# ============================================================
# Lifestyle analysis: ARG, BGC, Defense, CRISPR, Spacers, Cas
# MASS::glm.nb + sandwich::vcovCL, species-clustered SEs
# IRR = Symbiotic vs Free-living
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(MASS)
  library(sandwich)
  library(lmtest)
  library(dplyr)
  library(tidyr)
})

# ---------- 0. Paths ----------
base_dir <- "D:/Projects/Gigantic/submit/Update 1/Ecothesis"

bgc_path    <- file.path(base_dir, "BGC.xlsx")
bac_path    <- file.path(base_dir, "bac data.xlsx")
arg_path    <- file.path(base_dir, "ARG_counts_per_bacteria.csv")
crispr_path <- file.path(base_dir, "CRISPR.xlsx")
matrix_path <- file.path(base_dir, "Matrix.csv")

for (p in c(bgc_path, bac_path, arg_path, crispr_path, matrix_path)) {
  if (!file.exists(p)) stop("File not found: ", p)
}
cat("All input files found.\n\n")

# ---------- 1. Read ----------
bac    <- read_excel(bac_path)
bgc    <- read_excel(bgc_path)
arg    <- read.csv(arg_path,    stringsAsFactors = FALSE, check.names = FALSE)
crispr <- read_excel(crispr_path)
mat    <- read.csv(matrix_path, stringsAsFactors = FALSE, check.names = FALSE)

cat("bac rows:   ", nrow(bac),    "\n")
cat("bgc rows:   ", nrow(bgc),    "\n")
cat("arg rows:   ", nrow(arg),    "\n")
cat("arg columns:", paste(names(arg), collapse = " | "), "\n")
cat("crispr rows:", nrow(crispr), "\n")
cat("crispr columns:", paste(names(crispr), collapse = " | "), "\n")
cat("matrix rows:", nrow(mat),    "\n")
cat("matrix first column:", names(mat)[1], "\n\n")

# ---------- 2. Normalize names ----------
bac[["Bacteria name"]] <- trimws(as.character(bac[["Bacteria name"]]))
bgc[["Bacteria name"]] <- trimws(as.character(bgc[["Bacteria name"]]))
arg[["Bacteria name"]] <- trimws(as.character(arg[["Bacteria name"]]))
crispr[["Bac name"]]   <- trimws(as.character(crispr[["Bac name"]]))
mat[["Bac name"]]      <- trimws(as.character(mat[["Bac name"]]))

check_names <- function(vec, source) {
  miss <- setdiff(vec, bac[["Bacteria name"]])
  if (length(miss) > 0) {
    cat("\n[WARNING]", source, "- names not found in bac data:\n")
    print(miss)
  } else {
    cat("[OK]", source, "- all names matched bac data\n")
  }
}
check_names(bgc[["Bacteria name"]], "BGC.xlsx")
check_names(arg[["Bacteria name"]], "ARG csv")
check_names(crispr[["Bac name"]],   "CRISPR.xlsx")
check_names(mat[["Bac name"]],      "Matrix.csv")
cat("\n")

# ---------- 3. BGC total ----------
bgc <- bgc[!is.na(bgc$Type) & bgc$Type != "" & bgc$Type != "-" &
           !grepl("Nothing found", bgc$Region, ignore.case = TRUE), ]
bgc_total <- bgc %>% count(`Bacteria name`, name = "Total_BGC")
cat("BGC total table:", nrow(bgc_total), "genomes\n")

# ---------- 4. CRISPR table ----------
crispr_df <- data.frame(
  "Bacteria name" = crispr[["Bac name"]],
  CRISPR_arrays   = as.numeric(crispr[["CRISPR"]]),
  Spacers         = as.numeric(crispr[["Spacer"]]),
  Cas_genes       = as.numeric(crispr[["Cas"]]),
  check.names     = FALSE
)
cat("CRISPR table:", nrow(crispr_df), "genomes\n")

# ---------- 5. Defense module count ----------
mat_numeric <- mat[, sapply(mat, is.numeric), drop = FALSE]
mat$Defense_modules <- rowSums(mat_numeric, na.rm = TRUE)
defense_df <- data.frame(
  "Bacteria name" = mat[["Bac name"]],
  Defense_modules = mat$Defense_modules,
  check.names     = FALSE
)
cat("Defense table:", nrow(defense_df), "genomes\n")
cat("Defense modules range:",
    min(defense_df$Defense_modules), "-",
    max(defense_df$Defense_modules), "\n\n")

# ---------- 6. Merge all ----------
df <- bac %>%
  left_join(bgc_total,                              by = "Bacteria name") %>%
  left_join(arg[, c("Bacteria name", "ARG_count")], by = "Bacteria name") %>%
  left_join(crispr_df,                              by = "Bacteria name") %>%
  left_join(defense_df,                             by = "Bacteria name")

response_vars <- c("Total_BGC", "ARG_count", "Defense_modules",
                   "CRISPR_arrays", "Spacers", "Cas_genes")
for (v in response_vars) {
  df[[v]][is.na(df[[v]])] <- 0
}
cat("Merged:", nrow(df), "rows\n")

# ---------- 7. Lifestyle normalization ----------
df$Lifestyle <- trimws(tolower(as.character(df$Lifestyle)))
df$Lifestyle <- ifelse(
  df$Lifestyle %in% c("free-living", "free living", "freeliving", "free"),
  "Free-living",
  ifelse(
    df$Lifestyle %in% c("symbiotic", "symbiont", "symbiosis", "sym"),
    "Symbiotic",
    NA
  )
)
df$Lifestyle <- factor(df$Lifestyle, levels = c("Free-living", "Symbiotic"))
cat("\nLifestyle levels:\n")
print(table(df$Lifestyle, useNA = "ifany"))

# ---------- 8. Species cluster ----------
assign_cluster <- function(name) {
  p <- strsplit(name, " ")[[1]]
  if (length(p) >= 2 && p[1] == "Candidatus") {
    if (length(p) >= 3 &&
        !(p[3] %in% c("sp.", "isolate", "strain")) &&
        !grepl("^sp\\.", p[3])) return(paste(p[1], p[2], p[3]))
    return(name)
  }
  if (length(p) >= 2 &&
      !(p[2] %in% c("sp.", "isolate", "strain")) &&
      !grepl("^sp\\.", p[2])) return(paste(p[1], p[2]))
  name
}
df$species_cluster <- vapply(df[["Bacteria name"]], assign_cluster, character(1))
cat("Unique species clusters:", length(unique(df$species_cluster)), "\n")

# ---------- 9. Model variables ----------
df$completeness_c <- as.numeric(scale(df[["CheckM completeness"]],
                                      center = TRUE, scale = FALSE)) / 10
df$log_size <- log(df[["Size (Mb)"]])

# ---------- 10. Fit NB GLM ----------
fit_nb_lifestyle <- function(data, response, label = response) {
  d <- data[!is.na(data$Lifestyle), ]
  d$y <- d[[response]]

  if (sum(d$y) == 0) {
    cat("  skip", label, "(all zero)\n"); return(NULL)
  }
  if (length(unique(d$Lifestyle)) < 2) {
    cat("  skip", label, "(one lifestyle level)\n"); return(NULL)
  }

  free_pos <- sum(d$y[d$Lifestyle == "Free-living"] > 0)
  symb_pos <- sum(d$y[d$Lifestyle == "Symbiotic"] > 0)

  if (free_pos == 0 || symb_pos == 0) {
    cat("  separation for", label, "- Fisher exact\n")
    d$present <- as.integer(d$y > 0)
    ft <- tryCatch(fisher.test(table(d$Lifestyle, d$present)),
                   error = function(e) NULL)
    p_fisher <- if (!is.null(ft)) ft$p.value else NA
    return(data.frame(
      Feature    = label,
      IRR        = NA,
      CI_low     = NA,
      CI_high    = NA,
      p_value    = p_fisher,
      theta      = NA,
      n_free     = sum(d$Lifestyle == "Free-living"),
      n_symb     = sum(d$Lifestyle == "Symbiotic"),
      n_pos      = sum(d$y > 0),
      n_pos_free = free_pos,
      n_pos_symb = symb_pos,
      Note       = "complete separation (Fisher exact)",
      row.names  = NULL
    ))
  }

  m <- tryCatch(
    suppressWarnings(
      glm.nb(y ~ Lifestyle + completeness_c + offset(log_size), data = d)
    ),
    error = function(e) {
      cat("  glm.nb failed for", label, ":", conditionMessage(e), "\n"); NULL
    }
  )
  if (is.null(m)) return(NULL)

  vc <- tryCatch(
    vcovCL(m, cluster = d$species_cluster, type = "HC1"),
    error = function(e) {
      cat("  vcovCL failed for", label, ":", conditionMessage(e), "\n"); NULL
    }
  )
  if (is.null(vc)) return(NULL)

  ct <- coeftest(m, vcov = vc)
  rn <- grep("LifestyleSymbiotic", rownames(ct), value = TRUE)
  if (length(rn) == 0) return(NULL)

  b  <- ct[rn[1], "Estimate"]
  se <- ct[rn[1], "Std. Error"]
  p  <- ct[rn[1], "Pr(>|z|)"]
  unreliable <- is.na(se) || is.nan(se) || se <= 0

  data.frame(
    Feature    = label,
    IRR        = ifelse(unreliable, NA, exp(b)),
    CI_low     = ifelse(unreliable, NA, exp(b - 1.96 * se)),
    CI_high    = ifelse(unreliable, NA, exp(b + 1.96 * se)),
    p_value    = ifelse(unreliable, NA, p),
    theta      = m$theta,
    n_free     = sum(d$Lifestyle == "Free-living"),
    n_symb     = sum(d$Lifestyle == "Symbiotic"),
    n_pos      = sum(d$y > 0),
    n_pos_free = free_pos,
    n_pos_symb = symb_pos,
    Note       = ifelse(unreliable, "non-convergent", ""),
    row.names  = NULL
  )
}

# ---------- 11. Run all six responses ----------
labels <- c(
  Total_BGC       = "BGCs",
  ARG_count       = "ARGs",
  Defense_modules = "Defense modules",
  CRISPR_arrays   = "CRISPR arrays",
  Spacers         = "Spacers",
  Cas_genes       = "Cas genes"
)

cat("\n===== Fitting lifestyle models =====\n")
results_list <- lapply(names(labels), function(v) {
  cat(" fitting:", labels[v], "\n")
  fit_nb_lifestyle(df, v, label = labels[v])
})
results <- do.call(rbind, Filter(Negate(is.null), results_list))

if (!is.null(results) && nrow(results) > 0) {
  ok <- !is.na(results$p_value)
  results$FDR <- NA
  if (any(ok)) {
    results$FDR[ok] <- p.adjust(results$p_value[ok], method = "BH")
  }
  results <- results[order(results$p_value, na.last = TRUE), ]
}

cat("\n============================================================\n")
cat("Lifestyle GLM results (IRR = Symbiotic vs Free-living)\n")
cat("============================================================\n")
print(results, row.names = FALSE)

out_results <- file.path(base_dir, "lifestyle_all_results.csv")
write.csv(results, out_results, row.names = FALSE)
cat("\nSaved:", out_results, "\n")

# ---------- 12. Summary by lifestyle ----------
summary_tbl <- df %>%
  filter(!is.na(Lifestyle)) %>%
  group_by(Lifestyle) %>%
  summarise(
    n = n(),
    across(all_of(response_vars), mean, .names = "{.col}_mean"),
    .groups = "drop"
  )
cat("\nSummary by lifestyle (means):\n")
print(as.data.frame(summary_tbl))

out_summary <- file.path(base_dir, "lifestyle_all_summary.csv")
write.csv(summary_tbl, out_summary, row.names = FALSE)
cat("\nSaved:", out_summary, "\n")
