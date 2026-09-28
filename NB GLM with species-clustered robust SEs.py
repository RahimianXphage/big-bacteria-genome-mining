library(readxl)
library(MASS)
library(sandwich)
library(lmtest)

excel_path <- "D:/Projects/Gigantic/submit/Update 1/New fig/analyse.xlsx"
df <- read_excel(excel_path)

# ---------- 1. Parse taxonomy ----------
parse_taxonomy <- function(name) {
  name <- trimws(as.character(name))
  name <- sub("^[Uu]ncultured\\s+", "", name)
  parts <- strsplit(name, "\\s+")[[1]]
  if (length(parts) == 0) return(list(genus = "Unknown", species = "Unknown"))
  if (parts[1] == "Candidatus") {
    genus <- paste(parts[1:min(2, length(parts))], collapse = " ")
    rest  <- if (length(parts) > 2) parts[3:length(parts)] else character(0)
  } else {
    genus <- parts[1]
    rest  <- if (length(parts) > 1) parts[2:length(parts)] else character(0)
  }
  if (length(rest) > 0 &&
      !tolower(sub("\\.$", "", rest[1])) %in% c("sp", "isolate", "strain")) {
    species <- paste(genus, rest[1])
  } else {
    species <- name
  }
  list(genus = genus, species = species)
}

genus_vec   <- character(nrow(df))
species_vec <- character(nrow(df))
for (i in seq_len(nrow(df))) {
  p <- parse_taxonomy(df[["Bacteria name"]][i])
  genus_vec[i]   <- p$genus
  species_vec[i] <- p$species
}
df$Genus   <- genus_vec
df$Species <- species_vec

cat("Genus distribution BEFORE merge:\n")
print(table(df$Genus))

# ---------- 2. Merge Candidatus Thiomargarita into Thiomargarita ----------
# NOTE: This is done ONLY for the statistical model.
# Elsewhere in the manuscript, "Candidatus Thiomargarita" and "Thiomargarita"
# are reported separately, consistent with the Supplementary Material.
df$Genus[df$Genus == "Candidatus Thiomargarita"] <- "Thiomargarita"

cat("\nGenus distribution AFTER merge:\n")
print(table(df$Genus))

# ---------- 3. Prepare predictors ----------
df <- df[!is.na(df[["Size (Mb)"]]) & !is.na(df[["CheckM completeness"]]), ]
df <- df[df[["Size (Mb)"]] > 0, ]
df$log_size        <- log(df[["Size (Mb)"]])
df$Completeness_10 <- (df[["CheckM completeness"]] -
                         mean(df[["CheckM completeness"]])) / 10

# Reference genus
reference_genus <- "Beggiatoa"
df$Genus <- relevel(factor(df$Genus), ref = reference_genus)

cat("\nReference genus:", reference_genus, "\n")
cat("Final genus levels:\n")
print(levels(df$Genus))
cat("\nSpecies replication (top 10):\n")
print(head(sort(table(df$Species), decreasing = TRUE), 10))

# ---------- 4. Features ----------
features <- c("ARG count", "BGC count", "CRISPR", "Spacer", "Cas",
              "Other defense and anti-defense modules")

# ---------- 5. Fit NB GLM + species-clustered robust SEs ----------
results   <- list()
diag_rows <- list()

for (f in features) {
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("NB GLM (cluster-robust by Species) --", f, "\n")
  cat("  Reference:", reference_genus, "\n")
  cat(strrep("=", 70), "\n", sep = "")

  df$y <- as.integer(df[[f]])

  m <- glm.nb(y ~ Genus + Completeness_10 + offset(log_size), data = df)

  vcov_cl <- vcovCL(m, cluster = ~ Species, type = "HC1")
  co      <- coeftest(m, vcov. = vcov_cl)

  res_df <- data.frame(
    Feature  = f,
    Term     = rownames(co),
    Coef     = co[, 1],
    SE       = co[, 2],
    z        = co[, 3],
    p        = co[, 4],
    IRR      = exp(co[, 1]),
    IRR_lo   = exp(co[, 1] - 1.96 * co[, 2]),
    IRR_hi   = exp(co[, 1] + 1.96 * co[, 2]),
    row.names = NULL
  )
  results[[f]] <- res_df
  print(res_df, digits = 4)

  cat("\n  theta (NB dispersion, larger = less overdispersed):",
      round(m$theta, 4), "\n")
  cat("  n =", nrow(df),
      "| n species clusters =", length(unique(df$Species)), "\n")
  cat("  AIC =", round(AIC(m), 2),
      "| LogLik =", round(as.numeric(logLik(m)), 2), "\n")

  diag_rows[[f]] <- data.frame(
    Feature   = f,
    n         = nrow(df),
    n_species = length(unique(df$Species)),
    theta     = round(m$theta, 4),
    AIC       = round(AIC(m), 2),
    LogLik    = round(as.numeric(logLik(m)), 2)
  )
}

# ---------- 6. Save ----------
combined      <- do.call(rbind, results)
diag_combined <- do.call(rbind, diag_rows)

out_dir <- "D:/Projects/Gigantic/submit/Update 1/New fig"
write.csv(combined,       file.path(out_dir, "NB_GLM_R_results.csv"),
          row.names = FALSE)
write.csv(diag_combined,  file.path(out_dir, "NB_GLM_R_diagnostics.csv"),
          row.names = FALSE)

cat("\nSaved:\n",
    file.path(out_dir, "NB_GLM_R_results.csv"), "\n",
    file.path(out_dir, "NB_GLM_R_diagnostics.csv"), "\n", sep = "")
