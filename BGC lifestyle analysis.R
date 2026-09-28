suppressPackageStartupMessages({
  library(readxl)
  library(MASS)
  library(sandwich)
  library(lmtest)
  library(dplyr)
  library(tidyr)
})

# ---------- 0. File paths ----------
bgc_path <- "D:/Projects/Gigantic/submit/Update 1/Ecothesis/BGC.xlsx"
bac_path <- "D:/Projects/Gigantic/submit/Update 1/Ecothesis/bac data.xlsx"
out_dir  <- "D:/Projects/Gigantic/submit/Update 1/Ecothesis"

if (!file.exists(bgc_path)) stop("BGC file not found: ", bgc_path)
if (!file.exists(bac_path)) stop("Bac file not found: ", bac_path)
cat("Both input files found.\n")

# ---------- 1. Read ----------
bgc <- read_excel(bgc_path)
bac <- read_excel(bac_path)
bgc$`Bacteria name` <- trimws(as.character(bgc$`Bacteria name`))
bac$`Bacteria name` <- trimws(as.character(bac$`Bacteria name`))
cat("BGC rows:", nrow(bgc), "\n")
cat("Bac rows:", nrow(bac), "\n")

# ---------- 2. Clean BGC rows ----------
bgc <- bgc[!is.na(bgc$Type) &
             bgc$Type != "" &
             bgc$Type != "-" &
             !grepl("Nothing found", bgc$Region, ignore.case = TRUE), ]
cat("BGC rows after cleaning:", nrow(bgc), "\n")

# ---------- 3. Classify BGC types into product classes ----------
classify_bgc <- function(t) {
  if (is.na(t)) return("Other")
  t <- tolower(t)
  classes <- character(0)
  if (grepl("nrps", t))        classes <- c(classes, "NRPS")
  if (grepl("pks", t))         classes <- c(classes, "PKS")
  if (grepl("terpene", t))     classes <- c(classes, "Terpene")
  if (grepl("ranthipeptide|lanthipeptide|lassopeptide|cyanobactin|azole-containing-ripp|triceptide", t))
                               classes <- c(classes, "RiPP")
  if (grepl("phosphonate", t)) classes <- c(classes, "Phosphonate")
  if (grepl("arylpolyene", t)) classes <- c(classes, "Arylpolyene")
  if (grepl("siderophore", t)) classes <- c(classes, "Siderophore")
  if (length(classes) == 0) classes <- "Other"
  paste(sort(unique(classes)), collapse = ";")
}
bgc$Classes <- vapply(bgc$Type, classify_bgc, character(1))

# ---------- 4. Count BGCs per genome per class ----------
bgc_long <- bgc %>%
  separate_rows(Classes, sep = ";") %>%
  filter(Classes != "") %>%
  count(`Bacteria name`, Classes, name = "n")

bgc_wide <- bgc_long %>%
  pivot_wider(names_from = Classes, values_from = n, values_fill = 0)

bgc_total <- bgc %>% count(`Bacteria name`, name = "Total_BGC")
bgc_wide  <- left_join(bgc_total, bgc_wide, by = "Bacteria name")
bgc_wide[is.na(bgc_wide)] <- 0

cat("\nGenomes with BGCs:", nrow(bgc_wide), "\n")
cat("Classes:", paste(setdiff(names(bgc_wide),
    c("Bacteria name", "Total_BGC")), collapse = ", "), "\n")

# ---------- 5. Merge with metadata (numeric-only NA fill) ----------
merged <- bac %>% left_join(bgc_wide, by = "Bacteria name")

bgc_cols <- setdiff(names(bgc_wide), "Bacteria name")
for (cc in bgc_cols) {
  merged[[cc]][is.na(merged[[cc]])] <- 0
}

cat("\nMerged:", nrow(merged), "rows\n")

# ---------- 6. Normalize lifestyle labels ----------
merged$Lifestyle <- trimws(tolower(as.character(merged$Lifestyle)))
merged$Lifestyle <- ifelse(
  merged$Lifestyle %in% c("free-living", "free living", "freeliving", "free"),
  "Free-living",
  ifelse(
    merged$Lifestyle %in% c("symbiotic", "symbiont", "symbiosis", "sym"),
    "Symbiotic",
    NA
  )
)
merged$Lifestyle <- factor(merged$Lifestyle,
                           levels = c("Free-living", "Symbiotic"))

cat("\nLifestyle factor levels after normalization:\n")
print(table(merged$Lifestyle, useNA = "ifany"))

# ---------- 7. Species cluster ----------
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
merged$species_cluster <- vapply(merged$`Bacteria name`,
                                 assign_cluster, character(1))
cat("\nUnique species clusters:", length(unique(merged$species_cluster)), "\n")

# ---------- 8. Model variables ----------
merged$completeness_c <- as.numeric(scale(merged$`CheckM completeness`,
                                          center = TRUE, scale = FALSE)) / 10
merged$log_size       <- log(merged$`Size (Mb)`)

# ---------- 9. NB GLM with species-clustered SEs ----------
fit_nb_lifestyle <- function(df, response) {
  df <- df[!is.na(df$Lifestyle), ]
  df$y <- df[[response]]

  if (sum(df$y) == 0) {
    cat("  skip", response, "(all zero)\n"); return(NULL)
  }
  if (length(unique(df$Lifestyle)) < 2) {
    cat("  skip", response, "(only one lifestyle level)\n"); return(NULL)
  }

  # Check for complete separation: all positives in one lifestyle
  free_pos <- sum(df$y[df$Lifestyle == "Free-living"] > 0)
  symb_pos <- sum(df$y[df$Lifestyle == "Symbiotic"] > 0)

  if (free_pos == 0 || symb_pos == 0) {
    cat("  separation detected for", response,
        "- running Fisher exact on presence/absence\n")

    df$present <- as.integer(df$y > 0)
    tab <- table(df$Lifestyle, df$present)

    ft <- tryCatch(fisher.test(tab), error = function(e) NULL)
    p_fisher <- if (!is.null(ft)) ft$p.value else NA

    return(data.frame(
      Feature   = response,
      IRR       = NA,
      CI_low    = NA,
      CI_high   = NA,
      p_value   = p_fisher,
      theta     = NA,
      n_free    = sum(df$Lifestyle == "Free-living"),
      n_symb    = sum(df$Lifestyle == "Symbiotic"),
      n_pos     = sum(df$y > 0),
      n_pos_free = free_pos,
      n_pos_symb = symb_pos,
      Note      = "complete separation (Fisher exact)",
      row.names = NULL
    ))
  }

  # Regular NB fit
  m <- tryCatch(
    suppressWarnings(
      glm.nb(y ~ Lifestyle + completeness_c + offset(log_size), data = df)
    ),
    error = function(e) {
      cat("  glm.nb failed for", response, ":", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(m)) return(NULL)

  vc <- tryCatch(
    vcovCL(m, cluster = df$species_cluster, type = "HC1"),
    error = function(e) {
      cat("  vcovCL failed for", response, ":", conditionMessage(e), "\n")
      NULL
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
    Feature   = response,
    IRR       = ifelse(unreliable, NA, exp(b)),
    CI_low    = ifelse(unreliable, NA, exp(b - 1.96 * se)),
    CI_high   = ifelse(unreliable, NA, exp(b + 1.96 * se)),
    p_value   = ifelse(unreliable, NA, p),
    theta     = m$theta,
    n_free    = sum(df$Lifestyle == "Free-living"),
    n_symb    = sum(df$Lifestyle == "Symbiotic"),
    n_pos     = sum(df$y > 0),
    n_pos_free = free_pos,
    n_pos_symb = symb_pos,
    Note      = ifelse(unreliable, "non-convergent", ""),
    row.names = NULL
  )
}

# Primary analysis
cat("\nFitting primary model (Total_BGC)...\n")
primary <- fit_nb_lifestyle(merged, "Total_BGC")
if (!is.null(primary)) {
  primary$Analysis <- "Primary"
  primary$FDR      <- NA
}

# Secondary analyses
class_cols <- setdiff(names(bgc_wide), c("Bacteria name", "Total_BGC"))
cat("\nFitting secondary models for", length(class_cols), "classes...\n")
secondary_list <- lapply(class_cols, function(cl) {
  cat("  fitting:", cl, "\n")
  r <- fit_nb_lifestyle(merged, cl)
  if (!is.null(r)) {
    r$Analysis <- "Secondary"
    r$FDR      <- NA
  }
  r
})
secondary <- do.call(rbind, Filter(Negate(is.null), secondary_list))

# FDR correction on secondary analyses only, on valid p-values
if (!is.null(secondary) && nrow(secondary) > 0) {
  ok <- !is.na(secondary$p_value)
  if (any(ok)) {
    secondary$FDR[ok] <- p.adjust(secondary$p_value[ok], method = "BH")
  }
}

# Combine
results <- rbind(primary, secondary)
if (!is.null(results) && nrow(results) > 0) {
  results <- results[order(results$p_value, na.last = TRUE), ]
  cat("\n============================================================\n")
  cat("BGC lifestyle GLM results (IRR = Symbiotic vs Free-living)\n")
  cat("============================================================\n")
  print(results, row.names = FALSE)

  out_results <- file.path(out_dir, "BGC_lifestyle_results.csv")
  write.csv(results, out_results, row.names = FALSE)
  cat("\nSaved:", out_results, "\n")
} else {
  cat("\nNo results to save.\n")
}

# ---------- 10. Summary by lifestyle ----------
summary_tbl <- merged %>%
  filter(!is.na(Lifestyle)) %>%
  group_by(Lifestyle) %>%
  summarise(
    n = n(),
    across(c(Total_BGC, all_of(class_cols)),
           mean, .names = "{.col}_mean"),
    .groups = "drop"
  )
cat("\nSummary by lifestyle:\n")
print(as.data.frame(summary_tbl))

out_summary <- file.path(out_dir, "BGC_lifestyle_summary.csv")
write.csv(summary_tbl, out_summary, row.names = FALSE)
cat("\nSaved:", out_summary, "\n")
