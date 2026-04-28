# Paper: Gypsy, Irish Traveller, and Roma young people

# - Most scripts were designed to run on an HPC cluster
#   (Linux/SLURM) using the parallel package with mclapply.
# - To run sections independently, each can be extracted and
#   run as a standalone script.
# - The orchestrator script grt_subgroup_v0.91.R coordinates
#   several of these scripts with modified paths/captions for
#   the subgroup analysis (syear==2025 + syear==2023 grade 10/11).
#

# ============================================================
# SCRIPT: Descriptive Table, Table 1 & Poisson Models
# ============================================================
# Part A: Pooled prevalence estimates (%) with 95% CIs across
#         multiply imputed datasets using Rubin's rules
# Part B: Table 1 — Sample characteristics by ethnicity
# Part C: Unadjusted & Adjusted Poisson models (RRs)
# - Exposure: eth_8cat (categorical, ref = "Gypsy/Traveller")
# - Outcomes: 16 binary outcomes grouped by domain
# - Model 1: Unadjusted (eth_8cat only)
# - Model 2: Adjusted for gender + grade + fas
# - Cluster-robust SEs by id2c (school cluster)
# - Output: CSV files + combined Word document with all tables
# ============================================================

# ── 1. Setup ──
user_lib <- file.path(Sys.getenv("HOME"), "R_libs")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))

pkgs <- c("parallel", "dplyr", "clubSandwich", "mice", "tidyr",
          "officer", "flextable")
new_pkgs <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(new_pkgs) > 0) {
  install.packages(new_pkgs, lib = user_lib, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(parallel)
  library(dplyr)
  library(clubSandwich)
  library(mice)
  library(tidyr)
  library(officer)
  library(flextable)
})

# --- Configuration ---
N_CORES <- 3L
BASEDIR <- "/shared/home1/c.wppjw/SHRN/Gypsy RT"
OUTDIR  <- "/shared/home1/c.wppjw/SHRN/Gypsy RT/output/v0.6"
if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

IMP_PATH <- file.path(BASEDIR, "output/v0.4/GRT_2023_2025_imp_m44_parallel.rds")

# --- Outcomes ---
outcomes <- c("lifesatneg", "sdqinternalhigh", "sdqexternalhigh",
              "cannabis30d", "alcoholuse", "smokeweekly", "regvape",
              "phys7", "vigexercise", "sit",
              "dailysoft", "dailyenergy",
              "beenbulliedever", "cbeenbulliedever",
              "truantever", "excludedever")

# --- Outcome labels ---
outcome_labels <- c(
  "lifesatneg"       = "Low life satisfaction",
  "sdqexternalhigh"  = "Externalizing difficulties\u1d43",
  "sdqinternalhigh"  = "Internalizing difficulties\u1d47",
  "smokeweekly"      = "Smoke weekly",
  "regvape"          = "Vape weekly",
  "alcoholuse"       = "\u22652 alcoholic drinks a session",
  "cannabis30d"      = "Cannabis in past 30 days",
  "phys7"            = "60 minutes of MVPA per day",
  "vigexercise"      = "\u22653 days VPA per week\u1d9c",
  "sit"              = "Sedentary behaviour",
  "dailysoft"        = "Daily soft drink consumption",
  "dailyenergy"      = "Daily energy drink consumption",
  "beenbulliedever"  = "Ever been bullied",
  "cbeenbulliedever" = "Ever been cyberbullied",
  "truantever"       = "Ever truanted",
  "excludedever"     = "Ever excluded from school"
)

# --- Outcome groups ---
outcome_groups <- list(
  "Mental health and life satisfaction" = c("lifesatneg", "sdqexternalhigh", "sdqinternalhigh"),
  "Substance use"     = c("smokeweekly", "regvape", "alcoholuse", "cannabis30d"),
  "Physical activity"  = c("phys7", "vigexercise", "sit"),
  "Diet"              = c("dailysoft", "dailyenergy"),
  "Bullying"          = c("beenbulliedever", "cbeenbulliedever"),
  "School problems"   = c("truantever", "excludedever")
)

# --- Ethnicity labels ---
eth_labels <- c(
  "1" = "White",
  "2" = "White Other",
  "3" = "Gypsy/Traveller",
  "4" = "Roma",
  "5" = "Black",
  "6" = "Mixed",
  "7" = "Asian",
  "8" = "Other"
)

eth_order <- c("White", "White Other", "Gypsy/Traveller", "Roma",
               "Black", "Mixed", "Asian", "Other")

gender_labels <- c("1" = "Male", "2" = "Female", "3" = "Other/prefer not to say")
grade_labels  <- c("7" = "Year 7", "8" = "Year 8", "9" = "Year 9",
                    "10" = "Year 10", "11" = "Year 11")

# ══════════════════════════════════════════════════════════════
# 2. LOAD DATA
# ══════════════════════════════════════════════════════════════
cat("Loading imputed data...\n")
t0 <- Sys.time()
imp <- readRDS(IMP_PATH)

# Get N per ethnicity from original data
eth_tab <- table(imp$data$eth_8cat, useNA = "no")
eth_n <- as.data.frame(eth_tab)
names(eth_n) <- c("Ethnicity", "n")
eth_n$pct <- round(100 * eth_n$n / sum(eth_n$n), 1)
eth_n$n_pct <- paste0(format(eth_n$n, big.mark = ","), " (", eth_n$pct, "%)")
eth_n$Ethnicity <- factor(eth_n$Ethnicity, levels = eth_order)
eth_n <- eth_n[order(eth_n$Ethnicity), ]

n_imp <- imp$m
n_obs <- nrow(imp$data)

cat(sprintf("Loaded mids object: %d imputations, n = %d (%.1f sec)\n",
            n_imp, n_obs,
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
cat(sprintf("  mids object size: %.1f GB\n", object.size(imp) / 1e9))
cat(sprintf("  Current R memory usage: %.1f GB\n\n", sum(gc()[, 2]) / 1024))

# ══════════════════════════════════════════════════════════════
# PART A: DESCRIPTIVE TABLE — POOLED PREVALENCES (TABLE 2)
# ══════════════════════════════════════════════════════════════

compute_prevalence <- function(imp_list, outcome_var, eth_var = "eth_8cat") {
  m <- length(imp_list)

  props_list <- lapply(seq_len(m), function(i) {
    df <- imp_list[[i]]
    y <- df[[outcome_var]]
    if (is.null(y)) stop(sprintf("Variable '%s' not found in imputed dataset %d. Available columns: %s",
                                  outcome_var, i, paste(names(df), collapse = ", ")))
    stopifnot(length(y) == nrow(df))

    tapply(y, df[[eth_var]], function(x) {
      n <- sum(!is.na(x))
      p <- mean(x, na.rm = TRUE)
      p_clamped <- max(min(p, 1 - 1e-8), 1e-8)
      logit_p <- log(p_clamped / (1 - p_clamped))
      var_logit <- 1 / (n * p_clamped * (1 - p_clamped))
      c(p = p, logit_p = logit_p, var_logit = var_logit, n = n)
    })
  })

  groups <- names(props_list[[1]])

  results <- data.frame(
    Ethnicity = character(), pct = numeric(),
    lower = numeric(), upper = numeric(), n = numeric(),
    stringsAsFactors = FALSE
  )

  for (grp in groups) {
    logit_ests <- sapply(props_list, function(x) x[[grp]]["logit_p"])
    var_withins <- sapply(props_list, function(x) x[[grp]]["var_logit"])
    ns <- sapply(props_list, function(x) x[[grp]]["n"])

    q_bar <- mean(logit_ests)
    u_bar <- mean(var_withins)
    b_m <- var(logit_ests)
    total_var <- u_bar + (1 + 1/m) * b_m
    total_se <- sqrt(total_var)

    r <- (1 + 1/m) * b_m / u_bar
    df_old <- (m - 1) * (1 + 1/r)^2
    if (is.nan(df_old) || is.na(df_old)) df_old <- Inf
    t_crit <- qt(0.975, df = df_old)

    pct <- 100 * plogis(q_bar)
    lower <- 100 * plogis(q_bar - t_crit * total_se)
    upper <- 100 * plogis(q_bar + t_crit * total_se)

    results <- rbind(results, data.frame(
      Ethnicity = grp, pct = pct, lower = lower,
      upper = upper, n = round(mean(ns)),
      stringsAsFactors = FALSE
    ))
  }

  results$Outcome <- outcome_var
  results
}

# Extract datasets for prevalence computation (need eth_8cat labelled)
cat("Computing pooled prevalences...\n")
t_start <- Sys.time()

prev_imp_list <- lapply(seq_len(n_imp), function(i) {
  df <- mice::complete(imp, action = i)
  df$eth_8cat <- factor(df$eth_8cat, levels = names(eth_labels), labels = eth_labels)
  for (oc in outcomes) {
    if (is.factor(df[[oc]])) df[[oc]] <- as.integer(df[[oc]]) - 1L
  }
  df
})

all_prev <- list()
for (out in outcomes) {
  cat("  ", out, "...\n")
  all_prev[[out]] <- compute_prevalence(prev_imp_list, out)
}
prev_df <- do.call(rbind, all_prev)
rownames(prev_df) <- NULL

elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
cat(sprintf("Prevalences computed in %.1f seconds.\n\n", elapsed))

prev_df$formatted <- sprintf("%.1f (%.1f, %.1f)", prev_df$pct, prev_df$lower, prev_df$upper)

rm(prev_imp_list); gc(verbose = FALSE)

# --- Build wide table ---
rows <- list()
row_types <- character()
row_codes <- character()

for (grp_name in names(outcome_groups)) {
  rows[[length(rows) + 1]] <- c(grp_name, rep("", length(eth_order)))
  row_types <- c(row_types, "header")
  row_codes <- c(row_codes, grp_name)

  for (oc in outcome_groups[[grp_name]]) {
    vals <- c(outcome_labels[oc])
    for (eth in eth_order) {
      match_row <- prev_df[prev_df$Ethnicity == eth & prev_df$Outcome == oc, ]
      if (nrow(match_row) > 0) {
        vals <- c(vals, match_row$formatted[1])
      } else {
        vals <- c(vals, "")
      }
    }
    rows[[length(rows) + 1]] <- vals
    row_types <- c(row_types, "outcome")
    row_codes <- c(row_codes, oc)
  }
}

wide_df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
names(wide_df) <- c("Outcome", eth_order)

n_row <- c("N", eth_n$n_pct[match(eth_order, eth_n$Ethnicity)])
wide_df <- rbind(n_row, wide_df)

write.csv(wide_df, file.path(OUTDIR, "Descriptive_Table_Wide.csv"), row.names = FALSE)
cat("Saved: Descriptive_Table_Wide.csv\n")

# ══════════════════════════════════════════════════════════════
# PART B: TABLE 1 — SAMPLE CHARACTERISTICS BY ETHNICITY
# ══════════════════════════════════════════════════════════════
cat("\n====== TABLE 1 ======\n")

# Pool N per ethnicity across imputations
cat("Pooling N per ethnicity...\n")
n_mat <- matrix(0, nrow = n_imp, ncol = length(eth_order))
colnames(n_mat) <- eth_order

for (i in seq_len(n_imp)) {
  df_i <- complete(imp, i)
  df_i$eth_8cat <- factor(df_i$eth_8cat, levels = names(eth_labels),
                          labels = eth_labels)
  tab <- table(df_i$eth_8cat)
  for (eth in eth_order) {
    n_mat[i, eth] <- tab[eth]
  }
}

n_pooled     <- round(colMeans(n_mat))
total_pooled <- sum(n_pooled)
pct_pooled   <- round(100 * n_pooled / total_pooled, 1)

# Pool categorical cross-tabulations
pool_crosstab <- function(imp_obj, row_var, row_levels, row_labels,
                          col_var = "eth_8cat", col_levels = names(eth_labels),
                          col_labels = eth_labels) {
  m <- imp_obj$m
  n_rows <- length(row_levels)
  n_cols <- length(col_levels)
  count_array <- array(0, dim = c(m, n_rows, n_cols))

  for (i in seq_len(m)) {
    df_i <- complete(imp_obj, i)
    df_i[[col_var]] <- factor(df_i[[col_var]], levels = col_levels)
    df_i[[row_var]] <- factor(df_i[[row_var]], levels = row_levels)
    tab <- table(df_i[[row_var]], df_i[[col_var]])
    count_array[i, , ] <- as.matrix(tab)
  }

  mean_counts <- apply(count_array, c(2, 3), mean)
  rownames(mean_counts) <- row_labels
  colnames(mean_counts) <- col_labels
  col_totals <- colSums(mean_counts)
  pct_mat <- sweep(mean_counts, 2, col_totals, "/") * 100

  list(counts = round(mean_counts), pcts = pct_mat)
}

cat("Pooling gender x ethnicity...\n")
gender_pooled <- pool_crosstab(
  imp, row_var = "gender",
  row_levels = names(gender_labels), row_labels = gender_labels
)

cat("Pooling school year x ethnicity...\n")
grade_pooled <- pool_crosstab(
  imp, row_var = "grade",
  row_levels = names(grade_labels), row_labels = grade_labels
)

# Pool continuous variables
cat("Pooling continuous variables...\n")

pool_continuous <- function(imp_obj, var_name, col_var = "eth_8cat",
                            col_levels = names(eth_labels),
                            col_labels = eth_labels) {
  m <- imp_obj$m
  n_cols <- length(col_levels)
  means_mat <- matrix(0, nrow = m, ncol = n_cols)
  vars_mat  <- matrix(0, nrow = m, ncol = n_cols)
  ns_mat    <- matrix(0, nrow = m, ncol = n_cols)
  colnames(means_mat) <- col_labels
  colnames(vars_mat)  <- col_labels
  colnames(ns_mat)    <- col_labels

  for (i in seq_len(m)) {
    df_i <- complete(imp_obj, i)
    df_i[[col_var]] <- factor(df_i[[col_var]], levels = col_levels,
                              labels = col_labels)
    for (j in seq_along(col_labels)) {
      eth <- col_labels[j]
      vals <- df_i[[var_name]][df_i[[col_var]] == eth]
      means_mat[i, j] <- mean(vals, na.rm = TRUE)
      vars_mat[i, j]  <- var(vals, na.rm = TRUE)
      ns_mat[i, j]    <- sum(!is.na(vals))
    }
  }

  q_bar <- colMeans(means_mat)
  avg_var <- colMeans(vars_mat)
  avg_sd  <- sqrt(avg_var)

  data.frame(
    Ethnicity = col_labels,
    Mean      = round(q_bar, 2),
    SD        = round(avg_sd, 2),
    stringsAsFactors = FALSE
  )
}

fas_pooled     <- pool_continuous(imp, "fas")
famsupp_pooled <- pool_continuous(imp, "famsupp")
peer_pooled    <- pool_continuous(imp, "peer")
teacher_pooled <- pool_continuous(imp, "teacher")

# Build Table 1 data frame
cat("Building Table 1...\n")

format_n_pct_vec <- function(counts_row, pcts_row) {
  sapply(seq_along(counts_row), function(j) {
    pct_val <- round(pcts_row[j], 1)
    pct_str <- sub("\\.0$", "", formatC(pct_val, format = "f", digits = 1))
    sprintf("%s (%s%%)", format(counts_row[j], big.mark = ",", trim = TRUE),
            pct_str)
  })
}

format_mean_sd <- function(pooled_df) {
  sapply(seq_len(nrow(pooled_df)), function(j) {
    sprintf("%.2f (%.2f)", pooled_df$Mean[j], pooled_df$SD[j])
  })
}

t1_rows <- list()

t1_rows[[1]] <- c("N", "",
                sapply(eth_order, function(e) format(n_pooled[e], big.mark = ",", trim = TRUE)))

t1_rows[[2]] <- c("Gender, n (%)", "Male",
                format_n_pct_vec(gender_pooled$counts["Male", ],
                                 gender_pooled$pcts["Male", ]))
t1_rows[[3]] <- c("", "Female",
                format_n_pct_vec(gender_pooled$counts["Female", ],
                                 gender_pooled$pcts["Female", ]))
t1_rows[[4]] <- c("", "Other/prefer not to say",
                format_n_pct_vec(gender_pooled$counts["Other/prefer not to say", ],
                                 gender_pooled$pcts["Other/prefer not to say", ]))

t1_rows[[5]] <- c("School year, n (%)", "Year 7",
                format_n_pct_vec(grade_pooled$counts["Year 7", ],
                                 grade_pooled$pcts["Year 7", ]))
t1_rows[[6]] <- c("", "Year 8",
                format_n_pct_vec(grade_pooled$counts["Year 8", ],
                                 grade_pooled$pcts["Year 8", ]))
t1_rows[[7]] <- c("", "Year 9",
                format_n_pct_vec(grade_pooled$counts["Year 9", ],
                                 grade_pooled$pcts["Year 9", ]))
t1_rows[[8]] <- c("", "Year 10",
                format_n_pct_vec(grade_pooled$counts["Year 10", ],
                                 grade_pooled$pcts["Year 10", ]))
t1_rows[[9]] <- c("", "Year 11",
                format_n_pct_vec(grade_pooled$counts["Year 11", ],
                                 grade_pooled$pcts["Year 11", ]))

t1_rows[[10]] <- c("Family Affluence Scale (FAS), Mean (SD)",
                "Family Affluence Scale (FAS), Mean (SD)",
                format_mean_sd(fas_pooled))
t1_rows[[11]] <- c("Family support, Mean (SD)",
                "Family support, Mean (SD)",
                format_mean_sd(famsupp_pooled))
t1_rows[[12]] <- c("Peer support, Mean (SD)",
                "Peer support, Mean (SD)",
                format_mean_sd(peer_pooled))
t1_rows[[13]] <- c("Teacher support, Mean (SD)",
                "Teacher support, Mean (SD)",
                format_mean_sd(teacher_pooled))

tab1_df <- as.data.frame(do.call(rbind, t1_rows), stringsAsFactors = FALSE)
names(tab1_df) <- c("Variable", "Level", eth_order)

write.csv(tab1_df, file.path(OUTDIR, "Table1_Sample_Characteristics.csv"),
          row.names = FALSE)
cat("Saved: Table1_Sample_Characteristics.csv\n")


# ══════════════════════════════════════════════════════════════
# PREPARE DATASETS FOR MODELS
# ══════════════════════════════════════════════════════════════

keep_cols <- unique(c("eth_8cat", "gender", "grade", "fas", "id2c", outcomes))

prepare_df <- function(df) {
  df <- df[, keep_cols, drop = FALSE]
  df$eth_8cat <- factor(df$eth_8cat, levels = names(eth_labels),
                        labels = eth_labels)
  df$eth_8cat <- relevel(df$eth_8cat, ref = "Gypsy/Traveller")
  df$gender   <- as.factor(df$gender)
  df$grade    <- as.factor(df$grade)
  df$fas      <- as.numeric(df$fas)
  df$id2c     <- as.integer(as.factor(df$id2c))
  for (oc in outcomes) {
    if (is.factor(df[[oc]])) df[[oc]] <- as.integer(df[[oc]]) - 1L
  }
  df
}

cat("\nExtracting all imputed datasets...\n")
imp_list <- lapply(seq_len(n_imp), function(i) {
  prepare_df(mice::complete(imp, action = i))
})
cat(sprintf("  Extracted %d datasets (%.1f GB total)\n",
            n_imp, object.size(imp_list) / 1e9))

# Free the mids object
rm(imp); gc(verbose = FALSE)
cat(sprintf("  Memory after freeing mids object: %.1f GB\n\n", sum(gc()[, 2]) / 1024))

# ── Diagnostics ──
cat("=== SYSTEM DIAGNOSTICS ===\n")
cat(sprintf("  R version: %s\n", R.version.string))
cat(sprintf("  Platform: %s\n", R.version$platform))
cat(sprintf("  N cores requested: %d\n", N_CORES))
cat(sprintf("  N cores detected: %d\n", parallel::detectCores()))
mem_info <- tryCatch(system("free -h", intern = TRUE), error = function(e) NULL)
if (!is.null(mem_info)) {
  cat("  System memory:\n")
  for (line in mem_info) cat("    ", line, "\n")
}
slurm_mem <- Sys.getenv("SLURM_MEM_PER_NODE", unset = NA)
if (!is.na(slurm_mem)) cat(sprintf("  SLURM allocated memory: %s MB\n", slurm_mem))
cat(sprintf("  R process RSS: %.1f GB\n", sum(gc()[, 2]) / 1024))
cat(sprintf("  imp_list object size: %.2f GB\n", object.size(imp_list) / 1e9))
cat(sprintf("  Single dataset size: %.1f MB\n", object.size(imp_list[[1]]) / 1e6))
cat(sprintf("  Estimated peak memory (parent + %d forks): %.1f GB\n",
            N_CORES,
            (sum(gc()[, 2]) / 1024) * (1 + N_CORES)))
cat("===========================\n\n")

# ══════════════════════════════════════════════════════════════
# PART C: POISSON REGRESSION — UNADJUSTED & ADJUSTED RRs
# (TABLE S1 uses the unadjusted results)
# ══════════════════════════════════════════════════════════════

fit_model <- function(df, outcome_var, formula_rhs) {
  tryCatch({
    fml <- as.formula(paste0(outcome_var, " ~ ", formula_rhs))
    mod <- glm(fml, data = df, family = poisson(link = "log"))
    cf  <- coef(mod)
    vc  <- clubSandwich::vcovCR(mod, cluster = df$id2c, type = "CR1S")
    rm(mod)
    return(list(coef = cf, vcov = vc, error = NULL))
  }, error = function(e) return(list(error = e$message)))
}

pool_rubin <- function(fit_list) {
  results_only <- Filter(function(x) !is.null(x) && is.null(x$error), fit_list)
  m <- length(results_only)
  if (m == 0) stop("All models failed (check for OOM kills or model errors).")
  if (m < 3) warning(sprintf("Only %d/%d imputations succeeded — estimates unreliable.", m, length(fit_list)))

  q_m <- do.call(rbind, lapply(results_only, `[[`, "coef"))
  u_m <- lapply(results_only, `[[`, "vcov")

  q_bar <- colMeans(q_m)
  u_bar <- Reduce("+", u_m) / m
  b_m   <- cov(q_m)
  total_vcov <- u_bar + (1 + 1/m) * b_m

  se <- sqrt(diag(total_vcov))

  r <- (1 + 1/m) * diag(b_m) / diag(u_bar)
  v_old <- (m - 1) * (1 + 1/r)^2
  v_old[is.nan(v_old) | is.na(v_old)] <- Inf

  # Monte Carlo error diagnostics
  mc_error_est <- sqrt(diag(b_m) / m)
  mc_error_se  <- se / sqrt(2 * (m - 1))
  total_var    <- diag(total_vcov)
  fmi          <- (diag(b_m) + diag(b_m) / m) / total_var
  fmi[is.nan(fmi) | is.na(fmi)] <- 0

  mc_ratio     <- mc_error_est / se
  mc_ratio[is.nan(mc_ratio) | is.na(mc_ratio)] <- 0
  mc_adequate  <- ifelse(mc_ratio < 0.10, "OK", "INCREASE m")

  data.frame(
    term = names(q_bar),
    est  = q_bar,
    se   = se,
    RR   = exp(q_bar),
    L    = exp(q_bar - 1.96 * se),
    H    = exp(q_bar + 1.96 * se),
    p    = 2 * pt(-abs(q_bar / se), df = v_old),
    fmi  = round(fmi, 3),
    mc_error_est = round(mc_error_est, 5),
    mc_error_se  = round(mc_error_se, 5),
    mc_pct_of_se = round(100 * mc_ratio, 1),
    mc_adequate  = mc_adequate,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

run_models <- function(imp_list, outcomes, formula_rhs, model_label) {
  cat(sprintf("\n====== %s MODELS ======\n", toupper(model_label)))
  cat(sprintf("Formula: outcome ~ %s\n", formula_rhs))
  cat(sprintf("Reference category: Gypsy/Traveller\n"))
  cat(sprintf("Using %d cores for parallel processing\n\n", N_CORES))

  all_res <- list()
  t_models <- Sys.time()
  for (idx in seq_along(outcomes)) {
    out <- outcomes[idx]
    t_out <- Sys.time()
    cat(sprintf("  [%d/%d] Model: %s ...", idx, length(outcomes), out))
    fits <- mclapply(imp_list, fit_model, outcome_var = out,
                     formula_rhs = formula_rhs,
                     mc.cores = N_CORES, mc.preschedule = FALSE)
    n_null <- sum(sapply(fits, is.null))
    n_err  <- sum(sapply(fits, function(x) !is.null(x$error)))
    n_ok   <- length(fits) - n_null - n_err
    elapsed_out <- as.numeric(difftime(Sys.time(), t_out, units = "secs"))
    elapsed_total <- as.numeric(difftime(Sys.time(), t_models, units = "secs"))
    avg_per <- elapsed_total / idx
    remaining <- avg_per * (length(outcomes) - idx)
    cat(sprintf(" %d/%d converged (%.1fs, ~%.0fs remaining)\n",
                n_ok, length(imp_list), elapsed_out, remaining))
    if (n_null > 0) {
      cat(sprintf("  *** WARNING: %d workers returned NULL (likely OOM killed) ***\n", n_null))
      cat(sprintf("  *** R memory: %.1f GB | Consider increasing --mem ***\n", sum(gc()[, 2]) / 1024))
    }
    if (n_err > 0) {
      errs <- sapply(Filter(function(x) !is.null(x$error), fits), `[[`, "error")
      cat(sprintf("  *** %d model errors: %s ***\n", n_err, paste(unique(errs), collapse = "; ")))
    }
    all_res[[out]] <- pool_rubin(fits) %>% mutate(Outcome = out)
    rm(fits); gc(verbose = FALSE)
  }
  cat(sprintf("  %s models completed in %.1f seconds.\n\n",
              model_label, as.numeric(difftime(Sys.time(), t_models, units = "secs"))))

  results <- bind_rows(all_res)

  eth_res <- results %>%
    filter(grepl("eth_8cat", term)) %>%
    mutate(
      term  = gsub("eth_8cat", "", term),
      sig   = ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))),
      Model = model_label
    )

  cat(sprintf("\n-- %s: Ethnicity RRs (ref = Gypsy/Traveller) --\n", model_label))
  for (out in outcomes) {
    cat(sprintf("\n--- %s ---\n", out))
    d <- eth_res %>% filter(Outcome == out)
    for (i in 1:nrow(d)) {
      cat(sprintf("  %-20s  RR = %5.2f (%4.2f\u2013%5.2f)  p = %.4f %s\n",
                  d$term[i], d$RR[i], d$L[i], d$H[i], d$p[i], d$sig[i]))
    }
  }

  # Monte Carlo error summary
  cat(sprintf("\n-- %s: Monte Carlo Error Diagnostics --\n", model_label))
  cat("  (MC error < 10%% of SE = adequate number of imputations)\n\n")
  mc_flags <- eth_res %>% filter(mc_adequate == "INCREASE m")
  if (nrow(mc_flags) == 0) {
    cat(sprintf("  ALL estimates have adequate MC error (< 10%% of SE). m = %d is sufficient.\n", n_imp))
  } else {
    cat(sprintf("  WARNING: %d estimate(s) have MC error >= 10%% of SE:\n", nrow(mc_flags)))
    for (i in 1:nrow(mc_flags)) {
      cat(sprintf("    %-25s %-20s  MC error = %.1f%% of SE, FMI = %.3f\n",
                  mc_flags$Outcome[i], mc_flags$term[i],
                  mc_flags$mc_pct_of_se[i], mc_flags$fmi[i]))
    }
    cat("  Consider increasing m for these estimates.\n")
  }

  list(all_results = results, eth_results = eth_res)
}

# --- Run models ---
unadj <- run_models(imp_list, outcomes,
                    formula_rhs = "eth_8cat",
                    model_label = "Unadjusted")

adj <- run_models(imp_list, outcomes,
                  formula_rhs = "eth_8cat + gender + grade + fas",
                  model_label = "Adjusted")

# --- Save CSVs ---
write.csv(unadj$all_results, file.path(OUTDIR, "Unadjusted_Results.csv"), row.names = FALSE)
write.csv(adj$all_results,   file.path(OUTDIR, "Adjusted_Results.csv"),   row.names = FALSE)

# ══════════════════════════════════════════════════════════════
# PART D: COMBINED WORD DOCUMENT — ALL TABLES
# ══════════════════════════════════════════════════════════════
cat("\n====== CREATING COMBINED WORD DOCUMENT ======\n")

# --- Helper: build wide results table from model output ---
build_results_wide <- function(results_df, eth_order, outcome_groups,
                               outcome_labels, eth_n, ref_group = "Gypsy/Traveller") {
  format_rr <- function(rr, lo, hi) {
    sprintf("%.2f (%.2f, %.2f)", rr, lo, hi)
  }

  eth_res <- results_df %>%
    filter(grepl("eth_8cat", term)) %>%
    mutate(ethnicity = gsub("eth_8cat", "", term))

  col_names <- c("Outcome", eth_order)
  rows_list <- list()
  row_types <- character()

  # N row
  n_row_vals <- c(paste0("N                                                        ",
                          eth_n$n_pct[match(ref_group, eth_n$Ethnicity)]))
  for (eth in eth_order) {
    if (eth == ref_group) {
      n_row_vals <- c(n_row_vals, "")
    } else {
      n_row_vals <- c(n_row_vals, eth_n$n_pct[match(eth, eth_n$Ethnicity)])
    }
  }
  rows_list[[1]] <- setNames(as.data.frame(t(n_row_vals), stringsAsFactors = FALSE), col_names)
  row_types <- c(row_types, "n")

  for (section_name in names(outcome_groups)) {
    sec_vals <- c(section_name, rep("", length(eth_order)))
    rows_list[[length(rows_list) + 1]] <- setNames(
      as.data.frame(t(sec_vals), stringsAsFactors = FALSE), col_names)
    row_types <- c(row_types, "section")

    for (oc in outcome_groups[[section_name]]) {
      label <- outcome_labels[oc]
      sub <- eth_res %>% filter(Outcome == oc)

      vals <- c(label)
      for (eth in eth_order) {
        if (eth == ref_group) {
          vals <- c(vals, "1.00 (ref)")
        } else {
          eth_row <- sub %>% filter(ethnicity == eth)
          if (nrow(eth_row) == 1) {
            vals <- c(vals, format_rr(eth_row$RR, eth_row$L, eth_row$H))
          } else {
            vals <- c(vals, "")
          }
        }
      }
      rows_list[[length(rows_list) + 1]] <- setNames(
        as.data.frame(t(vals), stringsAsFactors = FALSE), col_names)
      row_types <- c(row_types, "data")
    }
  }

  list(df = bind_rows(rows_list), row_types = row_types)
}

# --- Build wide tables ---
unadj_wide <- build_results_wide(unadj$all_results, eth_order, outcome_groups,
                                  outcome_labels, eth_n, ref_group = "Gypsy/Traveller")
adj_wide   <- build_results_wide(adj$all_results, eth_order, outcome_groups,
                                  outcome_labels, eth_n, ref_group = "Gypsy/Traveller")

# --- Helper: apply csv-to-table formatting to a flextable ---
format_results_ft <- function(tbl_df, row_types) {
  n_row_idx      <- which(row_types == "n")
  section_idx    <- which(row_types == "section")
  data_idx       <- which(row_types == "data")
  ncols          <- ncol(tbl_df)

  ft <- flextable(tbl_df)

  ft <- width(ft, j = 1, width = 2.3)
  for (j in 2:ncols) ft <- width(ft, j = j, width = 1.15)

  ft <- font(ft, fontname = "Aptos", part = "all")
  ft <- fontsize(ft, size = 11, part = "header")
  ft <- fontsize(ft, size = 10, part = "body")

  ft <- bold(ft, part = "header")

  ft <- fontsize(ft, i = data_idx, j = 1, size = 11)

  if (length(n_row_idx) > 0) {
    ft <- fontsize(ft, i = n_row_idx, size = 11)
    ft <- bold(ft, i = n_row_idx, j = 1)
  }

  if (length(section_idx) > 0) {
    ft <- fontsize(ft, i = section_idx, size = 11)
    ft <- bold(ft, i = section_idx, j = 1)
    ft <- italic(ft, i = section_idx, j = 1)
  }

  ft <- align(ft, j = 1, align = "left", part = "all")
  if (ncols > 1) ft <- align(ft, j = 2:ncols, align = "center", part = "all")

  ft <- padding(ft, padding.top = 2, padding.bottom = 2,
                padding.left = 3, padding.right = 3, part = "all")

  ft <- border_remove(ft)
  ft <- hline_top(ft, part = "header",
                  border = fp_border(color = "black", width = 1.5))
  ft <- hline_bottom(ft, part = "header",
                     border = fp_border(color = "black", width = 1))
  if (length(n_row_idx) > 0) {
    ft <- hline(ft, i = n_row_idx, border = fp_border(color = "black", width = 0.5),
                part = "body")
  }
  ft <- hline_bottom(ft, part = "body",
                     border = fp_border(color = "black", width = 1.5))

  if (length(n_row_idx) > 0) {
    ft <- merge_at(ft, i = n_row_idx, j = 1:2, part = "body")
  }
  for (si in section_idx) {
    ft <- merge_at(ft, i = si, j = 1:ncols, part = "body")
  }

  ft
}

# --- Helper: apply csv-to-table formatting to Table 1 ---
format_table1_ft <- function(tab1_df) {
  ft <- flextable(tab1_df)
  ft <- set_header_labels(ft, Variable = "Variable", Level = "Level")

  for (i in 10:13) ft <- merge_at(ft, i = i, j = 1:2, part = "body")
  ft <- merge_at(ft, i = 2:4, j = 1, part = "body")
  ft <- merge_at(ft, i = 5:9, j = 1, part = "body")

  ft <- font(ft, fontname = "Aptos", part = "all")
  ft <- fontsize(ft, size = 11, part = "header")
  ft <- fontsize(ft, size = 10, part = "body")
  ft <- bold(ft, part = "header")
  ft <- bold(ft, j = 1, part = "body")
  ft <- bold(ft, i = 1, part = "body")
  ft <- align(ft, j = 1:2, align = "left", part = "all")
  ft <- align(ft, j = 3:10, align = "center", part = "all")
  ft <- valign(ft, valign = "top", part = "body")
  ft <- width(ft, j = 1, width = 1.5)
  ft <- width(ft, j = 2, width = 1.3)
  ft <- width(ft, j = 3:10, width = 1.05)
  ft <- padding(ft, padding.top = 2, padding.bottom = 2,
                padding.left = 3, padding.right = 3, part = "all")
  ft <- border_remove(ft)
  ft <- hline_top(ft, part = "header",
                  border = fp_border(color = "black", width = 1.5))
  ft <- hline_bottom(ft, part = "header",
                     border = fp_border(color = "black", width = 1))
  ft <- hline_bottom(ft, part = "body",
                     border = fp_border(color = "black", width = 1.5))
  ft
}

# --- Helper: apply csv-to-table formatting to descriptive table ---
format_descriptive_ft <- function(wide_df) {
  ncols <- ncol(wide_df)
  ft <- flextable(wide_df)

  ft <- font(ft, fontname = "Aptos", part = "all")
  ft <- fontsize(ft, size = 11, part = "header")
  ft <- fontsize(ft, size = 10, part = "body")
  ft <- bold(ft, part = "header")

  ft <- bold(ft, i = 1, part = "body")
  ft <- fontsize(ft, i = 1, size = 11)

  for (i in 2:nrow(wide_df)) {
    if (all(wide_df[i, 2:ncols] == "")) {
      ft <- bold(ft, i = i, j = 1)
      ft <- italic(ft, i = i, j = 1)
      ft <- fontsize(ft, i = i, j = 1, size = 11)
      ft <- merge_at(ft, i = i, j = 1:ncols, part = "body")
    }
  }

  ft <- align(ft, j = 1, align = "left", part = "all")
  if (ncols > 1) ft <- align(ft, j = 2:ncols, align = "center", part = "all")
  ft <- width(ft, j = 1, width = 2.3)
  for (j in 2:ncols) ft <- width(ft, j = j, width = 1.15)
  ft <- padding(ft, padding.top = 2, padding.bottom = 2,
                padding.left = 3, padding.right = 3, part = "all")
  ft <- border_remove(ft)
  ft <- hline_top(ft, part = "header",
                  border = fp_border(color = "black", width = 1.5))
  ft <- hline_bottom(ft, part = "header",
                     border = fp_border(color = "black", width = 1))
  ft <- hline(ft, i = 1, border = fp_border(color = "black", width = 0.5),
              part = "body")
  ft <- hline_bottom(ft, part = "body",
                     border = fp_border(color = "black", width = 1.5))
  ft
}

# --- Build flextables ---
cat("Formatting Table 1...\n")
ft_table1 <- format_table1_ft(tab1_df)

cat("Formatting Descriptive Table...\n")
ft_desc <- format_descriptive_ft(wide_df)

cat("Formatting Unadjusted Results Table...\n")
ft_unadj <- format_results_ft(unadj_wide$df, unadj_wide$row_types)

cat("Formatting Adjusted Results Table...\n")
ft_adj <- format_results_ft(adj_wide$df, adj_wide$row_types)

# --- Assemble Word document ---
cat("Assembling combined Word document...\n")

doc <- read_docx()

# Table 1
doc <- body_end_section_landscape(doc)
doc <- body_add_par(doc, "Table 1. Sample characteristics by ethnicity",
                    style = "Normal")
doc <- body_add_flextable(doc, ft_table1)
doc <- body_end_section_landscape(doc)

# Descriptive Table
doc <- body_add_par(doc,
  "Table 2. Pooled prevalence estimates (%, 95% CI) for all outcomes by ethnic group",
  style = "Normal")
doc <- body_add_flextable(doc, ft_desc)
doc <- body_end_section_landscape(doc)

# Unadjusted Results
doc <- body_add_par(doc,
  "Table S1. Unadjusted risk ratios (95% confidence interval) for all outcomes by ethnic group (ref: Gypsy/Traveller)",
  style = "Normal")
doc <- body_add_flextable(doc, ft_unadj)
doc <- body_end_section_landscape(doc)

# Adjusted Results
doc <- body_add_par(doc,
  "Table S2. Adjusted risk ratios (95% confidence interval) for all outcomes by ethnic group (ref: Gypsy/Traveller)",
  style = "Normal")
doc <- body_add_flextable(doc, ft_adj)
doc <- body_end_section_landscape(doc)

doc_path <- file.path(OUTDIR, "All_Tables.docx")
print(doc, target = doc_path)
cat(sprintf("Saved: %s\n", doc_path))

# ══════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════
cat(sprintf("\n\nTotal script time: %.1f minutes.\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("\nAll outputs saved to:", OUTDIR, "\n")
cat("  1.  Descriptive_Table_Wide.csv\n")
cat("  2.  Table1_Sample_Characteristics.csv\n")
cat("  3.  Unadjusted_Results.csv\n")
cat("  4.  Adjusted_Results.csv\n")
cat("  5.  All_Tables.docx (Table 1 + Descriptive + Unadjusted + Adjusted)\n")
cat("\nScript completed.\n")


# ************************************************************
# ************************************************************
#
# SECTION 2: TABLE S2 (adjusted for gender) and
#            TABLE S3 (adjusted for grade)
# Source: Sophie_adjusted_separate_v0.8.R
#
# Table S2 -- Risk ratios adjusted for gender identity only
# Table S3 -- Risk ratios adjusted for school year only
#
# ************************************************************
# ************************************************************

# ============================================================
# SCRIPT: Poisson Models — eth_8cat adjusted for gender only
#         AND eth_8cat adjusted for grade only (separate models)
# ============================================================
# - Exposure: eth_8cat (categorical, ref = "White")
# - Outcomes: 16 binary outcomes grouped by domain
# - Model 1: eth_8cat + gender
# - Model 2: eth_8cat + grade
# - Cluster-robust SEs by id2c (school cluster)
# - Output: Separate CSVs + combined Word table
# ============================================================

# ── 1. Setup ──
user_lib <- file.path(Sys.getenv("HOME"), "R_libs")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))

pkgs <- c("parallel", "dplyr", "clubSandwich", "mice", "tidyr",
           "flextable", "officer")
new_pkgs <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(new_pkgs) > 0) {
  install.packages(new_pkgs, lib = user_lib, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(parallel)
  library(dplyr)
  library(clubSandwich)
  library(mice)
  library(tidyr)
  library(flextable)
  library(officer)
})

# --- Configuration ---
N_CORES <- 3L
BASEDIR <- "/shared/home1/c.wppjw/SHRN/Gypsy RT"
OUTDIR  <- "/shared/home1/c.wppjw/SHRN/Gypsy RT/output/v0.8"
if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

IMP_PATH <- file.path(BASEDIR, "output/v0.4/GRT_2023_2025_imp_m44_parallel.rds")

# --- Outcomes ---
outcomes <- c("lifesatneg", "sdqinternalhigh", "sdqexternalhigh",
              "cannabis30d", "alcoholuse", "smokeweekly", "regvape",
              "phys7", "vigexercise", "sit",
              "dailysoft", "dailyenergy",
              "beenbulliedever", "cbeenbulliedever",
              "truantever", "excludedever")

# --- Outcome labels ---
outcome_labels <- c(
  "lifesatneg"       = "Low life satisfaction",
  "sdqexternalhigh"  = "Externalizing difficulties\u1d43",
  "sdqinternalhigh"  = "Internalizing difficulties\u1d47",
  "smokeweekly"      = "Smoke weekly",
  "regvape"          = "Vape weekly",
  "alcoholuse"       = "\u22652 alcoholic drinks a session",
  "cannabis30d"      = "Cannabis in past 30 days",
  "phys7"            = "60 minutes of MVPA per day",
  "vigexercise"      = "\u22653 days VPA per week\u1d9c",
  "sit"              = "Sedentary behaviour",
  "dailysoft"        = "Daily soft drink consumption",
  "dailyenergy"      = "Daily energy drink consumption",
  "beenbulliedever"  = "Ever been bullied",
  "cbeenbulliedever" = "Ever been cyberbullied",
  "truantever"       = "Ever truanted",
  "excludedever"     = "Ever excluded from school"
)

# --- Outcome groups ---
outcome_groups <- list(
  "Mental health and life satisfaction" = c("lifesatneg", "sdqexternalhigh", "sdqinternalhigh"),
  "Substance use"     = c("smokeweekly", "regvape", "alcoholuse", "cannabis30d"),
  "Physical activity"  = c("phys7", "vigexercise", "sit"),
  "Diet"              = c("dailysoft", "dailyenergy"),
  "Bullying"          = c("beenbulliedever", "cbeenbulliedever"),
  "School problems"   = c("truantever", "excludedever")
)

# --- Ethnicity labels ---
eth_labels <- c(
  "1" = "White",
  "2" = "White Other",
  "3" = "Gypsy/Traveller",
  "4" = "Roma",
  "5" = "Black",
  "6" = "Mixed",
  "7" = "Asian",
  "8" = "Other"
)

# Column order for output tables (excluding White reference)
eth_col_order <- c("Gypsy/Traveller", "Roma", "White Other",
                   "Black", "Mixed", "Asian", "Other")

# ══════════════════════════════════════════════════════════════
# 2. LOAD DATA
# ══════════════════════════════════════════════════════════════
cat("Loading imputed data...\n")
t0 <- Sys.time()
imp <- readRDS(IMP_PATH)

n_imp <- imp$m
n_obs <- nrow(imp$data)

cat(sprintf("Loaded mids object: %d imputations, n = %d (%.1f sec)\n",
            n_imp, n_obs,
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
cat(sprintf("  mids object size: %.1f GB\n", object.size(imp) / 1e9))
cat(sprintf("  Current R memory usage: %.1f GB\n\n", sum(gc()[, 2]) / 1024))

# ── Prepare datasets ──
keep_cols <- unique(c("eth_8cat", "gender", "grade", "id2c", outcomes))

prepare_df <- function(df) {
  df <- df[, keep_cols, drop = FALSE]
  df$eth_8cat <- factor(df$eth_8cat, levels = names(eth_labels),
                        labels = eth_labels)
  df$eth_8cat <- relevel(df$eth_8cat, ref = "White")
  df$gender   <- as.factor(df$gender)
  df$grade    <- as.factor(df$grade)
  df$id2c     <- as.integer(as.factor(df$id2c))
  for (oc in outcomes) {
    if (is.factor(df[[oc]])) df[[oc]] <- as.integer(df[[oc]]) - 1L
  }
  df
}

# Extract all imputed datasets upfront for parallel processing
cat("Extracting all imputed datasets...\n")
imp_list <- lapply(seq_len(n_imp), function(i) {
  prepare_df(mice::complete(imp, action = i))
})
cat(sprintf("  Extracted %d datasets (%.1f GB total)\n",
            n_imp, object.size(imp_list) / 1e9))

# Free the mids object
rm(imp); gc(verbose = FALSE)
cat(sprintf("  Memory after freeing mids object: %.1f GB\n\n", sum(gc()[, 2]) / 1024))

# ── Diagnostics ──
cat("=== SYSTEM DIAGNOSTICS ===\n")
cat(sprintf("  R version: %s\n", R.version.string))
cat(sprintf("  Platform: %s\n", R.version$platform))
cat(sprintf("  N cores requested: %d\n", N_CORES))
cat(sprintf("  N cores detected: %d\n", parallel::detectCores()))
mem_info <- tryCatch(system("free -h", intern = TRUE), error = function(e) NULL)
if (!is.null(mem_info)) {
  cat("  System memory:\n")
  for (line in mem_info) cat("    ", line, "\n")
}
slurm_mem <- Sys.getenv("SLURM_MEM_PER_NODE", unset = NA)
if (!is.na(slurm_mem)) cat(sprintf("  SLURM allocated memory: %s MB\n", slurm_mem))
cat(sprintf("  R process RSS: %.1f GB\n", sum(gc()[, 2]) / 1024))
cat(sprintf("  imp_list object size: %.2f GB\n", object.size(imp_list) / 1e9))
cat(sprintf("  Single dataset size: %.1f MB\n", object.size(imp_list[[1]]) / 1e6))
cat(sprintf("  Estimated peak memory (parent + %d forks): %.1f GB\n",
            N_CORES,
            (sum(gc()[, 2]) / 1024) * (1 + N_CORES)))
cat("===========================\n\n")

# ══════════════════════════════════════════════════════════════
# 3. MODEL FUNCTIONS
# ══════════════════════════════════════════════════════════════

fit_model <- function(df, outcome_var, formula_rhs) {
  tryCatch({
    fml <- as.formula(paste0(outcome_var, " ~ ", formula_rhs))
    mod <- glm(fml, data = df, family = poisson(link = "log"))
    cf  <- coef(mod)
    vc  <- clubSandwich::vcovCR(mod, cluster = df$id2c, type = "CR1S")
    rm(mod)
    return(list(coef = cf, vcov = vc, error = NULL))
  }, error = function(e) return(list(error = e$message)))
}

pool_rubin <- function(fit_list) {
  results_only <- Filter(function(x) !is.null(x) && is.null(x$error), fit_list)
  m <- length(results_only)
  if (m == 0) stop("All models failed (check for OOM kills or model errors).")
  if (m < 3) warning(sprintf("Only %d/%d imputations succeeded -- estimates unreliable.", m, length(fit_list)))

  q_m <- do.call(rbind, lapply(results_only, `[[`, "coef"))
  u_m <- lapply(results_only, `[[`, "vcov")

  q_bar <- colMeans(q_m)
  u_bar <- Reduce("+", u_m) / m
  b_m   <- cov(q_m)
  total_vcov <- u_bar + (1 + 1/m) * b_m

  se <- sqrt(diag(total_vcov))

  r <- (1 + 1/m) * diag(b_m) / diag(u_bar)
  v_old <- (m - 1) * (1 + 1/r)^2
  v_old[is.nan(v_old) | is.na(v_old)] <- Inf

  data.frame(
    term = names(q_bar),
    est  = q_bar,
    se   = se,
    RR   = exp(q_bar),
    L    = exp(q_bar - 1.96 * se),
    H    = exp(q_bar + 1.96 * se),
    p    = 2 * pt(-abs(q_bar / se), df = v_old),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

run_models <- function(formula_rhs, model_label) {
  cat(sprintf("\n====== %s ======\n", model_label))
  cat(sprintf("Formula: outcome ~ %s\n", formula_rhs))
  cat(sprintf("Using %d cores for parallel processing\n\n", N_CORES))

  all_res <- list()
  t_models <- Sys.time()
  for (idx in seq_along(outcomes)) {
    out <- outcomes[idx]
    t_out <- Sys.time()
    cat(sprintf("  [%d/%d] Model: %s ...", idx, length(outcomes), out))
    fits <- mclapply(imp_list, fit_model, outcome_var = out,
                     formula_rhs = formula_rhs,
                     mc.cores = N_CORES, mc.preschedule = FALSE)
    n_null <- sum(sapply(fits, is.null))
    n_err  <- sum(sapply(fits, function(x) !is.null(x$error)))
    n_ok   <- length(fits) - n_null - n_err
    elapsed_out <- as.numeric(difftime(Sys.time(), t_out, units = "secs"))
    elapsed_total <- as.numeric(difftime(Sys.time(), t_models, units = "secs"))
    avg_per <- elapsed_total / idx
    remaining <- avg_per * (length(outcomes) - idx)
    cat(sprintf(" %d/%d converged (%.1fs, ~%.0fs remaining)\n",
                n_ok, length(imp_list), elapsed_out, remaining))
    if (n_null > 0) {
      cat(sprintf("  *** WARNING: %d workers returned NULL (likely OOM killed) ***\n", n_null))
      cat(sprintf("  *** R memory: %.1f GB | Consider increasing --mem ***\n", sum(gc()[, 2]) / 1024))
    }
    if (n_err > 0) {
      errs <- sapply(Filter(function(x) !is.null(x$error), fits), `[[`, "error")
      cat(sprintf("  *** %d model errors: %s ***\n", n_err, paste(unique(errs), collapse = "; ")))
    }
    all_res[[out]] <- pool_rubin(fits) %>% mutate(Outcome = out)
    rm(fits); gc(verbose = FALSE)
  }
  cat(sprintf("  Models completed in %.1f seconds.\n\n",
              as.numeric(difftime(Sys.time(), t_models, units = "secs"))))

  results <- bind_rows(all_res)

  # Extract ethnicity results
  eth_res <- results %>%
    filter(grepl("eth_8cat", term)) %>%
    mutate(term = gsub("eth_8cat", "", term))

  # Print results
  cat(sprintf("\n-- Ethnicity RRs (ref = White), %s --\n", model_label))
  for (out in outcomes) {
    cat(sprintf("\n--- %s ---\n", out))
    d <- eth_res %>% filter(Outcome == out)
    for (i in 1:nrow(d)) {
      cat(sprintf("  %-20s  RR = %5.2f (%4.2f\u2013%5.2f)  p = %.4f\n",
                  d$term[i], d$RR[i], d$L[i], d$H[i], d$p[i]))
    }
  }

  eth_res
}

# ══════════════════════════════════════════════════════════════
# 4. RUN MODELS
# ══════════════════════════════════════════════════════════════

# Model 1: adjusted for gender only
eth_res_gender <- run_models("eth_8cat + gender", "ADJUSTED FOR GENDER ONLY")

# Model 2: adjusted for grade only
eth_res_grade <- run_models("eth_8cat + grade", "ADJUSTED FOR GRADE ONLY")

# ══════════════════════════════════════════════════════════════
# 5. SAVE CSVs
# ══════════════════════════════════════════════════════════════

write.csv(eth_res_gender, file.path(OUTDIR, "Adjusted_Gender_Only_Results.csv"), row.names = FALSE)
cat("\nSaved: Adjusted_Gender_Only_Results.csv\n")

write.csv(eth_res_grade, file.path(OUTDIR, "Adjusted_Grade_Only_Results.csv"), row.names = FALSE)
cat("Saved: Adjusted_Grade_Only_Results.csv\n")

# ══════════════════════════════════════════════════════════════
# 6. COMBINED WORD TABLE
# ══════════════════════════════════════════════════════════════

build_wide_table <- function(eth_res) {
  eth_res$formatted <- sprintf("%.2f (%.2f\u2013%.2f)", eth_res$RR, eth_res$L, eth_res$H)

  rows <- list()
  row_types <- character()

  for (grp_name in names(outcome_groups)) {
    # Group header row
    rows[[length(rows) + 1]] <- c(grp_name, rep("", length(eth_col_order)))
    row_types <- c(row_types, "header")

    for (oc in outcome_groups[[grp_name]]) {
      vals <- c(outcome_labels[oc])
      for (eth in eth_col_order) {
        match_row <- eth_res[eth_res$term == eth & eth_res$Outcome == oc, ]
        if (nrow(match_row) > 0) {
          vals <- c(vals, match_row$formatted[1])
        } else {
          vals <- c(vals, "")
        }
      }
      rows[[length(rows) + 1]] <- vals
      row_types <- c(row_types, "outcome")
    }
  }

  wide_df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(wide_df) <- c("Outcome", eth_col_order)

  list(df = wide_df, row_types = row_types)
}

make_flextable <- function(wide_df, row_types, caption) {
  ft <- flextable(wide_df) %>%
    set_caption(caption = caption) %>%
    fontsize(size = 8, part = "all") %>%
    font(fontname = "Arial", part = "all") %>%
    autofit() %>%
    set_table_properties(layout = "autofit")

  header_idx <- which(row_types == "header")
  if (length(header_idx) > 0) {
    ft <- bold(ft, i = header_idx, part = "body")
  }
  ft
}

# Build tables
tbl_gender <- build_wide_table(eth_res_gender)
tbl_grade  <- build_wide_table(eth_res_grade)

ft_gender <- make_flextable(
  tbl_gender$df, tbl_gender$row_types,
  "Risk ratios (95% CI) adjusted for gender only, reference: White British"
)
ft_grade <- make_flextable(
  tbl_grade$df, tbl_grade$row_types,
  "Risk ratios (95% CI) adjusted for school year only, reference: White British"
)

# Save combined Word document
doc <- read_docx() %>%
  body_add_par("Table S_. Risk ratios adjusted for gender only", style = "heading 1") %>%
  body_add_par("Poisson regression with cluster-robust SEs. Reference group: White British.", style = "Normal") %>%
  body_add_flextable(ft_gender) %>%
  body_add_par("", style = "Normal") %>%
  body_add_par("Table S_. Risk ratios adjusted for school year only", style = "heading 1") %>%
  body_add_par("Poisson regression with cluster-robust SEs. Reference group: White British.", style = "Normal") %>%
  body_add_flextable(ft_grade)

# Set landscape orientation
doc <- body_end_block_section(doc, block_section(
  prop_section(page_size = page_size(orient = "landscape"))
))

print(doc, target = file.path(OUTDIR, "Adjusted_Separate_Tables.docx"))
cat("Saved: Adjusted_Separate_Tables.docx\n")

# ══════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════
cat(sprintf("\n\nTotal script time: %.1f minutes.\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("\nAll outputs saved to:", OUTDIR, "\n")
cat("  1.  Adjusted_Gender_Only_Results.csv\n")
cat("  2.  Adjusted_Grade_Only_Results.csv\n")
cat("  3.  Adjusted_Separate_Tables.docx\n")
cat("\nScript completed.\n")


# ************************************************************
# ************************************************************
#
# SECTION 3: TABLE S4
# Source: Sophie_adjusted_gender_grade_v0.7.R
#
# Table S4 -- Risk ratios adjusted for gender identity and grade
# Model: eth_8cat + gender + grade
#
# ************************************************************
# ************************************************************

# ============================================================
# SCRIPT: Poisson Models — eth_8cat adjusted for gender + grade
# ============================================================
# - Exposure: eth_8cat (categorical, ref = "White")
# - Outcomes: 16 binary outcomes grouped by domain
# - Model: eth_8cat + gender + grade
# - Cluster-robust SEs by id2c (school cluster)
# - Output: Results CSV + Word table
# ============================================================

# ── 1. Setup ──
user_lib <- file.path(Sys.getenv("HOME"), "R_libs")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))

pkgs <- c("parallel", "dplyr", "clubSandwich", "mice", "tidyr",
           "flextable", "officer")
new_pkgs <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(new_pkgs) > 0) {
  install.packages(new_pkgs, lib = user_lib, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(parallel)
  library(dplyr)
  library(clubSandwich)
  library(mice)
  library(tidyr)
  library(flextable)
  library(officer)
})

# --- Configuration ---
N_CORES <- 3L
BASEDIR <- "/shared/home1/c.wppjw/SHRN/Gypsy RT"
OUTDIR  <- "/shared/home1/c.wppjw/SHRN/Gypsy RT/output/v0.7"
if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

IMP_PATH <- file.path(BASEDIR, "output/v0.4/GRT_2023_2025_imp_m44_parallel.rds")

# --- Outcomes ---
outcomes <- c("lifesatneg", "sdqinternalhigh", "sdqexternalhigh",
              "cannabis30d", "alcoholuse", "smokeweekly", "regvape",
              "phys7", "vigexercise", "sit",
              "dailysoft", "dailyenergy",
              "beenbulliedever", "cbeenbulliedever",
              "truantever", "excludedever")

# --- Outcome labels ---
outcome_labels <- c(
  "lifesatneg"       = "Low life satisfaction",
  "sdqexternalhigh"  = "Externalizing difficulties\u1d43",
  "sdqinternalhigh"  = "Internalizing difficulties\u1d47",
  "smokeweekly"      = "Smoke weekly",
  "regvape"          = "Vape weekly",
  "alcoholuse"       = "\u22652 alcoholic drinks a session",
  "cannabis30d"      = "Cannabis in past 30 days",
  "phys7"            = "60 minutes of MVPA per day",
  "vigexercise"      = "\u22653 days VPA per week\u1d9c",
  "sit"              = "Sedentary behaviour",
  "dailysoft"        = "Daily soft drink consumption",
  "dailyenergy"      = "Daily energy drink consumption",
  "beenbulliedever"  = "Ever been bullied",
  "cbeenbulliedever" = "Ever been cyberbullied",
  "truantever"       = "Ever truanted",
  "excludedever"     = "Ever excluded from school"
)

# --- Outcome groups ---
outcome_groups <- list(
  "Mental health and life satisfaction" = c("lifesatneg", "sdqexternalhigh", "sdqinternalhigh"),
  "Substance use"     = c("smokeweekly", "regvape", "alcoholuse", "cannabis30d"),
  "Physical activity"  = c("phys7", "vigexercise", "sit"),
  "Diet"              = c("dailysoft", "dailyenergy"),
  "Bullying"          = c("beenbulliedever", "cbeenbulliedever"),
  "School problems"   = c("truantever", "excludedever")
)

# --- Ethnicity labels ---
eth_labels <- c(
  "1" = "White",
  "2" = "White Other",
  "3" = "Gypsy/Traveller",
  "4" = "Roma",
  "5" = "Black",
  "6" = "Mixed",
  "7" = "Asian",
  "8" = "Other"
)

eth_order <- c("White", "White Other", "Gypsy/Traveller", "Roma",
               "Black", "Mixed", "Asian", "Other")

# ══════════════════════════════════════════════════════════════
# 2. LOAD DATA
# ══════════════════════════════════════════════════════════════
cat("Loading imputed data...\n")
t0 <- Sys.time()
imp <- readRDS(IMP_PATH)

n_imp <- imp$m
n_obs <- nrow(imp$data)

cat(sprintf("Loaded mids object: %d imputations, n = %d (%.1f sec)\n",
            n_imp, n_obs,
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
cat(sprintf("  mids object size: %.1f GB\n", object.size(imp) / 1e9))
cat(sprintf("  Current R memory usage: %.1f GB\n\n", sum(gc()[, 2]) / 1024))

# ── Prepare datasets ──
keep_cols <- unique(c("eth_8cat", "gender", "grade", "id2c", outcomes))

prepare_df <- function(df) {
  df <- df[, keep_cols, drop = FALSE]
  df$eth_8cat <- factor(df$eth_8cat, levels = names(eth_labels),
                        labels = eth_labels)
  df$eth_8cat <- relevel(df$eth_8cat, ref = "White")
  df$gender   <- as.factor(df$gender)
  df$grade    <- as.factor(df$grade)
  df$id2c     <- as.integer(as.factor(df$id2c))
  for (oc in outcomes) {
    if (is.factor(df[[oc]])) df[[oc]] <- as.integer(df[[oc]]) - 1L
  }
  df
}

# Extract all imputed datasets upfront for parallel processing
cat("Extracting all imputed datasets...\n")
imp_list <- lapply(seq_len(n_imp), function(i) {
  prepare_df(mice::complete(imp, action = i))
})
cat(sprintf("  Extracted %d datasets (%.1f GB total)\n",
            n_imp, object.size(imp_list) / 1e9))

# Free the mids object
rm(imp); gc(verbose = FALSE)
cat(sprintf("  Memory after freeing mids object: %.1f GB\n\n", sum(gc()[, 2]) / 1024))

# ── Diagnostics ──
cat("=== SYSTEM DIAGNOSTICS ===\n")
cat(sprintf("  R version: %s\n", R.version.string))
cat(sprintf("  Platform: %s\n", R.version$platform))
cat(sprintf("  N cores requested: %d\n", N_CORES))
cat(sprintf("  N cores detected: %d\n", parallel::detectCores()))
mem_info <- tryCatch(system("free -h", intern = TRUE), error = function(e) NULL)
if (!is.null(mem_info)) {
  cat("  System memory:\n")
  for (line in mem_info) cat("    ", line, "\n")
}
slurm_mem <- Sys.getenv("SLURM_MEM_PER_NODE", unset = NA)
if (!is.na(slurm_mem)) cat(sprintf("  SLURM allocated memory: %s MB\n", slurm_mem))
cat(sprintf("  R process RSS: %.1f GB\n", sum(gc()[, 2]) / 1024))
cat(sprintf("  imp_list object size: %.2f GB\n", object.size(imp_list) / 1e9))
cat(sprintf("  Single dataset size: %.1f MB\n", object.size(imp_list[[1]]) / 1e6))
cat(sprintf("  Estimated peak memory (parent + %d forks): %.1f GB\n",
            N_CORES,
            (sum(gc()[, 2]) / 1024) * (1 + N_CORES)))
cat("===========================\n\n")

# ══════════════════════════════════════════════════════════════
# POISSON REGRESSION — eth_8cat + gender + grade
# ══════════════════════════════════════════════════════════════

fit_model <- function(df, outcome_var, formula_rhs) {
  tryCatch({
    fml <- as.formula(paste0(outcome_var, " ~ ", formula_rhs))
    mod <- glm(fml, data = df, family = poisson(link = "log"))
    cf  <- coef(mod)
    vc  <- clubSandwich::vcovCR(mod, cluster = df$id2c, type = "CR1S")
    rm(mod)
    return(list(coef = cf, vcov = vc, error = NULL))
  }, error = function(e) return(list(error = e$message)))
}

pool_rubin <- function(fit_list) {
  results_only <- Filter(function(x) !is.null(x) && is.null(x$error), fit_list)
  m <- length(results_only)
  if (m == 0) stop("All models failed (check for OOM kills or model errors).")
  if (m < 3) warning(sprintf("Only %d/%d imputations succeeded — estimates unreliable.", m, length(fit_list)))

  q_m <- do.call(rbind, lapply(results_only, `[[`, "coef"))
  u_m <- lapply(results_only, `[[`, "vcov")

  q_bar <- colMeans(q_m)
  u_bar <- Reduce("+", u_m) / m
  b_m   <- cov(q_m)
  total_vcov <- u_bar + (1 + 1/m) * b_m

  se <- sqrt(diag(total_vcov))

  r <- (1 + 1/m) * diag(b_m) / diag(u_bar)
  v_old <- (m - 1) * (1 + 1/r)^2
  v_old[is.nan(v_old) | is.na(v_old)] <- Inf

  data.frame(
    term = names(q_bar),
    est  = q_bar,
    se   = se,
    RR   = exp(q_bar),
    L    = exp(q_bar - 1.96 * se),
    H    = exp(q_bar + 1.96 * se),
    p    = 2 * pt(-abs(q_bar / se), df = v_old),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

FORMULA_RHS <- "eth_8cat + gender + grade"

cat(sprintf("\n====== ADJUSTED MODELS (gender + grade) ======\n"))
cat(sprintf("Formula: outcome ~ %s\n", FORMULA_RHS))
cat(sprintf("Using %d cores for parallel processing\n\n", N_CORES))

all_res <- list()
t_models <- Sys.time()
for (idx in seq_along(outcomes)) {
  out <- outcomes[idx]
  t_out <- Sys.time()
  cat(sprintf("  [%d/%d] Model: %s ...", idx, length(outcomes), out))
  fits <- mclapply(imp_list, fit_model, outcome_var = out,
                   formula_rhs = FORMULA_RHS,
                   mc.cores = N_CORES, mc.preschedule = FALSE)
  n_null <- sum(sapply(fits, is.null))
  n_err  <- sum(sapply(fits, function(x) !is.null(x$error)))
  n_ok   <- length(fits) - n_null - n_err
  elapsed_out <- as.numeric(difftime(Sys.time(), t_out, units = "secs"))
  elapsed_total <- as.numeric(difftime(Sys.time(), t_models, units = "secs"))
  avg_per <- elapsed_total / idx
  remaining <- avg_per * (length(outcomes) - idx)
  cat(sprintf(" %d/%d converged (%.1fs, ~%.0fs remaining)\n",
              n_ok, length(imp_list), elapsed_out, remaining))
  if (n_null > 0) {
    cat(sprintf("  *** WARNING: %d workers returned NULL (likely OOM killed) ***\n", n_null))
    cat(sprintf("  *** R memory: %.1f GB | Consider increasing --mem ***\n", sum(gc()[, 2]) / 1024))
  }
  if (n_err > 0) {
    errs <- sapply(Filter(function(x) !is.null(x$error), fits), `[[`, "error")
    cat(sprintf("  *** %d model errors: %s ***\n", n_err, paste(unique(errs), collapse = "; ")))
  }
  all_res[[out]] <- pool_rubin(fits) %>% mutate(Outcome = out)
  rm(fits); gc(verbose = FALSE)
}
cat(sprintf("  Models completed in %.1f seconds.\n\n",
            as.numeric(difftime(Sys.time(), t_models, units = "secs"))))

results <- bind_rows(all_res)

# Extract ethnicity results
eth_res <- results %>%
  filter(grepl("eth_8cat", term)) %>%
  mutate(
    term = gsub("eth_8cat", "", term)
  )

# Print results
cat("\n-- Ethnicity RRs (ref = White), adjusted for gender + grade --\n")
for (out in outcomes) {
  cat(sprintf("\n--- %s ---\n", out))
  d <- eth_res %>% filter(Outcome == out)
  for (i in 1:nrow(d)) {
    cat(sprintf("  %-20s  RR = %5.2f (%4.2f\u2013%5.2f)  p = %.4f\n",
                d$term[i], d$RR[i], d$L[i], d$H[i], d$p[i]))
  }
}

# ── Save results CSV ──
write.csv(eth_res, file.path(OUTDIR, "Adjusted_Gender_Grade_Results.csv"), row.names = FALSE)
cat("\nSaved: Adjusted_Gender_Grade_Results.csv\n")

# ══════════════════════════════════════════════════════════════
# WORD TABLE
# ══════════════════════════════════════════════════════════════

# Build wide table: rows = outcomes (grouped), columns = ethnicity groups
eth_res$formatted <- sprintf("%.2f (%.2f\u2013%.2f)", eth_res$RR, eth_res$L, eth_res$H)

rows <- list()
row_types <- character()

for (grp_name in names(outcome_groups)) {
  # Group header row
  rows[[length(rows) + 1]] <- c(grp_name, rep("", length(eth_order) - 1))
  row_types <- c(row_types, "header")

  for (oc in outcome_groups[[grp_name]]) {
    vals <- c(outcome_labels[oc])
    for (eth in eth_order[-1]) {  # skip White (reference)
      match_row <- eth_res[eth_res$term == eth & eth_res$Outcome == oc, ]
      if (nrow(match_row) > 0) {
        vals <- c(vals, match_row$formatted[1])
      } else {
        vals <- c(vals, "")
      }
    }
    rows[[length(rows) + 1]] <- vals
    row_types <- c(row_types, "outcome")
  }
}

wide_df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
names(wide_df) <- c("Outcome", eth_order[-1])  # exclude White (reference)

# Create flextable
ft <- flextable(wide_df) %>%
  set_caption(caption = "Risk ratios (95% CI) adjusted for gender and grade, reference: White") %>%
  fontsize(size = 8, part = "all") %>%
  font(fontname = "Arial", part = "all") %>%
  autofit() %>%
  set_table_properties(layout = "autofit")

# Bold header rows
header_idx <- which(row_types == "header")
if (length(header_idx) > 0) {
  ft <- bold(ft, i = header_idx, part = "body")
}

# Save to Word
doc <- read_docx() %>%
  body_add_par("Adjusted Risk Ratios (gender + grade)", style = "heading 1") %>%
  body_add_par("Poisson regression with cluster-robust SEs. Reference group: White.", style = "Normal") %>%
  body_add_flextable(ft)

# Set landscape orientation for the wide table
doc <- body_end_block_section(doc, block_section(
  prop_section(page_size = page_size(orient = "landscape"))
))

print(doc, target = file.path(OUTDIR, "Adjusted_Gender_Grade_Table.docx"))
cat("Saved: Adjusted_Gender_Grade_Table.docx\n")

# ══════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════
cat(sprintf("\n\nTotal script time: %.1f minutes.\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("\nAll outputs saved to:", OUTDIR, "\n")
cat("  1.  Adjusted_Gender_Grade_Results.csv\n")
cat("  2.  Adjusted_Gender_Grade_Table.docx\n")
cat("\nScript completed.\n")


# ************************************************************
# ************************************************************
#
# SECTION 4: TABLE S5 (fully adjusted) and
#            CSV used by FIGURE 1
# Source: grt_subgroup_adj_only_v0.1.R
#
# Table S5 -- Adjusted risk ratios for all outcomes by ethnic
#             group, adjusted for gender, grade, and FAS
#
# Also produces: Adjusted_Results_exact_p.csv which is read
# by make_figures.R to create Figure 1.
#
# ************************************************************
# ************************************************************

#!/usr/bin/env Rscript
# ============================================================
# SCRIPT: GRT subgroup - fully adjusted model only
# ------------------------------------------------------------
# Refits the single fully-adjusted Poisson model:
#     outcome ~ eth_8cat + gender + grade + fas
# on the subgroup mids object produced by grt_subgroup_v0.91.R
# (syear == 2025 + (syear == 2023 & grade in 10,11)).
#
# Differences from the descriptives/estimation scripts:
#   - Only the fully adjusted model is fit (no unadjusted pass,
#     no gender/grade-separate, no descriptives)
#   - Reference category: White (not Gypsy/Traveller)
#   - CSV contains ALL coefficient estimates (every term, not
#     just ethnicity), with EXACT p-values (no rounding).
#
# Output:
#   output/v0.91/Adjusted_Results_exact_p.csv
# ============================================================

# ── Library path ──
user_lib <- file.path(Sys.getenv("HOME"), "R_libs")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))

suppressPackageStartupMessages({
  library(parallel)
  library(dplyr)
  library(clubSandwich)
  library(mice)
})

# ── Configuration ──
N_CORES  <- 3L
BASEDIR  <- "/shared/home1/c.wppjw/SHRN/Gypsy RT"
OUTDIR   <- file.path(BASEDIR, "output/v0.91")
IMP_PATH <- file.path(OUTDIR, "GRT_subgroup_imp.rds")
CSV_OUT  <- file.path(OUTDIR, "Adjusted_Results_exact_p.csv")

if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)
if (!file.exists(IMP_PATH)) {
  stop("Subgroup mids not found at ", IMP_PATH,
       "\n  Run grt_subgroup_v0.91.R first (step 1 writes this file).")
}

# ── Outcomes ──
outcomes <- c("lifesatneg", "sdqinternalhigh", "sdqexternalhigh",
              "cannabis30d", "alcoholuse", "smokeweekly", "regvape",
              "phys7", "vigexercise", "sit",
              "dailysoft", "dailyenergy",
              "beenbulliedever", "cbeenbulliedever",
              "truantever", "excludedever")

eth_labels <- c("1" = "White", "2" = "White Other", "3" = "Gypsy/Traveller",
                "4" = "Roma",  "5" = "Black",       "6" = "Mixed",
                "7" = "Asian", "8" = "Other")

# ============================================================
# 1. LOAD SUBGROUP MIDS
# ============================================================
cat("Loading subgroup mids from ", IMP_PATH, "\n", sep = "")
t0  <- Sys.time()
imp <- readRDS(IMP_PATH)
n_imp <- imp$m
cat(sprintf("  Loaded: %d imputations, n = %d (%.1f sec)\n",
            n_imp, nrow(imp$data),
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# ============================================================
# 2. EXTRACT AND PREPARE IMPUTED DATASETS
# ============================================================
keep_cols <- unique(c("eth_8cat", "gender", "grade", "fas", "id2c", outcomes))

prepare_df <- function(df) {
  df <- df[, keep_cols, drop = FALSE]
  df$eth_8cat <- factor(df$eth_8cat, levels = names(eth_labels),
                        labels = eth_labels)
  df$eth_8cat <- relevel(df$eth_8cat, ref = "White")
  df$gender   <- as.factor(df$gender)
  df$grade    <- as.factor(df$grade)
  df$fas      <- as.numeric(df$fas)
  df$id2c     <- as.integer(as.factor(df$id2c))
  for (oc in outcomes) {
    if (is.factor(df[[oc]])) df[[oc]] <- as.integer(df[[oc]]) - 1L
  }
  df
}

cat("\nExtracting imputed datasets...\n")
imp_list <- lapply(seq_len(n_imp), function(i) prepare_df(mice::complete(imp, i)))
rm(imp); gc(verbose = FALSE)
cat(sprintf("  Extracted %d datasets (%.2f GB total)\n",
            n_imp, object.size(imp_list) / 1e9))

# ============================================================
# 3. FIT + POOL (Rubin's rules, cluster-robust SEs)
# ============================================================
fit_model <- function(df, outcome_var, formula_rhs) {
  tryCatch({
    fml <- as.formula(paste0(outcome_var, " ~ ", formula_rhs))
    mod <- glm(fml, data = df, family = poisson(link = "log"))
    cf  <- coef(mod)
    vc  <- clubSandwich::vcovCR(mod, cluster = df$id2c, type = "CR1S")
    rm(mod)
    list(coef = cf, vcov = vc, error = NULL)
  }, error = function(e) list(error = e$message))
}

pool_rubin <- function(fit_list) {
  ok <- Filter(function(x) !is.null(x) && is.null(x$error), fit_list)
  m <- length(ok)
  if (m == 0) stop("All models failed.")
  if (m < 3) warning(sprintf("Only %d/%d imputations succeeded.", m, length(fit_list)))

  q_m   <- do.call(rbind, lapply(ok, `[[`, "coef"))
  u_m   <- lapply(ok, `[[`, "vcov")
  q_bar <- colMeans(q_m)
  u_bar <- Reduce("+", u_m) / m
  b_m   <- cov(q_m)
  total_vcov <- u_bar + (1 + 1/m) * b_m
  se    <- sqrt(diag(total_vcov))

  r     <- (1 + 1/m) * diag(b_m) / diag(u_bar)
  v_old <- (m - 1) * (1 + 1/r)^2
  v_old[is.nan(v_old) | is.na(v_old)] <- Inf

  data.frame(
    term = names(q_bar),
    estimate = q_bar,
    std.error = se,
    RR   = exp(q_bar),
    CI_lower = exp(q_bar - 1.96 * se),
    CI_upper = exp(q_bar + 1.96 * se),
    statistic = q_bar / se,
    df = v_old,
    p.value = 2 * pt(-abs(q_bar / se), df = v_old),  # full precision
    n_imp_converged = m,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

FORMULA_RHS <- "eth_8cat + gender + grade + fas"
cat("\n====== FULLY ADJUSTED MODELS ======\n")
cat("Formula: outcome ~ ", FORMULA_RHS, "\n", sep = "")
cat("Reference category: White\n")
cat(sprintf("Using %d cores for parallel processing\n\n", N_CORES))

all_res  <- list()
t_models <- Sys.time()
for (idx in seq_along(outcomes)) {
  out   <- outcomes[idx]
  t_out <- Sys.time()
  cat(sprintf("  [%d/%d] %s ...", idx, length(outcomes), out))
  fits  <- mclapply(imp_list, fit_model,
                    outcome_var = out, formula_rhs = FORMULA_RHS,
                    mc.cores = N_CORES, mc.preschedule = FALSE)
  n_ok  <- sum(sapply(fits, function(x) !is.null(x) && is.null(x$error)))
  elap  <- as.numeric(difftime(Sys.time(), t_out, units = "secs"))
  tot   <- as.numeric(difftime(Sys.time(), t_models, units = "secs"))
  rem   <- (tot / idx) * (length(outcomes) - idx)
  cat(sprintf(" %d/%d converged (%.1fs, ~%.0fs remaining)\n",
              n_ok, length(fits), elap, rem))
  pooled <- pool_rubin(fits)
  pooled$outcome <- out
  all_res[[out]] <- pooled
  rm(fits); gc(verbose = FALSE)
}
cat(sprintf("\n  Completed in %.1f seconds.\n",
            as.numeric(difftime(Sys.time(), t_models, units = "secs"))))

results <- bind_rows(all_res)
# Tidy column order: outcome first
results <- results[, c("outcome", "term", "estimate", "std.error",
                       "RR", "CI_lower", "CI_upper", "statistic",
                       "df", "p.value", "n_imp_converged")]

# ============================================================
# 4. WRITE CSV WITH EXACT P-VALUES
# ============================================================
# Use maximum precision so exact p-values survive write.csv rounding.
op <- options(scipen = 999, digits = 17)
write.csv(results, CSV_OUT, row.names = FALSE)
options(op)

cat(sprintf("\nSaved: %s  (%d rows)\n", CSV_OUT, nrow(results)))
cat("\n====== DONE ======\n")


# ************************************************************
# ************************************************************
#
# SECTION 5: FIGURE 1
# Source: Final R script for figures/make_figures.R
#
# Figure 1 -- Annotated heatmap of adjusted risk ratios for
#             health and behavioural outcomes for different
#             ethnic groups compared with White British
#
# Reads: Adjusted_Results_exact_p.csv (from Section 4)
#
# ************************************************************
# ************************************************************

## ============================================================
## Publication figure: Gypsy/Traveller health disparities
## Annotated heatmap of adjusted RRs
## ============================================================

library(tidyverse)
library(scales)

# ── Paths ──
results_dir <- "C:/Users/wppjw/OneDrive - Cardiff University/Papers/SHRN/Gyspy Irish Roma/Results/v0.91"
base_dir <- tryCatch(
  dirname(rstudioapi::getSourceEditorContext()$path),
  error = function(e) results_dir
)

adj <- read_csv(file.path(results_dir, "Adjusted_Results_exact_p.csv")) %>%
  rename(Outcome = outcome, p = p.value)

# ── Outcome labels & domain groupings ──
outcome_info <- tribble(
  ~Outcome,            ~label,                            ~domain,                        ~order,
  "lifesatneg",        "Low life satisfaction",            "Mental health & life satisfaction", 1,
  "sdqexternalhigh",   "Externalising difficulties",      "Mental health & life satisfaction", 2,
  "sdqinternalhigh",   "Internalising difficulties",      "Mental health & life satisfaction", 3,
  "smokeweekly",       "Weekly smoking",                  "Substance use",                    4,
  "regvape",           "Weekly vaping",                   "Substance use",                    5,
  "alcoholuse",        "\u22652 alcoholic drinks/session", "Substance use",                    6,
  "cannabis30d",       "Cannabis in past 30 days",        "Substance use",                    7,
  "phys7",             "<60 min MVPA per week",            "Physical activity",                8,
  "vigexercise",       "<3 days of VPA/week",              "Physical activity",                9,
  "sit",               "Sitting \u22657 hrs/day (weekdays)","Physical activity",               10,
  "dailysoft",         "Daily soft drink consumption",    "Diet",                            11,
  "dailyenergy",       "Daily energy drink consumption",  "Diet",                            12,
  "beenbulliedever",   "Ever been bullied",               "Bullying",                        13,
  "cbeenbulliedever",  "Ever been cyberbullied",          "Bullying",                        14,
  "truantever",        "Ever truanted",                   "School problems",                 15,
  "excludedever",      "Ever excluded from school",       "School problems",                 16
)

# Ethnicity mapping
eth_levels <- c("Gypsy/Traveller", "Roma", "Black", "Mixed",
                "White Other", "Asian", "Other")

eth_map <- c(
  "eth_8catGypsy/Traveller" = "Gypsy/Traveller",
  "eth_8catRoma"            = "Roma",
  "eth_8catBlack"           = "Black",
  "eth_8catMixed"           = "Mixed",
  "eth_8catWhite Other"     = "White Other",
  "eth_8catAsian"           = "Asian",
  "eth_8catOther"           = "Other"
)

# ── Prepare data ──
adj_eth <- adj %>%
  filter(str_starts(term, "eth_8cat")) %>%
  mutate(
    ethnicity = recode(term, !!!eth_map),
    ethnicity = factor(ethnicity, levels = eth_levels)
  ) %>%
  left_join(outcome_info, by = "Outcome") %>%
  filter(!is.na(label)) %>%
  mutate(
    label = factor(label, levels = rev(outcome_info$label))
  )

# ── Bonferroni correction ──
n_tests <- length(eth_levels) * nrow(outcome_info)
bonf_alpha <- 0.05 / n_tests
adj_eth <- adj_eth %>% mutate(sig_bonf = p < bonf_alpha)

# ── Colour palette ──
gt_colour <- "#D62728"

## ════════════════════════════════════════════════════════════
## Annotated heatmap of adjusted RRs
## ════════════════════════════════════════════════════════════

heat_data <- adj_eth %>%
  select(ethnicity, label, RR, p, sig_bonf, Outcome) %>%
  mutate(
    log_rr = log(RR),
    annot = case_when(
      sig_bonf & RR >= 10 ~ paste0(sprintf("%.1f", RR), "\u2020"),
      sig_bonf             ~ paste0(sprintf("%.2f", RR), "\u2020"),
      RR >= 10             ~ sprintf("%.1f", RR),
      TRUE                 ~ sprintf("%.2f", RR)
    ),
    fontface_val = if_else(ethnicity == "Gypsy/Traveller", "bold", "plain"),
    text_colour = if_else(abs(log_rr) > 1.2, "white", "black")
  )

# Cap log_rr for colour scale
cap_hi <- log(7)
cap_lo <- log(0.25)
heat_data <- heat_data %>%
  mutate(log_rr_capped = pmax(pmin(log_rr, cap_hi), cap_lo))

fig <- ggplot(heat_data, aes(x = ethnicity, y = label, fill = log_rr_capped)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = annot, colour = text_colour, fontface = fontface_val),
            size = 3.2, show.legend = FALSE) +
  scale_colour_identity() +
  scale_fill_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
    limits = c(cap_lo, cap_hi),
    name = "Risk Ratio",
    breaks = log(c(0.25, 0.5, 1, 2, 4, 7)),
    labels = sprintf("%.2f", c(0.25, 0.5, 1, 2, 4, 7))
  ) +
  annotate("rect",
           xmin = 0.5, xmax = 1.5,
           ymin = 0.5, ymax = n_distinct(heat_data$label) + 0.5,
           fill = NA, colour = gt_colour, linewidth = 1.5) +
  labs(x = NULL, y = NULL,
       caption = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.key.height = unit(1.1, "cm"),
    legend.key.width = unit(0.35, "cm"),
    plot.caption = element_text(size = 8, colour = "#444444", hjust = 0),
    plot.margin = margin(6, 6, 6, 6)
  )

ggsave(file.path(base_dir, "fig_heatmap.png"), fig,
       width = 6.5, height = 8.5, dpi = 300)
ggsave(file.path(base_dir, "fig_heatmap.svg"), fig,
       width = 6.5, height = 8.5)
ggsave(file.path(base_dir, "fig_heatmap.pdf"), fig,
       width = 6.5, height = 8.5)
cat("Heatmap saved.\n")


# ************************************************************
# ************************************************************
#
# SECTION 6: TABLE S6
# Source: Sophie_adjusted_cc_v0.1.R
#
# Table S6 -- Adjusted risk ratios in a complete case sample
# Model: eth_8cat + gender + grade + fas (no imputation)
#
# ************************************************************
# ************************************************************

# ============================================================
# SCRIPT: Adjusted Poisson Models — Complete Case Dataset
# ============================================================
# - Exposure: eth_8cat (categorical, ref = "Gypsy/Traveller")
# - Outcomes: 16 binary outcomes grouped by domain
# - Adjusted for gender + grade + fas
# - Cluster-robust SEs by id2c (school cluster)
# - Output: Word document with adjusted results table
# ============================================================

# ── 1. Setup ──
library(dplyr)
library(clubSandwich)
library(officer)
library(flextable)

OUTDIR <- "C:/Users/wppjw/OneDrive - Cardiff University/Papers/SHRN/Gyspy Irish Roma/Data/Complete case"
if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

CC_PATH <- file.path(OUTDIR, "GRT_2023_2025_complete_case.rds")

# --- Outcomes ---
outcomes <- c("lifesatneg", "sdqinternalhigh", "sdqexternalhigh",
              "cannabis30d", "alcoholuse", "smokeweekly", "regvape",
              "phys7", "vigexercise", "sit",
              "dailysoft", "dailyenergy",
              "beenbulliedever", "cbeenbulliedever",
              "truantever", "excludedever")

outcome_labels <- c(
  "lifesatneg"       = "Low life satisfaction",
  "sdqexternalhigh"  = "Externalizing difficulties\u1d43",
  "sdqinternalhigh"  = "Internalizing difficulties\u1d47",
  "smokeweekly"      = "Smoke weekly",
  "regvape"          = "Vape weekly",
  "alcoholuse"       = "\u22652 alcoholic drinks a session",
  "cannabis30d"      = "Cannabis in past 30 days",
  "phys7"            = "60 minutes of MVPA per day",
  "vigexercise"      = "\u22653 days VPA per week\u1d9c",
  "sit"              = "Sedentary behaviour",
  "dailysoft"        = "Daily soft drink consumption",
  "dailyenergy"      = "Daily energy drink consumption",
  "beenbulliedever"  = "Ever been bullied",
  "cbeenbulliedever" = "Ever been cyberbullied",
  "truantever"       = "Ever truanted",
  "excludedever"     = "Ever excluded from school"
)

outcome_groups <- list(
  "Mental health and life satisfaction" = c("lifesatneg", "sdqexternalhigh", "sdqinternalhigh"),
  "Substance use"     = c("smokeweekly", "regvape", "alcoholuse", "cannabis30d"),
  "Physical activity"  = c("phys7", "vigexercise", "sit"),
  "Diet"              = c("dailysoft", "dailyenergy"),
  "Bullying"          = c("beenbulliedever", "cbeenbulliedever"),
  "School problems"   = c("truantever", "excludedever")
)

eth_labels <- c(
  "1" = "White",
  "2" = "White Other",
  "3" = "Gypsy/Traveller",
  "4" = "Roma",
  "5" = "Black",
  "6" = "Mixed",
  "7" = "Asian",
  "8" = "Other"
)

eth_order <- c("Gypsy/Traveller", "Roma", "White Other",
               "Black", "Mixed", "Asian", "Other")

# ══════════════════════════════════════════════════════════════
# 2. LOAD DATA
# ══════════════════════════════════════════════════════════════
cat("Loading complete case data...\n")
t0 <- Sys.time()
shrndata <- readRDS(CC_PATH)

n_obs <- nrow(shrndata)
cat(sprintf("Complete case dataset: n = %d\n", n_obs))

# --- Prepare data ---
keep_cols <- unique(c("eth_8cat", "gender", "grade", "fas", "id2c", outcomes))
shrndata <- shrndata[, keep_cols, drop = FALSE]

shrndata$eth_8cat <- factor(shrndata$eth_8cat, levels = names(eth_labels),
                            labels = eth_labels)
shrndata$eth_8cat <- relevel(shrndata$eth_8cat, ref = "White")
shrndata$gender   <- as.factor(shrndata$gender)
shrndata$grade    <- as.factor(shrndata$grade)
shrndata$fas      <- as.numeric(shrndata$fas)
shrndata$id2c     <- as.integer(as.factor(shrndata$id2c))

for (oc in outcomes) {
  if (is.factor(shrndata[[oc]])) shrndata[[oc]] <- as.integer(shrndata[[oc]]) - 1L
}

# N per ethnicity
eth_tab <- table(shrndata$eth_8cat, useNA = "no")
eth_n <- as.data.frame(eth_tab)
names(eth_n) <- c("Ethnicity", "n")
eth_n$pct <- round(100 * eth_n$n / sum(eth_n$n), 1)
eth_n$n_pct <- paste0(format(eth_n$n, big.mark = ","), " (", eth_n$pct, "%)")
eth_n$Ethnicity <- factor(eth_n$Ethnicity, levels = eth_order)
eth_n <- eth_n[order(eth_n$Ethnicity), ]

# ══════════════════════════════════════════════════════════════
# 3. ADJUSTED POISSON REGRESSION
# ══════════════════════════════════════════════════════════════
cat("\n====== ADJUSTED MODELS ======\n")
cat("Formula: outcome ~ eth_8cat + gender + grade + fas\n")
cat("Reference category: White\n\n")

formula_rhs <- "eth_8cat + gender + grade + fas"

all_res <- list()
for (idx in seq_along(outcomes)) {
  out <- outcomes[idx]
  cat(sprintf("  [%d/%d] Model: %s ...", idx, length(outcomes), out))

  fml <- as.formula(paste0(out, " ~ ", formula_rhs))
  mod <- glm(fml, data = shrndata, family = poisson(link = "log"))
  cf  <- coef(mod)
  vc  <- clubSandwich::vcovCR(mod, cluster = shrndata$id2c, type = "CR1S")
  se  <- sqrt(diag(vc))

  res <- data.frame(
    term    = names(cf),
    est     = cf,
    se      = se,
    RR      = exp(cf),
    L       = exp(cf - 1.96 * se),
    H       = exp(cf + 1.96 * se),
    p       = 2 * pnorm(-abs(cf / se)),
    Outcome = out,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  all_res[[out]] <- res
  cat(" done\n")
  rm(mod)
}

results <- bind_rows(all_res)

eth_res <- results %>%
  filter(grepl("eth_8cat", term)) %>%
  mutate(
    term = gsub("eth_8cat", "", term),
    sig  = ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
  )

# Print results to console
cat("\n-- Adjusted: Ethnicity RRs (ref = Gypsy/Traveller) --\n")
for (out in outcomes) {
  cat(sprintf("\n--- %s ---\n", out))
  d <- eth_res %>% filter(Outcome == out)
  for (i in 1:nrow(d)) {
    cat(sprintf("  %-20s  RR = %5.2f (%4.2f\u2013%5.2f)  p = %.4f %s\n",
                d$term[i], d$RR[i], d$L[i], d$H[i], d$p[i], d$sig[i]))
  }
}

# ══════════════════════════════════════════════════════════════
# 4. BUILD WIDE TABLE & WORD DOCUMENT
# ══════════════════════════════════════════════════════════════
cat("\n====== CREATING WORD DOCUMENT ======\n")

format_rr <- function(rr, lo, hi) {
  sprintf("%.2f (%.2f, %.2f)", rr, lo, hi)
}

ref_group <- "White"
col_names <- c("Outcome", eth_order)
rows_list <- list()
row_types <- character()

# N row
n_row_vals <- c(paste0("N                                                        ",
                        eth_n$n_pct[match(ref_group, eth_n$Ethnicity)]))
for (eth in eth_order) {
  if (eth == ref_group) {
    n_row_vals <- c(n_row_vals, "")
  } else {
    n_row_vals <- c(n_row_vals, eth_n$n_pct[match(eth, eth_n$Ethnicity)])
  }
}
rows_list[[1]] <- setNames(as.data.frame(t(n_row_vals), stringsAsFactors = FALSE), col_names)
row_types <- c(row_types, "n")

for (section_name in names(outcome_groups)) {
  sec_vals <- c(section_name, rep("", length(eth_order)))
  rows_list[[length(rows_list) + 1]] <- setNames(
    as.data.frame(t(sec_vals), stringsAsFactors = FALSE), col_names)
  row_types <- c(row_types, "section")

  for (oc in outcome_groups[[section_name]]) {
    label <- outcome_labels[oc]
    sub <- eth_res %>% filter(Outcome == oc)

    vals <- c(label)
    for (eth in eth_order) {
      if (eth == ref_group) {
        vals <- c(vals, "1.00 (ref)")
      } else {
        eth_row <- sub %>% filter(term == eth)
        if (nrow(eth_row) == 1) {
          vals <- c(vals, format_rr(eth_row$RR, eth_row$L, eth_row$H))
        } else {
          vals <- c(vals, "")
        }
      }
    }
    rows_list[[length(rows_list) + 1]] <- setNames(
      as.data.frame(t(vals), stringsAsFactors = FALSE), col_names)
    row_types <- c(row_types, "data")
  }
}

wide_df <- bind_rows(rows_list)

# --- Format flextable ---
n_row_idx   <- which(row_types == "n")
section_idx <- which(row_types == "section")
data_idx    <- which(row_types == "data")
ncols       <- ncol(wide_df)

ft <- flextable(wide_df)

ft <- width(ft, j = 1, width = 2.3)
for (j in 2:ncols) ft <- width(ft, j = j, width = 1.15)

ft <- font(ft, fontname = "Aptos", part = "all")
ft <- fontsize(ft, size = 11, part = "header")
ft <- fontsize(ft, size = 10, part = "body")
ft <- bold(ft, part = "header")

ft <- fontsize(ft, i = data_idx, j = 1, size = 11)

if (length(n_row_idx) > 0) {
  ft <- fontsize(ft, i = n_row_idx, size = 11)
  ft <- bold(ft, i = n_row_idx, j = 1)
}

if (length(section_idx) > 0) {
  ft <- fontsize(ft, i = section_idx, size = 11)
  ft <- bold(ft, i = section_idx, j = 1)
  ft <- italic(ft, i = section_idx, j = 1)
}

ft <- align(ft, j = 1, align = "left", part = "all")
if (ncols > 1) ft <- align(ft, j = 2:ncols, align = "center", part = "all")

ft <- padding(ft, padding.top = 2, padding.bottom = 2,
              padding.left = 3, padding.right = 3, part = "all")

ft <- border_remove(ft)
ft <- hline_top(ft, part = "header",
                border = fp_border(color = "black", width = 1.5))
ft <- hline_bottom(ft, part = "header",
                   border = fp_border(color = "black", width = 1))
if (length(n_row_idx) > 0) {
  ft <- hline(ft, i = n_row_idx, border = fp_border(color = "black", width = 0.5),
              part = "body")
}
ft <- hline_bottom(ft, part = "body",
                   border = fp_border(color = "black", width = 1.5))

if (length(n_row_idx) > 0) {
  ft <- merge_at(ft, i = n_row_idx, j = 1:2, part = "body")
}
for (si in section_idx) {
  ft <- merge_at(ft, i = si, j = 1:ncols, part = "body")
}

# --- Assemble Word document ---
doc <- read_docx()
doc <- body_end_section_landscape(doc)
doc <- body_add_par(doc,
  "Table. Adjusted risk ratios (95% confidence interval) for all outcomes by ethnic group (ref: White) \u2014 Complete case analysis",
  style = "Normal")
doc <- body_add_flextable(doc, ft)
doc <- body_end_section_landscape(doc)

doc_path <- file.path(OUTDIR, "GRT_2023_2025_adjusted_results_cc.docx")
print(doc, target = doc_path)
cat(sprintf("Saved: %s\n", doc_path))

# ══════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════
cat(sprintf("\nTotal script time: %.1f seconds.\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
cat(sprintf("Output: %s\n", doc_path))
cat("Script completed.\n")


# ************************************************************
# ************************************************************
#
# SECTION 7: TABLE S9 (2023 data) and TABLE S10 (2025 data)
# Source: Sophie_adjusted_by_year_v0.1.R
#
# Table S9  -- Adjusted risk ratios using 2023 survey data
# Table S10 -- Adjusted risk ratios using 2025 survey data
#
# ************************************************************
# ************************************************************

# ============================================================
# SCRIPT: Adjusted Poisson Models by Survey Year (2023 & 2025)
# ============================================================
# - Exposure: eth_8cat (categorical, ref = "White")
# - Outcomes: 16 binary outcomes grouped by domain
# - Model: Adjusted for gender + grade + fas
# - Cluster-robust SEs by id2c (school cluster)
# - Stratified by syear (2023 and 2025)
# - Output: One combined CSV + Word document with two tables
# ============================================================

# ── 1. Setup ──
user_lib <- file.path(Sys.getenv("HOME"), "R_libs")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))

pkgs <- c("parallel", "dplyr", "clubSandwich", "mice", "tidyr",
          "officer", "flextable")
new_pkgs <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(new_pkgs) > 0) {
  install.packages(new_pkgs, lib = user_lib, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(parallel)
  library(dplyr)
  library(clubSandwich)
  library(mice)
  library(tidyr)
  library(officer)
  library(flextable)
})

# --- Configuration ---
N_CORES <- 3L
BASEDIR <- "/shared/home1/c.wppjw/SHRN/Gypsy RT"
OUTDIR  <- "/shared/home1/c.wppjw/SHRN/Gypsy RT/output"
if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

IMP_PATH <- file.path(BASEDIR, "output/v0.4/GRT_2023_2025_imp_m44_parallel.rds")

# --- Outcomes ---
outcomes <- c("lifesatneg", "sdqinternalhigh", "sdqexternalhigh",
              "cannabis30d", "alcoholuse", "smokeweekly", "regvape",
              "phys7", "vigexercise", "sit",
              "dailysoft", "dailyenergy",
              "beenbulliedever", "cbeenbulliedever",
              "truantever", "excludedever")

# --- Outcome labels ---
outcome_labels <- c(
  "lifesatneg"       = "Low life satisfaction",
  "sdqexternalhigh"  = "Externalizing difficulties\u1d43",
  "sdqinternalhigh"  = "Internalizing difficulties\u1d47",
  "smokeweekly"      = "Smoke weekly",
  "regvape"          = "Vape weekly",
  "alcoholuse"       = "\u22652 alcoholic drinks a session",
  "cannabis30d"      = "Cannabis in past 30 days",
  "phys7"            = "60 minutes of MVPA per day",
  "vigexercise"      = "\u22653 days VPA per week\u1d9c",
  "sit"              = "Sedentary behaviour",
  "dailysoft"        = "Daily soft drink consumption",
  "dailyenergy"      = "Daily energy drink consumption",
  "beenbulliedever"  = "Ever been bullied",
  "cbeenbulliedever" = "Ever been cyberbullied",
  "truantever"       = "Ever truanted",
  "excludedever"     = "Ever excluded from school"
)

# --- Outcome groups ---
outcome_groups <- list(
  "Mental health and life satisfaction" = c("lifesatneg", "sdqexternalhigh", "sdqinternalhigh"),
  "Substance use"     = c("smokeweekly", "regvape", "alcoholuse", "cannabis30d"),
  "Physical activity"  = c("phys7", "vigexercise", "sit"),
  "Diet"              = c("dailysoft", "dailyenergy"),
  "Bullying"          = c("beenbulliedever", "cbeenbulliedever"),
  "School problems"   = c("truantever", "excludedever")
)

# --- Ethnicity labels ---
eth_labels <- c(
  "1" = "White",
  "2" = "White Other",
  "3" = "Gypsy/Traveller",
  "4" = "Roma",
  "5" = "Black",
  "6" = "Mixed",
  "7" = "Asian",
  "8" = "Other"
)

eth_order <- c("White", "White Other", "Gypsy/Traveller", "Roma",
               "Black", "Mixed", "Asian", "Other")

# ══════════════════════════════════════════════════════════════
# 2. LOAD DATA
# ══════════════════════════════════════════════════════════════
cat("Loading imputed data...\n")
t0 <- Sys.time()
imp <- readRDS(IMP_PATH)

n_imp <- imp$m
n_obs <- nrow(imp$data)

cat(sprintf("Loaded mids object: %d imputations, n = %d (%.1f sec)\n",
            n_imp, n_obs,
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
cat(sprintf("  mids object size: %.1f GB\n", object.size(imp) / 1e9))
cat(sprintf("  Current R memory usage: %.1f GB\n\n", sum(gc()[, 2]) / 1024))

# ── Prepare datasets ──
keep_cols <- unique(c("eth_8cat", "gender", "grade", "fas", "id2c", "syear", outcomes))

prepare_df <- function(df) {
  df <- df[, keep_cols, drop = FALSE]
  df$eth_8cat <- factor(df$eth_8cat, levels = names(eth_labels),
                        labels = eth_labels)
  df$eth_8cat <- relevel(df$eth_8cat, ref = "White")
  df$gender   <- as.factor(df$gender)
  df$grade    <- as.factor(df$grade)
  df$fas      <- as.numeric(df$fas)
  df$id2c     <- as.integer(as.factor(df$id2c))
  for (oc in outcomes) {
    if (is.factor(df[[oc]])) df[[oc]] <- as.integer(df[[oc]]) - 1L
  }
  df
}

# Extract all imputed datasets
cat("Extracting all imputed datasets...\n")
imp_list <- lapply(seq_len(n_imp), function(i) {
  prepare_df(mice::complete(imp, action = i))
})
cat(sprintf("  Extracted %d datasets (%.1f GB total)\n",
            n_imp, object.size(imp_list) / 1e9))

# Free the mids object
rm(imp); gc(verbose = FALSE)
cat(sprintf("  Memory after freeing mids object: %.1f GB\n\n", sum(gc()[, 2]) / 1024))

# ── Diagnostics ──
cat("=== SYSTEM DIAGNOSTICS ===\n")
cat(sprintf("  R version: %s\n", R.version.string))
cat(sprintf("  Platform: %s\n", R.version$platform))
cat(sprintf("  N cores requested: %d\n", N_CORES))
cat(sprintf("  N cores detected: %d\n", parallel::detectCores()))
mem_info <- tryCatch(system("free -h", intern = TRUE), error = function(e) NULL)
if (!is.null(mem_info)) {
  cat("  System memory:\n")
  for (line in mem_info) cat("    ", line, "\n")
}
slurm_mem <- Sys.getenv("SLURM_MEM_PER_NODE", unset = NA)
if (!is.na(slurm_mem)) cat(sprintf("  SLURM allocated memory: %s MB\n", slurm_mem))
cat(sprintf("  R process RSS: %.1f GB\n", sum(gc()[, 2]) / 1024))
cat(sprintf("  imp_list object size: %.2f GB\n", object.size(imp_list) / 1e9))
cat(sprintf("  Single dataset size: %.1f MB\n", object.size(imp_list[[1]]) / 1e6))
cat(sprintf("  Estimated peak memory (parent + %d forks): %.1f GB\n",
            N_CORES,
            (sum(gc()[, 2]) / 1024) * (1 + N_CORES)))
cat("===========================\n\n")

# ══════════════════════════════════════════════════════════════
# 3. POISSON REGRESSION — ADJUSTED RRs BY SURVEY YEAR
# ══════════════════════════════════════════════════════════════

fit_model <- function(df, outcome_var, formula_rhs) {
  tryCatch({
    fml <- as.formula(paste0(outcome_var, " ~ ", formula_rhs))
    mod <- glm(fml, data = df, family = poisson(link = "log"))
    cf  <- coef(mod)
    vc  <- clubSandwich::vcovCR(mod, cluster = df$id2c, type = "CR1S")
    rm(mod)
    return(list(coef = cf, vcov = vc, error = NULL))
  }, error = function(e) return(list(error = e$message)))
}

pool_rubin <- function(fit_list) {
  results_only <- Filter(function(x) !is.null(x) && is.null(x$error), fit_list)
  m <- length(results_only)
  if (m == 0) stop("All models failed (check for OOM kills or model errors).")
  if (m < 3) warning(sprintf("Only %d/%d imputations succeeded — estimates unreliable.", m, length(fit_list)))

  q_m <- do.call(rbind, lapply(results_only, `[[`, "coef"))
  u_m <- lapply(results_only, `[[`, "vcov")

  q_bar <- colMeans(q_m)
  u_bar <- Reduce("+", u_m) / m
  b_m   <- cov(q_m)
  total_vcov <- u_bar + (1 + 1/m) * b_m

  se <- sqrt(diag(total_vcov))

  r <- (1 + 1/m) * diag(b_m) / diag(u_bar)
  v_old <- (m - 1) * (1 + 1/r)^2
  v_old[is.nan(v_old) | is.na(v_old)] <- Inf

  # Monte Carlo error diagnostics
  mc_error_est <- sqrt(diag(b_m) / m)
  mc_error_se  <- se / sqrt(2 * (m - 1))
  total_var    <- diag(total_vcov)
  fmi          <- (diag(b_m) + diag(b_m) / m) / total_var
  fmi[is.nan(fmi) | is.na(fmi)] <- 0

  mc_ratio     <- mc_error_est / se
  mc_ratio[is.nan(mc_ratio) | is.na(mc_ratio)] <- 0
  mc_adequate  <- ifelse(mc_ratio < 0.10, "OK", "INCREASE m")

  data.frame(
    term = names(q_bar),
    est  = q_bar,
    se   = se,
    RR   = exp(q_bar),
    L    = exp(q_bar - 1.96 * se),
    H    = exp(q_bar + 1.96 * se),
    p    = 2 * pt(-abs(q_bar / se), df = v_old),
    fmi  = round(fmi, 3),
    mc_error_est = round(mc_error_est, 5),
    mc_error_se  = round(mc_error_se, 5),
    mc_pct_of_se = round(100 * mc_ratio, 1),
    mc_adequate  = mc_adequate,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

run_adjusted_models <- function(imp_list, outcomes, year_label) {
  cat(sprintf("\n====== ADJUSTED MODELS — %s ======\n", year_label))
  formula_rhs <- "eth_8cat + gender + grade + fas"
  cat(sprintf("Formula: outcome ~ %s\n", formula_rhs))
  cat(sprintf("Using %d cores for parallel processing\n\n", N_CORES))

  # Subset imputed datasets to this survey year
  imp_year <- lapply(imp_list, function(df) df[df$syear == year_label, , drop = FALSE])
  n_obs_year <- nrow(imp_year[[1]])
  cat(sprintf("  N observations for %s: %d\n\n", year_label, n_obs_year))

  all_res <- list()
  t_models <- Sys.time()
  for (idx in seq_along(outcomes)) {
    out <- outcomes[idx]
    t_out <- Sys.time()
    cat(sprintf("  [%d/%d] Model: %s ...", idx, length(outcomes), out))
    fits <- mclapply(imp_year, fit_model, outcome_var = out,
                     formula_rhs = formula_rhs,
                     mc.cores = N_CORES, mc.preschedule = FALSE)
    n_null <- sum(sapply(fits, is.null))
    n_err  <- sum(sapply(fits, function(x) !is.null(x$error)))
    n_ok   <- length(fits) - n_null - n_err
    elapsed_out <- as.numeric(difftime(Sys.time(), t_out, units = "secs"))
    elapsed_total <- as.numeric(difftime(Sys.time(), t_models, units = "secs"))
    avg_per <- elapsed_total / idx
    remaining <- avg_per * (length(outcomes) - idx)
    cat(sprintf(" %d/%d converged (%.1fs, ~%.0fs remaining)\n",
                n_ok, length(imp_year), elapsed_out, remaining))
    if (n_null > 0) {
      cat(sprintf("  *** WARNING: %d workers returned NULL (likely OOM killed) ***\n", n_null))
      cat(sprintf("  *** R memory: %.1f GB | Consider increasing --mem ***\n", sum(gc()[, 2]) / 1024))
    }
    if (n_err > 0) {
      errs <- sapply(Filter(function(x) !is.null(x$error), fits), `[[`, "error")
      cat(sprintf("  *** %d model errors: %s ***\n", n_err, paste(unique(errs), collapse = "; ")))
    }
    all_res[[out]] <- pool_rubin(fits) %>% mutate(Outcome = out)
    rm(fits); gc(verbose = FALSE)
  }
  cat(sprintf("  Adjusted models for %s completed in %.1f seconds.\n\n",
              year_label, as.numeric(difftime(Sys.time(), t_models, units = "secs"))))

  results <- bind_rows(all_res)

  eth_res <- results %>%
    filter(grepl("eth_8cat", term)) %>%
    mutate(
      term  = gsub("eth_8cat", "", term),
      sig   = ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))),
      Model = "Adjusted",
      Year  = year_label
    )

  cat(sprintf("\n-- %s Adjusted: Ethnicity RRs (ref = White) --\n", year_label))
  for (out in outcomes) {
    cat(sprintf("\n--- %s ---\n", out))
    d <- eth_res %>% filter(Outcome == out)
    for (i in 1:nrow(d)) {
      cat(sprintf("  %-20s  RR = %5.2f (%4.2f\u2013%5.2f)  p = %.4f %s\n",
                  d$term[i], d$RR[i], d$L[i], d$H[i], d$p[i], d$sig[i]))
    }
  }

  # Monte Carlo error summary
  cat(sprintf("\n-- %s Adjusted: Monte Carlo Error Diagnostics --\n", year_label))
  cat("  (MC error < 10%% of SE = adequate number of imputations)\n\n")
  mc_flags <- eth_res %>% filter(mc_adequate == "INCREASE m")
  if (nrow(mc_flags) == 0) {
    cat(sprintf("  ALL estimates have adequate MC error (< 10%% of SE). m = %d is sufficient.\n", length(imp_list)))
  } else {
    cat(sprintf("  WARNING: %d estimate(s) have MC error >= 10%% of SE:\n", nrow(mc_flags)))
    for (i in 1:nrow(mc_flags)) {
      cat(sprintf("    %-25s %-20s  MC error = %.1f%% of SE, FMI = %.3f\n",
                  mc_flags$Outcome[i], mc_flags$term[i],
                  mc_flags$mc_pct_of_se[i], mc_flags$fmi[i]))
    }
    cat("  Consider increasing m for these estimates.\n")
  }

  rm(imp_year); gc(verbose = FALSE)
  list(all_results = results, eth_results = eth_res)
}

# --- Run adjusted models separately by year ---
adj_2023 <- run_adjusted_models(imp_list, outcomes, year_label = "2023")
adj_2025 <- run_adjusted_models(imp_list, outcomes, year_label = "2025")

# ══════════════════════════════════════════════════════════════
# 4. SAVE CSV — Combined results for both years
# ══════════════════════════════════════════════════════════════
combined_eth <- bind_rows(adj_2023$eth_results, adj_2025$eth_results)
write.csv(combined_eth, file.path(OUTDIR, "Adjusted_Results_by_Year.csv"), row.names = FALSE)
cat("Saved: Adjusted_Results_by_Year.csv\n")

# ══════════════════════════════════════════════════════════════
# 5. WORD DOCUMENT — Two tables (one per year)
# ══════════════════════════════════════════════════════════════

build_wide_table <- function(eth_res, year_label) {
  # Create formatted RR (95% CI) strings
  eth_res$formatted <- sprintf("%.2f (%.2f\u2013%.2f)",
                               eth_res$RR, eth_res$L, eth_res$H)

  # Pivot to wide: rows = outcomes (grouped), columns = ethnic groups (excl. ref)
  eth_groups <- setdiff(eth_order, "White")

  rows <- list()
  for (grp_name in names(outcome_groups)) {
    # Group header row
    header_row <- data.frame(Outcome = grp_name, stringsAsFactors = FALSE)
    for (eg in eth_groups) header_row[[eg]] <- ""
    header_row$is_header <- TRUE
    rows[[length(rows) + 1]] <- header_row

    for (oc in outcome_groups[[grp_name]]) {
      data_row <- data.frame(Outcome = outcome_labels[oc], stringsAsFactors = FALSE)
      for (eg in eth_groups) {
        match_val <- eth_res %>%
          filter(Outcome == oc & term == eg) %>%
          pull(formatted)
        data_row[[eg]] <- if (length(match_val) > 0) match_val[1] else ""
      }
      data_row$is_header <- FALSE
      rows[[length(rows) + 1]] <- data_row
    }
  }

  wide <- bind_rows(rows)
  wide
}

format_flextable <- function(wide_df, year_label) {
  is_header <- wide_df$is_header
  wide_df$is_header <- NULL

  ft <- flextable(wide_df)
  ft <- set_header_labels(ft, Outcome = "Outcome")
  ft <- set_caption(ft, caption = sprintf(
    "Adjusted risk ratios (95%% CI) for health outcomes by ethnicity (ref = White), %s. Adjusted for gender, grade, and FAS.",
    year_label))

  # Bold group header rows
  for (i in seq_along(is_header)) {
    if (is_header[i]) {
      ft <- bold(ft, i = i, part = "body")
    }
  }

  ft <- bold(ft, part = "header")
  ft <- fontsize(ft, size = 10.5, part = "all")
  ft <- font(ft, fontname = "Aptos", part = "all")
  ft <- autofit(ft)
  ft <- set_table_properties(ft, layout = "autofit")

  ft
}

cat("Building Word document...\n")

wide_2023 <- build_wide_table(adj_2023$eth_results, "2023")
wide_2025 <- build_wide_table(adj_2025$eth_results, "2025")

ft_2023 <- format_flextable(wide_2023, "2023")
ft_2025 <- format_flextable(wide_2025, "2025")

doc <- read_docx()
doc <- body_add_par(doc, "Adjusted Risk Ratios by Ethnicity — Stratified by Survey Year",
                    style = "heading 1")
doc <- body_add_par(doc, "")

doc <- body_add_par(doc, "Table 1: Survey Year 2023", style = "heading 2")
doc <- body_add_flextable(doc, ft_2023)
doc <- body_add_break(doc)

doc <- body_add_par(doc, "Table 2: Survey Year 2025", style = "heading 2")
doc <- body_add_flextable(doc, ft_2025)

print(doc, target = file.path(OUTDIR, "Adjusted_RR_by_Year.docx"))
cat("Saved: Adjusted_RR_by_Year.docx\n")

# ══════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════
cat(sprintf("\n\nTotal script time: %.1f minutes.\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("\nAll outputs saved to:", OUTDIR, "\n")
cat("  1.  Adjusted_Results_by_Year.csv\n")
cat("  2.  Adjusted_RR_by_Year.docx\n")
cat("\nScript completed.\n")


# ************************************************************
# ************************************************************
#
# SECTION 8: TABLE S11 and TABLE S12
# Source: grt_stratified_by_ethnicity_v0.1.R
#
# Table S11 -- Adjusted risk ratios for the difference in
#              outcomes between Gypsy/Traveller young people
#              according to gender identity
# Table S12 -- Adjusted risk ratios for the difference in
#              outcomes between Roma young people according
#              to gender identity
#
# ************************************************************
# ************************************************************

#!/usr/bin/env Rscript
# ============================================================
# SCRIPT: Ethnicity-stratified gender comparison
# ------------------------------------------------------------
# Uses the subgroup mids object:
#   syear == 2025 + (syear == 2023 & grade in 10,11)
#
# For Gypsy/Traveller and Roma subsets separately, fits:
#   Adjusted:  outcome ~ gender + grade + fas
#
# Gender reference = Boy. The genderGirl coefficient gives the
# exact within-group girls-vs-boys RR and p-value (i.e. the
# direct comparison of GT girls vs GT boys, and Roma girls vs
# Roma boys).
#
# Cluster-robust SEs by id2c (school). Pooled across imputations
# via Rubin's rules.
#
# Outputs (all in output/v0.93/):
#   CSVs:
#     - Stratified_GT_Adjusted.csv
#     - Stratified_Roma_Adjusted.csv
#     - Ethnicity_Stratified_Gender_Combined.csv
#   Word:
#     - Ethnicity_Stratified_Gender_Tables.docx
# ============================================================

# -- Library path --
user_lib <- file.path(Sys.getenv("HOME"), "R_libs")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))

suppressPackageStartupMessages({
  library(parallel)
  library(dplyr)
  library(tidyr)
  library(clubSandwich)
  library(mice)
  library(officer)
  library(flextable)
})

# -- Configuration --
N_CORES  <- 4L
BASEDIR  <- "/shared/home1/c.wppjw/SHRN/Gypsy RT"
OUTDIR   <- file.path(BASEDIR, "output/v0.93")
IMP_PATH <- file.path(BASEDIR, "output/v0.91/GRT_subgroup_imp.rds")

if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)
if (!file.exists(IMP_PATH)) {
  stop("Subgroup mids not found at ", IMP_PATH,
       "\n  Run grt_subgroup_v0.91.R first.")
}

# -- Outcomes --
outcomes <- c("lifesatneg", "sdqinternalhigh", "sdqexternalhigh",
              "cannabis30d", "alcoholuse", "smokeweekly", "regvape",
              "phys7", "vigexercise", "sit",
              "dailysoft", "dailyenergy",
              "beenbulliedever", "cbeenbulliedever",
              "truantever", "excludedever")

outcome_labels <- c(
  "lifesatneg"       = "Low life satisfaction",
  "sdqexternalhigh"  = "Externalising difficulties",
  "sdqinternalhigh"  = "Internalising difficulties",
  "smokeweekly"      = "Weekly smoking",
  "regvape"          = "Weekly vaping",
  "alcoholuse"       = "\u22652 alcoholic drinks/session",
  "cannabis30d"      = "Cannabis past 30 days",
  "phys7"            = "<60 min MVPA per day",
  "vigexercise"      = "\u22653 days VPA per week",
  "sit"              = "Sedentary behaviour",
  "dailysoft"        = "Daily soft drink consumption",
  "dailyenergy"      = "Daily energy drink consumption",
  "beenbulliedever"  = "Ever been bullied",
  "cbeenbulliedever" = "Ever been cyberbullied",
  "truantever"       = "Ever truanted",
  "excludedever"     = "Ever excluded from school"
)

eth_labels <- c("1" = "White", "2" = "White Other", "3" = "Gypsy/Traveller",
                "4" = "Roma",  "5" = "Black",       "6" = "Mixed",
                "7" = "Asian", "8" = "Other")

# Ethnic groups to stratify on
target_eth <- c("GT" = "Gypsy/Traveller", "Roma" = "Roma")

# ============================================================
# 1. LOAD SUBGROUP MIDS
# ============================================================
cat("Loading subgroup mids from ", IMP_PATH, "\n", sep = "")
t0  <- Sys.time()
imp <- readRDS(IMP_PATH)
n_imp <- imp$m
cat(sprintf("  Loaded: %d imputations, n = %d (%.1f sec)\n",
            n_imp, nrow(imp$data),
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# ============================================================
# 2. EXTRACT AND PREPARE IMPUTED DATASETS
# ============================================================
keep_cols <- unique(c("eth_8cat", "gender", "grade", "fas", "id2c", outcomes))

prepare_df <- function(df) {
  df <- df[, keep_cols, drop = FALSE]
  df$eth_8cat <- factor(df$eth_8cat, levels = names(eth_labels),
                        labels = eth_labels)
  df$eth_8cat <- relevel(df$eth_8cat, ref = "White")

  df$gender <- factor(df$gender, levels = c("1", "2", "3"),
                      labels = c("Boy", "Girl", "Other"))
  df$gender <- relevel(df$gender, ref = "Boy")

  df$grade <- as.factor(df$grade)
  df$fas   <- as.numeric(df$fas)
  df$id2c  <- as.integer(as.factor(df$id2c))
  for (oc in outcomes) {
    if (is.factor(df[[oc]])) df[[oc]] <- as.integer(df[[oc]]) - 1L
  }
  df
}

cat("\nExtracting imputed datasets...\n")
imp_list <- lapply(seq_len(n_imp), function(i) prepare_df(mice::complete(imp, i)))
rm(imp); gc(verbose = FALSE)
cat(sprintf("  Extracted %d datasets (%.2f GB total)\n",
            n_imp, object.size(imp_list) / 1e9))

# ============================================================
# 3. CORE FUNCTIONS
# ============================================================
fit_model <- function(df, outcome_var, formula_rhs) {
  tryCatch({
    fml <- as.formula(paste0(outcome_var, " ~ ", formula_rhs))
    mod <- glm(fml, data = df, family = poisson(link = "log"))
    cf  <- coef(mod)
    vc  <- clubSandwich::vcovCR(mod, cluster = df$id2c, type = "CR1S")
    rm(mod)
    list(coef = cf, vcov = vc, error = NULL)
  }, error = function(e) list(error = e$message))
}

pool_rubin <- function(fit_list) {
  ok <- Filter(function(x) !is.null(x) && is.null(x$error), fit_list)
  m <- length(ok)
  if (m == 0) stop("All models failed.")
  if (m < 3) warning(sprintf("Only %d/%d imputations succeeded.", m, length(fit_list)))

  q_m   <- do.call(rbind, lapply(ok, `[[`, "coef"))
  u_m   <- lapply(ok, `[[`, "vcov")
  q_bar <- colMeans(q_m)
  u_bar <- Reduce("+", u_m) / m
  b_m   <- cov(q_m)
  total_vcov <- u_bar + (1 + 1/m) * b_m
  se    <- sqrt(diag(total_vcov))

  r     <- (1 + 1/m) * diag(b_m) / diag(u_bar)
  v_old <- (m - 1) * (1 + 1/r)^2
  v_old[is.nan(v_old) | is.na(v_old)] <- Inf

  data.frame(
    term = names(q_bar),
    estimate = q_bar,
    std.error = se,
    RR   = exp(q_bar),
    CI_lower = exp(q_bar - 1.96 * se),
    CI_upper = exp(q_bar + 1.96 * se),
    statistic = q_bar / se,
    df = v_old,
    p.value = 2 * pt(-abs(q_bar / se), df = v_old),
    n_imp_converged = m,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

run_models <- function(data_list, outcomes, formula_rhs, model_label) {
  cat(sprintf("\n====== %s ======\n", toupper(model_label)))
  cat(sprintf("Formula: outcome ~ %s\n", formula_rhs))
  cat(sprintf("n = %d per imputation\n\n", nrow(data_list[[1]])))

  all_res  <- list()
  t_models <- Sys.time()
  for (idx in seq_along(outcomes)) {
    out   <- outcomes[idx]
    t_out <- Sys.time()
    cat(sprintf("  [%d/%d] %s ...", idx, length(outcomes), out))
    fits  <- mclapply(data_list, fit_model,
                      outcome_var = out, formula_rhs = formula_rhs,
                      mc.cores = N_CORES, mc.preschedule = FALSE)
    n_ok  <- sum(sapply(fits, function(x) !is.null(x) && is.null(x$error)))
    elap  <- as.numeric(difftime(Sys.time(), t_out, units = "secs"))
    tot   <- as.numeric(difftime(Sys.time(), t_models, units = "secs"))
    rem   <- if (idx > 0) (tot / idx) * (length(outcomes) - idx) else 0
    cat(sprintf(" %d/%d converged (%.1fs, ~%.0fs remaining)\n",
                n_ok, length(fits), elap, rem))
    pooled <- pool_rubin(fits)
    pooled$outcome <- out
    all_res[[out]] <- pooled
    rm(fits); gc(verbose = FALSE)
  }
  cat(sprintf("\n  Completed in %.1f seconds.\n",
              as.numeric(difftime(Sys.time(), t_models, units = "secs"))))

  bind_rows(all_res)
}

# ============================================================
# 4. FIT STRATIFIED MODELS (one ethnicity at a time)
# ============================================================
cat("\n\n########################################\n")
cat("# ETHNICITY-STRATIFIED GENDER MODELS\n")
cat("########################################\n")

strat_results <- list()

for (eth_short in names(target_eth)) {
  eth_full <- target_eth[[eth_short]]
  cat(sprintf("\n--- Ethnicity stratum: %s ---\n", eth_full))

  eth_list <- lapply(imp_list, function(df) df[df$eth_8cat == eth_full, , drop = FALSE])
  cat(sprintf("  n per imputation: %d\n", nrow(eth_list[[1]])))

  # Adjusted
  label_a <- paste0(eth_short, " - Adjusted (grade + FAS)")
  res_a <- run_models(eth_list, outcomes, "gender + grade + fas", label_a)
  res_a$ethnicity_stratum <- eth_full
  res_a$model <- "Adjusted"
  strat_results[[paste0(eth_short, "_adj")]] <- res_a

  rm(eth_list); gc(verbose = FALSE)
}

# ============================================================
# 5. WRITE INDIVIDUAL CSVs
# ============================================================
cat("\n\nSaving CSVs...\n")
op <- options(scipen = 999, digits = 17)

for (nm in names(strat_results)) {
  res <- strat_results[[nm]]
  eth <- res$ethnicity_stratum[1]
  mod <- res$model[1]
  eth_short <- names(target_eth)[target_eth == eth]
  fname <- paste0("Stratified_", eth_short, "_", mod, ".csv")
  write.csv(res, file.path(OUTDIR, fname), row.names = FALSE)
  cat(sprintf("  %s (%d rows)\n", fname, nrow(res)))
}

# ============================================================
# 6. COMBINED CSV: girls-vs-boys contrast for GT and Roma
# ============================================================
format_gender_rr <- function(results, eth_label) {
  results %>%
    filter(term %in% c("genderGirl", "genderOther")) %>%
    mutate(
      ethnicity = eth_label,
      gender_contrast = ifelse(term == "genderGirl",
                               "Girls vs Boys", "Other vs Boys"),
      outcome_label = outcome_labels[outcome],
      rr_ci = sprintf("%.2f (%.2f\u2013%.2f)", RR, CI_lower, CI_upper),
      p_fmt = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value)),
      p_exact = format.pval(p.value, digits = 4, eps = .Machine$double.eps)
    ) %>%
    select(ethnicity, model, outcome, outcome_label, gender_contrast,
           estimate, std.error, RR, CI_lower, CI_upper,
           statistic, df, p.value, p_exact, rr_ci, p_fmt,
           n_imp_converged)
}

combined <- bind_rows(
  format_gender_rr(strat_results$GT_adj,   "Gypsy/Traveller"),
  format_gender_rr(strat_results$Roma_adj, "Roma")
)

write.csv(combined,
          file.path(OUTDIR, "Ethnicity_Stratified_Gender_Combined.csv"),
          row.names = FALSE)
cat(sprintf("  Ethnicity_Stratified_Gender_Combined.csv (%d rows)\n",
            nrow(combined)))

options(op)

# ============================================================
# 7. BUILD WORD TABLES
# ============================================================
cat("\nBuilding Word tables...\n")

build_table <- function(combined, eth_label, model_type) {
  combined %>%
    filter(ethnicity == eth_label, model == model_type) %>%
    select(Outcome = outcome_label, Contrast = gender_contrast,
           `RR (95% CI)` = rr_ci, `p-value` = p_fmt,
           `Exact p` = p_exact) %>%
    arrange(Outcome, Contrast)
}

doc <- read_docx()

section_titles <- list(
  c("Gypsy/Traveller", "Adjusted",
    "Table 1. Gender contrasts within Gypsy/Traveller pupils (adjusted for grade + FAS)",
    paste0("Modified Poisson regression fit on Gypsy/Traveller pupils only, ",
           "adjusted for grade and family affluence (FAS). ",
           "Reference: Gypsy/Traveller boys. ",
           "\u2018Girls vs Boys\u2019 row = GT girls vs GT boys (direct within-group contrast). ",
           "Cluster-robust SEs by school. Pooled across imputations using Rubin\u2019s rules.")),
  c("Roma", "Adjusted",
    "Table 2. Gender contrasts within Roma pupils (adjusted for grade + FAS)",
    paste0("Modified Poisson regression fit on Roma pupils only, ",
           "adjusted for grade and family affluence (FAS). ",
           "Reference: Roma boys. ",
           "\u2018Girls vs Boys\u2019 row = Roma girls vs Roma boys (direct within-group contrast). ",
           "Cluster-robust SEs by school. Pooled across imputations using Rubin\u2019s rules."))
)

for (i in seq_along(section_titles)) {
  s <- section_titles[[i]]
  eth_label <- s[1]; model_type <- s[2]; title <- s[3]; caption <- s[4]

  tab <- build_table(combined, eth_label, model_type)

  doc <- body_add_par(doc, title, style = "heading 1")
  doc <- body_add_par(doc, caption, style = "Normal")

  ft <- flextable(tab) %>%
    theme_vanilla() %>%
    autofit() %>%
    font(fontname = "Arial", part = "all") %>%
    fontsize(size = 9, part = "body") %>%
    fontsize(size = 10, part = "header") %>%
    bold(part = "header") %>%
    merge_v(j = "Outcome") %>%
    valign(j = "Outcome", valign = "top")

  doc <- body_add_flextable(doc, ft)
  if (i < length(section_titles)) doc <- body_add_break(doc)
}

doc_path <- file.path(OUTDIR, "Ethnicity_Stratified_Gender_Tables.docx")
print(doc, target = doc_path)
cat(sprintf("\nSaved: %s\n", doc_path))

# ============================================================
cat("\n====== ALL ANALYSES COMPLETE ======\n")
cat(sprintf("All outputs under: %s\n", OUTDIR))
cat("\nCSVs:\n")
for (eth_short in names(target_eth)) {
  cat(sprintf("  - Stratified_%s_Adjusted.csv\n", eth_short))
}
cat("  - Ethnicity_Stratified_Gender_Combined.csv\n")
cat("\nWord:\n")
cat("  - Ethnicity_Stratified_Gender_Tables.docx (2 tables)\n")
cat("    T1: GT adjusted gender contrasts\n")
cat("    T2: Roma adjusted gender contrasts\n")


# ************************************************************
# ************************************************************
#
# SECTION 9: FIGURE 2
# Source: Paper/Figures 10.04.2026/redo_attenuation_ladder.R
#
# Figure 2 -- Unadjusted and adjusted risk ratios for health
#             and behavioural outcomes for Gypsy/Traveller and
#             Roma young people (attenuation ladder plot)
#
# Reads: Results/v0.7/Percent_Attenuation_GT_Roma.csv
#        (produced by Section 10 below)
#
# ************************************************************
# ************************************************************

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(patchwork)

results_dir <- file.path("C:/Users/wppjw",
                          "OneDrive - Cardiff University", "Papers", "SHRN",
                          "Gyspy Irish Roma", "Results")
out_dir <- file.path("C:/Users/wppjw",
                      "OneDrive - Cardiff University", "Papers", "SHRN",
                      "Gyspy Irish Roma", "Paper", "Figures 10.04.2026")

attenuation <- read_csv(file.path(results_dir, "v0.7", "Percent_Attenuation_GT_Roma.csv"),
                        show_col_types = FALSE)

# Clean outcome names
attenuation <- attenuation %>%
  mutate(Outcome = gsub("\u1D43|\u1D47|\u1D9C", "", Outcome),
         Outcome = trimws(Outcome))

att_long <- attenuation %>%
  rename(
    Unadjusted = RR_Unadj,
    `+ Gender` = RR_Gender,
    `+ School year` = RR_Grade,
    `+ Gender + Year` = RR_GenderGrade,
    `Fully adjusted` = RR_Full
  ) %>%
  pivot_longer(
    cols = c(Unadjusted, `+ Gender`, `+ School year`, `+ Gender + Year`, `Fully adjusted`),
    names_to = "model",
    values_to = "RR"
  ) %>%
  mutate(
    model = factor(model, levels = c("Unadjusted", "+ Gender", "+ School year",
                                     "+ Gender + Year", "Fully adjusted")),
    model_num = as.numeric(model)
  )

make_attenuation <- function(data, eth, panel_label) {
  d <- data %>% filter(Ethnicity == eth)

  ggplot(d, aes(x = model_num, y = RR, group = Outcome, colour = Outcome)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    geom_line(linewidth = 0.55, alpha = 0.8) +
    geom_point(size = 1.5) +
    scale_x_continuous(
      breaks = 1:5,
      labels = c("Unadjusted", "+ Gender", "+ School\nyear", "+ Gender\n+ Year", "Fully\nadjusted")
    ) +
    scale_y_log10(breaks = c(0.5, 0.75, 1, 1.5, 2, 3, 5, 7, 10, 13),
                  labels = c("0.5", "0.75", "1", "1.5", "2", "3", "5", "7", "10", "13")) +
    labs(x = NULL, y = "Risk ratio (log scale)",
         subtitle = panel_label, colour = "Outcome") +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 6.5),
      legend.title = element_text(size = 7.5, face = "bold"),
      legend.key.size = unit(0.35, "cm"),
      legend.key.spacing.y = unit(1, "pt"),
      plot.subtitle = element_text(face = "bold", size = 9.5),
      axis.text.x = element_text(size = 7),
      axis.text.y = element_text(size = 7),
      axis.title.y = element_text(size = 8),
      plot.margin = margin(t = 2, r = 2, b = 2, l = 4)
    ) +
    guides(colour = guide_legend(ncol = 1))
}

fig_gt   <- make_attenuation(att_long, "Gypsy/Traveller", "a  Gypsy/Traveller")
fig_roma <- make_attenuation(att_long, "Roma", "b  Roma")

fig <- fig_gt / fig_roma + plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave(file.path(out_dir, "Figure_attenuation.png"), fig,
       width = 9.5, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "Figure_attenuation.pdf"), fig,
       width = 9.5, height = 5.5, bg = "white")
cat("Figure (attenuation ladder) saved.\n")


# ************************************************************
# ************************************************************
#
# SECTION 10: TABLE S7 and TABLE S8
# Source: Results/v0.7/percent_attenuation.R
#     and Results/v0.7/percent_attenuation_from_docx.R
#
# Table S7 -- Percentage attenuation of relative risks across
#             adjustment models for Gypsy/Traveller
# Table S8 -- Percentage attenuation of relative risks across
#             adjustment models for Roma
#
# Two scripts contribute:
#   percent_attenuation.R -- computes attenuation from CSVs
#     of all ethnic groups (full table)
#   percent_attenuation_from_docx.R -- extracts RRs from Word
#     tables and computes GT/Roma-specific attenuation
#     (produces Percent_Attenuation_GT_Roma.csv used by Fig 2)
#
# Both are included below.
#
# ************************************************************
# ************************************************************

# ---------- percent_attenuation.R ----------

library(tidyverse)
library(flextable)
library(officer)

# ── 1. Read all files ──────────────────────────────────────────────────────

unadj_raw <- read.csv("Unadjusted_Results.csv", stringsAsFactors = FALSE)
gender_raw <- read.csv("Adjusted_Gender_Only_Results.csv", stringsAsFactors = FALSE)
grade_raw <- read.csv("Adjusted_Grade_Only_Results.csv", stringsAsFactors = FALSE)
gendergrade_raw <- read.csv("Adjusted_Gender_Grade_Results.csv", stringsAsFactors = FALSE)
full_raw <- read.csv("Adjusted_Results for gender and year and fas.csv", stringsAsFactors = FALSE)

# ── 2. Process unadjusted (reference = Gypsy/Traveller, re-ref to White) ──

unadj <- unadj_raw %>%
  filter(grepl("^eth_8cat", term)) %>%
  mutate(ethnicity = gsub("eth_8cat", "", term)) %>%
  select(ethnicity, est_raw = est, RR_raw = RR, Outcome)

# Get White's log-RR for each outcome (to re-reference)
white_ref_unadj <- unadj %>%
  filter(ethnicity == "White") %>%
  select(Outcome, est_white = est_raw)

unadj_reref <- unadj %>%
  filter(ethnicity != "White") %>%
  left_join(white_ref_unadj, by = "Outcome") %>%
  mutate(
    est_vs_white = est_raw - est_white,
    RR_unadj = exp(est_vs_white)
  ) %>%
  # Add Gypsy/Traveller (reference in original = 0, so vs White = 0 - est_white)
  bind_rows(
    white_ref_unadj %>%
      mutate(
        ethnicity = "Gypsy/Traveller",
        est_vs_white = 0 - est_white,
        RR_unadj = exp(est_vs_white)
      )
  ) %>%
  select(ethnicity, Outcome, RR_unadj)

# ── 3. Process gender-only adjusted (reference = White British) ──────────

gender_adj <- gender_raw %>%
  rename(ethnicity = term) %>%
  select(ethnicity, Outcome, RR_gender = RR)

# ── 4. Process grade-only adjusted (reference = White British) ────────────

grade_adj <- grade_raw %>%
  rename(ethnicity = term) %>%
  select(ethnicity, Outcome, RR_grade = RR)

# ── 5. Process gender+grade adjusted (wide format, parse RR from strings) ─

gendergrade <- gendergrade_raw %>%
  rename(ethnicity = term) %>%
  select(ethnicity, Outcome, RR_gendergrade = RR)

# ── 6. Process fully adjusted (reference = Gypsy/Traveller, re-ref to White)

full <- full_raw %>%
  filter(grepl("^eth_8cat", term)) %>%
  mutate(ethnicity = gsub("eth_8cat", "", term)) %>%
  select(ethnicity, est_raw = est, Outcome)

white_ref_full <- full %>%
  filter(ethnicity == "White") %>%
  select(Outcome, est_white = est_raw)

full_reref <- full %>%
  filter(ethnicity != "White") %>%
  left_join(white_ref_full, by = "Outcome") %>%
  mutate(
    est_vs_white = est_raw - est_white,
    RR_full = exp(est_vs_white)
  ) %>%
  bind_rows(
    white_ref_full %>%
      mutate(
        ethnicity = "Gypsy/Traveller",
        est_vs_white = 0 - est_white,
        RR_full = exp(est_vs_white)
      )
  ) %>%
  select(ethnicity, Outcome, RR_full)

# ── 7. Merge all ─────────────────────────────────────────────────────────

merged <- unadj_reref %>%
  left_join(gender_adj, by = c("ethnicity", "Outcome")) %>%
  left_join(grade_adj, by = c("ethnicity", "Outcome")) %>%
  left_join(gendergrade, by = c("ethnicity", "Outcome")) %>%
  left_join(full_reref, by = c("ethnicity", "Outcome"))

# ── 8. Calculate percent attenuation ─────────────────────────────────────
# Formula: ((RR_basic - 1) - (RR_adj - 1)) / (RR_basic - 1) * 100
# Simplifies to: (RR_basic - RR_adj) / (RR_basic - 1) * 100

merged <- merged %>%
  mutate(
    pct_atten_gender = ifelse(abs(RR_unadj - 1) < 0.001, NA,
      (RR_unadj - RR_gender) / (RR_unadj - 1) * 100),
    pct_atten_grade = ifelse(abs(RR_unadj - 1) < 0.001, NA,
      (RR_unadj - RR_grade) / (RR_unadj - 1) * 100),
    pct_atten_gendergrade = ifelse(abs(RR_unadj - 1) < 0.001, NA,
      (RR_unadj - RR_gendergrade) / (RR_unadj - 1) * 100),
    pct_atten_full = ifelse(abs(RR_unadj - 1) < 0.001, NA,
      (RR_unadj - RR_full) / (RR_unadj - 1) * 100)
  )

# ── 9. Create nice output table ──────────────────────────────────────────

# Outcome labels
outcome_labels <- c(
  lifesatneg = "Low life satisfaction",
  sdqinternalhigh = "High internalising",
  sdqexternalhigh = "High externalising",
  cannabis30d = "Cannabis use (30d)",
  alcoholuse = "Alcohol use",
  smokeweekly = "Weekly smoking",
  regvape = "Regular vaping",
  phys7 = "Daily physical activity",
  vigexercise = "Vigorous exercise",
  sit = "Sedentary time",
  dailysoft = "Daily soft drinks",
  dailyenergy = "Daily energy drinks",
  beenbulliedever = "Ever been bullied",
  cbeenbulliedever = "Ever cyberbullied",
  truantever = "Ever truanted",
  excludedever = "Ever excluded",
  breakfast = "Breakfast",
  lowdailyfruit = "Low daily fruit"
)

output <- merged %>%
  mutate(
    Outcome_label = ifelse(Outcome %in% names(outcome_labels),
      outcome_labels[Outcome], Outcome)
  ) %>%
  select(
    Ethnicity = ethnicity,
    Outcome = Outcome_label,
    RR_Unadjusted = RR_unadj,
    RR_Gender = RR_gender,
    RR_Grade = RR_grade,
    RR_Gender_Grade = RR_gendergrade,
    RR_Full = RR_full,
    `% Atten Gender` = pct_atten_gender,
    `% Atten Grade` = pct_atten_grade,
    `% Atten Gender+Grade` = pct_atten_gendergrade,
    `% Atten Full` = pct_atten_full
  ) %>%
  mutate(
    across(starts_with("RR"), ~ round(., 2)),
    across(starts_with("% Atten"), ~ round(., 1))
  ) %>%
  arrange(Ethnicity, Outcome)

# ── 10. Save CSV ──────────────────────────────────────────────────────────

write.csv(output, "Percent_Attenuation_Results.csv", row.names = FALSE)
cat("CSV saved.\n")

# ── 11. Create Word document ─────────────────────────────────────────────

# Focus table on Gypsy/Traveller and Roma (the key groups)
# But include all groups for completeness

# First, create a summary table for all ethnic groups
ft_data <- output %>%
  select(
    Ethnicity, Outcome,
    `RR (Unadj)`  = RR_Unadjusted,
    `RR (Gender)` = RR_Gender,
    `% Atten`     = `% Atten Gender`,
    `RR (Grade)`  = RR_Grade,
    `% Atten `    = `% Atten Grade`,
    `RR (G+Gr)`   = RR_Gender_Grade,
    `% Atten  `   = `% Atten Gender+Grade`,
    `RR (Full)`   = RR_Full,
    `% Atten   `  = `% Atten Full`
  )

# Split by ethnicity for cleaner tables
doc <- read_docx()
doc <- body_add_par(doc, "Percentage Attenuation of Relative Risks Across Adjustment Models",
  style = "heading 1")
doc <- body_add_par(doc,
  paste0("Formula: ((RR_unadjusted - 1) - (RR_adjusted - 1)) / (RR_unadjusted - 1) x 100. ",
         "Positive values indicate attenuation (reduction in excess risk); ",
         "negative values indicate amplification (increase in excess risk) after adjustment."),
  style = "Normal")
doc <- body_add_par(doc, "", style = "Normal")

for (eth in unique(ft_data$Ethnicity)) {
  eth_data <- ft_data %>% filter(Ethnicity == eth) %>% select(-Ethnicity)

  doc <- body_add_par(doc, eth, style = "heading 2")

  ft <- flextable(eth_data) %>%
    set_header_labels(
      Outcome = "Outcome",
      `RR (Unadj)` = "RR",
      `RR (Gender)` = "RR",
      `% Atten` = "% Atten",
      `RR (Grade)` = "RR",
      `% Atten ` = "% Atten",
      `RR (G+Gr)` = "RR",
      `% Atten  ` = "% Atten",
      `RR (Full)` = "RR",
      `% Atten   ` = "% Atten"
    ) %>%
    add_header_row(
      values = c("", "Unadjusted", "Adjusted: Gender", "Adjusted: Grade",
                  "Adjusted: Gender + Grade", "Adjusted: Gender + Grade + FAS"),
      colwidths = c(1, 1, 2, 2, 2, 2)
    ) %>%
    fontsize(size = 8, part = "all") %>%
    autofit() %>%
    theme_booktabs() %>%
    align(align = "center", part = "all") %>%
    align(j = 1, align = "left", part = "body") %>%
    bold(part = "header") %>%
    # Colour negative attenuation (amplification) in red
    color(j = c(3, 5, 7, 9), color = "red",
      i = ~ !is.na(`% Atten`) & `% Atten` < 0) %>%
    color(j = c(3, 5, 7, 9), color = "red",
      i = ~ !is.na(`% Atten `) & `% Atten ` < 0) %>%
    color(j = c(3, 5, 7, 9), color = "red",
      i = ~ !is.na(`% Atten  `) & `% Atten  ` < 0) %>%
    color(j = c(3, 5, 7, 9), color = "red",
      i = ~ !is.na(`% Atten   `) & `% Atten   ` < 0)

  doc <- body_add_flextable(doc, ft)
  doc <- body_add_par(doc, "", style = "Normal")
}

print(doc, target = "Percent_Attenuation_Results.docx")
cat("Word document saved.\n")
cat("\nDone!\n")

# ---------- percent_attenuation_from_docx.R ----------
# This script extracts RRs directly from the Word document
# tables and computes GT/Roma-specific attenuation.
# It produces Percent_Attenuation_GT_Roma.csv which is
# used by Figure 2 (Section 9 above).

library(officer)
library(flextable)
library(dplyr)
library(tidyr)

docx_path <- "C:/Users/wppjw/OneDrive - Cardiff University/Papers/SHRN/Gyspy Irish Roma/Paper/All_Tables 09.04.2026.docx"
out_path  <- "C:/Users/wppjw/OneDrive - Cardiff University/Papers/SHRN/Gyspy Irish Roma/Results/v0.7/Percent_Attenuation_GT_Roma.docx"

d <- read_docx(docx_path)
s <- docx_summary(d)

# Tables S1..S5 correspond to table_index 3..7 in this file
tbl_idx <- c(Unadj = 3, Gender = 4, Grade = 5, GenderGrade = 6, Full = 7)

# Columns in each table: 1=Outcome, 2=Gypsy/Traveller, 3=Roma, ...
extract_rr <- function(ti) {
  tc <- s[s$content_type == "table cell" & s$table_index == ti, ]
  w <- tc %>%
    select(row_id, cell_id, text) %>%
    pivot_wider(names_from = cell_id, values_from = text)
  names(w)[2:8] <- c("Outcome", "GT", "Roma", "WhiteOther", "Black", "Mixed", "Asian")
  # keep data rows only (those where GT has a numeric with parenthesis)
  w <- w[grepl("^[0-9]", w$GT %||% ""), ]
  parse_rr <- function(x) as.numeric(sub("^([0-9.]+).*", "\\1", x))
  data.frame(
    Outcome = trimws(w$Outcome),
    GT = parse_rr(w$GT),
    Roma = parse_rr(w$Roma),
    stringsAsFactors = FALSE
  )
}
`%||%` <- function(a, b) if (is.null(a)) b else a

rrs <- lapply(tbl_idx, extract_rr)

# Merge by Outcome
merge_one <- function(eth) {
  out <- rrs$Unadj[, c("Outcome", eth)]
  names(out)[2] <- "RR_Unadj"
  for (m in c("Gender", "Grade", "GenderGrade", "Full")) {
    tmp <- rrs[[m]][, c("Outcome", eth)]
    names(tmp)[2] <- paste0("RR_", m)
    out <- merge(out, tmp, by = "Outcome", all.x = TRUE, sort = FALSE)
  }
  out$Ethnicity <- ifelse(eth == "GT", "Gypsy/Traveller", "Roma")
  out
}

combined <- rbind(merge_one("GT"), merge_one("Roma"))

pa <- function(ref, adj) ifelse(abs(ref - 1) < 0.001, NA,
                                (ref - adj) / (ref - 1) * 100)

combined <- combined %>%
  mutate(
    `% Atten Gender`        = round(pa(RR_Unadj, RR_Gender), 1),
    `% Atten Grade`         = round(pa(RR_Unadj, RR_Grade), 1),
    `% Atten Gender+Grade`  = round(pa(RR_Unadj, RR_GenderGrade), 1),
    `% Atten Full`          = round(pa(RR_Unadj, RR_Full), 1)
  ) %>%
  mutate(across(starts_with("RR_"), ~ round(., 2)))

# Build output table: one row per outcome, columns paired RR + %Atten
ft_data <- combined %>%
  transmute(
    Ethnicity, Outcome,
    `RR (Unadj)`  = RR_Unadj,
    `RR (Gender)` = RR_Gender,
    `% Atten`     = `% Atten Gender`,
    `RR (Grade)`  = RR_Grade,
    `% Atten `    = `% Atten Grade`,
    `RR (G+Gr)`   = RR_GenderGrade,
    `% Atten  `   = `% Atten Gender+Grade`,
    `RR (Full)`   = RR_Full,
    `% Atten   `  = `% Atten Full`
  )

doc <- read_docx()
doc <- body_add_par(doc,
  "Percentage Attenuation of Relative Risks Across Adjustment Models (Gypsy/Traveller and Roma)",
  style = "heading 1")
doc <- body_add_par(doc,
  paste0("Source: RR point estimates extracted from All_Tables 09.04.2026.docx (Tables S1-S5). ",
         "Formula: ((RR_unadjusted - 1) - (RR_adjusted - 1)) / (RR_unadjusted - 1) x 100. ",
         "Positive values indicate attenuation (reduction in excess risk); ",
         "negative values indicate amplification after adjustment."),
  style = "Normal")
doc <- body_add_par(doc, "", style = "Normal")

for (eth in c("Gypsy/Traveller", "Roma")) {
  ed <- ft_data %>% filter(Ethnicity == eth) %>% select(-Ethnicity)
  doc <- body_add_par(doc, eth, style = "heading 2")
  ft <- flextable(ed) %>%
    set_header_labels(
      Outcome = "Outcome",
      `RR (Unadj)` = "RR", `RR (Gender)` = "RR", `% Atten` = "% Atten",
      `RR (Grade)` = "RR", `% Atten ` = "% Atten",
      `RR (G+Gr)` = "RR", `% Atten  ` = "% Atten",
      `RR (Full)` = "RR", `% Atten   ` = "% Atten"
    ) %>%
    add_header_row(
      values = c("", "Unadjusted", "Adjusted: Gender", "Adjusted: Grade",
                 "Adjusted: Gender + Grade", "Adjusted: Gender + Grade + FAS"),
      colwidths = c(1, 1, 2, 2, 2, 2)
    ) %>%
    fontsize(size = 8, part = "all") %>%
    autofit() %>%
    theme_booktabs() %>%
    align(align = "center", part = "all") %>%
    align(j = 1, align = "left", part = "body") %>%
    bold(part = "header") %>%
    color(j = c(3, 5, 7, 9), color = "red", i = ~ !is.na(`% Atten`)   & `% Atten`   < 0) %>%
    color(j = c(3, 5, 7, 9), color = "red", i = ~ !is.na(`% Atten `)  & `% Atten `  < 0) %>%
    color(j = c(3, 5, 7, 9), color = "red", i = ~ !is.na(`% Atten  `) & `% Atten  ` < 0) %>%
    color(j = c(3, 5, 7, 9), color = "red", i = ~ !is.na(`% Atten   `)& `% Atten   `< 0)
  doc <- body_add_flextable(doc, ft)
  doc <- body_add_par(doc, "", style = "Normal")
}

print(doc, target = out_path)
write.csv(combined, sub("\\.docx$", ".csv", out_path), row.names = FALSE)
cat("Saved:", out_path, "\n")


# ************************************************************
# ************************************************************
#
# ORCHESTRATOR SCRIPT (for reference)
# Source: grt_subgroup_v0.91.R
#
# This script coordinates the subgroup analysis by:
# 1. Filtering the imputed data to the subgroup
#    (syear==2025 + syear==2023 & grade 10/11)
# 2. Running Sophie_descriptives_v0.14.R (T1, T2, S1)
# 3. Running Sophie_adjusted_separate_v0.8.R (S2, S3)
# 4. Running Sophie_adjusted_gender_grade_v0.7.R (S4)
# 5. Running grt_estimation_2023_2025_v0.2.R (S5 + CSV)
#
# It rewrites paths/captions and stubs CSV writes so all
# outputs land in a single output directory.
#
# ************************************************************
# ************************************************************

cat("\n[ORCHESTRATOR: See grt_subgroup_v0.91.R for full code]\n")

# ============================================================
