# GRT script: 44 imputations - PARALLEL (Stata .dta input)
# Dataset: GRT Dataset SHW 2023-2025 v0.2.dta
# ====================================================================================
# Uses 'futuremice' for scientifically correct, reproducible parallel imputation.
# Reads Stata .dta file via haven package.
# ====================================================================================

# ── 1. Setup & Packages ──
user_lib <- file.path(Sys.getenv("HOME"), "R_libs")
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)
.libPaths(c(user_lib, .libPaths()))

pkgs <- c("mice", "miceadds", "estimatr", "future", "furrr", "haven")
new_pkgs <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(new_pkgs) > 0) {
  install.packages(new_pkgs, lib = user_lib, repos = "https://cloud.r-project.org")
}

library(mice)
library(miceadds)
library(estimatr)
library(future)
library(haven)

# ── 2. Resource Management ──
Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1")
Sys.setenv(OPENBLAS_NUM_THREADS = "1")

N_CORES <- 44L
plan(multisession, workers = N_CORES)

# ── 3. Data Loading & Prep ──
BASEDIR <- "/shared/home1/c.wppjw/SHRN/Gypsy RT"

INFILE <- file.path(BASEDIR, "GRT Dataset SHW 2023-2025 v0.2.dta")

cat("Loading data...\n")
shrndata <- read_dta(INFILE)
shrndata <- as.data.frame(shrndata)

# Define types
factor_vars <- c("gender", "grade", "eth_8cat", "syear",
                 "regvape", "lifesatneg", "sdqinternalhigh", "sdqexternalhigh",
                 "cannabis30d", "alcoholuse", "smokeweekly",
                 "phys7", "vigexercise", "dailysoft", "dailyenergy", "sit",
                 "beenbulliedever", "cbeenbulliedever",
                 "truantever", "excludedever")

numeric_vars <- c("fas", "peer", "teacher", "famsupp")

for (v in factor_vars) shrndata[[v]] <- as.factor(shrndata[[v]])
for (v in numeric_vars) shrndata[[v]] <- as.numeric(shrndata[[v]])

# ── 4. Setup Predictor Matrix ──
ini  <- mice(shrndata, maxit = 0, printFlag = FALSE)
meth <- ini$method
pred <- ini$predictorMatrix

# id4: participant ID - exclude from imputation and prediction entirely
meth["id4"] <- ""
pred["id4", ] <- 0
pred[, "id4"] <- 0

# id2c: clustering variable - do not impute, but use as predictor
meth["id2c"] <- ""
pred["id2c", ] <- 0
# id2c remains as a column predictor (pred[, "id2c"] stays as-is)

# syear: auxiliary variable with complete data - do not impute, use as predictor
meth["syear"] <- ""
pred["syear", ] <- 0
# syear remains as a column predictor (pred[, "syear"] stays as-is)

# ── 5. Run Parallel Imputation ──
cat(sprintf("Starting parallel imputation (m=44) on %d cores...\n", N_CORES))
t_start <- proc.time()

imp <- futuremice(
  shrndata,
  m               = 44,
  maxit           = 25,
  method          = meth,
  predictorMatrix = pred,
  nnet.MaxNWts    = 5000,
  parallelseed    = 54321,
  n.core          = N_CORES
)

t_total <- (proc.time() - t_start)[3]
cat(sprintf("Done! Total time: %.2f minutes\n", t_total / 60))

# ── 6. Save Results ──
OUTDIR <- "/shared/home1/c.wppjw/SHRN/Gypsy RT/output/v0.4"
if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

saveRDS(imp, file.path(OUTDIR, "GRT_2023_2025_imp_m44_parallel.rds"))

imputed <- complete(imp, action = "long", include = TRUE)
write.csv(imputed, file.path(OUTDIR, "GRT_2023_2025_m44_imputed_dataset.csv"), row.names = FALSE)

# ── 7. Diagnostics ──
pdf(file.path(OUTDIR, "GRT_2023_2025_convergence_m44.pdf"))
plot(imp)
dev.off()

cat("Script completed successfully.\n")
