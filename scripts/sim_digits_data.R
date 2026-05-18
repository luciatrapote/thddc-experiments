################################################################################
#######################  ALGORITHM SIMULATIONS AND BENCHMARKING ################
################################################################################

# -------------------------------------------------------------------------
# 1. DATA LOADING AND PREPARATION
# -------------------------------------------------------------------------

# Load the dataset from the CSV file
# ZZZ_df <- read.csv("data_ZZZ_df.csv")
ZZZ_df <- read.csv("~/Simulaciones: Extending tclust to higher dimensions/Simulaciones: Digits Data/data_ZZZ_df.csv")

# Extract numeric features for clustering (columns 1 to 256)
# High-dimensional input matrix
X_matrix <- as.matrix(ZZZ_df[, 1:256])


# -------------------------------------------------------------------------
# 2. PROPOSED MODEL EXECUTION: tHDDC
# -------------------------------------------------------------------------

set.seed(113) # For reproducibility

start_tHDDC <- Sys.time()

# Run the parallelized tHDDC algorithm
# Parameters: q (intrinsic dimensions), alpha (trimming percentage)
uX <- tHDDC_paralelizado_foreach(
  X_matrix, q = c(1, 1, 1), alpha = 0.20, 
  nstart = 200, cstep1 = 10, cstep2 = 150, nkeep = 5,
  equal.weights = TRUE, restr.fact = c(5, 1.1), zero.tol = 1e-16,
  trace = 0, opt = "HARD", sol_ini_p = FALSE, sol_ini = NA,
  catell_Y = TRUE, catell_threshold = 0.2, dim_max = 20)

end_tHDDC <- Sys.time()

#saveRDS(uX, file = "tHDDC_results_uX.rds")

cat("tHDDC Execution Time:", end_tHDDC - start_tHDDC, "\n")

# --- Ground Truth Label Preparation ---
# We assess if outliers are correctly trimmed and if digits 3, 5, and 8 
# are properly identified. 
# Labels 0, 2, 4, and 9 are merged into 'out1' (natural noise/other digits).
# Label 999 is renamed to 'out2' (artificially generated outlier patterns).

ZZZ_df$cls_fusionada <- ZZZ_df$cls
ZZZ_df$cls_fusionada[ZZZ_df$cls_fusionada %in% c(0, 2, 4, 9)] <- "out1"
ZZZ_df$cls_fusionada[ZZZ_df$cls_fusionada == 999] <- "out2"

# Confusion Matrix: tHDDC Assignments vs. Merged Ground Truth
cat("Confusion Matrix for tHDDC:\n")
table(uX$assig, ZZZ_df$cls_fusionada)


# -------------------------------------------------------------------------
# 3. BENCHMARK: HDDC (High-Dimensional Data Clustering)
# -------------------------------------------------------------------------

library(HDclassif)

set.seed(113) # For reproducibility

start_HDDC <- Sys.time()

res_hddc <- hddc(
  X_matrix,
  K = 1:3,
  model = "AkjBkQkDk",
  threshold = 0.2,
  criterion = "bic",
  itermax = 150,
  eps = 1e-16,
  algo = "EM",
  d_select = "Cattell",
  init = "mini-em",
  show = getHDclassif.show(),
  mini.nb = c(200, 10),
  scaling = FALSE,
  min.individuals = 2,
  noise.ctrl = 0,
  mc.cores = 80,
  nb.rep = 1,
  keepAllRes = TRUE,
  kmeans.control = list(),
  d_max = 20,
  subset = Inf
)

end_HDDC <- Sys.time()
cat("HDDC Execution Time:", end_HDDC - start_HDDC, "\n")

# Confusion Matrix for HDDC
cat("Confusion Matrix for HDDC:\n")
table(res_hddc$class, ZZZ_df$cls_fusionada)


# -------------------------------------------------------------------------
# 4. BENCHMARK: TCLUST (Robust Clustering)
# -------------------------------------------------------------------------

set.seed(113) # For reproducibility

start_TCLUST <- Sys.time()

uX_tclust <- tclust(
  X_matrix, K = 3, alpha = 0.2, 
  nstart = 200, cstep1 = 10, cstep2 = 150, nkeep = 5, 
  equal.weights = TRUE, restr.fact = 12, zero.tol = 1e-6, 
  trace = 0, opt = "HARD", sol_ini_p = FALSE, sol_ini = NA
)

end_TCLUST <- Sys.time()
cat("TCLUST Execution Time:", end_TCLUST - start_TCLUST, "\n")

# Confusion Matrix for TCLUST
cat("Confusion Matrix for TCLUST:\n")
table(uX_tclust$assig, ZZZ_df$cls_fusionada)