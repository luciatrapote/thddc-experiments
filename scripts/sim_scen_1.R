################################################################################
# SIMULATION EXPERIMENT: tHDDC vs tclust vs rlg
# Scenario 1: High-Dimensional Clusters (p=200) with Contamination
            # Intrinsic dimensions for clusters: 3 and 1
            # Each group has an othogonal matrix
################################################################################

library(clue)
library(tclust)
library(ggplot2)
library(patchwork)
library(doParallel)
library(doRNG)

# --- Global Settings ---
set.seed(123)
N <- 20                # Internal iterations
D <- 7                 # Number of delta values (separation)
p <- 200               # Dimensionality
delta_vals <- seq(-0.3, 0.3, length.out = D)
alpha_trim <- 0.05     # Trimming level

# Accuracy function (Diagonal agreement)
accuracy_clustering <- function(x, y) {
  part_x <- as.cl_partition(x)
  part_y <- as.cl_partition(y)
  cl_agreement(part_x, part_y, method = "diag")
}

# --- Data Structures for Results ---

# Initialize matrices for each method
accuracy_catell <- matrix(0, nrow = D, ncol = N)
accuracy_catell1  <- matrix(0, nrow = D, ncol = N)
accuracy_catell2  <- matrix(0, nrow = D, ncol = N)
accuracy_catell3  <- matrix(0, nrow = D, ncol = N)
accuracy_tclust  <- matrix(0, nrow = D, ncol = N)
accuracy_tclust1  <- matrix(0, nrow = D, ncol = N)
accuracy_tclust2  <- matrix(0, nrow = D, ncol = N)
accuracy_rlg     <- matrix(0, nrow = D, ncol = N)


time_catell <- matrix(0, nrow = D, ncol = N)
time_catell1  <- matrix(0, nrow = D, ncol = N)
time_catell2  <- matrix(0, nrow = D, ncol = N)
time_catell3  <- matrix(0, nrow = D, ncol = N)
time_tclust  <- matrix(0, nrow = D, ncol = N)
time_tclust1  <- matrix(0, nrow = D, ncol = N)
time_tclust2  <- matrix(0, nrow = D, ncol = N)
time_rlg     <- matrix(0, nrow = D, ncol = N)

#Matrices para almacenar las dimensiones estimadas
q_catell <- array(0, dim = c(D, N, 2))
q_catell1  <- array(0, dim = c(D, N, 2))
q_catell2  <- array(0, dim = c(D, N, 2))
q_catell3  <- array(0, dim = c(D, N, 2))

# --- Simulation Parameters ---
n1 <- 570
n2 <- 380
n_cont <- 50

# First group
q1=3
mu1 <- rep(0, p)
Delta1 <- diag(c(seq(9,7,length=q1), rep(0.15, p - q1)))

# Second group
q2=1
Delta2 <- diag(c(rep(5,q2), rep(0.45, p - q2)))

# Vector of intrinsic true dimensions
qg <- c(q1, q2)


# --- MAIN SIMULATION LOOP ---
for (d in 1:D)
{
  delta <- delta_vals[d]
  
  for (i in 1:N)
  {
    # 1. Data Generation
    # Group 1
    A1 <- matrix(rnorm((p + 1) * p), ncol = p)
    B1 <- cov(A1)
    U1 <- eigen(B1)$vector
    X1 <- rep(1, n1) %*% t(mu1) +
      matrix(rnorm(n1 * p), nrow = n1) %*% diag(sqrt(diag(Delta1))) %*% U1
    # Group 2
    mu2 <- rep(delta, p)
    A2 <- matrix(rnorm((p + 1) * p), ncol = p)
    B2 <- cov(A2)
    U2 <- eigen(B2)$vector
    X2 <- rep(1, n2) %*% t(mu2) +
      matrix(rnorm(n2 * p), nrow = n2) %*% diag(sqrt(diag(Delta2))) %*% U2
    # Contamination
    X.cont <- matrix(runif(n_cont*p,-2,2),ncol=p)
    
    # Final data set
    X <- rbind(X1, X2, X.cont)
    true_labels <- c(rep(1,n1), rep(2,n2), rep(3,50))
    
    # ---------------------------------------------------------
    # 1) tHDDC CON CATELL (3,1)
    # ---------------------------------------------------------
    ini <- Sys.time()
    uX <- tHDDC_paralelizado_foreach(
      X, q=qg, alpha=0.05, nstart = 250, cstep1=2, cstep2=25, nkeep=5,
      equal.weights=TRUE, restr.fact=c(5,3), zero.tol=1e-16,
      trace=0, opt="HARD", sol_ini_p=FALSE, sol_ini=NA,
      catell_Y=TRUE, catell_threshold=0.3, dim_max=20
    )
    fin <- Sys.time()
    
    time_catell[d,i]     <- as.numeric(fin-ini, units="secs")
    accuracy_catell[d,i] <- accuracy_clustering(uX$assig, true_labels)
    q_catell[d, i, ] <- uX$q
    
    # ---------------------------------------------------------
    # 2) tHDDC CON CATELL (1,1)
    # ---------------------------------------------------------
    ini <- Sys.time()
    uX <- tHDDC_paralelizado_foreach(
      X, q=c(1,1), alpha=0.05, nstart = 250, cstep1=2, cstep2=25, nkeep=5,
      equal.weights=TRUE, restr.fact=c(5,3), zero.tol=1e-16,
      trace=0, opt="HARD", sol_ini_p=FALSE, sol_ini=NA,
      catell_Y=TRUE, catell_threshold=0.3, dim_max=20
    )
    fin <- Sys.time()
    
    time_catell1[d,i]     <- as.numeric(fin-ini, units="secs")
    accuracy_catell1[d,i] <- accuracy_clustering(uX$assig, true_labels)
    q_catell1[d, i, ] <- uX$q
    
    # ---------------------------------------------------------
    # 3) tHDDC CON CATELL (3,3)
    # ---------------------------------------------------------
    ini <- Sys.time()
    uX <- tHDDC_paralelizado_foreach(
      X, q=c(3,3), alpha=0.05, nstart = 250, cstep1=2, cstep2=25, nkeep=5,
      equal.weights=TRUE, restr.fact=c(5,3), zero.tol=1e-16,
      trace=0, opt="HARD", sol_ini_p=FALSE, sol_ini=NA,
      catell_Y=TRUE, catell_threshold=0.3, dim_max=20
    )
    fin <- Sys.time()
    
    time_catell2[d,i]     <- as.numeric(fin-ini, units="secs")
    accuracy_catell2[d,i] <- accuracy_clustering(uX$assig, true_labels)
    q_catell2[d, i, ] <- uX$q
    
    # ---------------------------------------------------------
    # 4) tHDDC CON CATELL (5,5)
    # ---------------------------------------------------------
    ini <- Sys.time()
    uX <- tHDDC_paralelizado_foreach(
      X, q=c(5,5), alpha=0.05, nstart = 250, cstep1=2, cstep2=25, nkeep=5,
      equal.weights=TRUE, restr.fact=c(5,3), zero.tol=1e-16,
      trace=0, opt="HARD", sol_ini_p=FALSE, sol_ini=NA,
      catell_Y=TRUE, catell_threshold=0.3, dim_max=20
    )
    fin <- Sys.time()
    
    time_catell3[d,i]     <- as.numeric(fin-ini, units="secs")
    accuracy_catell3[d,i] <- accuracy_clustering(uX$assig, true_labels)
    q_catell3[d, i, ] <- uX$q
    
    # ---------------------------------------------------------
    # 5) TCLUST (2 clusters; accurate trimming percentage; c=12 (default) used to ensure flexibility)
    # ---------------------------------------------------------
    ini <- Sys.time()
    tc <- tclust(X,K=2,alpha = 0.05, nstart = 250, cstep1=2, cstep2=25, nkeep=5,  equal.weights = FALSE, 
                 restr.fact=12,zero.tol = 1e-16,  trace = 0,   opt="HARD",
                 sol_ini_p = FALSE, sol_ini = NA ) 
    fin <- Sys.time()
    
    time_tclust[d,i]     <- as.numeric(fin-ini, units="secs")
    accuracy_tclust[d,i] <- accuracy_clustering(tc$assig, true_labels)
    
    # ---------------------------------------------------------
    # 5) TCLUST (2 clusters; accurate trimming percentage; c=50 used to ensure flexibility)
    # ---------------------------------------------------------
    ini <- Sys.time()
    tc <- tclust(X,K=2,alpha = 0.05, nstart = 250, cstep1=2, cstep2=25, nkeep=5,  equal.weights = FALSE, 
                 restr.fact=50,zero.tol = 1e-16,  trace = 0,   opt="HARD",
                 sol_ini_p = FALSE, sol_ini = NA ) 
    fin <- Sys.time()
    
    time_tclust1[d,i]     <- as.numeric(fin-ini, units="secs")
    accuracy_tclust1[d,i] <- accuracy_clustering(tc$assig, true_labels)
    
    # ---------------------------------------------------------
    # 5) TCLUST (2 clusters; accurate trimming percentage; c=150 used to ensure flexibility)
    # ---------------------------------------------------------
    ini <- Sys.time()
    tc <- tclust(X,K=2,alpha = 0.05, nstart = 250, cstep1=2, cstep2=25, nkeep=5,  equal.weights = FALSE, 
                 restr.fact=150,zero.tol = 1e-16,  trace = 0,   opt="HARD",
                 sol_ini_p = FALSE, sol_ini = NA ) 
    fin <- Sys.time()
    
    time_tclust2[d,i]     <- as.numeric(fin-ini, units="secs")
    accuracy_tclust2[d,i] <- accuracy_clustering(tc$assig, true_labels)
    
    
    # ---------------------------------------------------------
    # 7) RLG (2 grupos, buenas dimensiones)
    # ---------------------------------------------------------
    Sys.setenv(OMP_NUM_THREADS="4")
    ini <- Sys.time()
    rg2 <- rlg(X, d=c(q1,q2), alpha=0.05, nstart=250, niter1=2, niter2=25, nkeep=5, parallel = TRUE, n.cores=20)
    fin <- Sys.time()
    
    time_rlg[d,i]     <- as.numeric(fin-ini, units="secs")
    accuracy_rlg[d,i] <- accuracy_clustering(rg2$cluster, true_labels)
    
    stopImplicitCluster()
  }
}

# --- VISUALIZATION HELPERS ---

load("~/Simulaciones: Extending tclust to higher dimensions/sim_scen_1.RData")

# -----------------------------------------------------------
# FIGURE 3 (a)
# -----------------------------------------------------------

# Boxplot of the accuracy obtained by TCLUST, RLG, and tHDDC (RLG and tHDDC 
#initialized with the true intrinsic dimensions) for different values of delta.

metodos <- c("tHDDC", "tclust", "rlg")
accuracy_list <- list(accuracy_catell, accuracy_tclust, accuracy_rlg)

df_accuracy <- do.call(rbind,
                       lapply(1:3, function(m) {
                         data.frame(
                           delta = rep(delta_vals, each=N),
                           metodo = metodos[m],
                           accuracy = as.vector(t(accuracy_list[[m]]))
                         )
                       })
)

df_accuracy$metodo <- factor(
  df_accuracy$metodo,
  levels = metodos
)

library(ggplot2)
library(patchwork)

p1 <- ggplot(
  df_accuracy,
  aes(
    x = factor(delta),
    y = accuracy,
    fill = metodo,
    colour = metodo
  )
) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(
    title = "",
    x = "Delta", 
    y = "Accuracy",
    fill = "Method",
    colour = "Method"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

p1 <- p1 +
  coord_cartesian(ylim = c(0.47, 1))

p1


# -----------------------------------------------------------
# FIGURE 4 (a)
# -----------------------------------------------------------

# Boxplot of the classification accuracy obtained by the tHDDC method when 
#using q1 and q2 equal to true values of the intrinsic dimensions, and also when
#using different values  q_{ini} for the  automatized determination of intrinsic
#dimensions based on the Cattell's approach, for each value of delta.

metodos <- c("tHDDC", "tHDDC q_{ini}=1", "tHDDC q_{ini}=3", "tHDDC q_{ini}=5")
accuracy_list <- list(accuracy_catell, accuracy_catell1, accuracy_catell2, accuracy_catell3)

labels_metodos <- c(
  expression(tHDDC~paste("(", 3, ",", 1, ")")),
  expression(tHDDC~q[ini]==1),
  expression(tHDDC~q[ini]==3),
  expression(tHDDC~q[ini]==5)
)


df_accuracy <- do.call(rbind,
                       lapply(1:4, function(m) {
                         data.frame(
                           delta = rep(delta_vals, each=N),
                           metodo = metodos[m],
                           accuracy = as.vector(t(accuracy_list[[m]]))
                         )
                       })
)

df_accuracy$metodo <- factor(
  df_accuracy$metodo,
  levels = metodos
)

library(ggplot2)
library(patchwork)

p1 <- ggplot(
  df_accuracy,
  aes(
    x = factor(delta),
    y = accuracy,
    fill = metodo,
    colour = metodo
  )
) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(
    title = "",
    x = "Delta", 
    y = "Accuracy",
    fill = "Method",
    colour = "Method"
  ) +
  scale_fill_discrete(labels = labels_metodos) +
  scale_colour_discrete(labels = labels_metodos) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

p1 <- p1 +
  coord_cartesian(ylim = c(0.5, 1))
p1


# -----------------------------------------------------------
# TABLE 1 (First row)
# -----------------------------------------------------------

#Ratio between the computation times of TCLUST and tHDDC for each value of the 
#separation parameter delta.

ratio <- rowMeans(time_tclust)/rowMeans(time_catell)

tabla_resumen <- data.frame(
  delta = delta_vals,
  ratio_time= ratio)

tabla_resumen

# -----------------------------------------------------------
# FIGURE 5 (a)
# -----------------------------------------------------------

#Boxplot of the computation time for tHDDC different initial q_{ini} values for
#automated determination of the intrinsic dimensions, for each value of delta.

metodos <- c("tHDDC", "tHDDC q_{ini}=1", "tHDDC q_{ini}=3", "tHDDC q_{ini}=5")
accuracy_list <- list(time_catell, time_catell1, time_catell2, time_catell3)

labels_metodos <- c(
  expression(tHDDC~paste("(", 3, ",", 1, ")")),
  expression(tHDDC~q[ini]==1),
  expression(tHDDC~q[ini]==3),
  expression(tHDDC~q[ini]==5)
)

df_time <- do.call(rbind,
                   lapply(1:4, function(m) {
                     data.frame(
                       delta = rep(delta_vals, each=N),
                       metodo = metodos[m],
                       accuracy = as.vector(t(accuracy_list[[m]]))
                     )
                   })
)

df_time$metodo <- factor(
  df_time$metodo,
  levels = metodos
)

library(ggplot2)
library(patchwork)

p1 <- ggplot(df_time, aes(x = factor(delta), y = accuracy, fill = metodo, colour=metodo)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(
    title = "",
    x = "Delta", 
    y = "Time",
    fill = "Method",
    colour = "Method"
  ) +
  scale_fill_discrete(labels = labels_metodos) +
  scale_colour_discrete(labels = labels_metodos) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

p1 <- p1 +
  coord_cartesian(ylim = c(6.5, 9.5))

p1


# -----------------------------------------------------------
# FIGURE 6 (a)
# -----------------------------------------------------------

# Boxplots of the estimated larger intrinsic dimension obtained by the tHDDC 
#method using the Cattell-based procedure for different values of the separation
#parameter delta.

k1=max(qg)

err_catell <- matrix(0, nrow = D, ncol = N)
err_catell1 <- matrix(0, nrow = D, ncol = N)
err_catell2 <- matrix(0, nrow = D, ncol = N)
err_catell3 <- matrix(0, nrow = D, ncol = N)


err_catell <- apply(q_catell, c(1, 2), max)-k1
err_catell1 <- apply(q_catell1, c(1, 2), max)-k1
err_catell2 <- apply(q_catell2, c(1, 2), max)-k1
err_catell3 <- apply(q_catell3, c(1, 2), max)-k1

err_list <- list(
  Catell  = err_catell,
  Catell1 = err_catell1,
  Catell2 = err_catell2,
  Catell3 = err_catell3
)

df <- data.frame()

for (d in 1:D) {
  for (name in names(err_list)) {
    df <- rbind(
      df,
      data.frame(
        delta  = delta_vals[d],
        metodo = name,
        error  = err_list[[name]][d, ]
      )
    )
  }
}


ggplot(df, aes(x = factor(delta), y = error,
               fill = metodo, colour = metodo)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5)+
  geom_boxplot(position = position_dodge(width = 0.8)) +
  scale_fill_discrete(
    name = "Method",
    labels = c(
      expression(tHDDC~~"(" * 3 * "," * 1 * ")"),
      expression(tHDDC~~q[ini]==1),
      expression(tHDDC~~q[ini]==3),
      expression(tHDDC~~q[ini]==5)
    )
  ) +
  scale_colour_discrete(
    name = "Method",
    labels = c(
      expression(tHDDC~~"(" * 3 * "," * 1 * ")"),
      expression(tHDDC~~q[ini]==1),
      expression(tHDDC~~q[ini]==3),
      expression(tHDDC~~q[ini]==5)
    )
  ) +
  labs(
    title = "",
    x = "Delta",
    y = "Dimension"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )+
  scale_y_continuous(
    breaks = seq(-2, 10, by = 2),
    labels = seq(1, 13, by = 2)
  ) +
  coord_cartesian(ylim = c(-2, 2))


# -----------------------------------------------------------
# FIGURE 7 (a)
# -----------------------------------------------------------

# Boxplots of the estimated smaller intrinsic dimension obtained by the tHDDC 
#method using the Cattell-based procedure for different values of the separation
#parameter delta.

k2=min(qg)

err_catell <- matrix(0, nrow = D, ncol = N)
err_catell1 <- matrix(0, nrow = D, ncol = N)
err_catell2 <- matrix(0, nrow = D, ncol = N)
err_catell3 <- matrix(0, nrow = D, ncol = N)


err_catell <- apply(q_catell, c(1, 2), min)-k2
err_catell1 <- apply(q_catell1, c(1, 2), min)-k2
err_catell2 <- apply(q_catell2, c(1, 2), min)-k2
err_catell3 <- apply(q_catell3, c(1, 2), min)-k2

err_list <- list(
  Catell  = err_catell,
  Catell1 = err_catell1,
  Catell2 = err_catell2,
  Catell3 = err_catell3
)

df <- data.frame()

for (d in 1:D) {
  for (name in names(err_list)) {
    df <- rbind(
      df,
      data.frame(
        delta  = delta_vals[d],
        metodo = name,
        error  = err_list[[name]][d, ]
      )
    )
  }
}


ggplot(df, aes(x = factor(delta), y = error,
               fill = metodo, colour = metodo)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5)+
  geom_boxplot(position = position_dodge(width = 0.8)) +
  scale_fill_discrete(
    name = "Method",
    labels = c(
      expression(tHDDC~~"(" * 3 * "," * 1 * ")"),
      expression(tHDDC~~q[ini]==1),
      expression(tHDDC~~q[ini]==2),
      expression(tHDDC~~q[ini]==3)
    )
  ) +
  scale_colour_discrete(
    name = "Method",
    labels = c(
      expression(tHDDC~~"(" * 3 * "," * 1 * ")"),
      expression(tHDDC~~q[ini]==1),
      expression(tHDDC~~q[ini]==2),
      expression(tHDDC~~q[ini]==3)
    )
  ) +
  labs(
    title = "",
    x = "Delta",
    y = "Dimension"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )+
  scale_y_continuous(
    breaks = seq(-2, 10, by = 2),
    labels = seq(-1, 11, by = 2)
  ) +
  coord_cartesian(ylim = c(-2, 2))



