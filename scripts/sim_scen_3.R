################################################################################
# SIMULATION EXPERIMENT: tHDDC vs tclust vs rlg
# Scenario 3: High-Dimensional Clusters (p=200) with Contamination
# Intrinsic dimensions for clusters: 10 and 5
# Both groups share the same othogonal matrix
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
acc1<-matrix(0, nrow = N, ncol = D)
acc2<-matrix(0, nrow = N, ncol = D)
acc3<-matrix(0, nrow = N, ncol = D)
acc4<-matrix(0, nrow = N, ncol = D)
acc5<-matrix(0, nrow = N, ncol = D)
acc6<-matrix(0, nrow = N, ncol = D)
acc7<-matrix(0, nrow = N, ncol = D)
acc8<-matrix(0, nrow = N, ncol = D)

time1<-matrix(0, nrow = N, ncol = D)
time2<-matrix(0, nrow = N, ncol = D)
time3<-matrix(0, nrow = N, ncol = D)
time4<-matrix(0, nrow = N, ncol = D)
time5<-matrix(0, nrow = N, ncol = D)
time6<-matrix(0, nrow = N, ncol = D)
time7<-matrix(0, nrow = N, ncol = D)
time8<-matrix(0, nrow = N, ncol = D)

q_catell1 <- array(0, dim = c(D, N, 2))
q_catell2  <- array(0, dim = c(D, N, 2))
q_catell3  <- array(0, dim = c(D, N, 2))
q_catell4  <- array(0, dim = c(D, N, 2))


# --- Simulation Parameters ---

# Observations in each cluster
n1 <- 570
n2 <- 380
n_cont <- 50

true_labels <- c(rep(1,n1), rep(2,n2), rep(3,n_cont))

# Dimensions of each cluster
q1 <- 10
q2 <- 5

# Vector of intrinsic true dimensions
qg <- c(q1, q2)

# Means of each group
mu1 <- rep(0,p)
# mean of second group changes with delta

# Covariance matrices of each group
Delta1 <- diag(c(seq(20,15,length=q1), rep(0.8, p - q1)))
Delta2 <- diag(c(seq(18,16,length=q2), rep(0.4, p - q2)))



# --- MAIN SIMULATION LOOP ---

for (d in 1:D) {
  for (n in 1:N) {
    # Data Generation

    # We generate the normal groups with the same rotation
    A <- matrix(rnorm((p + 1) * p), ncol = p)
    B <- cov(A)
    Ushared <- eigen(B)$vector
    
    X1 <- rep(1, n1) %*% t(mu1) +
      matrix(rnorm(n1 * p), nrow = n1) %*% diag(sqrt(diag(Delta1))) %*% Ushared
    
    mu2<-rep(delta[d],p)
    X2 <- rep(1, n2) %*% t(mu2) +
      matrix(rnorm(n2 * p), nrow = n2) %*% diag(sqrt(diag(Delta2))) %*% Ushared
    
    # Contamination
    X.cont <- matrix(runif(n_cont*p, -2, 2), ncol = p)
    
    X <- rbind(X1, X2, X.cont)
    
    # ---------------------------------------------------------
    # 1) tHDDC CON CATELL (10,5)
    # ---------------------------------------------------------
    ini <- Sys.time()
    uX <- tHDDC_paralelizado_foreach(
      X, q=qg, alpha=0.05, nstart = 250, cstep1=2, cstep2=25, nkeep=5,
      equal.weights=TRUE, restr.fact=c(5,3), zero.tol=1e-16,
      trace=0, opt="HARD", sol_ini_p=FALSE, sol_ini=NA,
      catell_Y=TRUE, catell_threshold=0.3, dim_max=40
    )
    fin <- Sys.time()
    time1[n,d]<-as.numeric(fin-ini, units="secs")
    acc1[n,d]<-accuracy_clustering(uX$assig, true_labels)
    q_catell1[d, n, ] <- uX$q
    
    # ---------------------------------------------------------
    # 2) tHDDC CON CATELL (1,1)
    # ---------------------------------------------------------
    ini <- Sys.time()
    uX <- tHDDC_paralelizado_foreach(
      X, q=c(1,1), alpha=0.05, nstart = 250, cstep1=2, cstep2=25, nkeep=5,
      equal.weights=TRUE, restr.fact=c(5,3), zero.tol=1e-16,
      trace=0, opt="HARD", sol_ini_p=FALSE, sol_ini=NA,
      catell_Y=TRUE, catell_threshold=0.3, dim_max=40
    )
    fin <- Sys.time()
    time2[n,d]<-as.numeric(fin-ini, units="secs")
    acc2[n,d]<-accuracy_clustering(uX$assig, true_labels)
    q_catell2[d, n, ] <- uX$q
    
    # ---------------------------------------------------------
    # 3) tHDDC CON CATELL (5,5)
    # ---------------------------------------------------------
    ini <- Sys.time()
    uX <- tHDDC_paralelizado_foreach(
      X, q=c(5,5), alpha=0.05, nstart = 250, cstep1=2, cstep2=25, nkeep=5,
      equal.weights=TRUE, restr.fact=c(5,3), zero.tol=1e-16,
      trace=0, opt="HARD", sol_ini_p=FALSE, sol_ini=NA,
      catell_Y=TRUE, catell_threshold=0.3, dim_max=40
    )
    fin <- Sys.time()
    time3[n,d]<-as.numeric(fin-ini, units="secs")
    acc3[n,d]<-accuracy_clustering(uX$assig, true_labels)
    q_catell3[d, n, ] <- uX$q
    
    # ---------------------------------------------------------
    # 4) tHDDC CON CATELL (15,15)
    # ---------------------------------------------------------
    ini <- Sys.time()
    uX <- tHDDC_paralelizado_foreach(
      X, q=c(15,15), alpha=0.05, nstart = 250, cstep1=2, cstep2=25, nkeep=5,
      equal.weights=TRUE, restr.fact=c(5,3), zero.tol=1e-16,
      trace=0, opt="HARD", sol_ini_p=FALSE, sol_ini=NA,
      catell_Y=TRUE, catell_threshold=0.3, dim_max=40
    )
    fin <- Sys.time()
    time4[n,d]<-as.numeric(fin-ini, units="secs")
    acc4[n,d]<-accuracy_clustering(uX$assig, true_labels)
    q_catell4[d, n, ] <- uX$q
    
    # ---------------------------------------------------------
    # 4) TCLUST (2 clusters; accurate trimming percentage; c=12 (default) used to ensure flexibility)
    # ---------------------------------------------------------
    ini <- Sys.time()
    tc <- tclust(X,K=2,alpha = 0.05, nstart = 250, cstep1=2, cstep2=25, nkeep=5, equal.weights = FALSE, 
                 restr.fact=12,zero.tol = 1e-16,  trace = 0,   opt="HARD",
                 sol_ini_p = FALSE, sol_ini = NA ) 
    fin <- Sys.time()
    time5[n,d]<-as.numeric(fin-ini, units="secs")
    acc5[n,d]<-accuracy_clustering(tc$assig, true_labels)
    
    # ---------------------------------------------------------
    # 5) TCLUST (2 clusters; accurate trimming percentage; c=50 used to ensure flexibility)
    # ---------------------------------------------------------
    ini <- Sys.time()
    tc <- tclust(X,K=2,alpha = 0.05, nstart = 250, cstep1=2, cstep2=25, nkeep=5,  equal.weights = FALSE, 
                 restr.fact=50, zero.tol = 1e-16,  trace = 0,   opt="HARD",
                 sol_ini_p = FALSE, sol_ini = NA ) 
    fin <- Sys.time()
    time6[n,d]<-as.numeric(fin-ini, units="secs")
    acc6[n,d]<-accuracy_clustering(tc$assig, true_labels)
    
    # ---------------------------------------------------------
    # 6) TCLUST (2 clusters; accurate trimming percentage; c=150 used to ensure flexibility)
    # ---------------------------------------------------------
    ini <- Sys.time()
    tc <- tclust(X,K=2,alpha = 0.05, nstart = 250, cstep1=2, cstep2=25, nkeep=5,  equal.weights = FALSE, 
                 restr.fact=150, zero.tol = 1e-16,  trace = 0,   opt="HARD",
                 sol_ini_p = FALSE, sol_ini = NA ) 
    fin <- Sys.time()
    time7[n,d]<-as.numeric(fin-ini, units="secs")
    acc7[n,d]<-accuracy_clustering(tc$assig, true_labels)
    
    # ---------------------------------------------------------
    # 7) RLG (2 grupos, buenas dimensiones)
    # ---------------------------------------------------------
    Sys.setenv(OMP_NUM_THREADS="4")
    ini <- Sys.time()
    rg2 <- rlg(X, d=c(q1,q2), alpha=0.05, nstart=250, niter1=2, niter2=25, nkeep=5, parallel = TRUE, n.cores=20)
    fin <- Sys.time()
    time8[n,d]     <- as.numeric(fin-ini, units="secs")
    acc8[n,d] <- accuracy_clustering(rg2$cluster, true_labels)
    
  }
  
}


# --- VISUALIZATION HELPERS ---

# -----------------------------------------------------------
# FIGURE 3 (c)
# -----------------------------------------------------------

load("~/Simulaciones: Extending tclust to higher dimensions/sim_scen_3.RData")

# Boxplot of the accuracy obtained by TCLUST, RLG, and tHDDC (RLG and tHDDC 
#initialized with the true intrinsic dimensions) for different values of delta.
metodos <- c("tHDDC", "tclust", "rlg")
accuracy_list <- list(acc1, acc5, acc8)

df_accuracy <- do.call(rbind,
                       lapply(1:3, function(m) {
                         data.frame(
                           delta = rep(delta, each=N),
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
# FIGURE 4 (c)
# -----------------------------------------------------------

# Boxplot of the classification accuracy obtained by the tHDDC method when 
#using q1 and q2 equal to true values of the intrinsic dimensions, and also when
#using different values  q_{ini} for the  automatized determination of intrinsic
#dimensions based on the Cattell's approach, for each value of delta.

metodos <- c("tHDDC", "tHDDC q_{ini}=1", "tHDDC q_{ini}=5", "tHDDC q_{ini}=15")
accuracy_list <- list(acc1, acc2, acc3, acc4)

labels_metodos <- c(
  expression(tHDDC~paste("(", 10, ",", 5, ")")),
  expression(tHDDC~q[ini]==1),
  expression(tHDDC~q[ini]==5),
  expression(tHDDC~q[ini]==15)
)


df_accuracy <- do.call(rbind,
                       lapply(1:4, function(m) {
                         data.frame(
                           delta = rep(delta, each=N),
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
# TABLE 1 (Third row)
# -----------------------------------------------------------

#Ratio between the computation times of TCLUST and tHDDC for each value of the 
#separation parameter delta.

ratio <- colMeans(time5)/colMeans(time1)

tabla_resumen <- data.frame(
  delta = delta_vals,
  ratio_time=ratio)
tabla_resumen


# -----------------------------------------------------------
# FIGURE 5 (c)
# -----------------------------------------------------------

#Boxplot of the computation time for tHDDC different initial q_{ini} values for
#automated determination of the intrinsic dimensions, for each value of delta.

metodos <- c("tHDDC", "tHDDC q_{ini}=1", "tHDDC q_{ini}=5", "tHDDC q_{ini}=15")
accuracy_list <- list(time1, time2, time3, time4)

labels_metodos <- c(
  expression(tHDDC~paste("(", 10, ",", 5, ")")),
  expression(tHDDC~q[ini]==1),
  expression(tHDDC~q[ini]==5),
  expression(tHDDC~q[ini]==15)
)

df_time <- do.call(rbind,
                   lapply(1:4, function(m) {
                     data.frame(
                       delta = rep(delta, each=N),
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

p1 <- ggplot(df_time, aes(x = factor(delta), y = accuracy, fill = metodo, colour = metodo)) +
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
# FIGURE 6 (c)
# -----------------------------------------------------------

# Boxplots of the estimated larger intrinsic dimension obtained by the tHDDC 
#method using the Cattell-based procedure for different values of the separation
#parameter delta.

k1=max(qg)

err_catell <- matrix(0, nrow = N, ncol = D)
err_catell1 <- matrix(0, nrow = N, ncol = D)
err_catell2 <- matrix(0, nrow = N, ncol = D)
err_catell3 <- matrix(0, nrow = N, ncol = D)


err_catell <- apply(q_catell1, c(1, 2), max)-k1
err_catell1 <- apply(q_catell2, c(1, 2), max)-k1
err_catell2 <- apply(q_catell3, c(1, 2), max)-k1
err_catell3 <- apply(q_catell4, c(1, 2), max)-k1

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
        delta  = delta[d],
        metodo = name,
        error  = err_list[[name]][ ,d]
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
      expression(tHDDC~~"(" * 10 * "," * 5 * ")"),
      expression(tHDDC~~q[ini]==1),
      expression(tHDDC~~q[ini]==5),
      expression(tHDDC~~q[ini]==15)
    )
  ) +
  scale_colour_discrete(
    name = "Method",
    labels = c(
      expression(tHDDC~~"(" * 10 * "," * 5 * ")"),
      expression(tHDDC~~q[ini]==1),
      expression(tHDDC~~q[ini]==5),
      expression(tHDDC~~q[ini]==15)
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
    labels = seq(8, 20, by = 2)
  )+
  coord_cartesian(ylim = c(-0.5, 4.5))


# -----------------------------------------------------------
# FIGURE 7 (b)
# -----------------------------------------------------------

# Boxplots of the estimated smaller intrinsic dimension obtained by the tHDDC 
#method using the Cattell-based procedure for different values of the separation
#parameter delta.

k2=min(qg)

err_catell <- matrix(0, nrow = N, ncol = D)
err_catell1 <- matrix(0, nrow = N, ncol = D)
err_catell2 <- matrix(0, nrow = N, ncol = D)
err_catell3 <- matrix(0, nrow = N, ncol = D)


err_catell <- apply(q_catell1, c(1, 2), min)-k2
err_catell1 <- apply(q_catell2, c(1, 2), min)-k2
err_catell2 <- apply(q_catell3, c(1, 2), min)-k2
err_catell3 <- apply(q_catell4, c(1, 2), min)-k2

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
        delta  = delta[d],
        metodo = name,
        error  = err_list[[name]][ ,d]
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
      expression(tHDDC~~"(" * 10 * "," * 5 * ")"),
      expression(tHDDC~~q[ini]==1),
      expression(tHDDC~~q[ini]==5),
      expression(tHDDC~~q[ini]==15)
    )
  ) +
  scale_colour_discrete(
    name = "Method",
    labels = c(
      expression(tHDDC~~"(" * 10 * "," * 5 * ")"),
      expression(tHDDC~~q[ini]==1),
      expression(tHDDC~~q[ini]==5),
      expression(tHDDC~~q[ini]==15)
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
    labels = seq(3, 15, by = 2)
  )+
  coord_cartesian(ylim = c(-0.5, 6.5))
