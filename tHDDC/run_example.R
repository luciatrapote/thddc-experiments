###################################
#### IMPLEMENTATION EXAMPLE ####
###################################

set.seed(123)

##### Generate example data #####

# Data space dimension
p <- 200

# Group sizes
n1 <- 570
n2 <- 380

# First group
q1 <- 3
mu1 <- rep(1, p)
Delta1 <- diag(c(seq(9, 7, length = q1), rep(0.15, p - q1)))

# Second group
q2 <- 1
mu2 <- rep(-1, p)
Delta2 <- diag(c(rep(5, q2), rep(0.45, p - q2)))

qg <- c(q1, q2)

# Group 1
A1 <- matrix(rnorm((p + 1) * p), ncol = p)
B1 <- cov(A1)
U1 <- eigen(B1)$vector

X1 <- rep(1, n1) %*% t(mu1) +
  matrix(rnorm(n1 * p), nrow = n1) %*%
  diag(sqrt(diag(Delta1))) %*% U1

# Group 2
A2 <- matrix(rnorm((p + 1) * p), ncol = p)
B2 <- cov(A2)
U2 <- eigen(B2)$vector

X2 <- rep(1, n2) %*% t(mu2) +
  matrix(rnorm(n2 * p), nrow = n2) %*%
  diag(sqrt(diag(Delta2))) %*% U2

# Contamination / outliers
X.cont <- matrix(runif(50 * p, -2, 2), ncol = p)

# Final dataset
X <- rbind(X1, X2, X.cont)

# True labels
true_labels <- c(rep(1, n1), rep(2, n2), rep(3, 50))


##### Run the tHDDC algorithm #####

uX <- tHDDC_paralelizado_foreach(
  X,
  q = qg,
  alpha = 0.05,
  nstart = 250,
  cstep1 = 2,
  cstep2 = 25,
  nkeep = 5,
  equal.weights = TRUE,
  restr.fact = c(5, 3),
  zero.tol = 1e-16,
  trace = 0,
  opt = "HARD",
  sol_ini_p = FALSE,
  sol_ini = NA,
  catell_Y = TRUE,
  catell_threshold = 0.3,
  dim_max = 20
)

# Results: confusion matrix
table(true_labels, uX$assig)

# Estimated intrinsic dimensions
uX$q