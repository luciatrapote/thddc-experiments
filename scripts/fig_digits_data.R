################################################################################
#           Visualization of Clustering Centers, Loading Vectors, and          #
#          Robustness Diagnostics (SD/OD Plots and Silhouette Analysis)        #
################################################################################

# -------------------------------------------------------------------------
# 1. LOAD REQUIRED LIBRERIES
# -------------------------------------------------------------------------

if (!require("robustbase")) install.packages("robustbase")
if (!require("mclust")) install.packages("mclust")
library(robustbase)
library(mclust)


# -------------------------------------------------------------------------
# 2. DATA LOADING AND PREPARATION
# -------------------------------------------------------------------------

# uX <- readRDS("tHDDC_results_uX.rds")
uX <- readRDS("~/Simulaciones: Extending tclust to higher dimensions/tHDDC_results_uX.rds")
# ZZZ_df <- read.csv("data_ZZZ_df.csv")
ZZZ_df <- read.csv("~/Simulaciones: Extending tclust to higher dimensions/Simulaciones: Digits Data/data_ZZZ_df.csv")

# --- Ground Truth Label Preparation ---
# We assess if outliers are correctly trimmed and if digits 3, 5, and 8 
# are properly identified. 
# Labels 0, 2, 4, and 9 are merged into 'out1' (natural noise/other digits).
# Label 999 is renamed to 'out2' (artificially generated outlier patterns).

ZZZ_df$cls_fusionada <- ZZZ_df$cls
ZZZ_df$cls_fusionada[ZZZ_df$cls_fusionada %in% c(0, 2, 4, 9)] <- "out1"
ZZZ_df$cls_fusionada[ZZZ_df$cls_fusionada == 999] <- "out2"


# Data dimensions
ZZZ <- as.matrix(ZZZ_df[, 1:256])
n <- dim(ZZZ)[1]
p <- dim(ZZZ)[2]
q <- uX$q


# -------------------------------------------------------------------------
# 3. FIGURE 8: CLUSTER CENTERS
# -------------------------------------------------------------------------

par(mfrow=c(1,3), mar = c(2.2, 2.2, 2.2, 2.2))
for(k in c(2,3,1)){
  image(matrix(uX$center[k,],ncol=16,byrow=TRUE)) 
}


# -------------------------------------------------------------------------
# 4. FIGURE 9: LOADING VECTORS
# -------------------------------------------------------------------------

par(mfrow=c(2,4), mar = c(2.2, 2.2, 2.2, 2.2))
for (k in c(2,3,1)){
  for(i in 1:8){
    m <- matrix(uX$evectors[[k]][,i], ncol = 16, byrow = TRUE)
    
    colores <- colorRampPalette(c("#2196F3", "white", "#FF9800"))(100)
    
    lim <- max(abs(uX$evectors[[k]]))
    breaks <- seq(-lim, lim, length.out = length(colores) + 1)
    
    image(m, col = colores, breaks = breaks, main = "")
  }
}


# -------------------------------------------------------------------------
# 5. FIGURE 10: TRIMMED OUTLIERS 
# -------------------------------------------------------------------------

par(mfrow=c(5,4), mar = c(2, 2, 2, 2))

indices<-c(1904,1887,1787,1939, 1762,1761,1780,1773,  1889,1786,1926,1828, 1923,1941,1836,1899, 1958,1968,1978,1988 )
for (i in 1:20){
  image(matrix(ZZZ[indices[i],],ncol=16,byrow=TRUE), axes=FALSE)
}


# -------------------------------------------------------------------------
#  6. FIGURE 11 (a), (b) and (c): SD vs OD PLOTS 
# -------------------------------------------------------------------------

# Pre-calculating cluster assignments
K <- length(q)
DK <- matrix(,nrow=n,ncol=K)

for (k in 1:K){
  A <- matrix(uX$evectors[[k]],ncol=q[k])
  Xcen <- ZZZ - rep(1,n) %*% t(uX$center[k,])
  TT <- (Xcen %*% A)
  A1 <- rowSums( (TT^2) %*% diag((uX$evalues[1:q[k],k])^(-1),q[k],q[k]) )
  A2 <- ( ( uX$evalues[q[k]+1,k] )^(-1) ) * rowSums( (Xcen %*% (diag(p) - A %*% t(A)))^2)
  A3 <- sum( log(uX$evalues[1:q[k],k]) )
  A4 <- ( p - q[k] ) * log( uX$evalues[q[k]+1,k] )
  A5 <- p * log(2*pi)
  DK[,k] <- log(uX$cw[k]*exp(- 0.5 * ( A1 + A2 + A3 + A4 + A5 )))
}

todos <- apply(DK,1,which.max)

# Plotting Loop
digit_colors <- c("#6B6B6B", "black", "#7B3294", "#1F77B4", "#FF8C00",
                  "#2CA02C", "black", "black", "#D62728", "#8C564B","#FF1493")

aux_num <- as.numeric(as.character(ZZZ_df$cls))

for (k in c(2,3,1)) {
  
  X.k <- ZZZ[todos==k,]
  n.k <- dim(X.k)[1]
  
  A <- matrix(uX$evectors[[k]], ncol=q[k])
  Xcen <- X.k - rep(1,n.k) %*% t(uX$center[k,])
  TT <- (Xcen %*% A)
  
  A1 <- rowSums((TT^2) %*% diag((uX$evalues[1:q[k],k])^(-1), q[k], q[k]))
  
  A2 <- (uX$evalues[q[k]+1,k])^(-1) *
    rowSums((Xcen %*% (diag(p) - A %*% t(A)))^2)
  
  A3 <- sum(log(uX$evalues[1:q[k],k]))
  A4 <- (p - q[k]) * log(uX$evalues[q[k]+1,k])
  A5 <- p * log(2*pi)
  
  ff <- uX$cw[k]*exp(-0.5 * (A1 + A2 + A3 + A4 + A5))
  
  SD <- sqrt(A1)
  OD <- sqrt(A2)
  
  idx_cluster <- which(todos == k)
  cols_k <- digit_colors[aux_num[idx_cluster] + 1]
  
  # Set transparency for overlapping points
  if(k == 3) cols_k[cols_k == "#2CA02C"] <- adjustcolor("#2CA02C", alpha.f = 0.4)
  if(k == 1) cols_k[cols_k == "#D62728"] <- adjustcolor("#D62728", alpha.f = 0.4)
  if(k == 2) cols_k[cols_k == "#1F77B4"] <- adjustcolor("#1F77B4", alpha.f = 0.4)
  
  # Visualization Setup
  par(mfrow=c(1,1), mar=c(5.1,4.1,4.1,2.1))
  plot(SD, OD, type="n", xlab="SD_i", ylab="OD_i")
  
  
  # Highlight points 
  if(k == 3){
    indices <- c(772,1205,1074,1030,712)
    highlight_color <- "#2CA02C"
  }
  
  if(k == 1){
    indices <- c(1310,1734)
    highlight_color <- "#D62728"
  }
  
  if(k == 2){
    indices <- c(44,617)
    highlight_color <- "#1F77B4"
  }
  
  pos_highlight <- which(idx_cluster %in% indices)
  highlight <- rep(FALSE, length(SD))
  highlight[pos_highlight] <- TRUE
  
  # Identify observations of out2
  is_star <- aux[idx_cluster] == 999
  
  # Plot observations (no out2, no highlight)
  text(SD[!highlight & !is_star],
       OD[!highlight & !is_star],
       aux[idx_cluster[!highlight & !is_star]],
       col = cols_k[!highlight & !is_star],
       cex = 1.4)
  
  # Plot highlight observations (no out2)
  if(length(pos_highlight) > 0){
    text(SD[highlight & !is_star],
         OD[highlight & !is_star],
         aux[idx_cluster[highlight & !is_star]],
         col = highlight_color,
         cex = 1.6)
  }
  
  # Mark artificial outliers (out2)
  if(any(is_star)){
    points(SD[is_star],
           OD[is_star],
           pch = 4,
           col = "black",
           cex = 1.2)
  }
  
  # Robust Cutoff (Hubert et al., 2005)
  abline(v = sqrt(qchisq(0.975, q[k])), lty=2) 
  
  A2.3.2 <- OD^(2/3)
  mcd_result <- covMcd(A2.3.2)
  
  abline(h = (mcd_result$center +
                sqrt(mcd_result$cov)*qnorm(0.975))^(3/2),
         lty=2)
}


# -------------------------------------------------------------------------
#  SD vs OD PLOTS  (FIGURE 11 (d))
# -------------------------------------------------------------------------

par(mfrow=c(3,4), mar = c(3, 2, 0.5, 2))

indices<-c(44,617,772,1205, 1074,1030,712,1438, 1310,1734,22,1928)
for (i in 1:12){
  image(matrix(ZZZ[indices[i],],ncol=16,byrow=TRUE), axes=FALSE)
}


# -------------------------------------------------------------------------
#  SILHOUETTE PLOT (FIGURE 12 (a))
# -------------------------------------------------------------------------

# Calculation of Discriminant Factors (DF)
threshold <- sort(apply(DK,1,max))[sum(uX$assig==0)+1]

DK.order <- t(apply(DK, 1, sort))
BF <- DK.order[,K-1]-DK.order[,K]
BF[uX$assig==0] <- DK.order[uX$assig==0,K] - threshold

BF.ordenados <- c()

vector <- c(1,3,2,0)
for (i in vector){
  BF.ordenados <-c(BF.ordenados,sort(BF[uX$assig==i]))
}

# Handling visualization limits
BF.ordenados[BF.ordenados==-Inf]=-600 # Capping infinite values

BF.ordenados=-BF.ordenados

# Plotting
par(mfrow=c(1,1))

# Highlight specific observations to analyze their classification confidence
# i0-i7: Selected indices for detailed visual inspection
i0 <- 392
i1 <- 104
i2 <- 1438
i3 <- 361
i4 <- 1939
#i5 <- 489
i5 <- 1899
i6 <- 1995
i7 <- 1958

# Initialize vectors to store ordered positions and original indices
posiciones <- c()
indices_originales <- c()



for (i in vector){
  idx_cluster <- which(uX$assig == i)
  BF_cluster <- BF[idx_cluster]
  orden_cluster <- order(BF_cluster)
  
  posiciones <- c(posiciones, BF_cluster[orden_cluster])
  indices_originales <- c(indices_originales, idx_cluster[orden_cluster])
}

idx_cluster_all <- c()

for (i in vector){
  
  idx_cluster <- which(uX$assig == i)
  idx_cluster_all <- c(idx_cluster_all, idx_cluster)
  
}
pos_i0 <- which(indices_originales == i0)
pos_i1 <- which(indices_originales == i1)
pos_i2 <- which(indices_originales == i2)
pos_i3 <- which(indices_originales == i3)
pos_i4 <- which(indices_originales == i4)
pos_i5 <- which(indices_originales == i5)
pos_i6 <- which(indices_originales == i6)
pos_i7 <- which(indices_originales == i7)



color_aux<-c(sum(uX$assig==1),sum(uX$assig==3),sum(uX$assig==2),sum(uX$assig==0))
colores <- rep(c(2,3,4,1), times = color_aux) 

bp <- barplot(BF.ordenados,
              horiz = TRUE,
              xlim = c(0,650),
              col = colores,
              border = NA)

# SILHOUETTE BARPLOT CUSTOMIZATION

# Redraw SPECIFIC bars with increased thickness for emphasis
# Labeling points with bold indices (1-8) for reference in the paper text

rect(xleft  = 0,
     xright = BF.ordenados[pos_i0],
     ybottom = bp[pos_i0] - 0.6,
     ytop    = bp[pos_i0] + 0.6,
     col = "red",
     border = "darkblue",
     lwd = 3)

points(BF.ordenados[pos_i0], bp[pos_i0], 
       pch = 19, cex = 1, col = "darkblue")

text(x = BF.ordenados[pos_i0]+1,
     y = bp[pos_i0],
     labels = expression(bold(1)),
     pos = 4,          
     offset = 0.4,
     col = "darkblue",
     cex = 0.9)


#########
rect(xleft  = 0,
     xright = BF.ordenados[pos_i1],
     ybottom = bp[pos_i1] - 0.6,
     ytop    = bp[pos_i1] + 0.6,
     col = "red",
     border = "darkblue",
     lwd = 3)

points(BF.ordenados[pos_i1], bp[pos_i1], 
       pch = 19, cex = 1, col = "darkblue")

text(x = BF.ordenados[pos_i1]+10,
     y = bp[pos_i1]+45,
     labels = expression(bold(5)),
     pos = 4, 
     offset = 0.4,
     col = "darkblue",
     cex = 0.9)

#########
rect(xleft  = 0,
     xright = BF.ordenados[pos_i2],
     ybottom = bp[pos_i2] - 0.6,
     ytop    = bp[pos_i2] + 0.6,
     col = "red",
     border = "darkgreen",
     lwd = 3)

points(BF.ordenados[pos_i2], bp[pos_i2], 
       pch = 19, cex = 1, col = "darkgreen")

text(x = BF.ordenados[pos_i2]-5,
     y = bp[pos_i2]+10,
     labels = expression(bold(2)),
     pos = 4, 
     offset = 0.4,
     col = "darkgreen",
     cex = 0.9)

#########
rect(xleft  = 0,
     xright = BF.ordenados[pos_i3],
     ybottom = bp[pos_i3] - 0.6,
     ytop    = bp[pos_i3] + 0.6,
     col = "red",
     border = "#8B0000",
     lwd = 3)

points(BF.ordenados[pos_i3], bp[pos_i3], 
       pch = 19, cex = 1, col = "#8B0000")

text(x = BF.ordenados[pos_i3]-2,
     y = bp[pos_i3],
     labels = expression(bold(3)),
     pos = 4, 
     offset = 0.4,
     col = "#8B0000",
     cex = 0.9)

######
rect(xleft  = 0,
     xright = BF.ordenados[pos_i4],
     ybottom = bp[pos_i4] - 0.6,
     ytop    = bp[pos_i4] + 0.6,
     col = "red",
     border = "orange",
     lwd = 3)

points(BF.ordenados[pos_i4], bp[pos_i4],
       pch = 19, cex = 1, col = "orange")

text(x = BF.ordenados[pos_i4]+10,
     y = bp[pos_i4]+50,
     labels = expression(bold(6)),
     pos = 4,          
     offset = 0.4,
     col = "orange",
     cex = 0.9)

######
rect(xleft  = 0,
     xright = BF.ordenados[pos_i5],
     ybottom = bp[pos_i5] - 0.6,
     ytop    = bp[pos_i5] + 0.6,
     col = "red",
     border = "orange",
     lwd = 3)

points(BF.ordenados[pos_i5], bp[pos_i5], 
       pch = 19, cex = 1, col =  "orange")

text(x = BF.ordenados[pos_i5]-2,
     y = bp[pos_i5]+20,
     labels = expression(bold(4)),
     pos = 4,          
     offset = 0.4,
     col =  "orange",
     cex = 0.9)

######
rect(xleft  = 0,
     xright = BF.ordenados[pos_i6],
     ybottom = bp[pos_i6] - 0.6,
     ytop    = bp[pos_i6] + 0.6,
     col = "red",
     border = "#FF1493",
     lwd = 3)

points(BF.ordenados[pos_i6], bp[pos_i6], 
       pch = 19, cex = 1, col =  "#FF1493")

text(x = BF.ordenados[pos_i6]-1,
     y = bp[pos_i6]+20,
     labels = expression(bold(7)),
     pos = 4,          
     offset = 0.4,
     col =  "#FF1493",
     cex = 0.9)


######
rect(xleft  = 0,
     xright = BF.ordenados[pos_i7],
     ybottom = bp[pos_i7] - 0.6,
     ytop    = bp[pos_i7] + 0.6,
     col = "red",
     border = "#FF1493",
     lwd = 3)

points(BF.ordenados[pos_i7], bp[pos_i7], 
       pch = 19, cex = 1, col =  "#FF1493")

text(x = BF.ordenados[pos_i7]-1,
     y = bp[pos_i7]-20,
     labels = expression(bold(8)),
     pos = 4,          
     offset = 0.4,
     col =  "#FF1493",
     cex = 0.9)


# -------------------------------------------------------------------------
#  SILHOUETTE PLOT (FIGURE 12 (b))
# -------------------------------------------------------------------------

par(mfrow=c(4,2), mar = c(3, 2, 0.5, 2))

indices <- c(392,1438,361,1899,104,1939,1995,1958)

for (i in 1:8){
  image(matrix(ZZZ[indices[i],], ncol=16, byrow=TRUE), axes=FALSE)
}

