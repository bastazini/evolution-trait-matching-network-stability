require(geiger)
require(plotrix)
require(picante)
require(bipartite)

# Complementarity values
complementarity_values <- c(0.25, 0.5, 0.75, 1)

# Function to generate synthetic data
geradados <- function(n_spe_H = 10, n_spe_L = 10, power_H = 0.001, power_L = 0.001, comp = 1) {
  tree_H <- sim.bdtree(b = 0.1, d = 0, stop = "taxa", n = n_spe_H, extinct = FALSE)
  tree_H$edge.length[tree_H$edge.length == 0] <- 0.01
  tree_L <- sim.bdtree(b = 0.1, d = 0, stop = "taxa", n = n_spe_L, extinct = FALSE)
  tree_L$edge.length[tree_L$edge.length == 0] <- 0.01
  
  tree_H$tip.label <- sprintf("H_%.3d", 1:n_spe_H)
  tree_L$tip.label <- sprintf("L_%.3d", 1:n_spe_L)
  
  trait_H <- matrix(rTraitCont(compute.brlen(tree_H, power = power_H), model = "BM"), n_spe_H, 1)
  trait_L <- matrix(rTraitCont(compute.brlen(tree_L, power = power_L), model = "BM"), n_spe_L, 1)
  trait_H <- apply(trait_H, 2, function(x) rescale(x, c(0, 1)))
  trait_L <- apply(trait_L, 2, function(x) rescale(x, c(0, 1)))
  rownames(trait_H) <- tree_H$tip.label
  rownames(trait_L) <- tree_L$tip.label
  
  d_H <- matrix(runif(n_spe_H, max = 0.25), n_spe_H, 1)
  d_L <- matrix(runif(n_spe_L, max = 0.25), n_spe_L, 1)
  
  web <- matrix(0, n_spe_L, n_spe_H)
  for (i in 1:n_spe_L) {
    for (j in 1:n_spe_H) {
      II <- abs(trait_H[j] - trait_L[i])
      III <- comp * 0.5 * (d_H[j] + d_L[i])
      web[i, j] <- ifelse(II < III, 1, 0)
    }
  }
  
  z_H <- which(colSums(web) == 0)
  if (length(z_H) > 0) {
    for (i in z_H) web[sample(1:n_spe_L, 1), i] <- 1
  }
  z_L <- which(rowSums(web) == 0)
  if (length(z_L) > 0) {
    for (i in z_L) web[i, sample(1:n_spe_H, 1)] <- 1
  }
  
  list(tree_H = tree_H, tree_L = tree_L, trait_H = trait_H, trait_L = trait_L, web = web)
}

# ---- Set up a full panel layout (3 rows x 3 cols per panel) x 4 panels = 12 rows x 3 cols ----
layout(matrix(rep(c(0, 0, 1, 0, 0, 2, 3, 4, 5,
                    0, 0, 6, 0, 0, 7, 8, 9,10,
                    0, 0,11, 0, 0,12,13,14,15,
                    0, 0,16, 0, 0,17,18,19,20), each = 1), 
              nrow = 12, byrow = TRUE),
       widths = c(0.4, 0.4, 1), heights = rep(c(0.4, 0.4, 1), 4))
par(mar = c(0.4, 0.4, 0.4, 0.4))

# Loop through complementarity values
for (comp in complementarity_values) {
  DADOS <- geradados(10, 10, power_H = 0.001, power_L = 0.001, comp = comp)
  
  # Tree H
  plot(DADOS$tree_H, show.tip.label = FALSE, direction = "downwards")
  
  # Trait H
  plot(DADOS$trait_H[, 1], rep(1, length(DADOS$trait_H[, 1])),
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", type = "n",
       ylim = c(0.9, 1.1), xlim = c(1, length(DADOS$trait_H[, 1])), bty = "n")
  points(1:length(DADOS$trait_H[, 1]), rep(1, length(DADOS$trait_H[, 1])),
         cex = DADOS$trait_H[, 1] + 1, pch = 19)
  
  # Tree L
  plot(DADOS$tree_L, show.tip.label = FALSE)
  
  # Trait L
  plot(DADOS$trait_L[, 1], rep(1, length(DADOS$trait_L[, 1])),
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", type = "n",
       xlim = c(0.9, 1.1), ylim = c(1, length(DADOS$trait_L[, 1])), bty = "n")
  points(rep(1, length(DADOS$trait_L[, 1])), 1:length(DADOS$trait_L[, 1]),
         cex = DADOS$trait_L[, 1] + 1, pch = 19)
  
  # Matrix
  color2D.matplot(DADOS$web, ylab = "", xlab = "", yrev = FALSE,
                  cs1 = c(1, 0), cs2 = c(1, 0), cs3 = c(1, 0),
                  border = "white", xaxt = "n", yaxt = "n", axes = FALSE)
  mtext(paste("Complementarity =", comp), side = 3, line = -1.5, adj = 0.05, cex = 0.8)
}
