# ===================================================================
# --- 1. LOAD REQUIRED PACKAGES ---
# ===================================================================
if (!require("geiger")) install.packages("geiger")
if (!require("ape")) install.packages("ape")
if (!require("phytools")) install.packages("phytools")
if (!require("RColorBrewer")) install.packages("RColorBrewer")

library(geiger)
library(ape)
library(phytools)
library(RColorBrewer)


# ===================================================================
# --- 2. HELPER FUNCTION ---
# ===================================================================
rescale_vector <- function(x, new_range) {
  if (length(unique(x)) == 1) return(rep(mean(new_range), length(x)))
  old_min <- min(x, na.rm = TRUE)
  old_max <- max(x, na.rm = TRUE)
  new_min <- new_range[1]
  new_max <- new_range[2]
  return(new_min + (x - old_min) * (new_max - new_min) / (old_max - old_min))
}


# ===================================================================
# --- 3. SET YOUR PARAMETERS FOR THE FIGURE HERE ---
# ===================================================================
set.seed(12) 

n_spe_L <- 10
n_spe_H <- 10 
alpha_L <- 0.001
alpha_H <- 0.001
p_match <- 0.2


# ===================================================================
# --- 4. SIMULATE DATA FOR A SINGLE RUN ---
# (This section is unchanged and correct)
# ===================================================================
tree_L <- sim.bdtree(b = 0.1, d = 0, stop = "taxa", n = n_spe_L)
tree_H <- sim.bdtree(b = 0.1, d = 0, stop = "taxa", n = n_spe_H)
tree_L$tip.label <- sprintf("L_%.3d", 1:n_spe_L)
tree_H$tip.label <- sprintf("H_%.3d", 1:n_spe_H)

trait_L_raw <- rTraitCont(tree_L, model = "OU", alpha = alpha_L)
trait_H_raw <- rTraitCont(tree_H, model = "OU", alpha = alpha_H)
trait_L <- rescale_vector(trait_L_raw, c(0, 1))
trait_H <- rescale_vector(trait_H_raw, c(0, 1))

d_matrix <- matrix(runif(n_spe_L * n_spe_H, min = 0, max = p_match), nrow = n_spe_L, ncol = n_spe_H)
web <- outer(trait_L, trait_H, FUN = function(x, y) abs(x - y) < d_matrix)
web <- ifelse(web, 1, 0)
dimnames(web) <- list(names(trait_L), names(trait_H))

empty_cols_H <- which(colSums(web) == 0)
if (length(empty_cols_H) > 0) {
  for (h_idx in empty_cols_H) web[sample(1:n_spe_L, 1), h_idx] <- 1
}
empty_rows_L <- which(rowSums(web) == 0)
if (length(empty_rows_L) > 0) {
  for (l_idx in empty_rows_L) web[l_idx, sample(1:n_spe_H, 1)] <- 1
}


# ===================================================================
# --- 5. PREPARE DATA FOR PLOTTING (ALIGNMENT STEP) ---
# (This section is also unchanged and correct)
# ===================================================================
tree_L_reordered <- reorder(tree_L, "postorder")
tip_indices_L <- tree_L_reordered$edge[tree_L_reordered$edge[, 2] <= n_spe_L, 2]
tip_order_L <- tree_L_reordered$tip.label[tip_indices_L]

tree_H_reordered <- reorder(tree_H, "postorder")
tip_indices_H <- tree_H_reordered$edge[tree_H_reordered$edge[, 2] <= n_spe_H, 2]
tip_order_H <- tree_H_reordered$tip.label[tip_indices_H]

web_reordered <- web[tip_order_L, tip_order_H]
trait_L_reordered <- trait_L[tip_order_L]
trait_H_reordered <- trait_H[tip_order_H]

color_palette <- colorRampPalette(c("blue", "yellow", "red"))(100)
trait_colors_L <- color_palette[cut(trait_L_reordered, breaks = 100, labels = FALSE)]
trait_colors_H <- color_palette[cut(trait_H_reordered, breaks = 100, labels = FALSE)]


# ===================================================================
# --- 6. CREATE THE FINAL PLOT (CORRECTED) ---
# This version fixes the image() dimension error.
# ===================================================================

# --- Define the layout ---
layout_matrix <- matrix(c(0, 1, 2, 3), nrow=2, ncol=2)
layout(layout_matrix, widths=c(1, 3), heights=c(1, 3))

# --- Plotting parameters ---
phylo_edge_width <- 2
trait_point_size <- 2.0

# --- Plot 1: Top Phylogeny (H) ---
par(mar = c(0, 0, 1, 1))
plot(tree_H_reordered, direction = "rightwards", show.tip.label = FALSE, edge.width = phylo_edge_width)
points(y = 1:n_spe_H,
       x = rep(max(node.depth.edgelength(tree_H_reordered)), n_spe_H),
       pch = 21, bg = trait_colors_H, col = "black", cex = trait_point_size)

# --- Plot 2: Left Phylogeny (L) and Text ---
par(mar = c(1, 1, 0, 0))
plot(tree_L_reordered, direction = "downwards", show.tip.label = FALSE, edge.width = phylo_edge_width)
points(x = 1:n_spe_L, y = rep(0, n_spe_L), pch = 21, bg = trait_colors_L,
       col = "black", cex = trait_point_size)
text(x = 1, y = max(node.depth.edgelength(tree_L_reordered)) * 0.9,
     labels = paste0("Trait matching = ", p_match, " | \u03B1(resource) = ", alpha_L, ", \u03B1(consumer) = ", alpha_H),
     pos = 4, cex = 1.2)


# --- Plot 3: Interaction Web (CORRECTED) ---
par(mar = c(1, 0, 0, 1))
#
# ** THE FIX IS HERE **
# We do NOT transpose web_reordered. It already has the correct dimensions:
# rows = n_spe_L (for x-axis), cols = n_spe_H (for y-axis).
#
image(x = 1:n_spe_L, 
      y = 1:n_spe_H, 
      z = web_reordered,
      col = c("white", "black"),
      axes = FALSE, xlab = "", ylab = "")
box()

# --- Reset the plotting device to default ---
layout(1)
par(mar = c(5.1, 4.1, 4.1, 2.1))

