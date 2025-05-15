require(geiger)
require(plotrix)
require(picante)
require(bipartite)
require(igraph)
require(ggplot2)
require(lme4)
require(car)
library(dplyr)

# Parameters
runs <- 10
power_H <- c(0.0001, 1, 2, 5)
power_L <- c(0.0001, 1, 2, 5)
complementarity <- c(0.25, 0.5, 0.75, 1)

RESULTADOS <- data.frame(real_eigenvalue = numeric(),
                         power_H = numeric(),
                         power_L = numeric(),
                         complementarity = numeric())

# Simulation Loop
for (p in complementarity) {
  for (o in power_H) {
    for (q in power_L) {
      for (m in 1:runs) {
        
        n_spe_H <- sample(5:8, 1)
        n_spe_L <- sample(5:8, 1)
        
        # Simulating Phylogenetic Trees
        tree_H <- sim.bdtree(b = 0.1, d = 0, stop = "taxa", n = n_spe_H, extinct = FALSE)
        tree_L <- sim.bdtree(b = 0.1, d = 0, stop = "taxa", n = n_spe_L, extinct = FALSE)
        
        tree_H$tip.label <- sprintf("H_%.3d", 1:length(tree_H$tip.label))
        tree_L$tip.label <- sprintf("L_%.3d", 1:length(tree_L$tip.label))
        
        # Trait Evolution
        trait_H <- rescale(rTraitCont(compute.brlen(tree_H, power = o), model = "BM"), c(0, 1))
        trait_L <- rescale(rTraitCont(compute.brlen(tree_L, power = q), model = "BM"), c(0, 1))
        
        # Constructing Interaction Matrix
        d_H <- matrix(runif(n_spe_H, max = p), nrow = n_spe_L, ncol = n_spe_H, byrow = TRUE)
        d_L <- matrix(runif(n_spe_L, max = p), nrow = n_spe_L, ncol = n_spe_H, byrow = FALSE)
        web <- outer(trait_L, trait_H, FUN = function(x, y) abs(x - y) < 0.5 * (d_H + d_L))
        web <- ifelse(web, 1, 0)
        
        # Ensuring No Empty Columns/Rows
        for (i in which(colSums(web) == 0)) web[sample(1:n_spe_L, 1), i] <- 1
        for (i in which(rowSums(web) == 0)) web[i, sample(1:n_spe_H, 1)] <- 1
        
        # Population Dynamics Simulation (Lotka-Volterra)
        n_spe <- n_spe_H + n_spe_L
        N <- runif(n_spe, 0.5, 1)  # Initial abundances
        
        interaction_matrix <- matrix(0, nrow = n_spe, ncol = n_spe)
        interaction_matrix[1:n_spe_L, (n_spe_L + 1):n_spe] <- exp(-abs(outer(trait_L, trait_H, "-")))  # Trait-based interactions
        interaction_matrix[(n_spe_L + 1):n_spe, 1:n_spe_L] <- t(interaction_matrix[1:n_spe_L, (n_spe_L + 1):n_spe])
        diag(interaction_matrix) <- -0.1
        
        time_steps <- 1000
        dt <- 0.01
        r <- runif(n_spe, 0.1, 0.2) * (1 + c(trait_L, trait_H))  # Trait-based growth rates
        
        for (t in 1:time_steps) {
          dN <- N * (r + interaction_matrix %*% N) * dt
          N <- pmax(0, N + dN)
        }
        
        # Stability Analysis (Eigenvalues of Jacobian)
        eigenvalues <- eigen(interaction_matrix)$values
        real_eigenvalue <- max(Re(eigenvalues))
        
        # Storing Results
        RESULTADOS <- rbind(RESULTADOS, data.frame(real_eigenvalue = real_eigenvalue,
                                                   power_H = o,
                                                   power_L = q,
                                                   complementarity = p))
      }
    }
  }
}
# Statistical test using linear mixed model (adjust as needed)
model <- lmer(real_eigenvalue ~ power_H * power_L * complementarity + (1|power_H) + (1|power_L), data = RESULTADOS)
anova_results <- Anova(model, type = "III")
print(anova_results)

# Aggregate data to avoid duplicate (x, y) points
RESULTADOS_grid <- RESULTADOS %>%
  group_by(power_H, power_L) %>%
  summarise(real_eigenvalue = mean(real_eigenvalue, na.rm = TRUE), .groups = "drop")

# Create the contour plot with the aggregated data
contour_plot <- ggplot(RESULTADOS_grid, aes(x = power_H, y = power_L, z = real_eigenvalue)) +
  geom_contour_filled() +
  labs(
    title = "Phylogenetic Signal and Stability",
    x = "Phylogenetic Signal (Power H)",
    y = "Phylogenetic Signal (Power L)",  # Updated y-axis label
    fill = expression(Re(lambda))
  ) +
  theme_minimal()

# Create the violin plot
violin_plot <- ggplot(RESULTADOS, aes(x = as.factor(power_H), y = real_eigenvalue, fill = as.factor(complementarity))) +
  geom_violin() +
  facet_wrap(~ as.factor(power_L)) +
  labs(title = "Stability",
       x = "Phylogenetic Signal Higher trophic level",
       y = expression(Re(lambda)),
       fill = "Complementarity") +
  theme_minimal() +
  coord_flip()

###Empirical Cumulative Distribution Functions
distributions=ggplot(RESULTADOS, aes(x = real_eigenvalue, color = interaction(as.factor(power_H), as.factor(power_L)))) +
  stat_ecdf() +
  labs(title = "Phylogenetic Signal and Stability",
       x = "Real Eigenvalue (Stability)",
       y = "Cumulative Probability",
       color = "Phylogenetic Signal (Power H) : Phylogenetic Signal (Power L)") +
  theme_minimal()

