rm(list = ls()) # Clear workspace
require(geiger)
require(picante)
require(bipartite)
require(igraph)
require(ggplot2)
library(dplyr)

# Function to rescale
rescale_vector <- function(x, new_range) {
  if (length(unique(x)) == 1) return(rep(mean(new_range), length(x)))
  old_min <- min(x, na.rm = TRUE)
  old_max <- max(x, na.rm = TRUE)
  new_min <- new_range[1]
  new_max <- new_range[2]
  return(new_min + (x - old_min) * (new_max - new_min) / (old_max - old_min))
}

# --- Simulation Parameters ---
runs <- 500# Number of runs per parameter combination
power_H <- c(0.0001, 1, 2, 5)
power_L <- c(0.0001, 1, 2, 5)
complementarity_levels <- c(0.25, 0.5, 0.75, 1)


time_steps_lv <- 2000 # Time steps for LV simulation
dt_lv <- 0.01         # Time step size for LV simulation
intraspecific_comp <- -0.5 # Diagonal elements for interaction_matrix_lv
initial_abundance_range <- c(0.1, 0.5) # Range for initial N
r_base_range <- c(0.1, 0.3) # Range for base growth rate
max_allowed_abundance <- 1e6 # To stop runaway simulations

set.seed(123)

RESULTADOS <- data.frame(max_Re_eigen_jacobian = numeric(),
                         power_H = numeric(),
                         power_L = numeric(),
                         complementarity = numeric(),
                         n_spe_H_val = numeric(),
                         n_spe_L_val = numeric(),
                         num_persistent_species = numeric(),
                         simulation_converged = logical()) # To track issues

contagem <- 0
total_iterations <- length(complementarity_levels) * length(power_H) * length(power_L) * runs
print(paste("Total iterations planned:", total_iterations))


for (p_comp in complementarity_levels) {
  for (o_pH in power_H) {
    for (q_pL in power_L) {
      for (m_run in 1:runs) {
        contagem <- contagem + 1
        if (contagem %% 10 == 0) {
          print(paste("Iteration:", contagem, "of", total_iterations,
                      "| p:", p_comp, "o:", o_pH, "q:", q_pL))
        }
        
        n_spe_H <- sample(10:50, 1) # Slightly smaller for speed
        n_spe_L <- sample(10:50, 1)
        
        tree_H <- sim.bdtree(b = 0.1, d = 0, stop = "taxa", n = n_spe_H, extinct = FALSE)
        tree_L <- sim.bdtree(b = 0.1, d = 0, stop = "taxa", n = n_spe_L, extinct = FALSE)
        tree_H$tip.label <- sprintf("H_%.3d", 1:n_spe_H)
        tree_L$tip.label <- sprintf("L_%.3d", 1:n_spe_L)
        
        trait_H_raw <- rTraitCont(compute.brlen(tree_H, power = o_pH), model = "BM")
        trait_H <- rescale_vector(trait_H_raw, c(0, 1))
        trait_L_raw <- rTraitCont(compute.brlen(tree_L, power = q_pL), model = "BM")
        trait_L <- rescale_vector(trait_L_raw, c(0, 1))
        
        d_H_matrix <- matrix(runif(n_spe_L * n_spe_H, min = 0, max = p_comp), nrow = n_spe_L, ncol = n_spe_H)
        d_L_matrix <- matrix(runif(n_spe_L * n_spe_H, min = 0, max = p_comp), nrow = n_spe_L, ncol = n_spe_H)
        
        web <- outer(trait_L, trait_H, FUN = function(x, y) abs(x - y) < 0.5 * (d_L_matrix + d_H_matrix))
        web <- ifelse(web, 1, 0)
        
        empty_cols_H <- which(colSums(web) == 0)
        if (length(empty_cols_H) > 0) {
          for (h_idx in empty_cols_H) web[sample(1:n_spe_L, 1), h_idx] <- 1
        }
        empty_rows_L <- which(rowSums(web) == 0)
        if (length(empty_rows_L) > 0) {
          for (l_idx in empty_rows_L) web[l_idx, sample(1:n_spe_H, 1)] <- 1
        }
        
        n_total_spe <- n_spe_L + n_spe_H
        A_interaction_matrix <- matrix(0, nrow = n_total_spe, ncol = n_total_spe)
        
        consumer_degrees_H <- colSums(web)
        consumer_degrees_H[consumer_degrees_H == 0] <- 1
        resource_degrees_L <- rowSums(web)
        resource_degrees_L[resource_degrees_L == 0] <- 1
        
        for (i_row_L in 1:n_spe_L) {
          for (j_col_H in 1:n_spe_H) {
            if (web[i_row_L, j_col_H] == 1) {
              lv_L_idx <- i_row_L
              lv_H_idx <- n_spe_L + j_col_H
              A_interaction_matrix[lv_L_idx, lv_H_idx] <- -1 / consumer_degrees_H[j_col_H]
              A_interaction_matrix[lv_H_idx, lv_L_idx] <- +1 / resource_degrees_L[i_row_L]
            }
          }
        }
        diag(A_interaction_matrix) <- intraspecific_comp
        
        # --- Lotka-Volterra Simulation ---
        N <- runif(n_total_spe, min = initial_abundance_range[1], max = initial_abundance_range[2])
        r_base_values <- runif(n_total_spe, min = r_base_range[1], max = r_base_range[2])
        # Assuming traits for r are ordered L then H
        all_traits_for_r <- c(trait_L, trait_H) # Make sure this order matches N and A_interaction_matrix
        r <- r_base_values * (1 + all_traits_for_r) 
        
        simulation_ok <- TRUE
        for (t_step in 1:time_steps_lv) {
          growth_rates_per_capita <- r + (A_interaction_matrix %*% N)
          dN <- N * growth_rates_per_capita * dt_lv
          N <- N + dN
          N[N < 1e-6] <- 0 # Threshold for extinction / numerical stability
          
          if(any(is.na(N)) || any(is.infinite(N)) || any(N > max_allowed_abundance)) {
            simulation_ok <- FALSE
            # print(paste("Run", contagem, "diverged or became NA/Inf at step", t_step))
            break 
          }
        }
        
        max_Re_val_jacobian <- NA # Default to NA if issues
        num_persistent <- 0
        
        if (simulation_ok && sum(N > 0) > 0) { # If simulation ran and at least one species persists
          persistent_species_indices <- which(N > 1e-6) # Consider species with N > 0 as persistent
          num_persistent <- length(persistent_species_indices)
          
          if (num_persistent > 0) {
            # Construct Jacobian for the sub-community of persistent species
            N_persistent <- N[persistent_species_indices]
            A_persistent <- A_interaction_matrix[persistent_species_indices, persistent_species_indices]
            r_persistent <- r[persistent_species_indices]
            
            # Jacobian J_ij = N_i * a_ij (for i != j) and J_ii = r_i + sum_k(a_ik N_k) + N_i * a_ii
            # Equivalent to: J = diag(N_persistent) %*% A_persistent + diag(r_persistent + A_persistent %*% N_persistent)
            # Simpler: J = diag(N_persistent) %*% A_persistent + diag(growth_rates_per_capita[persistent_species_indices])
            # where growth_rates_per_capita are near zero at equilibrium.
            # The standard Jacobian is: J_ik = delta_ik * (r_i + sum_j a_ij N_j) + N_i * a_ik
            # delta_ik is Kronecker delta (1 if i=k, 0 otherwise)
            
            Jacobian <- matrix(0, nrow = num_persistent, ncol = num_persistent)
            for (i_idx in 1:num_persistent) {
              for (k_idx in 1:num_persistent) {
                val_A_ik = A_persistent[i_idx, k_idx]
                if (i_idx == k_idx) { # Diagonal element J_ii
                  Jacobian[i_idx, i_idx] = (r_persistent[i_idx] + sum(A_persistent[i_idx,] * N_persistent)) + N_persistent[i_idx] * val_A_ik
                } else { # Off-diagonal element J_ik
                  Jacobian[i_idx, k_idx] = N_persistent[i_idx] * val_A_ik
                }
              }
            }
            
            if (nrow(Jacobian) > 0) {
              eigenvalues_jacobian <- eigen(Jacobian, symmetric = FALSE, only.values = TRUE)$values
              max_Re_val_jacobian <- max(Re(eigenvalues_jacobian), na.rm = TRUE)
            } else {
              max_Re_val_jacobian <- NA # No persistent species, or only 1
            }
          } else { # All species went extinct
            max_Re_val_jacobian <- NA # Or assign a value indicating total collapse
            simulation_ok <- FALSE # Mark as not converged to a community
          }
        } else { # Simulation had issues or all went extinct very early
          max_Re_val_jacobian <- NA # Or some indicator of instability
          simulation_ok <- FALSE
        }
        
        RESULTADOS <- rbind(RESULTADOS, data.frame(max_Re_eigen_jacobian = max_Re_val_jacobian,
                                                   power_H = o_pH,
                                                   power_L = q_pL,
                                                   complementarity = p_comp,
                                                   n_spe_H_val = n_spe_H,
                                                   n_spe_L_val = n_spe_L,
                                                   num_persistent_species = num_persistent,
                                                   simulation_converged = simulation_ok))
      }
    }
  }
}

print("Simulation finished.")
write.csv(RESULTADOS, "RESULTADOS_simulation.csv", row.names = FALSE)
min(RESULTADOS$power_L, na.rm = TRUE)
max(RESULTADOS$power_L, na.rm = TRUE)


# Filter out rows where Jacobian could not be calculated (e.g., due to divergence)
RESULTADOS_filtered <- RESULTADOS[!is.na(RESULTADOS$max_Re_eigen_jacobian), ]
print(paste("Number of successful runs (Jacobian calculated):", nrow(RESULTADOS_filtered)))
if (nrow(RESULTADOS_filtered) == 0) {
  print("Warning: No successful runs to plot or analyze. Check simulation parameters.")
} else {
  print(head(RESULTADOS_filtered))
}
