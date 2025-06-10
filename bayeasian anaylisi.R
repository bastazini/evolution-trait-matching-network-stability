# ---------------------------------------------------
# 1. LOAD REQUIRED PACKAGES
# ---------------------------------------------------
library(rstanarm)      # Bayesian mixed models
library(tidyverse)     # ggplot2, dplyr, tidyr, etc.
library(tidybayes)     # Posterior distribution visualization
library(bayesplot)     # MCMC diagnostics and PPC
library(ggeffects)     # Marginal effects (predicted values)
library(forcats)
library(bayesplot)    #For mcmc_rhat_hist

library(brms)


# ---------------------------------------------------
# 2. FIT THE MODEL 
# ---------------------------------------------------
model_stan1 <- stan_lmer(
  max_Re_eigen_jacobian ~ alpha_H * alpha_L * trait_matching+
  (1 | alpha_H) + (1 | alpha_L),
  data = RESULTADOS,
  chains = 4, cores = 4, iter = 2000
)

model_stan <- stan_glm(
   max_Re_eigen_jacobian ~ alpha_H * alpha_L * trait_matching ,
     #(1 | power_H) + (1 | power_L),
   data = RESULTADOS,
   chains = 4, cores = 4, iter = 2000
 )

 
 #calculate the posterior summary and Rhat values
 posterior_summary <- summary(model_stan, probs = c(0.025, 0.5, 0.975))
 print(posterior_summary) 
 rhat_values <- posterior_summary[, "Rhat"]
 print(rhat_values)
 #Plot Rhat values
 mcmc_rhat_hist(rhat_values) 
# ---------------------------------------------------
# 3. POSTERIOR DISTRIBUTIONS OF FIXED EFFECTS (excluding intercept)
# ---------------------------------------------------
fixed_effects <- names(fixef(model_stan))
fixed_effects <- fixed_effects[fixed_effects != "(Intercept)"]


model_stan %>%
  gather_draws(!!!rlang::syms(fixed_effects)) %>%
  ggplot(aes(x = .value, y = reorder(.variable, .value))) +
  stat_halfeye(.width = c(0.66, 0.95), fill = "steelblue", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
  labs(
    title = "Posterior Distributions of Fixed Effects (Intercept Removed)",
    subtitle = "Shaded areas show 66% and 95% credible intervals\nRed dashed line = null effect (zero)",
    x = "Posterior estimate",
    y = "Parameter"
  ) +
  theme_minimal(base_size = 14)

# ---------------------------------------------------
# 4. MARGINAL EFFECTS WITH INTERACTIONS
# ---------------------------------------------------
marginal_effects <- ggpredict(model_stan, terms = c("power_H", "power_L", "complementarity"))

plot(marginal_effects) +
  labs(
    title = "Predicted Marginal Effects with Credible Intervals",
    x = "Consumer Phylogenetic Signal",
    y = expression("Predicted " * Re(lambda) * "")
  ) +
  theme_minimal(base_size = 14)

# ---------------------------------------------------
# 5. POSTERIOR PREDICTIVE CHECK (PPC)
# ---------------------------------------------------
pp_check(model_stan, plotfun = "dens_overlay", nreps = 50) +
  ggtitle("Posterior Predictive Check: Density Overlay") +
  theme_minimal(base_size = 14)

# ---------------------------------------------------
# 6. MCMC TRACEPLOTS
# ---------------------------------------------------
mcmc_trace(
  as.array(model_stan),
  pars = fixed_effects
) +
  ggtitle("MCMC Traceplots for Fixed Effects") +
  theme_minimal(base_size = 14)

# ---------------------------------------------------
# 7. SAVE PLOTS 
# ---------------------------------------------------
ggsave("fixed_effects_posterior.png", width = 10, height = 6)
ggsave("marginal_effects.png", width = 8, height = 6)
ggsave("ppc_density.png", width = 8, height = 5)

 
 
 # ---------------------------------------------------
 # Plot Fixed effect  with the correct names
 # ---------------------------------------------------
 
 model_stan %>%
   gather_draws(!!!rlang::syms(fixed_effects)) %>%
   mutate(
     .variable = fct_recode(.variable,
                            "Trait matching" = "trait_matching",
                            "Phylo signal (C)" = "alpha_H",
                            "Phylo signal (R)" = "alpha_L",
                            "Phylo signal (C):Phylo signal (R)" = "alpha_H:alpha_L",
                            "Phylo signal (C):Trait matching" = "alpha_H:trait_matching",
                            "Phylo signal (R):Trait matching" = "alpha_L:trait_matching",
                            "Phylo signal (C):Phylo signal (R):Trait matching" = "alpha_H:alpha_L:trait_matching"
     )
   ) %>%
   ggplot(aes(x = .value, y = reorder(.variable, .value))) +
   stat_halfeye(.width = c(0.66, 0.95), fill = "steelblue", color = "black") +
   geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
   labs(
     title = "Posterior Distributions",
     subtitle = "",
     x = "Posterior estimate",
     y = "Parameter"
   ) +
   theme_minimal(base_size = 14)
 ggsave("posterior_plot.pdf", width = 10, height = 6)
 
 # ---------------------------------------------------
 # calculate the effect of terms
 # ---------------------------------------------------
 
 posterior <- as.data.frame(model_stan)
 mean(posterior$`trait_matching` >0)

 ###function to automatize extraction fo probs
 get_probabilities <- function(model) {
   # Extract posterior draws as a data frame
   posterior <- as.data.frame(model)
   
   # Remove intercept if present (optional)
   posterior <- posterior[, !grepl("(Intercept)", colnames(posterior))]
   
   # Calculate probabilities
   probs <- sapply(posterior, function(x) {
     c(
       prob_less_than_zero = mean(x < 0),
       prob_greater_than_zero = mean(x > 0)
     )
   })
   
   # Transpose and convert to data frame
   probs_df <- as.data.frame(t(probs))
   
   # Round for readability
   probs_df <- round(probs_df, 4)
   
   # Return with term names as a column
   probs_df$term <- rownames(probs_df)
   rownames(probs_df) <- NULL
   
   return(probs_df[, c("term", "prob_less_than_zero", "prob_greater_than_zero")])
 }
 
 probabilities <- get_probabilities(model_stan)
 print(probabilities)
 

# ---------------------------------------------------
# model prediction OU
# ---------------------------------------------------
 library(tidybayes)
 
 # Create new data grid to predict
 new_data <- expand.grid(
   alpha_H = c(0.01, 0.1, 0.5, 1),
   alpha_L = c(0.01, 0.1, 0.5, 1),
   trait_matching = seq(0.05, 0.95, length.out = 20)
 )
 
 # Get posterior predictions
 preds <- add_epred_draws(new_data, model_stan)
 
library(tidybayes)

 library(ggplot2)
 library(dplyr)
 
 # 1. Create quartiles for trait matching
 RESULTADOS <- RESULTADOS %>%
   mutate(trait_quartile = ntile(trait_matching, 4)) %>%
   mutate(trait_quartile = factor(trait_quartile, 
                                  labels = c("Trait matching (Q1)", "Trait matching (Q2)", "Trait matching (Q3)", "Trait matching (Q4)")))
 
 # 2. Aggregate results (e.g. mean eigenvalue for each alpha_H/L pair)
 agg_data <- RESULTADOS %>%
   group_by(alpha_H, alpha_L, trait_quartile) %>%
   summarize(mean_eigen = mean(max_Re_eigen_jacobian, na.rm = TRUE), .groups = "drop")
 
 # 3. Plot heatmaps faceted by trait_matching quartiles
 ggplot(agg_data, aes(x = factor(alpha_H), y = factor(alpha_L), fill = mean_eigen)) +
   geom_tile(color = "white") +
   scale_fill_viridis_c(option = "plasma", direction = -1, name = "Mean max Re") +
   facet_wrap(~ trait_quartile, ncol = 2) +
   labs(
     title = "Effect of Trait Evolution Constraints on Stability",
     x = expression(alpha[H]),
     y = expression(alpha[L])
   ) +
   theme_minimal(base_size = 14) +
   theme(
     strip.text = element_text(face = "bold"),
     panel.grid = element_blank(),
     axis.text = element_text(size = 12),
     axis.title = element_text(size = 14),
     legend.position = "bottom"
   )
 
 #### Code for random draws in OU
 # Load libraries
 library(dplyr)
 library(ggplot2)
 library(tidybayes)  # Not used here, but part of original script
 library(viridis)
 
 # Assume RESULTADOS_filtered is already created:
 # If not, filter it here:
 RESULTADOS_filtered <- RESULTADOS[!is.na(RESULTADOS$max_Re_eigen_jacobian), ]
 
 # STEP 1: Bin alpha_H and alpha_L into categories
 RESULTADOS_binned <- RESULTADOS_filtered %>%
   mutate(
     alpha_H_bin = cut(alpha_H, breaks = c(0, 0.5, 1, 2, 3, 5),
                       labels = c("[0–0.5]", "(0.5–1]", "(1–2]", "(2–3]", "(3–5]"),
                       include.lowest = TRUE),
     alpha_L_bin = cut(alpha_L, breaks = c(0, 0.5, 1, 2, 3, 5),
                       labels = c("[0–0.5]", "(0.5–1]", "(1–2]", "(2–3]", "(3–5]"),
                       include.lowest = TRUE),
     trait_quartile = ntile(trait_matching, 4),
     trait_quartile = factor(trait_quartile,
                             labels = c("Trait matching (Q1)", "Trait matching (Q2)",
                                        "Trait matching (Q3)", "Trait matching (Q4)"))
   )
 
 
 # STEP 2: Aggregate data by alpha_H_bin, alpha_L_bin, and trait quartile
 agg_data <- RESULTADOS_binned %>%
   group_by(alpha_H_bin, alpha_L_bin, trait_quartile) %>%
   summarize(mean_eigen = mean(max_Re_eigen_jacobian, na.rm = TRUE), .groups = "drop")
 
 # STEP 3: Plot heatmap faceted by trait quartile
 ggplot(agg_data, aes(x = alpha_H_bin, y = alpha_L_bin, fill = mean_eigen)) +
   geom_tile(color = "white") +
   scale_fill_viridis_c(option = "plasma", direction = -1, name = "Mean max Re") +
   facet_wrap(~ trait_quartile, ncol = 2) +
   labs(
     title = "Effect of Trait Evolution Constraints on Stability",
     x = expression(alpha[C]),
     y = expression(alpha[R])
   ) +
   theme_minimal(base_size = 14) +
   theme(
     strip.text = element_text(face = "bold"),
     panel.grid = element_blank(),
     axis.text = element_text(size = 12),
     axis.title = element_text(size = 14),
     legend.position = "bottom"
   )
 
 #Ploting distributions as alphas
 library(ggplot2)
 library(tidybayes)
 library(forcats)
 library(dplyr)
 
 # Create named vector with plotmath-style labels
 label_map <- c(
   "trait_matching" = "italic('Trait matching')",
   "alpha_H" = "alpha[C]",
   "alpha_L" = "alpha[R]",
   "alpha_H:alpha_L" = "alpha[C] * ':' * alpha[R]",
   "alpha_H:trait_matching" = "alpha[C] * ':' * italic('Trait matching')",
   "alpha_L:trait_matching" = "alpha[R] * ':' * italic('Trait matching')",
   "alpha_H:alpha_L:trait_matching" = "alpha[C] * ':' * alpha[R] * ':' * italic('Trait matching')"
 )
 
 # Plot
 model_stan %>%
   gather_draws(!!!rlang::syms(fixed_effects)) %>%
   mutate(.variable = factor(.variable, levels = names(label_map))) %>%
   ggplot(aes(x = .value, y = .variable)) +
   stat_halfeye(.width = c(0.66, 0.95), fill = "steelblue", color = "black") +
   geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
   scale_y_discrete(labels = label_map, guide = guide_axis(check.overlap = TRUE)) +
   labs(
     title = "Posterior Distributions",
     x = "Posterior estimate",
     y = "Parameter"
   ) +
   theme_minimal(base_size = 14) +
   theme(axis.text.y = element_text(size = 12)) +
   scale_y_discrete(labels = function(x) parse(text = label_map[x]))
 
 # Save to PDF
 ggsave("posterior_plot.pdf", width = 10, height = 6)
 
 
 ##########Plotting predicte dvalues for eigen values
 library(tidybayes)
 
 # For overall predicted distribution
 posterior_predict(model_stan) %>%
   as.data.frame() %>%
   pivot_longer(cols = everything()) %>%
   ggplot(aes(x = value)) +
   stat_halfeye(fill = "steelblue", alpha = 0.6) +
   labs(
     title = "Posterior Predictive Distribution",
     x = "Predicted Re(λ)"
   ) +
   theme_minimal()
 
 
 
 