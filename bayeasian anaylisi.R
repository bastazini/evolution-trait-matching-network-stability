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


# ---------------------------------------------------
# 2. FIT THE MODEL 
# ---------------------------------------------------
 model_stan <- stan_lmer(
  real_eigenvalue ~ power_H * power_L * complementarity +
     (1 | power_H) + (1 | power_L),
   data = RESULTADOS,
   chains = 4, cores = 4, iter = 2000
 )

 
 #calculate the posterior summary and Rhat values
 posterior_summary <- summary(model_stan, probs = c(0.025, 0.5, 0.975))
 print(posterior_summary) 
 rhat_values <- posterior_summary[, "Rhat"]
 print(rhat_values)
 #Plot Rhat values
 mcmc_rhat_hist(rhat(model_stan))
 
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
                            "Trait matching" = "complementarity",
                            "Phylo signal (C)" = "power_H",
                            "Phylo signal (R)" = "power_L",
                            "Phylo signal (C):Phylo signal (R)" = "power_H:power_L",
                            "Phylo signal (C):Trait matching" = "power_H:complementarity",
                            "Phylo signal (R):Trait matching" = "power_L:complementarity",
                            "Phylo signal (C):Phylo signal (R):Trait matching" = "power_H:power_L:complementarity"
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
 # calculate the effect of the interaction
 # ---------------------------------------------------
 
 posterior <- as.data.frame(model_stan)
 mean(posterior$`power_H:power_L` > 0)
 

# ---------------------------------------------------
# CONTRASTS
# ---------------------------------------------------
 
