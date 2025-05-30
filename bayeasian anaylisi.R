# ---------------------------------------------------
# 1. LOAD REQUIRED PACKAGES
# ---------------------------------------------------
library(rstanarm)      # Bayesian mixed models
library(tidyverse)     # ggplot2, dplyr, tidyr, etc.
library(tidybayes)     # Posterior distribution visualization
library(bayesplot)     # MCMC diagnostics and PPC
library(ggeffects)     # Marginal effects (predicted values)

# ---------------------------------------------------
# 2. FIT THE MODEL 
# ---------------------------------------------------
 model_stan <- stan_lmer(
  real_eigenvalue ~ power_H * power_L * complementarity +
     (1 | power_H) + (1 | power_L),
   data = RESULTADOS,
   chains = 4, cores = 4, iter = 2000
 )

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
    x = "power_H",
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
# 7. OPTIONAL: SAVE PLOTS TO FILES
# ---------------------------------------------------
ggsave("fixed_effects_posterior.png", width = 10, height = 6)
ggsave("marginal_effects.png", width = 8, height = 6)
 ggsave("ppc_density.png", width = 8, height = 5)

# ---------------------------------------------------
# END OF SCRIPT
# ---------------------------------------------------
 
