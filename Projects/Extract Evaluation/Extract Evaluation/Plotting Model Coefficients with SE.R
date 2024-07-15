# Load required packages
library(broom.mixed)
library(dplyr)

# Tidy the model results for both the conditional and zero-inflation parts
tidy_conditional <- tidy(zigam_21_exp1_int_adisPT_ziPT_ID, component = "cond")
tidy_zi <- tidy(zigam_21_exp1_int_adisPT_ziPT_ID, component = "zi")

# Add a column to indicate the model type
tidy_conditional <- tidy_conditional %>%
  mutate(model = "Conditional")

tidy_zi <- tidy_zi %>%
  mutate(model = "Zero Inflation")

# Combine the results into one data frame
tidy_model <- bind_rows(tidy_conditional, tidy_zi)




# Extract summary information
model_summary <- summary(zigam_21_exp1_int_adisPT_ziPT_ID)

# Extract coefficient names in the order they appear in the summary
cond_terms <- rownames(model_summary$coefficients$cond)
zi_terms <- rownames(model_summary$coefficients$zi)

# Create a combined order vector
ordered_terms <- c(cond_terms, zi_terms)

# View the order
print(ordered_terms)





# Filter the fixed effects (excluding intercept if not needed)
fixed_effects <- tidy_model %>%
  filter(effect == "fixed")

# Factorize the term column to ensure the order matches the summary order
fixed_effects <- fixed_effects %>%
  mutate(term = factor(term, levels = ordered_terms))

# View fixed effects
print(fixed_effects)

# Load required packages
library(ggplot2)

# Create the plot
ggplot(fixed_effects, aes(x = term, y = estimate, ymin = estimate - std.error, ymax = estimate + std.error, color = model)) +
  geom_pointrange(position = position_dodge(width = 0.5)) +
  coord_flip() +
  theme_update() +
  labs(
    x = "Model Coefficients",
    y = "Estimate",
    title = "Model Coefficients with Standard Errors",
    color = "Model"
  ) +
  theme(legend.position = "bottom")

