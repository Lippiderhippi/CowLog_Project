# New data for prediction
new_data <- data.frame(Phase = factor("Phase"),
                       Treatment = factor("Treatment"),
                       Individual_Total = factor("Individual_Total"))

# Number of simulations
n_sim <- 300  # Adjust as needed

# Simulate predictions
preds <- simulate(zigam_21_exp1_int_adisPT_ziPT_ID, 
                  newdata = new_data, 
                  nsim = n_sim)
view(preds)

# Extract mean and quantiles for prediction intervals
mean_pred <- apply(preds, 2, mean)
lower_quantile <- apply(preds, 2, function(x) quantile(x, 0.025))
upper_quantile <- apply(preds, 2, function(x) quantile(x, 0.975))

# Combine into a data frame
Pred_df_zi_exp1 <- data.frame(
  Phase = "Phase",  # Update with actual values
  Treatment = "Treatment",  # Update with actual values
  mean_zi_prob = mean_pred,
  se_zi_prob = (upper_quantile - lower_quantile) / 4  # Adjust as needed for your interval
)

# Plot with ggplot
Zi_exp1 <- ggplot(Pred_df_zi_exp1, aes(x = Phase, y = mean_zi_prob, group = Treatment, color = Treatment)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_zi_prob - se_zi_prob, ymax = mean_zi_prob + se_zi_prob), width = 0.1) +
  facet_grid(~ Treatment) +
  theme_bw() +
  labs(
    x = "Phase",
    y = "Zero Inflation Probabilities",
    title = "Effects of Phase and Treatment on Zero Inflation Probabilities"
  ) +
  theme(legend.position = "bottom")

# Display the plot
print(Zi_exp1)
















library(glmmTMB)
library(tidyverse)

# Assuming zigam_21_exp1_int_adisPT_ziPT_ID is your fitted glmmTMB model

# Number of simulations
n_sim <- 300  # Adjust as needed

# New data for prediction
new_data_conditional <- expand.grid(
  Phase = levels(df_Rana_exp1$Phase),
  Treatment = levels(df_Rana_exp1$Treatment),
  Individual_Total = unique(df_Rana_exp1$Individual_Total)
)

# Predict mean Active_seconds and confidence intervals from conditional model
pred_conditional <- predict(zigam_21_exp1_int_adisPT_ziPT_ID, type = "conditional", se.fit = TRUE)

# Extracting mean predictions and confidence intervals
pred_conditional <- data.frame(
  Phase = new_data_conditional$Phase,
  Treatment = new_data_conditional$Treatment,
  mean_Active_seconds = pred_conditional$fit,
  se_Active_seconds = pred_conditional$se.fit
)

# Calculate confidence intervals (CI)
pred_conditional <- pred_conditional %>%
  mutate(
    lower_CI_Active_seconds = mean_Active_seconds - 1.96 * se_Active_seconds,  # 95% CI
    upper_CI_Active_seconds = mean_Active_seconds + 1.96 * se_Active_seconds   # 95% CI
  )

# Plot mean Active_seconds with confidence intervals
ggplot(pred_conditional, aes(x = Phase, y = mean_Active_seconds,group = Treatment, color = Treatment)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower_CI_Active_seconds, ymax = upper_CI_Active_seconds), width = 0.1) +
  facet_grid(~ Treatment) +
  labs(
    x = "Phase",
    y = "Mean Active Seconds",
    title = "Conditional Model: Mean Active Seconds with 95% CI"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")



# Simulate predictions for probabilities of zeros
#pred_zero_inflated <- predict(zigam_21_exp1_int_adisPT_ziPT_ID, newdata = new_data_conditional, type = "zero_inflated", nsim = 300)
# Zero-inflation probabilities predictions
pred_zi <- predict(zigam_21_exp1_int_adisPT_ziPT_ID, type = "zprob", se.fit = TRUE)
# Calculate mean probabilities of zeros
mean_probs_zeros <- apply(pred_zi, 2, mean)


# Calculate 95% confidence intervals for probabilities of zeros
lower_CI_probs_zeros <- apply(pred_zero_inflated, 2, function(x) quantile(x, 0.025))
upper_CI_probs_zeros <- apply(pred_zero_inflated, 2, function(x) quantile(x, 0.975))

# Create data frame for plotting
Pred_df_zero_inflated <- data.frame(
  Phase = new_data_conditional$Phase,
  Treatment = new_data_conditional$Treatment,
  mean_prob_zeros = mean_probs_zeros,
  lower_CI_prob_zeros = lower_CI_probs_zeros,
  upper_CI_prob_zeros = upper_CI_probs_zeros
)

# Plot mean probabilities of zeros with confidence intervals
ggplot(Pred_df_zero_inflated, aes(x = Phase, y = mean_prob_zeros, color = Treatment)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower_CI_prob_zeros, ymax = upper_CI_prob_zeros), width = 0.1) +
  facet_grid(~ Treatment) +
  labs(
    x = "Phase",
    y = "Mean Probability of Zeros",
    title = "Zero-Inflated Model: Mean Probability of Zeros with 95% CI"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")














