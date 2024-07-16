# Load necessary libraries
library(ggplot2)
library(dplyr)
library(patchwork)
library(broom)

# Sample data frame creation (Assuming df_Bufo_exp1 is already available)
# df_Bufo_exp1 <- read.csv("your_data.csv")

# Boxplots #
Box_exp1 <- ggplot(df_Bufo_exp1, aes(x = Phase, y = Active_seconds, fill = Treatment)) +
  geom_boxplot(size = 0.75, alpha = 0.5) +  # Adjust alpha for overlay visibility
  facet_grid(rows = ~ Treatment, shrink = T, scales = "free_y", switch = "x", margins = F) +
  scale_y_continuous(limits = c(0, 300)) +
  labs(title = "Active Seconds Across Phases and Treatments",
       x = "",
       y = "Active Seconds [s]") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))  # Rotate x-axis labels to 0 degrees (horizontal)

# Extracting Values from Model predictions #
# Conditional model value predictions
pred_cond_exp1 <- predict(zigam_21_exp1_int_adisPT_zi1_ID, type = "conditional", se.fit = TRUE)
# Zero-inflation probabilities value predictions
pred_zi_exp1 <- predict(zigam_21_exp1_int_adisPT_zi1_ID, type = "zprob", se.fit = TRUE)

#Create a new dataframe for calculations and subsequent plotting
df_pred_exp1 <- data.frame(
  Phase = df_Bufo_exp1$Phase,
  Treatment = df_Bufo_exp1$Treatment,
  Active_seconds = pred_cond_exp1$fit,
  Active_seconds_se = pred_cond_exp1$se.fit,
  ZeroInflation = pred_zi_exp1$fit,
  ZeroInflation_se = pred_zi_exp1$se.fit)

df_pred_exp1$Phase <- factor(df_pred_exp1$Phase, ordered = TRUE)

df_pred_act_exp1 <- df_pred_exp1 %>% 
  select(Phase, Treatment, Active_seconds, Active_seconds_se)

df_pred_zi_exp1 <- df_pred_exp1 %>% 
  select(Phase, Treatment, ZeroInflation, ZeroInflation_se)

# Calculate means and standard errors for each combination of Treatment and Phase
df_calpred_act_exp1 <- df_pred_act_exp1 %>%
  group_by(Treatment, Phase) %>%
  summarise(Active_seconds = mean(Active_seconds, na.rm = TRUE),
            Active_seconds_se = mean(Active_seconds_se, na.rm = TRUE))

df_calpred_zi_exp1 <- df_pred_zi_exp1 %>%
  group_by(Treatment, Phase) %>%
  summarise(mean_zi_prob = mean(ZeroInflation, na.rm = TRUE),
            mean_zi_prob_se = mean(ZeroInflation_se, na.rm = TRUE))

# Combine Box_exp1 and Int_exp1 by overlaying them
Box_Int_combined <- Box_exp1 +
  geom_smooth(method = "lm", se = FALSE, data = df_calpred_act_exp1, aes(x = Phase, y = Active_seconds, group = Treatment, color = Treatment), size = 1) +
  geom_point(data = df_calpred_act_exp1, aes(x = Phase, y = Active_seconds), size = 1) +
  geom_errorbar(data = df_calpred_act_exp1, aes(x = Phase, ymin = Active_seconds - Active_seconds_se, ymax = Active_seconds + Active_seconds_se, color = Treatment), width = 0.1) +
  labs(title = "Active Seconds and Interaction Plot",
       y = "Active Seconds [s]",
       color = "Treatment",
       fill = "Treatment") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

# Plot for Zero Inflation Probabilities
Zi_exp1 <- ggplot(df_calpred_zi_exp1, aes(x = Phase, y = mean_zi_prob, group = Treatment, color = Treatment)) +
  geom_smooth(method = "lm", se = FALSE, size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_zi_prob - mean_zi_prob_se, ymax = mean_zi_prob + mean_zi_prob_se), width = 0.1) +
  facet_grid(rows = ~ Treatment, shrink = T, scales = "free_y", switch = "x", margins = F) +
  scale_y_continuous(limits = c(-0.05, 1)) +
  labs(
    x = "Phase",
    y = "Zero Activity Probability \n [mean +SE]",
    title = "Effects of Phase and Treatment on Zero Activity Probabilities"
  ) +
  theme_bw() +
  theme(legend.position = "none")

# Extract coefficients and confidence intervals from the model
model_coefs <- tidy(zigam_21_exp1_int_adisPT_zi1_ID, conf.int = TRUE)

# Plot for Model Coefficients
Coef_plot <- ggplot(model_coefs, aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_pointrange(color = Treatment) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Term",
    y = "Estimate",
    title = "Model Coefficients"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine plots using patchwork
combined_plot <- (Box_Int_combined / Zi_exp1) | Coef_plot + plot_layout(ncol = 2, widths = c(2, 1))

# Display combined plot
print(combined_plot)
