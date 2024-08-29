### Experiment 1 - Boxplot & EMMs + Zero-Inflation ---------------------------

# Calculate Estimated Marginal Means (EMMs) on Log Scale and Original Scale
emmeans_result_exp1_Rana <- emmeans(zigam_21_exp1_int_adisPT_ziPT_ID_Rana, ~ Phase * Treatment)
original_scale_result_exp1_Rana <- summary(emmeans_result_exp1_Rana, type = "response")

# Convert to data frame and rename columns
df_emm_exp1_Rana <- as.data.frame(original_scale_result_exp1_Rana)
colnames(df_emm_exp1_Rana) <- c("Phase", "Treatment", "Active_seconds", "SE", "df", "Lower_CI", "Upper_CI")
df_emm_exp1_Rana$Phase <- factor(df_emm_exp1_Rana$Phase, ordered = TRUE)

# Boxplot with EMMs
Box_Emm_exp1_Rana <- ggplot(df_Rana_exp1, aes(x = Phase, y = Active_seconds, fill = Treatment)) +
  geom_boxplot(size = 0.1, alpha = 0.3, outlier.shape = 16, outlier.color = "black", outlier.size = 1) +
  geom_point(data = df_emm_exp1_Rana, aes(x = Phase, y = Active_seconds, color = Treatment), size = 1, shape = 21) +
  geom_smooth(linewidth = 0.5, method = "lm", se = FALSE, aes(group = Treatment, color = Treatment)) +
  geom_errorbar(data = df_emm_exp1_Rana, aes(x = Phase, ymin = Lower_CI, ymax = Upper_CI, color = Treatment), width = 0.2) +
  facet_grid(rows = ~ Treatment, scales = "free_y", switch = "x", margins = F) +
  scale_y_continuous(limits = c(0, 200)) +
  labs(y = "Active Seconds \n [s]") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 0.5))

# Zero-Inflation Probabilities Predictions
pred_zi_exp1_Rana <- predict(zigam_21_exp1_int_adisPT_ziPT_ID_Rana, type = "zprob", se.fit = TRUE)
df_pred_zi_exp1_Rana <- data.frame(Phase = factor(df_Rana_exp1$Phase, ordered = TRUE),
                                   Treatment = df_Rana_exp1$Treatment,
                                   ZeroInflation = pred_zi_exp1_Rana$fit,
                                   ZeroInflation_se = pred_zi_exp1_Rana$se.fit)

# Calculate means and standard errors for Zero-Inflation Probabilities
df_calpred_zi_exp1_Rana <- df_pred_zi_exp1_Rana %>%
  group_by(Treatment, Phase) %>%
  summarise(mean_zi_prob = mean(ZeroInflation, na.rm = TRUE),
            mean_zi_prob_se = mean(ZeroInflation_se, na.rm = TRUE))

# Plot for Zero-Inflation Probabilities
Zi_exp1_Rana <- ggplot(df_calpred_zi_exp1_Rana, aes(x = Phase, y = mean_zi_prob, group = Treatment, color = Treatment)) +
  geom_smooth(linewidth = 0.5, method = "lm", se = FALSE) +
  geom_point(size = 1) +
  geom_errorbar(linewidth = 0.5, width = 0.2, aes(ymin = mean_zi_prob - mean_zi_prob_se, ymax = mean_zi_prob + mean_zi_prob_se)) +
  facet_grid(rows = ~ Treatment, scales = "free_y", switch = "x", margins = F) +
  scale_y_continuous(limits = c(-0.05, 1)) +
  labs(x = "Phase", y = "Probability of Freezing \n [mean + SE]") +
  theme_bw() +
  theme(legend.position = "none")

# Combine and display plots
combined_plot_exp1_Rana <- Box_Emm_exp1_Rana / Zi_exp1_Rana + plot_layout(ncol = 1)
print(combined_plot_exp1_Rana)

### Experiment 2 - Boxplot & EMMs + Zero-Inflation ---------------------------

# Replicate the steps for Experiment 2 using the same structure as Experiment 1
emmeans_result_exp2_Rana <- emmeans(zigam_21_exp2_int_adisPT_ziPT_ID_Rana, ~ Phase * Treatment)
original_scale_result_exp2_Rana <- summary(emmeans_result_exp2_Rana, type = "response")
df_emm_exp2_Rana <- as.data.frame(original_scale_result_exp2_Rana)
colnames(df_emm_exp2_Rana) <- c("Phase", "Treatment", "Active_seconds", "SE", "df", "Lower_CI", "Upper_CI")
df_emm_exp2_Rana$Phase <- factor(df_emm_exp2_Rana$Phase, ordered = TRUE)

Box_Emm_exp2_Rana <- ggplot(df_Rana_exp2, aes(x = Phase, y = Active_seconds, fill = Treatment)) +
  geom_boxplot(size = 0.1, alpha = 0.3, outlier.shape = 16, outlier.color = "black", outlier.size = 1) +
  geom_point(data = df_emm_exp2_Rana, aes(x = Phase, y = Active_seconds, color = Treatment), size = 1, shape = 21) +
  geom_smooth(linewidth = 0.5, method = "lm", se = FALSE, aes(group = Treatment, color = Treatment)) +
  geom_errorbar(data = df_emm_exp2_Rana, aes(x = Phase, ymin = Lower_CI, ymax = Upper_CI, color = Treatment), width = 0.2) +
  facet_grid(rows = ~ Treatment, scales = "free_y", switch = "x", margins = F) +
  scale_y_continuous(limits = c(0, 200)) +
  labs(y = "Active Seconds \n [s]") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 0.5))

pred_zi_exp2_Rana <- predict(zigam_21_exp2_int_adisPT_ziPT_ID_Rana, type = "zprob", se.fit = TRUE)
df_pred_zi_exp2_Rana <- data.frame(Phase = factor(df_Rana_exp2$Phase, ordered = TRUE),
                                   Treatment = df_Rana_exp2$Treatment,
                                   ZeroInflation = pred_zi_exp2_Rana$fit,
                                   ZeroInflation_se = pred_zi_exp2_Rana$se.fit)

df_calpred_zi_exp2_Rana <- df_pred_zi_exp2_Rana %>%
  group_by(Treatment, Phase) %>%
  summarise(mean_zi_prob = mean(ZeroInflation, na.rm = TRUE),
            mean_zi_prob_se = mean(ZeroInflation_se, na.rm = TRUE))

Zi_exp2_Rana <- ggplot(df_calpred_zi_exp2_Rana, aes(x = Phase, y = mean_zi_prob, group = Treatment, color = Treatment)) +
  geom_smooth(linewidth = 0.5, method = "lm", se = FALSE) +
  geom_point(size = 1) +
  geom_errorbar(linewidth = 0.5, width = 0.2, aes(ymin = mean_zi_prob - mean_zi_prob_se, ymax = mean_zi_prob + mean_zi_prob_se)) +
  facet_grid(rows = ~ Treatment, scales = "free_y", switch = "x", margins = F) +
  scale_y_continuous(limits = c(-0.05, 1)) +
  labs(x = "Phase", y = "Probability of Freezing \n [mean + SE]") +
  theme_bw() +
  theme(legend.position = "none")

combined_plot_exp2_Rana <- Box_Emm_exp2_Rana / Zi_exp2_Rana + plot_layout(ncol = 1)
print(combined_plot_exp2_Rana)


















#############################
##### RANA TEMPORARIA #######
#############################

### Experiment 1 - Boxplot & Interaction + Zero-Inflation ---------------------------

# Boxplots
Box_exp1_Rana <- ggplot(df_Rana_exp1, aes(x = Phase, y = Active_seconds, fill = Treatment)) +
  geom_boxplot(size = 0.1, alpha = 0.3) +
  facet_grid(rows = ~ Treatment, scales = "free_y", switch = "x", margins = F) +
  scale_y_continuous(limits = c(0, 200)) +
  labs(y = "Active Seconds \n [s]") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 0.5))

# Extracting Values from Model predictions
# Conditional model value predictions
pred_cond_exp1_Rana <- predict(zigam_21_exp1_int_adisPT_ziPT_ID_Rana, type = "conditional", se.fit = TRUE)
# Zero-inflation probabilities value predictions
pred_zi_exp1_Rana <- predict(zigam_21_exp1_int_adisPT_ziPT_ID_Rana, type = "zprob", se.fit = TRUE)

# Create a new dataframe for calculations and subsequent plotting
df_pred_exp1_Rana <- data.frame(
  Phase = factor(df_Rana_exp1$Phase, ordered = TRUE),
  Treatment = df_Rana_exp1$Treatment,
  Active_seconds = pred_cond_exp1_Rana$fit,
  Active_seconds_se = pred_cond_exp1_Rana$se.fit,
  ZeroInflation = pred_zi_exp1_Rana$fit,
  ZeroInflation_se = pred_zi_exp1_Rana$se.fit
)

# Calculate means and standard errors for each combination of Treatment and Phase
df_calpred_act_exp1_Rana <- df_pred_exp1_Rana %>%
  group_by(Treatment, Phase) %>%
  summarise(Active_seconds = mean(Active_seconds, na.rm = TRUE),
            Active_seconds_se = mean(Active_seconds_se, na.rm = TRUE))

df_calpred_zi_exp1_Rana <- df_pred_exp1_Rana %>%
  group_by(Treatment, Phase) %>%
  summarise(mean_zi_prob = mean(ZeroInflation, na.rm = TRUE),
            mean_zi_prob_se = mean(ZeroInflation_se, na.rm = TRUE))

# Combine Boxplots and Interaction Plot
Box_Int_combined_exp1_Rana <- Box_exp1_Rana +
  geom_smooth(linewidth = 0.5, method = "lm", se = FALSE, data = df_calpred_act_exp1, aes(x = Phase, y = Active_seconds, group = Treatment, color = Treatment)) +
  geom_point(size = 1, data = df_calpred_act_exp1, aes(x = Phase, y = Active_seconds)) +
  geom_errorbar(linewidth = 0.5, width = 0.2, data = df_calpred_act_exp1, aes(x = Phase, ymin = Active_seconds - Active_seconds_se, ymax = Active_seconds + Active_seconds_se, color = Treatment)) +
  labs(y = "Active Seconds \n [s]", color = "Treatment", fill = "Treatment")

# Plot for Zero-Inflation Probabilities
Zi_exp1_Rana <- ggplot(df_calpred_zi_exp1_Rana, aes(x = Phase, y = mean_zi_prob, group = Treatment, color = Treatment)) +
  geom_smooth(linewidth = 0.5, method = "lm", se = FALSE) +
  geom_point(size = 1) +
  geom_errorbar(linewidth = 0.5, width = 0.2, aes(ymin = mean_zi_prob - mean_zi_prob_se, ymax = mean_zi_prob + mean_zi_prob_se)) +
  facet_grid(rows = ~ Treatment, scales = "free_y", switch = "x", margins = F) +
  scale_y_continuous(limits = c(-0.05, 1)) +
  labs(x = "Phase", y = "Probability of Freezing \n [mean + SE]") +
  theme_bw() +
  theme(legend.position = "none")

# Combine plots using patchwork
combined_plot_exp1_Rana <- Box_Int_combined_exp1_Rana / Zi_exp1_Rana + plot_layout(ncol = 1, heights = c(1, 1))

# Display combined plot
print(combined_plot_exp1_Rana)

### Experiment 2 - Boxplot & Interaction + Zero-Inflation ---------------------------

# Repeat the same process for Experiment 2
Box_exp2_Rana <- ggplot(df_Rana_exp2, aes(x = Phase, y = Active_seconds, fill = Treatment)) +
  geom_boxplot(size = 0.1, alpha = 0.3) +
  facet_grid(rows = ~ Treatment, scales = "free_y", switch = "x", margins = F) +
  scale_y_continuous(limits = c(0, 200)) +
  labs(y = "Active Seconds \n [s]") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 0.5))

pred_cond_exp2_Rana <- predict(zigam_21_exp2_int_adisPT_ziPT_ID_Rana, type = "conditional", se.fit = TRUE)
pred_zi_exp2_Rana <- predict(zigam_21_exp2_int_adisPT_ziPT_ID_Rana, type = "zprob", se.fit = TRUE)

df_pred_exp2_Rana <- data.frame(
  Phase = factor(df_pred_exp2_Rana$Phase, ordered = TRUE),
  Treatment = df_pred_exp2_Rana$Treatment,
  Active_seconds = pred_cond_exp2_Rana$fit,
  Active_seconds_se = pred_cond_exp2_Rana$se.fit,
  ZeroInflation = pred_zi_exp2_Rana$fit,
  ZeroInflation_se = pred_zi_exp2_Rana$se.fit
)

df_calpred_act_exp2_Rana <- df_pred_exp2_Rana %>%
  group_by(Treatment, Phase) %>%
  summarise(Active_seconds = mean(Active_seconds, na.rm = TRUE),
            Active_seconds_se = mean(Active_seconds_se, na.rm = TRUE))

df_calpred_zi_exp2_Rana <- df_pred_exp2 %>%
  group_by(Treatment, Phase) %>%
  summarise(mean_zi_prob = mean(ZeroInflation, na.rm = TRUE),
            mean_zi_prob_se = mean(ZeroInflation_se, na.rm = TRUE))

Box_Int_combined_exp2_Rana <- Box_exp2_Rana +
  geom_smooth(linewidth = 0.5, method = "lm", se = FALSE, data = df_calpred_act_exp2, aes(x = Phase, y = Active_seconds, group = Treatment, color = Treatment)) +
  geom_point(size = 1, data = df_calpred_act_exp2, aes(x = Phase, y = Active_seconds)) +
  geom_errorbar(linewidth = 0.5, width = 0.2, data = df_calpred_act_exp2, aes(x = Phase, ymin = Active_seconds - Active_seconds_se, ymax = Active_seconds + Active_seconds_se, color = Treatment)) +
  labs(y = "Active Seconds \n [s]", color = "Treatment", fill = "Treatment")

Zi_exp2_Rana <- ggplot(df_calpred_zi_exp2_Rana, aes(x = Phase, y = mean_zi_prob, group = Treatment, color = Treatment)) +
  geom_smooth(linewidth = 0.5, method = "lm", se = FALSE) +
  geom_point(size = 1) +
  geom_errorbar(linewidth = 0.5, width = 0.2, aes(ymin = mean_zi_prob - mean_zi_prob_se, ymax = mean_zi_prob + mean_zi_prob_se)) +
  facet_grid(rows = ~ Treatment, scales = "free_y", switch = "x", margins = F) +
  scale_y_continuous(limits = c(-0.05, 1)) +
  labs(x = "Phase", y = "Probability of Freezing \n [mean + SE]") +
  theme_bw() +
  theme(legend.position = "none")

combined_plot_exp2_Rana <- Box_Int_combined_exp2_Rana / Zi_exp2_Rana + plot_layout(ncol = 1, heights = c(1, 1))

print(combined_plot_exp2_Rana)

# Combine two combined plots side by side
final_combined_plot_exp1and2_Rana <- combined_plot_exp1_Rana | combined_plot_exp2_Rana
print(final_combined_plot_Rana)



#######################
##### BUFO BUFO #######
#######################


### Experiment 1 - Boxplot & Interaction + Zero-Inflation for Bufo ---------------------------

# Boxplots
Box_exp1_Bufo <- ggplot(df_Bufo_exp1, aes(x = Phase, y = Active_seconds, fill = Treatment)) +
  geom_boxplot(size = 0.1, alpha = 0.3) +
  facet_grid(rows = ~ Treatment, scales = "free_y", switch = "x", margins = F) +
  scale_y_continuous(limits = c(0, 300)) +
  labs(y = "Active Seconds \n [s]") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 0.5))

# Extracting Values from Model predictions
# Conditional model value predictions
pred_cond_exp1_Bufo <- predict(zigam_21_exp1_int_adisPT_zi1_ID_Bufo, type = "conditional", se.fit = TRUE)
# Zero-inflation probabilities value predictions
pred_zi_exp1_Bufo <- predict(zigam_21_exp1_int_adisPT_zi1_ID_Bufo, type = "zprob", se.fit = TRUE)

# Create a new dataframe for calculations and subsequent plotting
df_pred_exp1_Bufo <- data.frame(
  Phase = factor(df_Bufo_exp1$Phase, ordered = TRUE),
  Treatment = df_Bufo_exp1$Treatment,
  Active_seconds = pred_cond_exp1_Bufo$fit,
  Active_seconds_se = pred_cond_exp1_Bufo$se.fit,
  ZeroInflation = pred_zi_exp1_Bufo$fit,
  ZeroInflation_se = pred_zi_exp1_Bufo$se.fit
)

# Calculate means and standard errors for each combination of Treatment and Phase
df_calpred_act_exp1_Bufo <- df_pred_exp1_Bufo %>%
  group_by(Treatment, Phase) %>%
  summarise(Active_seconds = mean(Active_seconds, na.rm = TRUE),
            Active_seconds_se = mean(Active_seconds_se, na.rm = TRUE))

df_calpred_zi_exp1_Bufo <- df_pred_exp1_Bufo %>%
  group_by(Treatment, Phase) %>%
  summarise(mean_zi_prob = mean(ZeroInflation, na.rm = TRUE),
            mean_zi_prob_se = mean(ZeroInflation_se, na.rm = TRUE))

# Combine Boxplots and Interaction Plot
Box_Int_combined_Bufo_exp1 <- Box_exp1_Bufo +
  geom_smooth(linewidth = 0.5, method = "lm", se = FALSE, data = df_calpred_act_exp1_Bufo, aes(x = Phase, y = Active_seconds, group = Treatment, color = Treatment)) +
  geom_point(size = 1, data = df_calpred_act_exp1_Bufo, aes(x = Phase, y = Active_seconds)) +
  geom_errorbar(linewidth = 0.5, width = 0.2, data = df_calpred_act_exp1_Bufo, aes(x = Phase, ymin = Active_seconds - Active_seconds_se, ymax = Active_seconds + Active_seconds_se, color = Treatment)) +
  labs(y = "Active Seconds \n [s]", color = "Treatment", fill = "Treatment")

# Plot for Zero-Inflation Probabilities
Zi_exp1_Bufo <- ggplot(df_calpred_zi_exp1_Bufo, aes(x = Phase, y = mean_zi_prob, group = Treatment, color = Treatment)) +
  geom_smooth(linewidth = 0.5, method = "lm", se = FALSE) +
  geom_point(size = 1) +
  geom_errorbar(linewidth = 0.5, width = 0.2, aes(ymin = mean_zi_prob - mean_zi_prob_se, ymax = mean_zi_prob + mean_zi_prob_se)) +
  facet_grid(rows = ~ Treatment, scales = "free_y", switch = "x", margins = F) +
  scale_y_continuous(limits = c(-0.05, 1)) +
  labs(x = "Phase", y = "Probability of Freezing \n [mean + SE]") +
  theme_bw() +
  theme(legend.position = "none")

# Combine plots using patchwork
combined_plot_Bufo_exp1 <- Box_Int_combined_Bufo_exp1 / Zi_exp1_Bufo + plot_layout(ncol = 1, heights = c(1, 1))

# Display combined plot
print(combined_plot_Bufo_exp1)

### Experiment 2 - Boxplot & Interaction + Zero-Inflation for Bufo ---------------------------

# Repeat the same process for Experiment 2
Box_exp2_Bufo <- ggplot(df_Bufo_exp2, aes(x = Phase, y = Active_seconds, fill = Treatment)) +
  geom_boxplot(size = 0.1, alpha = 0.3) +
  facet_grid(rows = ~ Treatment, scales = "free_y", switch = "x", margins = F) +
  scale_y_continuous(limits = c(0, 300)) +
  labs(y = "Active Seconds \n [s]") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 0.5))

pred_cond_exp2_Bufo <- predict(zigam_21_exp2_int_adisPT_zi1_ID_Bufo, type = "conditional", se.fit = TRUE)
pred_zi_exp2_Bufo <- predict(zigam_21_exp2_int_adisPT_zi1_ID_Bufo, type = "zprob", se.fit = TRUE)

df_pred_exp2_Bufo <- data.frame(
  Phase = factor(df_Bufo_exp2$Phase, ordered = TRUE),
  Treatment = df_Bufo_exp2$Treatment,
  Active_seconds = pred_cond_exp2_Bufo$fit,
  Active_seconds_se = pred_cond_exp2_Bufo$se.fit,
  ZeroInflation = pred_zi_exp2_Bufo$fit,
  ZeroInflation_se = pred_zi_exp2_Bufo$se.fit
)

df_calpred_act_exp2_Bufo <- df_pred_exp2_Bufo %>%
  group_by(Treatment, Phase) %>%
  summarise(Active_seconds = mean(Active_seconds, na.rm = TRUE),
            Active_seconds_se = mean(Active_seconds_se, na.rm = TRUE))

df_calpred_zi_exp2_Bufo <- df_pred_exp2_Bufo %>%
  group_by(Treatment, Phase) %>%
  summarise(mean_zi_prob = mean(ZeroInflation, na.rm = TRUE),
            mean_zi_prob_se = mean(ZeroInflation_se, na.rm = TRUE))

Box_Int_combined_Bufo_exp2 <- Box_exp2_Bufo +
  geom_smooth(linewidth = 0.5, method = "lm", se = FALSE, data = df_calpred_act_exp2_Bufo, aes(x = Phase, y = Active_seconds, group = Treatment, color = Treatment)) +
  geom_point(size = 1, data = df_calpred_act_exp2_Bufo, aes(x = Phase, y = Active_seconds)) +
  geom_errorbar(linewidth = 0.5, width = 0.2, data = df_calpred_act_exp2_Bufo, aes(x = Phase, ymin = Active_seconds - Active_seconds_se, ymax = Active_seconds + Active_seconds_se, color = Treatment)) +
  labs(y = "Active Seconds \n [s]", color = "Treatment", fill = "Treatment")

Zi_exp2_Bufo <- ggplot(df_calpred_zi_exp2_Bufo, aes(x = Phase, y = mean_zi_prob, group = Treatment, color = Treatment)) +
  geom_smooth(linewidth = 0.5, method = "lm", se = FALSE) +
  geom_point(size = 1) +
  geom_errorbar(linewidth = 0.5, width = 0.2, aes(ymin = mean_zi_prob - mean_zi_prob_se, ymax = mean_zi_prob + mean_zi_prob_se)) +
  facet_grid(rows = ~ Treatment, scales = "free_y", switch = "x", margins = F) +
  scale_y_continuous(limits = c(-0.05, 1)) +
  labs(x = "Phase", y = "Probability of Freezing \n [mean + SE]") +
  theme_bw() +
  theme(legend.position = "none")

combined_plot_Bufo_exp2 <- Box_Int_combined_Bufo_exp2 / Zi_exp2_Bufo + plot_layout(ncol = 1, heights = c(1, 1))

print(combined_plot_Bufo_exp2)

# Combine two combined plots side by side
final_combined_plot_Bufo <- combined_plot_Bufo_exp1 | combined_plot_Bufo_exp2
print(final_combined_plot_Bufo)




# Combine all plots in a grid with 2 rows and 2 columns
final_combined_plot_BufoRana <- (combined_plot_Bufo_exp1 | combined_plot_Bufo_exp2) / (combined_plot_Rana_exp1 | combined_plot_Rana_exp2)
print(final_combined_plot_BufoRana)
















#############################
##### RANA TEMPORARIA #######
#############################

### Experiment 1 - Boxplot & Interaction + Zero-Inflation ---------------------------

# Boxplots
Box_exp1_Rana <- ggplot(df_Rana_exp1, aes(x = Phase, y = Active_seconds, fill = Treatment)) +
  geom_boxplot(size = 0.1, alpha = 0.3) +
  facet_grid(rows = ~ Treatment, scales = "free_y", switch = "x", margins = F) +
  scale_y_continuous(limits = c(0, 200)) +
  labs(y = "Active Seconds \n [s]") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 0.5))

# Extracting Values from Model predictions
pred_cond_exp1_Rana <- predict(zigam_21_exp1_int_adisPT_ziPT_ID_Rana, type = "conditional", se.fit = TRUE)
pred_zi_exp1_Rana <- predict(zigam_21_exp1_int_adisPT_ziPT_ID_Rana, type = "zprob", se.fit = TRUE)

# Create a new dataframe for calculations and subsequent plotting
df_pred_exp1_Rana <- data.frame(
  Phase = factor(df_Rana_exp1$Phase, ordered = TRUE),
  Treatment = df_Rana_exp1$Treatment,
  Active_seconds = pred_cond_exp1_Rana$fit,
  Active_seconds_se = pred_cond_exp1_Rana$se.fit,
  ZeroInflation = pred_zi_exp1_Rana$fit,
  ZeroInflation_se = pred_zi_exp1_Rana$se.fit
)

df_calpred_act_exp1_Rana <- df_pred_exp1_Rana %>%
  group_by(Treatment, Phase) %>%
  summarise(Active_seconds = mean(Active_seconds, na.rm = TRUE),
            Active_seconds_se = mean(Active_seconds_se, na.rm = TRUE))

df_calpred_zi_exp1_Rana <- df_pred_exp1_Rana %>%
  group_by(Treatment, Phase) %>%
  summarise(mean_zi_prob = mean(ZeroInflation, na.rm = TRUE),
            mean_zi_prob_se = mean(ZeroInflation_se, na.rm = TRUE))

# Combine Boxplots and Interaction Plot
Box_Int_combined_exp1_Rana <- Box_exp1_Rana +
  geom_smooth(linewidth = 0.5, method = "lm", se = FALSE, data = df_calpred_act_exp1_Rana, aes(x = Phase, y = Active_seconds, group = Treatment, color = Treatment)) +
  geom_point(size = 1, data = df_calpred_act_exp1_Rana, aes(x = Phase, y = Active_seconds)) +
  geom_errorbar(linewidth = 0.5, width = 0.2, data = df_calpred_act_exp1_Rana, aes(x = Phase, ymin = Active_seconds - Active_seconds_se, ymax = Active_seconds + Active_seconds_se, color = Treatment)) +
  labs(y = "Active Seconds \n [s]", color = "Treatment", fill = "Treatment")

# Plot for Zero-Inflation Probabilities
Zi_exp1_Rana <- ggplot(df_calpred_zi_exp1_Rana, aes(x = Phase, y = mean_zi_prob, group = Treatment, color = Treatment)) +
  geom_smooth(linewidth = 0.5, method = "lm", se = FALSE) +
  geom_point(size = 1) +
  geom_errorbar(linewidth = 0.5, width = 0.2, aes(ymin = mean_zi_prob - mean_zi_prob_se, ymax = mean_zi_prob + mean_zi_prob_se)) +
  facet_grid(rows = ~ Treatment, scales = "free_y", switch = "x", margins = F) +
  scale_y_continuous(limits = c(-0.05, 1)) +
  labs(x = "Phase", y = "Probability of Freezing \n [mean + SE]") +
  theme_bw() +
  theme(legend.position = "none")

combined_plot_exp1_Rana <- Box_Int_combined_exp1_Rana / Zi_exp1_Rana + plot_layout(ncol = 1, heights = c(1, 1))
print(combined_plot_exp1_Rana)

### Experiment 2 - Boxplot & Interaction + Zero-Inflation ---------------------------

Box_exp2_Rana <- ggplot(df_Rana_exp2, aes(x = Phase, y = Active_seconds, fill = Treatment)) +
  geom_boxplot(size = 0.1, alpha = 0.3) +
  facet_grid(rows = ~ Treatment, scales = "free_y", switch = "x", margins = F) +
  scale_y_continuous(limits = c(0, 200)) +
  labs(y = "Active Seconds \n [s]") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 0.5))

pred_cond_exp2_Rana <- predict(zigam_21_exp2_int_adisPT_ziPT_ID_Rana, type = "conditional", se.fit = TRUE)
pred_zi_exp2_Rana <- predict(zigam_21_exp2_int_adisPT_ziPT_ID_Rana, type = "zprob", se.fit = TRUE)

df_pred_exp2_Rana <- data.frame(
  Phase = factor(df_Rana_exp2$Phase, ordered = TRUE),
  Treatment = df_Rana_exp2$Treatment,
  Active_seconds = pred_cond_exp2_Rana$fit,
  Active_seconds_se = pred_cond_exp2_Rana$se.fit,
  ZeroInflation = pred_zi_exp2_Rana$fit,
  ZeroInflation_se = pred_zi_exp2_Rana$se.fit
)

df_calpred_act_exp2_Rana <- df_pred_exp2_Rana %>%
  group_by(Treatment, Phase) %>%
  summarise(Active_seconds = mean(Active_seconds, na.rm = TRUE),
            Active_seconds_se = mean(Active_seconds_se, na.rm = TRUE))

df_calpred_zi_exp2_Rana <- df_pred_exp2_Rana %>%
  group_by(Treatment, Phase) %>%
  summarise(mean_zi_prob = mean(ZeroInflation, na.rm = TRUE),
            mean_zi_prob_se = mean(ZeroInflation_se, na.rm = TRUE))

Box_Int_combined_exp2_Rana <- Box_exp2_Rana +
  geom_smooth(linewidth = 0.5, method = "lm", se = FALSE, data = df_calpred_act_exp2_Rana, aes(x = Phase, y = Active_seconds, group = Treatment, color = Treatment)) +
  geom_point(size = 1, data = df_calpred_act_exp2_Rana, aes(x = Phase, y = Active_seconds)) +
  geom_errorbar(linewidth = 0.5, width = 0.2, data = df_calpred_act_exp2_Rana, aes(x = Phase, ymin = Active_seconds - Active_seconds_se, ymax = Active_seconds + Active_seconds_se, color = Treatment)) +
  labs(y = "Active Seconds \n [s]", color = "Treatment", fill = "Treatment")

Zi_exp2_Rana <- ggplot(df_calpred_zi_exp2_Rana, aes(x = Phase, y = mean_zi_prob, group = Treatment, color = Treatment)) +
  geom_smooth(linewidth = 0.5, method = "lm", se = FALSE) +
  geom_point(size = 1) +
  geom_errorbar(linewidth = 0.5, width = 0.2, aes(ymin = mean_zi_prob - mean_zi_prob_se, ymax = mean_zi_prob + mean_zi_prob_se)) +
  facet_grid(rows = ~ Treatment, scales = "free_y", switch = "x", margins = F) +
  scale_y_continuous(limits = c(-0.05, 1)) +
  labs(x = "Phase", y = "Probability of Freezing \n [mean + SE]") +
  theme_bw() +
  theme(legend.position = "none")

combined_plot_exp2_Rana <- Box_Int_combined_exp2_Rana / Zi_exp2_Rana + plot_layout(ncol = 1, heights = c(1, 1))
print(combined_plot_exp2_Rana)

# Combine two combined plots side by side
final_combined_plot_exp1and2_Rana <- combined_plot_exp1_Rana | combined_plot_exp2_Rana
print(final_combined_plot_exp1and2_Rana)

#######################
##### BUFO BUFO #######
#######################

### Experiment 1 - Boxplot & Interaction + Zero-Inflation for Bufo ---------------------------

Box_exp1_Bufo <- ggplot(df_Bufo_exp1, aes(x = Phase, y = Active_seconds, fill = Treatment)) +
  geom_boxplot(size = 0.1, alpha = 0.3) +
  facet_grid(rows = ~ Treatment, scales = "free_y", switch = "x", margins = F) +
  scale_y_continuous(limits = c(0, 300)) +
  labs(y = "Active Seconds \n [s]") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 0.5))

pred_cond_exp1_Bufo <- predict(zigam_21_exp1_int_adisPT_zi1_ID_Bufo, type = "conditional", se.fit = TRUE)
pred_zi_exp1_Bufo <- predict(zigam_21_exp1_int_adisPT_zi1_ID_Bufo, type = "zprob", se.fit = TRUE)

df_pred_exp1_Bufo <- data.frame(
  Phase = factor(df_Bufo_exp1$Phase, ordered = TRUE),
  Treatment = df_Bufo_exp1$Treatment,
  Active_seconds = pred_cond_exp1_Bufo$fit,
  Active_seconds_se = pred_cond_exp1_Bufo$se.fit,
  ZeroInflation = pred_zi_exp1_Bufo$fit,
  ZeroInflation_se = pred_zi_exp1_Bufo$se.fit
)

df_calpred_act_exp1_Bufo <- df_pred_exp1_Bufo %>%
  group_by(Treatment, Phase) %>%
  summarise(Active_seconds = mean(Active_seconds, na.rm = TRUE),
            Active_seconds_se = mean(Active_seconds_se, na.rm = TRUE))

df_calpred_zi_exp1_Bufo <- df_pred_exp1_Bufo %>%
  group_by(Treatment, Phase) %>%
  summarise(mean_zi_prob = mean(ZeroInflation, na.rm = TRUE),
            mean_zi_prob_se = mean(ZeroInflation_se, na.rm = TRUE))

Box_Int_combined_Bufo_exp1 <- Box_exp1_Bufo +
  geom_smooth(linewidth = 0.5, method = "lm", se = FALSE, data = df_calpred_act_exp1_Bufo, aes(x = Phase, y = Active_seconds, group = Treatment, color = Treatment)) +
  geom_point(size = 1, data = df_calpred_act_exp1_Bufo, aes(x = Phase, y = Active_seconds)) +
  geom_errorbar(linewidth = 0.5, width = 0.2, data = df_calpred_act_exp1_Bufo, aes(x = Phase, ymin = Active_seconds - Active_seconds_se, ymax = Active_seconds + Active_seconds_se, color = Treatment)) +
  labs(y = "Active Seconds \n [s]", color = "Treatment", fill = "Treatment")

Zi_exp1_Bufo <- ggplot(df_calpred_zi_exp1_Bufo, aes(x = Phase, y = mean_zi_prob, group = Treatment, color = Treatment)) +
  geom_smooth(linewidth = 0.5, method = "lm", se = FALSE) +
  geom_point(size = 1) +
  geom_errorbar(linewidth = 0.5, width = 0.2, aes(ymin = mean_zi_prob - mean_zi_prob_se, ymax = mean_zi_prob + mean_zi_prob_se)) +
  facet_grid(rows = ~ Treatment, scales = "free_y", switch = "x", margins = F) +
  scale_y_continuous(limits = c(-0.05, 1)) +
  labs(x = "Phase", y = "Probability of Freezing \n [mean + SE]") +
  theme_bw() +
  theme(legend.position = "none")

combined_plot_Bufo_exp1 <- Box_Int_combined_Bufo_exp1 / Zi_exp1_Bufo + plot_layout(ncol = 1, heights = c(1, 1))
print(combined_plot_Bufo_exp1)

### Experiment 2 - Boxplot & Interaction + Zero-Inflation for Bufo ---------------------------

Box_exp2_Bufo <- ggplot(df_Bufo_exp2, aes(x = Phase, y = Active_seconds, fill = Treatment)) +
  geom_boxplot(size = 0.1, alpha = 0.3) +
  facet_grid(rows = ~ Treatment, scales = "free_y", switch = "x", margins = F) +
  scale_y_continuous(limits = c(0, 300)) +
  labs(y = "Active Seconds \n [s]") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 0.5))

pred_cond_exp2_Bufo <- predict(zigam_21_exp2_int_adisPT_zi1_ID_Bufo, type = "conditional", se.fit = TRUE)
pred_zi_exp2_Bufo <- predict(zigam_21_exp2_int_adisPT_zi1_ID_Bufo, type = "zprob", se.fit = TRUE)

df_pred_exp2_Bufo <- data.frame(
  Phase = factor(df_Bufo_exp2$Phase, ordered = TRUE),
  Treatment = df_Bufo_exp2$Treatment,
  Active_seconds = pred_cond_exp2_Bufo$fit,
  Active_seconds_se = pred_cond_exp2_Bufo$se.fit,
  ZeroInflation = pred_zi_exp2_Bufo$fit,
  ZeroInflation_se = pred_zi_exp2_Bufo$se.fit
)

df_calpred_act_exp2_Bufo <- df_pred_exp2_Bufo %>%
  group_by(Treatment, Phase) %>%
  summarise(Active_seconds = mean(Active_seconds, na.rm = TRUE),
            Active_seconds_se = mean(Active_seconds_se, na.rm = TRUE))

df_calpred_zi_exp2_Bufo <- df_pred_exp2_Bufo %>%
  group_by(Treatment, Phase) %>%
  summarise(mean_zi_prob = mean(ZeroInflation, na.rm = TRUE),
            mean_zi_prob_se = mean(ZeroInflation_se, na.rm = TRUE))

Box_Int_combined_Bufo_exp2 <- Box_exp2_Bufo +
  geom_smooth(linewidth = 0.5, method = "lm", se = FALSE, data = df_calpred_act_exp2_Bufo, aes(x = Phase, y = Active_seconds, group = Treatment, color = Treatment)) +
  geom_point(size = 1, data = df_calpred_act_exp2_Bufo, aes(x = Phase, y = Active_seconds)) +
  geom_errorbar(linewidth = 0.5, width = 0.2, data = df_calpred_act_exp2_Bufo, aes(x = Phase, ymin = Active_seconds - Active_seconds_se, ymax = Active_seconds + Active_seconds_se, color = Treatment)) +
  labs(y = "Active Seconds \n [s]", color = "Treatment", fill = "Treatment")

Zi_exp2_Bufo <- ggplot(df_calpred_zi_exp2_Bufo, aes(x = Phase, y = mean_zi_prob, group = Treatment, color = Treatment)) +
  geom_smooth(linewidth = 0.5, method = "lm", se = FALSE) +
  geom_point(size = 1) +
  geom_errorbar(linewidth = 0.5, width = 0.2, aes(ymin = mean_zi_prob - mean_zi_prob_se, ymax = mean_zi_prob + mean_zi_prob_se)) +
  facet_grid(rows = ~ Treatment, scales = "free_y", switch = "x", margins = F) +
  scale_y_continuous(limits = c(-0.05, 1)) +
  labs(x = "Phase", y = "Probability of Freezing \n [mean + SE]") +
  theme_bw() +
  theme(legend.position = "none")

combined_plot_Bufo_exp2 <- Box_Int_combined_Bufo_exp2 / Zi_exp2_Bufo + plot_layout(ncol = 1, heights = c(1, 1))
print(combined_plot_Bufo_exp2)

# Combine two combined plots side by side
final_combined_plot_Bufo <- combined_plot_Bufo_exp1 | combined_plot_Bufo_exp2
print(final_combined_plot_Bufo)

# Combine all plots in a grid with 2 rows and 2 columns
final_combined_plot_BufoRana <- (combined_plot_Bufo_exp1 | combined_plot_Bufo_exp2) / (combined_plot_exp1_Rana | combined_plot_exp2_Rana)
print(final_combined_plot_BufoRana)

























