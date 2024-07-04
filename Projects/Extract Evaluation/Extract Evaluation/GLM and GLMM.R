### GLM & GLMM ###

install.packages("lme4")
install.packages("lmerTest")
install.packages("emmeans")

library(lme4)
library(lmerTest)
library(emmeans)

df_L <- read.csv("C:/Users/Konrad Lipkowski/Desktop/GitHub/ExtractEvaluation/0. Extract Evaluation_All Rana_v3_L.csv")
str(df_L)

library(dplyr)
library(tidyverse)

filtered_df_L <- df_L %>%
  group_by(Individual_Total) %>%
  filter(!(all(Active_seconds[Phase == 1] == 0) & 
             all(Active_seconds[Phase == 2] == 0) & 
             all(Active_seconds[Phase == 3] == 0))) %>%
  ungroup()
str(filtered_df_L)

# Filtering data based on Filter_2 == 1 (Experiments 1)
# alternative ====> df_exp1 <- df_noNA %>%  filter(Filter_2 == 1) #"more" elegant

df_exp1 <- filter(filtered_df_L, Filter_1 == 1)     #"less" elegant
df_exp2 <- filter(filtered_df_L, Filter_1 == 2)     #"less" elegant

str(df_exp1)
str(df_exp2)

# Convert Phase to a factor
df_exp1$Phase <- factor(df_exp1$Phase, levels = c(1, 2, 3), labels = c("Phase1", "Phase2", "Phase3"))
df_exp2$Phase <- factor(df_exp2$Phase, levels = c(1, 2, 3), labels = c("Phase1", "Phase2", "Phase3"))

#Do not know if necessary/harmful but here I set "C1 = Control1" as the reference for subsequent analysis
# "$" is used to refer to access/refer to a specific column within df
df_exp1$Treatment <- relevel(factor(df_exp1$Treatment), ref = "C1")
df_exp2$Treatment <- relevel(factor(df_exp2$Treatment), ref = "C2")

str(df_exp1)
str(df_exp2)

trans_df_exp1 <- df_exp1
trans_df_exp2 <- df_exp2

trans_df_exp1$Phase <- as.character(trans_df_exp1$Phase)
trans_df_exp2$Phase <- as.character(trans_df_exp2$Phase)

str(trans_df_exp1)
str(trans_df_exp2)

# Assuming 'df' is your dataframe with columns 'Activity_seconds', 'Phase', 'Treatment', and 'Individual_Total'

# Fit the GLMM
lmer_model <- lmer(Active_seconds ~ Phase * Treatment + (1 | Individual_Total), data = df_L, REML = FALSE)
glmer_model <- glmer(Active_seconds ~ Phase * Treatment + (1 | Individual_Total), data = df_L, family = tweedie(link = "log"))

# Check the model summary
summary(model)

# Obtain p-values using likelihood ratio test
anova(model)


install.packages("glmmTMB") # it might happen that the TMB package is new and not working with this glmmTMB package. If so re install this package
library(glmmTMB)

# Example Model with Gamma Distribution #
glmmTMB_lognormal <- glmmTMB(Active_seconds ~ Phase * Treatment, data = df_L, family = lognormal(link = "log"))
glmmTMB_Gamma     <- glmmTMB(Active_seconds ~ Phase * Treatment, data = df_L, family = Gamma(link = "log"))



zigam_exp1 <- glmmTMB(Active_seconds ~ Phase + Treatment,
                          ziformula = ~ 1,  # Model for zero-inflation (intercept only)
                          data = df_exp1,
                          family = ziGamma(link = "log"))

zigam_exp1_int <- glmmTMB(Active_seconds ~ Phase * Treatment,
                          ziformula = ~ 1,  # Model for zero-inflation (intercept only)
                          data = df_exp1,
                          family = ziGamma(link = "log"))

zigam_exp1_zi <- glmmTMB(Active_seconds ~ Phase + Treatment,
                                 ziformula = ~ Phase + Treatment,  # This specifies a zero-inflation model where the probability of zeros (Active_seconds == 0) is modeled using Phase and Treatment as predictors.
                                 data = df_exp1,
                                 family = ziGamma(link = "log"))

zigam_exp1_int_zi <- glmmTMB(Active_seconds ~ Phase * Treatment,
                              ziformula = ~ Phase + Treatment,  # Model for zero-inflation (intercept only)
                              data = df_exp1,
                              family = ziGamma(link = "log"))

zigam_exp1_int_zi_rd <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Batch) + (1 | Individual_Total),
                         ziformula = ~ Phase + Treatment,  # This specifies a zero-inflation model where the probability of zeros (Active_seconds == 0) is modeled using Phase and Treatment as predictors.
                         data = df_exp2,
                         family = ziGamma(link = "log"))

zigam_exp1_int_zi_rd_phaseischr <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Batch) + (1 | Individual_Total),
                                ziformula = ~ Phase + Treatment,  # This specifies a zero-inflation model where the probability of zeros (Active_seconds == 0) is modeled using Phase and Treatment as predictors.
                                data = df_exp1,
                                family = ziGamma(link = "log"))


# Model summary
summary(zigam_exp1)
summary(zigam_exp1_int)
summary(zigam_exp1_zi)
summary(zigam_exp1_int_zi)
summary(zigam_exp1_int_zi_rd)
summary()

library(performance) # https://easystats.github.io/see/articles/performance.html # https://easystats.github.io/performance/reference/check_model.html #
library(lme4)
library(patchwork)
library(DHARMa)

model_performance(zigam_exp1)
model_performance(zigam_exp1_int)
model_performance(zigam_exp1_zi)
model_performance(zigam_exp1_int_zi)
model_performance(zigam_exp1_int_zi_rd)
model_performance(zigam_exp1_int_zi_rd_phaseischr)

compare_performance(zigam_exp1, zigam_exp1_int, zigam_exp1_zi, zigam_exp1_int_zi, zigam_exp1_int_zi_rd, rank = T)
plot(compare_performance(zigam_exp1, zigam_exp1_int, zigam_exp1_zi, zigam_exp1_int_zi, zigam_exp1_int_zi_rd))

check_model(zigam_exp1)
check_model(zigam_exp1_int)
check_model(zigam_exp1_zi)
check_model(zigam_exp1_int_zi)
check_model(zigam_exp1_int_zi_rd)
check_model(zigam_exp1_int_zi_rd_phaseischr)

# Troubleshooting #
# 1. "Error: `check_model()` returned following error: input string 3 is invalid in this locale"
# helped get rid of this error and is somewhat linked to the locale settings
Sys.getlocale()
Sys.setlocale("LC_ALL", "C")

str(df_exp2)
transformed_df <- df_L
transformed_df$Phase <- as.character(transformed_df$Phase)
str(transformed_df)

library(ggplot2)

# Predicted values from the count model
predicted_counts <- predict(zigam_exp1_int_zi_rd, type = "response")

# Predicted probability of structural zeros from the zero-inflation model
predicted_zeros <- predict(zigam_exp1_int_zi_rd, type = "zprob")