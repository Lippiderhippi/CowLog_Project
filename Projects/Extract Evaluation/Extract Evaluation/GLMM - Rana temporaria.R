### GLMM - Rana temporaria ###

# This R script is designed to analyze tadpole response behavior using a Generalized Linear Mixed Model (GLMM) with a zero-inflated gamma distribution. 
# The zero-inflated gamma distribution is used to account for the excess zeros in the data and the continuous nature of the non-zero responses.
# The model incorporates both fixed effects (treatment and phases) and random effects (ID of individuals and batches of trials).


# ### Libraries required #### ---------------------------------------------
# pacman      = For simplifying installing, loading and managing R packages
# dplyr       = For data manipulation
# tidyverse   = For data manipulation
# glmmTMB     = For fitting the GLMM with zero-inflated distributions
# performance = For model evaluation and comparison
# ggplot2     = 
# patchwork   = For combining multiple ggplot2 plots into a single cohesive layout
# DHARMa      = For residual diagnostics for hierarchical (multi-level/mixed) regression models
# lme4        = 

# the following package "Pacman" is for package management and makes checking, installing and loading of packages more efficient.
if (!require("pacman")) install.packages("pacman")

# p_load() from {pacman} checks to see if a package is installed.
# If not it attempts to install the package and then loads it. 
# It can also be applied to several packages at once (see below)
pacman::p_load(dplyr, tidyverse, glmmTMB, performance, lme4, patchwork, DHARMa)


# ### 1. Load the data #### ------------------------------------------------

df_Rana_L <- read.csv("C:/Users/Konrad Lipkowski/Desktop/GitHub/CowLog_Project/Projects/Extract Evaluation/Extract Evaluation/0. Extract Evaluation_All Rana_v3_L.csv")

# ### 2. Data Structure & Modifications #### ---------------------------------------
# The respective Variables should have the proper type e.g. int (whole numbers), chr (string), num (continuous)
str(df_Rana_L)
# e.g. Phase should not be "int" because is might screw up the analysis as the model might think it is a measurement


# ### 2.1 Converting Variable to Factor #### ------------------------------
# Here I needed to convert Phase to a factor
df_Rana_L$Phase <- factor(df_Rana_L$Phase, levels = c(1, 2, 3), labels = c("Pre", "Stim", "Post"))
levels(df_Rana_L$Phase)
# this changes the variable type from whatever it is to a factor and relabels them as "1", "2" etc.
# I did this because "phase" as "int" caused a problem during the analysis not recognizing "phase" as an ordinal variable
str(df_Rana_L)


# ### 2.2 Excluding "lazy" individuals #### -------------------------------
# Since some animals were inactive in all three phases, these individuals do not hold any useful information for the analysis
filtered_df_Rana_L <- df_Rana_L %>%
  group_by(Individual_Total) %>%
  filter(!(all(Active_seconds[Phase == "Pre"] == 0) & 
            all(Active_seconds[Phase == "Stim"] == 0) & 
             all(Active_seconds[Phase == "Post"] == 0))) %>%
  ungroup()
# this excludes all individuals that were inactive in all phases 
str(filtered_df_Rana_L)

# This excludes all data for Phase "Stim" or any other
filtered_df_Rana_L <- filtered_df_Rana_L %>%
  filter(Phase != "Stim")
# This excludes all data for Phase "Post" or any other
filtered_df_Rana_L <- filtered_df_Rana_L %>%
  filter(Phase != "Post")

# Combine Variables to reduce collinearity - Phases
filtered_df_Rana_L$Combined_Phase <- ifelse(filtered_df_Rana_L$Phase %in% c("Stim", "Post"), "Stim_Post", "Pre")
# Combine Variables to reduce collinearity - Treatments
filtered_df_Rana_L$Combined_Treatment <- ifelse(data$Treatment %in% c("boiled", "frozen", "aged"), "processed", "unprocessed")
  
# This essentially create a binary or categorical variable that distinguishes between two conditions:
  # "Pre"-Phase: Represents the period before any stimulus is introduced.
  # "Stim_Post"-Phase: Represents both the "Stim" and "Post" phases combined, indicating the period during and after the stimulus introduction.
# Simplification: Combining simplifies the model by reducing the number of distinct phases from three to two. 
# Simplification: This can make it easier to interpret the effects of these phases on tadpole activity.
# Enhanced Interpretability: Instead of separately estimating the effects of "Stim" and "Post," the combined variable allows for analyzing the overall effect of stimulus introduction (both during and after) compared to the period before the stimulus
# Reduced Collinearity: By creating a combined variable, this potentially reduces the collinearity issue between "Stim" and "Post," 
# Reduced Collinearity: This is because they are now treated as a single predictor in the model. 
# Reduced Collinearity: This can lead to more stable and interpretable model estimates.
# Note: Combining categorical variables like "Phase" into a single variable does not inherently change the duration of observation associated with each phase.
# Note: Given that all phases are equal in duration (5 minutes), there is no inherent imbalance in observation times between phases after combining them into "stim_post."
# Note: This ensures that each phase contributes equally to the combined variable, reflecting the consistent observation protocol of your study.

 
# ### 2.3 Split up the data #### ------------------------------------------
df_Rana_exp1 <- filter(filtered_df_Rana_L, Filter_1 == 1)
df_Rana_exp2 <- filter(filtered_df_Rana_L, Filter_1 == 2)
# more elegant way of filtering is e.g. ===> df_Rana_exp1 <- df_Rana_exp1 %>%  filter(Filter_1 == 1)
str(df_Rana_exp1)
str(df_Rana_exp2)


# Centering Phase and Treatment to reduce collinearity
df_Rana_exp1$Phase_centered <- scale(df_Rana_exp1$Phase, center = TRUE, scale = FALSE)
df_Rana_exp1$Treatment_centered <- scale(df_Rana_exp1$Treatment, center = TRUE, scale = FALSE)




# ### 2.4 Re-level the data  ####  ----------------------------------------
df_Rana_exp1$Treatment <- relevel(factor(df_Rana_exp1$Treatment), ref = "C1")
df_Rana_exp2$Treatment <- relevel(factor(df_Rana_exp2$Treatment), ref = "C2")
# This will set "C1 = Control1" as the reference treatment for any subsequent analysis
# I do not know however if this is necessary or harmful
# It seemed to be reasonable from what the outcome is when using this
levels(df_Rana_exp1$Treatment) 
levels(df_Rana_exp2$Treatment) 
# Check if the desired reference treatment is at first place. That is important for subsequent analysis.
# The other treatments can be in any order. R seems to order them according to values and then letters.
str(df_Rana_exp1)
str(df_Rana_exp2)





# ### 3. Data Analysis #### -----------------------------------------------
# Assuming 'df' is your data frame with columns 'Activity_seconds', 'Phase', 'Treatment', and 'Individual_Total'
# Activity_seconds  = dependent variable
# Phase             = fixed factor (independent variable)
# Treatment         = fixed factor (independent variable)
# Individual_Total  = random factor
# Batch             = random factor


# Example Model #
# Structure Example:  model <- glmer(dependent variable ~ Fixed factors * Fixed factor + (1| Random factor) + (1|Random factor), data = dataframe, family = dependent variable distribution(link = "linkfuntion")
# Data example:       model <- glmer(Active_seconds ~ Phase * Treatment + (1 | Individual_Total) + (1 | Batches), data = df_Rana_L, family = lognormal(link = "log"))

# However, my data has a lot of "zeros" which are inherently meaningful for the interpretation
# Most distribution families however cannot deal with "zeros"
# Because of this I need a special kind of analysis that accounts for zero inflation.
# For this I need a separate package

# Load the following package for zero-inflated data sets
library(glmmTMB)
# it might happen that the TMB package (or another) is too new and not working with "glmmTMB".
# If so, just re install the package using... 
install.packages("glmmTMB") 

# Example Models with using glmmTMB #
#1 glmmTMB_lognormal  <- glmmTMB(Active_seconds ~ Phase * Treatment, data = df_Rana_L, family = lognormal(link = "log"))
#2 glmmTMB_Gamma      <- glmmTMB(Active_seconds ~ Phase * Treatment, data = df_Rana_L, family = Gamma(link = "log"))
#3 glmmTMB_ziGamma    <- glmmTMB(Active_seconds ~ Phase + Treatment, ziformula = ~ 1, data = df_Rana_exp1, family = ziGamma(link = "log"))

# Important for #3 is that the model includes "ziformula = ~ 1" 
# This argument allows for modelling different probabilities of "zeros" 
# "ziformula = ~ 1" = assumes the probability of "zeros" being the same across all factors
# Since however the probability of "zeros" is expected to be different in respect to the Phase and Treatments this argument has to be changed
# "ziformula = ~ Phase + Treatment" This specifies a zero-inflation model where the probability of zeros (Active_seconds == 0) is modeled using Phase and Treatment as predictors.

# Note: When including (1 | Individual_Total) this tells the model to allow each subject to have their own baseline level of the dependent variable (e.g., seconds). 
# Note: This random intercept varies for each subject, capturing the idea that measurements within the same subject are more similar to each other than to measurements from other subjects.



# ### 3.1.1 Modelling                   - Experiment 1 #### -----------------

### 1st Models ###
# Models I think, are appropriate #
# probability of zeros (Active_seconds == 0) is modeled using Phase and Treatment as predictors #
zigam_1_exp1_int_ziPT_Ba_ID <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Batch) + (1 | Individual_Total),
                                     ziformula = ~ Phase + Treatment,
                                     #dispformula = ~ Phase + Treatment,  # Allows for heteroscedasticity
                                     data = df_Rana_exp1,
                                     family = ziGamma(link = "log"))
                                    

zigam_2_exp1_int_ziPT_ID    <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Individual_Total),
                                     ziformula = ~ Phase + Treatment,
                                     #dispformula = ~ Phase + Treatment,  # Allows for heteroscedasticity
                                     data = df_Rana_exp1,
                                     family = ziGamma(link = "log"))

zigam_3_exp1_int_ziPT_Ba    <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Batch),
                                     ziformula = ~ Phase + Treatment,
                                     data = df_Rana_exp1,
                                     family = ziGamma(link = "log"))

zigam_4_exp1_int_ziPT       <- glmmTMB(Active_seconds ~ Phase * Treatment,
                                     ziformula = ~ Phase + Treatment,
                                     data = df_Rana_exp1,
                                     family = ziGamma(link = "log"))


### 2nd Models ###
# Models I think, could be appropriate #
# probability of zeros (Active_seconds == 0) is modeled using Phase and Treatment as predictors #
# No interaction between Phase and Treatment #
zigam_5_exp1_ziPT_Ba_ID     <- glmmTMB(Active_seconds ~ Phase + Treatment + (1 | Batch) + (1 | Individual_Total),
                                     ziformula = ~ Phase + Treatment,
                                     data = df_Rana_exp1,
                                     family = ziGamma(link = "log"))

zigam_6_exp1_ziPT_ID        <- glmmTMB(Active_seconds ~ Phase + Treatment + (1 | Individual_Total),
                                     ziformula = ~ Phase + Treatment,
                                     #dispformula = ~ Phase + Treatment,  # Allows for heteroscedasticity
                                     data = df_Rana_exp1,
                                     family = ziGamma(link = "log"))

zigam_7_exp1_ziPT_Ba        <- glmmTMB(Active_seconds ~ Phase + Treatment + (1 | Batch),
                                     ziformula = ~ Phase + Treatment,
                                     data = df_Rana_exp1,
                                     family = ziGamma(link = "log"))

zigam_8_exp1_ziPT           <- glmmTMB(Active_seconds ~ Phase + Treatment,
                                     ziformula = ~ Phase + Treatment,
                                     data = df_Rana_exp1,
                                     family = ziGamma(link = "log"))


### 3rd Models ###
# Models I think, are not appropriate #
# probability of zeros (Active_seconds == 0) is the same across Phase and Treatment as predictors #

zigam_9_exp1_int_zi1_Ba_ID    <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Batch) + (1 | Individual_Total),
                                       ziformula = ~ 1,
                                       data = df_Rana_exp1,
                                       family = ziGamma(link = "log"))

zigam_10_exp1_int_zi1_ID       <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Individual_Total),
                                       ziformula = ~ 1,
                                       data = df_Rana_exp1,
                                       family = ziGamma(link = "log"))

zigam_11_exp1_int_zi1_Ba       <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Batch),
                                       ziformula = ~ 1,
                                       data = df_Rana_exp1,
                                       family = ziGamma(link = "log"))

zigam_12_exp1_int_zi1          <- glmmTMB(Active_seconds ~ Phase * Treatment,
                                       ziformula = ~ 1,
                                       data = df_Rana_exp1,
                                       family = ziGamma(link = "log"))


### 4th Models ###
# Models I think, are not appropriate #
# probability of zeros (Active_seconds == 0) is the same across Phase and Treatment as predictors #
# No interaction between Phase and Treatment #

zigam_13_exp1_zi1_Ba_ID      <- glmmTMB(Active_seconds ~ Phase + Treatment + (1 | Batch) + (1 | Individual_Total),
                                     ziformula = ~ 1,
                                     data = df_Rana_exp1,
                                     family = ziGamma(link = "log"))

zigam_14_exp1_zi1_ID         <- glmmTMB(Active_seconds ~ Phase + Treatment + (1 | Individual_Total),
                                     ziformula = ~ 1,
                                     data = df_Rana_exp1,
                                     family = ziGamma(link = "log"))

zigam_15_exp1_zi1_Ba         <- glmmTMB(Active_seconds ~ Phase + Treatment + (1 | Batch),
                                     ziformula = ~ 1,
                                     data = df_Rana_exp1,
                                     family = ziGamma(link = "log"))

zigam_16_exp1_zi1            <- glmmTMB(Active_seconds ~ Phase + Treatment,
                                     ziformula = ~ 1,
                                     data = df_Rana_exp1,
                                     family = ziGamma(link = "log"))






# ### 3.2.1 Compare Model Performance   - Experiment 1 #### -----------------

https://www.youtube.com/watch?v=EPIxQ5i5oxs
https://easystats.github.io/see/articles/performance.html 
https://easystats.github.io/performance/reference/check_model.html 

# Troubleshooting #
# 1. "Error: `check_model()` returned following error: input string 3 is invalid in this locale"
# helped get rid of this error and is somewhat linked to the locale settings
Sys.getlocale()
Sys.setlocale("LC_ALL", "C")


# The following command creates a table and ranks the models for the overall performance
# you should still look at AIC and BIC values and consider models that do not perform the best if applicable for good reasons
compare_performance(zigam_1_exp1_int_ziPT_Ba_ID, zigam_2_exp1_int_ziPT_ID, zigam_3_exp1_int_ziPT_Ba, zigam_4_exp1_int_ziPT,
                    zigam_5_exp1_ziPT_Ba_ID, zigam_6_exp1_ziPT_ID, zigam_7_exp1_ziPT_Ba, zigam_8_exp1_ziPT,
                    zigam_9_exp1_int_zi1_Ba_ID, zigam_10_exp1_int_zi1_ID, zigam_11_exp1_int_zi1_Ba, zigam_12_exp1_int_zi1,
                    zigam_13_exp1_zi1_Ba_ID, zigam_14_exp1_zi1_ID, zigam_15_exp1_zi1_Ba, zigam_16_exp1_zi1,
                    rank = T)
# the following plot visualizes the performance differences
plot(compare_performance(zigam_1_exp1_int_ziPT_Ba_ID, zigam_2_exp1_int_ziPT_ID, zigam_3_exp1_int_ziPT_Ba, zigam_4_exp1_int_ziPT,
                         zigam_5_exp1_ziPT_Ba_ID, zigam_6_exp1_ziPT_ID, zigam_7_exp1_ziPT_Ba, zigam_8_exp1_ziPT,
                         zigam_9_exp1_int_zi1_Ba_ID, zigam_10_exp1_int_zi1_ID, zigam_11_exp1_int_zi1_Ba, zigam_12_exp1_int_zi1,
                         zigam_13_exp1_zi1_Ba_ID, zigam_14_exp1_zi1_ID, zigam_15_exp1_zi1_Ba, zigam_16_exp1_zi1))

# values in the center of the spider web plot indicate low performance while values on the edge reflect high performance
# some models might be good in some values but not other (e.g. over fitted models could have a high BIC  but low AIC)

# Nagel: Nagelkerke's R² is a pseudo R² statistic used for logistic regression models. 
# Nagel: This value indicates the proportion of variance in the dependent variable that is predictable from the independent variables. 
# Nagel: Higher values suggest a better model fit.

# RMSE (Root Mean Square Error): RMSE is a measure of the differences between predicted values by the model and the actual values.
# RMSE (Root Mean Square Error): Lower RMSE values indicate a better fit of the model to the data.

# Sigma: Sigma typically refers to the standard deviation of the residuals (errors) in the model. 
# Sigma: It indicates the typical error size between observed and predicted values. 
# Sigma: Smaller sigma values suggest a better model fit.

# AIC weights: AIC weights are derived from the Akaike Information Criterion, which assesses the quality of each model relative to the others.
# AIC weights: Higher AIC weights suggest a better model.
# AICc weights: AICc is a version of AIC corrected for small sample sizes. Higher AICc weights indicate a better model.

# BIC: BIC weights are derived from the Bayesian Information Criterion. 
# BIC: It is another criterion for model selection that includes a penalty for the number of parameters in the model to prevent overfitting. 
# BIC: Like AIC weights, BIC weights indicate the relative likelihood of each model being the best. 
# BIC: Higher BIC weights suggest a better model fit.
# Values using the compare_performance() have to be interpreted differently than values using model_performance!!!
# compare_performance:  ~ Higher values = Better fit
# model_performance:    ~ Lower values = Better fit

# Likelihood Ratio Test
anova(zigam_2_exp1_int_ziPT_ID,
      zigam_6_exp1_ziPT_ID)
#If the LRT p-value is significant, it suggests that the model with the interaction term fits the data better than the model without the interaction term.

# ### 3.3.1 Model Performance Indices   - Experiment 1 #### -----------------
### 1st Models ###
model_performance(zigam_1_exp1_int_ziPT_Ba_ID)
model_performance(zigam_2_exp1_int_ziPT_ID)
model_performance(zigam_3_exp1_int_ziPT_Ba)
model_performance(zigam_4_exp1_int_ziPT)

### 2nd Models ###
model_performance(zigam_5_exp1_ziPT_Ba_ID)
model_performance(zigam_6_exp1_ziPT_ID)
model_performance(zigam_7_exp1_ziPT_Ba)
model_performance(zigam_8_exp1_ziPT)

### 3rd Models ###
model_performance(zigam_9_exp1_int_zi1_Ba_ID)
model_performance(zigam_10_exp1_int_zi1_ID)
model_performance(zigam_11_exp1_int_zi1_Ba)
model_performance(zigam_12_exp1_int_zi1)

### 4th Models ###
model_performance(zigam_13_exp1_zi1_Ba_ID)
model_performance(zigam_14_exp1_zi1_ID)
model_performance(zigam_15_exp1_zi1_Ba)
model_performance(zigam_16_exp1_zi1)

# AIC: Relative quality of a statistical model. Lower AIC values indicate a better model

# AICc: A version of AIC corrected for small sample sizes. lower values of AICc indicate a better model. 
# AICc: It is more reliable than AIC when dealing with small datasets.

# BIC: Similar to AIC but with a different penalty for the number of parameters in the model. Lower BIC values indicate a better model

# R2 (cond): Proportion of variance explained by both the fixed and random effects in mixed models. Higher values indicating a better fit
# R2 (cond): It shows how well the model explains the variance in the data, considering both fixed and random effects.
# R2 (marg): Represents the proportion of variance explained by only the fixed effects in mixed models.
# R2 (marg): It only considers the fixed effects part of the model. Higher values indicating a better fit 

# ICC: Measures the proportion of total variance that is attributable to the grouping structure in the data (i.e., variance explained by the random effects).
# ICC: Higher ICC values indicate that a larger proportion of the total variance is due to differences between groups.

# RMSE: The square root of the average squared differences between observed and predicted values. 

# RMSE: Lower RMSE values indicate a better fit, as they suggest that the model's predictions are close to the observed data.
# Sigma: The standard deviation of the residuals (errors) of the model.
# Sigma: Lower sigma values indicate that the residuals are smaller, which means the model fits the data better.


# ### 3.4.1 Check Model Assumptions     - Experiment 1 #### ----------------
# Visual check of model assumptions can be done with several models and model types

### 1st Models ###
check_model(zigam_1_exp1_int_ziPT_Ba_ID)
check_model(zigam_2_exp1_int_ziPT_ID)
check_model(zigam_3_exp1_int_ziPT_Ba)
check_model(zigam_4_exp1_int_ziPT)


### 2nd Models ###
check_model(zigam_5_exp1_ziPT_Ba_ID)
check_model(zigam_6_exp1_ziPT_ID)
check_model(zigam_7_exp1_ziPT_Ba)
check_model(zigam_8_exp1_ziPT)

### 3rd Models ###
check_model(zigam_9_exp1_int_zi1_Ba_ID)
check_model(zigam_10_exp1_int_zi1_ID)
check_model(zigam_11_exp1_int_zi1_Ba)
check_model(zigam_12_exp1_int_zi1)

### 4th Models ###
check_model(zigam_13_exp1_zi1_Ba_ID)
check_model(zigam_14_exp1_zi1_ID)
check_model(zigam_15_exp1_zi1_Ba)
check_model(zigam_16_exp1_zi1)

# Collinearity
# Collinearity occurs when two or more predictor variables in a regression model are highly correlated.
# This means they contain similar information about the variance in the outcome variable.
# This is represented by high "VIF" values
# High Collinearity might inflate parameter uncertainty

# High VIF in Conditional Model
  # Inflated Standard Errors: High collinearity increases the standard errors of the estimated coefficients.
  # Unstable Estimates: The coefficients might change significantly with small changes in the model.
  # Reduced Significance: High collinearity can make it harder to detect significant relationships between predictors and the outcome.

# Low VIF in Zero-Inflated Model
# Low collinearity means the predictors in the zero-inflated component of the model are not highly correlated.
# Stable Estimates: The coefficients are more reliable and stable.
# More Reliable Significance Testing: The significance tests for the predictors are more trustworthy.


# ### 3.5.1 Model Summary               - Experiment 1 #### -----------------
### Checking Model Summary will give you parameter coefficients ###
#keep in mind that since re leveling the results should be in comparison to "C1" (the Control)

### 1st Models ###
summary(zigam_1_exp1_int_ziPT_Ba_ID)
summary(zigam_2_exp1_int_ziPT_ID)
summary(zigam_3_exp1_int_ziPT_Ba)
summary(zigam_4_exp1_int_ziPT)

### 2nd Models ###
summary(zigam_5_exp1_ziPT_Ba_ID)
summary(zigam_6_exp1_ziPT_ID)
summary(zigam_7_exp1_ziPT_Ba)
summary(zigam_8_exp1_ziPT)

### 3rd Models ###
summary(zigam_9_exp1_int_zi1_Ba_ID)
summary(zigam_10_exp1_int_zi1_ID)
summary(zigam_11_exp1_int_zi1_Ba)
summary(zigam_12_exp1_int_zi1)

### 4th Models ###
summary(zigam_13_exp1_zi1_Ba_ID)
summary(zigam_14_exp1_zi1_ID)
summary(zigam_15_exp1_zi1_Ba)
summary(zigam_16_exp1_zi1)






# ### 3.1.2 Modelling                   - Experiment 2 #### -----------------

### 1st Models ###
# Models I think, are appropriate #
# probability of zeros (Active_seconds == 0) is modeled using Phase and Treatment as predictors #
zigam_1_exp2_int_ziPT_Ba_ID <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Batch) + (1 | Individual_Total),
                                     ziformula = ~ Phase + Treatment,
                                     data = df_Rana_exp2,
                                     family = ziGamma(link = "log"))

zigam_2_exp2_int_ziPT_ID    <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Individual_Total),
                                     ziformula = ~ Phase + Treatment,
                                     data = df_Rana_exp2,
                                     family = ziGamma(link = "log"))

zigam_3_exp2_int_ziPT_Ba    <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Batch),
                                     ziformula = ~ Phase + Treatment,
                                     data = df_Rana_exp2,
                                     family = ziGamma(link = "log"))

zigam_4_exp2_int_ziPT      <- glmmTMB(Active_seconds ~ Phase * Treatment,
                                     ziformula = ~ Phase + Treatment,
                                     data = df_Rana_exp2,
                                     family = ziGamma(link = "log"))


### 2nd Models ###
# Models I think, could be appropriate #
# probability of zeros (Active_seconds == 0) is modeled using Phase and Treatment as predictors #
# No interaction between Phase and Treatment #
zigam_5_exp2_ziPT_Ba_ID     <- glmmTMB(Active_seconds ~ Phase + Treatment + (1 | Batch) + (1 | Individual_Total),
                                     ziformula = ~ Phase + Treatment,
                                     data = df_Rana_exp2,
                                     family = ziGamma(link = "log"))

zigam_6_exp2_ziPT_ID        <- glmmTMB(Active_seconds ~ Phase + Treatment + (1 | Individual_Total),
                                     ziformula = ~ Phase + Treatment,
                                     data = df_Rana_exp2,
                                     family = ziGamma(link = "log"))

zigam_7_exp2_ziPT_Ba        <- glmmTMB(Active_seconds ~ Phase + Treatment + (1 | Batch),
                                     ziformula = ~ Phase + Treatment,
                                     data = df_Rana_exp2,
                                     family = ziGamma(link = "log"))

zigam_8_exp2_ziPT           <- glmmTMB(Active_seconds ~ Phase + Treatment,
                                     ziformula = ~ Phase + Treatment,
                                     data = df_Rana_exp2,
                                     family = ziGamma(link = "log"))


### 3rd Models ###
# Models I think, are not appropriate #
# probability of zeros (Active_seconds == 0) is the same across Phase and Treatment as predictors #

zigam_9_exp2_int_zi1_Ba_ID    <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Batch) + (1 | Individual_Total),
                                       ziformula = ~ 1,
                                       data = df_Rana_exp2,
                                       family = ziGamma(link = "log"))

zigam_10_exp2_int_zi1_ID       <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Individual_Total),
                                       ziformula = ~ 1,
                                       data = df_Rana_exp2,
                                       family = ziGamma(link = "log"))

zigam_11_exp2_int_zi1_Ba       <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Batch),
                                       ziformula = ~ 1,
                                       data = df_Rana_exp2,
                                       family = ziGamma(link = "log"))

zigam_12_exp2_int_zi1          <- glmmTMB(Active_seconds ~ Phase * Treatment,
                                       ziformula = ~ 1,
                                       data = df_Rana_exp2,
                                       family = ziGamma(link = "log"))


### 4th Models ###
# Models I think, are not appropriate #
# probability of zeros (Active_seconds == 0) is the same across Phase and Treatment as predictors #
# No interaction between Phase and Treatment #

zigam_13_exp2_zi1_Ba_ID      <- glmmTMB(Active_seconds ~ Phase + Treatment + (1 | Batch) + (1 | Individual_Total),
                                     ziformula = ~ 1,
                                     data = df_Rana_exp2,
                                     family = ziGamma(link = "log"))

zigam_14_exp2_zi1_ID         <- glmmTMB(Active_seconds ~ Phase + Treatment + (1 | Individual_Total),
                                     ziformula = ~ 1,
                                     data = df_Rana_exp2,
                                     family = ziGamma(link = "log"))

zigam_15_exp2_zi1_Ba         <- glmmTMB(Active_seconds ~ Phase + Treatment + (1 | Batch),
                                     ziformula = ~ 1,
                                     data = df_Rana_exp2,
                                     family = ziGamma(link = "log"))

zigam_16_exp2_zi1            <- glmmTMB(Active_seconds ~ Phase + Treatment,
                                     ziformula = ~ 1,
                                     data = df_Rana_exp2,
                                     family = ziGamma(link = "log"))

# ### 3.2.2 Compare Model Performance   - Experiment 2 #### -----------------

# Troubleshooting #
# 1. "Error: `check_model()` returned following error: input string 3 is invalid in this locale"
# helped get rid of this error and is somewhat linked to the locale settings
Sys.getlocale()
Sys.setlocale("LC_ALL", "C")

# The following command creates a table and ranks the models for the overall performance
# you should still look at AIC and BIC values and consider models that do not perform the best if applicable for good reasons
compare_performance(zigam_1_exp2_int_ziPT_Ba_ID, zigam_2_exp2_int_ziPT_ID, zigam_3_exp2_int_ziPT_Ba, zigam_4_exp2_int_ziPT,
                    zigam_5_exp2_ziPT_Ba_ID, zigam_6_exp2_ziPT_ID, zigam_7_exp2_ziPT_Ba, zigam_8_exp2_ziPT,
                    zigam_9_exp2_int_zi1_Ba_ID, zigam_10_exp2_int_zi1_ID, zigam_11_exp2_int_zi1_Ba, zigam_12_exp2_int_zi1,
                    zigam_13_exp2_zi1_Ba_ID, zigam_14_exp2_zi1_ID, zigam_15_exp2_zi1_Ba, zigam_16_exp2_zi1,
                    rank = T)

# the following plot visualizes the performance differences
plot(compare_performance(zigam_1_exp2_int_ziPT_Ba_ID, zigam_2_exp2_int_ziPT_ID, zigam_3_exp2_int_ziPT_Ba, zigam_4_exp2_int_ziPT,
                         zigam_5_exp2_ziPT_Ba_ID, zigam_6_exp2_ziPT_ID, zigam_7_exp2_ziPT_Ba, zigam_8_exp2_ziPT,
                         zigam_9_exp2_int_zi1_Ba_ID, zigam_10_exp2_int_zi1_ID, zigam_11_exp2_int_zi1_Ba, zigam_12_exp2_int_zi1,
                         zigam_13_exp2_zi1_Ba_ID, zigam_14_exp2_zi1_ID, zigam_15_exp2_zi1_Ba, zigam_16_exp2_zi1))

# values in the center of the spider web plot indicate low performance while values on the edge reflect high performance
# some models might be good in some values but not other (e.g. over fitted models could have a high BIC  but low AIC)

# ### 3.3.2 Model Performance Indices   - Experiment 2 #### -----------------

### 1st Models ###
model_performance(zigam_1_exp2_int_ziPT_Ba_ID)
model_performance(zigam_2_exp2_int_ziPT_ID)
model_performance(zigam_3_exp2_int_ziPT_Ba)
model_performance(zigam_4_exp2_int_ziPT)

### 2nd Models ###
model_performance(zigam_5_exp2_ziPT_Ba_ID)
model_performance(zigam_6_exp2_ziPT_ID)
model_performance(zigam_7_exp2_ziPT_Ba)
model_performance(zigam_8_exp2_ziPT)

### 3rd Models ###
model_performance(zigam_9_exp2_int_zi1_Ba_ID)
model_performance(zigam_10_exp2_int_zi1_ID)
model_performance(zigam_11_exp2_int_zi1_Ba)
model_performance(zigam_12_exp2_int_zi1)

### 4th Models ###
model_performance(zigam_13_exp2_zi1_Ba_ID)
model_performance(zigam_14_exp2_zi1_ID)
model_performance(zigam_15_exp2_zi1_Ba)
model_performance(zigam_16_exp2_zi1)


# ### 3.4.2 Check Model Assumptions     - Experiment 2 #### -----------------

# Visual check of model assumptions can be done with several models and model types
### 1st Models ###
check_model(zigam_1_exp2_int_ziPT_Ba_ID)
check_model(zigam_2_exp2_int_ziPT_ID)
check_model(zigam_3_exp2_int_ziPT_Ba)
check_model(zigam_4_exp2_int_ziPT)

### 2nd Models ###
check_model(zigam_5_exp2_ziPT_Ba_ID)
check_model(zigam_6_exp2_ziPT_ID)
check_model(zigam_7_exp2_ziPT_Ba)
check_model(zigam_8_exp2_ziPT)

### 3rd Models ###
check_model(zigam_9_exp2_int_zi1_Ba_ID)
check_model(zigam_10_exp2_int_zi1_ID)
check_model(zigam_11_exp2_int_zi1_Ba)
check_model(zigam_12_exp2_int_zi1)

### 4th Models ###
check_model(zigam_13_exp2_zi1_Ba_ID)
check_model(zigam_14_exp2_zi1_ID)
check_model(zigam_15_exp2_zi1_Ba)
check_model(zigam_16_exp2_zi1)


# ### 3.5.2 Model Summary               - Experiment 2 #### -----------------
### Checking Model Summary will give you parameter coefficients ###
#keep in mind that since releveling the resulst should be in comparison to "C2" (the Control)

### 1st Models ###
summary(zigam_1_exp2_int_ziPT_Ba_ID)
summary(zigam_2_exp2_int_ziPT_ID)
summary(zigam_3_exp2_int_ziPT_Ba)
summary(zigam_4_exp2_int_ziPT)

### 2nd Models ###
summary(zigam_5_exp2_ziPT_Ba_ID)
summary(zigam_6_exp2_ziPT_ID)
summary(zigam_7_exp2_ziPT_Ba)
summary(zigam_8_exp2_ziPT)

### 3rd Models ###
summary(zigam_9_exp2_int_zi1_Ba_ID)
summary(zigam_10_exp2_int_zi1_ID)
summary(zigam_11_exp2_int_zi1_Ba)
summary(zigam_12_exp2_int_zi1)

### 4th Models ###
summary(zigam_13_exp2_zi1_Ba_ID)
summary(zigam_14_exp2_zi1_ID)
summary(zigam_15_exp2_zi1_Ba)
summary(zigam_16_exp2_zi1)









# #### Abstellgleis #### --------------------------------------------------




library(ggplot2)

# Predicted values from the count model
predicted_counts <- predict(zigam_exp1_int_zi_rd, type = "response")

# Predicted probability of structural zeros from the zero-inflation model
predicted_zeros <- predict(zigam_exp1_int_zi_rd, type = "zprob")




trans_df_Rana_exp1 <- df_Rana_exp1
trans_df_Rana_exp2 <- df_Rana_exp2

trans_df_Rana_exp1$Phase <- as.character(trans_df_Rana_exp1$Phase)
trans_df_Rana_exp2$Phase <- as.character(trans_df_Rana_exp2$Phase)

str(trans_df_Rana_exp1)
str(trans_df_Rana_exp2)

lmer_model <- lmer(Active_seconds ~ Phase * Treatment + (1 | Individual_Total), data = df_Rana_L, REML = FALSE)

# Obtain p-values using likelihood ratio test after models
anova(model)


str(df_Rana_exp2)
transformed_df <- df_Rana_L
transformed_df$Phase <- as.character(transformed_df$Phase)
str(transformed_df)

# Load the needed Packages #
library(dplyr)
library(tidyverse)
library(glmmTMB) # package for zero-inflated gamma
library(performance) # package for model perfomance comparisons
library(lme4)
library(patchwork)
library(DHARMa)