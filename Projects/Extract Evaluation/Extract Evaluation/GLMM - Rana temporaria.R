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
# FactoMineR  = 

# the following package "Pacman" is for package management and makes checking, installing and loading of packages more efficient.
if (!require("pacman")) install.packages("pacman")

# p_load() from {pacman} checks to see if a package is installed.
# If not it attempts to install the package and then loads it. 
# It can also be applied to several packages at once (see below)
pacman::p_load(dplyr, tidyverse, glmmTMB, performance, lme4, patchwork, DHARMa, FactoMineR, interactions, sjPlot, reshape2, broom.mixed)


# ### 1. Load the data #### ------------------------------------------------

df_Rana_L <- read.csv("0. Extract Evaluation_All Rana_v3_L.csv")

# ### 2. Data Structure & Modifications #### ---------------------------------------
# The respective Variables should have the proper type e.g. int (whole numbers), chr (string), num (continuous)
str(df_Rana_L)
# e.g. Phase should not be "int" because is might screw up the analysis as the model might think it is a measurement

# Also check for non-ASCII characters in relevant varaible like - Phase or Treatment
# This can lead to errors in plotting graphs etc. 
any(grepl("[^\x01-\x7F]", df_Rana_L$Phase))
any(grepl("[^\x01-\x7F]", df_Rana_L$Treatment))


# ### 2.1 Converting Variable to Factor #### ------------------------------
# Here I needed to convert Phase to a factor
df_Rana_L$Phase <- factor(df_Rana_L$Phase, levels = c(1, 2, 3), labels = c("Pre", "Stim", "Post"), ordered = TRUE)
#df_Rana_L$Phase <- factor(df_Rana_L$Phase, levels = c(1, 2, 3), labels = c("Pre", "Stim", "Post"))
str(df_Rana_L)

# Convert the ordered factor to numeric
df_Rana_L$Phase_numeric <- as.numeric(df_Rana_L$Phase)
str(df_Rana_L)

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

### Exluding certain Phases or Treatments ###

# This excludes all data for Phase "Stim" or any other
filtered_df_Rana_L <- filtered_df_Rana_L %>% filter(Treatment != "95C")
# This excludes all data for Phase "Post" or any other
filtered_df_Rana_L <- filtered_df_Rana_L %>% filter(Phase != "Post")

#### Combine Variables to reduce collinearity - Phases ###
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
#df_Rana_exp1_1 <- filter(filtered_df_Rana_L, Filter_1_1 == 1)
#df_Rana_exp1_2 <- filter(filtered_df_Rana_L, Filter_1_2 == 1)
# more elegant way of filtering is e.g. ===> df_Rana_exp1 <- df_Rana_exp1 %>%  filter(Filter_1 == 1)
str(df_Rana_exp1)
str(df_Rana_exp2)


#### Centering Phase and Treatment to reduce collinearity ###
#df_Rana_exp1$Phase_centered <- scale(df_Rana_exp1$Phase, center = TRUE, scale = FALSE)
#df_Rana_exp1$Treatment_centered <- scale(df_Rana_exp1$Treatment, center = TRUE, scale = FALSE)


# ### 2.4 Re-level the data  ####  ----------------------------------------
# First set the order of Treatments in the right order so they will be displayed by the plots as you would like them to appear in the plots 
# Assuming Treatment is a factor variable in your data frame df_Rana_exp1
df_Rana_exp1$Treatment <- factor(df_Rana_exp1$Treatment, levels = c("C1", "BFT", "LN", "MS222", "20C", "24h", "65C", "95C", "Prot-K"))
df_Rana_exp2$Treatment <- factor(df_Rana_exp2$Treatment, levels = c("C2", "Arg.2", "Arg.02", "Arg.002", "ArgTric", "Tric.02", "Tric.002"))
df_Rana_exp1$Treatment <- relevel(factor(df_Rana_exp1$Treatment), ref = "C1")
df_Rana_exp2$Treatment <- relevel(factor(df_Rana_exp2$Treatment), ref = "C2")
#df_Rana_exp1_1$Treatment <- relevel(factor(df_Rana_exp1_1$Treatment), ref = "C1")
#df_Rana_exp1_2$Treatment <- relevel(factor(df_Rana_exp1_2$Treatment), ref = "C1")
# This will set "C1 = Control1" as the reference treatment for any subsequent analysis
# I do not know however if this is necessary or harmful
# It seemed to be reasonable from what the outcome is when using this
levels(df_Rana_exp1$Treatment)
levels(df_Rana_exp1$Phase)
levels(df_Rana_exp2$Treatment) 
levels(df_Rana_exp2$Phase)
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
#1 glmmTMB_lognormal  <- glmmTMB(Active_seconds ~ Phase * Treatment, family = lognormal(link = "log"),data = df_Rana_L)
#2 glmmTMB_Gamma      <- glmmTMB(Active_seconds ~ Phase * Treatment, data = df_Rana_L, family = Gamma(link = "log"))
#3 glmmTMB_ziGamma    <- glmmTMB(Active_seconds ~ Phase + Treatment, ziformula = ~ 1, data = df_Rana_exp1, family = ziGamma(link = "log"))

# Important for #3 is that the model includes "ziformula = ~ 1" 
# This argument allows for modelling different probabilities of "zeros" 
# "ziformula = ~ 1" = assumes the probability of "zeros" being the same across all factors
# Since however the probability of "zeros" is expected to be different in respect to the Phase and Treatments this argument has to be changed
# "ziformula = ~ Phase + Treatment" This specifies a zero-inflation model where the probability of zeros (Active_seconds == 0) is modeled using Phase and Treatment as predictors.
# "dispformula = ~ Phase * Treatment" 

# These p-values are relatively high, indicating that these specific interaction terms are not statistically significant at conventional levels (like 0.05). This suggests that the interaction between Phase and Treatment levels does not significantly contribute to explaining the variation in the dispersion of Active_seconds.

# Note: When including (1 | Individual_Total) this tells the model to allow each subject to have their own baseline level of the dependent variable (e.g., seconds). 
# Note: This random intercept varies for each subject, capturing the idea that measurements within the same subject are more similar to each other than to measurements from other subjects.

# ?Generally? # If the majority of interaction terms are not statistically significant (have p-values > 0.05), it suggests that including these interactions in your dispersion model may not be justified.




# ### 3.1.1 Modelling                   - Experiment 1 #### -----------------

### 1st initial Models ###
# Models I think, are appropriate #
# probability of zeros (Active_seconds == 0) is modeled using Phase and Treatment as predictors #

zigam_1_exp1_int_ziPT_Ba_ID <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Individual_Total) + (1 | Batch),
                                          ziformula = ~ Phase + Treatment + (1 | Individual_Total) + (1 | Batch),
                                          data = df_Rana_exp1,
                                          family = ziGamma(link = "log"))
                                    
zigam_2_exp1_int_ziPT_ID    <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Individual_Total),
                                          ziformula = ~ Phase + Treatment + (1 | Individual_Total),
                                          data = df_Rana_exp1,
                                          family = ziGamma(link = "log"))

zigam_3_exp1_int_ziPT_Ba    <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Batch),
                                          ziformula = ~ Phase + Treatment + (1 | Batch),
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
                                          ziformula = ~ Phase + Treatment + (1 | Batch) + (1 | Individual_Total),
                                          data = df_Rana_exp1,
                                          family = ziGamma(link = "log"))

zigam_6_exp1_ziPT_ID        <- glmmTMB(Active_seconds ~ Phase + Treatment + (1 | Individual_Total),
                                          ziformula = ~ Phase + Treatment + (1 | Individual_Total),
                                          data = df_Rana_exp1,
                                          family = ziGamma(link = "log"))

zigam_7_exp1_ziPT_Ba        <- glmmTMB(Active_seconds ~ Phase + Treatment + (1 | Batch),
                                          ziformula = ~ Phase + Treatment + (1 | Batch),
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

### 5th Models ###
# Models I think, are most appropriate #
# probability of zeros (Active_seconds == 0) is modeled using Phase and Treatment as predictors #
# Models that account for Heteroscedasticity #

zigam_17_exp1_int_mdisPT_ziPT_Ba_ID      <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Individual_Total) + (1 | Batch),
                                                    dispformula = ~ Phase * Treatment,
                                                      ziformula = ~ Phase + Treatment + (1 | Individual_Total) + (1 | Batch),
                                                   family = ziGamma(link = "log"),
                                                   data = df_Rana_exp1)
 
zigam_18_exp1_int_mdisPT_ziPT_ID         <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Individual_Total),
                                                    dispformula = ~ Phase * Treatment,
                                                      ziformula = ~ Phase + Treatment + (1 | Individual_Total),
                                                   family = ziGamma(link = "log"),
                                                   data = df_Rana_exp1)

zigam_19_exp1_int_mdisPT_ziPT            <- glmmTMB(Active_seconds ~ Phase * Treatment,
                                                    dispformula = ~ Phase * Treatment,
                                                      ziformula = ~ Phase + Treatment,
                                                   family = ziGamma(link = "log"),
                                                   data = df_Rana_exp1)

zigam_20_exp1_int_adisPT_ziPT_Ba_ID      <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Individual_Total) + (1 | Batch),
                                                    dispformula = ~ Phase + Treatment,
                                                      ziformula = ~ Phase + Treatment + (1 | Individual_Total) + (1 | Batch),
                                                   family = ziGamma(link = "log"),
                                                   data = df_Rana_exp1)

zigam_21_exp1_int_adisPT_ziPT_ID         <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Individual_Total),
                                                    dispformula = ~ Phase + Treatment,
                                                      ziformula = ~ Phase + Treatment + (1 | Individual_Total),
                                                   family = ziGamma(link = "log"),
                                                   data = df_Rana_exp1)

zigam_22_exp1_int_adisPT_ziPT            <- glmmTMB(Active_seconds ~ Phase * Treatment,
                                                    dispformula = ~ Phase + Treatment,
                                                      ziformula = ~ Phase + Treatment,
                                                   family = ziGamma(link = "log"),
                                                   data = df_Rana_exp1)

zigam_23_exp1_mdisPT_ziPT_Ba_ID          <- glmmTMB(Active_seconds ~ Phase + Treatment + (1 | Individual_Total) + (1 | Batch),
                                                    dispformula = ~ Phase * Treatment,
                                                      ziformula = ~ Phase + Treatment + (1 | Individual_Total) + (1 | Batch),
                                                   family = ziGamma(link = "log"),
                                                   data = df_Rana_exp1)

zigam_24_exp1_mdisPT_ziPT_ID             <- glmmTMB(Active_seconds ~ Phase + Treatment + (1 | Individual_Total),
                                                    dispformula = ~ Phase * Treatment,
                                                      ziformula = ~ Phase + Treatment + (1 | Individual_Total),
                                                   family = ziGamma(link = "log"),
                                                   data = df_Rana_exp1)

zigam_25_exp1_mdisPT_ziPT                <- glmmTMB(Active_seconds ~ Phase + Treatment,
                                                    dispformula = ~ Phase * Treatment,
                                                      ziformula = ~ Phase + Treatment,
                                                   family = ziGamma(link = "log"),
                                                   data = df_Rana_exp1)

zigam_26_exp1_adisPT_ziPT_Ba_ID         <- glmmTMB(Active_seconds ~ Phase + Treatment + (1 | Individual_Total) + (1 | Batch),
                                                    dispformula = ~ Phase + Treatment,
                                                      ziformula = ~ Phase + Treatment + (1 | Individual_Total) + (1 | Batch),
                                                   family = ziGamma(link = "log"),
                                                   data = df_Rana_exp1)

zigam_27_exp1_adisPT_ziPT_ID            <- glmmTMB(Active_seconds ~ Phase + Treatment + (1 | Individual_Total),
                                                    dispformula = ~ Phase + Treatment,
                                                      ziformula = ~ Phase + Treatment + (1 | Individual_Total),
                                                   family = ziGamma(link = "log"),
                                                   data = df_Rana_exp1)

zigam_28_exp1_adisPT_ziPT               <- glmmTMB(Active_seconds ~ Phase + Treatment,
                                                    dispformula = ~ Phase + Treatment,
                                                      ziformula = ~ Phase + Treatment,
                                                   family = ziGamma(link = "log"),
                                                   data = df_Rana_exp1)





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

# Models initial 1st - 4th #
compare_performance(zigam_1_exp1_int_ziPT_Ba_ID, zigam_2_exp1_int_ziPT_ID, zigam_3_exp1_int_ziPT_Ba, zigam_4_exp1_int_ziPT,
                    zigam_5_exp1_ziPT_Ba_ID, zigam_6_exp1_ziPT_ID, zigam_7_exp1_ziPT_Ba, zigam_8_exp1_ziPT,
                    zigam_9_exp1_int_zi1_Ba_ID, zigam_10_exp1_int_zi1_ID, zigam_11_exp1_int_zi1_Ba, zigam_12_exp1_int_zi1,
                    zigam_13_exp1_zi1_Ba_ID, zigam_14_exp1_zi1_ID, zigam_15_exp1_zi1_Ba, zigam_16_exp1_zi1,
                    rank = T)

# Models 5th with dispersion #
compare_performance(zigam_17_exp1_int_mdisPT_ziPT_Ba_ID, zigam_18_exp1_int_mdisPT_ziPT_ID, zigam_19_exp1_int_mdisPT_ziPT,
                    zigam_20_exp1_int_adisPT_ziPT_Ba_ID, zigam_21_exp1_int_adisPT_ziPT_ID, zigam_22_exp1_int_adisPT_ziPT,
                    zigam_23_exp1_mdisPT_ziPT_Ba_ID, zigam_24_exp1_mdisPT_ziPT_ID, zigam_25_exp1_mdisPT_ziPT,
                    zigam_26_exp1_adisPT_ziPT_Ba_ID, zigam_27_exp1_adisPT_ziPT_ID, zigam_28_exp1_adisPT_ziPT,
                    rank = T)


# the following plot visualizes the performance differences
# Models initial 1st - 4th #
plot(compare_performance(zigam_1_exp1_int_ziPT_Ba_ID, zigam_2_exp1_int_ziPT_ID, zigam_3_exp1_int_ziPT_Ba, zigam_4_exp1_int_ziPT,
                         zigam_5_exp1_ziPT_Ba_ID, zigam_6_exp1_ziPT_ID, zigam_7_exp1_ziPT_Ba, zigam_8_exp1_ziPT,
                         zigam_9_exp1_int_zi1_Ba_ID, zigam_10_exp1_int_zi1_ID, zigam_11_exp1_int_zi1_Ba, zigam_12_exp1_int_zi1,
                         zigam_13_exp1_zi1_Ba_ID, zigam_14_exp1_zi1_ID, zigam_15_exp1_zi1_Ba, zigam_16_exp1_zi1))

# Models 5th with dispersion #
plot(compare_performance(zigam_17_exp1_int_mdisPT_ziPT_Ba_ID, zigam_18_exp1_int_mdisPT_ziPT_ID, zigam_19_exp1_int_mdisPT_ziPT,
                         zigam_20_exp1_int_adisPT_ziPT_Ba_ID, zigam_21_exp1_int_adisPT_ziPT_ID, zigam_22_exp1_int_adisPT_ziPT,
                         zigam_23_exp1_mdisPT_ziPT_Ba_ID, zigam_24_exp1_mdisPT_ziPT_ID, zigam_25_exp1_mdisPT_ziPT,
                         zigam_26_exp1_adisPT_ziPT_Ba_ID, zigam_27_exp1_adisPT_ziPT_ID, zigam_28_exp1_adisPT_ziPT))



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

# If the LRT p-value is significant, it suggests that the model with the interaction term fits the data better than the model without the interaction term.

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

### 5th Models ###
model_performance(zigam_17_exp1_int_mdisPT_ziPT_Ba_ID)
model_performance(zigam_18_exp1_int_mdisPT_ziPT_ID)
model_performance(zigam_19_exp1_int_mdisPT_ziPT)
model_performance(zigam_20_exp1_int_adisPT_ziPT_Ba_ID)
model_performance(zigam_21_exp1_int_adisPT_ziPT_ID)
model_performance(zigam_22_exp1_int_adisPT_ziPT)
model_performance(zigam_23_exp1_mdisPT_ziPT_Ba_ID)
model_performance(zigam_24_exp1_mdisPT_ziPT_ID)
model_performance(zigam_25_exp1_mdisPT_ziPT)
model_performance(zigam_26_exp1_adisPT_ziPT_Ba_ID)
model_performance(zigam_27_exp1_adisPT_ziPT_ID)
model_performance(zigam_28_exp1_adisPT_ziPT)

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

### 5th Models ###
check_model(zigam_17_exp1_int_mdisPT_ziPT_Ba_ID)
check_model(zigam_18_exp1_int_mdisPT_ziPT_ID)
check_model(zigam_19_exp1_int_mdisPT_ziPT)
check_model(zigam_20_exp1_int_adisPT_ziPT_Ba_ID)
check_model(zigam_21_exp1_int_adisPT_ziPT_ID)
check_model(zigam_22_exp1_int_adisPT_ziPT)
check_model(zigam_23_exp1_mdisPT_ziPT_Ba_ID)
check_model(zigam_24_exp1_mdisPT_ziPT_ID)
check_model(zigam_25_exp1_mdisPT_ziPT)
check_model(zigam_26_exp1_adisPT_ziPT_Ba_ID)
check_model(zigam_27_exp1_adisPT_ziPT_ID)
check_model(zigam_28_exp1_adisPT_ziPT)

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
summary(zigam_optimized)

### 5th Models ###
summary(zigam_17_exp1_int_mdisPT_ziPT_Ba_ID)
summary(zigam_18_exp1_int_mdisPT_ziPT_ID)
summary(zigam_19_exp1_int_mdisPT_ziPT)
summary(zigam_20_exp1_int_adisPT_ziPT_Ba_ID)
summary(zigam_21_exp1_int_adisPT_ziPT_ID)
summary(zigam_22_exp1_int_adisPT_ziPT)
summary(zigam_23_exp1_mdisPT_ziPT_Ba_ID)
summary(zigam_24_exp1_mdisPT_ziPT_ID)
summary(zigam_25_exp1_mdisPT_ziPT)
summary(zigam_26_exp1_adisPT_ziPT_Ba_ID)
summary(zigam_27_exp1_adisPT_ziPT_ID)
summary(zigam_28_exp1_adisPT_ziPT)

# ### Heteroscedasticity checks ###  --------------------------------------

# DHARMa residuals
plot(s1 <- simulateResiduals(zigam_21_exp1_int_adisPT_ziPT_ID))
par(mfrow=c(1,2))
plotResiduals(s1, df_Rana_exp1$Phase, rank = FALSE)
plotResiduals(s1, df_Rana_exp1$Treatment, rank  = FALSE)
# under H0 (perfect model), we would expect those boxes to range homogeneously from 0.25-0.75. 
# To see whether there are deviations from this expectation, the plot calculates a test for uniformity per box, and a test for homgeneity of variances between boxes. 
# A positive test will be highlighted in red.



# ### 3.1.2 Modelling                   - Experiment 2 #### -----------------

### 1st Models ###
# Models I think, are appropriate #
# probability of zeros (Active_seconds == 0) is modeled using Phase and Treatment as predictors #
zigam_1_exp2_int_ziPT_Ba_ID <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Batch) + (1 | Individual_Total),
                                     ziformula = ~ Phase + Treatment + (1 | Batch) + (1 | Individual_Total),
                                     data = df_Rana_exp2,
                                     family = ziGamma(link = "log"))

zigam_2_exp2_int_ziPT_ID    <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Individual_Total),
                                     ziformula = ~ Phase + Treatment + (1 | Individual_Total),
                                     data = df_Rana_exp2,
                                     family = ziGamma(link = "log"))

zigam_3_exp2_int_ziPT_Ba    <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Batch),
                                     ziformula = ~ Phase + Treatment + (1 | Batch),
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
                                     ziformula = ~ Phase + Treatment + (1 | Batch) + (1 | Individual_Total),
                                     data = df_Rana_exp2,
                                     family = ziGamma(link = "log"))

zigam_6_exp2_ziPT_ID        <- glmmTMB(Active_seconds ~ Phase + Treatment + (1 | Individual_Total),
                                     ziformula = ~ Phase + Treatment + (1 | Individual_Total),
                                     data = df_Rana_exp2,
                                     family = ziGamma(link = "log"))

zigam_7_exp2_ziPT_Ba        <- glmmTMB(Active_seconds ~ Phase + Treatment + (1 | Batch),
                                     ziformula = ~ Phase + Treatment + (1 | Batch),
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

### 5th Models ###
# Models I think, are most appropriate #
# probability of zeros (Active_seconds == 0) is modeled using Phase and Treatment as predictors #
# Models that account for Heteroscedasticity #

zigam_17_exp2_int_mdisPT_ziPT_Ba_ID      <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Individual_Total) + (1 | Batch),
                                                   dispformula = ~ Phase * Treatment,
                                                   ziformula = ~ Phase + Treatment + (1 | Individual_Total) + (1 | Batch),
                                                   family = ziGamma(link = "log"),
                                                   data = df_Rana_exp2)

zigam_18_exp2_int_mdisPT_ziPT_ID         <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Individual_Total),
                                                   dispformula = ~ Phase * Treatment,
                                                   ziformula = ~ Phase + Treatment + (1 | Individual_Total),
                                                   family = ziGamma(link = "log"),
                                                   data = df_Rana_exp2)

zigam_19_exp2_int_mdisPT_ziPT            <- glmmTMB(Active_seconds ~ Phase * Treatment,
                                                   dispformula = ~ Phase * Treatment,
                                                   ziformula = ~ Phase + Treatment,
                                                   family = ziGamma(link = "log"),
                                                   data = df_Rana_exp2)

zigam_20_exp2_int_adisPT_ziPT_Ba_ID      <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Individual_Total) + (1 | Batch),
                                                   dispformula = ~ Phase + Treatment,
                                                   ziformula = ~ Phase + Treatment + (1 | Individual_Total) + (1 | Batch),
                                                   family = ziGamma(link = "log"),
                                                   data = df_Rana_exp2)

zigam_21_exp2_int_adisPT_ziPT_ID         <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Individual_Total),
                                                   dispformula = ~ Phase + Treatment,
                                                   ziformula = ~ Phase + Treatment + (1 | Individual_Total),
                                                   family = ziGamma(link = "log"),
                                                   data = df_Rana_exp2)

zigam_22_exp2_int_adisPT_ziPT            <- glmmTMB(Active_seconds ~ Phase * Treatment,
                                                   dispformula = ~ Phase + Treatment,
                                                   ziformula = ~ Phase + Treatment,
                                                   family = ziGamma(link = "log"),
                                                   data = df_Rana_exp2)

zigam_23_exp2_mdisPT_ziPT_Ba_ID          <- glmmTMB(Active_seconds ~ Phase + Treatment + (1 | Individual_Total) + (1 | Batch),
                                                   dispformula = ~ Phase * Treatment,
                                                   ziformula = ~ Phase + Treatment + (1 | Individual_Total) + (1 | Batch),
                                                   family = ziGamma(link = "log"),
                                                   data = df_Rana_exp2)

zigam_24_exp2_mdisPT_ziPT_ID             <- glmmTMB(Active_seconds ~ Phase + Treatment + (1 | Individual_Total),
                                                   dispformula = ~ Phase * Treatment,
                                                   ziformula = ~ Phase + Treatment + (1 | Individual_Total),
                                                   family = ziGamma(link = "log"),
                                                   data = df_Rana_exp2)

zigam_25_exp2_mdisPT_ziPT                <- glmmTMB(Active_seconds ~ Phase + Treatment,
                                                   dispformula = ~ Phase * Treatment,
                                                   ziformula = ~ Phase + Treatment,
                                                   family = ziGamma(link = "log"),
                                                   data = df_Rana_exp2)

zigam_26_exp2_adisPT_ziPT_Ba_ID         <- glmmTMB(Active_seconds ~ Phase + Treatment + (1 | Individual_Total) + (1 | Batch),
                                                   dispformula = ~ Phase + Treatment,
                                                   ziformula = ~ Phase + Treatment + (1 | Individual_Total) + (1 | Batch),
                                                   family = ziGamma(link = "log"),
                                                   data = df_Rana_exp2)

zigam_27_exp2_adisPT_ziPT_ID            <- glmmTMB(Active_seconds ~ Phase + Treatment + (1 | Individual_Total),
                                                   dispformula = ~ Phase + Treatment,
                                                   ziformula = ~ Phase + Treatment + (1 | Individual_Total),
                                                   family = ziGamma(link = "log"),
                                                   data = df_Rana_exp2)

zigam_28_exp2_adisPT_ziPT               <- glmmTMB(Active_seconds ~ Phase + Treatment,
                                                   dispformula = ~ Phase + Treatment,
                                                   ziformula = ~ Phase + Treatment,
                                                   family = ziGamma(link = "log"),
                                                   data = df_Rana_exp2)


# ### 3.2.2 Compare Model Performance   - Experiment 2 #### -----------------

# Troubleshooting #
# 1. "Error: `check_model()` returned following error: input string 3 is invalid in this locale"
# helped get rid of this error and is somewhat linked to the locale settings
Sys.getlocale()
Sys.setlocale("LC_ALL", "C")

# The following command creates a table and ranks the models for the overall performance
# you should still look at AIC and BIC values and consider models that do not perform the best if applicable for good reasons
# Models initial 1st - 4th #
compare_performance(zigam_1_exp2_int_ziPT_Ba_ID, zigam_2_exp2_int_ziPT_ID, zigam_3_exp2_int_ziPT_Ba, zigam_4_exp2_int_ziPT,
                    zigam_5_exp2_ziPT_Ba_ID, zigam_6_exp2_ziPT_ID, zigam_7_exp2_ziPT_Ba, zigam_8_exp2_ziPT,
                    zigam_9_exp2_int_zi1_Ba_ID, zigam_10_exp2_int_zi1_ID, zigam_11_exp2_int_zi1_Ba, zigam_12_exp2_int_zi1,
                    zigam_13_exp2_zi1_Ba_ID, zigam_14_exp2_zi1_ID, zigam_15_exp2_zi1_Ba, zigam_16_exp2_zi1,
                    rank = T)

# Models 5th with dispersion #
compare_performance(zigam_17_exp2_int_mdisPT_ziPT_Ba_ID, zigam_18_exp2_int_mdisPT_ziPT_ID, zigam_19_exp2_int_mdisPT_ziPT,
                    #zigam_20_exp2_int_adisPT_ziPT_Ba_ID, 
                    zigam_21_exp2_int_adisPT_ziPT_ID, zigam_22_exp2_int_adisPT_ziPT,
                    zigam_23_exp2_mdisPT_ziPT_Ba_ID, zigam_24_exp2_mdisPT_ziPT_ID, zigam_25_exp2_mdisPT_ziPT,
                    zigam_26_exp2_adisPT_ziPT_Ba_ID, zigam_27_exp2_adisPT_ziPT_ID, zigam_28_exp2_adisPT_ziPT,
                    rank = T)

# the following plot visualizes the performance differences
# Models initial 1st - 4th #
plot(compare_performance(zigam_1_exp2_int_ziPT_Ba_ID, zigam_2_exp2_int_ziPT_ID, zigam_3_exp2_int_ziPT_Ba, zigam_4_exp2_int_ziPT,
                         zigam_5_exp2_ziPT_Ba_ID, zigam_6_exp2_ziPT_ID, zigam_7_exp2_ziPT_Ba, zigam_8_exp2_ziPT,
                         zigam_9_exp2_int_zi1_Ba_ID, zigam_10_exp2_int_zi1_ID, zigam_11_exp2_int_zi1_Ba, zigam_12_exp2_int_zi1,
                         zigam_13_exp2_zi1_Ba_ID, zigam_14_exp2_zi1_ID, zigam_15_exp2_zi1_Ba, zigam_16_exp2_zi1))

# Models 5th with dispersion #
plot(compare_performance(zigam_17_exp2_int_mdisPT_ziPT_Ba_ID, zigam_18_exp2_int_mdisPT_ziPT_ID, zigam_19_exp2_int_mdisPT_ziPT,
                         zigam_20_exp2_int_adisPT_ziPT_Ba_ID, zigam_21_exp2_int_adisPT_ziPT_ID, zigam_22_exp2_int_adisPT_ziPT,
                         zigam_23_exp2_mdisPT_ziPT_Ba_ID, zigam_24_exp2_mdisPT_ziPT_ID, zigam_25_exp2_mdisPT_ziPT,
                         zigam_26_exp2_adisPT_ziPT_Ba_ID, zigam_27_exp2_adisPT_ziPT_ID, zigam_28_exp2_adisPT_ziPT))



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

### 5th Models ###
model_performance(zigam_17_exp2_int_mdisPT_ziPT_Ba_ID)
model_performance(zigam_18_exp2_int_mdisPT_ziPT_ID)
model_performance(zigam_19_exp2_int_mdisPT_ziPT)
model_performance(zigam_20_exp2_int_adisPT_ziPT_Ba_ID)
model_performance(zigam_21_exp2_int_adisPT_ziPT_ID)
model_performance(zigam_22_exp2_int_adisPT_ziPT)
model_performance(zigam_23_exp2_mdisPT_ziPT_Ba_ID)
model_performance(zigam_24_exp2_mdisPT_ziPT_ID)
model_performance(zigam_25_exp2_mdisPT_ziPT)
model_performance(zigam_26_exp2_adisPT_ziPT_Ba_ID)
model_performance(zigam_27_exp2_adisPT_ziPT_ID)
model_performance(zigam_28_exp2_adisPT_ziPT)



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

### 5th Models ###
check_model(zigam_17_exp2_int_mdisPT_ziPT_Ba_ID)
check_model(zigam_18_exp2_int_mdisPT_ziPT_ID)
check_model(zigam_19_exp2_int_mdisPT_ziPT)
check_model(zigam_20_exp2_int_adisPT_ziPT_Ba_ID)
check_model(zigam_21_exp2_int_adisPT_ziPT_ID)
check_model(zigam_22_exp2_int_adisPT_ziPT)
check_model(zigam_23_exp2_mdisPT_ziPT_Ba_ID)
check_model(zigam_24_exp2_mdisPT_ziPT_ID)
check_model(zigam_25_exp2_mdisPT_ziPT)
check_model(zigam_26_exp2_adisPT_ziPT_Ba_ID)
check_model(zigam_27_exp2_adisPT_ziPT_ID)
check_model(zigam_28_exp2_adisPT_ziPT)


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

### 5th Models ###
summary(zigam_17_exp2_int_mdisPT_ziPT_Ba_ID)
summary(zigam_18_exp2_int_mdisPT_ziPT_ID)
summary(zigam_19_exp2_int_mdisPT_ziPT)
summary(zigam_20_exp2_int_adisPT_ziPT_Ba_ID)
summary(zigam_21_exp2_int_adisPT_ziPT_ID)
summary(zigam_22_exp2_int_adisPT_ziPT)
summary(zigam_23_exp2_mdisPT_ziPT_Ba_ID)
summary(zigam_24_exp2_mdisPT_ziPT_ID)
summary(zigam_25_exp2_mdisPT_ziPT)
summary(zigam_26_exp2_adisPT_ziPT_Ba_ID)
summary(zigam_27_exp2_adisPT_ziPT_ID)
summary(zigam_28_exp2_adisPT_ziPT)


# ### Heteroscedasticity checks ###  --------------------------------------

# DHARMa residuals
plot(s1 <- simulateResiduals(zigam_21_exp1_int_adisPT_ziPT_ID),refit = T)
par(mfrow=c(1,2))
plotResiduals(s1, df_Rana_exp2$Phase, rank = FALSE)
plotResiduals(s1, df_Rana_exp2$Treatment, rank  = FALSE)

testZeroInflation(zigam_27_exp2_adisPT_ziPT_ID)

plotConventionalResiduals(zigam_21_exp1_int_adisPT_ziPT_ID)

# under H0 (perfect model), we would expect those boxes to range homogeneously from 0.25-0.75. 
# To see whether there are deviations from this expectation, the plot calculates a test for uniformity per box, and a test for homgeneity of variances between boxes. 
# A positive test will be highlighted in red.


# ### Data Plots ### ------------------------------------------------------

### Boxplot - horizontal #

Box_exp1 <- ggplot(df_Rana_exp1, aes(x = Phase, y = Active_seconds, fill = Treatment)) +
  geom_boxplot(size = 1) +
  facet_grid(rows = ~ Treatment, shrink = T, scales = "fixed", switch = "x", margins = F) +
  labs(title = "Active Seconds Across Phases for Each Treatment",
       x = "Phase",
       y = "Active Seconds [s]") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))  # Rotate x-axis labels to 0 degrees (horizontal)

# Calculate means and standard errors for each combination of Treatment and Phase
interaction_means_exp1 <- df_Rana_exp1 %>%
  group_by(Treatment, Phase) %>%
  summarise(mean_active_seconds = mean(Active_seconds, na.rm = TRUE),
            se_active_seconds = sd(Active_seconds, na.rm = TRUE) / sqrt(n()))

# Interaction plot using ggplot2
Int_exp1 <- ggplot(interaction_means_exp1, aes(x = Phase, y = mean_active_seconds, group = Treatment, color = Treatment)) +
  geom_line(size=1) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin = mean_active_seconds - se_active_seconds, ymax = mean_active_seconds + se_active_seconds), width = 0.1) +
  facet_grid(rows = ~ Treatment, shrink = T, scales = "fixed", switch = "x", margins = F) +
  labs(title = "Interaction Plot of Treatment and Phase on Active Seconds",
       x = "Treatment",
       y = "Mean Active Seconds",
       color = "Treatment") +  # Ensure the legend title is correctly specified
  theme_bw()


# Conditional model predictions
pred_conditional <- predict(zigam_21_exp1_int_adisPT_ziPT_ID, type = "conditional", se.fit = TRUE)
# Zero-inflation probabilities predictions
pred_zi <- predict(zigam_21_exp1_int_adisPT_ziPT_ID, type = "zprob", se.fit = TRUE)


pred_df <- data.frame(
  Phase = df_Rana_exp1$Phase,
  Treatment = df_Rana_exp1$Treatment,
  ActiveSeconds = pred_conditional$fit,
  ActiveSeconds_se = pred_conditional$se.fit,
  ZeroInflation = pred_zi$fit,
  ZeroInflation_se = pred_zi$se.fit
)
pred_df$Phase <- factor(pred_df$Phase, ordered = TRUE)
str(pred_df)

pred_df_zi <- pred_df %>% 
  select(Phase, Treatment, ZeroInflation, ZeroInflation_se) %>% 
  rename(Value = ZeroInflation, se.fit = ZeroInflation_se)


# Calculate means and standard errors for each combination of Treatment and Phase

Pred_df_zi_exp1 <- pred_df_zi %>%
  group_by(Treatment, Phase) %>%
  summarise(mean_zi_prob = mean(Value, na.rm = TRUE),
            se_zi_prob = mean(se.fit, na.rm = TRUE))

# Plot for Zero Inflation Probabilities
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

# Combine plots using patchwork
combined_plot <- Box_exp1 / Int_exp1 / Zi_exp1 + plot_layout(ncol = 1, heights = c(1, 1, 1))

# Display combined plot
combined_plot




Box_exp2 <- ggplot(df_Rana_exp2, aes(x = Phase, y = Active_seconds, fill = Treatment)) +
  geom_boxplot(size = 1) +
  facet_grid(rows = ~ Treatment, shrink = T, scales = "fixed", switch = "x", margins = F) +
  labs(title = "Active Seconds Across Phases for Each Treatment",
       x = "Phase",
       y = "Active Seconds [s]") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))  # Rotate x-axis labels to 0 degrees (horizontal)

# Calculate means and standard errors for each combination of Treatment and Phase
interaction_means_exp2 <- df_Rana_exp2 %>%
  group_by(Treatment, Phase) %>%
  summarise(mean_active_seconds = mean(Active_seconds, na.rm = TRUE),
            se_active_seconds = sd(Active_seconds, na.rm = TRUE) / sqrt(n()))

# Interaction plot using ggplot2
Int_exp2 <- ggplot(interaction_means_exp2, aes(x = Phase, y = mean_active_seconds, group = Treatment, color = Treatment)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_active_seconds - se_active_seconds, ymax = mean_active_seconds + se_active_seconds), width = 0.1) +
  facet_grid(rows = ~ Treatment, shrink = T, scales = "fixed", switch = "x", margins = F) +
  labs(title = "Interaction Plot of Treatment and Phase on Active Seconds",
       x = "Treatment",
       y = "Mean Active Seconds",
       color = "Treatment") +  # Ensure the legend title is correctly specified
  theme_bw()


# Combine plots using patchwork
combined_plot <- Box_exp2 / Int_exp2 + plot_layout(ncol = 1, heights = c(1, 1))

# Display combined plot
combined_plot


--------------






# theme_update()# Predict zero-inflation probabilities
zi_probs <- predict(zigam_21_exp1_int_adisPT_ziPT_ID, newdata = new_data, type = "zprob")

# Add zero-inflation probabilities to the new data
new_data$ZI_Prob <- zi_probs

# Plot zero-inflation probabilities
ggplot(new_data, aes(x = Phase, y = ZI_Prob, color = Treatment, group = Treatment)) +
  geom_line() +
  geom_point() +
  facet_grid(rows = ~ Treatment, shrink = T, scales = "fixed", switch = "x", margins = F)
  labs(title = "Zero-Inflation Probabilities by Phase and Treatment",
       x = "Phase",
       y = "Probability of Zero Activity") +
  theme_bw()


# Conditional model predictions
pred_conditional <- predict(zigam_21_exp1_int_adisPT_ziPT_ID, type = "conditional", se.fit = TRUE)
# Zero-inflation probabilities predictions
pred_zi <- predict(zigam_21_exp1_int_adisPT_ziPT_ID, type = "zprob", se.fit = TRUE)


pred_df <- data.frame(
  Phase = df_Rana_exp1$Phase,
  Treatment = df_Rana_exp1$Treatment,
  ActiveSeconds = pred_conditional$fit,
  ActiveSeconds_se = pred_conditional$se.fit,
  ZeroInflation = pred_zi$fit,
  ZeroInflation_se = pred_zi$se.fit
)

# Separate data frames for Active Seconds and Zero Inflation
pred_df_active <- pred_df %>% 
  select(Phase, Treatment, ActiveSeconds, ActiveSeconds_se) %>% 
  rename(Value = ActiveSeconds, se.fit = ActiveSeconds_se)

pred_df_zi <- pred_df %>% 
  select(Phase, Treatment, ZeroInflation, ZeroInflation_se) %>% 
  rename(Value = ZeroInflation, se.fit = ZeroInflation_se)

# Calculate means and standard errors for each combination of Treatment and Phase
Pred_df_zi_exp1 <- pred_df_zi %>%
  group_by(Treatment, Phase) %>%
  summarise(mean_zi_prob = mean(Value, na.rm = TRUE),
            se_zi_prob = mean(se.fit, na.rm = TRUE))

library(ggplot2)
library(patchwork)

# Plot for Active Seconds
p1 <- ggplot(pred_df_active, aes(x = Phase, y = Value, color = Treatment)) +
  geom_line() +
  geom_ribbon(aes(ymin = Value - se.fit, ymax = Value + se.fit), alpha = 0.2) +
  facet_grid(~ Treatment) +
  theme_bw() +
  labs(
    x = "Phase",
    y = "Active Seconds",
    title = "Effects of Phase and Treatment on Active Seconds"
  ) +
  theme(legend.position = "bottom")

# Plot for Zero Inflation Probabilities
p2 <- ggplot(Pred_df_zi_exp1, aes(x = Phase, y = mean_zi_prob, color = Treatment)) +
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

# Combine plots using patchwork
combined_plot <- p1 / p2 + plot_layout(ncol = 1, heights = c(1, 1))

# Display combined plot
combined_plot


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


BFGS_glmm <- glmmTMB(
  Active_seconds ~ Phase * Treatment + (1 | Individual_Total),
  ziformula = ~ Phase + Treatment,
  data = df_Rana_exp1,
  family = ziGamma(link = "log"),
  control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))

L_BFGS_B_glmm <- glmmTMB(
  Active_seconds ~ Phase * Treatment + (1 | Individual_Total),
  ziformula = ~ Phase + Treatment,
  data = df_Rana_exp1,
  family = ziGamma(link = "log"),
  control = glmmTMBControl(optimizer = optim, optArgs = list(method = "L-BFGS-B")))

CG_glmm <- glmmTMB(
  Active_seconds ~ Phase * Treatment + (1 | Individual_Total),
  ziformula = ~ Phase + Treatment,
  data = df_Rana_exp1,
  family = ziGamma(link = "log"),
  control = glmmTMBControl(optimizer = optim, optArgs = list(method = "CG")))




zigam_more_iter <- glmmTMB(
  Active_seconds ~ Phase * Treatment + (1 | Individual_Total),
  ziformula = ~ Phase + Treatment,
  data = df_Rana_exp1,
  family = ziGamma(link = "log"),
  control = glmmTMBControl(optimizer = optim, optArgs = list(maxit = 1000000))
)

check_model(optimized_glmm)
summary(optimized_glmm)
check_model(zigam_more_iter)


# Generate a new data frame for prediction
new_data <- expand.grid(
  Phase = unique(df_Rana_exp1$Phase),
  Treatment = unique(df_Rana_exp1$Treatment),
  Individual_Total = NA)

# Predict values
new_data$Predicted <- predict(zigam_21_exp1_int_adisPT_ziPT_ID, newdata = new_data, type = "response")

# Plot the effects of Phase and Treatment
ggplot(new_data, aes(x = Phase, y = Predicted, color = Treatment, group = Treatment)) +
  geom_line() +
  geom_point() +
  facet_grid(rows = ~ Treatment, shrink = T, scales = "fixed", switch = "x", margins = F) +
  labs(title = "Effects of Phase and Treatment on Active Seconds",
       x = "Phase",
       y = "Predicted Active Seconds") +
  theme_bw()



# Zero-inflation probabilities predictions
pred_zi <- predict(zigam_21_exp1_int_adisPT_ziPT_ID, type = "zprob", se.fit = TRUE)

# theme_update()# Predict zero-inflation probabilities
zi_probs <- predict(zigam_21_exp1_int_adisPT_ziPT_ID, newdata = new_data, type = "zprob")

# Add zero-inflation probabilities to the new data
new_data$ZI_Prob <- zi_probs

# Plot zero-inflation probabilities
Zi_exp1 <- ggplot(new_data, aes(x = Phase, y = ZI_Prob, color = Treatment, group = Treatment)) +
  geom_line(size=1) +
  geom_point(size=2) +
  facet_grid(rows = ~ Treatment, shrink = T, scales = "fixed", switch = "x", margins = F) +
  labs(title = "Zero-Inflation Probabilities by Phase and Treatment",
       x = "Phase",
       y = "Probability of Zero Activity") +
  theme_bw()

# Combine plots using patchwork
combined_plot <- Box_exp1 / Int_exp1 / Zi_exp1 + plot_layout(ncol = 1, heights = c(1, 1, 1))

# Display combined plot
combined_plot