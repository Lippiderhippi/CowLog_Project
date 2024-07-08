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
df_Rana_L$Phase <- factor(df_Rana_L$Phase, levels = c(1, 2, 3), labels = c("1", "2", "3"))
levels(df_Rana_L$Phase)
# this changes the variable type from whatever it is to a factor and relabels them as "1", "2" etc.
# I did this because "phase" as "int" caused a problem during the analysis not recognizing "phase" as an ordinal variable
str(df_Rana_L)


# ### 2.2 Excluding "lazy" individuals #### -------------------------------
# Since some animals were inactive in all three phases, these individuals do not hold any useful information for the analysis
filtered_df_Rana_L <- df_Rana_L %>%
  group_by(Individual_Total) %>%
  filter(!(all(Active_seconds[Phase == "1"] == 0) & 
            all(Active_seconds[Phase == "2"] == 0) & 
             all(Active_seconds[Phase == "3"] == 0))) %>%
  ungroup()
# this excludes all individuals that were inactive in all phases 
str(filtered_df_Rana_L)


# ### 2.3 Split up the data #### ------------------------------------------
df_Rana_exp1 <- filter(filtered_df_Rana_L, Filter_1 == 1)
df_Rana_exp2 <- filter(filtered_df_Rana_L, Filter_1 == 2)
# more elegant way of filtering is e.g. ===> df_Rana_exp1 <- df_Rana_exp1 %>%  filter(Filter_1 == 1)
str(df_Rana_exp1)
str(df_Rana_exp2)


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

# Now some models from very complex to less complex #



# ### 3.1.1 Modelling                   - Experiment 1 #### -----------------

### 1st Models ###
# Models I think, are appropriate #
# probability of zeros (Active_seconds == 0) is modeled using Phase and Treatment as predictors #
zigam_exp1_int_ziPT_Ba_ID <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Batch) + (1 | Individual_Total),
                                     ziformula = ~ Phase + Treatment,
                                     #dispformula = ~ Phase + Treatment,  # Allows for heteroscedasticity
                                     data = df_Rana_exp1,
                                     family = ziGamma(link = "log"))

zigam_exp1_int_ziPT_ID    <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Individual_Total),
                                     ziformula = ~ Phase + Treatment,
                                     #dispformula = ~ Phase + Treatment,  # Allows for heteroscedasticity
                                     data = df_Rana_exp1,
                                     family = ziGamma(link = "log"))

zigam_exp1_int_ziPT_Ba    <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Batch),
                                     ziformula = ~ Phase + Treatment,
                                     data = df_Rana_exp1,
                                     family = ziGamma(link = "log"))

zigam_exp1_int_ziPT       <- glmmTMB(Active_seconds ~ Phase * Treatment,
                                     ziformula = ~ Phase + Treatment,
                                     data = df_Rana_exp1,
                                     family = ziGamma(link = "log"))


### 2nd Models ###
# Models I think, could be appropriate #
# probability of zeros (Active_seconds == 0) is modeled using Phase and Treatment as predictors #
# No interaction between Phase and Treatment #
zigam_exp1_ziPT_Ba_ID     <- glmmTMB(Active_seconds ~ Phase + Treatment + (1 | Batch) + (1 | Individual_Total),
                                     ziformula = ~ Phase + Treatment,
                                     data = df_Rana_exp1,
                                     family = ziGamma(link = "log"))

zigam_exp1_ziPT_ID        <- glmmTMB(Active_seconds ~ Phase + Treatment + (1 | Individual_Total),
                                     ziformula = ~ Phase + Treatment,
                                     data = df_Rana_exp1,
                                     family = ziGamma(link = "log"))

zigam_exp1_ziPT_Ba        <- glmmTMB(Active_seconds ~ Phase + Treatment + (1 | Batch),
                                     ziformula = ~ Phase + Treatment,
                                     data = df_Rana_exp1,
                                     family = ziGamma(link = "log"))

zigam_exp1_ziPT           <- glmmTMB(Active_seconds ~ Phase + Treatment,
                                     ziformula = ~ Phase + Treatment,
                                     data = df_Rana_exp1,
                                     family = ziGamma(link = "log"))


### 3rd Models ###
# Models I think, are not appropriate #
# probability of zeros (Active_seconds == 0) is the same across Phase and Treatment as predictors #

zigam_exp1_int_zi1_Ba_ID    <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Batch) + (1 | Individual_Total),
                                       ziformula = ~ 1,
                                       data = df_Rana_exp1,
                                       family = ziGamma(link = "log"))

zigam_exp1_int_zi1_ID       <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Individual_Total),
                                       ziformula = ~ 1,
                                       data = df_Rana_exp1,
                                       family = ziGamma(link = "log"))

zigam_exp1_int_zi1_Ba       <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Batch),
                                       ziformula = ~ 1,
                                       data = df_Rana_exp1,
                                       family = ziGamma(link = "log"))

zigam_exp1_int_zi1          <- glmmTMB(Active_seconds ~ Phase * Treatment,
                                       ziformula = ~ 1,
                                       data = df_Rana_exp1,
                                       family = ziGamma(link = "log"))


### 4th Models ###
# Models I think, are not appropriate #
# probability of zeros (Active_seconds == 0) is the same across Phase and Treatment as predictors #
# No interaction between Phase and Treatment #

zigam_exp1_zi1_Ba_ID      <- glmmTMB(Active_seconds ~ Phase + Treatment + (1 | Batch) + (1 | Individual_Total),
                                     ziformula = ~ 1,
                                     data = df_Rana_exp1,
                                     family = ziGamma(link = "log"))

zigam_exp1_zi1_ID         <- glmmTMB(Active_seconds ~ Phase + Treatment + (1 | Individual_Total),
                                     ziformula = ~ 1,
                                     data = df_Rana_exp1,
                                     family = ziGamma(link = "log"))

zigam_exp1_zi1_Ba         <- glmmTMB(Active_seconds ~ Phase + Treatment + (1 | Batch),
                                     ziformula = ~ 1,
                                     data = df_Rana_exp1,
                                     family = ziGamma(link = "log"))

zigam_exp1_zi1            <- glmmTMB(Active_seconds ~ Phase + Treatment,
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
compare_performance(zigam_exp1_int_ziPT_Ba_ID, zigam_exp1_int_ziPT_ID, zigam_exp1_int_ziPT_Ba, zigam_exp1_int_ziPT,
                    zigam_exp1_ziPT_Ba_ID, zigam_exp1_ziPT_ID, zigam_exp1_ziPT_Ba, zigam_exp1_ziPT,
                    zigam_exp1_int_zi1_Ba_ID, zigam_exp1_int_zi1_ID, zigam_exp1_int_zi1_Ba, zigam_exp1_int_zi1,
                    zigam_exp1_zi1_Ba_ID, zigam_exp1_zi1_ID, zigam_exp1_zi1_Ba, zigam_exp1_zi1,
                    rank = T)
# the following plot visualizes the performance differences
plot(compare_performance(zigam_exp1_int_ziPT_Ba_ID, zigam_exp1_int_ziPT_ID, zigam_exp1_int_ziPT_Ba, zigam_exp1_int_ziPT,
                         zigam_exp1_ziPT_Ba_ID, zigam_exp1_ziPT_ID, zigam_exp1_ziPT_Ba, zigam_exp1_ziPT,
                         zigam_exp1_int_zi1_Ba_ID, zigam_exp1_int_zi1_ID, zigam_exp1_int_zi1_Ba, zigam_exp1_int_zi1,
                         zigam_exp1_zi1_Ba_ID, zigam_exp1_zi1_ID, zigam_exp1_zi1_Ba, zigam_exp1_zi1))

# values in the center of the spider web plot indicate low performance while values on the edge reflect high performance
# some models might be good in some values but not other (e.g. over fitted models could have a high BIC  but low AIC)


# ### 3.3.1 Model Performance Indices   - Experiment 1 #### -----------------
### 1st Models ###
model_performance(zigam_exp1_int_ziPT_Ba_ID)
model_performance(zigam_exp1_int_ziPT_ID)
model_performance(zigam_exp1_int_ziPT_Ba)
model_performance(zigam_exp1_int_ziPT)

### 2nd Models ###
model_performance(zigam_exp1_ziPT_Ba_ID)
model_performance(zigam_exp1_ziPT_ID)
model_performance(zigam_exp1_ziPT_Ba)
model_performance(zigam_exp1_ziPT)

### 3rd Models ###
model_performance(zigam_exp1_int_zi1_Ba_ID)
model_performance(zigam_exp1_int_zi1_ID)
model_performance(zigam_exp1_int_zi1_Ba)
model_performance(zigam_exp1_int_zi1)

### 4th Models ###
model_performance(zigam_exp1_zi1_Ba_ID)
model_performance(zigam_exp1_zi1_ID)
model_performance(zigam_exp1_zi1_Ba)
model_performance(zigam_exp1_zi1)



# ### 3.4.1 Check Model Assumptions     - Experiment 1 #### ----------------
# Visual check of model assumptions can be done with several models and model types

### 1st Models ###
check_model(zigam_exp1_int_ziPT_Ba_ID)
check_model(zigam_exp1_int_ziPT_ID)
check_model(zigam_exp1_int_ziPT_Ba)
check_model(zigam_exp1_int_ziPT)

### 2nd Models ###
check_model(zigam_exp1_ziPT_Ba_ID)
check_model(zigam_exp1_ziPT_ID)
check_model(zigam_exp1_ziPT_Ba)
check_model(zigam_exp1_ziPT)

### 3rd Models ###
check_model(zigam_exp1_int_zi1_Ba_ID)
check_model(zigam_exp1_int_zi1_ID)
check_model(zigam_exp1_int_zi1_Ba)
check_model(zigam_exp1_int_zi1)

### 4th Models ###
check_model(zigam_exp1_zi1_Ba_ID)
check_model(zigam_exp1_zi1_ID)
check_model(zigam_exp1_zi1_Ba)
check_model(zigam_exp1_zi1)


# ### 3.5.1 Model Summary               - Experiment 1 #### -----------------
### Checking Model Summary will give you parameter coefficients ###
#keep in mind that since releveling the resulst should be in comparison to "C1" (the Control)

### 1st Models ###
summary(zigam_exp1_int_ziPT_Ba_ID)
summary(zigam_exp1_int_ziPT_ID)
summary(zigam_exp1_int_ziPT_Ba)
summary(zigam_exp1_int_ziPT)

### 2nd Models ###
summary(zigam_exp1_ziPT_Ba_ID)
summary(zigam_exp1_ziPT_ID)
summary(zigam_exp1_ziPT_Ba)
summary(zigam_exp1_ziPT)

### 3rd Models ###
summary(zigam_exp1_int_zi1_Ba_ID)
summary(zigam_exp1_int_zi1_ID)
summary(zigam_exp1_int_zi1_Ba)
summary(zigam_exp1_int_zi1)

### 4th Models ###
summary(zigam_exp1_zi1_Ba_ID)
summary(zigam_exp1_zi1_ID)
summary(zigam_exp1_zi1_Ba)
summary(zigam_exp1_zi1)






# ### 3.1.2 Modelling                   - Experiment 2 #### -----------------

### 1st Models ###
# Models I think, are appropriate #
# probability of zeros (Active_seconds == 0) is modeled using Phase and Treatment as predictors #
zigam_exp2_int_ziPT_Ba_ID <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Batch) + (1 | Individual_Total),
                                     ziformula = ~ Phase + Treatment,
                                     data = df_Rana_exp2,
                                     family = ziGamma(link = "log"))

zigam_exp2_int_ziPT_ID    <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Individual_Total),
                                     ziformula = ~ Phase + Treatment,
                                     data = df_Rana_exp2,
                                     family = ziGamma(link = "log"))

zigam_exp2_int_ziPT_Ba    <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Batch),
                                     ziformula = ~ Phase + Treatment,
                                     data = df_Rana_exp2,
                                     family = ziGamma(link = "log"))

zigam_exp2_int_ziPT       <- glmmTMB(Active_seconds ~ Phase * Treatment,
                                     ziformula = ~ Phase + Treatment,
                                     data = df_Rana_exp2,
                                     family = ziGamma(link = "log"))


### 2nd Models ###
# Models I think, could be appropriate #
# probability of zeros (Active_seconds == 0) is modeled using Phase and Treatment as predictors #
# No interaction between Phase and Treatment #
zigam_exp2_ziPT_Ba_ID     <- glmmTMB(Active_seconds ~ Phase + Treatment + (1 | Batch) + (1 | Individual_Total),
                                     ziformula = ~ Phase + Treatment,
                                     data = df_Rana_exp2,
                                     family = ziGamma(link = "log"))

zigam_exp2_ziPT_ID        <- glmmTMB(Active_seconds ~ Phase + Treatment + (1 | Individual_Total),
                                     ziformula = ~ Phase + Treatment,
                                     data = df_Rana_exp2,
                                     family = ziGamma(link = "log"))

zigam_exp2_ziPT_Ba        <- glmmTMB(Active_seconds ~ Phase + Treatment + (1 | Batch),
                                     ziformula = ~ Phase + Treatment,
                                     data = df_Rana_exp2,
                                     family = ziGamma(link = "log"))

zigam_exp2_ziPT           <- glmmTMB(Active_seconds ~ Phase + Treatment,
                                     ziformula = ~ Phase + Treatment,
                                     data = df_Rana_exp2,
                                     family = ziGamma(link = "log"))


### 3rd Models ###
# Models I think, are not appropriate #
# probability of zeros (Active_seconds == 0) is the same across Phase and Treatment as predictors #

zigam_exp2_int_zi1_Ba_ID    <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Batch) + (1 | Individual_Total),
                                       ziformula = ~ 1,
                                       data = df_Rana_exp2,
                                       family = ziGamma(link = "log"))

zigam_exp2_int_zi1_ID       <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Individual_Total),
                                       ziformula = ~ 1,
                                       data = df_Rana_exp2,
                                       family = ziGamma(link = "log"))

zigam_exp2_int_zi1_Ba       <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Batch),
                                       ziformula = ~ 1,
                                       data = df_Rana_exp2,
                                       family = ziGamma(link = "log"))

zigam_exp2_int_zi1          <- glmmTMB(Active_seconds ~ Phase * Treatment,
                                       ziformula = ~ 1,
                                       data = df_Rana_exp2,
                                       family = ziGamma(link = "log"))


### 4th Models ###
# Models I think, are not appropriate #
# probability of zeros (Active_seconds == 0) is the same across Phase and Treatment as predictors #
# No interaction between Phase and Treatment #

zigam_exp2_zi1_Ba_ID      <- glmmTMB(Active_seconds ~ Phase + Treatment + (1 | Batch) + (1 | Individual_Total),
                                     ziformula = ~ 1,
                                     data = df_Rana_exp2,
                                     family = ziGamma(link = "log"))

zigam_exp2_zi1_ID         <- glmmTMB(Active_seconds ~ Phase + Treatment + (1 | Individual_Total),
                                     ziformula = ~ 1,
                                     data = df_Rana_exp2,
                                     family = ziGamma(link = "log"))

zigam_exp2_zi1_Ba         <- glmmTMB(Active_seconds ~ Phase + Treatment + (1 | Batch),
                                     ziformula = ~ 1,
                                     data = df_Rana_exp2,
                                     family = ziGamma(link = "log"))

zigam_exp2_zi1            <- glmmTMB(Active_seconds ~ Phase + Treatment,
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
compare_performance(zigam_exp2_int_ziPT_Ba_ID, zigam_exp2_int_ziPT_ID, zigam_exp2_int_ziPT_Ba, zigam_exp2_int_ziPT,
                    zigam_exp2_ziPT_Ba_ID, zigam_exp2_ziPT_ID, zigam_exp2_ziPT_Ba, zigam_exp2_ziPT,
                    zigam_exp2_int_zi1_Ba_ID, zigam_exp2_int_zi1_ID, zigam_exp2_int_zi1_Ba, zigam_exp2_int_zi1,
                    zigam_exp2_zi1_Ba_ID, zigam_exp2_zi1_ID, zigam_exp2_zi1_Ba, zigam_exp2_zi1,
                    rank = T)

# the following plot visualizes the performance differences
plot(compare_performance(zigam_exp2_int_ziPT_Ba_ID, zigam_exp2_int_ziPT_ID, zigam_exp2_int_ziPT_Ba, zigam_exp2_int_ziPT,
                         zigam_exp2_ziPT_Ba_ID, zigam_exp2_ziPT_ID, zigam_exp2_ziPT_Ba, zigam_exp2_ziPT,
                         zigam_exp2_int_zi1_Ba_ID, zigam_exp2_int_zi1_ID, zigam_exp2_int_zi1_Ba, zigam_exp2_int_zi1,
                         zigam_exp2_zi1_Ba_ID, zigam_exp2_zi1_ID, zigam_exp2_zi1_Ba, zigam_exp2_zi1))

# values in the center of the spider web plot indicate low performance while values on the edge reflect high performance
# some models might be good in some values but not other (e.g. over fitted models could have a high BIC  but low AIC)

# ### 3.3.2 Model Performance Indices   - Experiment 2 #### -----------------

### 1st Models ###
model_performance(zigam_exp2_int_ziPT_Ba_ID)
model_performance(zigam_exp2_int_ziPT_ID)
model_performance(zigam_exp2_int_ziPT_Ba)
model_performance(zigam_exp2_int_ziPT)

### 2nd Models ###
model_performance(zigam_exp2_ziPT_Ba_ID)
model_performance(zigam_exp2_ziPT_ID)
model_performance(zigam_exp2_ziPT_Ba)
model_performance(zigam_exp2_ziPT)

### 3rd Models ###
model_performance(zigam_exp2_int_zi1_Ba_ID)
model_performance(zigam_exp2_int_zi1_ID)
model_performance(zigam_exp2_int_zi1_Ba)
model_performance(zigam_exp2_int_zi1)

### 4th Models ###
model_performance(zigam_exp2_zi1_Ba_ID)
model_performance(zigam_exp2_zi1_ID)
model_performance(zigam_exp2_zi1_Ba)
model_performance(zigam_exp2_zi1)


# ### 3.4.2 Check Model Assumptions     - Experiment 2 #### -----------------

# Visual check of model assumptions can be done with several models and model types
### 1st Models ###
check_model(zigam_exp2_int_ziPT_Ba_ID)
check_model(zigam_exp2_int_ziPT_ID)
check_model(zigam_exp2_int_ziPT_Ba)
check_model(zigam_exp2_int_ziPT)

### 2nd Models ###
check_model(zigam_exp2_ziPT_Ba_ID)
check_model(zigam_exp2_ziPT_ID)
check_model(zigam_exp2_ziPT_Ba)
check_model(zigam_exp2_ziPT)

### 3rd Models ###
check_model(zigam_exp2_int_zi1_Ba_ID)
check_model(zigam_exp2_int_zi1_ID)
check_model(zigam_exp2_int_zi1_Ba)
check_model(zigam_exp2_int_zi1)

### 4th Models ###
check_model(zigam_exp2_zi1_Ba_ID)
check_model(zigam_exp2_zi1_ID)
check_model(zigam_exp2_zi1_Ba)
check_model(zigam_exp2_zi1)


# ### 3.5.2 Model Summary               - Experiment 2 #### -----------------
### Checking Model Summary will give you parameter coefficients ###
#keep in mind that since releveling the resulst should be in comparison to "C2" (the Control)

### 1st Models ###
summary(zigam_exp2_int_ziPT_Ba_ID)
summary(zigam_exp2_int_ziPT_ID)
summary(zigam_exp2_int_ziPT_Ba)
summary(zigam_exp2_int_ziPT)

### 2nd Models ###
summary(zigam_exp2_ziPT_Ba_ID)
summary(zigam_exp2_ziPT_ID)
summary(zigam_exp2_ziPT_Ba)
summary(zigam_exp2_ziPT)

### 3rd Models ###
summary(zigam_exp2_int_zi1_Ba_ID)
summary(zigam_exp2_int_zi1_ID)
summary(zigam_exp2_int_zi1_Ba)
summary(zigam_exp2_int_zi1)

### 4th Models ###
summary(zigam_exp2_zi1_Ba_ID)
summary(zigam_exp2_zi1_ID)
summary(zigam_exp2_zi1_Ba)
summary(zigam_exp2_zi1)









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
