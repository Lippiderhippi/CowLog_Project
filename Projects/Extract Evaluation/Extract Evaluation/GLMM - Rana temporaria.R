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
# emmeans     = For calculating/extracting estimated marginal means
# ggplot2     = For plotting data
# patchwork   = For combining multiple ggplot2 plots into a single cohesive layout
# DHARMa      = For residual diagnostics for hierarchical (multi-level/mixed) regression models
# lme4        = 
# FactoMineR  = 

# the following package "Pacman" is for package management and makes checking, installing and loading of packages more efficient.
if (!require("pacman")) install.packages("pacman")

# p_load() from {pacman} checks to see if a package is installed.
# If not it attempts to install the package and then loads it. 
# It can also be applied to several packages at once (see below)
pacman::p_load(dplyr, tidyverse, emmeans, glmmTMB, performance, lme4, patchwork, DHARMa, FactoMineR, interactions, sjPlot, reshape2, broom.mixed, ggplot2)


# ### 1. Load the data #### ------------------------------------------------
getwd() # Check what the working directory is (should be the one where the excelfiles are in)
setwd("C:/Users/Lippi/Desktop/GitRepository/Projects/Extract Evaluation/Extract Evaluation") # Changes the working directory to where the excel files ar in

df_Rana_L <- read.csv("0. Extract Evaluation_All Rana_v3_L.csv")

# ### 2. Data Structure & Modifications #### ---------------------------------------
# The respective Variables should have the proper type e.g. int (whole numbers), chr (string), num (continuous)
# e.g. Phase should not be "int" because is might screw up the analysis as the model might think it is a measurement
# Also check for non-ASCII characters in relevant variable like - Phase or Treatment
# This can lead to errors in plotting graphs etc. 
str(df_Rana_L)
any(grepl("[^\x01-\x7F]", df_Rana_L$Phase))
any(grepl("[^\x01-\x7F]", df_Rana_L$Treatment))


# ### 2.1 Converting Variable to Factor #### ------------------------------
# Here I needed to convert Phase to a factor
df_Rana_L$Phase <- factor(df_Rana_L$Phase, levels = c(1, 2, 3), labels = c("Pre", "Stim", "Post"), ordered = TRUE)
str(df_Rana_L)


# ### 2.2 Excluding "lazy" individuals #### -------------------------------
# Since some animals were inactive in all three phases, these individuals do not hold any useful information for the analysis
# this excludes all individuals that were inactive in all phases 
filtered_df_Rana_L <- df_Rana_L %>%
  group_by(Individual_Total) %>%
  filter(!(all(Active_seconds[Phase == "Pre"] == 0) & 
            all(Active_seconds[Phase == "Stim"] == 0) & 
             all(Active_seconds[Phase == "Post"] == 0))) %>%
  ungroup()

# The following then counts the number of excluded individuals per treatment and summarizes how many remain

initial_individuals_per_treatment <- df_Rana_L %>%
  group_by(Treatment) %>%
  summarise(Initial_Individuals = n_distinct(Individual_Total))

remaining_count <- filtered_df_Rana_L %>%
  distinct(Individual_Total, Treatment) %>%
  group_by(Treatment) %>%
  summarise(Remaining = n())

summary_per_treatment <- initial_individuals_per_treatment %>%
  left_join(remaining_count, by = "Treatment") %>%
  mutate(Remaining = coalesce(Remaining, 0),  # Replace NA with 0 if no exclusions for a treatment
         Excluded = Initial_Individuals - Remaining) %>%
  select(Treatment, Initial_Individuals, Remaining, Excluded)

# Display the summary of how many Individuals were included/excluded
summary_per_treatment



# ### OPTIONAL - Excluding certain Phases or Treatments ### --------------------------
# This is OPTIONAL and was not done for the analysis of this data set

# This excludes all data for Treatment "95C" or any other
#filtered_df_Rana_L <- filtered_df_Rana_L %>% filter(Treatment != "95C")
# This excludes all data for Phase "Post" or any other
#filtered_df_Rana_L <- filtered_df_Rana_L %>% filter(Phase != "Post")



# ### OPTIONAL - Combining to reduce collinearity ### ------------------------------------------
# This is OPTIONAL and was not done for the analysis of this data set

# Combining Phases #
#filtered_df_Rana_L$Combined_Phase <- ifelse(filtered_df_Rana_L$Phase %in% c("Stim", "Post"), "Stim_Post", "Pre")

# Combining Treatments #
#filtered_df_Rana_L$Combined_Treatment <- ifelse(data$Treatment %in% c("boiled", "frozen", "aged"), "processed", "unprocessed")

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
# I split up the data according to the respective available controls 
# i.e. I have a control for Experiments 1 and another control for Experiments 2
df_Rana_exp1 <- filter(filtered_df_Rana_L, Filter_1 == 1)
df_Rana_exp2 <- filter(filtered_df_Rana_L, Filter_1 == 2)
#df_Rana_exp1_1 <- filter(filtered_df_Rana_L, Filter_1_1 == 1)
#df_Rana_exp1_2 <- filter(filtered_df_Rana_L, Filter_1_2 == 1)
str(df_Rana_exp1)
str(df_Rana_exp2)


# ### OPTIONAL - Centering to reducing collinearity ### ------------------------------
#df_Rana_exp1$Phase_centered <- scale(df_Rana_exp1$Phase, center = TRUE, scale = FALSE)
#df_Rana_exp1$Treatment_centered <- scale(df_Rana_exp1$Treatment, center = TRUE, scale = FALSE)


# ### Counting abundances of "Zero" ### ------------------------------------
# This reports the number of "zeros" per treatment and Phase
# This gives an idea of how to model the zero-inflation part of of the model
# e.g. because species might have different "freezing" behaviors or abundances of zeros.
# In this particular case Rana had a lot of "zeros" while Bufo did not.
# Therefore I choose different "ziformula", modeling different "zero" occurrences probabilities for Rana and Bufo 

zero_counts <- df_Rana_exp1 %>%
  mutate(Zero_Active_seconds = ifelse(Active_seconds == 0, 1, 0)) %>%
  group_by(Phase, Treatment) %>%
  summarise(
    Zero_Count = sum(Zero_Active_seconds),
    Total_Count = n(),
    Zero_Proportion = Zero_Count / Total_Count)
print(n = 50, zero_counts) # Might have to be adjusted to show all Treatments and Phases

zero_counts <- df_Rana_exp2 %>%
  mutate(Zero_Active_seconds = ifelse(Active_seconds == 0, 1, 0)) %>%
  group_by(Phase, Treatment) %>%
  summarise(
    Zero_Count = sum(Zero_Active_seconds),
    Total_Count = n(),
    Zero_Proportion = Zero_Count / Total_Count)
print(n = 50, zero_counts) # Might have to be adjusted to show all Treatments and Phases



# ### 2.4 Re-level the data  ####  ----------------------------------------
# First set the order of Treatments in the right order so they will be displayed by the plots as you would like them to appear in the plots 
# Assuming Treatment is a factor variable in your data frame df_Rana_exp1
# This will set "C1 = Control1" or "C2 = Control2" as the reference treatment for any subsequent analysis
# I do not know however if this is necessary or harmful
# It seemed to be reasonable from what the outcome is when using this
df_Rana_exp1$Treatment <- factor(df_Rana_exp1$Treatment, levels = c("C1", "BFT", "LN2", "MS222", "-20C", "24h", "65C", "95C", "Prot-K"))
df_Rana_exp2$Treatment <- factor(df_Rana_exp2$Treatment, levels = c("C2", "Arg.2", "Arg.02", "Arg.002", "ArgTric", "Tric.02", "Tric.002"))
df_Rana_exp1$Treatment <- relevel(factor(df_Rana_exp1$Treatment), ref = "C1")
df_Rana_exp2$Treatment <- relevel(factor(df_Rana_exp2$Treatment), ref = "C2")

# This shows the levels per Factor and the order in which they will be reported
levels(df_Rana_exp1$Treatment)
levels(df_Rana_exp1$Phase)
levels(df_Rana_exp2$Treatment) 
levels(df_Rana_exp2$Phase)
# Check if the desired reference treatment is at first place.
# That is important for subsequent analysis because the first place is always the reference to which the other levels will be compared.
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
# probability of zeros (Active_seconds == 0) is modeled using Phase and Treatment as predictors and random effects
# That means that the probability of "zeros" is different per Phase and Treatments

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
# 1. Interaction Term ~ Phase*Treatment: Is accounting for any baseline imbalances 
# 2. Ziformula: Probability of "zeros" (Active_seconds == 0) is modeled using Phase and Treatment as predictors #
# 2. Ziformula: This means the probability of "zeros" occurring is different across Phases and Treatments #
# 3. Dispformula: Models account for Heteroscedasticity #
# 3. Dispformula: By doing this, I am essentially allowing the dispersion to vary according to one or more covariates.
# 3. Dispformula: This can effectively account for heteroscedasticity in my data
# 3. Dispformula: Note: "dispformula" does not allow the inclusion of random factors

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

zigam_21_exp1_int_adisPT_ziPT_ID_Rana         <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Individual_Total),
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





# ### 3.2.1 Compare Model Performance   - Experiments 1 #### -----------------

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
                    zigam_20_exp1_int_adisPT_ziPT_Ba_ID, zigam_21_exp1_int_adisPT_ziPT_ID_Rana, zigam_22_exp1_int_adisPT_ziPT,
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
                         zigam_20_exp1_int_adisPT_ziPT_Ba_ID, zigam_21_exp1_int_adisPT_ziPT_ID_Rana, zigam_22_exp1_int_adisPT_ziPT,
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
model_performance(zigam_21_exp1_int_adisPT_ziPT_ID_Rana)
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
check_model(zigam_21_exp1_int_adisPT_ziPT_ID_Rana)
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
summary(zigam_21_exp1_int_adisPT_ziPT_ID_Rana)
summary(zigam_22_exp1_int_adisPT_ziPT)
summary(zigam_23_exp1_mdisPT_ziPT_Ba_ID)
summary(zigam_24_exp1_mdisPT_ziPT_ID)
summary(zigam_25_exp1_mdisPT_ziPT)
summary(zigam_26_exp1_adisPT_ziPT_Ba_ID)
summary(zigam_27_exp1_adisPT_ziPT_ID)
summary(zigam_28_exp1_adisPT_ziPT)


# ### Extracting values from the model to report ### ---------------------------------------
# Extract coefficients



### 1st version ###

model_summary <- summary(zigam_21_exp1_int_adisPT_ziPT_ID_Rana)

# Extract coefficients and standard errors for the conditional model
coef_values_cond <- model_summary$coefficients$cond[, "Estimate"]
se_values_cond <- model_summary$coefficients$cond[, "Std. Error"]

# Extract coefficients and standard errors for the zero-inflation model
coef_values_zi <- model_summary$coefficients$zi[, "Estimate"]
se_values_zi <- model_summary$coefficients$zi[, "Std. Error"]

# Extract coefficients and standard errors for the dispersion model
coef_values_disp <- model_summary$coefficients$disp[, "Estimate"]
se_values_disp <- model_summary$coefficients$disp[, "Std. Error"]

# Function to calculate exponentiated coefficient and its standard error
calculate_exp_and_se <- function(coef, se) {
  exp_coef <- exp(coef)
  se_exp_coef <- exp_coef * se
  return(list(exp_coef = exp_coef, se_exp_coef = se_exp_coef))
}

# Calculate exponentiated values and standard errors for conditional model
exp_and_se_cond <- mapply(calculate_exp_and_se, coef_values_cond, se_values_cond, SIMPLIFY = FALSE)

# Calculate exponentiated values and standard errors for zero-inflation model
exp_and_se_zi <- mapply(calculate_exp_and_se, coef_values_zi, se_values_zi, SIMPLIFY = FALSE)

# Calculate exponentiated values and standard errors for dispersion model
exp_and_se_disp <- mapply(calculate_exp_and_se, coef_values_disp, se_values_disp, SIMPLIFY = FALSE)

# Print results for conditional model
cat("Conditional Model Results:\n")
for (i in seq_along(exp_and_se_cond)) {
  cat("Term:", names(exp_and_se_cond)[i], "\n")
  cat("Exponentiated Coefficient:", exp_and_se_cond[[i]]$exp_coef, "\n")
  cat("Standard Error of Exponentiated Coefficient:", exp_and_se_cond[[i]]$se_exp_coef, "\n\n")
}

# Print results for zero-inflation model
cat("Zero-Inflation Model Results:\n")
for (i in seq_along(exp_and_se_zi)) {
  cat("Term:", names(exp_and_se_zi)[i], "\n")
  cat("Exponentiated Coefficient:", exp_and_se_zi[[i]]$exp_coef, "\n")
  cat("Standard Error of Exponentiated Coefficient:", exp_and_se_zi[[i]]$se_exp_coef, "\n\n")
}

# Print results for dispersion model
cat("Dispersion Model Results:\n")
for (i in seq_along(exp_and_se_disp)) {
  cat("Term:", names(exp_and_se_disp)[i], "\n")
  cat("Exponentiated Coefficient:", exp_and_se_disp[[i]]$exp_coef, "\n")
  cat("Standard Error of Exponentiated Coefficient:", exp_and_se_disp[[i]]$se_exp_coef, "\n\n")
}

### 2nd version ###

# Extract coefficients and standard errors for the conditional model
coef_values_cond <- model_summary$coefficients$cond[, "Estimate"]
se_values_cond <- model_summary$coefficients$cond[, "Std. Error"]

# Extract coefficients and standard errors for the zero-inflation model
coef_values_zi <- model_summary$coefficients$zi[, "Estimate"]
se_values_zi <- model_summary$coefficients$zi[, "Std. Error"]

# Extract coefficients and standard errors for the dispersion model
coef_values_disp <- model_summary$coefficients$disp[, "Estimate"]
se_values_disp <- model_summary$coefficients$disp[, "Std. Error"]

# Function to calculate exponentiated coefficient and its standard error
calculate_exp_and_se <- function(coef, se) {
  exp_coef <- exp(coef)
  se_exp_coef <- exp_coef * se
  return(list(exp_coef = exp_coef, se_exp_coef = se_exp_coef))
}

# Calculate exponentiated values and standard errors for conditional model
exp_and_se_cond <- mapply(calculate_exp_and_se, coef_values_cond, se_values_cond, SIMPLIFY = FALSE)

# Calculate exponentiated values and standard errors for zero-inflation model
exp_and_se_zi <- mapply(calculate_exp_and_se, coef_values_zi, se_values_zi, SIMPLIFY = FALSE)

# Calculate exponentiated values and standard errors for dispersion model
exp_and_se_disp <- mapply(calculate_exp_and_se, coef_values_disp, se_values_disp, SIMPLIFY = FALSE)

# Print results for conditional model
cat("Conditional Model Results:\n")
for (i in seq_along(exp_and_se_cond)) {
  cat("Term:", names(exp_and_se_cond)[i], "\n")
  cat("Exponentiated Coefficient:", exp_and_se_cond[[i]]$exp_coef, "\n")
  cat("Standard Error of Exponentiated Coefficient:", exp_and_se_cond[[i]]$se_exp_coef, "\n\n")
}

# Print results for zero-inflation model
cat("Zero-Inflation Model Results:\n")
for (i in seq_along(exp_and_se_zi)) {
  cat("Term:", names(exp_and_se_zi)[i], "\n")
  cat("Exponentiated Coefficient:", exp_and_se_zi[[i]]$exp_coef, "\n")
  cat("Standard Error of Exponentiated Coefficient:", exp_and_se_zi[[i]]$se_exp_coef, "\n\n")
}

# Print results for dispersion model
cat("Dispersion Model Results:\n")
for (i in seq_along(exp_and_se_disp)) {
  cat("Term:", names(exp_and_se_disp)[i], "\n")
  cat("Exponentiated Coefficient:", exp_and_se_disp[[i]]$exp_coef, "\n")
  cat("Standard Error of Exponentiated Coefficient:", exp_and_se_disp[[i]]$se_exp_coef, "\n\n")
}


### version 3 ###

model_summary <- summary(zigam_21_exp1_int_adisPT_ziPT_ID_Rana)

# Function to calculate exponentiated coefficient and its standard error
calculate_exp_and_se <- function(coef, se) {
  exp_coef <- exp(coef)
  se_exp_coef <- exp_coef * se  # Approximate SE for exponentiated coefficient
  return(list(exp_coef = exp_coef, se_exp_coef = se_exp_coef))
}

# Extract coefficients and standard errors for the conditional model
coef_values_cond <- model_summary$coefficients$cond[, "Estimate"]
se_values_cond <- model_summary$coefficients$cond[, "Std. Error"]
exp_and_se_cond <- mapply(calculate_exp_and_se, coef_values_cond, se_values_cond, SIMPLIFY = FALSE)

# Extract coefficients and standard errors for the dispersion model
coef_values_disp <- model_summary$coefficients$disp[, "Estimate"]
se_values_disp <- model_summary$coefficients$disp[, "Std. Error"]
exp_and_se_disp <- mapply(calculate_exp_and_se, coef_values_disp, se_values_disp, SIMPLIFY = FALSE)

# Extract coefficients and standard errors for the zero-inflation model
coef_values_zero <- model_summary$coefficients$zi[, "Estimate"]
se_values_zero <- model_summary$coefficients$zi[, "Std. Error"]
exp_and_se_zero <- mapply(calculate_exp_and_se, coef_values_zero, se_values_zero, SIMPLIFY = FALSE)

# Print results for all model components
print_results <- function(exp_and_se, model_type) {
  cat(paste("\nExponentiated Coefficients and Standard Errors for", model_type, "Model:\n"))
  for (term in names(exp_and_se)) {
    exp_coef <- exp_and_se[[term]]$exp_coef
    se_exp_coef <- exp_and_se[[term]]$se_exp_coef
    cat("Term:", term, "\n")
    cat("Exponentiated Coefficient:", exp_coef, "\n")
    cat("Standard Error of Exponentiated Coefficient:", se_exp_coef, "\n")
    cat("\n")
  }
}

print_results(exp_and_se_cond, "Conditional")
print_results(exp_and_se_disp, "Dispersion")
print_results(exp_and_se_zero, "Zero-Inflation")

# Function to calculate and print combined effects for interaction terms
calculate_interaction_effects <- function(coef_values, se_values, exp_and_se, model_type) {
  interaction_terms <- grep(":", names(exp_and_se), value = TRUE)
  cat(paste("\nCombined Effects for Interaction Terms in", model_type, "Model:\n"))
  for (interaction in interaction_terms) {
    # Split interaction into its components
    parts <- strsplit(interaction, ":")[[1]]
    term1 <- parts[1]
    term2 <- parts[2]
    
    # Check if the terms are in the main effects list
    if (term1 %in% names(exp_and_se) && term2 %in% names(exp_and_se)) {
      term1_coef <- coef_values[term1]
      term1_se <- se_values[term1]
      
      term2_coef <- coef_values[term2]
      term2_se <- se_values[term2]
      
      interaction_coef <- coef_values[interaction]
      interaction_se <- se_values[interaction]
      
      # Calculate combined effect
      combined_coef <- term1_coef + term2_coef + interaction_coef       # cutting out "term2_coef" ???
      exp_combined_coef <- exp(combined_coef)
      
      # Calculate variance for combined effect (assuming errors are independent)
      combined_var <- (term1_se^2) + (term2_se^2) + (interaction_se^2)  # cutting out "term2_coef" ???
      combined_se <- sqrt(combined_var)
      
      # Exponentiate the standard error
      exp_combined_se <- exp(combined_coef) * combined_se
      
      cat("Interaction Term:", interaction, "\n")
      cat("Exponentiated Combined Coefficient:", exp_combined_coef, "\n")
      cat("Standard Error of Exponentiated Combined Coefficient:", exp_combined_se, "\n")
      cat("\n")
    }
  }
}

# Calculate and print combined effects for interaction terms in each model
calculate_interaction_effects(coef_values_cond, se_values_cond, exp_and_se_cond, "Conditional")
calculate_interaction_effects(coef_values_disp, se_values_disp, exp_and_se_disp, "Dispersion")
calculate_interaction_effects(coef_values_zero, se_values_zero, exp_and_se_zero, "Zero-Inflation")



exp(-0.149 - 0.671) # Phase:BFT expected activity after Treatment effect is ~44 seconds
100*exp(-0.149 - 0.174) # Phase:MS222 expected activity after Treatment effect is ~72 seconds

sqrt(0.05809^2 + 0.17498^2)
exp(-0.149)*0.058
exp(0.149 - 0.225)


exp(-0.149 - 0.527)

# Estimated marginal means #
# Log Scale #
emmeans_result_exp1 <- emmeans(zigam_21_exp1_int_adisPT_ziPT_ID_Rana, ~ Phase * Treatment)
summary(emmeans_result_exp1) # Results are given on the log (not the response) scale

# Original Scale #
original_scale_result_exp1 <- summary(emmeans_result_exp1, type = "response") # The type argument specifies the scale on which to return the results.
# Setting type = "response" tells the summary function to transform the estimated marginal means from the link function scale (linear predictor) back to the original response scale.
# type = link           :returns the estimated marginal means (EMMs) on the scale of the linear predictor
# type = response       :returns the EMMs on the original response scale, which means the results are transformed back from the link scale to the original scale of the response variable
# type = prob           :provides the estimated probabilities when the model involves a binary outcome
# type = lp or linpred  :similar to "link" 
# type = cumulative     :Used in ordinal models to provide cumulative probabilities

summary(original_scale_result_exp1) # Results are back-transformed from the log scale


# ### Heteroscedasticity checks ###  --------------------------------------

# DHARMa residuals
plot(s1 <- simulateResiduals(zigam_21_exp1_int_adisPT_ziPT_ID_Rana))
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
# By doing this, I am essentially allowing the dispersion to vary according to one or more covariates.
# This can effectively account for heteroscedasticity in my data
# Note: "dispformula" does not allow the inclusion of random factors

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

zigam_21_exp2_int_adisPT_ziPT_ID_Rana         <- glmmTMB(Active_seconds ~ Phase * Treatment + (1 | Individual_Total),
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
                    zigam_21_exp2_int_adisPT_ziPT_ID_Rana, zigam_22_exp2_int_adisPT_ziPT,
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
                         zigam_20_exp2_int_adisPT_ziPT_Ba_ID, zigam_21_exp2_int_adisPT_ziPT_ID_Rana, zigam_22_exp2_int_adisPT_ziPT,
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
model_performance(zigam_21_exp2_int_adisPT_ziPT_ID_Rana)
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
check_model(zigam_21_exp2_int_adisPT_ziPT_ID_Rana)
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
summary(zigam_21_exp2_int_adisPT_ziPT_ID_Rana)
summary(zigam_22_exp2_int_adisPT_ziPT)
summary(zigam_23_exp2_mdisPT_ziPT_Ba_ID)
summary(zigam_24_exp2_mdisPT_ziPT_ID)
summary(zigam_25_exp2_mdisPT_ziPT)
summary(zigam_26_exp2_adisPT_ziPT_Ba_ID)
summary(zigam_27_exp2_adisPT_ziPT_ID)
summary(zigam_28_exp2_adisPT_ziPT)


# Estimated marginal means #
# Log Scale #
emmeans_result_exp2 <- emmeans(zigam_21_exp2_int_adisPT_ziPT_ID_Rana, ~ Phase * Treatment)
summary(emmeans_result_exp2) # Results are given on the log (not the response) scale

# Original Scale #
original_scale_result_exp2 <- summary(emmeans_result_exp2, type = "response") # The type argument specifies the scale on which to return the results.
# Setting type = "response" tells the summary function to transform the estimated marginal means from the link function scale (linear predictor) back to the original response scale.
# type = link           :returns the estimated marginal means (EMMs) on the scale of the linear predictor
# type = response       :returns the EMMs on the original response scale, which means the results are transformed back from the link scale to the original scale of the response variable
# type = prob           :provides the estimated probabilities when the model involves a binary outcome
# type = lp or linpred  :similar to "link" 
# type = cumulative     :Used in ordinal models to provide cumulative probabilities

summary(original_scale_result_exp2) # Results are back-transformed from the log scale

# ### Heteroscedasticity checks ###  --------------------------------------

# DHARMa residuals
plot(s1 <- simulateResiduals(zigam_21_exp1_int_adisPT_ziPT_ID_Rana),refit = T)
par(mfrow=c(1,2))
plotResiduals(s1, df_Rana_exp2$Phase, rank = FALSE)
plotResiduals(s1, df_Rana_exp2$Treatment, rank  = FALSE)

testZeroInflation(zigam_27_exp2_adisPT_ziPT_ID)

plotConventionalResiduals(zigam_21_exp1_int_adisPT_ziPT_ID_Rana)

# under H0 (perfect model), we would expect those boxes to range homogeneously from 0.25-0.75. 
# To see whether there are deviations from this expectation, the plot calculates a test for uniformity per box, and a test for homgeneity of variances between boxes. 
# A positive test will be highlighted in red.


# ### Data Plots ### ------------------------------------------------------


# ### Experiment 1 - Boxplot & EMMs + Zero-Inflation ### -----------

# Calculate Estimated Marginal Means #
# Log Scale #
emmeans_result_exp1_Rana <- emmeans(zigam_21_exp1_int_adisPT_ziPT_ID_Rana, ~ Phase * Treatment)
summary(emmeans_result_exp1_Rana)
# Original Scale #
original_scale_result_exp1_Rana <- summary(emmeans_result_exp1_Rana, type = "response") # The type argument specifies the scale on which to return the results.
# Setting type = "response" tells the summary function to transform the estimated marginal means from the link function scale (linear predictor) back to the original response scale.
summary(original_scale_result_exp1_Rana)


# Convert to data frame for easier plotting
df_emm_exp1_Rana <- as.data.frame(original_scale_result_exp1_Rana)
# Check colum names and rename them (see below). 
# I did this because I got an error saying "Active_seconds" was not found.
# In theroy this should not cause any problems because each layer "df_rana_exp1" and "df_emm_exp1" were given the correct colum names for the values.
# For "df_Rana_exp1"  = Active_seconds
# For "df_emm_exp1"   = response
# Because of this annoying error I just renamed the colums in "df_emm_exp1" to include "Active_seconds"
# Now the code stops bitching around...
str(df_emm_exp1_Rana) 
names(df_emm_exp1_Rana)
colnames(df_emm_exp1_Rana) <- c("Phase", "Treatment", "Active_seconds", "SE", "df", "Lower_CI", "Upper_CI") # Renaming the colums
df_emm_exp1_Rana$Phase <- factor(df_emm_exp1_Rana$Phase, ordered = TRUE)

# Plot for Boxplots with Emms
Box_Emm_exp1_Rana <- ggplot(df_Rana_exp1, aes(x = Phase, y = Active_seconds, fill = Treatment)) +
  geom_boxplot(size = 0.1, alpha = 0.3, outlier.shape = 16, outlier.color = "black", outlier.size = 1) +  # Box plot for original data
  geom_point(data = df_emm_exp1_Rana, aes(x = Phase, y = Active_seconds, color = Treatment), size = 1, shape = 21) +  # EMM points
  geom_smooth(linewidth = 0.5, alpha = 1, method = "lm", se = FALSE, data = df_emm_exp1_Rana, aes(x = Phase, y = Active_seconds, group = Treatment, color = Treatment)) +
  geom_errorbar(data = df_emm_exp1_Rana, aes(x = Phase, ymin = Lower_CI, ymax = Upper_CI, color = Treatment), width = 0.2) +  # Error bars
  facet_grid(rows = ~ Treatment, shrink = T, scales = "free_y", switch = "x", margins = F) +
  scale_y_continuous(limits = c(0, 200)) +
  labs(title = "",
       x = "",
       y = "Active Seconds \n [s]") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0, hjust = 0.5))  # Rotate x-axis labels to 0 degrees (horizontal)

# Extracting Values from Model predictions #
# Conditional model value predictions
#pred_cond_exp1 <- predict(zigam_21_exp1_int_adisPT_ziPT_ID_Rana, type = "conditional", se.fit = TRUE)
# Zero-inflation probabilities value predictions
pred_zi_exp1_Rana <- predict(zigam_21_exp1_int_adisPT_ziPT_ID_Rana, type = "zprob", se.fit = TRUE)

#Create a new dataframe for calculations and subsequent plotting
df_pred_exp1_Rana <- data.frame(
  Phase = df_Rana_exp1$Phase,
  Treatment = df_Rana_exp1$Treatment,
  #Active_seconds = pred_cond_exp1$fit,
  #Active_seconds_se = pred_cond_exp1$se.fit,
  ZeroInflation = pred_zi_exp1_Rana$fit,
  ZeroInflation_se = pred_zi_exp1_Rana$se.fit)

df_pred_exp1_Rana$Phase <- factor(df_pred_exp1_Rana$Phase, ordered = TRUE)

#df_pred_act_exp1 <- df_pred_exp1 %>% 
#select(Phase, Treatment, Active_seconds, Active_seconds_se)

df_pred_zi_exp1_Rana <- df_pred_exp1_Rana %>% 
  select(Phase, Treatment, ZeroInflation, ZeroInflation_se)

# Calculate means and standard errors for each combination of Treatment and Phase
#df_calpred_act_exp1 <- df_pred_act_exp1 %>%
# group_by(Treatment, Phase) %>%
# summarise(Active_seconds = mean(Active_seconds, na.rm = TRUE),
#  Active_seconds_se = mean(Active_seconds_se, na.rm = TRUE))

df_calpred_zi_exp1_Rana <- df_pred_zi_exp1_Rana %>%
  group_by(Treatment, Phase) %>%
  summarise(mean_zi_prob = mean(ZeroInflation, na.rm = TRUE),
            mean_zi_prob_se = mean(ZeroInflation_se, na.rm = TRUE))

# Plot for Zero Inflation Probabilities
Zi_exp1_Rana <- ggplot(df_calpred_zi_exp1_Rana, aes(x = Phase, y = mean_zi_prob, group = Treatment, color = Treatment)) +
  geom_smooth(linewidth = 0.5, method = "lm", se = FALSE) +
  geom_point(size = 1) +
  geom_errorbar(linewidth = 0.5, width = 0.2, aes(ymin = mean_zi_prob - mean_zi_prob_se, ymax = mean_zi_prob + mean_zi_prob_se)) +
  facet_grid(rows = ~ Treatment, shrink = T, scales = "free_y", switch = "x", margins = F) +
  scale_y_continuous(limits = c(-0.05, 1)) +
  labs(
    x = "Phase",
    y = "Propability of freezing \n [mean +SE]",
    title = ""
  ) +
  theme_bw() +
  theme(legend.position = "none")

# Combine plots using patchwork
combined_plot_exp1_Rana <- Box_Emm_exp1_Rana / Zi_exp1_Rana + plot_layout(ncol = 1, heights = c(1, 1))

# Display combined plot
print(combined_plot_exp1_Rana)

# Print individual the plots
print(Box_Emm_exp1_Rana)
print(Zi_exp1_Rana)


# ### Experiment 2 - Boxplot & EMMs + Zero-Inflation ### -----------

# Calculate Estimated Marginal Means #
# Log Scale #
emmeans_result_exp2_Rana <- emmeans(zigam_21_exp2_int_adisPT_ziPT_ID_Rana, ~ Phase * Treatment)
summary(emmeans_result_exp2_Rana)
# Original Scale #
original_scale_result_exp2_Rana <- summary(emmeans_result_exp2_Rana, type = "response") # The type argument specifies the scale on which to return the results.
# Setting type = "response" tells the summary function to transform the estimated marginal means from the link function scale (linear predictor) back to the original response scale.
summary(original_scale_result_exp2_Rana)


# Convert to data frame for easier plotting
df_emm_exp2_Rana <- as.data.frame(original_scale_result_exp2_Rana)
# Check colum names and rename them (see below). 
# I did this because I got an error saying "Active_seconds" was not found.
# In theroy this should not cause any problems because each layer "df_rana_exp2" and "df_emm_exp2" were given the correct colum names for the values.
# For "df_Rana_exp2"  = Active_seconds
# For "df_emm_exp2"   = response
# Because of this annoying error I just renamed the colums in "df_emm_exp2" to include "Active_seconds"
# Now the code stops bitching around...
str(df_emm_exp2_Rana) 
names(df_emm_exp2_Rana)
colnames(df_emm_exp2_Rana) <- c("Phase", "Treatment", "Active_seconds", "SE", "df", "Lower_CI", "Upper_CI") # Renaming the colums
df_emm_exp2_Rana$Phase <- factor(df_emm_exp2_Rana$Phase, ordered = TRUE)

# Plot for Boxplots with Emms
Box_Emm_exp2_Rana <- ggplot(df_Rana_exp2, aes(x = Phase, y = Active_seconds, fill = Treatment)) +
  geom_boxplot(size = 0.1, alpha = 0.3, outlier.shape = 16, outlier.color = "black", outlier.size = 1) +  # Box plot for original data
  geom_point(data = df_emm_exp2_Rana, aes(x = Phase, y = Active_seconds, color = Treatment), size = 1, shape = 21) +  # EMM points
  geom_smooth(linewidth = 0.5, alpha = 1, method = "lm", se = FALSE, data = df_emm_exp2_Rana, aes(x = Phase, y = Active_seconds, group = Treatment, color = Treatment)) +
  geom_errorbar(data = df_emm_exp2_Rana, aes(x = Phase, ymin = Lower_CI, ymax = Upper_CI, color = Treatment), width = 0.2) +  # Error bars
  facet_grid(rows = ~ Treatment, shrink = T, scales = "free_y", switch = "x", margins = F) +
  scale_y_continuous(limits = c(0, 200)) +
  labs(title = "",
       x = "",
       y = "Active Seconds \n [s]") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0, hjust = 0.5))  # Rotate x-axis labels to 0 degrees (horizontal)

# Extracting Values from Model predictions #
# Conditional model value predictions
#pred_cond_exp2 <- predict(zigam_21_exp2_int_adisPT_ziPT_ID_Rana, type = "conditional", se.fit = TRUE)
# Zero-inflation probabilities value predictions
pred_zi_exp2_Rana <- predict(zigam_21_exp2_int_adisPT_ziPT_ID_Rana, type = "zprob", se.fit = TRUE)

#Create a new dataframe for calculations and subsequent plotting
df_pred_exp2_Rana <- data.frame(
  Phase = df_Rana_exp2$Phase,
  Treatment = df_Rana_exp2$Treatment,
  #Active_seconds = pred_cond_exp2$fit,
  #Active_seconds_se = pred_cond_exp2$se.fit,
  ZeroInflation = pred_zi_exp2_Rana$fit,
  ZeroInflation_se = pred_zi_exp2_Rana$se.fit)

df_pred_exp2_Rana$Phase <- factor(df_pred_exp2_Rana$Phase, ordered = TRUE)

#df_pred_act_exp2 <- df_pred_exp2 %>% 
#select(Phase, Treatment, Active_seconds, Active_seconds_se)

df_pred_zi_exp2_Rana <- df_pred_exp2_Rana %>% 
  select(Phase, Treatment, ZeroInflation, ZeroInflation_se)

# Calculate means and standard errors for each combination of Treatment and Phase
#df_calpred_act_exp2 <- df_pred_act_exp2 %>%
# group_by(Treatment, Phase) %>%
# summarise(Active_seconds = mean(Active_seconds, na.rm = TRUE),
#  Active_seconds_se = mean(Active_seconds_se, na.rm = TRUE))

df_calpred_zi_exp2_Rana <- df_pred_zi_exp2_Rana %>%
  group_by(Treatment, Phase) %>%
  summarise(mean_zi_prob = mean(ZeroInflation, na.rm = TRUE),
            mean_zi_prob_se = mean(ZeroInflation_se, na.rm = TRUE))

# Plot for Zero Inflation Probabilities
Zi_exp2_Rana <- ggplot(df_calpred_zi_exp2_Rana, aes(x = Phase, y = mean_zi_prob, group = Treatment, color = Treatment)) +
  geom_smooth(linewidth = 0.5, method = "lm", se = FALSE) +
  geom_point(size = 1) +
  geom_errorbar(linewidth = 0.5, width = 0.2, aes(ymin = mean_zi_prob - mean_zi_prob_se, ymax = mean_zi_prob + mean_zi_prob_se)) +
  facet_grid(rows = ~ Treatment, shrink = T, scales = "free_y", switch = "x", margins = F) +
  scale_y_continuous(limits = c(-0.05, 1)) +
  labs(
    x = "Phase",
    y = "Propability of freezing \n [mean +SE]",
    title = ""
  ) +
  theme_bw() +
  theme(legend.position = "none")

# Combine plots using patchwork
combined_plot_exp2_Rana <- Box_Emm_exp2_Rana / Zi_exp2_Rana + plot_layout(ncol = 1, heights = c(1, 1))

# Display combined plot
print(combined_plot_exp2_Rana)

# Print individual the plots
print(Box_Emm_exp2_Rana)
print(Zi_exp2_Rana)

### Combine two combined plots ###

# Combine two combined plots side by side
final_combined_plot_Rana <- combined_plot_exp1_Rana | combined_plot_exp2_Rana

# Display the final combined plot
print(final_combined_plot_Rana)


# Combine all plots in a grid with 2 rows and 2 columns
final_combined_plot_BufoRana <- (combined_plot_exp1_Bufo | combined_plot_exp2_Bufo) / (combined_plot_exp1_Rana | combined_plot_exp2_Rana)

# Display the final combined plot
print(final_combined_plot_BufoRana)



# ### Experiment 1 - Boxplot & Interaction + Zero-Inflation ### -----------

# Boxplots #
Box_exp1 <- ggplot(df_Rana_exp1, aes(x = Phase, y = Active_seconds, fill = Treatment)) +
  geom_boxplot(size = 0.1, alpha = 0.3) +  # Adjust alpha for overlay visibility
  facet_grid(rows = ~ Treatment, shrink = T, scales = "free_y", switch = "x", margins = F) +
  scale_y_continuous(limits = c(0, 200)) +
  labs(title = "",
       x = "",
       y = "Active Seconds \n [s]") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))  # Rotate x-axis labels to 0 degrees (horizontal)

# Extracting Values from Model predictions #
# Conditional model value predictions
pred_cond_exp1 <- predict(zigam_21_exp1_int_adisPT_ziPT_ID_Rana, type = "conditional", se.fit = TRUE)
# Zero-inflation probabilities value predictions
pred_zi_exp1 <- predict(zigam_21_exp1_int_adisPT_ziPT_ID_Rana, type = "zprob", se.fit = TRUE)

#Create a new dataframe for calculations and subsequent plotting
df_pred_exp1 <- data.frame(
  Phase = df_Rana_exp1$Phase,
  Treatment = df_Rana_exp1$Treatment,
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
  geom_smooth(linewidth = 2, alpha = 2, method = "lm", se = FALSE, data = df_calpred_act_exp1, aes(x = Phase, y = Active_seconds, group = Treatment, color = Treatment)) +
  geom_point(size = 2, alpha = 2, data = df_calpred_act_exp1, aes(x = Phase, y = Active_seconds)) +
  geom_errorbar(linewidth = 1, width = 0.2, data = df_calpred_act_exp1, aes(x = Phase, ymin = Active_seconds - Active_seconds_se, ymax = Active_seconds + Active_seconds_se, color = Treatment)) +
  labs(title = "",
       y = "Active Seconds \n [s]",
       color = "Treatment",
       fill = "Treatment") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

# Plot for Zero Inflation Probabilities
Zi_exp1 <- ggplot(df_calpred_zi_exp1, aes(x = Phase, y = mean_zi_prob, group = Treatment, color = Treatment)) +
  geom_smooth(linewidth = 1, method = "lm", se = FALSE) +
  geom_point(size = 2) +
  geom_errorbar(linewidth = 1, width = 0.2, aes(ymin = mean_zi_prob - mean_zi_prob_se, ymax = mean_zi_prob + mean_zi_prob_se)) +
  facet_grid(rows = ~ Treatment, shrink = T, scales = "free_y", switch = "x", margins = F) +
  scale_y_continuous(limits = c(-0.05, 1)) +
  labs(
    x = "Phase",
    y = "Propability of freezing \n [mean +SE]",
    title = ""
  ) +
  theme_bw() +
  theme(legend.position = "none")

# Combine plots using patchwork
combined_plot <- Box_Int_combined / Zi_exp1 + plot_layout(ncol = 1, heights = c(1, 1))

# Display combined plot
print(combined_plot)



library(dplyr)

# Generate predictions with error handling
tryCatch({
  pred_cond_exp1 <- predict(zigam_21_exp1_int_adisPT_ziPT_ID_Rana, type = "conditional", se.fit = TRUE)
  pred_zi_exp1 <- predict(zigam_21_exp1_int_adisPT_ziPT_ID_Rana, type = "zprob", se.fit = TRUE)
}, error = function(e) {
  stop("Prediction failed: ", e$message)
})

# Create a dataframe with individual predictions
individual_pred_df <- data.frame(
  Phase = df_Rana_exp1$Phase,
  Treatment = df_Rana_exp1$Treatment,
  Active_seconds = pred_cond_exp1$fit,  # Keep on log scale
  Active_seconds_se = pred_cond_exp1$se.fit,  # Keep SE on log scale
  ZeroInflation = pred_zi_exp1$fit,  # Assuming this is on the probability scale
  ZeroInflation_se = pred_zi_exp1$se.fit
)

# Calculate mean and standard deviation for each Phase within each Treatment
summary_pred_df <- individual_pred_df %>%
  group_by(Phase, Treatment) %>%
  summarize(
    Mean_Active_seconds = mean(Active_seconds, na.rm = TRUE),
    SD_Active_seconds = sd(Active_seconds, na.rm = TRUE),
    Mean_ZeroInflation = mean(ZeroInflation, na.rm = TRUE),
    SD_ZeroInflation = sd(ZeroInflation, na.rm = TRUE)
  )

# Print all rows of the summarized dataframe with mean and standard deviation
print(summary_pred_df, n = Inf)



### Plot the Marginal Effects ###########################################

# Calculate the marginal means
emmeans_result <- emmeans(zigam_21_exp1_int_adisPT_ziPT_ID_Rana, ~ Phase * Treatment, type = "response")
# Extract the marginal means and their confidence intervals
df_emm <- as.data.frame(summary(emmeans_result))
# Rename columns for clarity (optional)
colnames(df_emm) <- c("Phase", "Treatment", "Marginal_Mean", "SE", "df", "Lower_CI", "Upper_CI")

# Create the plot for marginal effects with faceting by Treatment
Plot_Marginal_Effects_Faceted <- ggplot(df_emm, aes(x = Phase, y = Marginal_Mean)) +
  geom_line(aes(group = Treatment, color = Treatment), size = 1) +  # Line for marginal means
  geom_point(size = 3, aes(color = Treatment)) +  # Points for marginal means
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI, color = Treatment), width = 0.2, size = 0.8) +  # Error bars
  facet_grid(rows = ~ Treatment, shrink = T, scales = "free_y", switch = "x", margins = F) +  # Facet by Treatment
  labs(title = "Marginal Effects of Treatment and Phase",
       y = "Predicted Active Seconds [s]",
       color = "Treatment") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),  # Adjust x-axis text
        legend.position = "none")  # Remove legend since treatments are faceted

# Print the plot
print(Plot_Marginal_Effects_Faceted)

# Print the plot
print(Plot_Marginal_Effects)



# ### Experiment 2 - Boxplot & Interaction + Zero-Inflation ### -----------

Box_exp2 <- ggplot(df_Rana_exp2, aes(x = Phase, y = Active_seconds, fill = Treatment)) +
  geom_boxplot(size = 0.1, alpha = 0) +  # Adjust alpha for overlay visibility
  facet_grid(rows = ~ Treatment, shrink = T, scales = "free_y", switch = "x", margins = F) +
  scale_y_continuous(limits = c(0, 150)) +
  labs(title = "",
       x = "",
       y = "Active Seconds \n [s]") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))  # Rotate x-axis labels to 0 degrees (horizontal)

# Extracting Values from Model predictions #
# Conditional model value predictions
pred_cond_exp2 <- predict(zigam_21_exp2_int_adisPT_ziPT_ID_Rana, type = "conditional", se.fit = TRUE)
# Zero-inflation probabilities value predictions
pred_zi_exp2 <- predict(zigam_21_exp2_int_adisPT_ziPT_ID_Rana, type = "zprob", se.fit = TRUE)

#Create a new dataframe for calculations and subsequent plotting
df_pred_exp2 <- data.frame(
  Phase = df_Rana_exp2$Phase,
  Treatment = df_Rana_exp2$Treatment,
  Active_seconds = pred_cond_exp2$fit,
  Active_seconds_se = pred_cond_exp2$se.fit,
  ZeroInflation = pred_zi_exp2$fit,
  ZeroInflation_se = pred_zi_exp2$se.fit)

df_pred_exp2$Phase <- factor(df_pred_exp2$Phase, ordered = TRUE)

df_pred_act_exp2 <- df_pred_exp2 %>% 
  select(Phase, Treatment, Active_seconds, Active_seconds_se)

df_pred_zi_exp2 <- df_pred_exp2 %>% 
  select(Phase, Treatment, ZeroInflation, ZeroInflation_se)

# Calculate means and standard errors for each combination of Treatment and Phase
df_calpred_act_exp2 <- df_pred_act_exp2 %>%
  group_by(Treatment, Phase) %>%
  summarise(Active_seconds = mean(Active_seconds, na.rm = TRUE),
            Active_seconds_se = mean(Active_seconds_se, na.rm = TRUE))

df_calpred_zi_exp2 <- df_pred_zi_exp2 %>%
  group_by(Treatment, Phase) %>%
  summarise(mean_zi_prob = mean(ZeroInflation, na.rm = TRUE),
            mean_zi_prob_se = mean(ZeroInflation_se, na.rm = TRUE))

# Combine Box_exp2 and Int_exp2 by overlaying them
Box_Int_combined <- Box_exp2 +
  geom_smooth(linewidth = 2, alpha = 2, method = "lm", se = FALSE, data = df_calpred_act_exp2, aes(x = Phase, y = Active_seconds, group = Treatment, color = Treatment)) +
  geom_point(size = 2, alpha = 2, data = df_calpred_act_exp2, aes(x = Phase, y = Active_seconds)) +
  geom_errorbar(linewidth = 1, width = 0.2, data = df_calpred_act_exp2, aes(x = Phase, ymin = Active_seconds - Active_seconds_se, ymax = Active_seconds + Active_seconds_se, color = Treatment)) +
  labs(title = "",
       y = "Active Seconds \n [s]",
       color = "Treatment",
       fill = "Treatment") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

# Plot for Zero Inflation Probabilities
Zi_exp2 <- ggplot(df_calpred_zi_exp2, aes(x = Phase, y = mean_zi_prob, group = Treatment, color = Treatment)) +
  geom_smooth(linewidth = 1, method = "lm", se = FALSE) +
  geom_point(size = 2) +
  geom_errorbar(linewidth = 1, width = 0.2, aes(ymin = mean_zi_prob - mean_zi_prob_se, ymax = mean_zi_prob + mean_zi_prob_se)) +
  facet_grid(rows = ~ Treatment, shrink = T, scales = "free_y", switch = "x", margins = F) +
  scale_y_continuous(limits = c(-0.05,1)) +
  labs(
    x = "Phase",
    y = "Propability of freezing \n [mean +SE]",
    title = ""
  ) +
  theme_bw() +
  theme(legend.position = "none")

# Combine plots using patchwork
combined_plot <- Box_Int_combined / Zi_exp2 + plot_layout(ncol = 1, heights = c(1, 1))

# Display combined plot
print(combined_plot)




# #### Abstellgleis #### --------------------------------------------------


-----------------------------------------------------------------------------------
### Oldschool Plots with three rows ###
### Boxplot - horizontal ###

Box_exp1 <- ggplot(df_Rana_exp1, aes(x = Phase, y = Active_seconds, fill = Treatment)) +
  geom_boxplot(size = 1) +
  facet_grid(rows = ~ Treatment, shrink = T, scales = "fixed", switch = "x", margins = F) +
  labs(title = "Active Seconds Across Phases and Treatments",
       x = "",
       y = "Active Seconds [s]") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))  # Rotate x-axis labels to 0 degrees (horizontal)


# Conditional model predictions
pred_cond_exp1    <- predict(zigam_21_exp1_int_adisPT_ziPT_ID_Rana, type = "conditional", se.fit = TRUE)
# Zero-inflation probabilities predictions
pred_zi_exp1      <- predict(zigam_21_exp1_int_adisPT_ziPT_ID_Rana, type = "zprob", se.fit = TRUE)


df_pred_exp1 <- data.frame(
  Phase = df_Rana_exp1$Phase,
  Treatment = df_Rana_exp1$Treatment,
  ActiveSeconds = pred_cond_exp1$fit,
  ActiveSeconds_se = pred_cond_exp1$se.fit,
  ZeroInflation = pred_zi_exp1$fit,
  ZeroInflation_se = pred_zi_exp1$se.fit)

df_pred_exp1$Phase <- factor(df_pred_exp1$Phase, ordered = TRUE)

df_pred_act_exp1 <- df_pred_exp1 %>% 
  select(Phase, Treatment, ActiveSeconds, ActiveSeconds_se)

df_pred_zi_exp1 <- df_pred_exp1 %>% 
  select(Phase, Treatment, ZeroInflation, ZeroInflation_se)

# Calculate means and standard errors for each combination of Treatment and Phase
# CALCULATED FROM MODEL PREDICTIONS!

df_calpred_act_exp1 <- df_pred_act_exp1 %>%
  group_by(Treatment, Phase) %>%
  summarise(mean_ActiveSeconds = mean(ActiveSeconds, na.rm = TRUE),
            mean_ActiveSeconds_se = mean(ActiveSeconds_se, na.rm = TRUE))

df_calpred_zi_exp1 <- df_pred_zi_exp1 %>%
  group_by(Treatment, Phase) %>%
  summarise(mean_zi_prob = mean(ZeroInflation, na.rm = TRUE),
            mean_zi_prob_se = mean(ZeroInflation_se, na.rm = TRUE))


# Interaction plot for Active Seconds #
Int_exp1 <- ggplot(df_calpred_act_exp1, aes(x = Phase, y = mean_ActiveSeconds, group = Treatment, color = Treatment)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_ActiveSeconds - mean_ActiveSeconds_se, ymax = mean_ActiveSeconds + mean_ActiveSeconds_se), width = 0.1) +
  facet_grid(rows = ~ Treatment, shrink = T, scales = "fixed", switch = "x", margins = F) +
  labs(title = "Effects of Treatment and Phase on Active Seconds",
       x = "",
       y = "Active Seconds [mean + SE]",
       color = "Treatment")
theme_bw() + 
  theme(legend.position = "none")


# Plot for Zero Inflation Probabilities #
Zi_exp1 <- ggplot(df_calpred_zi_exp1, aes(x = Phase, y = mean_zi_prob, group = Treatment, color = Treatment)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_zi_prob - mean_zi_prob_se, ymax = mean_zi_prob + mean_zi_prob_se), width = 0.1) +
  facet_grid(rows = ~ Treatment, shrink = T, scales = "fixed", switch = "x", margins = F) +
  labs(
    x = "Treatments",
    y = "Probability of Zero Activity",
    title = "Effects of Phase and Treatment on Zero Activity Probabilities") +
  theme_bw() +
  theme(legend.position = "none")

# Combine plots using patchwork
combined_plot <- Box_exp1 / Int_exp1 / Zi_exp1 + plot_layout(ncol = 1, heights = c(1, 1, 1))

# Display combined plot
combined_plot 


Box_exp2 <- ggplot(df_Rana_exp2, aes(x = Phase, y = Active_seconds, fill = Treatment)) +
  geom_boxplot(size = 1) +
  facet_grid(rows = ~ Treatment, shrink = T, scales = "fixed", switch = "x", margins = F) +
  labs(title = "Active Seconds Across Phases and Treatments",
       x = "",
       y = "Active Seconds [s]") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))  # Rotate x-axis labels to 0 degrees (horizontal)


# Conditional model predictions
pred_cond_exp2    <- predict(zigam_21_exp2_int_adisPT_ziPT_ID_Rana, type = "conditional", se.fit = TRUE)
# Zero-inflation probabilities predictions
pred_zi_exp2      <- predict(zigam_21_exp2_int_adisPT_ziPT_ID_Rana, type = "zprob", se.fit = TRUE)


df_pred_exp2 <- data.frame(
  Phase = df_Rana_exp2$Phase,
  Treatment = df_Rana_exp2$Treatment,
  ActiveSeconds = pred_cond_exp2$fit,
  ActiveSeconds_se = pred_cond_exp2$se.fit,
  ZeroInflation = pred_zi_exp2$fit,
  ZeroInflation_se = pred_zi_exp2$se.fit)

df_pred_exp2$Phase <- factor(df_pred_exp2$Phase, ordered = TRUE)

df_pred_act_exp2 <- df_pred_exp2 %>% 
  select(Phase, Treatment, ActiveSeconds, ActiveSeconds_se)

df_pred_zi_exp2 <- df_pred_exp2 %>% 
  select(Phase, Treatment, ZeroInflation, ZeroInflation_se)

# Calculate means and standard errors for each combination of Treatment and Phase
# CALCULATED FROM MODEL PREDICTIONS!

df_calpred_act_exp2 <- df_pred_act_exp2 %>%
  group_by(Treatment, Phase) %>%
  summarise(mean_ActiveSeconds = mean(ActiveSeconds, na.rm = TRUE),
            mean_ActiveSeconds_se = mean(ActiveSeconds_se, na.rm = TRUE))

df_calpred_zi_exp2 <- df_pred_zi_exp2 %>%
  group_by(Treatment, Phase) %>%
  summarise(mean_zi_prob = mean(ZeroInflation, na.rm = TRUE),
            mean_zi_prob_se = mean(ZeroInflation_se, na.rm = TRUE))


# Interaction plot for Active Seconds #
Int_exp2 <- ggplot(df_calpred_act_exp2, aes(x = Phase, y = mean_ActiveSeconds, group = Treatment, color = Treatment)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_ActiveSeconds - mean_ActiveSeconds_se, ymax = mean_ActiveSeconds + mean_ActiveSeconds_se), width = 0.1) +
  facet_grid(rows = ~ Treatment, shrink = T, scales = "fixed", switch = "x", margins = F) +
  labs(title = "Effects of Treatment and Phase on Active Seconds",
       x = "",
       y = "Active Seconds [mean + SE]",
       color = "Treatment") +
  theme_bw() + 
  theme(legend.position = "none")


# Plot for Zero Inflation Probabilities #
Zi_exp2 <- ggplot(df_calpred_zi_exp2, aes(x = Phase, y = mean_zi_prob, group = Treatment, color = Treatment)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_zi_prob - mean_zi_prob_se, ymax = mean_zi_prob + mean_zi_prob_se), width = 0.1) +
  facet_grid(rows = ~ Treatment, shrink = T, scales = "fixed", switch = "x", margins = F) +
  labs(
    x = "Treatments",
    y = "Probability of Zero Activity",
    title = "Effects of Phase and Treatment on Zero Activity Probabilities") +
  theme_bw() +
  theme(legend.position = "none")

# Combine plots using patchwork
combined_plot <- Box_exp2 / Int_exp2 / Zi_exp2 + plot_layout(ncol = 1, heights = c(1, 1, 1))

# Display combined plot
combined_plot


### Means and Standard Errors from Raw Data ###
# Calculate means and standard errors for each combination of Treatment and Phase
# CALCULATD FROM RAW DATA!
interaction_means_exp1 <- df_Rana_exp1 %>%
  group_by(Treatment, Phase) %>%
  summarise(mean_active_seconds = mean(Active_seconds, na.rm = TRUE),
            se_active_seconds = sd(Active_seconds, na.rm = TRUE) / sqrt(n()))
# Calculate means and standard errors for each combination of Treatment and Phase
interaction_means_exp2 <- df_Rana_exp2 %>%
  group_by(Treatment, Phase) %>%
  summarise(mean_active_seconds = mean(Active_seconds, na.rm = TRUE),
            se_active_seconds = sd(Active_seconds, na.rm = TRUE) / sqrt(n()))
# Note these values are different from the model predictions as the purely describe the data and do not predict

#---------------------------------------------------------------------------------




### Alternative to Boxplot & Interaction + Zero inflation ###
### Boxplot + Interaction & Zeroinflation ###

### Boxplot - Active Seconds ###
Box_exp1 <- ggplot(df_Rana_exp1, aes(x = Phase, y = Active_seconds, fill = Treatment)) +
  geom_boxplot(size = 1, alpha = 0.5) +
  facet_grid(rows = ~ Treatment, shrink = T, scales = "fixed", switch = "x", margins = F) +
  labs(title = "Active Seconds Across Phases and Treatments",
       x = "",
       y = "Active Seconds [s]") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))  # Rotate x-axis labels to 0 degrees (horizontal)

# Conditional model predictions
pred_cond_exp1 <- predict(zigam_21_exp1_int_adisPT_ziPT_ID_Rana, type = "conditional", se.fit = TRUE)
# Zero-inflation probabilities predictions
pred_zi_exp1 <- predict(zigam_21_exp1_int_adisPT_ziPT_ID_Rana, type = "zprob", se.fit = TRUE)

df_pred_exp1 <- data.frame(
  Phase = df_Rana_exp1$Phase,
  Treatment = df_Rana_exp1$Treatment,
  ActiveSeconds = pred_cond_exp1$fit,
  ActiveSeconds_se = pred_cond_exp1$se.fit,
  ZeroInflation = pred_zi_exp1$fit,
  ZeroInflation_se = pred_zi_exp1$se.fit
)

df_pred_exp1$Phase <- factor(df_pred_exp1$Phase, ordered = TRUE)

# Summarize data
df_calpred_exp1 <- df_pred_exp1 %>%
  group_by(Treatment, Phase) %>%
  summarise(mean_ActiveSeconds = mean(ActiveSeconds, na.rm = TRUE),
            mean_ActiveSeconds_se = mean(ActiveSeconds_se, na.rm = TRUE),
            mean_zi_prob = mean(ZeroInflation, na.rm = TRUE),
            mean_zi_prob_se = mean(ZeroInflation_se, na.rm = TRUE))

### Combined plot for Active Seconds and Zero Inflation Probabilities ###
combined_exp1 <- ggplot(df_calpred_exp1, aes(x = Phase, group = Treatment, color = Treatment)) +
  geom_line(aes(y = mean_ActiveSeconds), size = 1) +
  geom_point(aes(y = mean_ActiveSeconds), size = 2) +
  geom_errorbar(aes(ymin = mean_ActiveSeconds - mean_ActiveSeconds_se, ymax = mean_ActiveSeconds + mean_ActiveSeconds_se), width = 0.1) +
  geom_line(aes(y = mean_zi_prob * 100), size = 1, linetype = "dashed") + # Scale probability for second axis
  geom_point(aes(y = mean_zi_prob * 100), size = 2, shape = 4) +
  geom_errorbar(aes(ymin = (mean_zi_prob - mean_zi_prob_se) * 100, ymax = (mean_zi_prob + mean_zi_prob_se) * 100), width = 0.1, linetype = "dashed") +
  facet_grid(rows = ~ Treatment, shrink = T, scales = "fixed", switch = "x", margins = F) +
  labs(
    title = "Effects of Treatment and Phase on Active Seconds and Zero Activity Probabilities",
    x = "",
    y = "Active Seconds [mean + SE]",
    color = "Treatment"
  ) +
  scale_y_continuous(
    sec.axis = sec_axis(~./100, name = "Probability of Zero Activity [%]")  # Scale secondary axis back to original
  ) +
  theme_bw() +
  theme(legend.position = "none")

# Combine plots using patchwork
combined_plot <- Box_exp1 / combined_exp1 + plot_layout(heights = c(1, 1))

# Display combined plot
print(combined_plot)
