### GLM & GLMM ###
df_L <- read.csv("C:/Users/Konrad Lipkowski/Desktop/GitHub/ExtractEvaluation/0. Extract Evaluation_All Rana_v3_L.csv")

install.packages("lme4")
install.packages("lmerTest")
install.packages("emmeans")

library(lme4)
library(lmerTest)
library(emmeans)



# Assuming 'df' is your dataframe with columns 'Activity', 'Phase', 'Group', and 'Individual'

# Fit the GLMM
lmer_model <- lmer(Active_seconds ~ Phase * Treatment + (1 | Individual_Total), data = df_L, REML = FALSE)
glmer_model <- glmer(Active_seconds ~ Phase * Treatment + (1 | Individual_Total), data = df_L, family = tweedie(link = "log"))

# Check the model summary
summary(model)

# Obtain p-values using likelihood ratio test
anova(model)


install.packages("glmmTMB")
library(glmmTMB)

glmmTMB_model <- glmmTMB(Active_seconds ~ Phase * Treatment, data = df_L, family = ziGamma(link = "log"))
zigamma

class(df_L$Active_seconds)
anyNA(df_L$Active_seconds)
any(df_L$Active_seconds <= 0)

zi_gamma_model <- glmmTMB(Active_seconds ~ Phase + Treatment,
                          ziformula = ~ 1,  # Model for zero-inflation (intercept only)
                          data = df_L,
                          family = ziGamma(link = "log"))

zi_gamma_model_int <- glmmTMB(Active_seconds ~ Phase * Treatment,
                          ziformula = ~ 1,  # Model for zero-inflation (intercept only)
                          data = df_L,
                          family = ziGamma(link = "log"))
zi_gamma_model_zi <- glmmTMB(Active_seconds ~ Phase + Treatment,
                                 ziformula = ~ Phase + Treatment,  # Model for zero-inflation (intercept only)
                                 data = df_L,
                                 family = ziGamma(link = "log"))
zi_gamma_model_int_zi <- glmmTMB(Active_seconds ~ Phase * Treatment,
                              ziformula = ~ Phase + Treatment,  # Model for zero-inflation (intercept only)
                              data = df_L,
                              family = ziGamma(link = "log"))


# Model summary
summary(zi_gamma_model)
summary(zi_gamma_model_int)
summary(zi_gamma_model_zi)
summary(zi_gamma_model_int_zi)
