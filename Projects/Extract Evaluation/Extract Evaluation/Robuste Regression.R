### Robust regression is a technique that can reduce the impact of outliers, violation of the distribution assumption and heterogeneity in variance. It should be noted that the linearity assumption is still needed for proper inference using robust regression.

# Install the necessary packages if you haven't already
install.packages("robustbase")
install.packages("multcomp")
install.packages("olsrr")
install.packages("sandwitch")
install.packages("lmtest")

# Load the packages
library(dplyr) #to filter data based on a variable "filter()"
library(olsrr) # to plot outliers "ols_plot_resid_lev()"
library(lmtest) # Heteroscedasticity-Consistent Covariance Matrix Estimation
library(sandwich) # Heteroscedasticity-Consistent Covariance Matrix Estimation
library(robustbase) # robust regression "MM-type Estimators for Linear Regression"
library(multcomp) # postHoc multiple comparisons
library(tidyverse)


library(ggplot2)


df_W <- read.csv("C:/Users/Konrad Lipkowski/Desktop/GitHub/ExtractEvaluation/0. Extract Evaluation_All Rana_v3_W.csv")

# Omit missing values by creating a new dataframe from original df
df_noNA <- na.omit(df_W)

# Filtering data based on Filter_2 == 1 (Experiments 1)
# alternative ====> df_exp1 <- df_noNA %>%  filter(Filter_2 == 1) #"more" elegant

df_exp1 <- filter(df_noNA, Filter_2 == 1)     #"less" elegant
df_exp2 <- filter(df_noNA, Filter_2 == 2)     #"less" elegant

#Do not know if necessary/harmful but here I set "C1 = Control1" as the reference for subsequent analysis
# "$" is used to refer to access/refer to a specific column within df
df_exp1$Treat <- relevel(factor(df_exp1$Treat), ref = "C1")
df_exp2$Treat <- relevel(factor(df_exp2$Treat), ref = "C2")

#Normal LM#
normal_model1 <- lm(Diff_PostPre ~ Pre + Treat, data = df_exp1)
summary(normal_model1)
coeftest(normal_model1, vcov=vcovHC(normal_model1, type=c("HC3")))
plot(normal_model1, which = 1)
#
normal_model2 <- lm(Diff_PostPre ~ Pre + Treat, data = df_exp2)
summary(normal_model2)
coeftest(normal_model2, vcov=vcovHC(normal_model2, type=c("HC3")))
plot(normal_model2, which = 1)



# Outlier Plots #
ols_plot_resid_lev(normal_model1) #only works with 'lm'-function-model
ols_plot_resid_lev(normal_model2) #only works with 'lm'-function-model

# Robust Regression #
robust_model1 <- lmrob(Diff_PostPre ~ Pre + Treat, data = df_exp1)
summary(robust_model1)
par(mfrow = c(2, 2))
plot(robust_model1)
# Adjusting parameters in the lmrob function
robust_model1_KS2011 <- lmrob(Diff_PostPre ~ Pre + Treat, data = df_exp1, setting = "KS2011")
summary(robust_model1_KS2011)
robust_model1_KS2014 <- lmrob(Diff_PostPre ~ Pre + Treat, data = df_exp1, setting = "KS2014")
summary(robust_model1_KS2014)
#
robust_model2 <- lmrob(Diff_PostPre ~ Pre + Treat, data = df_exp2)
summary(robust_model2)
par(mfrow = c(2, 2))
plot(robust_model2)
# Adjusting parameters in the lmrob function 
robust_model2_KS2011 <- lmrob(Diff_PostPre ~ Pre + Treat, data = df_exp2, setting = "KS2011")
summary(robust_model2_KS2011)
robust_model2_KS2014 <- lmrob(Diff_PostPre ~ Pre + Treat, data = df_exp2, setting = "KS2014")
summary(robust_model2_KS2014)

#
robust_post_hoc1 <- glht(robust_model1, linfct = mcp(Treat = "Tukey"))
summary(robust_post_hoc1)
df_exp1$predicted <- predict(robust_model2, newdata = df_exp1)
ggplot(df_exp1, aes(x = Pre, y = Diff_PostPre, color = Treat)) +
  geom_point() +
  geom_smooth(method = "lmrob", formula = y ~ x, se = FALSE, aes(group = 1)) +
  labs(title = "Scatterplot with Robust Regression Lines",
       x = "Pre",
       y = "Diff_PostPre",
       color = "Treatment")
#
robust_post_hoc2 <- glht(robust_model2, linfct = mcp(Treat = "Tukey"))
summary(robust_post_hoc2)
df_exp2$predicted <- predict(robust_model2, newdata = df_exp2)
ggplot(df_exp2, aes(x = Pre, y = Diff_PostPre, color = Treat)) +
  geom_point() +
  geom_smooth(method = "lmrob", formula = y ~ x, se = FALSE, aes(group = 1)) +
  labs(title = "Scatterplot with Robust Regression Lines",
       x = "Pre",
       y = "Diff_PostPre",
       color = "Treatment")


par(mfrow = c(2, 2))
plot(robust_model2)

