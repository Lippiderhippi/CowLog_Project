#Robust ANCOVA#

# Install the necessary package if you haven't already
install.packages("WRS2")

# Load the package
library(WRS2)
library(tidyverse)
library(dplyr)

df_W <- read.csv("C:/Users/Konrad Lipkowski/Desktop/GitHub/ExtractEvaluation/0. Extract Evaluation_All Rana_v3_W.csv")
#
df_L <- read.csv("C:/Users/Konrad Lipkowski/Desktop/GitHub/ExtractEvaluation/0. Extract Evaluation_All Rana_v3_L.csv")

# Inspect your data
head(df_L)
str(df_L)

# Omit missing values by creating a new dataframe from original df
df_noNA <- na.omit(df_W)

# Filtering data based on Filter_2 == 1 (Experiments 1)
# alternative ====> df_exp1 <- df_noNA %>%  filter(Filter_2 == 1) #"more" elegant
df_exp1 <- filter(df_noNA, Filter_2 == 1)     #"less" elegant
df_exp2 <- filter(df_noNA, Filter_2 == 2)     #"less" elegant

#Dont know if necessary/harmful but here I set "C1 = Control1" as the reference for subsequent analysis
# "$" is used to refere to access/refere to a specific column within df
df_exp1$Treat <- relevel(factor(df_exp1$Treat), ref = "C1")
df_exp2$Treat <- relevel(factor(df_exp2$Treat), ref = "C2")

# Perform robust ANCOVA
robust_ancova <- rmanova(Diff_PostPre ~ Treat + Pre, groups = df_exp1$Treat, data = df_exp1)
robust_ancova <- rmanova(Diff_PostPre ~ Treat + Pre, groups = df_exp2$Treat, data = df_exp2)

# Summarize the results
summary(robust_ancova)
