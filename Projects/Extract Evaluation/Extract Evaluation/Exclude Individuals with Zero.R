# Exclude Rows with inactive Individuals in all Phases #

df_L <- read.csv("C:/Users/Konrad Lipkowski/Desktop/GitHub/ExtractEvaluation/0. Extract Evaluation_All Rana_v3_L.csv")

library(dplyr)
library(tidyverse)

filtered_df_L <- df_L %>%
  group_by(Individual_Total) %>%
  filter(!(all(Active_seconds[Phase == 1] == 0) & 
             all(Active_seconds[Phase == 2] == 0) & 
             all(Active_seconds[Phase == 3] == 0))) %>%
  ungroup()

# Filtering data based on Filter_2 == 1 (Experiments 1)
# alternative ====> df_exp1 <- df_noNA %>%  filter(Filter_2 == 1) #"more" elegant

df_exp1 <- filter(filtered_df_L, Filter_1 == 1)     #"less" elegant
df_exp2 <- filter(filtered_df_L, Filter_1 == 2)     #"less" elegant

#Do not know if necessary/harmful but here I set "C1 = Control1" as the reference for subsequent analysis
# "$" is used to refer to access/refer to a specific column within df
df_exp1$Treatment <- relevel(factor(df_exp1$Treatment), ref = "C1")
df_exp2$Treatment <- relevel(factor(df_exp2$Treatment), ref = "C2")





