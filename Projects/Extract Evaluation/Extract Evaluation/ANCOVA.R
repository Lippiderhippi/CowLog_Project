# ANCOVA #

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

# Testing ANCOVA models #
normal_model1 <- lm(Diff_PostPre ~ Pre + Treat, data = df_exp1)
anova(normal_model1)
normal_model1_int <- lm(Diff_PostPre ~ Pre + Treat + Pre:Treat, data = df_exp1)
anova(normal_model1_int)

# General linear F-TEst #
anova(normal_model1, normal_model1_int)

# ===> If Model with Interaction is sig then use Interaction Term Model #