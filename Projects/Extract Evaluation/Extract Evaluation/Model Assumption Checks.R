### Model Assumption Checks ###
https://www.youtube.com/watch?v=EPIxQ5i5oxs

install.packages("performance")
library(performance) # https://rdrr.io/cran/performance/man/check_model.html # https://easystats.github.io/performance/reference/check_model.html #
library(lme4)
library(patchwork)


# Model Performance #
#no AIC & BIC for lmrob()-functions; AIC should be close/even to adjR2
model_performance(normal_model1)
model_performance(normal_model2) 
model_performance(robust_model1)
model_performance(robust_model2)
model_performance(robust_model1_KS2011)
model_performance(robust_model2_KS2011)
model_performance(robust_model1_KS2014)
model_performance(robust_model2_KS2014)

# MOdel Comparison #
# Models will be compared; should be of same/similar type
compare_performance(normal_model1, normal_model1_int, robust_model1, robust_model1_KS2011, robust_model1_KS2014, rank = T)
compare_performance(normal_model2, normal_model2_int, robust_model2, robust_model2_KS2011, robust_model2_KS2014, rank = T)

plot(compare_performance(normal_model1, robust_model1, robust_model1_KS2011, robust_model1_KS2014))
plot(compare_performance(normal_model2, robust_model2, robust_model2_KS2011, robust_model2_KS2014))
# in the plot points closer to the center indicate worse fit indices #

check_model(normal_model1)
check_model(normal_model1_int)

check_model(robust_model1)
check_model(robust_model1_KS2011)
check_model(robust_model1_KS2014)
check_predictions(normal_model1)


check_model(normal_model2)
check_model(normal_model2_int)

check_model(robust_model2)
check_model(robust_model2_KS2011)
check_model(robust_model2_KS2014)
check_predictions(normal_model2)

# Troubleshooting #
# 1. "Error: `check_model()` returned following error: input string 3 is invalid in this locale"
# helped get rid of this error and is somewhat linked to the locale settings
Sys.getlocale()
Sys.setlocale("LC_ALL", "C")

library(car)
vif(robust_model1)
par(mfrow = c(2, 2))
plot(robust_model1)

influence.measures(robust_model1)
