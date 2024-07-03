### Model Assumption Checks ###
https://www.youtube.com/watch?v=EPIxQ5i5oxs

install.packages("performance")
library(performance)
library(lme4)


model_performance(normal_model1)
model_performance(normal_model2)



model_performance(normal_model1)
model_performance(normal_model2)

compare_performance(normal_model1, robust_model1, rank = T)
plot(compare_performance(normal_model1, robust_model1))

r2(robust_model1)
r2(normal_model1)
AIC(robust_model1)
AIC(normal_model1)

check_model(normal_model1)
check_model(robust_model1)
