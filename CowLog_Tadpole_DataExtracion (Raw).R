if (!"devtools" %in% installed.packages()) install.packages("devtools"); library(devtools)

install_github("kellyjwallace/cowlogdata")


library(cowlogdata)
library(dplyr)
library(ggplot2)
library(viridis)
library(stringr)
library(broom)


clflag(pathtofile = "C:/Users/Lippi/Desktop/Daten-Extraktion")

summary_df = cldata(pathtofile = "C:/Users/Lippi/Desktop/Daten-Extraktion",
       outputdataname = "dataframe_round",
       outputzonename = "list_of_zones",
       factor = T, factorindex = 2, factorname = "Treat")

write.csv(summary_df, "summary_df.csv")
View(list_of_zones)
View(summary_df)
View(dataframe_round)


clseries(pathtofile = "C:/Users/Lippi/Desktop/Daten-Extraktion",
         zonename = list_of_zones,
         seglength = 30,
         factor = T, factorindex = 2, factorname = "Treat")


clpie(dataname = dataframe_round, zonename = list_of_zones, factor = F, factorname = "Treat")

clboxplot(dataname = dataframe_round, factor = T, factorname = "Treat")

clreg(data = dataframe_round, zonename = list_of_zones, factor = T, factorname = "Treat")

