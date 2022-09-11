if (!"devtools" %in% installed.packages()) install.packages("devtools"); library(devtools)

install_github("kellyjwallace/cowlogdata")

library(cowlogdata)
library(dplyr)
library(ggplot2)
library(viridis)
library(stringr)
library(broom)

clflag(pathtofile = "E:/Balkan 2022/Alarm Cue/Data sheets")

cldata(pathtofile = "E:/Balkan 2022/Alarm Cue/Data sheets",
       outputdataname = "dataframe_round",
       outputzonename = "list_of_zones",
       factor = T, factorindex = 2, factorname = "Motu")

clseries(pathtofile = "E:/Balkan 2022/Alarm Cue/Data sheets",
         zonename = list_of_zones,
         seglength = 90,
         factor = T, factorindex = 2, factorname = "Motu")

clpie(dataname = dataframe_round, zonename = list_of_zones, factor = F)

clboxplot(dataname = dataframe_round, factor = T, factorname = "Motu")

clreg(data = dataframe_round, zonename = list_of_zones, factor = T, factorname = "Motu")
