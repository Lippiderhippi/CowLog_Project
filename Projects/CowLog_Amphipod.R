if (!"devtools" %in% installed.packages()) install.packages("devtools"); library(devtools)

install_github("kellyjwallace/cowlogdata") #htafjvv

#### Cow Log script starts ####
#### Load the library packages first #### 

library(cowlogdata)
library(dplyr)
library(ggplot2)
library(viridis)
library(stringr)
library(broom)

#### then check data-tables for errors with "clflag" ####


clflag(pathtofile = "E:/Balkan 2022/Alarm Cue/Data sheets")

summary_df = cldata(pathtofile = "E:/Balkan 2022/Alarm Cue/Data sheets",
       outputdataname = "dataframe_round",
       outputzonename = "list_of_zones",
       factor = T, factorindex = 2, factorname = "Motu")

write.csv(summary_df, "summary_df.csv")

#### "clflag" will show tables that have invalid times and in which row of the table those occur ####
#### correct false times and check again with "clflag" ####

#### clseries will calculate the total time spend in each zone by subtracting time stamps from cow log ####
#### clseries will also show the total time spend in each zone per given time segment (seglength) ####

clseries(pathtofile = "E:/Balkan 2022/Alarm Cue/Data sheets",
         zonename = list_of_zones,
         seglength = 90,
         factor = T, factorindex = 2, factorname = "Motu")


clpie(dataname = dataframe_round, zonename = list_of_zones, factor = F)

clboxplot(dataname = dataframe_round, factor = T, factorname = "Motu")

clreg(data = dataframe_round, zonename = list_of_zones, factor = T, factorname = "Motu")

