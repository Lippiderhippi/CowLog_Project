if (!"devtools" %in% installed.packages()) install.packages("devtools"); library(devtools)

install_github("kellyjwallace/cowlogdata")

#### Cow Log script starts ####
#### Load the library packages first #### 

library(cowlogdata)
library(dplyr)
library(ggplot2)
library(viridis)
library(stringr)
library(broom)

#### then check data-tables for errors with "clflag" ####

clflag(pathtofile = "C:/Users/Konrad Lipkowski/Desktop/Cowlog-Data Extraction")

summary_df = cldata(pathtofile = "C:/Users/Konrad Lipkowski/Desktop/Cowlog-Data Extraction",
       outputdataname = "dataframe_round",
       outputzonename = "list_of_zones",
       factor = T, factorindex = 2, factorname = "Treat")

write.csv(summary_df, "summary_df.csv")

#### "clflag" will show tables that have invalid times and in which row of the table those occur ####
#### correct false times and check again with "clflag" ####

#### clseries will calculate the total time spend in each zone by subtracting time stamps from cow log ####
#### clseries will also show the total time spend in each zone per given time segment (seglength) ####

clseries(pathtofile = "C:/Users/Konrad Lipkowski/Desktop/Cowlog-Data Extraction",
         zonename = list_of_zones,
         seglength = 30,
         factor = T, factorindex = 2, factorname = "Treat")


clpie(dataname = dataframe_round, zonename = list_of_zones, factor = F)

clboxplot(dataname = dataframe_round, factor = T, factorname = "Treat")

clreg(data = dataframe_round, zonename = list_of_zones, factor = T, factorname = "Treat")

