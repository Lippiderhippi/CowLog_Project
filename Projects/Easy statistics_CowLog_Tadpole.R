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


clflag(pathtofile = "C:/Users/Lippi/Desktop/Daten-Extraktion")

summary_df = cldata(pathtofile = "C:/Users/Lippi/Desktop/Daten-Extraktion",
       outputdataname = "dataframe_round",
       outputzonename = "list_of_zones",
       factor = T, factorindex = 2, factorname = "Treat")

write.csv(summary_df, "summary_df.csv")
View(summary_df)
View(list_of_zones)
View(dataframe_round)


#### "clflag" will show tables that have invalid times and in which row of the table those occur ####
#### correct false times and check again with "clflag" ####

#### clseries will calculate the total time spend in each zone by subtracting time stamps from cow log ####
#### clseries will also show the total time spend in each zone per given time segment (seglength) ####


clseries(pathtofile = "C:/Users/Lippi/Desktop/Daten-Extraktion",
         zonename = list_of_zones,
         seglength = 30,
         factor = T, factorindex = 2, factorname = "Treat")


clpie(dataname = dataframe_round, zonename = list_of_zones, factor = F, factorname = "Treat")

clboxplot(dataname = dataframe_round, factor = T, factorname = "Treat")

clreg(data = dataframe_round, zonename = list_of_zones, factor = T, factorname = "Treat")


#### Einfache Statistic ####

# Installiere das Paket "psych" falls noch nicht installiert
install.packages("psych")

# Lade das Paket "psych"
library(psych)

# Lade die CSV-Datei in ein dataframe. Die Datei muss im "CowLog R-Project" Ordner sein.
df <- read.csv("RT-I-summary_df_v2.csv")
View(df)

# Berechne die Statistiken nach einer Gruppe (hier am Beispiel der Spalte "Gruppe")
stats_by_group <- describeBy(df[, c("Active_seconds", "Active_entries", "Active_prop")], group = df$Treat, mat = TRUE)

# Gib die berechneten Statistiken aus
stats_by_group
View(stats_by_group)
write.csv(stats_by_group, file = "output.csv", row.names = TRUE)

# Subtrahiere die Aktivitätswerte zwischen Gruppen
diff_by_group <- with(df, tapply(Active_seconds, Treat, function(x) x - x[1]))

# Output anzeigen
print(diff_by_group)
View(diff_by_group)


