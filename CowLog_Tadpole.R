if (!"devtools" %in% installed.packages()) install.packages("devtools"); library(devtools)

install_github("kellyjwallace/cowlogdata")
install.packages("doBy")
install.packages("plyr")

#### Cow Log script starts ####
#### Load the library packages first #### 

library(cowlogdata)
library(plyr)       # Load plyr first ! #
library(dplyr)      # Load dplyr second ! #
library(ggplot2)
library(viridis)
library(stringr)
library(broom)

library(doBy)

######## "clflag" - Screening Data-tables  For Invalid Datapoints ########

    # "clflag" will scan data-tables for invalid datapoints and the respective row in which the invalid datapoint occures.
    # Example: [1] "18-EA_Trea(xRS)Pre.csv WARNING: LESS THAN THREE ENTRIES"

clflag(pathtofile = "C:/Users/Konrad Lipkowski/Desktop/Cowlog-Data Extraction")

# correct invalid datapoints and rerun "clflag" to check again if changes were applied. 
# Note: Depending on your usecase the error: "WARNING: LESS THAN THREE ENTRIES" can be ignored.  
# When everything is correct proceed with extracting activity-data.


######## Behaviour Data Extraction ########

summary_df = cldata(pathtofile = "C:/Users/Konrad Lipkowski/Desktop/Cowlog-Data Extraction",
       outputdataname = "dataframe_round",
       outputzonename = "list_of_zones",
       factor = T, factorindex = 2, factorname = "Treat")

write.csv(summary_df, "summary_df.csv")   # creates a data-table with calculated activity-data as a.csv.
View(summary_df)                          # to check if table and data is structured correctly.

########"clseries" - Visualizing toal time spend ########

clseries(pathtofile = "C:/Users/Konrad Lipkowski/Desktop/Cowlog-Data Extraction",
         zonename = list_of_zones,
         seglength = 30,
         factor = T, factorindex = 2, factorname = "Treat")

# clseries calculates and visualizes the total time spend in each "zone" or "status" by subtracting time stamps
# clseries visualizes time spend in 10 time segments with (x=seglenght) seconds each


clpie(dataname = dataframe_round, zonename = list_of_zones, factor = F)

clboxplot(dataname = dataframe_round, factor = T, factorname = "Treat")

clreg(data = dataframe_round, zonename = list_of_zones, factor = T, factorname = "Treat")



######## Simple Descriptive Across All Factors ########

# first loads your data into the working memory of R from your summary.csv into a dataframe.

df <- read.csv("summary_df.csv",stringsAsFactors = F,header=T) 

# to get a better understanding of the structure of your data-table and data itself run "dim()" and "str()"
      # "dim()" shows how much rows and colums your data-table has. Not necessary but nice to get an overview over the structure of you data-table
      # "str()" displays the structur of the object in R (not only for table but almost everything!)

dim(df)
str(df)

# There are three kinds of objects in data frame. i.e character (chr), integer (int) and number (num)

summary(df)   # summary of your data with means, medians etc. BUT across all factors!                                      
summary(df[c("Active_seconds","Resting_seconds")]) # summary of "specific variables" (columns) BUT across all factors!




######## Simple Descriptive Of One Variable by one specific Factor ########

cdata <- ddply(df, c("Treat"), summarise,
               N    = length(Active_seconds),
               mean = mean(Active_seconds, na.rm=T), 
               sd   = sd(Active_seconds, na.rm =T),
               se   = sd / sqrt(N))
cdata
write.csv(cdata, "desriptive.csv")





## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

summarySE(df, measurevar= "Active_seconds", groupvars=c("Treat"), na.rm=TRUE)
