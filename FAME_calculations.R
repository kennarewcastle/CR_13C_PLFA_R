# Calculations for the concentrations of FAME compounds in each sample.
# 18 April 2018
# KER


# Read in data frame that was created with append functions in FAME_data_prep.R
data<-read.csv("MASTER_DAT_W_SOIL.csv") 
data<-data[,2:9] # Get rid of extra numerical ID column

