#read in the data file
LCA_variables <- read.csv2("")

library(psych)

# identifying classes
library(mclust)
TSM <- data.frame(LCA_variables$LG1x1,LCA_variables$LG1y1,LCA_variables$LG1z1) #LG1x1, LG1y1, and LG1z1 are example variable names.
TSM_lc<-Mclust(TSM,G=1:5)

summary(TSM_lc)
show(TSM_lc$BIC)
summary(TSM_lc$BIC)
plot(TSM_lc$BIC)

# plotting the classes
clPairs(TSM,TSM_lc$class)
table(TSM_lc$class)

#adding the class assignment for each subject to the data file
#library(psych)
TSM_class <-data.frame(LCA_variables$ID,TSM_lc$class, TSM_lc$z)
TSM_class$ID <- TSM_class$LCA_variables.ID #renaming variable
TSM_class$LCA_variables.ID <- NULL #deleting redundant
LCA_variables_wclass <- merge(LCA_variables,TSM_class, by="ID")
describeBy(LCA_variables_wclass$LG1x1,LCA_variables_wclass$TSM_lc.class)
describeBy(LCA_variables_wclass$LG1y1,LCA_variables_wclass$TSM_lc.class)
describeBy(LCA_variables_wclass$LG1z1,LCA_variables_wclass$TSM_lc.class)

#save data frame to excel
library(writexl)
write_xlsx(LCA_variables_wclass,".xlsx")