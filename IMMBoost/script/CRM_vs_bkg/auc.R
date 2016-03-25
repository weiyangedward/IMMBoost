# this script computes AUC given pred and label

suppressMessages(library(ROCR))

args <- commandArgs(trailingOnly=TRUE)

data <- read.table(args[1], sep=" ")
pred <- prediction( data[,2], data[,1])     
auc <- attr(performance(pred ,"auc"), "y.values")
write.table(auc,args[2])
