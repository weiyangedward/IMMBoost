# This script trains a randomforest model given
# training data and make prediction on test data
# Feature importance will also be computed.

suppressMessages(library(randomForest))
suppressMessages(library(ROCR))
args <- commandArgs(trailingOnly=TRUE)

##====== import train and test data ========
trainData <- read.table(args[1], header = T, sep="\t", row.names = 1)
testData <- read.table(args[3], header = T, sep="\t", row.names = 1)
trainLab <- read.table(args[2])
testLab <- read.table(args[4])

##====== parse train and test data ========
trainData  = data.matrix(trainData)
col = dim(trainData)[2]
trainLab = as.factor(trainLab[,1])
trainData = trainData[,1:col-1]

testData  = data.matrix(testData)
testLab = as.factor(testLab[,1])
testData = testData[,1:col-1]

##======= train and prediction =======
train.rf <- randomForest(trainData, trainLab, ntree = 1000, nodesize=1, importance = TRUE)
testPred = predict(train.rf, testData, type=c("prob"))

##======= output prediction =======##
labAndPred = cbind(testLab, testPred[,2])
labAndPred  = round(labAndPred, digits=3)
write.table(labAndPred,args[5],row.names=F, col.names=F)

##====== output importance ====##
featureImportance = round(importance(train.rf, type=1),3)
write.table(featureImportance,args[6],row.names=T, col.names=F)
