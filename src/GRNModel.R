#####GRN MODELING USING LINEAR REGRESSION Author: Kavya Kannan Email: kkannan2@illinois.edu####

library(dplyr)
library(data.table)
library(igraph)
library(BioNet)

#####GRN MODELING USING LINEAR REGRESSION####

GraphEx <- loadNetwork.sif("Input/LeakeyCSoybeanMRRanks0.8Corr0.6.sif", format=c("igraph"), directed = TRUE) ##Load the network .sif file
NodeList <- V(GraphEx)$name #List the network nodes
NodeList <- NodeList[grepl("^G", NodeList)]##only picking genes coding for enzymes
LMData <- read.table("Input/TrainingDataLeafSoy.txt", header = TRUE, sep = "\t") #Load the training data used for GRNmodel
TrainingData <- LMData[c(1:213),]
### The for loop below generates a list having training data for every enzyme coding gene and basically subsets the training data such that only columns for the enzyme coding gene and the TFs regulaing it are selected. 
###So total of 11 datasets are subsetted this way for 11 enzyme coding genes
for (k in NodeList){
  NodeInDegList <- V(GraphEx)$name[neighbors(GraphEx, k, mode = "in")]
  if(length(NodeInDegList)!=0){
    NodeInDegList <- c(NodeInDegList, k)
    NodeTrainData <- TrainingData[, NodeInDegList]
    assign(paste("NodeTrainData", k, sep = "_"), NodeTrainData)
  }
}
##Creating a list of linear models for the 11 enzyme coding genes.
outlist=list() 
for (i in 1:length(NodeList)){
  f <- paste(colnames(get(paste("NodeTrainData", NodeList[i], sep = "_")))[ncol(get(paste("NodeTrainData", NodeList[i], sep = "_")))], "~", paste(colnames(get(paste("NodeTrainData", NodeList[i], sep = "_")))[1: ncol(get(paste("NodeTrainData", NodeList[i], sep = "_")))-1], collapse=" + "))
  outlist[[NodeList[i]]] <- lm(f, data = get(paste("NodeTrainData", NodeList[i], sep = "_")))
}
outListDirected <- outlist
###LMTest function predicts gene expression values for an enzyme coding gene based on linear model weights and test data TF expression values.
LMTest <- function(TestDataN, outlistN){
  predictions <- list()
  for (i in 1:length(outlistN)){
    tempList <- outlistN[[i]]$coefficients
    if(length(tempList) != 0){
      predictions[[names(outlistN[[i]]$model[1])]] <- predict(outlistN[[i]], TestDataN, se.fit = TRUE)
      TestDataN[,names(outlistN[[i]]$model[1])] <- predictions[[names(outlistN[[i]]$model[1])]]$fit
    }
  }
  return(TestDataN)
}


##PercDiff function calculates difference in expression values for TF perturbation results with two (WT and treatment specific) linear model prediction data as inputs
PercDiff <- function(TestData_m,TestData_n){
  PercentageChange <- c("Averages_Elevated","Average_ambient")
  for (i in names(TestData_m)){
    tempmax <- ((TestData_m[,i] - TestData_n[,i])*100/abs(TestData_n[,i]))
    tempmax <- data.frame(tempmax)
    colnames(tempmax)[1] <- i
    PercentageChange <- cbind(PercentageChange,tempmax)
  }
  return(PercentageChange)
}

###InSilico Perturbation example for bHLH TF knockout###
### bHLHB TF Ko ###
TestData <- LMData[c(214:229), c(2:495)]
TestData[17,] <- colMeans(TestData[c(1:8),])
TestData[18,] <- colMeans(TestData[c(9:16),])
TestData <- TestData[-c(1:16),]
TestDataStatic <- TestData
TestDataControlLM <- LMTest(TestData, outListDirected)
TFGlyma.18G115700Ko <-  TestDataStatic
TFGlyma.18G115700Ko$TFGlyma.18G115700 <- TFGlyma.18G115700Ko$TFGlyma.18G115700*0
Testbhlhb1Ko <- LMTest(TFGlyma.18G115700Ko, outListDirected)
bhlhbTFKoPercDiff <- PercDiff(Testbhlhb1Ko,TestDataControlLM)
write.table(bhlhbTFKoPercDiff,"Output/GRN_Output.txt",quote=F,row.names = F, sep="\t")


#### END OF GRN MODELING ####

##Crossvalidation and overfitting test for LM model results##

##Predicted Residual Sum of Squares calculation function 
PRESS <- function(linear.model) {
  #calculate the predictive residuals
  pr <- residuals(linear.model)/(1-lm.influence(linear.model)$hat)
  #calculate the PRESS
  PRESS <- sum(pr^2)
  
  return(PRESS)
}
##Function to calculate predicted R square using PRESS value and linear model result
pred_r_squared <- function(linear.model) {
  #Use anova() to get the sum of squares for the linear model
  lm.anova <- anova(linear.model)
  #Calculate the total sum of squares
  tss <- sum(lm.anova$'Sum Sq')
  #Calculate the predictive R^2
  pred.r.squared <- 1-PRESS(linear.model)/(tss)
  
  return(pred.r.squared)
}
##Function to display linear model related results
model_fit_stats <- function(linear.model) {
  r.sqr <- summary(linear.model)$r.squared
  adj.r.sqr <- summary(linear.model)$adj.r.squared
  ratio.adjr2.to.r2 <- (adj.r.sqr/r.sqr)
  pre.r.sqr <- pred_r_squared(linear.model)
  press <- PRESS(linear.model)
  return.df <- data.frame("R-squared" = r.sqr, "Adj R-squared" = adj.r.sqr, 
                          "Ratio Adj.R2 to R2" = ratio.adjr2.to.r2, "Pred R-squared" = pre.r.sqr, PRESS = press)
  return(round(return.df,3))
}
