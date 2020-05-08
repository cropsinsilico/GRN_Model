#####GRN MODELING USING LINEAR REGRESSION Author: Kavya Kannan Email: kkannan2@illinois.edu####

library(dplyr)
library(data.table)
library(igraph)
library(BioNet)
library(grid)
library(gridExtra)
library(lattice)

##Deatch and then reattach yggdrasil/zeallot to prevent use of igraph %<-% operator
if (Sys.getenv("YGG_SUBPROCESS") != "") {
  detach("package:yggdrasil")
  detach("package:zeallot")
  library(yggdrasil)
  library(zeallot)
}

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

#####GRN MODELING USING LINEAR REGRESSION####
GRNModel <- function(LMData, OptionsFile) {
  # Convert lists to data frames
  if (Sys.getenv("YGG_SUBPROCESS") != "") {
    OptionsFile <- t(data.frame(unlist(OptionsFile)))
    }
if(OptionsFile[1,2] == 1){
  genes <- as.character(OptionsFile[1,3])
  genevals <- OptionsFile[1,4]
} else {
  genes <- unlist(strsplit(as.character(OptionsFile[1,3]), ",", fixed = TRUE))
  genevals <- unlist(strsplit(as.character(OptionsFile[1,4]), ",", fixed = TRUE)) 
}
 
  GraphEx <- loadNetwork.sif("Input/LeakeyCSoybeanMRRanks0.8Corr0.6.sif", format=c("igraph"), directed = TRUE) ##Load the network .sif file
  NodeList <- V(GraphEx)$name #List the network nodes
  NodeList <- NodeList[grepl("^G", NodeList)]##only picking genes coding for enzymes
  TrainingData <- LMData[c(1:213),]
  ### The for loop below generates a list having training data for every enzyme coding gene and basically subsets the training data such that   only columns for the enzyme coding gene and the TFs regulaing it are selected. 
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
  
  PredictedRSquaredTargets <- list(model_fit_stats(outListDirected[1]$Glyma.10G268500), model_fit_stats(outListDirected[2]$Glyma.19G022900), model_fit_stats(outListDirected[3]$Glyma.07G142700), model_fit_stats(outListDirected[4]$Glyma.08G165500), model_fit_stats(outListDirected[5]$Glyma.10G293500), model_fit_stats(outListDirected[6]$Glyma.16G168000), model_fit_stats(outListDirected[7]$Glyma.19G046800), model_fit_stats(outListDirected[8]$Glyma.11G226900), model_fit_stats(outListDirected[9]$Glyma.19G089100), model_fit_stats(outListDirected[10]$Glyma.17G015600, model_fit_stats(outListDirected[11]$Glyma.08G165400)))

  ##Regression plots of target genes in GRN
  LMDataWithoutSampleInfo <- LMData[,c(2:495)]
  LMDataPredictedAll <- LMTest(LMDataWithoutSampleInfo, outListDirected)
  colnames(LMDataPredictedAll) <- paste(colnames(LMDataPredictedAll), "1", sep = ".")
  AllDataTrainAndPredicted <- cbind(LMDataWithoutSampleInfo, LMDataPredictedAll)

  RegPlot1 = ggplot(AllDataTrainAndPredicted, aes(x = Glyma.10G268500, y = Glyma.10G268500.1)) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red")+labs(title = paste("Adj R2 = ",signif(summary(lm(AllDataTrainAndPredicted$Glyma.10G268500.1~AllDataTrainAndPredicted$Glyma.10G268500))$adj.r.squared, 5)), x = "Glyma.10G268500", y = "Glyma.10G268500_predicted")

  RegPlot2 = ggplot(AllDataTrainAndPredicted, aes(x = Glyma.19G022900, y = Glyma.19G022900.1)) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red")+labs(title = paste("Adj R2 = ",signif(summary(lm(AllDataTrainAndPredicted$Glyma.19G022900.1~AllDataTrainAndPredicted$Glyma.19G022900))$adj.r.squared, 5)), x = "Glyma.19G022900", y = "Glyma.19G022900_predicted")

  RegPlot3 = ggplot(AllDataTrainAndPredicted, aes(x = Glyma.07G142700, y = Glyma.07G142700.1)) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red")+labs(title = paste("Adj R2 = ",signif(summary(lm(AllDataTrainAndPredicted$Glyma.07G142700.1~AllDataTrainAndPredicted$Glyma.07G142700))$adj.r.squared, 5)), x = "Glyma.07G142700", y = "Glyma.07G142700_predicted")

  RegPlot4 = ggplot(AllDataTrainAndPredicted, aes(x = Glyma.08G165500, y = Glyma.08G165500.1)) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red")+labs(title = paste("Adj R2 = ",signif(summary(lm(AllDataTrainAndPredicted$Glyma.08G165500.1~AllDataTrainAndPredicted$Glyma.08G165500))$adj.r.squared, 5)), x = "Glyma.08G165500", y = "Glyma.08G165500_predicted")

  RegPlot5 = ggplot(AllDataTrainAndPredicted, aes(x = Glyma.10G293500, y = Glyma.10G293500.1)) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red")+labs(title = paste("Adj R2 = ",signif(summary(lm(AllDataTrainAndPredicted$Glyma.10G293500.1~AllDataTrainAndPredicted$Glyma.10G293500))$adj.r.squared, 5)), x = "Glyma.10G293500", y = "Glyma.10G293500_predicted")

  RegPlot6 = ggplot(AllDataTrainAndPredicted, aes(x = Glyma.16G168000, y = Glyma.16G168000.1)) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red")+labs(title = paste("Adj R2 = ",signif(summary(lm(AllDataTrainAndPredicted$Glyma.16G168000.1~AllDataTrainAndPredicted$Glyma.16G168000))$adj.r.squared, 5)), x = "Glyma.16G168000", y = "Glyma.16G168000_predicted")

  RegPlot7 = ggplot(AllDataTrainAndPredicted, aes(x = Glyma.19G046800, y = Glyma.19G046800.1)) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red")+labs(title = paste("Adj R2 = ",signif(summary(lm(AllDataTrainAndPredicted$Glyma.19G046800.1~AllDataTrainAndPredicted$Glyma.19G046800))$adj.r.squared, 5)), x = "Glyma.19G046800", y = "Glyma.19G046800_predicted")

  RegPlot8 = ggplot(AllDataTrainAndPredicted, aes(x = Glyma.11G226900, y = Glyma.11G226900.1)) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red")+labs(title = paste("Adj R2 = ",signif(summary(lm(AllDataTrainAndPredicted$Glyma.11G226900.1~AllDataTrainAndPredicted$Glyma.11G226900))$adj.r.squared, 5)), x = "Glyma.11G226900", y = "Glyma.11G226900_predicted")

  RegPlot9 = ggplot(AllDataTrainAndPredicted, aes(x = Glyma.19G089100, y = Glyma.19G089100.1)) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red")+labs(title = paste("Adj R2 = ",signif(summary(lm(AllDataTrainAndPredicted$Glyma.19G089100.1~AllDataTrainAndPredicted$Glyma.19G089100))$adj.r.squared, 5)), x = "Glyma.19G089100", y = "Glyma.19G089100_predicted")

  RegPlot10 = ggplot(AllDataTrainAndPredicted, aes(x = Glyma.17G015600, y = Glyma.17G015600.1)) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red")+labs(title = paste("Adj R2 = ",signif(summary(lm(AllDataTrainAndPredicted$Glyma.17G015600.1~AllDataTrainAndPredicted$Glyma.17G015600))$adj.r.squared, 5)), x = "Glyma.17G015600", y = "Glyma.17G015600_predicted")

  RegPlot11 = ggplot(AllDataTrainAndPredicted, aes(x = Glyma.08G165400, y = Glyma.08G165400.1)) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red")+labs(title = paste("Adj R2 = ",signif(summary(lm(AllDataTrainAndPredicted$Glyma.08G165400.1~AllDataTrainAndPredicted$Glyma.08G165400))$adj.r.squared, 5)), x = "Glyma.08G165400", y = "Glyma.08G165400_predicted")


  Figure4 <- grid.arrange(RegPlot1, RegPlot2, RegPlot3, RegPlot4, RegPlot5, RegPlot6, RegPlot7, RegPlot8, RegPlot9, RegPlot10, RegPlot11)

  ###InSilico Perturbation example for bHLH TF knockout###
  ### bHLHB TF Ko ###
  TestData <- LMData[c(214:229), c(2:495)]
  TestData[17,] <- colMeans(TestData[c(1:8),])
  TestData[18,] <- colMeans(TestData[c(9:16),])
  TestData <- TestData[-c(1:16),]
  TestDataStatic <- TestData
  TestDataControlLM <- LMTest(TestData, outListDirected)
  TFPerturbation <-  TestDataStatic
  if (as.character(OptionsFile[1,1]) == "Mutant") {
    for (i in 1:length(genes)){
    TFPerturbation[,as.character(genes[i])] <- TFPerturbation[,as.character(genes[i])]*as.numeric(genevals[i])
    }
    TFPerturbationLM <- LMTest(TFPerturbation, outListDirected)
    ExpPercDiff <- PercDiff(TFPerturbationLM,TestDataControlLM)
    return(ExpPercDiff)
  } else if (as.character(OptionsFile[1,1]) == "WildType") {
    ExpPercDiff <- PercDiff(TestDataControlLM,TestDataControlLM)
    return(ExpPercDiff)
  }

}
#### END OF GRN MODELING ####

if (Sys.getenv("YGG_SUBPROCESS") == "") {
  LMData <- read.table("Input/TrainingDataLeafSoy.txt", header = TRUE, sep = "\t") #Load the training data used for GRNmodel
  OptionsFile <- read.table("Input/Input_options.txt",nrow=1)
  results <- GRNModel(LMData, OptionsFile)
  write.table(results,"Output/GRN_Output.txt",quote=F,row.names = F, sep="\t")
}
