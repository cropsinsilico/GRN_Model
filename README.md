# Gene regulatory network model

## Model Description

A Gene Regulatory Network (GRN) model was built using a linear regression algorithm to infer relationships between a dependent variable (in this case the expression of the putative target gene) and one or more independent variables (or predictors; in this case TFs). The resulting linear model was used to predict the response variable based on the states of the dependent variables. The regression algorithm was run in R using the LM function (Team, 2013), which optimizes variables of the linear model using a least squares fit between the response and dependent variables on training data (Bjorck, 1996). 
For example:
	mRNA_Target=g0+(W1*mRNA_TF1)+(W2*mRNA_TF2)+⋯ 	        
where, g0, W1, W2 are least squares optimized parameters for the linear model. mRNA_TF1, mRNA_TF2, etc. are expression values of Transcription factors (TFs) predicted to regulated target genes of interest in the static GRN. Parameters were optimized using training data, which ultimately resulted in a weight (Wx) that corresponds to the level of influence that a TF exerts on a predicted target gene’s expression. A linear model was generated for every gene in the static CO2-responsive GRN, and used to simulate the expression of genes of interest in both ambient and elevated CO2 environments in mature soybean leaf. 
In order to perform cross validation and check for overfitting, a predicted r-square value was calculated for every linear regression model (Quan, 1988). The predictive r-square helps determine how well the proposed regression model predicts responses for new observations. It is calculated by systematically removing each observation from the dataset, estimating the regression equation, and determining how well the model predicts the removed observations. Thus, it is a leave-one-out cross validation approach. The CO2-responsive dataset that is used to build the static GRN was used as a test dataset for the linear model, to predict expression of the target gene using the optimized weight associated with every TF. 

## References

Bjorck A. 1996. Numerical methods for least squares problems: Siam.
Team RC. 2013. R: A language and environment for statistical computing.
Quan NT. 1988. The prediction sum of squares as a general measure for regression diagnostics. Journal of Business & Economic Statistics, 6: 501-504.

## Running the model

```
R src/GRNModel.R
##Inputs used in the R script are in 'Input' folder
##Percentage change in expression levels due to TF perturbation generated by the R script provided in 'Output' folder. 
```

## Model Inputs/Outputs

### Inputs

Name | Units | Data Type | Description
---- | ----- | --------- | -----------
Gene regulatory network interactions | NA | NA | Interactions that denote equations for linear model (ExpGene1 = w1*ExpGene2 +w2*ExpGene3...)
Training data | NA | float | Normalized microarray mRNA expression values of soybean leaf tissues under different treatments for training linear model.

### Outputs

Name | Units | Data Type | Description
---- | ----- | --------- | -----------
Percentage difference in predicted gene expression | NA | float | Change in predicted gene expression due to a TF knockout or overexpression as compared to control.


