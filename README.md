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

## Installing the dependencies

Issue the following commands from the R command prompt (make sure you run them one at a time):

```
install.packages("dplyr")
install.packages("data.table")
install.packages("igraph")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("BioNet")
```

## Running the model

### Without ``yggdrasil``

```
R src/GRNModel.R
##Inputs used in the R script are in 'Input' folder
##Percentage change in expression levels due to TF perturbation generated by the R script provided in 'Output' folder. 
```

### With ``yggdrasil``

```
yggrun GRNModel.yml GRNModel_tofile.yml
##Inputs used in the R script are in 'Input' folder
##Percentage change in expression levels due to TF perturbation generated by the R script provided in 'Output' folder. 
```

## Model Inputs/Outputs

### Inputs

Name | Units | Data Type | Description
---- | ----- | --------- | -----------
LeakeyCSoybeanMRRanks0.8Corr0.6.sif | NA | NA | Interactions that denote equations for linear model (ExpGene1 = w1*ExpGene2 +w2*ExpGene3...)
TrainingDataLeafSoy.txt | NA | float | Normalized microarray mRNA expression values of soybean leaf tissues under different treatments for training linear model.
Input_options.txt | NA | NA | A text file with four charaters specifying input options for GRN model as follows - 'Mutant' or 'WildType', Number of TFs perturbed if 'Mutant', list of TF IDs to be perturbed, factor by which TF expression needs to be increased (0, 1, 2 etc). Input options are given in that order with tab separation. List of TFs and list of factors are provided as comma separated values if more than 1 TF perturbed.

### Output

Name | Units | Data Type | Description
---- | ----- | --------- | -----------
bhlhbTFKoPercDiff.txt | NA | float | Percentage change in predicted gene expression due to a TF knockout or overexpression as compared to control.

## Additional notes regarding reproducibility of some tables and figures from GRN model output in Kannan et al., 2019 (doi: 'https://doi.org/10.1093/insilicoplants/diz008')

Figure 4 data provided as an OriginPro (https://www.originlab.com/) data analysis and graphing software file 'Figure4.opju' and uses 'Output/GRN_Ouput.txt' results

Figure 5 and Supplemental Figure S5 are a Gene Regulatory Networks (GRN) created in 'Cytoscape' (https://cytoscape.org/download.html) visualization tool using interactions/edges between TFs and their predicted Target genes in 'Input/LeakeyCSoybeanMRRanks0.8Corr0.6.sif'

Supplemental Figure S3 consists of linear regression fit between measured mRNA expression of Target genes and predicted mRNA expression of Target genes using the linear model with their adjusted R-squared value. Script used to generate all graphs in this figure are provided in 'Regression plots of target genes in GRN' section of 'src/GRNModel.R' R script.

Supplemental table S1 consists of predicted R-squared values and least squares linear regression optimized weights for individual 'Target' genes. 'outListDirected' and 'PredictedRSquaredTargets' R objects generated in the 'src/GRNModel.R' R script consists of the least squares optimized weights for individual TFs connected to the targets in the GRN and predicted R-squared values for the linear model of the particular target respectively.



