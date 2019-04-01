# Neural_Networks
This work performs two different types of Neural Networks training and 5-fold CV, each one on three datasets that differ in three of the attributes. It uses 100 iterations with random seeds and after all iterations are completed, it calculates the average of different performance metrics


#Usage: 
The only things that needs to be changed are the input and output files directories and WEKA directory. Everything else is done automatically.

##Data:
This is a real data of polymorphisms and some of their conservation and biochemical properties. The names of the proteins and the polymorphisms were masked by indexes (2nd and 1st columns, respectively) and random noise was added to each of the columns.

#Goal:
Neural Networks will eliminate most of the columns. Each of two files has a common SubstitutionMatrix, SIFT and PolyPhen DIV columns (at least originally, before addition of noise). Columns called "entropy" and "pssm native" are different and come from computation of conservation patterns in homologous or orthologous sequences (non_congruent_known_homologs.arff and non_congruent_known_orthologs.arff, respectively). Both file contain "tag" column in which "1" corresponds to "Pathogenic Mutation" label and "0" to "Neutral Mutation"

The goal is to see whether using the conservation information based on orthologs will yield a better pathogenicity predictor than using information based on homologs. In each iteration, NN are training using:
1. Only SIFT and PolyPhen DIV columns
2. SIFT, PolyPhen DIV, substitutionmatrix,  entropy and pssm-native based on homologs sequences
3. SIFT, PolyPhen DIV, substitutionmatrix,  entropy and pssm-native based on orthologues sequences
(1-3 are repeated using NN with 0 hidden layers and 1 hidden layer and 2 nodes)

#Methods:
Training with only values of SIFT and PolyPhen (well-known pathogenicity predictors) is established as baseline. If next predictors will have lower performance than this one, the predictors are worthless, as they are less useful than just using the values of currently known predictors (the goal is to improve the pathogenicity predictions).

Because there are >3x more pathogenic variants than neutral ones, SMOTE is performed during training and neutral variants are oversampled. The exact way the SMOTE is performed is described in commented-out section in the beginning of the script.

This script doesn't include test phase, only training and 5-fold cross-validation. The performance is assessed directly from .rbf files generated during cross-validation, by calculating true and postitives/negatives from which (in the python script) sensitivities, MCC, specificities and accuracies are computed for each neural network condition and for each dataset. Because there are 100 of these values for each condition, average of each performance value is outputed. The reason the performance is calculated directly from .rbf instead of the usual training/test partition is because the dataset is too small and dividing it would reduce the training dataset even more. I didn't have an "outside" data to test my model so this is the reason the "Build Model" step was also skipped.

##Example of output and interpretation (it will vary because seeds are generated randomly):

homologs_siftpoly_ohl averages
MCC:
0.161
Sensitivity
0.517
Specificity
0.688
ACC
0.55

homologs_siftpoly_2hl averages
MCC:
0.1615
Sensitivity
0.6005
Specificity
0.604
ACC
0.6015

homologs_siftpoly_pssmnatentropy_ohl averages
MCC:
0.2545
Sensitivity
0.6425
Specificity
0.6775
ACC
0.6495

homologs_siftpoly_pssmnatentropy_2hl averages
MCC:
0.2455
Sensitivity
0.603
Specificity
0.709
ACC
0.6235

orthologs_siftpoly_pssmnatentropy_0hl averages
MCC:
0.204
Sensitivity
0.569
Specificity
0.6915
ACC
0.592

orthologs_siftpoly_pssmnatentropy_2hl averages
MCC:
0.2135
Sensitivity
0.71
Specificity
0.547
ACC
0.679


In both cases of 0 and 1 hidden layers, adding information about substitutionmatrix, entropy and pssm-nat improves the performance (except for specificities in some cases). Using ortholgs instead of homologs impoves overall sensitivity in case of NN with 1 hidden layer and 2 nodes. All other performance metrics are better when using homologs dataset instead of orthologs, in both 0 and 1 hidden layers conditions
