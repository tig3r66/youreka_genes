# Deep learning transcriptomic model for prediction of pan-drug chemotherapeutic sensitivity

## Abstract
Emerging precision oncology studies have yet to generate a predictive biomarker that utilizes gene expression profiles to stratify tumours into similar pan-drug sensitivity profiles. This development would allow for the identification of candidate drugs for treatments that maximize therapeutic response and minimize cytotoxic burden. As such, this study utilized cell line sensitivity and molecular profiling data to generate a combinatorial gene expression predictive biomarker via feature selection and a deep learning model. A cohort of cell line gene expression data from Genomics of Drug Sensitivity in Cancer (GDSC) was clustered into two response groups. Cell line response groups showed a significant difference in pan-drug chemotherapeutic sensitivity. Due to the high dimensional nature of the microarray data, the Boruta feature selection algorithm was used to identify genes with the highest predictive value. The feature space was reduced to 300 genes; functional profiling was enriched primarily for the focal adhesion, ECM-receptor and proteoglycan interaction pathways. Using the selected genes, a deep learning neural network architecture was developed to predict response groups. It was determined that a 10 hidden layer neural network was optimal for the dataset; following hyperparameter tuning and 5-fold cross-validation, the model showed a predictive accuracy of 95.4\% (AUC = 0.913 $\pm$ 0.049).  This validates the postulate that cell lines with similar gene expression profiles present similar pan-drug chemotherapeutic sensitivity, and it suggests the potential utility of similar combinatorial biomarkers for the selection of potent candidate drugs.

## Source Code

Within [src/clustering](src/clustering/) is the code for the clustering analysis. Here is a brief description of important files:

- [gdsc2.r](src/clustering/gdsc2.r): code to create the heatmap, volcano plot, and KEGG functional annotation plot.
- [mutations.r](src/clustering/mutations.r): creating KRAS and TP53 PCA plots.
- [pca.r](src/clustering/pca.r): visualization of response group clusters. Also combines the heatmap, volcano plot, and response group cluster plot into one panelled figure.
- [tissueid.r](src/clustering/tissueid.r): PCA plot of tissue ID and tumour type plots.

Within [src/neural_net](src/neural_net/) is the code for feature selection, hyperparameter optimization, and cross-validation. Here is a brief description of important files:

- [boruta_trials.ipynb](src/neural_net/boruta_trials.ipynb): the Boruta algorithm to select statistically relevant genes.
- [hyperparameterized_model.ipynb](src/neural_net/hyperparameterized_model.ipynb): grid search for optimal hyperparameters for the neural network architecture.
- [cv.ipynb](src/neural_net/cv.ipynb): K-folds cross-validation for the neural network architectures analyzed.
- [layer_eval.ipynb](src/neural_net/layer_eval.ipynb): for model creation and quick investigation of parameters.
