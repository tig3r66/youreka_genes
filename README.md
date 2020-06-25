# Deep learning transcriptomic model for prediction of pan-drug chemotherapeutic sensitivity

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
