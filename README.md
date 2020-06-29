# Deep learning transcriptomic model for prediction of pan-drug chemotherapeutic sensitivity

## Source Code

Within [src/clustering](src/clustering/) is the code for the clustering analysis. Here is a brief description of important files:

- [gdsc2.r](src/clustering/gdsc2.r): code to create the heatmap, volcano plot, and KEGG functional annotation plot.
- [mutations.r](src/clustering/mutations.r): creating KRAS and TP53 PCA plots.
- [pca.r](src/clustering/pca.r): visualization of response group clusters. Also combines the heatmap, volcano plot, and response group cluster plot into one panelled figure.
- [tissueid.r](src/clustering/tissueid.r): PCA plot of tissue ID and tumour type plots.

Within [src/neural_net](src/neural_net/) is the code for feature selection, hyperparameter optimization, and cross-validation. Here is a brief description of important files in order of usage:

- [boruta_trials.ipynb](src/neural_net/boruta_trials.ipynb): the Boruta algorithm to select statistically relevant genes.
- [data_split.ipynb](src/neural_net/data_split.ipynb): splits data into X_train, X_test, y_train, and y_test datasets (80% training, 20% testing).
- [build_fns.py](src/neural_net/build_fns.py): build functions for neural networks with 1, 5, 10, and 15 hidden layers.
- [grid_params.ipynb](src/neural_net/grid_params.ipynb): grid search for optimal hyperparameters for the neural network architecture. **Warning:** 10 and 15 layer grid search is extremely computationally expensive. Ensure that you have at least 32 gb of RAM before running.
- [cv_stratified.ipynb](src/neural_net/cv_stratified.ipynb): K-folds stratified cross-validation for the neural network architectures analyzed.

## Neural Network Models
Within [src/neural_net/models](src/neural_net/models) are the trained Keras neural networks with 1, 5, 10, and 15 hiudden layers. To load them, use `model = keras.models.load_model('path/to/location')`, where `model` is the object you wish to load the model into.
