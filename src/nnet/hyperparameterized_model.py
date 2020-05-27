import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# for data scaling and splitting
from sklearn.preprocessing import MinMaxScaler 
from sklearn.model_selection import train_test_split
# for neural net
from tensorflow.keras.models import Sequential, load_model
from tensorflow.keras.layers import Dense, Activation, Dropout
from tensorflow.keras.wrappers.scikit_learn import KerasClassifier
from tensorflow.keras.callbacks import EarlyStopping
# for evaluation
from sklearn.model_selection import KFold, cross_val_score, GridSearchCV
from sklearn.metrics import classification_report, confusion_matrix, roc_curve, auc

data = pd.read_csv("data/combined_expression.csv")
data.head()

selected_genes = pd.read_csv('cleaned/boruta-99-25-0.01.csv')
selected_genes = selected_genes.values.tolist()
selected_genes = list(itertools.chain(*selected_genes))

# retrieving proper columns
X = data.loc[:, selected_genes]
y = data['classification'].values

# scaling the data
scalar = MinMaxScaler()
x_scaled = scalar.fit_transform(X)

# splitting data (20% test, 80% train)
X_train, X_test, y_train, y_test = train_test_split(x_scaled, y, test_size=0.2, random_state=0)


# #### Gridsearch for Input and Output Layer (4 hidden layers)

# ## Optimizing Epochs and Batches

def create_model(optimizer='',init=''):
    model = Sequential()
    # adding layers and adding droplayers to avoid overfitting
    hidden_layers = len(selected_genes)
    model.add(Dense(hidden_layers, activation='relu'))
    model.add(Dropout(0.3))
    
    model.add(Dense((hidden_layers*0.5), activation='relu'))
    model.add(Dropout(0.3))
    
    model.add(Dense((hidden_layers*0.25), activation='relu'))
    model.add(Dropout(0.3))
    
    model.add(Dense((hidden_layers*0.25), activation='relu'))
    model.add(Dropout(0.3))
    
    model.add(Dense(1, activation='sigmoid'))
    # compiling
    model.compile(optimizer=optimizer, loss='binary_crossentropy', metrics=['accuracy'])
    return model


model = KerasClassifier(build_fn=create_model)
epochs = [50, 75, 100, 150]
batches = [16, 32, 64, 128]
optimizers = ['SGD', 'RMSprop', 'Adagrad', 'Adadelta', 'Adam', 'Adamax', 'Nadam']
init = ['glorot_uniform', 'normal', 'uniform']
param_grid = dict(epochs=epochs, batch_size=batches,optimizer=optimizers,init=init)
grid = GridSearchCV(estimator=model, param_grid=param_grid, cv=2, verbose=1, n_jobs=-1)
grid_result = grid.fit(X_train, y_train)

#Output the best parameters
print("Best: %f using %s" % (grid_result.best_score_, grid_result.best_params_))
grid_result.cv_results_

## Testing the Model

print(grid_result.best_params_['optimizer'])
model = KerasClassifier(build_fn=create_model, epochs=grid_result.best_params_['epochs'], batch_size=grid_result.best_params_['batch_size'],optimizer=grid_result.best_params_['optimizer'],init=grid_result.best_params_['init'])
kfold = KFold(n_splits=5, shuffle=True)
results = cross_val_score(model, X_train, y_train, cv=kfold)
print("Baseline Accuracy: %.2f%% (%.2f%%)" % (results.mean()*100, results.std()*100))


history = model.fit(x=X_train, y=y_train)
test_predictions = model.predict(X_test)

model_loss = pd.DataFrame(history.history)
model_loss.plot()

#visualize the model evaluation
print(classification_report(y_test, test_predictions))
print(confusion_matrix(y_test,test_predictions))

#save the model
model.model.save('model_1.h5')

