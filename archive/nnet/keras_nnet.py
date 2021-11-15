import sklearn
import os
import tensorflow
import numpy as np
import pandas as pd
import seaborn as sns
import itertools
import matplotlib.pyplot as plt
from os import walk
from tensorflow.keras.wrappers.scikit_learn import KerasClassifier
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import cross_val_score
from sklearn import datasets
from tensorflow.keras.utils import to_categorical
from tensorflow.keras.models import Sequential, load_model
from tensorflow.keras.layers import Dense,Dropout
from sklearn.model_selection import KFold
#from tensorflow.keras.models import load_model
from tensorflow.keras.callbacks import EarlyStopping
from sklearn.metrics import classification_report, confusion_matrix

#navigate to the directory and load the data into a panda dataframe
os.chdir(r'B:\Research\STEM Fellowship Hackathon\youreka_genes-master\src\nnet\data')
selected_genes = pd.read_csv('selected_genes.csv')
selected_genes = selected_genes.values.tolist()
selected_genes = list(itertools.chain(*selected_genes))
data_id = pd.read_csv('combined_expression.csv')
data_id1 = data_id.loc[:,selected_genes]

#data segregration and partition into the test and training 
X = data_id1
y = data_id['classification'].values
y = to_categorical(y)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20, random_state=42)

#normalize the expression levels
scaler = MinMaxScaler()
X_train = scaler.fit_transform(X_train)
X_test = scaler.transform(X_test)

#building the model

def baseline_model():
    model = Sequential()
        #adding hidden layers and dropout layers to prevent overfitting
    model.add(Dense(900,activation="relu"))
    model.add(Dropout(0.2))

    #model.add(Dense(450,activation="relu"))
    #model.add(Dropout(0.2))

    #model.add(Dense(450,activation="relu"))
    #model.add(Dropout(0.3))

    model.add(Dense(8,activation="softmax"))
    model.compile(optimizer='adam',loss='categorical_crossentropy',metrics=['accuracy'])
    return model

#compiling the model and taking measure to decrease overfitting with early_stopping and dropout neurons
    #early_stop = EarlyStopping(monitor='val_loss',mode='min',verbose=1, patience=50)
estimator = KerasClassifier(build_fn=baseline_model, epochs=100, batch_size=5, verbose=1)

#plotting the training and validation losses
model_loss = pd.DataFrame(model.history.history)
model_loss.plot()

#K-fold cross validation and model evaluation
kfold = KFold(n_splits=2, shuffle=True)
results = cross_val_score(estimator, X_train, y_train, cv=kfold)
print("Baseline: %.2f%% (%.2f%%)" % (results.mean()*100, results.std()*100))

#predicting based on the test set and visualizing the report on how well the model makes predictions
test_predictions = model.predict(X_test)
#print(classification_report(y_test,test_predictions)) #doesn't work yet
#print(confusion_matrix(y_test,test_predictions))  #doesn't work yet
model.save('classification_model_1.h5')
