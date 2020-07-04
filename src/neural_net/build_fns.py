from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, BatchNormalization
from tensorflow.keras import regularizers


def hidden1(optimizer='sgd', init='normal', dropout=0.3, neurons=100, kernel_regularizer=1e-3):
    model = Sequential()
    # hidden layer
    model.add(Dense(neurons, activation='relu', kernel_regularizer=regularizers.l2(kernel_regularizer)))
    model.add(BatchNormalization())
    model.add(Dropout(dropout))
    # output layer
    model.add(Dense(1, activation='sigmoid'))
    # compiling
    model.compile(optimizer=optimizer, loss='binary_crossentropy', metrics=['accuracy'])
    return model


def hidden5(optimizer='sgd', init='normal', dropout=0.3, neurons=100, kernel_regularizer=1e-3):
    model = Sequential()
    # add 5 hidden layers
    for i in range(5):
        model.add(Dense(neurons, activation='relu', kernel_regularizer=regularizers.l2(kernel_regularizer)))
        model.add(BatchNormalization())
        model.add(Dropout(dropout))

    model.add(Dense(1, activation='sigmoid'))
    # compiling
    model.compile(optimizer=optimizer, loss='binary_crossentropy', metrics=['accuracy'])
    return model


def hidden10(optimizer='sgd', init='normal', dropout=0.3, neurons=100, kernel_regularizer=1e-3):
    model = Sequential()
    # add 10 hidden layers
    for i in range(10):
        model.add(Dense(neurons, activation='relu', kernel_regularizer=regularizers.l2(kernel_regularizer)))
        model.add(BatchNormalization())
        model.add(Dropout(dropout))

    model.add(Dense(1, activation='sigmoid'))
    # compiling
    model.compile(optimizer=optimizer, loss='binary_crossentropy', metrics=['accuracy'])
    return model


def hidden15(optimizer='sgd', init='normal', dropout=0.3, neurons=100, kernel_regularizer=1e-3):
    model = Sequential()
    # add 15 hidden layers
    for i in range(15):
        model.add(Dense(neurons, activation='relu',kernel_regularizer=regularizers.l2(kernel_regularizer)))
        model.add(BatchNormalization())
        model.add(Dropout(dropout))

    model.add(Dense(1, activation='sigmoid'))
    # compiling
    model.compile(optimizer=optimizer, loss='binary_crossentropy', metrics=['accuracy'])
    return model
