library(keras)

# Reading and dropping unnecessary column
data <- read.csv("data/combined_expression.csv", header=T)
data <- data.frame(data)
data <- subset(data, select = -c(CELL_LINE_NAME))


# Normalizing the input matrix
mat_data <- as.matrix(data)
dimnames(mat_data) <- NULL
mat_data[, 2:ncol(mat_data)] <- normalize(mat_data[, 2:ncol(mat_data)])
mat_data[, 1] <- as.numeric(data[, 1]) - 1


# partitioning (80:20 training-test split)
set.seed(1234)
ind <- sample(2, nrow(mat_data), replace=T, prob=c(0.7, 0.3))
## testing and training sets
train <- mat_data[ind==1, 2:ncol(mat_data)]
test <- mat_data[ind==2, 2:ncol(mat_data)]
## testing and training labels to one hot encoding
train_target <- mat_data[ind==1, 1]
test_target <- mat_data[ind==2, 1]
train_labels <- to_categorical(train_target)
test_labels <- to_categorical(test_target)


# creating sequential model
model <- keras_model_sequential()
model %>%
    layer_dense(units=10920, activation="relu", input_shape=c(16381)) %>%
    layer_dense(units=7, activation="softmax")
summary(model)
model %>%
    compile(loss="categorical_crossentropy",
            optimizer="adam",
            metrics="accuracy")
history <- model %>%
    fit(train,
        train_labels,
        epoch=5,
        batch_size=64,
        validation_split=0.2)


# evaluating
## evaluating accuracy
model2 <- model %>%
    evaluate(test, test_labels)
## creating confusion matrix
prob <- model %>%
    predict_proba(test)
pred <- model %>%
    predict_classes(test)
table <- table(Predicted=pred, Actual=test_target)
## results with prediction
results <- cbind(prob, pred, test_target)

