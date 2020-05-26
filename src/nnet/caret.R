library(caret)
library(RANN)

train <- read.csv("data/combined_expression.csv")

# imputing missing values using KNN
id <- train$classification
preProcValues <- preProcess(train)
train_processed <- predict(preProcValues, train)
train_processed$CELL_LINE_NAME <- NULL
train_processed$classification <- id


# splitting data using caret
index <- createDataPartition(train_processed$classification, p=0.8, list=F)
trainSet <- train_processed[index,]
testSet <- train_processed[-index,]

# (recursive) feature selection w/ caret
control <- rfeControl(functions=rfFuncs, method="repeatedcv", repeats=5, verbose=T)
outcomeName <- "classification"
predictors <- names(trainSet)[!names(trainSet) %in% outcomeName]
Pred_Profile <- rfe(trainSet[, predictors], trainSet[, outcomeName], rfeControl=control)
Pred_Profile

# taking top 5 predictors
predictors <- c("Credit_History", "LoanAmount", "Loan_Amount_Term", "ApplicantIncome", "CoapplicantIncome")

model_nnet <- train(trainSet[, predictors], trainSet[, outcomeName], method="nnet")
fitControl <- trainControl(method="repeatedcv", number=5, repeats=5)
plot(model_nnet)
plot(varImp(object=model_nnet))

# predictions
predictions <- predict.train(object=model_nnet, testSet[,predictors], type="raw")
table(predictions)
confusionMatrix(predictions, testSet[,outcomeName])

# plotting the neural net
#library(NeuralNetTools)
#plotnet(model_nnet$finalModel)
