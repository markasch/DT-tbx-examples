# tune NN with caret
library(car)
library(caret)
set.seed(123)
trainIndex <- createDataPartition(Prestige$income, p=.7, list=F)
prestige.train <- Prestige[trainIndex, ]
prestige.test <- Prestige[-trainIndex, ]
# inspect the dataset
head(prestige.train)
# tune single-layer MLP:
#  size = number of neurons
#  decay = decay rate for weights
my.grid <- expand.grid(.decay = c(0.5, 0.1), .size = c(5, 6, 7))
prestige.fit <- train(income ~ prestige + education, 
                      data = prestige.train,
                      method = "nnet", 
                      maxit = 1000, 
                      tuneGrid = my.grid, 
                      trace = F, 
                      linout = 1)
print(prestige.fit)
# predict on test data
prestige.predict <- predict(prestige.fit, newdata = prestige.test)
prestige.rmse <- sqrt(mean((prestige.predict - prestige.test$income)^2)) 
##
## tune NN with 3 hidden layers
##
tune.grid.neuralnet <- expand.grid(
  layer1 = c(10,15),
  layer2 = c(5,10),
  layer3 = 5
)
model.neuralnet.caret <- train(
  income ~ prestige + education,
  data = prestige.train,
  method = "neuralnet",
  linear.output = TRUE, 
  tuneGrid = tune.grid.neuralnet, # cannot pass parameter hidden directly!!
  metric = "RMSE",
  trControl = trainControl(method = "cv")
  )
print(model.neuralnet.caret)
# predict on test data
neural.predict <- predict(model.neuralnet.caret, newdata = prestige.test)
neural.rmse <- sqrt(mean((neural.predict - prestige.test$income)^2)) 
# Compare the two RMSEs
paste("RMSE for nnet      =",prestige.rmse)
paste("RMSE for neuralnet =",neural.rmse)