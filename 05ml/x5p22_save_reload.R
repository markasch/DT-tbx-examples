library(neuralnet)
set.seed(123)
# Split data
train_idx  <- sample(nrow(iris), 2/3 * nrow(iris))
iris_train <- iris[train_idx, ]
iris_test  <- iris[-train_idx, ]

# Binary classification
nn <- neuralnet(Species == "setosa" ~ Petal.Length + Petal.Width, 
                iris_train,
                linear.output = FALSE)
pred <- predict(nn, iris_test)
table(iris_test$Species == "setosa", pred[, 1] > 0.5)

# Multiclass classification
nn <- neuralnet((Species == "setosa") + (Species == "versicolor") 
           + (Species == "virginica") ~ Petal.Length + Petal.Width, 
            iris_train, 
            linear.output = FALSE)
pred <- predict(nn, iris_test)
table(iris_test$Species, apply(pred, 1, which.max))

# Save and reload
saveRDS(nn, "nn_model.rds")            # save the MLP model
save(iris_test, file = "data.RData")   # save the test data
# the next day...
rm(list = ls()) # start with a clean slate
saved_nn <- readRDS("nn_model.rds")    # load the NN model
load("data.RData")                     # load the test data
# perform inference on the new data
new_data <- iris_test
new_predict <- predict(saved_nn, new_data)
table(iris_test$Species, apply(new_predict, 1, which.max))