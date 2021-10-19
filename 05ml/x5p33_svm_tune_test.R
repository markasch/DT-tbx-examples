#' ---
#' title: "Test tune options on svm"
#' author: 
#' date: 
#' output:
#'    pdf_document: 
#'       keep_tex: yes
#' ---
# tuning of SVM
library(e1071)      # for SVM and tune functions
set.seed(41)        # to make it reproducible
# load data
data(iris)
# tune 
obj <- tune.svm(Species~., data = iris, sampling = "fix",
                gamma = 2^c(-8,-6,-4,-2,0,2,4), cost = 2^c(-8,-6,-4,-2,-1,0,2,4,6,8))
best0 <- obj$best.model
table(Prediction =  predict(best0, iris), Truth = iris$Species)
cat(sprintf("prediction accuracy for best model = %.3f \n",
            sum(iris$Species ==  predict(best0, iris))/nrow(iris)))
summary(best0)
# a simpler example, tuned on parameter cost only 
# (much shorter CPU time)
set.seed(41)        # to make it reproducible
data(iris)
obj1 <- tune.svm(Species~., data = iris, 
                cost = 2^(2:8), 
                kernel = "linear") 
plot(obj1) # to show error evolution
best1 <- obj1$best.model
table(Prediction =  predict(best1, iris), Truth = iris$Species)
cat(sprintf("prediction accuracy for best model = %.3f \n",
            sum(iris$Species ==  predict(best1, iris))/nrow(iris)))
summary(best1)
# we can use best.svm directly...
set.seed(41)        # to make it reproducible
best2 <- best.svm(Species~., data = iris, 
                 cost = 2^(2:8), 
                 kernel = "linear") 
table(Prediction =  predict(best2, iris), Truth = iris$Species)
cat(sprintf("prediction accuracy for best model = %.3f \n",
            sum(iris$Species ==  predict(best2, iris))/nrow(iris)))
summary(best2)

