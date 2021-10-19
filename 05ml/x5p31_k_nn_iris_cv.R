# We consider the IRIS dataset - a classic in statistical learning
library(class)
data(iris)
head(iris)
summary(iris)
# Begin by separateing the data into 2 parts: 
#  - a training sample
#  - a test sampole.
train = iris[c(1:30,51:80,101:130),1:5]
test = iris[c(31:50,81:100,131:150),1:5]
# We want to determine the iris species from measurments of sepals and petals.
# We will use the knn() function to obtain the predictions.
pred = knn(train[,1:4], test[,1:4], train[,5], k = 2)
# Print the confusion matrix:
table(pred,test[,5])
# The number of correctly predicted values are on the diagonal.
# These results show that the k-NN classifier with k=2 leads to 
# a small classification error of only 3.33%.

# But, to choose k, it is recommended to perform cross-validation:
  
# 5-fold cross-validation for choosing k
# in the set {1,...,10}
fold = sample(rep(1:5,each=18)) # creation of groups B_v
cvpred = matrix(NA,nrow=90,ncol=10) # initialization of the matrix of predictors
for (k in 1:10)
  for (v in 1:5)
  {
    sample1 = train[which(fold!=v),1:4]
    sample2 = train[which(fold==v),1:4]
    class1 = train[which(fold!=v),5]
    cvpred[which(fold==v),k] = knn(sample1,sample2,class1,k=k)
  }
class = as.numeric(train[,5])
# Display the misclassification rates for k=1:10
apply(cvpred,2,function(x) sum(class!=x)) # calculate classif. error

# The best classification is obtained for k=6... 
# But this result then needs to be confronted with the knowledge of a domain expert 
# i.e. a botanist.
