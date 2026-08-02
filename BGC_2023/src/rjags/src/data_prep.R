TrainingData <- read.table("../data/TrainingData.txt", header=F)

N <- nrow(TrainingData)

Y <- as.numeric(TrainingData[,1])
X <- as.matrix(TrainingData[,2:7], ncol=6, byrow=T)

#-----------------------------
#Add an intercept column to X
#-----------------------------
X <- cbind(rep(1, times=N), X)

#--------------------------
#Store number of covariates
#--------------------------
p <- ncol(X)

#-------------------
#Read the predictors
#-------------------
X.tilde <- as.matrix(read.table("../data/PredictionData.txt", header=F), ncol=6, byrow=T)

N.pred <- nrow(X.tilde)

#----------------------------------
#Add an intercept column to X.tilde
#----------------------------------
X.tilde <- cbind(rep(1, times=N.pred), X.tilde)

jags.data = list(N=N, p=p, N.pred=N.pred, Y=Y, X=X, X.tilde=X.tilde)


