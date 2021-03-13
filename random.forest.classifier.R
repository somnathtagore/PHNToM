install.packages(c("Hmisc","caret","rpart","tree","e1071","ggplot2","readr","randomForest"))
install.packages("randomForestExplainer")
library(Hmisc)
library("caret")
library("rpart")
library("tree")
library("e1071")
library(ggplot2) # Data visualization
library(readr) # CSV file I/O, e.g. the read_csv function
library(randomForest)
library(randomForestExplainer)
set.seed(1)

# Input data files are available in the "../input/" directory.
# For example, running this (by clicking run or pressing Shift+Enter) will list the files in the input directory

#system("ls ../input")

# Read data
egfr.domains.data<-readRDS("~/Documents/viper-docs/benchmarking_1/signature_with_gtex/luad-gtex/egfr_luad/egfr.tf.cotf.top100.domain.wise.rds")
dim(egfr.domains.data)
pik3r1.train<-temp.plot.data.3.1.c.s.high.conf.lof.pi3kr1.ucec.analysis
dim(pik3r1.train)
colnames(pik3r1.train)
pik3r1.test<-temp.plot.data.3.1.c.s.unknown.mutations.pik3r1.ucec.samples
dim(pik3r1.test)
colnames(pik3r1.test)
colnames(egfr.domains.data)<-c("ECD.I","ECD.II","ECD.III","TM",
                               "Exon18","Exon19","Exon20","Exon21",
                               "CTT")
train <- read.csv("train.csv" , stringsAsFactors = FALSE)
test  <- read.csv("test.csv",   stringsAsFactors = FALSE)

egfr.domains.df<-as.data.frame(egfr.domains.data)
train<-cbind(egfr.domains.df$Exon19,egfr.domains.df$Exon21)
#train<-cbind(egfr.domains.df$Exon19)
#test<-cbind(egfr.domains.df$Exon20)
  test<-cbind(egfr.domains.df$ECD.I,egfr.domains.df$ECD.II,
              egfr.domains.df$ECD.III,egfr.domains.df$TM,
              egfr.domains.df$Exon18,egfr.domains.df$Exon20,
              egfr.domains.df$CTT)
train<-as.data.frame(train)
colnames(train)<-c("Exon19","Exon21")
#colnames(train)<-c("Exon19")
rownames(train)<-rownames(egfr.domains.df)
test<-as.data.frame(test)
colnames(test)<-c("ECD.I","ECD.II","ECD.III","TM","Exon18","Exon20","CTT")
#colnames(test)<-c("Exon20")
rownames(test)<-rownames(egfr.domains.df)

head(test)

str(train)

describe(train)

class(train)


#A function to extract features
# Missing data imputation
 # extractFeatures <- function(data) {
 #   features <- c("Pclass",
 #                 "Age",
 #                 "Sex",
 #                 "Parch",
 #                 "SibSp",
 #                 "Fare",
 #                 "Embarked")
 #   fea <- data[,features]
 #   fea$Age[is.na(fea$Age)] <- -1
 #   fea$Fare[is.na(fea$Fare)] <- median(fea$Fare, na.rm=TRUE)
 #   fea$Embarked[fea$Embarked==""] = "S"
 #   fea$Sex      <- as.factor(fea$Sex)
 #   fea$Embarked <- as.factor(fea$Embarked)
 #   #fea <- cbind(fea, fea$Age * fea$Age)
 #   return(fea)
 # }

#summary(extractFeatures(train))

#summary(extractFeatures(test))

#rf <- randomForest(extractFeatures(train), as.factor(train$Survived), ntree=100, importance=TRUE)
pik3r1.train<-as.data.frame(t(pik3r1.train))

train$group <- ifelse(train$V1>0,1,0)

rf <- randomForest(train, as.factor(train$group), ntree=1000, localImp=TRUE)
rf
rf <- randomForest(pik3r1.train, as.factor(pik3r1.train$TAF8), ntree=10000, localImp=TRUE)
rf
rf <- randomForest(formula=TAF8 ~ ., data=pik3r1.train, ntree= 10000, localImp = TRUE)
rf
rf$oob.times

# # create submission file
# submission <- data.frame( PassengerId= test$PassengerId )  # create a dataframe

# using model rf fit on training data to predict test data

test$group <- ifelse(test$V1>0,1,0)
rf.test <- predict( rf, test)  

pik3r1.test<-as.data.frame(t(pik3r1.test))
any(is.na(pik3r1.test))
rf.test <- predict( rf, pik3r1.test)
table(rf.test)
plot(rf.test)
class(rf.test)
rf$predicted

install.packages("pROC")
library(pROC)
rf.test.1 <- predict( rf, pik3r1.test,type="prob")
rf.test.1.roc <- multiclass.roc(pik3r1.test$HEXIM1,rf.test.1[,2])
plot(rf.test.1.roc)
rf.test.1.roc$response

auc(rf.test.1.roc)
class(rf.test.1)
rf.test.1

require(ROCR)
install.packages("party")
library(party)
x.cf <- cforest(HEXIM1 ~ ., data=pik3r1.test, control = cforest_unbiased(mtry = ncol(pik3r1.test)-2))
#x.cf.pred <-predict(x.cf,pik3r1.test)
x.cf.pred <- predict(x.cf, newdata=pik3r1.test)
#x.cf.prob <- 1-unlist(treeresponse(rf,pik3r1.test))
x.cf.prob <-  1- unlist(treeresponse(x.cf, pik3r1.test), use.names=F)[seq(1,nrow(pik3r1.test)*2,2)]

x.ct.prob.rocr <- prediction(x.cf.prob, pik3r1.test)
x.ct.perf <- performance(x.ct.prob.rocr, "tpr","fpr")
# add=TRUE draws on the existing chart 
plot(x.ct.perf, col=4, main="ROC curves of different machine learning classifier")
#roc_rf<-prediction(x.cf.pred,rf.test.1[,2])

# write results to CSV file
#write.csv(submission, file = "1_random_forest_r_submission.csv", row.names=FALSE)
# plot importance of preditors

min_depth<-min_depth_distribution(rf)
head(min_depth,n=10)
plot_min_depth_distribution(min_depth)

imp <- importance(rf, type=1 )
imp
print(row.names(imp))
featureImportance <- data.frame(Feature = row.names(imp), Importance = imp[, 1])
#plot_multi_way_importance(featureImportance,size_measure = "no_of_nodes",min_no_of_trees = 30)

# Use ggplot to plot the importance
p <- ggplot(featureImportance, aes(x= reorder(Feature, Importance) , y = Importance) ) +
  geom_bar(stat = "identity", fill = "#52cfff") +
  coord_flip() +
  theme_light(base_size = 20) +
  xlab("") + 
  ylab("Importance")+
  ggtitle("Random Forest Feature Importance\n") +
  theme(plot.title= element_text(size=18))

ggsave("2_feature_importance.png", p)
p

#various variable importance measures
importance_frame<-measure_importance(rf)
head(importance_frame,n=10)

#multi-way importance plot
plot_multi_way_importance(importance_frame,size_measure = "no_of_nodes",min_no_of_trees = 30)
plot_multi_way_importance(importance_frame,x_measure = "accuracy_decrease",y_measure = "gini_decrease",size_measure = "p_value")

explain_forest(rf,interactions = TRUE,pik3r1.train)
rf$predicted

#compare measures using ggpairs
plot_importance_ggpairs(importance_frame)

#compare different rankings
plot_importance_rankings(importance_frame)
# Classification Tree with rpart
library(rpart)

# grow tree 
fol= formula( as.factor(Survived) ~ Pclass + Age + Sex + Parch + SibSp + Fare + Embarked)
fit <- rpart( fol, data=train, method= "class")

print(fit)     #print results  
printcp(fit)   #display cp table  
plotcp(fit)    #plot cross-validation results 
#rsq.rpart(fit) #plot approximate R-squared and relative error for different splits (2 plots). labels are only appropriate for the "anova" method.  
summary(fit)   #detailed results including surrogate splits

# plot tree 
plot(fit, uniform=TRUE, main="Classification Tree")      #plot decision tree  
text(fit, use.n=FALSE, all=TRUE, cex=.8 )      #label the decision tree plot


####### new

library(caret)
## a vector of your selected features
fset = unique(as.vector(fmat[1:fi,]))
fset = unique(as.vector(ref.sig))
## convert to feature formula:
fset_formula = as.formula(paste("resp ~ ",paste(fset, collapse="+"),sep=""))
## train random forest, select
mfit <- train(fset_formula,
              method = "rf",
              trControl=trainControl(method="oob"),## optimal mtry parameters can be chosen by oob validation
              metric = "Kappa",
              tuneGrid=data.frame(mtry=floor(sqrt(length(fset)))), ## rf parameter: mtry, conventionally sqrt(length(features)) is used for classification
              data = data.frame(t(vp.train.bal),
                                resp=resp.train.bal)) ## data consists of  e.g. viper and response vector(classes)

require(randomForest)
install.packages("mlbench")
library(mlbench)
data(Sonar)
dim(Sonar)
Sonar[1:10,60:61]
str(Sonar[, 1:10])
write.csv(Sonar,file="Sonar.csv")
Sonar<-read.csv("plot.data.brca.pik3ca.top50.bottom50.csv",header=TRUE, sep=",")
Sonar<-read.csv("pik3ca.brca.gof.lof.neo.plot.data.gof.consensus.new.copy.csv",header=TRUE, sep=",")
rownames(Sonar) <- Sonar[,1]
#Sonar<-Sonar[,-1]
rownames(Sonar)
dim(Sonar)

mix.pik3ca.brca.cl.1.2.3.4.wt<-read.csv("~/Documents/random_forest/mix.pik3ca.brca.cl.1.2.3.4.wt.csv",header=TRUE, sep=",")
mix.pik3ca.brca.cl.1.2.3.4.wt<-mix.pik3ca.brca.cl.1.2.3.4.wt[,intersect(colnames(mix.pik3ca.brca.cl.1.2.3.4.wt),colnames(Sonar))]
#mix.pik3ca.brca.cl.1.2.3.4.wt<-mix.pik3ca.brca.cl.1.2.3.4.wt[unique(rownames(mix.pik3ca.brca.cl.1.2.3.4.wt)),]
rownames(mix.pik3ca.brca.cl.1.2.3.4.wt) <- mix.pik3ca.brca.cl.1.2.3.4.wt[,1]
mix.pik3ca.brca.cl.1.2.3.4.wt<-mix.pik3ca.brca.cl.1.2.3.4.wt[,-1]

Sonar<-Sonar[,intersect(colnames(Sonar),colnames(mix.pik3ca.brca.cl.1.2.3.4.wt))]
dim(Sonar)
# Sonar<-t(Sonar)
# Sonar<-as.data.frame(Sonar)
#Sonar[,41]
library(caret)
set.seed(998)
inTraining <- createDataPartition(Sonar$Class, p = .5, list = FALSE)
#inTraining <- createDataPartition(Sonar[,101], p = .5, list = FALSE)
training <- Sonar[ inTraining,]
testing  <- Sonar[-inTraining,]



testing <- mix.pik3ca.brca.cl.1.2.3.4.wt
fitControl <- trainControl(## 5-fold CV
  method = "cv",
  summaryFunction=twoClassSummary, 
  classProbs=T,
  savePredictions = T,
  number = 5)
#,
  ## repeated ten times
 # repeats = 5)

set.seed(825)
rfFit1 <- train(Class ~ ., data = testing, 
                 method = "rf", 
                 trControl = fitControl,
                 metric = "ROC",
                 verbose = TRUE)
rfFit1

fm <- rfFit1$finalModel


varImp(rfFit1)

#importance(fm)
# estimate variable importance
importance <- varImp(rfFit1, scale=FALSE)
# summarize importance
print(importance,20)
# plot importance
plot(importance,20)

plot(rfFit1)
ggplot(rfFit1)  

###svm
set.seed(123)
svmFit1 <- train(
  # diabetes ~., data = train.data, method = "svmLinear",
  # trControl = trainControl("cv", number = 10),
  # preProcess = c("center","scale")
  
  Class ~ ., data = testing, 
  method = "svmLinear", 
  trControl = fitControl,
  metric = "ROC",
  verbose = TRUE
)
importance <- varImp(svmFit1, scale=FALSE)
# summarize importance
print(importance)
# plot importance
plot(importance,20)
plot(svmFit1)

#knn
set.seed(123)
knnFit1 <- train(
  # diabetes ~., data = train.data, method = "svmLinear",
  # trControl = trainControl("cv", number = 10),
  # preProcess = c("center","scale")
  
  Class ~ ., data = testing, 
  method = "knn", 
  trControl = fitControl,
  metric = "ROC"
)
importance <- varImp(knnFit1, scale=FALSE)
# summarize importance
print(importance)
# plot importance
plot(importance,20)
plot(knnFit1)


#DT
set.seed(123)
dtFit1 <- train(
  # diabetes ~., data = train.data, method = "svmLinear",
  # trControl = trainControl("cv", number = 10),
  # preProcess = c("center","scale")
  
  Class ~ ., data = testing, 
  method = "rpart", 
  trControl = fitControl,
  metric = "ROC"
)
importance <- varImp(dtFit1, scale=FALSE)
# summarize importance
print(importance)
# plot importance
plot(importance,10)
plot(dtFit1)

# Make predictions on the test data
predicted.classes <- model %>% predict(test.data)
head(predicted.classes)


install.packages("MLeval")
library(MLeval)

predict.rf=predict(rfFit1,testing,type="prob")
res <- evalm(rfFit1)
res <- evalm(svmFit1)
res <- evalm(knnFit1)
res <- evalm(dtFit1)

######## saturation plot
#v <- read.csv(file="egfr.reg.improv.luad.best.90.percent.best.reg.avg.nes.csv",sep=",")
v <- c(0.741,0.752,0.763,0.772,0.779,0.781,0.784,0.791,0.793,
       0.80,0.802,0.803,0.806,0.807,0.808,0.8081,0.8082,0.809,
       0.81,0.82,0.83,0.835,0.837,0.84,0.841,0.845,0.847,0.85)
#v <- read.csv(file="avg.nes.plot.csv",sep=",")
#best.reg.1.avg.nes.copy
v<-as.data.frame(v)
v$n<-c(1:28)
#row.names(v) <- v[,1]
#v<-v[,-1]
# which(v$best.reg.1.avg.nes.copy!='')
# v$Gen<-c(1:length(which(v$best.reg.1.avg.nes.copy!='')))
# #plot(v[,1],type="o",col="blue",xlab="Generations",ylab="Best regulon avg.NES")
# v$best.reg.1.avg.nes.copy[v$best.reg.1.avg.nes.copy!='']
# v1<-as.data.frame(as.numeric(v[v!=""]))
# a<-length(v1$`as.numeric(v[v != ""])`)
# b<-v1$`as.numeric(v[v != ""])`
#ggplot(v1,aes(x = as.numeric(c(1:a)), y = b)) 
ggplot(v,aes(x = as.numeric(n), 
             y = as.numeric(v))) +
  #geom_hline(yintercept = 0, colour = "black") +
 # geom_point(fill = "#A4A4A4",colour = "darkred",  alpha = 1) + 
#  scale_shape_manual(values=c(17))+
  geom_point(shape = 21, colour = "black", fill = "white", size = 3, stroke = 3)+
 # geom_smooth()+
  #geom_point()+
  #geom_smooth(method=lm)
  scale_y_continuous(limits = c(0,1), breaks = seq(from = 0, to = 1, by = 0.1)) +theme_bw()+
  theme(axis.text = element_text(size=10,angle = -270)) +
  #geom_density_2d()
  ggtitle("Saturation Plot")
