library(data.table)
data = read.csv("dataset/2_walk/Perosn0101/hasc-130528-110507-acc.csv", header=F)
str(data) #     

colnames(data) = c("sec", "x", "y", "z") #

summary(data)

hist(diff(data$sec), col='blue',breaks=100);grid()

path = "dataset"
df = lapply(dir(path)[3:8], function(act){ 
  print(act)
  path = paste(path,act,sep="/")
  data = lapply(dir(path)[2:10], function(person){ 
    print(person)
    path = paste(path,person,sep="/")
    data = lapply(dir(path,pattern="-acc.csv")[1], function(file){ 
      print(file)
      path = paste(path, file, sep='/')
      data = read.csv(path, header = F)
      colnames(data) = c('sec', 'x', 'y', 'z')
      return(data.table(data, act, person, file))
    })
    return(rbindlist(data)) #data.table
  })
  return(rbindlist(data)) #data.table
})
df = rbindlist(df)
df$act = factor(df$act)
str(df)

plot(  rbind(data.frame(axis = 'x', value = df$x),
             data.frame(axis = 'y', value = df$y),
             data.frame(axis = 'z', value = df$z)
), col='orange'); grid()

acts = unique(df$act)
persons = unique(df$person) 
plot(df$x,df$y,col=df$act, cex=0.01);grid() 
legend('topleft', legend=acts, lwd=1, col=seq(acts))

calc_features = function(win){
  win = win[,.(x,y,z)] 
  ### intensity
  intensity = sqrt(win$x^2 + win$y^2 + win$z^2)
  
  ### mean/var
  means = colMeans(win, na.rm=T)
  vars = apply(win,2, function(x)var(x,na.rm=T))
  
  ### eigen covariance
  library(stats)
  eigens<-c(NA,NA)
  try({eigens <- eigen(cov(win))$values[1:2]},T)
  
  ### polar coordinates
  theta = with(win, sum(z/sqrt(x^2+y^2+z^2),na.rm=T))
  phy = with(win, sum(x/sqrt(x^2+y^2),na.rm=T))
  
  
  ### correlation between axis 
  mymean = function(x)mean(x, na.rm=T)
  mycov <- function(x, y, z){
    ret = mymean((x-mymean(x))*(y-mymean(y))/ (z-mymean(z))^2)
    if (is.na(ret)) ret = 0
    return(ret)
  }
  cor_x <- with(win, mycov(y, z, x))
  cor_y <- with(win, mycov(z, x, y))
  cor_z <- with(win, mycov(x, y, z))
  
  mydifcov <- function(x, y,z){
    mycov(diff(x),diff(y),diff(z))
  }
  difcov_x <- with(win, mydifcov(y,z,x))
  difcov_y <- with(win, mydifcov(z,x,y))
  difcov_z <- with(win, mydifcov(x,y,z))
  
  ### averaged frequency abs
  fft_x <- tryCatch(abs(fft(win$x))[-1],error=function()NA)
  fft_y <- tryCatch(abs(fft(win$y))[-1],error=function()NA)
  fft_z <- tryCatch(abs(fft(win$z))[-1],error=function()NA)
  
  fft_intensity <-  tryCatch(abs(fft(intensity))[-1],error=function()NA)
  
  mean_fft_x <- mymean(fft_x)
  mean_fft_y <- mymean(fft_y)
  mean_fft_z <- mymean(fft_z)
  mean_fft_intensity = mymean(fft_intensity)

  ### entropy of acceleration energy
  na.zero =function(x){if(is.na(x))0 else x}
  library(entropy, quietly=T)
  entropy_x <- na.zero(entropy(fft_x))
  entropy_y <- na.zero(entropy(fft_y))
  entropy_z <- na.zero(entropy(fft_z))
  entropy_intensity = na.omit(entropy(fft_intensity))
  
  
  ### original
  dev = intensity - mymean(intensity)
  
  calc_meancross = function(dev)sum(abs(diff(dev)>0))
  meancross <- calc_meancross(dev)
  
  outzero <-sum(abs(dev)>0.1)
  
  zonecross <- calc_meancross(dev[abs(dev)>0.1]) 
  
  return(
    data.table(eigen1 = eigens[1], eigen2= eigens[2], theta, phy, cor_x, cor_y, cor_z, 
               difcov_x, difcov_y, difcov_z, 
               mean_fft_x, mean_fft_y, mean_fft_z, mean_fft_intensity, 
               entropy_x, entropy_y, entropy_z, entropy_intensity, 
               meancross, outzero, zonecross))
}

twidth = 3000
tshift = 1500
system.time({
  feats = lapply(unique(df$file),function(file){ 
    print(file)
    ix = df$file == file
    data = df[ix,]
    act = data$act[1]; person = data$person[1]
    
    feats = lapply(seq(min(data$sec),max(data$sec), by = tshift), function(t){ 
      win = data[data$sec>=t & data$sec<t+twidth,]
      return(calc_features(win))
    })
    feats = rbindlist(feats)
    feats$act = act; feats$person = person
    return(feats)
  })
})

feats = rbindlist(feats)

featcols = colnames(feats)[-(22:23)]

save(feats, acts, persons, featcols, file = 'features.rdata')

load('features.rdata')
summary(feats)

tmp = lapply(featcols, function(col){
  hist(unlist(feats[,col,with=F]), col='blue', main=col)
})

tmp = lapply(featcols, function(col){
  plot(feats$act, unlist(feats[,col,with=F]), col='orange', main=col)
})

featcols = c('eigen1', 'mean_fft_intensity', 'outzero')

tmp = lapply(featcols, function(col){
  plot(feats$act, unlist(feats[,col,with=F]), col='orange', main=col)
})

do_eval = function(label, pre){
  if(length(label)!=length(pre))return()
  tab = table(pre, label)
  
  accuracy = sum(diag(tab))/sum(tab)
  cat("accuracy:", round(accuracy*100, 2), "%, ")
  
  precision = mean(diag(tab) / rowSums(tab))
  cat("precision:", round(precision*100, 2), "%, ")
  
  recall = mean(diag(tab / colSums(tab)))
  cat("recall:", round(recall*100, 2), "%, ")
  
  f_measure = 2*precision*recall / (precision + recall) 
  cat("F-measure", round(f_measure*100,2), "%\n")
  
  print(table(pre, label))
}

table(feats$person, feats$act)

train = feats[,c(featcols,"act"),with=F]
test = train

summary(train)

library(rpart)
cat('Tree:?n')
model = rpart(act~., train, method="class")
pre = predict(model, test, type="class")
do_eval(test$act, pre)

#NaiveBayes
cat('NaiveBayes:?n')
library(e1071)
model = naiveBayes(act~., train)
pre = predict(model, test)
do_eval(test$act,pre)


#SVM
cat('SVM:?n')
model = svm(act~., train)
pre = predict(model, test)
do_eval(test$act,pre)

#RandomForest
cat('RandomForest:?n')
library(randomForest)
model = randomForest(act~., train)
pre = predict(model, test)
do_eval(test$act,pre)

#KNN
cat('1NN:?n')
library(FNN)
pre = knn(train[,featcols,with=F],test[,featcols,with=F], train$act)
do_eval(test$act,pre)

#logistic
cat('logistic  :?n')
library(glmnet)
model = glmnet(as.matrix(train[,featcols, with=F]), train$act, family="multinomial")
pre = predict(model, as.matrix(test[,featcols,with=F]), type="class", s=0.01)
do_eval(test$act, pre)

#CV
table(feats$person, feats$act)

train = feats[feats$person %in% persons[1:5],c(featcols,"act"),with=F]
test = feats[feats$person %in% persons[6:10],c(featcols,"act"),with=F]


summary(train)
summary(test)

#1-person-leave-out CV

library(rpart)
res = rbindlist(lapply(persons, function(person){
  ix = feats$person == person
  train = feats[!ix, c(featcols,"act"),with=F]; test = feats[ix, c(featcols,"act"),with=F]
  model = rpart(act~., train, method="class")
  pre = predict(model, test, type="class")
  return(data.frame(act = test$act,pre))
}))
do_eval(res$act,res$pre)

#NaiveBayes
cat('NaiveBayes:?n')
library(e1071)
res = rbindlist(lapply(persons, function(person){
  ix = feats$person == person
  train = feats[!ix, c(featcols,"act"),with=F]; test = feats[ix, c(featcols,"act"),with=F]
  model = naiveBayes(act~., train)
  pre = predict(model, test)
  return(data.frame(act = test$act,pre))
}))
do_eval(res$act,res$pre)


#SVM
cat('SVM:?n')
res = rbindlist(lapply(persons, function(person){
  ix = feats$person == person
  train = feats[!ix, c(featcols,"act"),with=F]; test = feats[ix, c(featcols,"act"),with=F]
  model = svm(act~., train)
  pre = predict(model, test)
  return(data.frame(act = test$act,pre))
}))
do_eval(res$act,res$pre)



#RandomForest
cat('RandomForest:\n')
library(randomForest)
res = rbindlist(lapply(persons, function(person){
  ix = feats$person == person
  train = feats[!ix, c(featcols,"act"),with=F]; test = feats[ix, c(featcols,"act"),with=F]
  
  model = randomForest(act~., train)
  pre = predict(model, test)
  return(data.frame(act = test$act,pre))
}))
do_eval(res$act,res$pre)


#KNN
cat('1NN:?n')
library(FNN)
res = rbindlist(lapply(persons, function(person){
  ix = feats$person == person
  train = feats[!ix, c(featcols,"act"),with=F]; test = feats[ix, c(featcols,"act"),with=F]
  
  pre = knn(train[,featcols,with=F],test[,featcols,with=F], train$act)
  return(data.frame(act = test$act,pre))
}))
do_eval(res$act,res$pre)


#logistic
cat('logistic:?n')
library(glmnet)
res = rbindlist(lapply(persons, function(person){
  ix = feats$person == person
  train = feats[!ix, c(featcols,"act"),with=F]; test = feats[ix, c(featcols,"act"),with=F]
  
  model = glmnet(as.matrix(train[,featcols, with=F]), train$act, family="multinomial")
  pre = c(predict(model, as.matrix(test[,featcols,with=F]), type="class", s=0.01))
  return(data.frame(act = test$act, pre))
}))
do_eval(res$act,res$pre)

