
act <- dir("/Users/kawac/Documents/behavior_recognition/HAR_study")  #HAR_study内のファイル格納

##データを入れる配列作り##
resultlist <- data.frame()  #学習データ
lab <- c("stay","walk","jog","skip","stUP","stDown")
data <- data.frame()
l <- 1

##学習データ作成##
for (j in 3:8){
  setwd(paste("/Users/kawac/Documents/behavior_recognition/HAR_study",act[j],sep="/"))   #ディレクトリを１降りる
  per <- dir(".")
  poz <- getwd()
  for (k in 2:11){
    setwd(paste(poz,per[k],sep="/"))
    csvfile <- dir(".")
    cnt_csv <- grep("csv",csvfile)
    for (i in 1:length(cnt_csv)){
      data <- read.csv(csvfile[cnt_csv[i]],header=F)
      time <- as.numeric(row.names(data))/100
      data <- cbind(time,data)
      
      #特徴量の計算#
      features <- data.frame()
      windowTime <- 2
      sift <- 0.5
      
      startcsv <- data$time[1]
      finishcsv <- data$time[nrow(data)]
      wendcsv <- startcsv + windowTime
      
      while(wendcsv <= finishcsv){
        wincsv <- data[(data$time >= startcsv)&(data$time< wendcsv),]
        v <- sqrt(var(wincsv$V2)^2 + var(wincsv$V3)^2 + var(wincsv$V4)^2)
        m <- sqrt(mean(wincsv$V2)^2 + mean(wincsv$V3)^2 + mean(wincsv$V4)^2)
        
        startcsv <- startcsv + sift
        wendcsv <- startcsv + windowTime
        
        features <-rbind(features,data.frame("var"=v, "mean"=m,"act"=lab[l]))  #1つの学習データ完成
      }
      resultlist <- rbind(resultlist,features) #1人分のデータ
    }
  }
  cat(j)
  l <- l+1
}

#テストデータの特長量計算#
person0 <- read.csv("/Users/kawac/Documents/behavior_recognition/HAR_study/0_sequence/person0101/hasc-690101-120131-acc.csv",header = F)
label_person0 <- read.csv("/Users/kawac/Documents/behavior_recognition/HAR_study/0_sequence/person0101/hasc-690101-120131-acc.label",skip = 1, header = F)	
test <- data.frame()

for(i in 1:nrow(label_person0)){
  t_data <- person0[(person0$V1 >= label_person0[i,1])&(person0$V1 <= label_person0[i,2]),]
  time <- c(1:nrow(t_data)/100)
  t_data <- cbind(time,t_data)
  
  windowTime <- 2
  sift <- 0.5
  
  startcsv <- t_data$time[1]
  finishcsv <- t_data$time[nrow(t_data)]
  wendcsv <- startcsv + windowTime
  
  while(wendcsv <= finishcsv){
    wincsv <- t_data[(t_data$time >= startcsv)&(t_data$time< wendcsv),]
    v <- sqrt(var(wincsv$V2)^2 + var(wincsv$V3)^2 + var(wincsv$V4)^2)
    m <- sqrt(mean(wincsv$V2)^2 + mean(wincsv$V3)^2 + mean(wincsv$V4)^2)
    
    startcsv <- startcsv + sift
    wendcsv <-startcsv + windowTime
    
    test <-rbind(test,data.frame("var"=v, "mean"=m, "act"=as.character(label_person0[i,3])))
  }
}

#機械学習#
L <- length(resultlist[,1])
result <- data.frame()

#SVM#
library(e1071)
model_act <- svm(act~var+mean,resultlist)
pre_act <- predict(model_act,test)
result <- rbind(result,data.frame(test$act,pre_act))

#精度#
result$test.act <- factor(as.character(result$test.act))　#ラベルの並び替え#
result$pre_act <- factor(as.character(result$pre_act))
tab <- table(result$pre_act, result$test.act)
print(tab)
acc <- sum(diag(tab))/sum(tab)
cat("精度：",acc*100,"%")

levels(result$test.act)
levels(result$pre_act)
str(result)
table(result)
tab <- table(result)
sum(diag(tab))
sum(diag(tab))/ sum(tab) *100
      
      
    


