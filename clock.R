#1. load libraries
library(dplyr)
library(readr)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(glmnet)
library(purrr)
library(randomForest)
library(optRF)
library(Boruta)
library(caret)
#install.packages("randomForest")
install.packages("optRF")
install.packages("Boruta")
install.packages("caret")

setwd("/home/danielruzzante/BrookTroutEpigenomics2024/methyldata")

#2. read in table of file names
methyl_list <- list.files(path = '/home/danielruzzante/BrookTroutEpigenomics2024/methyldata')

#3. read in the first file to begin it 
allmethyl <- read.table(methyl_list[1], header=TRUE, stringsAsFactors = FALSE, sep = "")
label2 <- methyl_list[1]
label3 <- substr(label2, 1, 13)
allmethyl$Individual <- label3[1] #makes a column with the individual name 

#4. Loop that joins a temporary file to the larger file 
for (i in 2:length(methyl_list)) {  ## a) starting on the 2nd file, read through every file on methyl_list
  tmp <- read.table (methyl_list[i], header=TRUE, stringsAsFactors = FALSE, sep = "")
  label <- methyl_list[i] 
  # b) select different cut points depending on the label (one is 2 characters, the others are 3 which messed things up)  
  if(grepl("CV", label)) {
    cut_point <- 13
  } else if (grepl("RCD", label)) {
    cut_point <- 14
  } else if (grepl("WWD", label)) {
    cut_point <- 14
  }
  
  label1 <- substr(label, 1, cut_point)  #c) cuts the label in the predetermined cutpoint
  tmp$Individual <- label1[1]  #d) makes a new column with the individual name
  allmethyl <- rbind(allmethyl, tmp) ##e) joins each file to the allmethyl dataframe
}

##5. filter allmethyls columns
longdf<- allmethyl  %>% filter(strand == "F") %>% na.omit()
longdf <- longdf[,c("chrBase","freqC","Individual")]

#df withonly 6 cpg sites seelcted by boruta
#with data set from kselect and bpruta selection

df <- subset(longdf, (longdf$chrBase %in% c('NC_074666.1.97009861', 'NC_074684.1.43247228', 'NC_074687.1.653304', 'NC_074689.1.41397363', 'NW_026600289.1.186696', 'NW_026602756.1.2848' )))


borutadf2<- subset(longdf, (longdf$chrBase %in% c('NC_074666.1.97004172', 'NC_074666.1.97009861', 'NC_074667.1.67273212', 'NC_074675.1.58512044', 'NC_074675.1.58512055',
                                                  'NC_074680.1.6219926', 'NC_074684.1.43247228', 'NC_074685.1.49691716', 'NC_074687.1.653304', 'NC_074689.1.41397363',
                                                  'NW_026600289.1.186696', 'NW_026601214.1.40052', 'NW_026602756.1.2848')))

borutadf<- subset(longdf, (longdf$chrBase %in% c('NC_074666.1.68036787', 'NC_074666.1.97000648', 'NC_074666.1.97275220', 'NC_074675.1.58512044',
                                                  'NC_074680.1.6219926', 'NC_074689.1.41397363', 'NW_026602756.1.2848' )))

streamdf <- subset(longdf, (longdf$chrBase %in% c( 'NC_074666.1.68022932', 'NC_074666.1.68035386', 'NC_074672.1.8566415', 'NC_074698.1.39784365',
      
                                                                                               'NW_026600879.1.102880', 'NC_074666.1.67608631',  'NC_074672.1.8566415', 'NC_074685.1.49804959', 'NC_074694.1.42670735'))) 

borutadf <- subset(longdf, (longdf$chrBase %in% c('NC_074666.1.97009861', 'NC_074687.1.653304', 'NC_074689.1.41397363', 'NW_026600289.1.186696', 'NW_026602756.1.2848', 'NC_074684.1.43247228')))

#longdfr <- subset(longdf, !(longdf$Individual %in% c('23_CV_SFO_005', '23_WWD_SFO_016')))
mergeddf <- merge(longdf, agedf) #combine cpg and age info into one dataset
#6. make a matrix for y with age 
age <- brooktrout...Sheet1 %>% select(Individual, clock_age) 
age <- brooktrout...Sheet1 
agedf <- age[,c("Individual", "clock_age")]
write.csv(agedf, file = "age_data")

# only two columns
#7. make a matrix for x spreading out the CpG sites using freqC as values
wide <- borutadf2 %>% pivot_wider(names_from = chrBase, values_from = freqC) 

#8. Remove NA valus without breaking the datset
sitesna <- colSums(is.na(wide))
#Histogram of Na number per CpG site
colSums(is.na(widefilt)) %>% hist( main = paste("Histogram of raw CpG's"), col = "lightblue", ylab = "number of NA")

number <- c(26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1)
count_row <- sapply(wide, function(x) sum(is.na(x))) #counts number of na values
cpg <- names(count_row[count_row %in% number]) #identifies rows with a specific number of NA values
widefilt <- wide[, !(names(wide) %in% cpg)] #removes cpg column in wide dataset with he requested number of NA values

#mergeddf <- merge(widefilt, agedf) #combine cpg and age info into one dataset

#removes individual, NA, and turns it into a matrix
clean_kselect <- widefilt %>%
  select(-Individual) %>% na.omit() #%>% as.matrix()

#export clean as a dataset
write.csv(clean, file = "clock_data")


#code for age data
x <- trainclean
y <- trainage$clock_age 

#important fitting code
fit <- glmnet(x, y)
plot(fit) 
#plots the coefficients
print(fit)  #Df = number of coefficients, %Dev = goodness of fit, lambda = lambda
cvfit <- cv.glmnet( #alpha= 0 ridge
  x, y,
  alpha = 0.95,
  nfolds = 7)



plot(cvfit)
print(cvfit)
#makes a dataset with all the coeficients and their values 
lamb <- cvfit$lambda.min


coef <- coef(cvfit, s = lamb) #get model coefficient at a specified lambda 
coef_df <- as.data.frame(as.matrix(coef))
coef_df <- cbind(variable = rownames(coef_df), coef_df) #dataset with all coefficients and their value
colnames(coef_df)[2] <- "value"

# Filter and print coefficients greater than zero
nonzero_coefs <- coef_df[coef_df$value != 0, ]
print(nonzero_coefs)
count(nonzero_coefs)

lambda <- cvfit$lambda.min

final_model <- glmnet(
  x, y,
  alpha = 0.5,
  lambda = lamb
)
plot(final_model)
print(coef(final_model))
print(final_model %>% filter(coef > 0))

final <- coef(final_model)
print(final)


predictions <- predict(final_model, newx = testclean, s = lamb)
print(predictions)

predictions <- as.vector(predictions)
plot <- data.frame(real = testage$clock_age,
                   predicted = predictions)

ggplot(plot, aes( x = real, y = predictions)) +
  geom_point(alpha = 0.5)


#calculate mse

mse <- postResample(pred = predictions, obs = testage$clock_age)
  print(mse)
  
plot(model, xvar = "lamb", label = TRUE)


### split data into testing and training: put the worse indivs into the testing set
##1. TESTING data x values
testing <- borutadf2[(borutadf2$Individual == '23_RCD_SFO_012') | (borutadf2$Individual == '23_RCD_SFO_021')  
                 | (borutadf2$Individual == '23_RCD_SFO_007') | (borutadf2$Individual == '23_CV_SFO_017') | (borutadf2$Individual =='23_RCD_SFO_009')]
borutadf2

testing <- longdf[(longdf$Individual == '23_RCD_SFO_012') | (longdf$Individual == '23_RCD_SFO_021') 
                     | (longdf$Individual == '23_RCD_SFO_007') | (longdf$Individual == '23_CV_SFO_017') | (longdf$Individual =='23_RCD_SFO_009')]
#use this the other one broke for some reason
testing <- subset(df, (df$Individual %in% c('23_RCD_SFO_012','23_RCD_SFO_021', '23_RCD_SFO_007',
                                                         '23_CV_SFO_017', '23_RCD_SFO_009' )))
#spread wider
widetest <- testing %>% pivot_wider(names_from = chrBase, values_from = freqC)
#remove NA 
number <- c(26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1)
count_row_te <- sapply(widetest, function(x) sum(is.na(x))) #counts number of na values
cpg_te <- names(count_row_te[count_row_te %in% number]) #identifies rows with a specific number of NA values
widefilt_te <- widetest[, !(names(widetest) %in% cpg)] #removes cpg column in wide dataset with he requested number of NA values
#remove individual, final clean up, convert to matrix
testclean <- widefilt_te %>%
  dplyr::select(-Individual) %>% na.omit() %>% as.matrix()


##2. TRAINING data x values
training1 <- subset(borutadf, !(borutadf$Individual %in% c('23_RCD_SFO_012','23_RCD_SFO_021', '23_RCD_SFO_007',
                                               '23_CV_SFO_017', '23_RCD_SFO_009' )))

training <- subset(longdf, !(longdf$Individual %in% c( '23_RCD_SFO_012','23_RCD_SFO_021', '23_RCD_SFO_007',
                                                         '23_CV_SFO_017', '23_RCD_SFO_009' )))

widetrain <- training1 %>% pivot_wider(names_from = chrBase, values_from = freqC)
#remove NA 
count_row_tr <- sapply(widetrain, function(x) sum(is.na(x))) #counts number of na values
cpg_tr <- names(count_row_tr[count_row_tr %in% number]) #identifies rows with a specific number of NA values
widefilt_tr <- widetrain[, !(names(widetrain) %in% cpg)] #removes cpg column in wide dataset with he requested number of NA values
#remove individual, final clean up, convert to matrix
trainclean <- widefilt_tr %>%
  dplyr::select(-Individual) %>% as.matrix()

### Make train and test y values: 
#3. TESTING y
testage <- subset(brooktrout...Sheet1, (brooktrout...Sheet1$Individual %in% c( '23_RCD_SFO_012','23_RCD_SFO_021', '23_RCD_SFO_007', '23_CV_SFO_017', '23_RCD_SFO_009'
                                                ))) 
#4. TRAINING y
trainage <- subset(brooktrout...Sheet1, !(brooktrout...Sheet1$Individual %in% c( '23_RCD_SFO_012','23_RCD_SFO_021', '23_RCD_SFO_007',  '23_CV_SFO_017', '23_RCD_SFO_009')))                                                                                                                                                             

1
predictions <- predict(final_model, newx = trainclean, s = lamb)
print(predictions)

 fullage <- subset(brooktrout...Sheet1, !(brooktrout...Sheet1$Individual %in% c('23_CV_SFO_005',
                                                                                '23_WWD_SFO_016'))) 
 
 #stream import
 stream <- `stream(Sheet1)` %>% select(-X)
 
 #Random forest
x <- trainclean
 y <- trainage$clock_age# %>% as.factor()
 xtest <- testclean %>% as.factor()
 
 ystream <- stream$stream %>% as.factor()
 
 #filtering clean dataset with Kselect best features
 kselect <- selcted_features_chi2$Selected
 #filter clean dataset based on those selected CPG'S
 selectdf <- longdf %>% filter(chrBase %in% selcted_features_chi2$Selected)
 
 #Boruta feature selection for variable importance (vi)
 sites <- Boruta( x = x,
   y = y,
   maxRuns = 200)
 
 sites_stream <- Boruta(x = x, y = ystream )
 sites_stream
 
 print(sites)
 plot(sites)
 vi <- attStats(sites)
vi$cpg = row.names(vi)
 filtered_vi = vi %>% filter(decision == 'Confirmed') #filter datafrane boruta producves to grab on,y important cpg's
 
 filtered_vi %>%  ggplot(x = cpg, y= meanImp)
 
 #model building
 #select optimal number of tree with optrf
 #have to reformat data to match package requirements, merge age and cpgs together by individual?
 optdf <- merge(clean, agedf) %>% select(-Individual)
 
 data(optdf)
  result <- opt_prediction(y = optdf$clock_age, x = optdf[,-8])
  
    #alpha = 0.5,
    num.trees_values = c(250, 500, 700, 1000, 2000, 3000, 10000, 100000, 50000)
  )
 
modelstream <- randomForest(x, ystream, ntree = 1000)
modelstream

 model1 <- randomForest(x, y, ntree = 5000)
 model1
 
 modeltrain <- randomForest(x, y, ntree = 2000)
 modeltrain
 
 modelfinal <- randomForest(x, y, ntree = 2000)
 modelfinal
 
 which.min(modeltrain$mse)
 sqrt(modeltrain$mse[which.min(modeltrain$mse)])
 
 plot(modelfinal)
 varImpPlot(modelfinal)
 
 
 predict(model1, newdata = testclean)
 
 write.csv(filtered_vi)
 
 
 

