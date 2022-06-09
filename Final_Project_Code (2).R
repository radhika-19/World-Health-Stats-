ibrary(dplyr)
library(tidyverse)
#install.packages("mice")
library(mice)
#install.packages("VIM")
library(VIM)
library(cluster)
#install.packages("psych")
library(psych)
#install.packages("corrplot")
library(corrplot)
#install.packages("factoextra")
library(factoextra)
#install.packages("sjmisc")
library(sjmisc)
library(readr)
library(rgl)
library(cluster)
library(car) #- lavene test


data_org <- read_csv("D:/Semester 4/DANA/Project/Final Project/data.csv")
data <- data_org

#Removing the columns Year and YearCode since they are constant.
#Removin the column Country Code, since its like primary key and unique for each row.
data <- data[,-c(1,2,3)]

#Dealing with Missing values
Null_Cnt <- sapply(data, function(x){ sum(is.na(x))})
Null_percnt <- sapply(data, function(x){ round((sum(is.na(x))/length(x))*100,2) })
Null_Smry <- cbind(Null_Cnt,Null_percnt)
Null_Smry <- as.data.frame(Null_Smry)

col_list <- row.names(Null_Smry)[Null_percnt > 25] #"ari" "gen_eq_rate"  "gov_debt" "gini" "ishare_low20" "hiv_fe15up

data <- data[,!(names(data) %in% col_list)]

# removing rows with more than 25% null value
data<- data[which(round(rowSums(is.na(data))/dim(data)[2]*100)<22),] #34 observations removed - only 1 NO-OBS in data
data


summary(data)

#######################################################
#Data Imputation using MICE
########################################################
md.pattern(data, color = c("orange","dark green"),rotate.names = TRUE) #161 observations without NA's

mice_plot <- aggr(data, col=c('brown','orange'),
                  numbers=TRUE, sortVars=TRUE,
                  labels=names(data), cex.axis=.7,
                  gap=3, ylab=c("Missing data","Pattern"))

imputed_Data <- mice(data[,-c(1,2)] , m=5, maxit = 50, method = 'cart', seed = 500)
summary(imputed_Data)

#complete_data <- merge_imputations(data,imputed_Data,summary="hist")

complete_data <- complete(imputed_Data,5)
nrow(complete_data[complete.cases(complete_data),]) #180

#corPlot(complete_data, numbers=FALSE, zlim = NULL, n.legend=5, scale=TRUE,stars=TRUE,  MAR=TRUE, cex.axis=0.6)
corrplot(cor(complete_data), type = 'upper', method = 'shade',col=colorRampPalette(c("dark green","lightblue","brown"))(100), tl.cex = 0.7)

#From the correlation matrix we see that, the Fertility rate and Mortality rate is highly correlated.

#Finding outliers using mahalanobis distances
# Finding the center point 
complete_data.center  = colMeans(complete_data)

# Finding the covariance matrix
complete_data.cov = cov(complete_data)

# Finding distances
distances <- mahalanobis(x = complete_data , center = complete_data.center , cov = complete_data.cov, tol=1e-40)


# Cutoff value for ditances from Chi-Sqaure Dist. 
# with p = 0.95 df = 6 which in ncol(df_num)
cutoff <- qchisq(p = 0.95 , df = ncol(complete_data))
complete_data$distances <- as.factor(ifelse(distances > cutoff, 0, 1))

complete_data_num <- complete_data[,-c(24)]

#complete_data_num = cbind(data[,c(1,2)],complete_data_num)
#complete_data_num<- complete_data_num[,-c(1,2)]
###################################################################################
#We will confine the number of variables used for the analysis using PCA.
###################################################################################
data_PCA <- prcomp(complete_data_num, scale = TRUE)
summary(data_PCA)
print(data_PCA)
#data_PCA$x

#dev.off()

# Selecting top 6 Principal components
comp <- data.frame(data_PCA$x[,1:6])


#Scree Plot
fviz_eig(data_PCA, barfill = "dark green", barcolor ="brown", linecolor  = "black" )

#Graph of individuals
fviz_pca_ind(data_PCA,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
###########################################################################################
#Variable-PCA
fviz_pca_var(data_PCA,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

#Biplot
fviz_pca_biplot(data_PCA, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)

# Eigenvalues
eig.val <- get_eigenvalue(data_PCA)
eig.val
###########################################################################################
# Determine number of clusters using WSS and Gap Stat
wss <- (nrow(complete_data_num)-1)*sum(apply(complete_data_num,2,var))
for (i in 2:8) wss[i] <- sum(kmeans(complete_data_num,
                                    centers=i)$withinss)
plot(1:8, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

fviz_nbclust(x = scale(complete_data_num),FUNcluster = kmeans, method = 'gap_stat' )

#Euclidean distance
ed = dist(complete_data_num)


####### k=3
km2 <- kmeans(comp, 3, nstart=25, iter.max=1000)
km2$size # 50 43 87

# Adding the region and country to the data
complete_data = cbind(data[,c(1,2)],complete_data)

complete_data$cluster <- km2$cluster
table(complete_data$WHO_Region, complete_data$cluster)
##          1  2  3
##AFR     41  0  5
##AMR      1  4 26
##EMR      5  3 17
##EUR      0 30 16
##NA- AMR  0  0  1
##SEAR     1  1  9
##WPR      2  5 13



fviz_cluster(km2, geom = "point", data =complete_data_num) + ggtitle(" K = 3")



names(complete_data)[1] <- "Country"

complete_data[,c(1,2,27)]<-lapply(complete_data[,c(1,2,27)], as.factor)

par(mfrow=c(3,4)) # define 2x5 multiframe graphic
for (i in 3:(ncol(complete_data)-2)) # make box plots for all columns except the cluster label
{
  if (is.numeric(complete_data[,i])) # if numeric -> boxplot
  {
    plot(complete_data$cluster, complete_data[,i], main= colnames(complete_data)[i], col= "blue")
  }
  else # if factor -> barplot
  {
    count<-table(complete_data[,i],complete_data$cluster)
    barplot(count, legend = rownames(count), main= colnames(complete_data)[i], col=rainbow(6))
  }
}

par(mfrow=c(1,1)) # define 2x5 multiframe graphic

count<-table(complete_data[,"cluster"],complete_data$cluster)
barplot(count, legend = rownames(count), main= colnames(complete_data)[27], col=rainbow(6))

barplot(count,
        main="Cluster Distribution",
        xlab="Cluster",
        ylab="Count",
        border="red",
        col="blue"
)

write.csv(complete_data, "D:/Semester 4/DANA/Project/Final Project/Complete_data_with_cluster.csv",row.names = FALSE)


##############################################
###Cluster Validation
##############################################
# Bartlett's test when data is normally distributed
# H(0) = There is no difference between the variances of 3 cluster groups
# H(a) = The 3 cluster groups have variance

bartlett.test(c_gdp ~ cluster, data = complete_data)

# p-value =  2.2e-16, means variance in c_gdp is significantlly different for the 3 cluster groups

# Lavene's test when data is not normally distributed

leveneTest(c_gdp ~ cluster, data = complete_data)

# p-value =  0.001226 , means variance in c_gdp is significantlly different for the 3 cluster groups

#########################Anova Test ###########################

# H(0) = There is no difference in the gdp means of 3 cluster groups
# H(a) = There is difference between gdp means of 3 cluster groups

oneway.test(c_gdp ~ cluster, data = complete_data, var.equal = TRUE)
# p-value =  0.000645 , means variance in c_gdp is significantlly different for the 3 cluster groups




#2.--------------------------------------------
# Health Spend

bartlett.test(healthexp_gdp ~ cluster, data = complete_data)   # pvalue - 0.04 , reject null

oneway.test(healthexp_gdp ~ cluster, data = complete_data, var.equal = TRUE)
# p-value =  2.2e-16 , means variance in c_gdp is significantlly different for the 3 cluster groups

#3.--------------------------------------------
# Mortality rate

bartlett.test(mortrate_inf ~ cluster, data = complete_data)   # pvalue - 1.07e-14 , reject null

oneway.test(mortrate_inf ~ cluster, data = complete_data,  var.equal = TRUE)
# p-value =  2.2e-16 , means variance in c_gdp is significantlly different for the 3 cluster groups

#4.--------------------------------------------
# Life Expectancy

bartlett.test(lifeexp ~ cluster, data = complete_data)   # pvalue - 3.927e-05 , reject null

oneway.test(lifeexp ~ cluster, data = complete_data,  var.equal = TRUE)
# p-value =  2.2e-16 , means variance in c_gdp is significantlly different for the 3 cluster groups

#5.--------------------------------------------
# Birth Rate


bartlett.test(brate ~ cluster, data = complete_data)   # pvalue - 4.006e-05 , reject null

oneway.test(brate ~ cluster, data = complete_data,  var.equal = TRUE)
# p-value =  2.2e-16 , means variance in c_gdp is significantlly different for the 3 cluster groups

