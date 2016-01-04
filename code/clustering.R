library(downloader)

first_der <- function(X){
   X_temp <- cbind(NA, X)
   X_temp <- X_temp[, -dim(X_temp)[2]]
   X_der <- X-X_temp
   X_der <- X_der[,-1]
   X <- X_der
   X
}    


#+ Data download ----------------------------------------------------------
# Create a "Data" folder in your current working directory

dat_dir <- "data/"

download("https://www.dropbox.com/s/1kdoupic0vmdz38/clustingdata.zip?dl=0", "./data/clusteringdata.zip", mode="wb")
unzip("./data/clusteringdata.zip", exdir="./data", junkpaths=T, overwrite=T)


htsxt_spec <- read.csv(paste(dat_dir, "mir_withchem.csv", sep=""))
htsxt_spec_ethiosis <- read.csv(paste(dat_dir, "mir_ethiosis_withchem.csv", sep=""))

htsxt_afsis <- htsxt_spec[, 29:dim(htsxt_spec)[2]]
htsxt_ethiosis <- htsxt_spec_ethiosis[,26:dim(htsxt_spec_ethiosis)[2]]

# only one row of ethiosis data (row 1905)
htsxt_all <- rbind(htsxt_afsis, htsxt_ethiosis[1,])
htsxt_all_der <- first_der(htsxt_all)

cor_htsxt <- cor(t(htsxt_all))
cor_htsxt_dist <- as.dist(1-cor_htsxt)


dist_htsxt <- dist(htsxt_all)
hclust_clusters <- hclust(cor_htsxt_dist, method="average")

plot(hclust_clusters)


# kmeans_clusters <- kmeans(cor_htsxt_dist, centers=2)

# # train, test sets 
#
library(randomForest)
train_index <- sample(1:dim(htsxt_all)[1], size=floor(0.8*dim(htsxt_all)[1]))

train_set <- htsxt_all[train_index,]

test_set <- htsxt_all[-train_index,]


# first derivative
train_set_der <- first_der(train_set)
test_set_der <- first_der(test_set)   

carbon <- c(htsxt_spec$Acidified.Carbon, htsxt_spec_ethiosis$SOC)

carbon_train <- carbon[train_index]

train_data <- cbind(carbon_train, train_set_der)
train_data <- na.omit(train_data)
rf_model <- randomForest(carbon_train~., data= train_data, xtest=test_set_der, ntree=50)
mse_vec <- sqrt(mean((carbon[-train_index]- rf_model$test$predicted)^2, na.rm=TRUE))

# using cross-validation to figure which level of clusters 

heights_all <- sort(hclust_clusters$heights)



for(h in heights_all[1:5]){
    mem_tree <- cuttree(hclust_clusters, h = h)
    resid <- NA
    for(k in 1:length(unique(mem_tree))){
        X_data <- htsxt_all[mem_tree==k,]
        carbon_data <- carbon[mem_tree==k]

        train_index <- sample(1:dim(X_data)[1], size=floor(0.8*dim(X_data)[1]))
        train_set <- X_data[train_index,]
        test_set <- X_data[-train_index,]

        train_set_der <- first_der(train_set)
        test_set_der <- first_der(test_set)   
        carbon_train <- carbon_data[train_index]

        train_data <- cbind(carbon_train, train_set_der)
        train_data <- na.omit(train_data)
        rf_model <- randomForest(carbon_train~., data= train_data, xtest=test_set_der)

        resid <- c(resid, (carbon_data[-train_index]- rf_model$test$predicted)) 
    }
    mse_vec <- c(mse_vec, mean(resid^2, na.rm=TRUE))
}


k_min <- heights_all[which.min(mse_vec)]

