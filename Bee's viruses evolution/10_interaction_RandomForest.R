library(rsample)      # data splitting 
library(randomForest) # basic implementation
library(caret)        # an aggregator package for performing many machine learning mode
library(dplyr)
library(randomForestExplainer)

rm(list = ls())

# loading data
data <- read.csv("new_data.csv")

# data exploration: summary(), str()=glimpse(), execute in console
# data preprocessing: remove or change variables 
data$Geocluster <- as.factor(data$Geocluster)

data <- data[ ,-c(2,3,6,7,8,9,11,12,13,14,15,16,17)]
rownames(data) <- data$ID
data <- data[ ,-1]


# using Random Forest models to see interaction is my goal, 
# so I did not split the data into train & test and not tuning trees

# ABPV
# for reproduciblity
set.seed(123)
# default RF model
m1 <- randomForest(
  formula = ABPV ~ .,
  data    = data,
  na.action = na.omit,
  localImp=T
)

m1
plot(m1)
which.min(m1$mse)
sqrt(m1$mse[which.min(m1$mse)])

#-------------------------------------------------------

min_depth_frame.1 <- min_depth_distribution(m1)
plot_min_depth_distribution(min_depth_frame.1)
plot_min_depth_distribution(min_depth_frame.1, mean_sample = "relevant_trees", k = 20)


importance_frame.1 <- measure_importance(m1)
plot_multi_way_importance(importance_frame.1, size_measure = 'no_of_nodes')
plot_multi_way_importance(importance_frame.1, x_measure = "mse_increase", 
                          y_measure = "node_purity_increase", size_measure = "p_value", 
                          no_of_labels = 5)
plot_importance_ggpairs(importance_frame.1)
plot_importance_rankings(importance_frame.1)


vars.1 <- important_variables(importance_frame.1, k=13, 
                            measures=c("mean_min_depth","no_of_trees") )
interactions_frame.1 <- min_depth_interactions(m1, vars.1)
plot_min_depth_interactions(interactions_frame.1)

###########################################################################

# AM
# for reproduciblity
set.seed(123)
# default RF model
m2 <- randomForest(
  formula = AM ~ .,
  data    = data,
  na.action = na.omit,
  localImp=T
)

m2
plot(m2)
which.min(m2$mse)
sqrt(m2$mse[which.min(m2$mse)])

#-------------------------------------------------------

min_depth_frame.2 <- min_depth_distribution(m2)
plot_min_depth_distribution(min_depth_frame.2)
plot_min_depth_distribution(min_depth_frame.2, mean_sample = "relevant_trees", k = 20)


importance_frame.2 <- measure_importance(m2)
plot_multi_way_importance(importance_frame.2, size_measure = 'no_of_nodes')
plot_multi_way_importance(importance_frame.2, x_measure = "mse_increase", 
                          y_measure = "node_purity_increase", size_measure = "p_value", 
                          no_of_labels = 5)
plot_importance_ggpairs(importance_frame.2)
plot_importance_rankings(importance_frame.2)


vars.2 <- important_variables(importance_frame.2, k=11, 
                            measures=c("mean_min_depth","no_of_trees") )
interactions_frame.2 <- min_depth_interactions(m2, vars.2)
plot_min_depth_interactions(interactions_frame.2)

##################################################################################

# BQCV
# for reproduciblity
set.seed(123)
# default RF model
m3 <- randomForest(
  formula = BQCV ~ .,
  data    = data,
  na.action = na.omit,
  localImp=T
)

m3
plot(m3)
which.min(m3$mse)
sqrt(m3$mse[which.min(m3$mse)])

#-------------------------------------------------------

min_depth_frame.3 <- min_depth_distribution(m3)
plot_min_depth_distribution(min_depth_frame.3)
plot_min_depth_distribution(min_depth_frame.3, mean_sample = "relevant_trees", k = 20)


importance_frame.3 <- measure_importance(m3)
plot_multi_way_importance(importance_frame.3, size_measure = 'no_of_nodes')
plot_multi_way_importance(importance_frame.3, x_measure = "mse_increase", 
                          y_measure = "node_purity_increase", size_measure = "p_value", 
                          no_of_labels = 5)
plot_importance_ggpairs(importance_frame.3)
plot_importance_rankings(importance_frame.3)


vars.3 <- important_variables(importance_frame.3, k=11, 
                            measures=c("mean_min_depth","no_of_trees") )
interactions_frame.3 <- min_depth_interactions(m3, vars.3)
plot_min_depth_interactions(interactions_frame.3)









