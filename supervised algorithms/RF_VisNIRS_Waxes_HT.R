# Load Required Libraries
library(readxl)
library(prospectr)
library(caret)
library(stringr)
library(dplyr)
library(ggplot2)
library(ranger)

figures_dir <- file.path(getwd(), "Figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
app_dir_candidates <- c(file.path(getwd(), "App"), file.path(getwd(), "..", "App"))
app_dir <- app_dir_candidates[file.exists(app_dir_candidates)][1]
if (is.na(app_dir)) {
  app_dir <- app_dir_candidates[1]
  dir.create(app_dir, recursive = TRUE, showWarnings = FALSE)
}

# Load Parallelization when available
cl <- NULL
if (requireNamespace("doParallel", quietly = TRUE)) {
  cl <- parallel::makePSOCKcluster(8)
  doParallel::registerDoParallel(cl)
}

# Load Data
pw_data <- read_excel("~/Documents/Doctorado/Tesis Doctoral/Investigación Cepsa/Vis-NIR/XDS-NIR_FOSS/Estudio según Tipo de Parafina e Hidrotratamiento/NIRS_HT_PW.xlsx",
                      sheet = "HT_Class")

# First derivative and Savitzky Golay Smoothing
pw_data$Sample <- as.factor(pw_data$Sample)

sgvec <- savitzkyGolay(X = pw_data[,-1], p = 3, w = 11, m = 1)
pw_sg <- cbind.data.frame(Sample = pw_data$Sample, sgvec)

# Data partition
set.seed(537)

intrain <- createDataPartition(y = pw_sg$Sample, 
                               p = 0.7, 
                               list = FALSE)
pw_train <- pw_sg[intrain,]
pw_test <- pw_sg[-intrain,]

# Ensure valid column names
colnames(pw_train) <- make.names(colnames(pw_train))
colnames(pw_test) <- make.names(colnames(pw_test))

# Define parameter grid for ranger
rf_grid <- expand.grid(
  mtry = c(sqrt(ncol(pw_train[,-1]))),             # Example values for mtry
  splitrule = c("gini"),            # Default splitrule for classification
  min.node.size = c(1, 5, 10)       # Example values for min.node.size
)

ntrees <- c(seq(2, 100, 2))         # Define a sequence for the number of trees
params <- expand.grid(ntrees = ntrees)

store_maxnode <- vector("list", nrow(params))

# Hyperparameter tuning
set.seed(537)
trctrl <- trainControl(method = "cv", number = 10)

start_time_1 <- Sys.time()

for(i in 1:nrow(params)) {
  ntree <- params[i, 1]
  set.seed(537)
  rf_model <- train(
    Sample ~ ., 
    data = pw_train,
    method = "ranger",
    metric = "Accuracy",
    tuneGrid = rf_grid,
    trControl = trctrl,
    num.trees = ntree,
    importance = "permutation"
  )
  store_maxnode[[i]] <- rf_model
}

names(store_maxnode) <- paste("ntrees:", params$ntrees)

rf_results <- resamples(store_maxnode)
rf_results

lapply(store_maxnode, 
       function(x) x$results[x$results$Accuracy == max(x$results$Accuracy),])

total_time_1 <- Sys.time() - start_time_1
print(total_time_1)

# Accuracy vs ntrees
rf_plot <- cbind.data.frame(Accuracy = colMeans(rf_results$values[, c(seq(2, 100, 2))]),
                            Ntrees = c(seq(1, 100, 2)))
rf_plot <- as.vector(rf_plot)
accuracy <- (rf_plot$Accuracy) * 100
ntrees <- rf_plot$Ntrees

plot(x = ntrees,
     y = accuracy,
     type = "b", 
     pch = 19, 
     lty = 2,
     col = "#287D8EFF",
     xlab = "Number of decision trees",
  ylab = "Accuracy (%) (10-Fold CV)")

# Final RF model with ranger
set.seed(537)

start_time_2 <- Sys.time()

best_rf <- ranger(
  formula = Sample ~ ., 
  data = pw_train, 
  num.trees = 100, 
  mtry = c(sqrt(ncol(pw_train[,-1]))),
  splitrule = "gini",
  min.node.size = 5,
  importance = 'permutation',
  classification = TRUE
)

total_time_2 <- Sys.time() - start_time_2
print(total_time_2)

saveRDS(best_rf, file.path(app_dir, "rf.rds"))

# Train set performance for RF
pred_train <- predict(best_rf, pw_train)$predictions
cmatrix_train <- confusionMatrix(as.factor(pred_train), pw_train$Sample)
print(cmatrix_train)

# Test set performance for RF
pred_test <- predict(best_rf, pw_test)$predictions
cmatrix_test <- confusionMatrix(as.factor(pred_test), pw_test$Sample)
print(cmatrix_test)

# Variable Importance
var_imp <- importance(best_rf)

# Convert to DataFrame for visualization
var_imp_df <- data.frame(
  Variable = names(var_imp),
  Importance = var_imp
)

# Sort by descending importance
var_imp_df <- var_imp_df[order(var_imp_df$Importance, decreasing = TRUE), ]

print(var_imp_df)

# Remove "X" prefix from labels on the vertical axis
var_imp_df$Variable <- gsub("^X", "", var_imp_df$Variable)

# Create the plot with a gradient of colors
importance_plot <- ggplot(var_imp_df[1:10, ], aes(x = reorder(Variable, Importance), y = Importance, fill = Importance)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradientn(colors = c("#FDE725FF", "#287D8EFF", "#440154FF", "#73D055FF", "#404788FF")) +
  labs(title = "Variable Importance (Top 10)",
       x = "Variables",
       y = "Importance") +
  theme_light() +
  theme(legend.position = "right")

importance_plot

ggsave(
  filename = file.path(figures_dir, "Fig.5.png"),
  plot = importance_plot,
  width = 10,
  height = 10,
  units = "in",
  dpi = 600,
  bg = "white"
)

# Stop Parallelization
if (!is.null(cl)) {
  parallel::stopCluster(cl)
}