# 03_fit_ML

# this script takes the processed recruitment and covariate data and fits ML models.

library(here)
library(tidyverse)
library(ncdf4)
library(sf)
library(viridis)
library(ggpubr)
library(paran)
library(corrplot)
library(tidymodels)
library(glmnet)
library(nnet)
library(mgcv)
library(ranger)
library(vip)
library(factoextra)

#### Load data ####

# Load covariates
feature_df <- read.csv(here::here("ML_inputs", "features_2025-01-14.csv"))
feature_df %>% dplyr::rename(Yr = X) -> feature_df

# correlation structure in our features
feature_cor_mat <- cor(dplyr::select(feature_df, -Yr))

as.data.frame(feature_cor_mat) %>% 
  rownames_to_column("Var1") %>% 
  pivot_longer(cols = -c("Var1"), names_to = "Var2", values_to = "Value") -> feature_cor_values

feature_cor_values %>%
  filter(!(Var1 == Var2)) %>% 
  arrange(desc(abs(Value)))
corrplot(feature_cor_mat)
# well, some are quite correlated but none are >= 0.8

# Load each of the recruitment dfs
sablefish_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "sablefish_recdevs.csv"))
black_rf_nca_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "black_rf_nca_recdevs.csv"))
black_rf_or_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "black_rf_or_recdevs.csv"))
black_rf_wa_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "black_rf_wa_recdevs.csv"))
copper_rf_north_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "copper_rf_north_recdevs.csv"))
copper_rf_south_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "copper_rf_south_recdevs.csv"))
canary_rf_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "canary_rf_recdevs.csv"))
dover_sole_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "dover_sole_recdevs.csv"))
lingcod_north_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "lingcod_north_recdevs.csv"))
lingcod_south_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "lingcod_south_recdevs.csv"))
petrale_sole_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "petrale_sole_recdevs.csv"))
rex_sole_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "rex_sole_recdevs.csv"))
shortspine_thornyhead_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "shortspine_thornyhead_recdevs.csv"))
vermilion_rf_sca_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "vermilion_rf_sca_recdevs.csv"))
vermilion_rf_nca_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "vermilion_rf_nca_recdevs.csv"))
vermilion_rf_or_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "vermilion_rf_or_recdevs.csv"))
vermilion_rf_wa_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "vermilion_rf_wa_recdevs.csv"))
quillback_rf_ca_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "quillback_rf_ca_recdevs.csv"))
rougheye_rf_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "rougheye_rf_recdevs.csv"))
yellowtail_rf_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "yellowtail_rf_recdevs.csv"))

# join features with recruitment data
join_features_recdevs <- function(recdevs){
  full_join(recdevs, feature_df, by = "Yr") %>% 
    dplyr::filter(!(is.na(bbv_PC1)) & !(is.na(Value))) -> output
  
  return(output)
}

join_features_recdevs(sablefish_cat_rec) -> sablefish
join_features_recdevs(black_rf_nca_cat_rec) -> black_rf_nca
join_features_recdevs(black_rf_or_cat_rec) -> black_rf_or
join_features_recdevs(black_rf_wa_cat_rec) -> black_rf_wa
join_features_recdevs(copper_rf_north_cat_rec) -> copper_rf_north
join_features_recdevs(copper_rf_south_cat_rec) -> copper_rf_south
join_features_recdevs(canary_rf_cat_rec) -> canary_rf
join_features_recdevs(dover_sole_cat_rec) -> dover_sole
join_features_recdevs(lingcod_north_cat_rec) -> lingcod_north
join_features_recdevs(lingcod_south_cat_rec) -> lingcod_south
join_features_recdevs(petrale_sole_cat_rec) -> petrale_sole
join_features_recdevs(rex_sole_cat_rec) -> rex_sole
join_features_recdevs(shortspine_thornyhead_cat_rec) -> shortspine_thornyhead
join_features_recdevs(vermilion_rf_sca_cat_rec) -> vermilion_rf_sca
join_features_recdevs(vermilion_rf_nca_cat_rec) -> vermilion_rf_nca
join_features_recdevs(vermilion_rf_or_cat_rec) -> vermilion_rf_or
join_features_recdevs(vermilion_rf_wa_cat_rec) -> vermilion_rf_wa
join_features_recdevs(quillback_rf_ca_cat_rec) -> quillback_rf_ca
join_features_recdevs(rougheye_rf_cat_rec) -> rougheye_rf
join_features_recdevs(yellowtail_rf_cat_rec) -> yellowtail_rf

#### Explore autocorrelation structure in the data ####

# plot all of them
acf(sablefish$Value)
acf(black_rf_nca$Value)
acf(black_rf_or$Value)
acf(black_rf_wa$Value) # autocorrelated at lag 1
acf(copper_rf_north$Value) # autocorrelated at lag 1
acf(copper_rf_south$Value) # autocorrelated at lag 1 and 2
acf(canary_rf$Value) # autocorrelated at lag 1, 2, 3
acf(dover_sole$Value) # autocorrelated at lag 1; negatively autocorrelated at lag 4, 5, 6
acf(lingcod_north$Value)
acf(lingcod_south$Value) # autocorrelated at lag 1 and 2
acf(petrale_sole$Value) # autocorrelated at lag 1
acf(rex_sole$Value) # autocorrelated at lag 1; negatively autocorrelated at lag 8, 9, 10, 11
acf(shortspine_thornyhead$Value) # completed autocorrelated - should be dropped. shortspine has known issues
# with aging that lead to recruitment not being estimated in the assessment model well at all.
acf(vermilion_rf_sca$Value) # autocorrelated at lag 1
acf(vermilion_rf_nca$Value)
acf(vermilion_rf_or$Value) # autocorrelated at lag 1
acf(vermilion_rf_wa$Value) # autocorrelated at lag 1
acf(quillback_rf_ca$Value)
acf(rougheye_rf$Value) # autocorrelated at lag 1
acf(yellowtail_rf$Value)





#### Explore clustering of years by environmental covariates ####

# Extract only environmental covariates
feature_df %>%
  column_to_rownames("Yr") -> features_only

# run a PCA
features_pca <- prcomp(features_only)
summary(features_pca)

# extra pca scores
as.data.frame(features_pca$x) %>%
  rownames_to_column("Yr") %>%
  mutate(Yr = as.integer(Yr)) %>%
  left_join(sablefish %>% dplyr::select(Yr, one_sd_boom_bust), by = "Yr") -> features_pca_scores_sablefish

# extract loadings
features_pca_loadings <- as.data.frame(features_pca$rotation) %>%
  rownames_to_column("variable") %>%
  dplyr::select(variable, PC1, PC2)

# Scale loadings for plotting (so arrows fit nicely in the score plot)
loading_scale <- 3
features_pca_loadings <- features_pca_loadings %>%
  mutate(PC1 = PC1 * loading_scale,
         PC2 = PC2 * loading_scale)

# Calculate variance explained
var_explained <- round(100 * features_pca$sdev^2 / sum(features_pca$sdev^2), 1)

# Plot
ggplot() +
  # Year points colored by recruitment category
  geom_point(data = features_pca_scores_sablefish,
             aes(x = PC1, y = PC2, color = one_sd_boom_bust),
             size = 3) +
  # Year labels
  geom_text(data = features_pca_scores_sablefish,
            aes(x = PC1, y = PC2, label = Yr, color = one_sd_boom_bust),
            nudge_y = 0.2, size = 3) +
  # Loading arrows
  geom_segment(data = features_pca_loadings,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "gray40") +
  # Loading labels
  geom_text(data = features_pca_loadings,
            aes(x = PC1, y = PC2, label = variable),
            color = "gray40", size = 3, nudge_y = 0.15) +
  # Axis labels with variance explained
  labs(x = paste0("PC1 (", var_explained[1], "%)"),
       y = paste0("PC2 (", var_explained[2], "%)"),
       color = "Recruitment") +
  scale_color_manual(values = c("Average" = "gray60",
                                "Boom" = "steelblue",
                                "Bust" = "tomato")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray80") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray80") +
  theme_minimal()


# what if you just plot two of the covariates at a time?
# as an example:
ggplot(sablefish, aes(x = BEUTI_north, y = sst_PC1, color = one_sd_boom_bust)) +
  geom_point() +
  scale_color_manual(values = c("Average" = "gray60",
                                "Boom" = "steelblue",
                                "Bust" = "tomato")) +
  theme_minimal()




#### Run models - use leave-future-out cross-validation, for last ten years, to compare performance ####

# drop all columns that aren't ID, features, or outcome
sablefish %>% 
  dplyr::select(-c("Value", "Parm_StDev", "recr_scaled")) -> sablefish

# code one_sd_boom_bust as a factor
sablefish$one_sd_boom_bust <- as.factor(sablefish$one_sd_boom_bust)

# get first ten years
first_10_yrs <- (min(sablefish$Yr):(min(sablefish$Yr)+9))

# get last ten years
last_10_yrs <- (max(sablefish$Yr)-9):max(sablefish$Yr)

# First, we will split our data into training and testing data. The testing
# data is the last ten years of our dataset; the training data is the rest

# Second, we will take our training data and use rolling_origin() to create folds
# within it; following the workflow of Allen Akselrud 2024, we will start with
# the first ten years as our first set of analysis data, and continue to predict
# one year ahead

# 

# Here, we are going to use essentially rolling LFO-CV twice: once for 

# create the resamples for ASSESSMENT from the TRAINING data
assessment_resamples <- rolling_origin(
  subset(sablefish, !(Yr %in% last_10_yrs)), # don't include the testing data!
  initial = 20,      # number of initial years
  assess  = 1,       # predict 1 future year
  cumulative = TRUE,
  skip = 0
)

# create the resamples for TESTING
lfo_resamples <- rolling_origin(
  sablefish,
  initial = nrow(subset(sablefish, !(Yr %in% last_10_yrs))),      # number of initial years
  assess  = 1,       # predict 1 future year
  cumulative = TRUE,
  skip = 0
)  

# create a table that has the truth (one_sd_boom_bust) for the left out years
lfo_df <- lfo_resamples %>%
  mutate(
    assessment_data = map(splits, assessment),
    .row = map_int(splits, ~ as.integer(assessment(.x)$Yr - min(sablefish$Yr) + 1))
  ) %>%
  dplyr::select(id, assessment_data, .row) %>%
  unnest(assessment_data) %>%
  dplyr::select(id, one_sd_boom_bust, .row)


#### Fit a null model as a baseline ####

# This model just uses the proportions from the training data

# Calculate class proportions from training data
null_probs <- subset(sablefish, !(Yr %in% last_10_yrs)) %>%
  count(one_sd_boom_bust) %>%
  mutate(prop = n / sum(n))

# Construct a null prediction tibble using training proportions
null_predictions <- lfo_df %>%
  mutate(
    .pred_Average = null_probs$prop[null_probs$one_sd_boom_bust == "Average"],
    .pred_Boom    = null_probs$prop[null_probs$one_sd_boom_bust == "Boom"],
    .pred_Bust    = null_probs$prop[null_probs$one_sd_boom_bust == "Bust"]
  )

# Calculate null log-loss
null_predictions %>%
  mn_log_loss(truth = one_sd_boom_bust,
              .pred_Average, .pred_Boom, .pred_Bust) -> null_log_loss

# Calculate null ROC-AUC
# calculate metrics
null_predictions %>% 
  roc_auc(truth = one_sd_boom_bust,
          .pred_Average, .pred_Boom, .pred_Bust) -> roc_auc_output

# overall accuracy
null_predictions %>% 
  mutate(.pred_class = factor("Average", levels = c("Average", "Boom", "Bust"))) %>% 
  accuracy(truth = one_sd_boom_bust, estimate = .pred_class) -> accuracy_output

null_log_loss %>% 
  bind_rows(roc_auc_output, accuracy_output) %>% 
  mutate(model = "null") -> null_metrics

# start a table to compare models with the null metrics table
model_comparison_table <- null_metrics


# fit each of our candidate models
# for (i in last_10_yrs){
  
  # get data for all years up to the one to be predicted
  dat <- subset(sablefish, Yr < i)

  #### Fit multinomial logistic regression ####
  
  
  #### Fit random forest ####
  # parallelize
  cores <- parallel::detectCores()
  
  # set up the model
  rf_mod <- 
    rand_forest(mtry = 2, min_n = 4, trees = 1000) %>% 
    set_engine("ranger", num.threads = cores) %>% 
    set_mode("classification")
  
  ## testing zone
  # set up the model
  rf_mod <- 
    rand_forest(mtry = tune(), min_n = tune(), trees = tune()) %>% 
    set_engine("ranger", num.threads = cores) %>% 
    set_mode("classification")
  
  # set up the recipe - formula + code columns
  rf_recipe <- recipe(one_sd_boom_bust ~., data = sablefish, !(Yr %in% last_10_yrs)) %>% 
    update_role(Yr, new_role = "ID")
  
  # package these together into a workflow for the model,
  workflow() %>%
    add_model(rf_mod) %>% 
    add_recipe(rf_recipe)  -> rf_workflow
  
  # Tune the model
  # create at uning grid
  rf_grid <- grid_regular(
    mtry(range = c(2, 10)),        # can't exceed number of predictors
    min_n(range = c(2, 10)), # keep well below number of training rows
    trees(range = c(50, 2500)),
    # levels = c(5, 5, 25) # try more trees than mtry or min_n - takes forever
    levels = c(5, 5, 5)
  )
  
  # tune the hyperparameters on the assessment resamples
  rf_res <- 
    rf_workflow %>% 
    tune_grid(assessment_resamples,
              grid = rf_grid,
              control = control_grid(save_pred = TRUE),
              metrics = metric_set(roc_auc))
  
  # pull out the top models
  rf_res %>% 
    show_best(metric = "roc_auc", n = 15)
  # terrible!!
  
  autoplot(rf_res)
  
  # select the best model according to the ROC AUC metric
  rf_best <- rf_res %>% 
    select_best(metric = "roc_auc")
  
  # run the model 100 times to see how random seed affects final predictions
  
  
  ## end testing zone
  
  # fit the model
  rf_mod %>% 
    fit(one_sd_boom_bust ~. -Yr, data = dat) -> rf_fit
  
  one_yr_pred <- predict(rf_fit, 
                         new_data = subset(sablefish, Yr == i),
                         type = "prob")
  
  rf_res <- rf_mod %>%
    fit_resamples(
      one_sd_boom_bust ~ . - Yr,
      resamples = lfo_resamples,
      control = control_resamples(save_pred = TRUE)
    )
  
  rf_one_yr_ahead_forecasts <- collect_predictions(rf_res)
  
  # calculate metrics
  rf_one_yr_ahead_forecasts %>% 
    roc_auc(truth = one_sd_boom_bust,
            .pred_Average, .pred_Boom, .pred_Bust) -> rf_roc_auc_output
  
  # overall accuracy
  rf_one_yr_ahead_forecasts %>% 
    accuracy(truth = one_sd_boom_bust, estimate = .pred_class) -> rf_accuracy_output
  
  # mean log-loss
  rf_one_yr_ahead_forecasts %>% 
    mn_log_loss(truth = one_sd_boom_bust,
              .pred_Average, .pred_Boom, .pred_Bust) -> rf_log_loss_output
  
  # combine metrics for rf into one table
  rf_log_loss_output %>% 
    bind_rows(rf_roc_auc_output, rf_accuracy_output) %>% 
    mutate(model = "rf") -> rf_metrics
  
  # add rf output to model comparison table
  model_comparison_table %>% 
    bind_rows(rf_metrics) -> model_comparison_table
  

# }




#### Recruitment as categorical ####


#### Split data into training (first 30 years) and testing (last 10 years) ####
# drop all columns that aren't ID, features, or outcome
sablefish %>% 
  dplyr::select(-c("Value", "Parm_StDev", "recr_scaled")) -> sablefish

# code one_sd_boom_bust as a factor
sablefish$one_sd_boom_bust <- as.factor(sablefish$one_sd_boom_bust)

sablefish_split <- initial_validation_time_split(sablefish, prop = c(0.6, 0.2))
sablefish_training <- training(sablefish_split)
sablefish_validation <- validation(sablefish_split)
sablefish_testing <- testing(sablefish_split)

# these train vs. test datasets are quite different, much higher proportion of boom (and no bust) in testing dataset
table(sablefish_training$one_sd_boom_bust)
table(sablefish_validation$one_sd_boom_bust)
table(sablefish_testing$one_sd_boom_bust)

#### Fit model on training set (random forest, but can change this) ####

# Split non-testing data into training and validation
# sablefish_val_set <- validation_split(sablefish_other, strata = one_sd_boom_bust, prop = 0.75)

#### Build the multinomial logistic regression model - version 2 - hyperparameter tuning ####

# set up the model
mlr <- multinom_reg(penalty = tune())

# set up the recipe - formula + code columns
mlr_recipe <- recipe(one_sd_boom_bust ~., data = sablefish_training) %>% 
  update_role(Yr, new_role = "ID")

# package these together into a workflow for the model,
workflow() %>%
  add_model(mlr) %>% 
  add_recipe(mlr_recipe)  -> mlr_workflow

# create the grid for tuning
mlr_reg_grid <- tibble(penalty = seq(0, 3, length.out = 30))

# Tune the model
mlr_res <- 
  mlr_workflow %>% 
  tune_grid(validation_set(sablefish_split),
            grid = mlr_reg_grid,
            control = control_grid(save_pred = TRUE),
            metrics = metric_set(roc_auc))

mlr_plot <- 
  mlr_res %>% 
  collect_metrics() %>% 
  ggplot(aes(x = penalty, y = mean)) + 
  geom_point() + 
  geom_line() + 
  ylab("Area under the ROC Curve")

mlr_plot 

# pull out the top models
top_models <-
  mlr_res %>% 
  show_best(metric = "roc_auc", n = 15) %>%
  arrange(penalty)

# pull out the model with equivalent ROC but highest penalty
top_models %>% 
  filter(mean == max(top_models$mean)) %>% 
  filter(penalty == max(penalty)) -> mlr_best_model

# visualize validation set AUC curve
mlr_auc <- 
  mlr_res %>% 
  collect_predictions(parameters = mlr_best_model) %>% 
  roc_curve(  truth = one_sd_boom_bust,
              c(".pred_Average", ".pred_Boom", ".pred_Bust")) %>% 
  mutate(model = "Multinomial Logistic Regression")

autoplot(mlr_auc)

### The last fit - fit the best model to combined training + validation data, evaluate on tesitng

# the last model
last_mlr_mod <- multinom_reg(penalty = mlr_best_model$penalty)

# the last workflow
last_mlr_workflow <- 
  mlr_workflow %>% 
  update_model(last_mlr_mod)

# the last fit
last_mlr_fit <- last_mlr_workflow %>% 
  last_fit(sablefish_split)
# here's a fun thing - our last time period has no "bust" recruitment.

# evaluate performance
last_mlr_fit %>% 
  collect_metrics()
# it's terrible!

# look at your parameters
last_mlr_fit %>% 
  extract_fit_parsnip()



#### Build the random forest model - version 2 - hyperparameter tuning ####

# parallelize
cores <- parallel::detectCores()

# set up the model
rf_mod <- 
  rand_forest(mtry = tune(), min_n = tune(), trees = 1000) %>% 
  set_engine("ranger", num.threads = cores) %>% 
  set_mode("classification")

# set up the recipe - formula + code columns
rf_recipe <- recipe(one_sd_boom_bust ~., data = sablefish_training) %>% 
  update_role(Yr, new_role = "ID")

# package these together into a workflow for the model,
workflow() %>%
  add_model(rf_mod) %>% 
  add_recipe(rf_recipe)  -> rf_workflow

# Tune the model
# first, create a custom grid that makes sense for my small dataset:
rf_grid <- grid_regular(
  mtry(range = c(2, 10)),        # can't exceed number of predictors
  min_n(range = c(2, 10)),       # keep well below your 26 training rows
  levels = 5
)

set.seed(123)
rf_res <- 
  rf_workflow %>% 
  tune_grid(validation_set(sablefish_split),
            grid = rf_grid,
            control = control_grid(save_pred = TRUE),
            metrics = metric_set(roc_auc))

# pull out the top models
rf_res %>% 
  show_best(metric = "roc_auc", n = 15)
# terrible!!

autoplot(rf_res)

# select the best model according to the ROC AUC metric
rf_best <- rf_res %>% 
  select_best(metric = "roc_auc")

# get the predictions for our best model
rf_auc <- 
  rf_res %>% 
  collect_predictions(parameters = rf_best) %>% 
  roc_curve(  truth = one_sd_boom_bust,
              c(".pred_Average", ".pred_Boom", ".pred_Bust")) %>% 
  mutate(model = "Random Forest")

# Plot ROC curve
autoplot(rf_auc)
# stinky stinky


### The last fit - fit the best model to combined training + validation data, evaluate on tesitng

# the last model
last_rf_mod <- rand_forest(mtry = rf_best$mtry, min_n = rf_best$min_n, trees = 1000) %>% 
  set_engine("ranger", num.threads = cores, importance = "impurity") %>% 
  set_mode("classification")

# the last workflow
last_rf_workflow <- 
  rf_workflow %>% 
  update_model(last_rf_mod)

# the last fit
last_rf_fit <- last_rf_workflow %>% 
  last_fit(sablefish_split)
# here's a fun thing - our last time period has no "bust" recruitment.

# evaluate performance
last_rf_fit %>% 
  collect_metrics()
# it's definitely better than the mlr

# look at your parameters
last_rf_fit %>% 
  extract_fit_parsnip() %>% 
  vip(num_features = 20)

#### Build the multinomial logistic regression model - version 1 - no tuning ####
mlr <- multinom_reg(penalty = 0.1)
# mlr <- multinom_reg()

# workflow() %>% 
#   add_model(mlr) -> mlr_workflow

# fit the model
mlr %>% 
  fit(one_sd_boom_bust ~., data = sablefish_other) -> mlr_sablefish_fit

# predict to the test data
bind_cols(
  predict(mlr_sablefish_fit, sablefish_test),
  predict(mlr_sablefish_fit, sablefish_test, type = "prob"),
  sablefish_test$one_sd_boom_bust
) -> mlr_sablefish_pred_testing

colnames(mlr_sablefish_pred_testing) <- c("predicted", "pred_average", "pred_boom", "pred_bust", "actual")

# compute log-loss
mn_log_loss(
  data = mlr_sablefish_pred_testing,
  truth = actual,
  c("pred_average", "pred_boom", "pred_bust")
) %>% 
  mutate(model = "mlr_sablefish") -> mlr_sablefish_log_loss




#### Tune the hyperparameters




#### Evaluate model on testing set using log-loss ####









#### Recruitment as continuous ####





#### Define functions to fit models ####

## Simple linear regression on RecDevs (stats)
fit_simple_regression <- function(data, predictor_variable, stock_name){
  formula <- as.formula(paste0("Value ~ ", predictor_variable))
  model <- lm(formula, data = data)
  as.data.frame(summary(model)$coefficients) %>% 
    mutate(model = paste0("Value ~ ", predictor_variable), stock = stock_name) -> model_fit
  
  return(model_fit)
}

## GAM on RecDevs (mgcv)
gam(Value ~ bbv_PC1, data = sablefish)


## Elastic net regression on RecDevs, to avoid p-hacking with multiple terms (glmnet)

# Create a matrix of models
sablefish_x <- model.matrix(Value ~ ., sablefish)[,-c(1:6)]
sablefish_y <- sablefish$Value

sablefish_glmnet_fit <- glmnet(sablefish_x, sablefish_y)

plot(sablefish_glmnet_fit)
print(sablefish_glmnet_fit)

sablefish_glmnet_cvfit <- cv.glmnet(sablefish_x, sablefish_y)
plot(sablefish_glmnet_cvfit)

# get the min lambda value (the value of lambda that gives the minimum mean cross-validated error)
# and get the coefficients at that value of lambda:
coef(sablefish_glmnet_cvfit, s = "lambda.min")

## One predictor multinomial logistic regression on boom/average/bust (nnet)

# parameter estimates are relative to the baseline, which is average recruitment
sablefish_multinom <- nnet:multinom(one_sd_boom_bust ~ bbv_PC1, data = sablefish)
summary(sablefish_multinom)



## One predictor multinomial logistic regression on boom/average/bust with GAM (mgcv)

# response must be numeric
sablefish$boom_bust_numeric <- as.numeric(factor(sablefish$one_sd_boom_bust, levels = c("Bust", "Average", "Boom")))-1

sablefish_cat_multinom_gam <- gam(list(boom_bust_numeric ~ s(bbv_PC1) + 
                                         s(bbv_PC2), ~s(bbv_PC1) + 
                                         s(bbv_PC2)), data = sablefish, family = mgcv::multinom(K = 2))

sablefish_cat_multinom_gam <- gam(list(boom_bust_numeric ~ s(bbv_PC1), ~s(bbv_PC1)), data = sablefish, family = mgcv::multinom(K = 2))


## Elastic net multiple predictor multinomial logistic regression on boom/average/bust (glmnet)



## Multiple predictor multinomial logistic regression on boom/average/bust, with GAMs



## Artificial Neural Network



## Random forest




#### Fit models to all stocks ####

# Sablefish
sablefish_lm_fits <- data.frame()
for (i in 1:(ncol(sablefish)-5)){
  model_fit <- fit_simple_regression(data = sablefish, predictor_variable = colnames(sablefish)[i+5], stock_name = "Sablefish")
  rbind(sablefish_lm_fits, model_fit) -> sablefish_lm_fits
}

