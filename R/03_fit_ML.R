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
library(brulee)

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

# get number of cores for parallelization
# parallelize
cores <- parallel::detectCores()

# drop all columns that aren't ID, features, or outcome
sablefish %>% 
  dplyr::select(-c("one_sd_boom_bust", "Parm_StDev", "recr_scaled")) -> sablefish

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
  initial = 10,      # number of initial years
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

# create a table that has the truth (Value) for the left out years
lfo_df <- lfo_resamples %>%
  mutate(
    assessment_data = map(splits, assessment),
    .row = map_int(splits, ~ as.integer(assessment(.x)$Yr - min(sablefish$Yr) + 1))
  ) %>%
  dplyr::select(id, assessment_data, .row) %>%
  unnest(assessment_data) %>%
  dplyr::select(id, Value, .row)


#### Fit a null model as a baseline ####

# Use same null model as Ward et al. 2024:
# The null model consisted of using the  historical mean preceding each holdout 
# point (equivalent to fitting  an intercept only model).
# this will be very close to assuming a recDev of zero, but not exactly

# Calculate class proportions from training data
null_model_pred <- data.frame(Yr = last_10_yrs, .pred_Value = rep(NA, 10))
for (i in last_10_yrs){
  pred_val <- mean(subset(sablefish, !(Yr %in% i:max(sablefish$Yr)))$Value)
  null_model_pred[null_model_pred$Yr == i,]$.pred_Value <- pred_val
}

# Construct a null prediction tibble using training proportions
null_model_pred <- cbind(lfo_df, null_model_pred)

# Calculate null RMSE + rsq
null_model_pred %>%
  rmse(truth = Value,
              .pred_Value) %>% 
  mutate(model = "null") -> null_rmse

null_model_pred %>%
  rsq(truth = Value,
      .pred_Value) %>% 
  mutate(model = "null") -> null_rsq

# start a table to compare models with the null metrics table
model_comparison_table <- bind_rows(null_rmse, null_rsq)


# fit each of our candidate models

#### Fit multiple linear regression ####

### Generate all possible combinations of 1, 2, and 3 predictors
# this will be used for linear regressions + GAMs

predictors <- colnames(dplyr::select(sablefish, -c("Yr", "Value")))

# use the map function
predictor_combos <- bind_rows(map(1:3, ~ combn(predictors, .x, simplify = FALSE) %>% 
  map(~ tibble(predictors = list(.x))))
) %>% 
  mutate(model_id = row_number())

# how many models?
nrow(predictor_combos)

  
## Function to fit a lm model for a given set of predictors
fit_one_lm_lfo <- function(pred_set, resamples){
  # Build the formula
  lm_f <- as.formula(paste0("Value ~ ", paste(pred_set[[1]], collapse = " + ")))
  
  # Set up the model
  lm_mod <- linear_reg() %>% 
    set_engine("lm", num.threads = cores)
  
  # Set up the recipe
  lm_recipe <- recipe(lm_f, data = subset(sablefish, !(Yr %in% last_10_yrs)))

  # package these together into a workflow for the model
  workflow() %>%
    add_model(lm_mod) %>% 
    add_recipe(lm_recipe)  -> lm_workflow
  
  # Fit and evaluate using LFO-CV
  # note that this is unable to compute rsq on a single left-out value,
  # so we will calculate our metrics after the function is run
  fit_resamples(
    lm_workflow,
    resamples = resamples,
    metrics = metric_set(rmse),
    control = control_resamples(save_pred = TRUE)
  ) -> lm_res
  
  return(lm_res)
  
  # fit_resamples(
  #   lm_workflow,
  #   resamples = resamples,
  #   metrics = metric_set(rmse),
  #   control = control_resamples(save_pred = TRUE)
  # ) -> lm_res
  #
  # lm_res %>% 
  #   collect_predictions() %>% 
  #   rsq(truth = Value, estimate = .pred) -> lm_rsq
  # 
  # lm_res %>% 
  #   collect_predictions() %>% 
  #   rmse(truth = Value, estimate = .pred) -> lm_rmse
  # 
  # # combine these all into one df and return
  # bind_rows(lm_rsq, lm_rmse) %>% 
  #   mutate(predictors = as.character(lm_f[3])) -> lm_metrics
  # 
  # return(lm_metrics)
}

# Fit all models
predictor_combos %>%
  mutate(lfo_cv_fit = map(predictor_combos$predictors, ~fit_one_lm_lfo(.x, assessment_resamples))) -> lm_results

lm_results %>% 
  mutate(
    # RMSE from fold-level evaluation (already computed per fold)
    # this is calculated as an average for each held out year
    rmse = map_dbl(lfo_cv_fit, ~ collect_metrics(.x) %>% 
                     filter(.metric == "rmse") %>% 
                     pull(mean)),
    # rsq computed across all pooled held-out predictions
    # this is calculated
    rsq = map_dbl(lfo_cv_fit, ~ collect_predictions(.x) %>% 
                    rsq(truth = Value, estimate = .pred) %>% 
                    pull(.estimate)),
    n_predictors = map_int(predictors, length),
    predictor_names = map_chr(predictors, ~ paste(.x, collapse = " + "))
  ) %>% 
  dplyr::select(model_id, n_predictors, predictor_names, rmse, rsq) -> lm_metrics

# extract the best model based on RMSE
lm_metrics %>% 
  arrange(rmse) %>% 
  head(1) -> best_fit_lm

### The last fit - fit the best model to combined training + validation data, evaluate on testing

# Build the last
lm_last_f <- as.formula(paste0("Value ~ ", paste(best_fit_lm$predictor_names, collapse = " + ")))

# the last recipe
lm_last_recipe <- recipe(lm_last_f, data = sablefish)

# the last workflow
# update the recipe with the new formula
last_lm_workflow <- 
  lm_workflow %>% 
  update_recipe(lm_last_recipe)

# the last fit
# first, we need to make a split that has all training data and all testing data separated
sablefish_split <- make_splits(x = list(analysis = 1:(nrow(sablefish)-10),
                                        assessment = (nrow(sablefish)-9):nrow(sablefish)),
                               data = sablefish)
last_lm_fit <- last_lm_workflow %>% 
  last_fit(sablefish_split)

# evaluate performance
last_lm_fit %>% 
  collect_metrics() %>% 
  mutate(model = paste0("lm(", deparse(lm_last_f), ")")) %>%
  dplyr::select(-.config) -> last_lm_fit_metrics

# add these to our table
model_comparison_table %>% 
  bind_rows(last_lm_fit_metrics) -> model_comparison_table


#### Fit GAM ####

## Function to fit a gam model for a given set of predictors
fit_one_gam_lfo <- function(pred_set, resamples){
  # Build the formulas - two in this case for GAMs
  gam_f_splines <- as.formula(paste0("Value ~ ", paste("s(", pred_set[[1]], ", k = 3)", collapse = " + ")))
  
  gam_f <- as.formula(paste0("Value ~ ", paste(pred_set[[1]], collapse = " + ")))

  # Set up the model
  gam_mod <- gen_additive_mod() %>% 
    set_engine("mgcv", num.threads = cores) %>% 
    set_mode("regression")
  
  # package these together into a workflow for the model
  workflow() %>%
    add_model(gam_mod, formula = gam_f_splines) %>% 
    add_formula(gam_f) %>% 
    fit(data = subset(sablefish, !(Yr %in% last_10_yrs)))  -> gam_workflow
  
  # Fit and evaluate using LFO-CV
  # note that this is unable to compute rsq on a single left-out value,
  # so we will calculate our metrics after the function is run
  fit_resamples(
    gam_workflow,
    resamples = resamples,
    metrics = metric_set(rmse),
    control = control_resamples(save_pred = TRUE)
  ) -> gam_res
  
  return(gam_res)
  
  # fit_resamples(
  #   gam_workflow,
  #   resamples = resamples,
  #   metrics = metric_set(rmse),
  #   control = control_resamples(save_pred = TRUE)
  # ) -> gam_res
  #
  # gam_res %>% 
  #   collect_predictions() %>% 
  #   rsq(truth = Value, estimate = .pred) -> gam_rsq
  # 
  # gam_res %>% 
  #   collect_predictions() %>% 
  #   rmse(truth = Value, estimate = .pred) -> gam_rmse
  # 
  # # combine these all into one df and return
  # bind_rows(gam_rsq, gam_rmse) %>% 
  #   mutate(predictors = as.character(gam_f[3])) -> gam_metrics
  # 
  # return(gam_metrics)
}

# Fit all models
predictor_combos %>%
  mutate(lfo_cv_fit = map(predictor_combos$predictors, ~fit_one_gam_lfo(.x, assessment_resamples))) -> gam_results

gam_results %>% 
  mutate(
    # RMSE from fold-level evaluation (already computed per fold)
    # this is calculated as an average for each held out year
    rmse = map_dbl(lfo_cv_fit, ~ collect_metrics(.x) %>% 
                     filter(.metric == "rmse") %>% 
                     pull(mean)),
    # rsq computed across all pooled held-out predictions
    # this is calculated
    rsq = map_dbl(lfo_cv_fit, ~ collect_predictions(.x) %>% 
                    rsq(truth = Value, estimate = .pred) %>% 
                    pull(.estimate)),
    n_predictors = map_int(predictors, length),
    predictor_names = map_chr(predictors, ~ paste(.x, collapse = " + "))
  ) %>% 
  dplyr::select(model_id, n_predictors, predictor_names, rmse, rsq) -> gam_metrics

# extract the best model based on RMSE
gam_metrics %>% 
  arrange(rmse) %>% 
  head(1) -> best_fit_gam

### The last fit - fit the best model to combined training + validation data, evaluate on testing

# Build the last model
gam_last_f_splines <- as.formula(paste0("Value ~ ", paste("s(", best_fit_gam$predictor_names, ", k = 3)", collapse = " + ")))
gam_last_f <- as.formula(paste0("Value ~ ", paste(best_fit_gam$predictor_names, collapse = " + ")))


# the last workflow
# update the workflow with the new formula
last_gam_workflow <- 
  gam_workflow %>% 
  update_model(gam_mod, formula = gam_last_f_splines) %>% 
  update_formula(gam_last_f)

# the last fit
# first, we need to make a split that has all training data and all testing data separated
sablefish_split <- make_splits(x = list(analysis = 1:(nrow(sablefish)-10),
                                        assessment = (nrow(sablefish)-9):nrow(sablefish)),
                               data = sablefish)
last_gam_fit <- last_gam_workflow %>% 
  last_fit(sablefish_split)

# evaluate performance
last_gam_fit %>% 
  collect_metrics() %>% 
  mutate(model = paste0("gam(", deparse(gam_last_f), ")")) %>%
  dplyr::select(-.config) -> last_gam_fit_metrics

# add these to our table
model_comparison_table %>% 
  bind_rows(last_gam_fit_metrics) -> model_comparison_table
  
  
  #### Fit random forest ####
  # parallelize
  cores <- parallel::detectCores()
  
  # set up the model
  rf_mod <- 
    rand_forest(mtry = 2, min_n = 4, trees = 1000) %>% 
    set_engine("ranger", num.threads = cores) %>% 
    set_mode("regression")
  
  ## testing zone
  # set up the model
  rf_mod <- 
    rand_forest(mtry = tune(), min_n = tune(), trees = tune()) %>% 
    set_engine("ranger", num.threads = cores) %>% 
    set_mode("regression")
  
  # set up the recipe - formula + code columns
  rf_recipe <- recipe(Value ~., data = subset(sablefish, !(Yr %in% last_10_yrs))) %>% 
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
              metrics = metric_set(rmse))
  
  # select the best model according to RMSE
  rf_best <- rf_res %>% 
    select_best(metric = "rmse")
  
  autoplot(rf_res)
  
  # evaluate model on the held-out testing data and store rmse
  ### The last fit - fit the best model to combined training + validation data, evaluate on testing
  
  # the last model
  last_rf_mod <- rand_forest(mtry = rf_best$mtry, min_n = rf_best$min_n, trees = rf_best$trees) %>% 
    set_engine("ranger", num.threads = cores, importance = "impurity") %>% 
    set_mode("regression")
  
  # the last workflow
  last_rf_workflow <- 
    rf_workflow %>% 
    update_model(last_rf_mod)
  
  # the last fit
  # first, we need to make a split that has all training data and all testing data separated
  sablefish_split <- make_splits(x = list(analysis = 1:(nrow(sablefish)-10),
                                          assessment = (nrow(sablefish)-9):nrow(sablefish)),
                                          data = sablefish)
  last_rf_fit <- last_rf_workflow %>% 
    last_fit(sablefish_split)
  
  # evaluate performance
  last_rf_fit %>% 
    collect_metrics() %>% 
    mutate(model = "random forest") %>% 
    dplyr::select(-.config) -> last_rf_fit_metrics
  
  # add these to our table
  model_comparison_table %>% 
    bind_rows(last_rf_fit_metrics) -> model_comparison_table
  
  #### random forest v2 - try this again with a different grid ####
  # set up the model
  rf2_mod <- 
    rand_forest(mtry = tune(), min_n = tune(), trees = tune()) %>% 
    set_engine("ranger", num.threads = cores) %>% 
    set_mode("regression")
  
  # set up the recipe - formula + code columns
  rf2_recipe <- recipe(Value ~., data = subset(sablefish, !(Yr %in% last_10_yrs))) %>% 
    update_role(Yr, new_role = "ID")
  
  # package these together into a workflow for the model,
  workflow() %>%
    add_model(rf2_mod) %>% 
    add_recipe(rf2_recipe)  -> rf2_workflow
  
  rf2_grid <- grid_regular(
    mtry(range = c(2, 10)),        # can't exceed number of predictors
    min_n(range = c(2, 10)), # keep well below number of training rows
    trees(range = c(500, 2500)),
    # levels = c(5, 5, 25) # try more trees than mtry or min_n - takes forever
    levels = c(5, 5, 5)
  )
  
  # tune the hyperparameters on the assessment resamples
  rf2_res <- 
    rf2_workflow %>% 
    tune_grid(assessment_resamples,
              grid = rf2_grid,
              control = control_grid(save_pred = TRUE),
              metrics = metric_set(rmse))
  
  # select the best model according to RMSE
  rf2_best <- rf2_res %>% 
    select_best(metric = "rmse")
  
  autoplot(rf2_res)
  
  # evaluate model on the held-out testing data and store rmse
  ### The last fit - fit the best model to combined training + validation data, evaluate on testing
  
  # the last model
  last_rf2_mod <- rand_forest(mtry = rf2_best$mtry, min_n = rf2_best$min_n, trees = rf2_best$trees) %>% 
    set_engine("ranger", num.threads = cores, importance = "impurity") %>% 
    set_mode("regression")
  
  # the last workflow
  last_rf2_workflow <- 
    rf2_workflow %>% 
    update_model(last_rf2_mod)
  
  # the last fit
  # first, we need to make a split that has all training data and all testing data separated
  sablefish_split <- make_splits(x = list(analysis = 1:(nrow(sablefish)-10),
                                          assessment = (nrow(sablefish)-9):nrow(sablefish)),
                                 data = sablefish)
  last_rf2_fit <- last_rf2_workflow %>% 
    last_fit(sablefish_split)
  
  # evaluate perf2ormance
  last_rf2_fit %>% 
    collect_metrics() %>% 
    mutate(model = "random forest 2") %>% 
    dplyr::select(-.config) -> last_rf2_fit_metrics
  
  # add these to our table
  model_comparison_table %>% 
    bind_rows(last_rf2_fit_metrics) -> model_comparison_table
  
  # look at your parameters
  last_rf2_fit_metrics %>% 
    extract_fit_parsnip() %>% 
    vip(num_features = 20)
  
#### Fit neural network ####
  
  # create combinations of 3 for initial testing
  NN_predictor_combos <- bind_rows(map(3, ~ combn(predictors, .x, simplify = FALSE) %>% 
                                         map(~ tibble(predictors = list(.x))))
  ) %>% 
    mutate(model_id = row_number())
  
  
  ## Function to fit a nn model for a given set of predictors
  fit_one_nn_lfo <- function(pred_set, resamples){
    
    # Build the formula
    nn_f <- as.formula(paste0("Value ~ ", paste(pred_set[[1]], collapse = " + ")))
    
    # Set up the recipe
    nn_recipe <- recipe(nn_f, data = subset(sablefish, !(Yr %in% last_10_yrs)))
    
    # set model specifications
    nn_spec <- # mlp(epochs = 100, hidden_units = 2, penalty = 0.01, learn_rate = 0.1) %>% 
      mlp(epochs = 100, penalty = 0) %>% 
      # set_engine("brulee", validation = 0) %>% 
      set_engine("nnet") %>% 
      set_mode("regression")
    
    # set up the workflow
    nn_workflow <- nn_recipe %>% 
      workflow(nn_spec)
    
    # Fit and evaluate using LFO-CV
    # note that this is unable to compute rsq on a single left-out value,
    # so we will calculate our metrics after the function is run
    # mlp_reg_fit <- nn_spec |> fit(nn_f, data = subset(sablefish, !(Yr %in% last_10_yrs)))
    # 3-5-1 network with 26 weights
    
    fit_resamples(
      nn_workflow,
      resamples = resamples,
      metrics = metric_set(rmse),
      control = control_resamples(save_pred = TRUE)
    ) -> nn_res
    
    return(nn_res)
  }
  
  # Fit all models
  NN_predictor_combos %>%
    mutate(lfo_cv_fit = map(NN_predictor_combos$predictors, ~fit_one_nn_lfo(.x, assessment_resamples))) -> nn_results
  
  nn_results %>% 
    mutate(
      # RMSE from fold-level evaluation (already computed per fold)
      # this is calculated as an average for each held out year
      rmse = map_dbl(lfo_cv_fit, ~ collect_metrics(.x) %>% 
                       filter(.metric == "rmse") %>% 
                       pull(mean)),
      # rsq computed across all pooled held-out predictions
      # this is calculated
      rsq = map_dbl(lfo_cv_fit, ~ collect_predictions(.x) %>% 
                      rsq(truth = Value, estimate = .pred) %>% 
                      pull(.estimate)),
      n_predictors = map_int(predictors, length),
      predictor_names = map_chr(predictors, ~ paste(.x, collapse = " + "))
    ) %>% 
    dplyr::select(model_id, n_predictors, predictor_names, rmse, rsq) -> nn_metrics
  
  # extract the best model based on RMSE
  nn_metrics %>% 
    arrange(rmse) %>% 
    head(1) -> best_fit_nn
  
  best_fit_nn <- subset(nn_metrics, model_id == 469)
  
  ### The last fit - fit the best model to combined training + validation data, evaluate on testing
  
  # Build the last
  nn_last_f <- as.formula(paste0("Value ~ ", paste(best_fit_nn$predictor_names, collapse = " + ")))
  
  # the last recipe
  nn_last_recipe <- recipe(nn_last_f, data = sablefish)
  
  # the last workflow
  # update the recipe with the new formula
  last_nn_workflow <- 
    nn_workflow %>% 
    update_recipe(nn_last_recipe)
  
  # the last fit
  # first, we need to make a split that has all training data and all testing data separated
  sablefish_split <- make_splits(x = list(analysis = 1:(nrow(sablefish)-10),
                                          assessment = (nrow(sablefish)-9):nrow(sablefish)),
                                 data = sablefish)
  last_nn_fit <- last_nn_workflow %>% 
    last_fit(sablefish_split)
  
  # evaluate performance
  last_nn_fit %>% 
    collect_metrics() %>% 
    mutate(model = paste0("nn(", deparse(nn_last_f), ")")) %>%
    dplyr::select(-.config) -> last_nn_fit_metrics
  
  # add these to our table
  model_comparison_table %>% 
    bind_rows(last_nn_fit_metrics) -> model_comparison_table
  
  
  
  