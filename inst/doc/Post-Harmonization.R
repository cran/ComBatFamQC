## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, eval = FALSE------------------------------------------------------
#  library(ComBatFamQC)
#  data(age_df)

## ----eval = FALSE-------------------------------------------------------------
#  age_df <- data.frame(age_df)
#  features <- colnames(age_df)[c(6:56)]
#  age <- "age"
#  sex <- "sex"
#  icv <- "ICV_baseline"
#  age_df[[sex]] <- as.factor(age_df[[sex]])

## ----eval = FALSE-------------------------------------------------------------
#  # Create sub_df for different features
#  sub_df_list <- lapply(seq_len(length(features)), function(i){
#      sub_df <- age_df[,c(features[i], age, sex, icv)] %>% na.omit()
#      colnames(sub_df) <- c(features[i], "age", "sex", "icv")
#      return(sub_df)
#    })

## ----eval = FALSE-------------------------------------------------------------
#  # For MAC users
#  library(parallel)
#  age_list <- mclapply(seq_len(length(features)), function(w){
#    age_sub <- age_list_gen (sub_df = sub_df_list[[w]],  lq = 0.25, hq = 0.75)
#    return(age_sub)
#  }, mc.cores = detectCores())
#  
#  # For Windows users
#  age_list <- mclapply(1:length(features), function(w){
#    age_sub <- age_list_gen (sub_df = sub_df_list[[w]],  lq = 0.25, hq = 0.75)
#    return(age_sub)
#  }, mc.cores = 1)
#  
#  names(age_list) <- features
#  
#  quantile_type <- c(paste0("quantile_", 100*0.25), "median", paste0("quantile_", 100*0.75))

## ----eval=FALSE---------------------------------------------------------------
#  # plotly: interactive plot
#  ComBatFamQC::age_shiny(age_list, features, quantile_type, use_plotly = TRUE)
#  # ggplot: static plot
#  ComBatFamQC::age_shiny(age_list, features, quantile_type, use_plotly = FALSE)

## ----eval=FALSE---------------------------------------------------------------
#  # Save age trend table
#  temp_dir <- tempfile()
#  dir.create(temp_dir)
#  age_save(path = temp_dir, age_list = age_list)
#  
#  # Save GAMLSS Model
#  gamlss_model <- lapply(seq_len(length(age_list)), function(i){
#          g_model <- age_list[[i]]$model
#          return(g_model)})
#  names(gamlss_model) <- names(age_list)
#  saveRDS(gamlss_model, file = file.path(temp_dir, "gamlss_model.rds"))

## ----eval=FALSE---------------------------------------------------------------
#  features <- colnames(adni)[c(43:104)]
#  covariates <- c("timedays", "AGE", "SEX", "DIAGNOSIS")
#  interaction <- c("timedays,DIAGNOSIS")
#  batch <- "manufac"
#  combat_model <- combat_harm(type = "lm", features = features, batch = batch, covariates = covariates, interaction = interaction, smooth = NULL, random = NULL, df = adni)
#  harmonized_df <- combat_model$harmonized_df

## ----eval=FALSE---------------------------------------------------------------
#  # generate residuals by removing timedays and DIAGNOSIS effects, while preserving AGE and SEX effects.
#  result_residual <- residual_gen(type = "lm", features = features, covariates = covariates, interaction = interaction, smooth = NULL, df = harmonized_df, rm = c("timedays", "DIAGNOSIS"))
#  
#  # save residual data set
#  write.csv(result_residual$residual, file.path(temp_dir, "residual.csv"))
#  
#  # save regression model
#  saveRDS(result_residual$model, file.path(temp_dir, "regression_model.rds"))

## ----eval=FALSE---------------------------------------------------------------
#  result_residual <- residual_gen(df = harmonized_df, rm = c("timedays", "DIAGNOSIS"), model = TRUE, model_path = file.path(temp_dir, "regression_model.rds"))
#  
#  # save residual data set
#  write.csv(result_residual$residual, file.path(temp_dir, "residual.csv"))
#  # Clean up the temporary file
#  unlink(temp_dir, recursive = TRUE)

