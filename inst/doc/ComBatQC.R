## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, eval=FALSE--------------------------------------------------------
#  library(ComBatFamQC)
#  library(dplyr)
#  data(adni)

## ----eval=FALSE---------------------------------------------------------------
#  features <- colnames(adni)[c(43:104)]
#  covariates <- c("timedays", "AGE", "SEX", "DIAGNOSIS")
#  interaction <- c("timedays,DIAGNOSIS")
#  batch <- "manufac"
#  result_orig <- visual_prep(type = "lm", features = features, batch = batch, covariates = covariates, interaction = interaction, smooth = NULL, random = NULL, df = adni)
#  comfam_shiny(result_orig)

## ----eval=FALSE---------------------------------------------------------------
#  result_gam <- visual_prep(type = "gam", features = features, batch = batch, covariates = covariates, interaction = interaction, smooth_int_type = "linear", smooth = "AGE", df = adni)
#  comfam_shiny(result_gam)

## ----eval=FALSE---------------------------------------------------------------
#  result_lmer <- visual_prep(type = "lmer", features = features, batch = batch, covariates = covariates, interaction = interaction, smooth = NULL, random = "subid", df = adni)
#  comfam_shiny(result_lmer)

## ----eval=FALSE---------------------------------------------------------------
#  #library(quarto)
#  temp_dir <- tempfile()
#  dir.create(temp_dir)
#  diag_save(path = temp_dir, result = result_lmer, use_quarto = TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  diag_save(path = temp_dir, result = result_lmer, use_quarto = FALSE)

## ----eval=FALSE---------------------------------------------------------------
#  features <- colnames(adni)[c(43:104)]
#  covariates <- c("timedays", "AGE", "SEX", "DIAGNOSIS")
#  interaction <- c("timedays,DIAGNOSIS")
#  batch <- "manufac"
#  ## Harmonize using evaluation results as the inputs
#  combat_model <- combat_harm(result = result_orig, type = "lm", interaction = interaction, smooth = NULL, random = NULL, df = adni)
#  ## Harmonize through specifying features, batch, covariates and df arguments
#  combat_model_copy <- combat_harm(type = "lm", features = features, batch = batch, covariates = covariates, interaction = interaction, smooth = NULL, random = NULL, df = adni)
#  ## Expect to get the same harmonization results
#  identical(combat_model$harmonized_df, combat_model_copy$harmonized_df)
#  
#  # save harmonized data
#  write.csv(combat_model$harmonized_df, file.path(temp_dir, "harmonized.csv"))
#  
#  # save combat model
#  saveRDS(combat_model$combat.object, file.path(temp_dir, "combat_model.rds"))
#  # Clean up the temporary file
#  unlink(temp_dir, recursive = TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  ## Harmonize using evaluation results as the inputs
#  combat_model_lmer <- combat_harm(result = result_lmer, type = "lmer", interaction = interaction, smooth = NULL, random = "subid", df = adni)
#  ## Harmonize through specifying features, batch, covariates and df arguments
#  combat_model_lmer_copy <- combat_harm(type = "lmer", features = features, batch = batch, covariates = covariates, interaction = interaction, smooth = NULL, random = "subid", df = adni)
#  ## Expect to get the same harmonization results
#  identical(combat_model_lmer$harmonized_df, combat_model_lmer_copy$harmonized_df)

## ----eval=FALSE---------------------------------------------------------------
#  ## Harmonize using evaluation results as the inputs
#  combat_model_gam <- combat_harm(result = result_gam, type = "gam", interaction = interaction, smooth = "AGE", smooth_int_type = "linear", df = adni)
#  ## Harmonize through specifying features, batch, covariates and df arguments
#  combat_model_gam_copy <- combat_harm(type = "gam", features = features, batch = batch, covariates = covariates, interaction = interaction, smooth = "AGE", smooth_int_type = "linear", df = adni)
#  ## Expect to get the same harmonization results
#  identical(combat_model_gam$harmonized_df, combat_model_gam_copy$harmonized_df)

## ----eval=FALSE---------------------------------------------------------------
#  ## Harmonize using evaluation results as the inputs
#  covbat_model <- combat_harm(result = result_gam, type = "gam", interaction = interaction, smooth = "AGE", smooth_int_type = "linear", df = adni, family = "covfam")
#  ## Harmonize through specifying features, batch, covariates and df arguments
#  covbat_model_copy <- combat_harm(type = "gam", features = features, batch = batch, covariates = covariates, interaction = interaction, smooth_int_type = "linear", smooth = "AGE", df = adni, family = "covfam")
#  ## Expect to get the same harmonization results
#  identical(covbat_model$harmonized_df, covbat_model_copy$harmonized_df)

## ----eval=FALSE---------------------------------------------------------------
#  saved_model <- combat_model_gam$combat.object
#  harm_predict <- combat_harm(df = adni %>% head(1000), predict = TRUE, object = saved_model)

## ----eval=FALSE---------------------------------------------------------------
#  # harmonize reference data
#  reference_site <- adni %>% group_by(site) %>% summarize(count = n()) %>% arrange(desc(count)) %>% pull(site) %>% head(30)
#  reference_df <- adni %>% filter(site %in% reference_site)
#  features <- colnames(reference_df)[c(43:104)]
#  covariates <- c("timedays", "AGE", "SEX", "DIAGNOSIS")
#  interaction <- c("timedays,DIAGNOSIS")
#  batch <- "site"
#  ref_model <- combat_harm(type = "lmer", features = features, batch = batch, covariates = covariates, interaction = interaction, smooth = NULL, random = "subid", df = reference_df)
#  
#  # harmonize new data to the reference data
#  harm_new <- combat_harm(type = "lmer", features = features, batch = batch, covariates = covariates, interaction = interaction, smooth = NULL, random = "subid", df = adni, reference = ref_model$harmonized_df)

