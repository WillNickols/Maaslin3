# Load Required Packages
for (lib in c(
  'dplyr',
  'pbapply',
  'lmerTest',
  'parallel',
  'lme4',
  'plyr'
)) {
  suppressPackageStartupMessages(require(lib, character.only = TRUE))
}

# Function to augment data for logistic fitting
augment_data <- function(formula, random_effects_formula, dat_sub){
  if (is.null(random_effects_formula)) { # No random effects
    formula <- formula(formula)
    
    # Get model matrix
    mm <- data.frame(model.matrix(formula, dat_sub), check.names = F)
    org_data_length <- nrow(mm)
    org_col_names <- colnames(mm)
    
    # Add in feature abundance
    mm$expr <- dat_sub$expr
    
    # Generate the new rows
    new_rows <- matrix(nrow = 0, ncol = ncol(mm))
    for (colname in colnames(mm)[!(colnames(mm) %in% c("(Intercept)", "expr"))]) {
      col = which(colnames(mm) == colname)
      new_rows_addition <- matrix(colMeans(mm), nrow = 4, ncol = ncol(mm), byrow = T)
      new_rows_addition[,col] <- new_rows_addition[,col] + c(-1, 1, -1, 1) * sd(mm[,col])
      new_rows_addition[,ncol(new_rows_addition)] <- c(1, 1, 0, 0)
      new_rows <- rbind(new_rows, new_rows_addition)
    }
    
    # Append these new rows to the model matrix
    mm_names <- colnames(mm)
    mm <- rbind(as.matrix(mm), new_rows)
    colnames(mm) <- mm_names
    mm <- data.frame(mm, check.names = F)
    if ("(Intercept)" %in% colnames(mm)) { mm$`(Intercept)` <- NULL }
    
    # Calculate the weights
    p <- ncol(mm) - 1
    weight_scheme <- c(rep(1, org_data_length), rep((p + 1) / (4 * p), nrow(mm) - org_data_length))
    
    # Get ready to return augmented data
    mm_input = mm
    new_formula <- paste0("expr ~ ", paste0("`", paste0(org_col_names[org_col_names != "(Intercept)"], collapse = "` + `"), "`"))
    return(list(mm_input = mm_input, weight_scheme = weight_scheme, new_formula = new_formula))
  } else { # With random effects
    input_formula <- formula(formula)
    
    # Generate the model matrix for fixed effects
    lmer_parsed <- lme4::lFormula(input_formula, dat_sub)
    mm <- data.frame(lmer_parsed$X, check.names = F)
    org_data_dim <- dim(mm)
    org_col_names <- colnames(mm)
    
    # Add on feature abundance
    mm$expr <- dat_sub$expr
    
    # Get random effects
    random_effects <- names(lmer_parsed$reTrms$cnms)
    
    # Function to duplicate a random effects variable matrix for the 4 data points added
    duplicate_rows <- function(mat) {
      mat <- as.matrix(mat)
      num_rows <- nrow(mat)
      duplicated_mat <- matrix(, ncol = ncol(mat), nrow = num_rows * 4)
      for (i in 1:num_rows) {
        duplicated_mat[((i - 1) * 4 + 1):(i * 4), ] <- matrix(rep(mat[i, ], each = 4), ncol = ncol(mat))
      }
      return(duplicated_mat)
    }
    
    # Add the new rows
    new_rows <- matrix(nrow = 0, ncol = ncol(mm) + length(random_effects))
    for (colname in colnames(mm)[!(colnames(mm) %in% c("(Intercept)", "expr"))]) {
      if (length(random_effects) > 1) { # Get the set of possible random effects combinations
        new_re <- expand.grid(lapply(dat_sub[,random_effects], unique))
      } else {
        new_re <- matrix(unique(dat_sub[,random_effects]))
      }
      col = which(colnames(mm) == colname)
      
      # Add the 4 new points per metadata and per random effect combination
      new_rows_addition <- matrix(colMeans(mm), nrow = 4 * nrow(new_re), ncol = ncol(mm), byrow = T)
      new_rows_addition[,col] <- new_rows_addition[,col] + rep(c(-1, 1, -1, 1) * sd(mm[,col]), nrow(new_re))
      new_rows_addition[,ncol(new_rows_addition)] <- rep(c(1, 1, 0, 0), nrow(new_re))
      new_rows_addition <- cbind(new_rows_addition, duplicate_rows(new_re))
      new_rows <- rbind(new_rows, new_rows_addition)
    }
    new_rows <- data.frame(new_rows)
    colnames(new_rows) <- c(colnames(mm), random_effects)
    
    # Add on other necessary columns for the random effect (e.g., when using random slopes)
    extra_cols <- colnames(dat_sub)[!(colnames(dat_sub) %in% colnames(new_rows))]
    if (length(extra_cols) > 0) {
      # Identify missed columns and build a model matrix
      mm_2 <- model.matrix(formula(paste0("expr ~ ", paste0(extra_cols, collapse = " + "))), dat_sub)
      mm_2 <- mm_2[,-1] # Remove intercept
      new_colnames <- colnames(mm_2)[!(colnames(mm_2) %in% colnames(new_rows))]
      
      # Use the model matrix to add missed columns
      if (length(new_colnames) > 0) {
        new_rows[,new_colnames] <- rep(colMeans(as.matrix(mm_2[,new_colnames])), each = nrow(new_rows))
        prev_colnames <- colnames(mm)
        mm <- cbind(mm, mm_2[,new_colnames])
        colnames(mm) <- c(prev_colnames, new_colnames)
      }
    }
    
    # Join the new rows with the original data
    mm[,random_effects] <- dat_sub[,random_effects]
    new_rows <- new_rows[,order(mapvalues(colnames(new_rows), colnames(mm), 1:ncol(mm)))]
    mm <- rbind(mm, new_rows)
    mm <- data.frame(mm, check.names = F)
    
    # Turn numeric columns to numeric
    for (colname in colnames(mm)) {
      mm[,colname] <- tryCatch({
        as.numeric(mm[,colname])
      }, warning = function(w) { 
        return(mm[,colname])
      })
    }
    
    if ("(Intercept)" %in% colnames(mm)) { mm$`(Intercept)` <- NULL }
    
    # Establish weighting scheme
    p <- ncol(mm) - 1
    weight_scheme <- c(rep(1, org_data_dim[1]), 
                       rep(org_data_dim[2] / (nrow(mm) - org_data_dim[1]), nrow(mm) - org_data_dim[1]))
    
    # Create new formula based on model matrix
    new_formula <- paste0("expr ~ ", paste0("`", paste0(org_col_names[org_col_names != "(Intercept)"], collapse = "` + `"), "`"), " + ",
                          paste0("(", paste0(unlist(findbars(formula(input_formula))), collapse = ") + ("), ")"))
    
    return(list(mm_input = mm, weight_scheme = weight_scheme, new_formula = new_formula))
  }
}

safe_deparse <- function(formula) {
  paste0(deparse(formula), collapse = "")
}

# fit the data using the model selected and applying the correction
fit.model <- function(
    features,
    metadata,
    model,
    formula = NULL,
    random_effects_formula = NULL,
    correction = "BH",
    save_models = FALSE,
    augment = FALSE,
    cores = 1) {
      
  # Check formulas are valid
  if (is.null(random_effects_formula)) {
    if (is.null(formula)) {
      logging::logerror(
        paste("Both formula and random_effects_formula are null")
      )
      stop()
    }
  }
  
  formula <- formula(formula)
  
  # Extract group components
  groups <- regmatches(safe_deparse(formula), gregexpr("group\\((.*?)\\)", safe_deparse(formula)))[[1]]
  groups <- gsub("\\)$", "", gsub("^group\\(", "", groups))
  formula <- formula(gsub("(\\s*\\+\\s*)?group\\(.*?\\)", "", safe_deparse(formula)))
  formula <- formula(gsub("~ \\+", "~", safe_deparse(formula)))
  
  #############################################################
  # Determine the function and summary for the model selected #
  #############################################################
  
  ################
  # Linear Model #
  ################
  
  if (model == "LM") {
      if (is.null(random_effects_formula)) { # Fixed effects only
          model_function <-
              function(formula, data, na.action) {
                  return(glm(formula(formula), 
                             data = data,
                             family = 'gaussian',
                             na.action = na.action))
              }
          summary_function <- function(fit) {
              lm_summary <- summary(fit)$coefficients
              if (nrow(lm_summary) < length(coef(fit))) { # If deficient rank, make sure all rownames are included
                store_names <- rownames(lm_summary)
                rows_to_add = names(coef(fit))[!(names(coef(fit)) %in% store_names)]
                lm_summary <- rbind(lm_summary, matrix(rep(NaN, 4 * length(rows_to_add)), nrow=length(rows_to_add)))
                rownames(lm_summary) <- c(store_names, rows_to_add)
              }
              para <- as.data.frame(lm_summary)[-1, -3]
              para$name <- rownames(lm_summary)[-1]
              return(para)
          }
          compare_function <- function(fit1, fit2) {
            return(anova(fit1, fit2, test = 'F')[2, 'Pr(>F)'])
          }
      } else { # Random effects
        ranef_function <- lme4::ranef
        model_function <-
              function(formula, data, na.action) {
                  return(lmerTest::lmer(
                      formula(formula), 
                      data = data, 
                      na.action = na.action))
              }
          summary_function <- function(fit) {
              lm_summary <- coef(summary(fit))
              if (nrow(lm_summary) < length(coef(fit))) { # If deficient rank, make sure all rownames are included
                store_names <- rownames(lm_summary)
                rows_to_add = names(coef(fit))[!(names(coef(fit)) %in% store_names)]
                lm_summary <- rbind(lm_summary, matrix(rep(NaN, 5 * length(rows_to_add)), nrow=length(rows_to_add)))
                rownames(lm_summary) <- c(store_names, rows_to_add)
              }
              para <- as.data.frame(lm_summary)[-1, -c(3:4)]
              para$name <- rownames(lm_summary)[-1]
              return(para)
          }
          compare_function <- function(fit1, fit2) {
            return(anova(fit1, fit2, test = 'LRT')[2, 'Pr(>Chisq)'])
          }
      }
  }

  ##################
  # Logistic Model #
  ##################

  if (model == "logistic") {
    if (is.null(random_effects_formula)) { # Fixed effects only
      if (augment) {
        model_function <- function(formula, mm, weight_scheme, na.action) {
            return(glm(
              formula = formula(formula),
              family = 'binomial',
              data = mm,
              weights = weight_scheme,
              na.action = na.action,
            ))
          }
      } else {
        model_function <-
          function(formula, data, na.action) {
            return(glm(
              formula(formula),
              data = data,
              family = 'binomial',
              na.action = na.action,
            ))
          }
      }
      summary_function <- function(fit) {
        lm_summary <- summary(fit)$coefficients
        if (nrow(lm_summary) < length(coef(fit))) { # If equal numbers of predictors and observations
          store_names <- rownames(lm_summary)
          rows_to_add = names(coef(fit))[!(names(coef(fit)) %in% store_names)]
          lm_summary <- rbind(lm_summary, matrix(rep(NaN, 4 * length(rows_to_add)), nrow=length(rows_to_add)))
          rownames(lm_summary) <- c(store_names, rows_to_add)
        }
        para <- as.data.frame(lm_summary)[-1, -3]
        para$name <- rownames(lm_summary)[-1]
        return(para)
      }
      compare_function <- function(fit1, fit2) {
        return(anova(fit1, fit2, test = 'Chisq')[2, 'Pr(>Chi)'])
      }
    } else { # Random effects
      ranef_function <- lme4::ranef
      if (augment) {
        model_function <-
          function(formula, mm, weight_scheme, na.action) {
            return(lme4::glmer(
              formula(formula), 
              data = mm, 
              family = 'binomial',
              na.action = na.action,
              weights = weight_scheme))
          }
      } else {
        model_function <-
          function(formula, data, na.action) {
            return(lme4::glmer(
              formula(formula), 
              data = data, 
              family = 'binomial',
              na.action = na.action))
          }
      }
      summary_function <- function(fit) {
        lm_summary <- coef(summary(fit))
        if (nrow(lm_summary) < length(coef(fit))) { # If deficient rank, make sure all rownames are included
          store_names <- rownames(lm_summary)
          rows_to_add = names(coef(fit))[!(names(coef(fit)) %in% store_names)]
          lm_summary <- rbind(lm_summary, matrix(rep(NaN, 4 * length(rows_to_add)), nrow=length(rows_to_add)))
          rownames(lm_summary) <- c(store_names, rows_to_add)
        }
        para <- as.data.frame(lm_summary)[-1, -3]
        para$name <- rownames(lm_summary)[-1]
        return(para)
      }
      compare_function <- function(fit1, fit2) {
        return(anova(fit1, fit2, test = 'LRT')[2, 'Pr(>Chisq)'])
      }
    }
  }
  
  #######################################
  # Init cluster for parallel computing #
  #######################################
  
  cluster <- NULL
  if (cores > 1) {
      logging::loginfo("Creating cluster of %s R processes", cores)
      cluster <- parallel::makeCluster(cores)
  }
  
  ##############################
  # Apply per-feature modeling #
  ##############################
  outputs <-
      pbapply::pblapply(seq_len(ncol(features)), cl = cluster, function(x) {
          # Extract Features One by One
          featuresVector <- features[, x]
          
          logging::loginfo(
              "Fitting model to feature number %d, %s",
              x,
              colnames(features)[x])
          
          # Make fitting matrix of features and metadata
          dat_sub <-
              data.frame(expr = as.numeric(featuresVector), metadata)
          
          # 0 or 1 observations
          if (length(unique(dat_sub$expr)) < 2) {
            output <- list()
            
            # List fixed effects that will be included
            names_to_include <- c()
            if (is.null(random_effects_formula)) { # Fixed and group effects only
              names_to_include <- colnames(model.matrix(formula(gsub("^expr ", "", safe_deparse(formula))), dat_sub))
              names_to_include <- names_to_include[names_to_include != "(Intercept)"]
            } else { # Random effects
              pattern <- paste0("\\b", paste(gsub("\\|", "\\\\|", findbars(formula(gsub("^expr ", "", safe_deparse(formula))))), collapse = "|"), "\\b")
              
              fixed_effects_only <- gsub(pattern, "", 
                                         paste0(trimws(safe_deparse(formula(gsub("^expr ", "", safe_deparse(formula))))), collapse = " "), 
                                         ignore.case = TRUE)
              fixed_effects_only <- trimws(gsub("\\(\\)", "", fixed_effects_only))
              fixed_effects_only <- trimws(gsub("\\+$|", "", fixed_effects_only))
              names_to_include <- colnames(model.matrix(formula(fixed_effects_only), dat_sub))
              names_to_include <- names_to_include[names_to_include != "(Intercept)"]
            }
            if (length(groups) > 0) {
              names_to_include <- c(names_to_include, groups)
            }
            
            output$para <-
              as.data.frame(matrix(NA,
                                   nrow = length(names_to_include), ncol = 3))
            output$para$name <- names_to_include
            
            output$residuals <- NA
            output$fitted <- NA
            if (!(is.null(random_effects_formula))) output$ranef <- NA
            output$fit <- NA
            
            colnames(output$para) <- c('coef', 'stderr' , 'pval', 'name')
            output$para$feature <- colnames(features)[x]
            output$para$error <- ifelse(model == "logistic", "All logistic values are the same",
                                        "All LM values are the same")
            return(output)
          }

          # Augment logistic fitting
          if (augment & model == "logistic" & length(unique(featuresVector)) >= 2) {
            fit_and_message <- tryCatch({
              withCallingHandlers({ # Catch non-integer # successes first
                
                if (length(groups) > 0) {
                  formula_new <- formula(paste0(c(safe_deparse(formula), groups), collapse = " + "))
                  
                  augmented_data <- augment_data(formula_new, random_effects_formula, dat_sub)
                  mm_input <- augmented_data[["mm_input"]]
                  weight_scheme <- augmented_data[["weight_scheme"]]
                  formula_new <- augmented_data[["new_formula"]]
                  
                  # Fit augmnented model
                  fit1 <-
                    model_function(
                      formula = formula_new,
                      mm = mm_input, 
                      weight_scheme = weight_scheme,
                      na.action = na.exclude)
                  
                  fit_list <- list(fit1)
                  
                  for (group in groups) {
                    # Remove terms that were added for the group
                    if (is.null(random_effects_formula)) {
                      input_formula <- formula(formula(paste0(c(safe_deparse(formula), groups[groups != group]), collapse = " + ")))
                      mm <- data.frame(model.matrix(input_formula, dat_sub), check.names = F)
                      org_col_names <- colnames(mm)
                      formula_new <- paste0("expr ~ ", paste0("`", paste0(org_col_names[org_col_names != "(Intercept)"], collapse = "` + `"), "`"))
                    } else {
                      input_formula <- formula(formula(paste0(c(safe_deparse(formula), groups[groups != group]), collapse = " + ")))
                      lmer_parsed <- lme4::lFormula(input_formula, dat_sub)
                      mm <- data.frame(lmer_parsed$X, check.names = F)
                      org_col_names <- colnames(mm)
                      # Create new formula based on model matrix
                      formula_new <- paste0("expr ~ ", paste0("`", paste0(org_col_names[org_col_names != "(Intercept)"], collapse = "` + `"), "`"), " + ",
                                            paste0("(", paste0(unlist(findbars(formula(input_formula))), collapse = ") + ("), ")"))
                    }

                    # Fit augmnented model
                    fit1 <-
                      model_function(
                        formula = formula_new,
                        mm = mm_input, 
                        weight_scheme = weight_scheme,
                        na.action = na.exclude)
                    
                    fit_list <- c(fit_list, list(fit1))
                  }
                  fit_list <- c(fit_list, NA)
                } else {
                  # Augment data
                  augmented_data <- augment_data(formula, random_effects_formula, dat_sub)
                  mm_input <- augmented_data[["mm_input"]]
                  weight_scheme <- augmented_data[["weight_scheme"]]
                  formula_new <- augmented_data[["new_formula"]]
                  
                  # Fit augmnented model
                  fit1 <-
                    model_function(
                      formula = formula_new,
                      mm = mm_input, 
                      weight_scheme = weight_scheme,
                      na.action = na.exclude)
                  
                  fit_list <- list(fit1, NA)
                }
              }, warning=function(w) {
                if (w$message == "non-integer #successes in a binomial glm!") {
                  # Still worked
                  invokeRestart("muffleWarning")
                }
                return(fit_list)
              })
            }, warning = function(w) {
              message(paste("Feature", colnames(features)[x], ":", w))
              logging::logwarn(paste(
                "Fitting problem for feature",
                x,
                "a warning was issued"))
              
              fit_list <-
                try({
                  # Data augmentation
                  if (length(groups) > 0) {
                    formula_new <- formula(paste0(c(safe_deparse(formula), groups), collapse = " + "))
                    
                    augmented_data <- augment_data(formula_new, random_effects_formula, dat_sub)
                    mm_input <- augmented_data[["mm_input"]]
                    weight_scheme <- augmented_data[["weight_scheme"]]
                    formula_new <- augmented_data[["new_formula"]]
                    
                    # Fit augmnented model
                    fit1 <-
                      model_function(
                        formula = formula_new,
                        mm = mm_input, 
                        weight_scheme = weight_scheme,
                        na.action = na.exclude)
                    
                    fit_list <- list(fit1)
                    
                    for (group in groups) {
                      if (is.null(random_effects_formula)) {
                        input_formula <- formula(formula(paste0(c(safe_deparse(formula), groups[groups != group]), collapse = " + ")))
                        mm <- data.frame(model.matrix(input_formula, dat_sub), check.names = F)
                        org_col_names <- colnames(mm)
                        formula_new <- paste0("expr ~ ", paste0("`", paste0(org_col_names[org_col_names != "(Intercept)"], collapse = "` + `"), "`"))
                      } else {
                        input_formula <- formula(formula(paste0(c(safe_deparse(formula), groups[groups != group]), collapse = " + ")))
                        lmer_parsed <- lme4::lFormula(input_formula, dat_sub)
                        mm <- data.frame(lmer_parsed$X, check.names = F)
                        org_col_names <- colnames(mm)
                        # Create new formula based on model matrix
                        formula_new <- paste0("expr ~ ", paste0("`", paste0(org_col_names[org_col_names != "(Intercept)"], collapse = "` + `"), "`"), " + ",
                                              paste0("(", paste0(unlist(findbars(formula(input_formula))), collapse = ") + ("), ")"))
                      }
                      
                      # Fit augmnented model
                      fit1 <-
                        model_function(
                          formula = formula_new,
                          mm = mm_input, 
                          weight_scheme = weight_scheme,
                          na.action = na.exclude)
                      
                      fit_list <- c(fit_list, list(fit1))
                    }
                    fit_list
                  } else {
                    # Augment data
                    augmented_data <- augment_data(formula, random_effects_formula, dat_sub)
                    mm_input <- augmented_data[["mm_input"]]
                    weight_scheme <- augmented_data[["weight_scheme"]]
                    formula_new <- augmented_data[["new_formula"]]
                    
                    # Fit augmnented model
                    fit1 <-
                      model_function(
                        formula = formula_new,
                        mm = mm_input, 
                        weight_scheme = weight_scheme,
                        na.action = na.exclude)
                    
                    list(fit1)
                  }
                })
              if (inherits(fit_list, "try-error")) { # Sometimes warning process triggers error that gives unlisted result
                return(c(list(fit_list), list(w$message)))
              }
              return(c(fit_list, list(w$message)))
            }, error = function(err) { # Warn on augmented ataa
              fit_list <-
                try({
                  # Data augmentation
                  if (length(groups) > 0) {
                    formula_new <- formula(paste0(c(safe_deparse(formula), groups), collapse = " + "))
                    
                    augmented_data <- augment_data(formula_new, random_effects_formula, dat_sub)
                    mm_input <- augmented_data[["mm_input"]]
                    weight_scheme <- augmented_data[["weight_scheme"]]
                    formula_new <- augmented_data[["new_formula"]]
                    
                    # Fit augmnented model
                    fit1 <-
                      model_function(
                        formula = formula_new,
                        mm = mm_input, 
                        weight_scheme = weight_scheme,
                        na.action = na.exclude)
                    
                    fit_list <- list(fit1)
                    
                    for (group in groups) {
                      if (is.null(random_effects_formula)) {
                        input_formula <- formula(formula(paste0(c(safe_deparse(formula), groups[groups != group]), collapse = " + ")))
                        mm <- data.frame(model.matrix(input_formula, dat_sub), check.names = F)
                        org_col_names <- colnames(mm)
                        formula_new <- paste0("expr ~ ", paste0("`", paste0(org_col_names[org_col_names != "(Intercept)"], collapse = "` + `"), "`"))
                      } else {
                        input_formula <- formula(formula(paste0(c(safe_deparse(formula), groups[groups != group]), collapse = " + ")))
                        lmer_parsed <- lme4::lFormula(input_formula, dat_sub)
                        mm <- data.frame(lmer_parsed$X, check.names = F)
                        org_col_names <- colnames(mm)
                        # Create new formula based on model matrix
                        formula_new <- paste0("expr ~ ", paste0("`", paste0(org_col_names[org_col_names != "(Intercept)"], collapse = "` + `"), "`"), " + ",
                                              paste0("(", paste0(unlist(findbars(formula(input_formula))), collapse = ") + ("), ")"))
                      }
                      
                      # Fit augmnented model
                      fit1 <-
                        model_function(
                          formula = formula_new,
                          mm = mm_input, 
                          weight_scheme = weight_scheme,
                          na.action = na.exclude)
                      
                      fit_list <- c(fit_list, list(fit1))
                    }
                    fit_list
                  } else {
                    # Augment data
                    augmented_data <- augment_data(formula, random_effects_formula, dat_sub)
                    mm_input <- augmented_data[["mm_input"]]
                    weight_scheme <- augmented_data[["weight_scheme"]]
                    formula_new <- augmented_data[["new_formula"]]
                    
                    # Fit augmnented model
                    fit1 <-
                      model_function(
                        formula = formula_new,
                        mm = mm_input, 
                        weight_scheme = weight_scheme,
                        na.action = na.exclude)
                    
                    list(fit1)
                  }
                })
              return(c(list(fit_list), list(err$message)))
            })
          } else { # LM or non-augmented logistic
            fit_and_message <- tryCatch({
              if (length(groups) > 0) {
                formula_new <- formula(paste0(c(safe_deparse(formula), groups), collapse = " + "))
                
                fit1 <-
                  model_function(
                    formula_new, 
                    data = dat_sub, 
                    na.action = na.exclude)
                
                fit_list <- list(fit1)
                
                for (group in groups) {
                  formula_new <- formula(paste0(c(safe_deparse(formula), groups[groups != group]), collapse = " + "))
                  
                  fit1 <-
                    model_function(
                      formula_new, 
                      data = dat_sub, 
                      na.action = na.exclude)
                  
                  fit_list <- c(fit_list, list(fit1))
                }
                fit_list <- c(fit_list, NA)
              } else {
                fit1 <-
                  model_function(
                    formula, 
                    data = dat_sub, 
                    na.action = na.exclude)
                fit_list <- list(fit1, NA)
              }
            }, warning = function(w) {
              message(paste("Feature", colnames(features)[x], ":", w))
              logging::logwarn(paste(
                "Fitting problem for feature", 
                x, 
                "a warning was issued"))
              
              fit_list <-
                try({
                  if (length(groups) > 0) {
                    formula_new <- formula(paste0(c(safe_deparse(formula), groups), collapse = " + "))
                    
                    fit1 <-
                      model_function(
                        formula_new, 
                        data = dat_sub, 
                        na.action = na.exclude)
                    
                    fit_list <- list(fit1)
                    
                    for (group in groups) {
                      formula_new <- formula(paste0(c(safe_deparse(formula), groups[groups != group]), collapse = " + "))
                      
                      fit1 <-
                        model_function(
                          formula_new, 
                          data = dat_sub, 
                          na.action = na.exclude)
                      
                      fit_list <- c(fit_list, list(fit1))
                    }
                    fit_list
                  } else {
                    fit1 <-
                      model_function(
                        formula, 
                        data = dat_sub, 
                        na.action = na.exclude)
                    list(fit1)
                  }
                })
              if (inherits(fit_list, "try-error")) { # Sometimes warning process triggers error that gives unlisted result
                return(c(list(fit_list), list(w$message)))
              }
              return(c(fit_list, list(w$message)))
            }, error = function(err) {
              fit_list <-
                try({
                  if (length(groups) > 0) {
                    formula_new <- formula(paste0(c(safe_deparse(formula), groups), collapse = " + "))
                    
                    fit1 <-
                      model_function(
                        formula_new, 
                        data = dat_sub, 
                        na.action = na.exclude)
                    
                    fit_list <- list(fit1)
                    
                    for (group in groups) {
                      formula_new <- formula(paste0(c(safe_deparse(formula), groups[groups != group]), collapse = " + "))
                      
                      fit1 <-
                        model_function(
                          formula_new, 
                          data = dat_sub, 
                          na.action = na.exclude)
                      
                      fit_list <- c(fit_list, list(fit1))
                    }
                    fit_list
                  } else {
                    fit1 <-
                      model_function(
                        formula, 
                        data = dat_sub, 
                        na.action = na.exclude)
                    list(fit1)
                  }
                })
              return(c(list(fit_list), list(err$message)))
            })
          }
          
          fit <- fit_and_message[[1]]
          
          # Gather Output
          output <- list()
          if (all(!inherits(fit, "try-error"))) {
            output$para <- summary_function(fit)
            if (length(groups) > 0) {
              i <- 2
              for (group in groups) {
                if (all(!inherits(fit_and_message[[i]], "try-error"))) {
                  output$para <- tryCatch({
                    rbind(output$para, list(NA, NA, compare_function(fit_and_message[[i]], fit_and_message[[1]]), group))
                  },
                  warning = function(w) {
                    rbind(output$para, list(NA, NA, NA, group))
                  },
                  error = function(err) {
                    rbind(output$para, list(NA, NA, NA, group))
                  })
                  rownames(output$para) <- c(rownames(output$para)[-nrow(output$para)], group)
                  
                } else {
                  output$para <- rbind(output$para, list(NA, NA, NA, NA))
                }
                
                i <- i + 1
              }
            }
            
            if (is.null(random_effects_formula)) { # If no random effects, the missing names are the model matrix names
              names_to_include <- colnames(model.matrix(formula(gsub("^expr ", "", safe_deparse(formula))), dat_sub))
              names_to_include <- names_to_include[names_to_include != "(Intercept)"]
            } else { # If there are random effects, get the fixed effects model matrix names
              # Random effects to remove
              pattern <- paste0("\\b", paste(gsub("\\|", "\\\\|", findbars(formula(gsub("^expr ", "", safe_deparse(formula))))), collapse = "|"), "\\b")
              
              # Fixed effects piece of formula
              fixed_effects_only <- gsub(pattern, "", 
                                         paste0(trimws(safe_deparse(formula(gsub("^expr ", "", safe_deparse(formula))))), collapse = " "), 
                                         ignore.case = TRUE)
              fixed_effects_only <- trimws(gsub("\\(\\)", "", fixed_effects_only))
              fixed_effects_only <- trimws(gsub("\\+$|", "", fixed_effects_only))
              
              # Fixed effects terms
              names_to_include <- colnames(model.matrix(formula(fixed_effects_only), dat_sub))
              names_to_include <- names_to_include[names_to_include != "(Intercept)"]
            }
            if (length(groups) > 0) {
              names_to_include <- c(names_to_include, groups)
            }
            
            if (any(!(names_to_include %in% rownames(output$para)))) {
              fit_properly <- FALSE
              fit_and_message[[2]] <- "Metadata dropped during fitting (rank deficient)"
            } else { # No errors, summaries are correct
              fit_properly <- TRUE
            }
            
            output$para <- output$para[rownames(output$para) %in% names_to_include,]
          } else { # Fit issue occurred
              fit_properly <- FALSE
          }
          
          if (fit_properly) {
            output$residuals <- residuals(fit)
            output$fitted <- fitted(fit)
            if (!(is.null(random_effects_formula))) {
              # Returns a list with a table for each random effect
              l <- ranef_function(fit)
              
              # Rename rows as random effect labels if only random intercepts
              if (length(l) == 1 & ncol(l[[1]]) == 1 & colnames(l[[1]])[1] == "(Intercept)") {
                d<-as.vector(unlist(l))
                names(d)<-unlist(lapply(l, row.names))
                d[setdiff(unique(metadata[,names(l)]), names(d))] <- NA
                d <- d[order(unique(metadata[,names(l)]))]
                output$ranef<-d
              } else {
                # Otherwise return the random effects list
                output$ranef <- l
              }
            }
            if (save_models) {
              output$fit <- fit
            } else {
              output$fit <- NA
            }
          } else {
            logging::logwarn(paste(
              "Fitting problem for feature", 
              x, 
              "returning NA"))
            
            names_to_include <- c()
            
            if (is.null(random_effects_formula)) { # If no random effects, the missing names are the model matrix names
              names_to_include <- colnames(model.matrix(formula(gsub("^expr ", "", safe_deparse(formula))), dat_sub))
              names_to_include <- names_to_include[names_to_include != "(Intercept)"]
            } else { # If there are random effects, get the fixed effects model matrix names
              # Random effects to remove
              pattern <- paste0("\\b", paste(gsub("\\|", "\\\\|", findbars(formula(gsub("^expr ", "", safe_deparse(formula))))), collapse = "|"), "\\b")
              
              # Fixed effects piece of formula
              fixed_effects_only <- gsub(pattern, "", 
                                         paste0(trimws(safe_deparse(formula(gsub("^expr ", "", safe_deparse(formula))))), collapse = " "), 
                                         ignore.case = TRUE)
              fixed_effects_only <- trimws(gsub("\\(\\)", "", fixed_effects_only))
              fixed_effects_only <- trimws(gsub("\\+$|", "", fixed_effects_only))
              
              # Fixed effects terms
              names_to_include <- colnames(model.matrix(formula(fixed_effects_only), dat_sub))
              names_to_include <- names_to_include[names_to_include != "(Intercept)"]
            }
            if (length(groups) > 0) {
              names_to_include <- c(names_to_include, groups)
            }
            
            # Store NA values for missing outputs
            output$para <-
              as.data.frame(matrix(NA,
                                   nrow = length(names_to_include), ncol = 3))
            output$para$name <- names_to_include
            
            output$residuals <- NA
            output$fitted <- NA
            if (!(is.null(random_effects_formula))) output$ranef <- NA
            output$fit <- NA
          }

          colnames(output$para) <- c('coef', 'stderr' , 'pval', 'name')
          output$para$feature <- colnames(features)[x]
          output$para$error <- fit_and_message[[length(fit_and_message)]]
          
          return(output)
          })
  
  
  # stop the cluster
  if (!is.null(cluster))
      parallel::stopCluster(cluster)
  
  # bind the results for each feature
  paras <-
      do.call(rbind, lapply(outputs, function(x) {
          return(x$para)
      }))
  residuals <-
      do.call(rbind, lapply(outputs, function(x) {
          return(x$residuals)
      }))
  
  row.names(residuals) <- colnames(features)

  fitted <-
    do.call(rbind, lapply(outputs, function(x) {
      return(x$fitted)
    }))
  row.names(fitted) <- colnames(features)   
  
  fits <-
    lapply(outputs, function(x) {
      return(x$fit)
    })
  names(fits) <- colnames(features)  
  
  # Return NULL rather than empty object if fits aren't saved
  if (all(is.na(fits))) {
    fits <- NULL
  }
  
  if (!(is.null(random_effects_formula))) {
    ranef <-
      do.call(rbind, lapply(outputs, function(x) {
        return(x$ranef)
      }))
    row.names(ranef) <- colnames(features) 
  }
    
  ################################
  # Apply correction to p-values #
  ################################
  
  paras$qval <- as.numeric(p.adjust(paras$pval, method = correction))
  
  #####################################################
  # Determine the metadata names from the model names #
  #####################################################
  
  metadata_names <- colnames(metadata)
  # order the metadata names by decreasing length
  metadata_names_ordered <-
      metadata_names[order(
          nchar(metadata_names), decreasing = TRUE)]
  # find the metadata name based on the match 
  # to the beginning of the string
  extract_metadata_name <- function(name) {
    tmp_val <- metadata_names_ordered[mapply(
      startsWith, 
      name, 
      metadata_names_ordered)][1]
    if (is.na(tmp_val)) {
      metadata_names_ordered[mapply(
        grepl, 
        metadata_names_ordered,
        name)][1]
    } else {
      return(tmp_val)
    }
  }
  paras$metadata <- unlist(lapply(paras$name, extract_metadata_name))
  # compute the value as the model contrast minus metadata
  paras$value <-
      mapply(function(x, y) {
          if (x == y)
              x
          else
              gsub("^\\:", "", gsub(x, "", y))
      }, paras$metadata, paras$name)
  
  ##############################
  # Sort by decreasing q-value #
  ##############################
  
  paras <- paras[order(paras$qval, decreasing = FALSE), ]
  paras <-
      dplyr::select(
          paras,
          c('feature', 'metadata', 'value'),
          dplyr::everything())
  paras$model <- model
  rownames(paras)<-NULL
  
  if (!(is.null(random_effects_formula))) {
    return(list("results" = paras, "residuals" = residuals, "fitted" = fitted, "ranef" = ranef, "fits" = fits))
  } else {
    return(list("results" = paras, "residuals" = residuals, "fitted" = fitted, "ranef" = NULL, "fits" = fits))
  }
}

# Get joint significance for zeros and non-zeros
add_joint_signif <- function(fit_data_non_zero, fit_data_binary, analysis_method, correction) {
  # Subset to shared columns
  fit_data_binary_signif <- fit_data_binary$results[,c("feature", "metadata", "value", "name", "pval", "error")]
  colnames(fit_data_binary_signif) <- c("feature", "metadata", "value", "name", "logistic", "logistic_error")
  fit_data_non_zero_signif <- fit_data_non_zero$results[,c("feature", "metadata", "value", "name", "pval", "error")]
  colnames(fit_data_non_zero_signif) <- c("feature", "metadata", "value", "name", analysis_method, "LM_error")
  
  # Join and check linear and logistic pieces
  merged_signif <- full_join(fit_data_binary_signif, fit_data_non_zero_signif, by=c("feature", "metadata", "value", "name"))
  if (nrow(merged_signif) != nrow(fit_data_binary_signif) | nrow(merged_signif) != nrow(fit_data_non_zero_signif)) {
    print(nrow(merged_signif))
    print(nrow(fit_data_binary_signif))
    print(nrow(fit_data_non_zero_signif))
    print(anti_join(fit_data_binary_signif, fit_data_non_zero_signif, by=c("feature", "metadata")))
    print(anti_join(fit_data_non_zero_signif, fit_data_binary_signif, by=c("feature", "metadata")))
    stop("Merged significance tables have different associations")
  }
  
  # Create a combined p-value
  merged_signif$pval_joint <- pbeta(pmin(merged_signif[,analysis_method], 
                                         merged_signif[,"logistic"]), 
                                    1, 2)
  
  # If NA or model errored, use the p-value of the non-NA
  merged_signif$pval_joint <- ifelse(is.na(merged_signif[,analysis_method]) | !is.na(merged_signif$LM_error), merged_signif[,"logistic"], merged_signif$pval_joint)
  merged_signif$pval_joint <- ifelse(is.na(merged_signif[,"logistic"]) | !is.na(merged_signif$logistic_error), merged_signif[,analysis_method], merged_signif$pval_joint)
  merged_signif$pval_joint <- ifelse((is.na(merged_signif[,"logistic"]) | !is.na(merged_signif$logistic_error)) & 
                                       (is.na(merged_signif[,analysis_method]) | !is.na(merged_signif$LM_error)), 
                                     NA, merged_signif$pval_joint)
  merged_signif$qval_joint <- as.numeric(p.adjust(merged_signif$pval_joint, method = correction))
  return(list(append_joint(fit_data_non_zero, merged_signif), append_joint(fit_data_binary, merged_signif)))
}

# Take logistic or LM component and add on the merged significance pieces
append_joint <- function(outputs, merged_signif) {
  merged_signif <- merged_signif[,c("feature", "metadata", "value", "name", "pval_joint", "qval_joint")]
  tmp_colnames <- colnames(outputs$results)
  tmp_colnames <- case_when(tmp_colnames == "pval" ~ "pval_single",
                            tmp_colnames == "qval" ~ "qval_single",
                            TRUE ~ tmp_colnames)
  colnames(outputs$results) <- tmp_colnames
  
  merged_signif <- merge(outputs$results, merged_signif, 
                         by=c("feature", "metadata", "value", "name"))
  
  merged_signif <- merged_signif[order(merged_signif$qval_joint),]
  
  return(merged_signif)
}

bind_and_reorder <- function(growing_resid_mat, outputs_tmp_residuals, params_and_data_and_formula_data) {
  tmp_out <- rbind(growing_resid_mat, outputs_tmp_residuals)
  tmp_out <- 
    tmp_out[order(as.numeric(mapvalues(rownames(tmp_out),
                                       colnames(params_and_data_and_formula_data),
                                       1:ncol(params_and_data_and_formula_data), warn_missing = F))),]
  return(tmp_out)
}
