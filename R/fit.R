# Load Required Packages
for (lib in c(
  'dplyr',
  'pbapply',
  'lmerTest',
  'parallel',
  'lme4'
)) {
  suppressPackageStartupMessages(require(lib, character.only = TRUE))
}

# fit the data using the model selected and applying the correction
fit.model <-
    function(
        features,
        metadata,
        model,
        formula = NULL,
        random_effects_formula = NULL,
        correction = "BH",
        save_models = FALSE,
        cores = 1) {
      
        if (is.null(random_effects_formula)) {
          if (is.null(formula)) {
            logging::logerror(
              paste("Both formula and random_effects_formula are null")
            )
            stop()
          }
        }

        #############################################################
        # Determine the function and summary for the model selected #
        #############################################################
        
        ################
        # Linear Model #
        ################
        
        if (model == "LM") {
            if (is.null(random_effects_formula)) {
                model_function <-
                    function(formula, data, na.action) {
                        return(glm(
                            formula,
                            data = data,
                            family = 'gaussian',
                            na.action = na.action
                        ))
                    }
                summary_function <- function(fit) {
                    lm_summary <- summary(fit)$coefficients
                    para <- as.data.frame(lm_summary)[-1, -3]
                    para$name <- rownames(lm_summary)[-1]
                    return(para)
                }
            } else {
              ranef_function <- lme4::ranef
              model_function <-
                    function(formula, data, na.action) {
                        return(lmerTest::lmer(
                            formula, 
                            data = data, 
                            na.action = na.action))
                    }
                summary_function <- function(fit) {
                    lm_summary <- coef(summary(fit))
                    para <- as.data.frame(lm_summary)[-1, -c(3:4)]
                    para$name <- rownames(lm_summary)[-1]
                    return(para)
                }
            }
        }
      
        ##################
        # Logistic Model #
        ##################
      
        if (model == "logistic") {
          if (is.null(random_effects_formula)) {
            model_function <-
              function(formula, data, na.action) {
                return(glm(
                  formula,
                  data = data,
                  family = 'binomial',
                  na.action = na.action
                ))
              }
            summary_function <- function(fit) {
              lm_summary <- summary(fit)$coefficients
              para <- as.data.frame(lm_summary)[-1, -3]
              para$name <- rownames(lm_summary)[-1]
              return(para)
            }
          } else {
            ranef_function <- lme4::ranef
            model_function <-
              function(formula, data, na.action) {
                return(lme4::glmer(
                  formula, 
                  data = data, 
                  family = 'binomial',
                  na.action = na.action))
              }
            summary_function <- function(fit) {
              lm_summary <- coef(summary(fit))
              para <- as.data.frame(lm_summary)[-1, -3]
              para$name <- rownames(lm_summary)[-1]
              return(para)
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
                
                # Fit Model
                logging::loginfo(
                    "Fitting model to feature number %d, %s",
                    x,
                    colnames(features)[x])
                
                dat_sub <-
                    data.frame(expr = as.numeric(featuresVector), metadata)
                
                fit_and_message <- tryCatch({
                    fit1 <-
                        model_function(
                            formula, 
                            data = dat_sub, 
                            na.action = na.exclude)
                    fit1 <- list(fit1, NA)
                }, warning = function(w) { 
                  message(paste("Feature", colnames(features)[x], ":", w))
                  logging::logwarn(paste(
                    "Fitting problem for feature", 
                    x, 
                    "a warning was issued"))
                  fit1 <-
                    try({
                      model_function(
                        formula, 
                        data = dat_sub, 
                        na.action = na.exclude)
                    })
                  return(list(fit1, w$message))
                  }, error = function(err) {
                    fit1 <-
                        try({
                            model_function(
                                formula, 
                                data = dat_sub, 
                                na.action = na.exclude)
                        })
                    return(list(fit1, err$message))
                })
                
                fit <- fit_and_message[[1]]
                
                # Gather Output
                output <- list()
                if (all(!inherits(fit, "try-error"))) {
                    output$para <- summary_function(fit)
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
                    output$para <-
                        as.data.frame(matrix(NA, 
                            nrow = ncol(metadata), ncol = 3))
                    output$para$name <- colnames(metadata)
                    output$residuals <- NA
                    output$fitted <- NA
                    if (!(is.null(random_effects_formula))) output$ranef <- NA
                    output$fit <- NA
                }

                colnames(output$para) <- c('coef', 'stderr' , 'pval', 'name')
                output$para$feature <- colnames(features)[x]
                output$para$error <- fit_and_message[[2]]
                
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
  
# Get joint significance for zeros and non-zeros
add_joint_signif <- function(fit_data_non_zero, fit_data_binary, analysis_method) {
  fit_data_binary_signif <- fit_data_binary$results[,c("feature", "metadata", "value", "name", "pval")]
  colnames(fit_data_binary_signif) <- c("feature", "metadata", "value", "name", analysis_method)
  fit_data_non_zero_signif <- fit_data_non_zero$results[,c("feature", "metadata", "value", "name", "pval")]
  colnames(fit_data_non_zero_signif) <- c("feature", "metadata", "value", "name", "logistic")
  merged_signif <- merge(fit_data_binary_signif, fit_data_non_zero_signif, by=c("feature", "metadata", "value", "name"))
  merged_signif$pval_joint <- pbeta(pmin(merged_signif[,analysis_method], 
                                         merged_signif[,"logistic"]), 
                                    1, 2)
  merged_signif$qval_joint <- as.numeric(p.adjust(merged_signif$pval_joint, method = correction))
  return(list(append_joint(fit_data_non_zero, merged_signif), append_joint(fit_data_binary, merged_signif)))
}



