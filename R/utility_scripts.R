# Load Required Packages
for (lib in c('vegan', 'chemometrics', 'car', 'metagenomeSeq', 'edgeR')) {
    suppressPackageStartupMessages(require(lib, character.only = TRUE))
}


###################
## Transformation #
###################

transformFeatures = function(features, transformation) {
    if (transformation == 'LOG')     {
        features <- apply(features, 2, LOG)
    }
    
    if (transformation == 'LOGIT')     {
        features <- apply(features, 2, LOGIT)
    }
    
    if (transformation == 'AST')     {
        features <- apply(features, 2, AST)
    }
    
    return(features)
}


##################
## Normalization #
##################

normalizeFeatures = function(features, normalization) {
    if (normalization == 'TSS')
    {
        features <- TSSnorm(features)
    }
    
    if (normalization == 'CLR')
    {
        features <- CLRnorm(features)
    }
    
    if (normalization == 'CSS')
    {
        features <- CSSnorm(features)
    }
    
    if (normalization == 'TMM')
    {
        features <- TMMnorm(features)
    }
    
    if (normalization == 'NONE')
    {
        features <- features
    }
    
    return(features)
}

######################
## TSS Normalization #
######################

# Apply TSS Normalization To A Dataset

TSSnorm = function(features) {
    # Convert to Matrix from Data Frame
    features_norm = as.matrix(features)
    dd <- colnames(features_norm)
    
    # TSS Normalizing the Data
    features_TSS <-
        vegan::decostand(
            features_norm,
            method = "total",
            MARGIN = 1,
            na.rm = TRUE)
    
    # Convert back to data frame
    features_TSS <- as.data.frame(features_TSS)
    
    # Rename the True Positive Features - Same Format as Before
    colnames(features_TSS) <- dd
    
    
    # Return
    return(features_TSS)
}


######################
## CLR Normalization #
######################

# Apply CLR Normalization To A Dataset

CLRnorm = function(features) {
    # Convert to Matrix from Data Frame
    features_norm = as.matrix(features)
    dd <- colnames(features_norm)
    
    # CLR Normalizing the Data
    features_CLR <- chemometrics::clr(features_norm + 1)
    
    # Convert back to data frame
    features_CLR <- as.data.frame(features_CLR)
    
    # Rename the True Positive Features - Same Format as Before
    colnames(features_CLR) <- dd
    
    
    # Return
    return(features_CLR)
}

######################
## CSS Normalization #
######################

<<<<<<< Updated upstream
# Apply CSS Normalization To A Dataset
=======
######################
# From metagenomeSeq #
######################

calcNormFactors <- function(x, p = cumNormStat(x)) {
  xx <- x
  xx[xx == 0] <- NA
  qs = matrixStats::colQuantiles(xx, probs = p, na.rm = TRUE)
  normFactors <- sapply(1:ncol(xx), function(i) {
    xx = (x[, i] - .Machine$double.eps)
    sum(xx[xx <= qs[i]])
  })
  names(normFactors) <- colnames(x)
  as.data.frame(normFactors)
}

cumNormStat <- function (counts, qFlag = TRUE, pFlag = FALSE, rel = 0.1, ...) {
  mat = counts
  if (any(colSums(mat) == 0)) 
    stop("Warning empty sample")
  smat = sapply(1:ncol(mat), function(i) {
    sort(mat[, i], decreasing = FALSE)
  })
  ref = rowMeans(smat)
  yy = mat
  yy[yy == 0] = NA
  ncols = ncol(mat)
  refS = sort(ref)
  k = which(refS > 0)[1]
  lo = (length(refS) - k + 1)
  if (qFlag == TRUE) {
    diffr = sapply(1:ncols, function(i) {
      refS[k:length(refS)] - quantile(yy[, i], p = seq(0, 
                                                       1, length.out = lo), na.rm = TRUE)
    })
  }
  if (qFlag == FALSE) {
    diffr = sapply(1:ncols, function(i) {
      refS[k:length(refS)] - approx(sort(yy[, i], decreasing = FALSE), 
                                    n = lo)$y
    })
  }
  diffr2 = matrixStats::rowMedians(abs(diffr), na.rm = TRUE)
  if (pFlag == TRUE) {
    plot(abs(diff(diffr2[diffr2 > 0]))/diffr2[diffr2 > 0][-1], 
         type = "h", ylab = "Relative difference for reference", 
         xaxt = "n", ...)
    abline(h = rel)
    axis(1, at = seq(0, length(diffr2), length.out = 5), 
         labels = seq(0, 1, length.out = 5))
  }
  x = which(abs(diff(diffr2))/diffr2[-1] > rel)[1]/length(diffr2)
  if (x <= 0.5) {
    message("Default value being used.")
    x = 0.5
  }
  return(x)
}

cumNormMat <- function (x, p = cumNormStat(x), sl = 1000) {
  xx <- x
  xx[xx == 0] <- NA
  qs = matrixStats::colQuantiles(xx, probs = p, na.rm = TRUE)
  newMat <- sapply(1:ncol(xx), function(i) {
    xx = (x[, i] - .Machine$double.eps)
    sum(xx[xx <= qs[i]])
  })
  nmat <- sweep(x, 2, newMat/sl, "/")
  return(nmat)
}

MRcounts <- function (counts, norm_factors, sl = 1000)  {
  if (any(is.na(norm_factors))) {
    x = cumNormMat(as.matrix(counts), sl = sl)
  } else {
    x = sweep(as.matrix(counts), 2, norm_factors/sl, "/")
  }
  return(x)
}
>>>>>>> Stashed changes

CSSnorm = function(features) {
    # Convert to Matrix from Data Frame
    features_norm = as.matrix(features)
    dd <- colnames(features_norm)
    
    # CSS Normalizing the Data
    # Create the metagenomeSeq object
    MGS = metagenomeSeq::newMRexperiment(
        t(features_norm),
        featureData = NULL,
        libSize = NULL,
        normFactors = NULL
    )
    # Trigger metagenomeSeq to calculate its Cumulative Sum scaling factor.
    MGS = metagenomeSeq::cumNorm(MGS, p = metagenomeSeq::cumNormStat(MGS))
    # Save the normalized data as data.frame
    features_CSS = as.data.frame(t(
        metagenomeSeq::MRcounts(MGS, norm = TRUE, log = FALSE)))
    
    # Rename the True Positive Features - Same Format as Before
    colnames(features_CSS) <- dd
    
    
    # Return as list
    return(features_CSS)
}

######################
## TMM Normalization #
######################

# Apply TMM Normalization To A Dataset

TMMnorm = function(features) {
    # Convert to Matrix from Data Frame
    features_norm = as.matrix(features)
    dd <- colnames(features_norm)
    
    # TMM Normalizing the Data
    X <- t(features_norm)
    
    libSize = edgeR::calcNormFactors(X, method = "TMM")
    eff.lib.size = colSums(X) * libSize
    
    ref.lib.size = mean(eff.lib.size)
    #Use the mean of the effective library sizes as a reference library size
    X.output = sweep(X, MARGIN = 2, eff.lib.size, "/") * ref.lib.size
    #Normalized read counts
    
    # Convert back to data frame
    features_TMM <- as.data.frame(t(X.output))
    
    # Rename the True Positive Features - Same Format as Before
    colnames(features_TMM) <- dd
    
    
    # Return as list
    return(features_TMM)
}

#######################################
# Arc-Sine Square Root Transformation #
#######################################

# Arc Sine Square Root Transformation
AST <- function(x) {
    y <- sign(x) * asin(sqrt(abs(x)))
    if(any(is.na(y))) {
        logging::logerror(
            paste0("AST transform is only valid for values between -1 and 1. ",
                   "Please select an appropriate normalization option or ",
                   "normalize your data prior to running."))
        stop()
    }
    return(y)
}

########################
# Logit Transformation #
########################

# Zero-inflated Logit Transformation (Does not work well for microbiome data)
LOGIT <- function(x) {
    y <- car::logit(x, adjust = 0)
    y[!is.finite(y)] <- 0
    return(y)
}

# Shifted Logit Transformation (Lukens et al, 2014, Nature)
# LOGIT_S<-function(x){
#     y<-0.5*log(x/(1-x)) + 10
#     y[!is.finite(y)]<-0
#     return(y)
# }

######################
# Log Transformation #
######################

# Log Transformation
LOG <- function(x) {
    y <- replace(x, x == 0, min(x[x>0]) / 2)
    return(log2(y))
}
<<<<<<< Updated upstream
=======

############################
# Write out the model fits #
############################

write_fits <- function(params_data_formula_fit) {
  param_list <- maaslin_parse_param_list(params_data_formula_fit[["param_list"]])
  output <- param_list[["output"]]
  fit_data <- params_data_formula_fit[["fit_data"]]

  fits_folder <- file.path(output, "fits")
  if (!file.exists(fits_folder)) {
    print("Creating output fits folder")
    dir.create(fits_folder)
  }
  
  ################################
  # Write out the raw model fits #
  ################################
  
  if (param_list[["save_models"]]) {
    model_file = file.path(fits_folder, "models.rds")
    # remove models file if already exists (since models append)
    if (file.exists(model_file)) {
      logging::logwarn(
        "Deleting existing model objects file: %s", model_file)
      unlink(model_file)
    }
    logging::loginfo("Writing model objects to file %s", model_file)
    saveRDS(fit_data$fits, file = model_file)   
  }
  
  ###########################
  # Write residuals to file #
  ###########################
  
  residuals_file = file.path(fits_folder, "residuals.rds")
  # remove residuals file if already exists (since residuals append)
  if (file.exists(residuals_file)) {
    logging::logwarn(
      "Deleting existing residuals file: %s", residuals_file)
    unlink(residuals_file)
  }
  logging::loginfo("Writing residuals to file %s", residuals_file)
  saveRDS(fit_data$residuals, file = residuals_file)
  
  ###############################
  # Write fitted values to file #
  ###############################
  
  fitted_file = file.path(fits_folder, "fitted.rds")
  # remove fitted file if already exists (since fitted append)
  if (file.exists(fitted_file)) {
    logging::logwarn(
      "Deleting existing fitted file: %s", fitted_file)
    unlink(fitted_file)
  }
  logging::loginfo("Writing fitted values to file %s", fitted_file)
  saveRDS(fit_data$fitted, file = fitted_file)
  
  #########################################################
  # Write extracted random effects to file (if specified) #
  #########################################################
  
  if (!is.null(param_list[["random_effects"]])) {
    ranef_file = file.path(fits_folder, "ranef.rds")
    # remove ranef file if already exists (since ranef append)
    if (file.exists(ranef_file)) {
      logging::logwarn(
        "Deleting existing ranef file: %s", ranef_file)
      unlink(ranef_file)
    }
    logging::loginfo("Writing extracted random effects to file %s", ranef_file)
    saveRDS(fit_data$ranef, file = ranef_file)
  }
}

write_results <- function(params_data_formula_fit) {
  param_list <- maaslin_parse_param_list(params_data_formula_fit[["param_list"]])
  output <- param_list[["output"]]
  max_significance <- param_list[["max_significance"]]
  fit_data <- params_data_formula_fit[["fit_data"]]

  #############################
  # Write all results to file #
  #############################
  
  results_file <- file.path(output, "all_results.tsv")
  logging::loginfo(
    "Writing all results to file (ordered by increasing q-values): %s",
    results_file)
  ordered_results <- fit_data$results[order(fit_data$results$qval), ]
  # Remove any that are NA for the q-value
  ordered_results <-
    ordered_results[!is.na(ordered_results$qval), ]
  write.table(
    ordered_results[c(
      "feature",
      "metadata",
      "value",
      "name",
      "coef",
      "stderr",
      "N",
      "N.not.zero",
      "pval",
      "qval")],
    file = results_file,
    sep = "\t",
    quote = FALSE,
    col.names = c(
      "feature",
      "metadata",
      "name",
      "value",
      "coef",
      "stderr",
      "N",
      "N.not.0",
      "pval",
      "qval"
    ),
    row.names = FALSE
  )
  
  ###########################################
  # Write results passing threshold to file #
  ###########################################
  
  significant_results <-
    ordered_results[ordered_results$qval <= max_significance, ]
  significant_results_file <-
    file.path(output, "significant_results.tsv")
  logging::loginfo(
    paste("Writing the significant results",
          "(those which are less than or equal to the threshold",
          "of %f ) to file (ordered by increasing q-values): %s"),
    max_significance,
    significant_results_file
  )
  write.table(
    significant_results[c(
      "feature",
      "metadata",
      "value",
      "name",
      "coef",
      "stderr",
      "N",
      "N.not.zero",
      "pval",
      "qval")],
    file = significant_results_file,
    sep = "\t",
    quote = FALSE,
    col.names = c(
      "feature",
      "metadata",
      "value",
      "name",
      "coef",
      "stderr",
      "N",
      "N.not.0",
      "pval",
      "qval"
    ),
    row.names = FALSE
  )
}



>>>>>>> Stashed changes
