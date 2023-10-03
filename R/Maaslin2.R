#!/usr/bin/env Rscript

###############################################################################
# MaAsLin2

# Copyright (c) 2018 Harvard School of Public Health

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
###############################################################################

# load in the required libraries, report an error if they are not installed

for (lib in c('optparse', 'logging', 'data.table', 'dplyr')) {
    suppressPackageStartupMessages(require(lib, character.only = TRUE))
}

###############################################################
# If running on the command line, load other Maaslin2 modules #
###############################################################

# this evaluates to true if script is being called directly as an executable
if (identical(environment(), globalenv()) &&
        !length(grep("^source\\(", sys.calls()))) {
    # source all R in Maaslin2 package, relative to this folder
    # same method as original maaslin
    script_options <- commandArgs(trailingOnly = FALSE)
    script_path <-
        sub("--file=", "", script_options[grep("--file=", script_options)])
    script_dir <- dirname(script_path)
    script_name <- basename(script_path)
    
    for (R_file in dir(script_dir, pattern = "*.R$"))
    {
        if (!(R_file == script_name))
            source(file.path(script_dir, R_file))
    }
}

###########################
# Set the default options #
###########################

normalization_choices <- c("TSS", "CLR", "CSS", "NONE", "TMM")
analysis_method_choices_names <-
    c("LM", "CPLM", "NEGBIN", "ZINB")
transform_choices <- c("LOG", "LOGIT", "AST", "NONE")
valid_choice_method_norm <- hash::hash()
valid_choice_method_norm[[analysis_method_choices_names[3]]] <-
    normalization_choices[3:5]
valid_choice_method_norm[[analysis_method_choices_names[4]]] <-
    normalization_choices[3:5]
valid_choice_method_transform <- analysis_method_choices_names[1:2]
valid_choice_transform_norm <- hash::hash()
valid_choice_transform_norm[[transform_choices[2]]] <-
    normalization_choices[c(1, 4)]
valid_choice_transform_norm[[transform_choices[3]]] <-
    normalization_choices[c(1, 4)]
correction_choices <-
    c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY")

# set the default run options
args <- list()
args$input_data <- NULL
args$input_metadata <- NULL
args$output <- NULL
args$min_abundance <- 0.0
args$min_prevalence <- 0.1
args$min_variance <- 0.0
args$max_significance <- 0.25
args$normalization <- normalization_choices[1]
args$transform <- transform_choices[1]
args$analysis_method <- analysis_method_choices_names[1]
args$random_effects <- NULL
args$fixed_effects <- NULL
args$formula <- NULL
args$correction <- correction_choices[1]
args$standardize <- TRUE
args$plot_heatmap <- TRUE
args$heatmap_first_n <- 50
args$plot_scatter <- TRUE
args$max_pngs <- 10
args$save_scatter <- FALSE
args$cores <- 1
args$save_models <- FALSE
args$reference <- NULL

##############################
# Add command line arguments #
##############################

options <-
    optparse::OptionParser(usage = paste(
        "%prog [options]",
        " <data.tsv> ",
        "<metadata.tsv> ",
        "<output_folder>" 
        )
    )
options <-
    optparse::add_option(
        options,
        c("-a", "--min_abundance"),
        type = "double",
        dest = "min_abundance",
        default = args$min_abundance,
        help = paste0("The minimum abundance for each feature",
            "[ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-p", "--min_prevalence"),
        type = "double",
        dest = "min_prevalence",
        default = args$min_prevalence,
        help = paste0("The minimum percent of samples for which",
            "a feature is detected at minimum abundance",
            "[ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-b", "--min_variance"),
        type = "double",
        dest = "min_variance",
        default = args$min_variance,
        help = paste0("Keep features with variances",
            "greater than value",
            "[ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-s", "--max_significance"),
        type = "double",
        dest = "max_significance",
        default = args$max_significance,
        help = paste0("The q-value threshold for significance",
            "[ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-n", "--normalization"),
        type = "character",
        dest = "normalization",
        default = args$normalization,
        help = paste(
            "The normalization method to apply",
            "[ Default: %default ] [ Choices:",
            toString(normalization_choices),
            "]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-t", "--transform"),
        type = "character",
        dest = "transform",
        default = args$transform,
        help = paste(
            "The transform to apply [ Default: %default ] [ Choices:",
            toString(transform_choices),
            "]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-m", "--analysis_method"),
        type = "character",
        dest = "analysis_method",
        default = args$analysis_method,
        help = paste(
            "The analysis method to apply [ Default: %default ] [ Choices:",
            toString(analysis_method_choices_names),
            "]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-r", "--random_effects"),
        type = "character",
        dest = "random_effects",
        default = args$random_effects,
        help = paste("The random effects for the model, ",
            "comma-delimited for multiple effects",
            "[ Default: none ]"
        )
    )
options <-
<<<<<<< Updated upstream
    optparse::add_option(
        options,
        c("-f", "--fixed_effects"),
        type = "character",
        dest = "fixed_effects",
        default = args$fixed_effects,
        help = paste("The fixed effects for the model,",
            "comma-delimited for multiple effects",
            "[ Default: all ]"
        )
=======
  optparse::add_option(
    options,
    c("--formula"),
    type = "character",
    dest = "formula",
    default = args$formula,
    help = paste("The formula for the model,",
                 "[ Default: all variables fixed ]"
    )
  )
options <-
  optparse::add_option(
    options,
    c("-c", "--correction"),
    type = "character",
    dest = "correction",
    default = args$correction,
    help = paste("The correction method for computing",
                 "the q-value [ Default: %default ]"
>>>>>>> Stashed changes
    )
options <-
    optparse::add_option(
        options,
        c("-c", "--correction"),
        type = "character",
        dest = "correction",
        default = args$correction,
        help = paste("The correction method for computing",
            "the q-value [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-z", "--standardize"),
        type = "logical",
        dest = "standardize",
        default = args$standardize,
        help = paste("Apply z-score so continuous metadata are on",
            "the same scale [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-l", "--plot_heatmap"),
        type = "logical",
        dest = "plot_heatmap",
        default = args$plot_heatmap,
        help = paste("Generate a heatmap for the significant",
            "associations [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-i", "--heatmap_first_n"),
        type = "double",
        dest = "heatmap_first_n",
        default = args$heatmap_first_n,
        help = paste("In heatmap, plot top N features with significant",
            "associations [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-o", "--plot_scatter"),
        type = "logical",
        dest = "plot_scatter",
        default = args$plot_scatter,
        help = paste("Generate scatter plots for the significant",
            "associations [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-g", "--max_pngs"),
        type = "double",
        dest = "max_pngs",
        default = args$max_pngs,
        help = paste("The maximum number of scatterplots for signficant",
                     "associations to save as png files [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-O", "--save_scatter"),
        type = "logical",
        dest = "save_scatter",
        default = args$save_scatter,
        help = paste("Save all scatter plot ggplot objects",
                     "to an RData file [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-e", "--cores"),
        type = "double",
        dest = "cores",
        default = args$cores,
        help = paste("The number of R processes to ",
            "run in parallel [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-j", "--save_models"),
        type = "logical",
        dest = "save_models",
        default = args$save_models,
        help = paste("Return the full model outputs ",
                     "and save to an RData file [ Default: %default ]"
        )
    )
options <-
    optparse::add_option(
        options,
        c("-d", "--reference"),
        type = "character",
        dest = "reference",
        default = args$reference,
        help = paste("The factor to use as a reference for",
            "a variable with more than two levels",
            "provided as a string of 'variable,reference'",
            "semi-colon delimited for multiple variables [ Default: NA ]"
        )
    )

option_not_valid_error <- function(message, valid_options) {
<<<<<<< Updated upstream
    logging::logerror(paste(message, ": %s"), toString(valid_options))
    stop("Option not valid", call. = FALSE)
=======
  logging::logerror(paste(message, ": %s"), toString(valid_options))
  stop("Option not valid", call. = FALSE)
}

##################################################
# Fill in parameter list with missing parameters #
##################################################

maaslin_parse_param_list <- function(param_list) {
  if (sum(c("input_data", "input_metadata", "output") %in% names(param_list)) != 3) {
    stop("param_list must include input_data, input_metadata, and output")
  }
  default_params = list(min_abundance = args$min_abundance,
                        min_prevalence = args$min_prevalence,
                        min_variance = args$min_variance,
                        normalization = args$normalization,
                        transform = args$transform,
                        analysis_method = args$analysis_method,
                        max_significance = args$max_significance,
                        random_effects = args$random_effects,
                        fixed_effects = args$fixed_effects,
                        formula = args$formula,
                        correction = args$correction,
                        standardize = args$standardize,
                        cores = args$cores,
                        plot_heatmap = args$plot_heatmap,
                        heatmap_first_n = args$heatmap_first_n,
                        plot_scatter = args$plot_scatter,
                        max_pngs = args$max_pngs,
                        save_scatter = args$save_scatter,
                        save_models = args$save_models,
                        reference = args$reference)
  
  missing_params <- setdiff(names(default_params), names(param_list))
  
  for (param in missing_params) {
    param_list[[param]] <- default_params[[param]]
  }
  
  # Allow for lower case variables
  param_list[["normalization"]] <- toupper(param_list[["normalization"]])
  param_list[["transform"]] <- toupper(param_list[["transform"]])
  param_list[["analysis_method"]] <- toupper(param_list[["analysis_method"]])
  
  # Match variable ignoring case then set correctly as required for p.adjust
  param_list[["correction"]] <- correction_choices[match(toupper(param_list[["correction"]]), 
                                                         toupper(correction_choices))]
  
  return(param_list)
}

####################################
# Check valid options are selected #
####################################

maaslin_check_arguments <- function(param_list) {
  param_list <- maaslin_parse_param_list(param_list)
  
  # Check valid normalization option selected
  logging::loginfo("Verifying options selected are valid")
  if (!param_list[["normalization"]] %in% normalization_choices) {
    option_not_valid_error(
      paste(
        "Please select a normalization",
        "from the list of available options"),
      toString(normalization_choices)
    )
  }
  
  # check valid transform option selected
  if (!param_list[["transform"]] %in% transform_choices) {
    option_not_valid_error(
      "Please select a transform from the list of available options",
      toString(transform_choices)
    )
  }
  
  # check valid method option selected
  if (!param_list[["analysis_method"]] %in% analysis_method_choices_names) {
    option_not_valid_error(
      paste(
        "Please select an analysis method",
        "from the list of available options"),
      toString(analysis_method_choices_names)
    )
  }
  
  # check valid correction method selected
  if (!param_list[["correction"]] %in% correction_choices) {
    option_not_valid_error(
      paste("Please select a correction method",
            "from the list of available options"),
      toString(correction_choices)
    )
  }
  
  # check a valid choice combination is selected
  for (limited_transform in hash::keys(valid_choice_transform_norm)) {
    if (param_list[["transform"]] == limited_transform) {
      if (!param_list[["normalization"]] %in% 
          valid_choice_transform_norm[[limited_transform]]) {
        option_not_valid_error(
          paste0("This transform can only be used",
                 " with a subset of normalizations. ",
                 "Please select from the following list"
          ),
          toString(
            valid_choice_transform_norm[[limited_transform]])
        )
      }
    }
  }
  
  # check that plots are generated if to be saved
  if (!param_list[["plot_scatter"]] && param_list[["save_scatter"]]) {
    logging::logerror("Scatter plots cannot be saved if they are not plotted")
    stop("Option not valid", call. = FALSE)
  }
  
  return(param_list)
}

#####################################
# Create log file and log arguments #
#####################################

maaslin_log_arguments <- function(param_list) {
  param_list <- maaslin_parse_param_list(param_list)
  output <- param_list[["output"]]
  
  # create an output folder
  if (!file.exists(output)) {
    print("Creating output folder")
    dir.create(output)
  }
  
  # create log file (write info to stdout and debug level to log file)
  # set level to finest so all log levels are reviewed
  log_file <- file.path(output, "maaslin2.log")
  
  # remove log file if already exists (to avoid append)
  if (file.exists(log_file)) {
    print(paste("Warning: Deleting existing log file:", log_file))
    unlink(log_file)
  }
  
  logging::basicConfig(level = 'FINEST')
  logging::addHandler(logging::writeToFile, 
                      file = log_file, level = "DEBUG")
  logging::setLevel(20, logging::getHandler('basic.stdout'))
  
  logging::loginfo("Writing function arguments to log file")
  logging::logdebug("Function arguments")
  if (is.character(param_list[["input_data"]])) {
    logging::logdebug("Input data file: %s", param_list[["input_data"]])
  }
  if (is.character(param_list[["input_metadata"]])) {
    logging::logdebug("Input metadata file: %s", param_list[["input_metadata"]])
  }
  logging::logdebug("Output folder: %s", param_list[["output"]])
  logging::logdebug("Min Abundance: %f", param_list[["min_abundance"]])
  logging::logdebug("Min Prevalence: %f", param_list[["min_prevalence"]])
  logging::logdebug("Normalization: %s", param_list[["normalization"]])
  logging::logdebug("Transform: %s", param_list[["transform"]])
  logging::logdebug("Analysis method: %s", param_list[["analysis_method"]])
  logging::logdebug("Max significance: %f", param_list[["max_significance"]])
  logging::logdebug("Random effects: %s", param_list[["random_effects"]])
  logging::logdebug("Fixed effects: %s", param_list[["fixed_effects"]])
  logging::logdebug("Formula: %s", param_list[["formula"]])
  logging::logdebug("Correction method: %s", param_list[["correction"]])
  logging::logdebug("Standardize: %s", param_list[["standardize"]])
  logging::logdebug("Cores: %d", param_list[["cores"]])
  
  maaslin_check_arguments(param_list)
  
  return(param_list)
}

#################################
# Read in the data and metadata #
#################################

maaslin_read_data <- function(param_list) {
  param_list <- maaslin_parse_param_list(param_list)
  input_data <- param_list[["input_data"]]
  input_metadata <- param_list[["input_metadata"]]
  
  # if a character string then this is a file name, else it 
  # is a data frame
  if (is.character(input_data) && file.exists(input_data)) {
    data <-
      data.frame(data.table::fread(
        input_data, header = TRUE, sep = "\t"),
        row.names = 1)
    if (nrow(data) == 1) {
      # read again to get row name
      data <- read.table(input_data, header = TRUE, row.names = 1)
    }
  } else if (is.data.frame(input_data)) {
    if (!tibble::has_rownames(input_data)) {
      stop("If supplying input_data as a data frame, it must have appropriate rownames!")
    }
    data <- as.data.frame(input_data) # in case it's a tibble or something
  } else if (is.matrix(input_data)) {
    logging::logwarn("Input is a matrix, passing through as.data.frame() .")
    data <- as.data.frame(input_data)
  } else {
    stop("input_data is neither a file nor a data frame!")
  }
  
  if (is.character(input_metadata) && file.exists(input_metadata)) {
    metadata <-
      data.frame(data.table::fread(
        input_metadata, header = TRUE, sep = "\t"),
        row.names = 1)
    if (nrow(metadata) == 1) {
      metadata <- read.table(input_metadata,
                             header = TRUE,
                             row.names = 1)
    }
  } else if (is.data.frame(input_metadata)) {
    if (!tibble::has_rownames(input_metadata)) {
      stop("If supplying input_metadata as a data frame, it must have appropriate rownames!")
    }
    metadata <- as.data.frame(input_metadata) # in case it's a tibble or something
  } else {
    stop("input_metadata is neither a file nor a data frame!")
  }
  return(list("param_list" = param_list, 
              "data" = data, 
              "metadata" = metadata))
}

###############################################################
# Determine orientation of data in input and reorder to match #
###############################################################

maaslin_reorder_data <- function(params_and_data) {
  param_list <- maaslin_parse_param_list(params_and_data[["param_list"]])
  data <- params_and_data[["data"]]
  metadata <- params_and_data[["metadata"]]
  
  logging::loginfo("Determining format of input files")
  samples_row_row <- intersect(rownames(data), rownames(metadata))
  if (length(samples_row_row) > 0) {
    # this is the expected formatting so do not modify data frames
    logging::loginfo(
      paste(
        "Input format is data samples",
        "as rows and metadata samples as rows"))
  } else {
    samples_column_row <- intersect(colnames(data), rownames(metadata))
    
    if (length(samples_column_row) == 0) {
      # modify possibly included special chars in sample names in metadata
      rownames(metadata) <- make.names(rownames(metadata))
      
      samples_column_row <- intersect(colnames(data), rownames(metadata))
    }
    
    if (length(samples_column_row) > 0) {
      logging::loginfo(
        paste(
          "Input format is data samples",
          "as columns and metadata samples as rows"))
      # transpose data frame so samples are rows
      data <- as.data.frame(t(data))
      logging::logdebug(
        "Transformed data so samples are rows")
    } else {
      samples_column_column <- 
        intersect(colnames(data), colnames(metadata))
      if (length(samples_column_column) > 0) {
        logging::loginfo(
          paste(
            "Input format is data samples",
            "as columns and metadata samples as columns"))
        data <- as.data.frame(t(data))
        metadata <- type.convert(as.data.frame(t(metadata)))
        logging::logdebug(
          "Transformed data and metadata so samples are rows")
      } else {
        samples_row_column <- 
          intersect(rownames(data), colnames(metadata))
        
        if (length(samples_row_column) == 0) {
          # modify possibly included special chars in sample names in data
          rownames(data) <- make.names(rownames(data))
          
          samples_row_column <- intersect(rownames(data), colnames(metadata))
        }
        
        if (length(samples_row_column) > 0) {
          logging::loginfo(
            paste(
              "Input format is data samples",
              "as rows and metadata samples as columns"))
          metadata <- type.convert(as.data.frame(t(metadata)))
          logging::logdebug(
            "Transformed metadata so samples are rows")
        } else {
          logging::logerror(
            paste("Unable to find samples in data and",
                  "metadata files.",
                  "Rows/columns do not match."))
          logging::logdebug(
            "Data rows: %s", 
            paste(rownames(data), collapse = ","))
          logging::logdebug(
            "Data columns: %s", 
            paste(colnames(data), collapse = ","))
          logging::logdebug(
            "Metadata rows: %s", 
            paste(rownames(metadata), collapse = ","))
          logging::logdebug(
            "Metadata columns: %s",
            paste(colnames(data), collapse = ","))
          stop()
        }
      }
    }
  }
  
  # replace unexpected characters in feature names
  colnames(data) <- make.names(colnames(data))
  
  # check for samples without metadata
  extra_feature_samples <-
    setdiff(rownames(data), rownames(metadata))
  if (length(extra_feature_samples) > 0)
    logging::logdebug(
      paste("The following samples were found",
            "to have features but no metadata.",
            "They will be removed. %s"),
      paste(extra_feature_samples, collapse = ",")
    )
  
  # check for metadata samples without features
  extra_metadata_samples <-
    setdiff(rownames(metadata), rownames(data))
  if (length(extra_metadata_samples) > 0)
    logging::logdebug(
      paste("The following samples were found",
            "to have metadata but no features.",
            "They will be removed. %s"),
      paste(extra_metadata_samples, collapse = ",")
    )
  
  # get a set of the samples with both metadata and features
  intersect_samples <- intersect(rownames(data), rownames(metadata))
  logging::logdebug(
    "A total of %s samples were found in both the data and metadata",
    length(intersect_samples)
  )
  
  # now order both data and metadata with the same sample ordering
  logging::logdebug(
    "Reordering data/metadata to use same sample ordering")
  data <- data[intersect_samples, , drop = FALSE]
  metadata <- metadata[intersect_samples, , drop = FALSE]
  
  return(list("param_list" = param_list, 
              "data" = data, 
              "metadata" = metadata))
}

###########################################
# Compute the formula based on user input #
###########################################

maaslin_compute_formula <- function(params_and_data) {
  param_list <- maaslin_parse_param_list(params_and_data[["param_list"]])
  data <- params_and_data[["data"]]
  metadata <- params_and_data[["metadata"]]
  
  fixed_effects <- param_list[["fixed_effects"]]
  random_effects <- param_list[["random_effects"]]
  
  if (!is.null(param_list[["formula"]])) {
    if (!is.null(param_list[["fixed_effects"]]) | 
         !is.null(param_list[["random_effects"]])) {
      logging::logwarn(
        paste("fixed_effects and random_effects provided in addition to formula,", 
              "using fixed_effects and random_effects"))
    } else {
      logging::logwarn(
        paste("maaslin_compute_formula called even though a formula is provided",
              "without fixed_effects or random_effects,",
              "creating new formula based on metadata"))
    }
  }
  
  random_effects_formula <- NULL
  # use all metadata if no fixed effects are provided
  if (is.null(fixed_effects)) {
    fixed_effects <- colnames(metadata)
  } else {
    fixed_effects <- unlist(strsplit(fixed_effects, ",", fixed = TRUE))
    # remove any fixed effects not found in metadata names
    to_remove <- setdiff(fixed_effects, colnames(metadata))
    if (length(to_remove) > 0)
      logging::logwarn(
        paste("Feature name not found in metadata",
              "so not applied to formula as fixed effect: %s"),
        paste(to_remove, collapse = " , ")
      )
    
    fixed_effects <- setdiff(fixed_effects, to_remove)
    if (length(fixed_effects) == 0) {
      logging::logerror("No fixed effects included in formula.")
      stop()
    }
  }
  
  if (!is.null(random_effects)) {
    random_effects <-
      unlist(strsplit(random_effects, ",", fixed = TRUE))
    # subtract random effects from fixed effects
    common_variables <- intersect(fixed_effects, random_effects)
    if (length(common_variables) > 0) {
      logging::logerror(
        paste("Feature name included as fixed and random effect,",
              "check that this is intended: %s"),
        paste(common_variables, collapse = " , ")
      )
      stop()
    }
    # remove any random effects not found in metadata
    to_remove <- setdiff(random_effects, colnames(metadata))
    if (length(to_remove) > 0) {
      logging::logerror(
        paste("Feature name not found in metadata",
              "so not applied to formula as random effect: %s"),
        paste(to_remove, collapse = " , ")
      )
      stop()
    }
      
    random_effects <- setdiff(random_effects, to_remove)
    
    # create formula
    if (length(random_effects) > 0) {
      random_effects_formula_text <-
        paste(
          "expr ~ (1 | ",
          paste(
            random_effects,
            ")",
            sep = '',
            collapse = " + (1 | "
          ),
          sep = '')
      logging::loginfo("Formula for random effects: %s",
                       random_effects_formula_text)
      random_effects_formula <-
        tryCatch(
          as.formula(random_effects_formula_text),
          error = function(e)
            stop(
              paste(
                "Invalid formula for random effects: ",
                random_effects_formula_text
              )
            )
        )
    }
  }
  
  # reduce metadata to only include fixed/random effects in formula
  effects_names <- union(fixed_effects, random_effects)
  metadata <- metadata[, effects_names, drop = FALSE]
  
  # create the fixed effects formula text
  formula_text <-
    paste("expr ~ ", paste(fixed_effects, collapse = " + "))
  logging::loginfo("Formula for fixed effects: %s", formula_text)
  formula <-
    tryCatch(
      as.formula(formula_text),
      error = function(e)
        stop(
          paste(
            "Invalid formula.",
            "Please provide a different formula: ",
            formula_text
          )
        )
    )
  
  if (!(is.null(random_effects_formula))) {
    formula <-
      paste(
        '. ~', 
        paste(all.vars(formula)[-1], collapse = ' + '), 
        '.', 
        sep = ' + ')
    formula <- update(random_effects_formula, formula)
  }
  
  
  return(list("param_list" = param_list, 
              "data" = data, 
              "metadata" = metadata, 
              "formula" = list("formula" = formula, 
                               "random_effects_formula" = random_effects_formula)))
}

##############################
# Check a user input formula #
##############################

maaslin_check_formula <- function(params_and_data) {
  param_list <- maaslin_parse_param_list(params_and_data[["param_list"]])
  data <- params_and_data[["data"]]
  metadata <- params_and_data[["metadata"]]
  
  input_formula <- param_list[["formula"]]
  
  random_effects_formula <- NULL
  # use all metadata if no fixed effects are provided
  if (is.null(input_formula)) {
    logging::logwarn(
      paste("No user formula provided,",
            "building one from fixed_effects and random_effects"))
    return(maaslin_compute_formula(params_and_data))
  }
  
  if (!is.null(param_list[["fixed_effects"]]) | !is.null(param_list[["random_effects"]])) {
    logging::logwarn(
      paste("fixed_effects and random_effects provided in addition to formula,", 
            "using only formula"))
  }
  
  # Remove anything before the tilde if necessary
  input_formula <- sub(".*~\\s*", "", input_formula)
  input_formula <- paste0("expr ~ ", input_formula)
  
  formula <-
    tryCatch(
      as.formula(input_formula),
      error = function(e)
        stop(
          paste(
            "Invalid formula: ",
            input_formula
          )
        )
    )
  
  formula_terms <- all.vars(formula)
  formula_terms <- formula_terms[formula_terms != "expr"]
  
  to_remove <- setdiff(formula_terms, colnames(metadata))
  if (length(to_remove) > 0) {
    logging::logerror(
      paste("Feature name not found in metadata: %s"),
      paste(to_remove, collapse = " , ")
    )
    stop()
  }
  
  term_labels <- attr(terms(formula), "term.labels")
  
  if (sum(!grepl("\\|", term_labels)) == 0) {
    logging::logerror("No fixed effects included in formula.")
    stop()
  }
  
  # create formula
  if (sum(grepl("\\|", term_labels)) > 0) {
    random_effects_formula_text <- deparse(formula)
    logging::loginfo("Formula for random effects: %s",
                     input_formula)
    random_effects_formula <- formula
  } else {
    # create the fixed effects formula text
    formula_text <- deparse(formula)
    logging::loginfo("Formula for fixed effects: %s", formula_text)
    random_effects_formula <- NULL
  }
  
  # reduce metadata to only include fixed/random effects in formula
  metadata <- metadata[, formula_terms, drop = FALSE]
  
  return(list("param_list" = param_list, 
              "data" = data, 
              "metadata" = metadata, 
              "formula" = list("formula" = formula, 
                               "random_effects_formula" = random_effects_formula)))
}

##############################################################################################
# Filter data based on min abundance, min prevalence, and min variance; standardize metadata #
##############################################################################################

maaslin_filter_and_standardize <- function(params_and_data_and_formula) {
  param_list <- maaslin_parse_param_list(params_and_data_and_formula[["param_list"]])
  data <- params_and_data_and_formula[["data"]]
  metadata <- params_and_data_and_formula[["metadata"]]
  
  #########################################################
  # Filter data based on min abundance and min prevalence #
  #########################################################
  
  reference <- param_list[["reference"]]
  if (is.null(reference)) {
    reference <- ","
  }
  split_reference <- unlist(strsplit(reference, "[,;]"))
  
  fixed_effects <- param_list[["fixed_effects"]]
  # for each fixed effect, check that a reference level has been set if necessary: number of levels > 2 and metadata isn't already an ordered factor
  for (i in fixed_effects) {
    # don't check for or require reference levels for numeric metadata
    if (is.numeric(metadata[,i])) {
      next
    }
    # respect ordering if a factor is explicitly passed in with no reference set
    if (is.factor(metadata[,i]) && !(i %in% split_reference)) {
      logging::loginfo(paste("Factor detected for categorial metadata '", 
                             i, "'. Provide a reference argument or manually set factor ordering to change reference level.", sep=""))
      next
    }
    
    # set metadata as a factor (ordered alphabetically)
    metadata[,i] <- as.factor(metadata[,i])
    mlevels <- levels(metadata[,i])
    
    # get reference level for variable being considered, returns NA if not found
    ref <- split_reference[match(i, split_reference)+1]
    
    # if metadata has 2 levels, allow but don't require setting reference level, otherwise require it
    if ((length(mlevels) == 2)) {
      if(!is.na(ref)) {
        metadata[,i] = relevel(metadata[,i], ref = ref)
      }
    } else if (length(mlevels) > 2) {
      if (!is.na(ref)) {
        metadata[,i] = relevel(metadata[,i], ref = ref)
      } else {
        stop(paste("Please provide the reference for the variable '",
                   i, "' which includes more than 2 levels: ",
                   paste(as.character(mlevels), collapse=", "), ".", sep=""))   
      } 
    } else {
      stop("Provided categorical metadata has fewer than 2 unique, non-NA values.")
    }
  }
  
  unfiltered_data <- data
  unfiltered_metadata <- metadata
  
  # require at least total samples * min prevalence values 
  # for each feature to be greater than min abundance
  logging::loginfo(
    "Filter data based on min abundance and min prevalence")
  total_samples <- nrow(unfiltered_data)
  logging::loginfo("Total samples in data: %d", total_samples)
  min_samples <- total_samples * param_list[["min_prevalence"]]
  logging::loginfo(
    paste("Min samples required with min abundance",
          "for a feature not to be filtered: %f"),
    min_samples
  )
  
  # Filter by abundance using zero as value for NAs
  data_zeros <- unfiltered_data
  data_zeros[is.na(data_zeros)] <- 0
  filtered_data <-
    unfiltered_data[, 
                    colSums(data_zeros > param_list[["min_abundance"]]) > min_samples,
                    drop = FALSE]
  total_filtered_features <-
    ncol(unfiltered_data) - ncol(filtered_data)
  logging::loginfo("Total filtered features: %d", total_filtered_features)
  filtered_feature_names <-
    setdiff(names(unfiltered_data), names(filtered_data))
  logging::loginfo("Filtered feature names from abundance and prevalence filtering: %s",
                   toString(filtered_feature_names))
  
  #################################
  # Filter data based on variance #
  #################################
  
  sds <- apply(filtered_data, 2, sd)
  variance_filtered_data <- filtered_data[, which(sds > param_list[["min_variance"]]), drop = FALSE]
  variance_filtered_features <- ncol(filtered_data) - ncol(variance_filtered_data)
  logging::loginfo("Total filtered features with variance filtering: %d", variance_filtered_features)
  variance_filtered_feature_names <- setdiff(names(filtered_data), names(variance_filtered_data))
  logging::loginfo("Filtered feature names from variance filtering: %s",
                   toString(variance_filtered_feature_names))
  filtered_data <- variance_filtered_data
  
  output <- param_list[["output"]]
  features_folder <- file.path(output, "features")
  if (!file.exists(features_folder)) {
    print("Creating output feature tables folder")
    dir.create(features_folder, recursive = T)
  }
  
  filtered_file = file.path(features_folder, "filtered_data.tsv")
  logging::loginfo("Writing filtered data to file %s", filtered_file)
  write.table(
    data.frame("feature" = rownames(filtered_data), filtered_data), 
    file = filtered_file, 
    sep = "\t", 
    quote = FALSE, 
    row.names = FALSE
  )
  
  ################################
  # Standardize metadata, if set #
  ################################
  
  if (param_list[["standardize"]]) {
    logging::loginfo(
      "Applying z-score to standardize continuous metadata")
    metadata <- metadata %>% dplyr::mutate_if(is.numeric, scale)
  } else {
    logging::loginfo("Bypass z-score application to metadata")
  }
  
  return(list("param_list" = params_and_data_and_formula[["param_list"]],
              "data" = data, 
              "metadata" = metadata,
              "filtered_data" = filtered_data, 
              "unfiltered_metadata" = unfiltered_metadata, 
              "formula" = params_and_data_and_formula[["formula"]]))
}

#################
# Normalization #
#################

maaslin_normalize = function(params_and_data_and_formula) {
  param_list <- maaslin_parse_param_list(params_and_data_and_formula[["param_list"]])
  features <- params_and_data_and_formula[["filtered_data"]]
  normalization <- param_list[["normalization"]]
  
  logging::loginfo(
    "Running selected normalization method: %s", normalization)
  
  if (normalization == 'TSS') {features <- TSSnorm(features)}
  if (normalization == 'CLR') {features <- CLRnorm(features)}
  if (normalization == 'CSS') {features <- CSSnorm(features)}
  if (normalization == 'TMM') {features <- TMMnorm(features)}
  if (normalization == 'NONE') {features <- features}
  
  output <- param_list[["output"]]
  features_folder <- file.path(output, "features")
  if (!file.exists(features_folder)) {
    print("Creating output feature tables folder")
    dir.create(features_folder, recursive = T)
  }
  
  filtered_data_norm_file = file.path(features_folder, "filtered_data_norm.tsv")
  logging::loginfo("Writing filtered, normalized data to file %s", filtered_data_norm_file)
  write.table(
    data.frame("feature" = rownames(features), features), 
    file = filtered_data_norm_file, 
    sep = "\t", 
    quote = FALSE, 
    row.names = FALSE
  )
  
  return(list("param_list" = params_and_data_and_formula[["param_list"]], 
              "data" = params_and_data_and_formula[["data"]], 
              "metadata" = params_and_data_and_formula[["metadata"]],
              "unfiltered_metadata" = params_and_data_and_formula[["unfiltered_metadata"]], 
              "filtered_data" = params_and_data_and_formula[["filtered_data"]], 
              "filtered_data_norm" = features,
              "formula" = params_and_data_and_formula[["formula"]]))
}

##################
# Transformation #
##################

maaslin_transform = function(params_and_data_and_formula) {
  param_list <- maaslin_parse_param_list(params_and_data_and_formula[["param_list"]])
  features <- params_and_data_and_formula[["filtered_data_norm"]]
  transformation <- param_list[["transform"]]
  
  logging::loginfo("Running selected transform method: %s", transformation)
  
  if (transformation == 'LOG') {features <- apply(features, 2, LOG)}
  if (transformation == 'LOGIT') {features <- apply(features, 2, LOGIT)}
  if (transformation == 'AST') {features <- apply(features, 2, AST)}
  
  output <- param_list[["output"]]
  features_folder <- file.path(output, "features")
  if (!file.exists(features_folder)) {
    print("Creating output feature tables folder")
    dir.create(features_folder, recursive = T)
  }
  
  filtered_data_norm_transformed_file = file.path(features_folder, "filtered_data_norm_transformed.tsv")
  logging::loginfo("Writing filtered, normalized, transformed data to file %s", filtered_data_norm_transformed_file)
  write.table(
    data.frame("feature" = rownames(features), features), 
    file = filtered_data_norm_transformed_file, 
    sep = "\t", 
    quote = FALSE, 
    row.names = FALSE
  )
  
  return(list("param_list" = params_and_data_and_formula[["param_list"]], 
              "data" = params_and_data_and_formula[["data"]], 
              "metadata" = params_and_data_and_formula[["metadata"]],
              "unfiltered_metadata" = params_and_data_and_formula[["unfiltered_metadata"]], 
              "filtered_data" = params_and_data_and_formula[["filtered_data"]], 
              "filtered_data_norm_transformed" = features,
              "formula" = params_and_data_and_formula[["formula"]]))
}

##########
# Fit LM #
##########

maaslin_fit = function(params_and_data_and_formula) {
  param_list <- maaslin_parse_param_list(params_and_data_and_formula[["param_list"]])
  analysis_method <- param_list[["analysis_method"]]
  
  logging::loginfo(
    "Running selected analysis method: %s", analysis_method)
  
  fit_data <-
    fit.data(
      params_and_data_and_formula[["filtered_data_norm_transformed"]],
      params_and_data_and_formula[["metadata"]],
      analysis_method,
      formula = params_and_data_and_formula[["formula"]][["formula"]],
      random_effects_formula = params_and_data_and_formula[["formula"]][["random_effects_formula"]],
      correction = param_list[["correction"]],
      save_models = param_list[["save_models"]],
      cores = param_list[["cores"]]
    )
  
  #################################################################
  # Count the total values for each feature (untransformed space) #
  #################################################################
  
  logging::loginfo("Counting total values for each feature")
  
  fit_data$results$N <-
    apply(
      fit_data$results,
      1,
      FUN = function(x)
        length(params_and_data_and_formula[["filtered_data_norm_transformed"]][, x[1]])
    )
  fit_data$results$N.not.zero <-
    apply(
      fit_data$results,
      1,
      FUN = function(x)
        length(which(params_and_data_and_formula[["filtered_data"]][, x[1]] > 0))
    )
  
  return(list("param_list" = params_and_data_and_formula[["param_list"]], 
              "data" = params_and_data_and_formula[["data"]], 
              "metadata" = params_and_data_and_formula[["metadata"]],
              "unfiltered_metadata" = params_and_data_and_formula[["unfiltered_metadata"]], 
              "filtered_data" = params_and_data_and_formula[["filtered_data"]], 
              "filtered_data_norm_transformed" = params_and_data_and_formula[["filtered_data_norm_transformed"]],
              "formula" = params_and_data_and_formula[["formula"]],
              "fit_data" = fit_data
              ))
}

###########################
# Write tables of results #
###########################

maaslin_write_results <- function(params_data_formula_fit) {
  param_list <- maaslin_parse_param_list(params_data_formula_fit[["param_list"]])
  output <- param_list[["output"]]
  
  # create an output folder if it does not exist
  if (!file.exists(output)) {
    print("Creating output folder")
    dir.create(output)
  }
  
  write_fits(params_data_formula_fit)
  write_results(params_data_formula_fit)
  
  return(params_data_formula_fit)
}

#######################################################
# Create visualizations for results passing threshold #
#######################################################

maaslin_plot_results <- function(params_data_formula_fit) {
  param_list <- maaslin_parse_param_list(params_data_formula_fit[["param_list"]])
  output <- param_list[["output"]]
  
  # create an output folder and figures folder if it does not exist
  if (!file.exists(output)) {
    print("Creating output folder")
    dir.create(output)
  }
  
  if (param_list[["plot_heatmap"]] || param_list[["plot_scatter"]]) {
    figures_folder <- file.path(output, "figures")
    if (!file.exists(figures_folder)) {
      print("Creating output figures folder")
      dir.create(figures_folder)
    }
  }
  
  if (param_list[["plot_heatmap"]]) {
    heatmap_file <- file.path(output, "heatmap.pdf")
    logging::loginfo(
      "Writing heatmap of significant results to file: %s",
      heatmap_file)
    significant_results_file <-
      file.path(output, "significant_results.tsv")
    save_heatmap(significant_results_file, heatmap_file, figures_folder,
                 first_n = param_list[["heatmap_first_n"]])
  }
  
  if (param_list[["plot_scatter"]]) {
    logging::loginfo(
      paste("Writing association plots",
            "(one for each significant association)",
            "to output folder: %s"),
      output
    )
    significant_results_file <-
      file.path(output, "significant_results.tsv")
    saved_plots <- maaslin2_association_plots(
      params_data_formula_fit[["unfiltered_metadata"]],
      params_data_formula_fit[["filtered_data"]],
      significant_results_file,
      output,
      figures_folder,
      param_list[["max_pngs"]],
      param_list[["save_scatter"]])
    if (param_list[["save_scatter"]]) {
      scatter_file <- file.path(figures_folder, "scatter_plots.rds")
      # remove plots file if already exists
      if (file.exists(scatter_file)) {
        logging::logwarn(
          "Deleting existing scatter plot objects file: %s", scatter_file)
        unlink(scatter_file)
      }
      logging::loginfo("Writing scatter plot objects to file %s", scatter_file)
      saveRDS(saved_plots, file = scatter_file)   
    }
  }
  
  return(params_data_formula_fit)
>>>>>>> Stashed changes
}

#######################################################
# Main maaslin2 function (defaults same command line) #
#######################################################

<<<<<<< Updated upstream
Maaslin2 <-
    function(
        input_data,
        input_metadata,
        output,
        min_abundance = 0.0,
        min_prevalence = 0.1,
        min_variance = 0.0,
        normalization = "TSS",
        transform = "LOG",
        analysis_method = "LM",
        max_significance = 0.25,
        random_effects = NULL,
        fixed_effects = NULL,
        correction = "BH",
        standardize = TRUE,
        cores = 1,
        plot_heatmap = TRUE,
        heatmap_first_n = 50,
        plot_scatter = TRUE,
        max_pngs = 10,
        save_scatter = FALSE,
        save_models = FALSE,
        reference = NULL)
    {
        # Allow for lower case variables
        normalization <- toupper(normalization)
        transform <- toupper(transform)
        analysis_method <- toupper(analysis_method)
        
        # Match variable ignoring case then set correctly as required for p.adjust
        correction <- correction_choices[match(toupper(correction), toupper(correction_choices))]

        #################################################################
        # Read in the data and metadata, create output folder, init log #
        #################################################################
        # if a character string then this is a file name, else it 
        # is a data frame
        if (is.character(input_data) && file.exists(input_data)) {
            data <-
                data.frame(data.table::fread(
                    input_data, header = TRUE, sep = "\t"),
                        row.names = 1)
            if (nrow(data) == 1) {
                # read again to get row name
                data <- read.table(input_data, header = TRUE, row.names = 1)
            }
        } else if (is.data.frame(input_data)) {
            if (!tibble::has_rownames(input_data)) {
              stop("If supplying input_data as a data frame, it must have appropriate rownames!")
            }
            data <- as.data.frame(input_data) # in case it's a tibble or something
        } else if (is.matrix(input_data)) {
            logging::logwarn("Input is a matrix, passing through as.data.frame() .")
            data <- as.data.frame(input_data)
        } else {
            stop("input_data is neither a file nor a data frame!")
        }

        if (is.character(input_metadata) && file.exists(input_metadata)) {
            metadata <-
                data.frame(data.table::fread(
                    input_metadata, header = TRUE, sep = "\t"),
                        row.names = 1)
            if (nrow(metadata) == 1) {
                metadata <- read.table(input_metadata,
                    header = TRUE,
                    row.names = 1)
            }
        } else if (is.data.frame(input_metadata)) {
            if (!tibble::has_rownames(input_metadata)) {
              stop("If supplying input_metadata as a data frame, it must have appropriate rownames!")
            }
            metadata <- as.data.frame(input_metadata) # in case it's a tibble or something
        } else {
          stop("input_metadata is neither a file nor a data frame!")
        } 
        # create an output folder and figures folder if it does not exist
        if (!file.exists(output)) {
            print("Creating output folder")
            dir.create(output)
        }
        
        features_folder <- file.path(output, "features")
        if (!file.exists(features_folder)) {
            print("Creating output feature tables folder")
            dir.create(features_folder)
        }
        
        fits_folder <- file.path(output, "fits")
        if (!file.exists(fits_folder)) {
            print("Creating output fits folder")
            dir.create(fits_folder)
        }

        if (plot_heatmap || plot_scatter) {
            figures_folder <- file.path(output, "figures")
            if (!file.exists(figures_folder)) {
                print("Creating output figures folder")
                dir.create(figures_folder)
            }
        }
    
        
        # create log file (write info to stdout and debug level to log file)
        # set level to finest so all log levels are reviewed
        log_file <- file.path(output, "maaslin2.log")
        # remove log file if already exists (to avoid append)
        if (file.exists(log_file)) {
            print(paste("Warning: Deleting existing log file:", log_file))
            unlink(log_file)
        }
        logging::basicConfig(level = 'FINEST')
        logging::addHandler(logging::writeToFile, 
            file = log_file, level = "DEBUG")
        logging::setLevel(20, logging::getHandler('basic.stdout'))
        
        #####################
        # Log the arguments #
        #####################
        
        logging::loginfo("Writing function arguments to log file")
        logging::logdebug("Function arguments")
        if (is.character(input_data)) {
            logging::logdebug("Input data file: %s", input_data)
        }
        if (is.character(input_metadata)) {
            logging::logdebug("Input metadata file: %s", input_metadata)
        }
        logging::logdebug("Output folder: %s", output)
        logging::logdebug("Min Abundance: %f", min_abundance)
        logging::logdebug("Min Prevalence: %f", min_prevalence)
        logging::logdebug("Normalization: %s", normalization)
        logging::logdebug("Transform: %s", transform)
        logging::logdebug("Analysis method: %s", analysis_method)
        logging::logdebug("Max significance: %f", max_significance)
        logging::logdebug("Random effects: %s", random_effects)
        logging::logdebug("Fixed effects: %s", fixed_effects)
        logging::logdebug("Correction method: %s", correction)
        logging::logdebug("Standardize: %s", standardize)
        logging::logdebug("Cores: %d", cores)
        
        ####################################
        # Check valid options are selected #
        ####################################
        
        # Check valid normalization option selected
        logging::loginfo("Verifying options selected are valid")
        if (!normalization %in% normalization_choices) {
            option_not_valid_error(
                paste(
                    "Please select a normalization",
                    "from the list of available options"),
                toString(normalization_choices)
            )
        }
        
        # check valid transform option selected
        if (!transform %in% transform_choices) {
            option_not_valid_error(
                "Please select a transform from the list of available options",
                toString(transform_choices)
            )
        }
        
        # check valid method option selected
        if (!analysis_method %in% analysis_method_choices_names) {
            option_not_valid_error(
                paste(
                    "Please select an analysis method",
                    "from the list of available options"),
                toString(analysis_method_choices_names)
            )
        }
        
        # check valid correction method selected
        if (!correction %in% correction_choices) {
            option_not_valid_error(
                paste("Please select a correction method",
                    "from the list of available options"),
                toString(correction_choices)
            )
        }
        
        # check a valid choice combination is selected
        for (limited_method in hash::keys(
            valid_choice_method_norm)) {
            if (analysis_method == limited_method) {
                if (!normalization %in% 
                    valid_choice_method_norm[[limited_method]]) {
                    option_not_valid_error(
                        paste0("This method can only be used ",
                            "with a subset of normalizations. ",
                            "Please select from the following list"
                        ),
                        toString(valid_choice_method_norm[[limited_method]])
                    )
                }
            }
        }
        for (limited_transform in hash::keys(valid_choice_transform_norm)) {
            if (transform == limited_transform) {
                if (!normalization %in% 
                    valid_choice_transform_norm[[limited_transform]]) {
                    option_not_valid_error(
                        paste0("This transform can only be used",
                            " with a subset of normalizations. ",
                            "Please select from the following list"
                        ),
                        toString(
                            valid_choice_transform_norm[[limited_transform]])
                    )
                }
            }
        }
        
        # check that the transform can be applied to the method selected
        if (transform != "NONE")
        {
            if (!analysis_method %in% valid_choice_method_transform) {
                option_not_valid_error(
                    paste0("The transform selected can only be used",
                        " with some methods. ",
                        "Please select from the following list"
                    ),
                    toString(valid_choice_method_transform)
                )
            }
        }
        # check that plots are generated if to be saved
        if (!plot_scatter && save_scatter) {
            logging::logerror("Scatter plots cannot be saved if they are not plotted")
            stop("Option not valid", call. = FALSE)
        }
        
        ###############################################################
        # Determine orientation of data in input and reorder to match #
        ###############################################################
        
        logging::loginfo("Determining format of input files")
        samples_row_row <- intersect(rownames(data), rownames(metadata))
        if (length(samples_row_row) > 0) {
            # this is the expected formatting so do not modify data frames
            logging::loginfo(
                paste(
                    "Input format is data samples",
                    "as rows and metadata samples as rows"))
        } else {
            samples_column_row <- intersect(colnames(data), rownames(metadata))

            if (length(samples_column_row) == 0) {
                # modify possibly included special chars in sample names in metadata
                rownames(metadata) <- make.names(rownames(metadata))
            
                samples_column_row <- intersect(colnames(data), rownames(metadata))
            }

            if (length(samples_column_row) > 0) {
                logging::loginfo(
                    paste(
                        "Input format is data samples",
                        "as columns and metadata samples as rows"))
                # transpose data frame so samples are rows
                data <- as.data.frame(t(data))
                logging::logdebug(
                    "Transformed data so samples are rows")
            } else {
                samples_column_column <- 
                    intersect(colnames(data), colnames(metadata))
                if (length(samples_column_column) > 0) {
                    logging::loginfo(
                        paste(
                            "Input format is data samples",
                            "as columns and metadata samples as columns"))
                    data <- as.data.frame(t(data))
                    metadata <- type.convert(as.data.frame(t(metadata)))
                    logging::logdebug(
                        "Transformed data and metadata so samples are rows")
                } else {
                    samples_row_column <- 
                        intersect(rownames(data), colnames(metadata))

                    if (length(samples_row_column) == 0) {
                        # modify possibly included special chars in sample names in data
                        rownames(data) <- make.names(rownames(data))
            
                        samples_row_column <- intersect(rownames(data), colnames(metadata))
                    }

                    if (length(samples_row_column) > 0) {
                        logging::loginfo(
                            paste(
                                "Input format is data samples",
                                "as rows and metadata samples as columns"))
                        metadata <- type.convert(as.data.frame(t(metadata)))
                        logging::logdebug(
                            "Transformed metadata so samples are rows")
                    } else {
                        logging::logerror(
                            paste("Unable to find samples in data and",
                                "metadata files.",
                                "Rows/columns do not match."))
                        logging::logdebug(
                            "Data rows: %s", 
                            paste(rownames(data), collapse = ","))
                        logging::logdebug(
                            "Data columns: %s", 
                            paste(colnames(data), collapse = ","))
                        logging::logdebug(
                            "Metadata rows: %s", 
                            paste(rownames(metadata), collapse = ","))
                        logging::logdebug(
                            "Metadata columns: %s",
                            paste(colnames(data), collapse = ","))
                        stop()
                    }
                }
            }
        }
       
        # replace unexpected characters in feature names
        colnames(data) <- make.names(colnames(data))
 
        # check for samples without metadata
        extra_feature_samples <-
            setdiff(rownames(data), rownames(metadata))
        if (length(extra_feature_samples) > 0)
            logging::logdebug(
                paste("The following samples were found",
                    "to have features but no metadata.",
                    "They will be removed. %s"),
                paste(extra_feature_samples, collapse = ",")
            )
        
        # check for metadata samples without features
        extra_metadata_samples <-
            setdiff(rownames(metadata), rownames(data))
        if (length(extra_metadata_samples) > 0)
            logging::logdebug(
                paste("The following samples were found",
                    "to have metadata but no features.",
                    "They will be removed. %s"),
                paste(extra_metadata_samples, collapse = ",")
            )
        
        # get a set of the samples with both metadata and features
        intersect_samples <- intersect(rownames(data), rownames(metadata))
        logging::logdebug(
            "A total of %s samples were found in both the data and metadata",
            length(intersect_samples)
        )
        
        # now order both data and metadata with the same sample ordering
        logging::logdebug(
            "Reordering data/metadata to use same sample ordering")
        data <- data[intersect_samples, , drop = FALSE]
        metadata <- metadata[intersect_samples, , drop = FALSE]
        
        ###########################################
        # Compute the formula based on user input #
        ###########################################
        
        random_effects_formula <- NULL
        # use all metadata if no fixed effects are provided
        if (is.null(fixed_effects)) {
            fixed_effects <- colnames(metadata)
        } else {
            fixed_effects <- unlist(strsplit(fixed_effects, ",", fixed = TRUE))
            # remove any fixed effects not found in metadata names
            to_remove <- setdiff(fixed_effects, colnames(metadata))
            if (length(to_remove) > 0)
                logging::logwarn(
                    paste("Feature name not found in metadata",
                        "so not applied to formula as fixed effect: %s"),
                    paste(to_remove, collapse = " , ")
                )
            fixed_effects <- setdiff(fixed_effects, to_remove)
            if (length(fixed_effects) == 0) {
                logging::logerror("No fixed effects included in formula.")
                stop()
            }
        }
        
        if (!is.null(random_effects)) {
            random_effects <-
                unlist(strsplit(random_effects, ",", fixed = TRUE))
            # subtract random effects from fixed effects
            common_variables <- intersect(fixed_effects, random_effects)
            if (length(common_variables) > 0) {
                logging::logwarn(
                    paste("Feature name included as fixed and random effect,",
                          "check that this is intended: %s"),
                    paste(common_variables, collapse = " , ")
                )
            }
            # remove any random effects not found in metadata
            to_remove <- setdiff(random_effects, colnames(metadata))
            if (length(to_remove) > 0)
                logging::logwarn(
                    paste("Feature name not found in metadata",
                        "so not applied to formula as random effect: %s"),
                    paste(to_remove, collapse = " , ")
                )
            random_effects <- setdiff(random_effects, to_remove)
            
            # create formula
            if (length(random_effects) > 0) {
                random_effects_formula_text <-
                    paste(
                        "expr ~ (1 | ",
                        paste(
                            random_effects,
                            ")",
                            sep = '',
                            collapse = " + (1 | "
                        ),
                        sep = '')
                logging::loginfo("Formula for random effects: %s",
                    random_effects_formula_text)
                random_effects_formula <-
                    tryCatch(
                        as.formula(random_effects_formula_text),
                        error = function(e)
                            stop(
                                paste(
                                    "Invalid formula for random effects: ",
                                    random_effects_formula_text
                                )
                            )
                    )
            }
        }
        
        # reduce metadata to only include fixed/random effects in formula
        effects_names <- union(fixed_effects, random_effects)
        metadata <- metadata[, effects_names, drop = FALSE]
        
        # create the fixed effects formula text
        formula_text <-
            paste("expr ~ ", paste(fixed_effects, collapse = " + "))
        logging::loginfo("Formula for fixed effects: %s", formula_text)
        formula <-
            tryCatch(
                as.formula(formula_text),
                error = function(e)
                    stop(
                        paste(
                            "Invalid formula.",
                            "Please provide a different formula: ",
                            formula_text
                        )
                    )
            )
        
        #########################################################
        # Filter data based on min abundance and min prevalence #
        #########################################################

        if (is.null(reference)) {
            reference <- ","
        }
        split_reference <- unlist(strsplit(reference, "[,;]"))
        
        # for each fixed effect, check that a reference level has been set if necessary: number of levels > 2 and metadata isn't already an ordered factor
        for (i in fixed_effects) {
            # don't check for or require reference levels for numeric metadata
            if (is.numeric(metadata[,i])) {
                next
            }
            # respect ordering if a factor is explicitly passed in with no reference set
            if (is.factor(metadata[,i]) && !(i %in% split_reference)) {
                logging::loginfo(paste("Factor detected for categorial metadata '", 
                                       i, "'. Provide a reference argument or manually set factor ordering to change reference level.", sep=""))
                next
            }
            
            # set metadata as a factor (ordered alphabetically)
            metadata[,i] <- as.factor(metadata[,i])
            mlevels <- levels(metadata[,i])
            
            # get reference level for variable being considered, returns NA if not found
            ref <- split_reference[match(i, split_reference)+1]
            
            # if metadata has 2 levels, allow but don't require setting reference level, otherwise require it
            if ((length(mlevels) == 2)) {
                if(!is.na(ref)) {
                    metadata[,i] = relevel(metadata[,i], ref = ref)
                }
            } else if (length(mlevels) > 2) {
                if (!is.na(ref)) {
                    metadata[,i] = relevel(metadata[,i], ref = ref)
                } else {
                    stop(paste("Please provide the reference for the variable '",
                               i, "' which includes more than 2 levels: ",
                               paste(as.character(mlevels), collapse=", "), ".", sep=""))   
                } 
            } else {
                stop("Provided categorical metadata has fewer than 2 unique, non-NA values.")
            }
        }
 
        unfiltered_data <- data
        unfiltered_metadata <- metadata
        
        # require at least total samples * min prevalence values 
        # for each feature to be greater than min abundance
        logging::loginfo(
            "Filter data based on min abundance and min prevalence")
        total_samples <- nrow(unfiltered_data)
        logging::loginfo("Total samples in data: %d", total_samples)
        min_samples <- total_samples * min_prevalence
        logging::loginfo(
            paste("Min samples required with min abundance",
                "for a feature not to be filtered: %f"),
            min_samples
        )
        
        # Filter by abundance using zero as value for NAs
        data_zeros <- unfiltered_data
        data_zeros[is.na(data_zeros)] <- 0
        filtered_data <-
            unfiltered_data[, 
                colSums(data_zeros > min_abundance) > min_samples,
                drop = FALSE]
        total_filtered_features <-
            ncol(unfiltered_data) - ncol(filtered_data)
        logging::loginfo("Total filtered features: %d", total_filtered_features)
        filtered_feature_names <-
            setdiff(names(unfiltered_data), names(filtered_data))
        logging::loginfo("Filtered feature names from abundance and prevalence filtering: %s",
            toString(filtered_feature_names))
        
        #################################
        # Filter data based on variance #
        #################################
        
        sds <- apply(filtered_data, 2, sd)
        variance_filtered_data <- filtered_data[, which(sds > min_variance), drop = FALSE]
        variance_filtered_features <- ncol(filtered_data) - ncol(variance_filtered_data)
        logging::loginfo("Total filtered features with variance filtering: %d", variance_filtered_features)
        variance_filtered_feature_names <- setdiff(names(filtered_data), names(variance_filtered_data))
        logging::loginfo("Filtered feature names from variance filtering: %s",
                         toString(variance_filtered_feature_names))
        filtered_data <- variance_filtered_data
       
        ######################
        # Normalize features #
        ######################
        
        logging::loginfo(
            "Running selected normalization method: %s", normalization)
        filtered_data_norm <-
            normalizeFeatures(filtered_data, normalization = normalization)
        
        ################################
        # Standardize metadata, if set #
        ################################
        
        if (standardize) {
            logging::loginfo(
                "Applying z-score to standardize continuous metadata")
            metadata <- metadata %>% dplyr::mutate_if(is.numeric, scale)
        } else {
            logging::loginfo("Bypass z-score application to metadata")
        }
        
        ############################
        # Transform and run method #
        ############################
       
        # transform features
        logging::loginfo("Running selected transform method: %s", transform)
        filtered_data_norm_transformed <-
            transformFeatures(filtered_data_norm, transformation = transform)
        
        # apply the method to the data with the correction
        logging::loginfo(
            "Running selected analysis method: %s", analysis_method)
        fit_data <-
            fit.data(
                filtered_data_norm_transformed,
                metadata,
                analysis_method,
                formula = formula,
                random_effects_formula = random_effects_formula,
                correction = correction,
                save_models = save_models,
                cores = cores
            )
        
        #################################################################
        # Count the total values for each feature (untransformed space) #
        #################################################################
        
        logging::loginfo("Counting total values for each feature")
        fit_data$results$N <-
            apply(
                fit_data$results,
                1,
                FUN = function(x)
                    length(filtered_data_norm[, x[1]])
            )
        fit_data$results$N.not.zero <-
            apply(
                fit_data$results,
                1,
                FUN = function(x)
                    length(which(filtered_data[, x[1]] > 0))
            )
        
        ################################
        # Write out the raw model fits #
        ################################
        
        if (save_models) {
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
        
        ##########################################
        # Write processed feature tables to file #
        ##########################################
        
        filtered_file = file.path(features_folder, "filtered_data.tsv")
        logging::loginfo("Writing filtered data to file %s", filtered_file)
        write.table(
            data.frame("feature" = rownames(filtered_data), filtered_data), 
            file = filtered_file, 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE
            )
        
        filtered_data_norm_file = file.path(features_folder, "filtered_data_norm.tsv")
        logging::loginfo("Writing filtered, normalized data to file %s", filtered_data_norm_file)
        write.table(
            data.frame("feature" = rownames(filtered_data_norm), filtered_data_norm), 
            file = filtered_data_norm_file, 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE
        )
        
        filtered_data_norm_transformed_file = file.path(features_folder, "filtered_data_norm_transformed.tsv")
        logging::loginfo("Writing filtered, normalized, transformed data to file %s", filtered_data_norm_transformed_file)
        write.table(
            data.frame("feature" = rownames(filtered_data_norm_transformed), filtered_data_norm_transformed), 
            file = filtered_data_norm_transformed_file, 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE
        )
        
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
        
        if (!is.null(random_effects)) {
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
                "coef",
                "stderr",
                "N",
                "N.not.0",
                "pval",
                "qval"
            ),
            row.names = FALSE
        )
        
        #######################################################
        # Create visualizations for results passing threshold #
        #######################################################
        
        if (plot_heatmap) {
            heatmap_file <- file.path(output, "heatmap.pdf")
            logging::loginfo(
                "Writing heatmap of significant results to file: %s",
                heatmap_file)
            save_heatmap(significant_results_file, heatmap_file, figures_folder,
                first_n = heatmap_first_n)
        }
        
        if (plot_scatter) {
            logging::loginfo(
                paste("Writing association plots",
                    "(one for each significant association)",
                    "to output folder: %s"),
                output
            )
            saved_plots <- maaslin2_association_plots(
                unfiltered_metadata,
                filtered_data,
                significant_results_file,
                output,
                figures_folder,
                max_pngs,
                save_scatter)
            if (save_scatter) {
                scatter_file <- file.path(figures_folder, "scatter_plots.rds")
                # remove plots file if already exists
                if (file.exists(scatter_file)) {
                    logging::logwarn(
                        "Deleting existing scatter plot objects file: %s", scatter_file)
                    unlink(scatter_file)
                }
                logging::loginfo("Writing scatter plot objects to file %s", scatter_file)
                saveRDS(saved_plots, file = scatter_file)   
            }
        }
        
        return(fit_data)
    }
=======
Maaslin2 <- function(param_list = list()) {
  # Create log file, log arguments, and check arguments
  params_tmp <- maaslin_log_arguments(param_list) %>%
    maaslin_read_data() %>%
    maaslin_reorder_data() 
  
  if (is.null(params_tmp[['param_list']][['formula']])) {
    params_tmp_2 <- params_tmp %>%
      maaslin_compute_formula()
  } else {
    params_tmp_2 <- params_tmp %>%
      maaslin_check_formula()
  }
  
  params_data_formula_fit <- params_tmp_2 %>%
    maaslin_filter_and_standardize() %>%
    maaslin_normalize() %>%
    maaslin_transform() %>%
    maaslin_fit() %>%
    maaslin_write_results() %>%
    maaslin_plot_results()
  
  return(params_data_formula_fit[["fit_data"]])
}
>>>>>>> Stashed changes

###########################################################################
# If running on the command line, get arguments and call maaslin function #
###########################################################################

# this evaluates to true if script is being called directly as an executable
if (identical(environment(), globalenv()) &&
<<<<<<< Updated upstream
        !length(grep("^source\\(", sys.calls()))) {
    # get command line options and positional arguments
    parsed_arguments = optparse::parse_args(options, 
        positional_arguments = TRUE)
    current_args <- parsed_arguments[["options"]]
    positional_args <- parsed_arguments[["args"]]
    # check three positional arguments are provided
    if (length(positional_args) != 3) {
        optparse::print_help(options)
        stop(
            paste("Please provide the required",
                "positional arguments",
                "<data.tsv> <metadata.tsv> <output_folder>")
        )
    }
    
    # call maaslin with the command line options
    fit_data <-
        Maaslin2(
            positional_args[1],
            positional_args[2],
            positional_args[3],
            current_args$min_abundance,
            current_args$min_prevalence,
            current_args$min_variance,
            current_args$normalization,
            current_args$transform,
            current_args$analysis_method,
            current_args$max_significance,
            current_args$random_effects,
            current_args$fixed_effects,
            current_args$correction,
            current_args$standardize,
            current_args$cores,
            current_args$plot_heatmap,
            current_args$heatmap_first_n,
            current_args$plot_scatter,
            current_args$max_pngs,
            current_args$save_scatter,
            current_args$save_models,
            current_args$reference
        )
=======
    !length(grep("^source\\(", sys.calls()))) {
  # get command line options and positional arguments
  parsed_arguments = optparse::parse_args(options, 
                                          positional_arguments = TRUE)
  current_args <- parsed_arguments[["options"]]
  positional_args <- parsed_arguments[["args"]]
  # check three positional arguments are provided
  if (length(positional_args) != 3) {
    optparse::print_help(options)
    stop(
      paste("Please provide the required",
            "positional arguments",
            "<data.tsv> <metadata.tsv> <output_folder>")
    )
  }
  
  # call maaslin with the command line options
  fit_data <-
    Maaslin2(list(
      positional_args[1],
      positional_args[2],
      positional_args[3],
      current_args$min_abundance,
      current_args$min_prevalence,
      current_args$min_variance,
      current_args$normalization,
      current_args$transform,
      current_args$analysis_method,
      current_args$max_significance,
      current_args$random_effects,
      current_args$fixed_effects,
      current_args$formula,
      current_args$correction,
      current_args$standardize,
      current_args$cores,
      current_args$plot_heatmap,
      current_args$heatmap_first_n,
      current_args$plot_scatter,
      current_args$max_pngs,
      current_args$save_scatter,
      current_args$save_models,
      current_args$reference
    ))
>>>>>>> Stashed changes
}
