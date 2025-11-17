#' DNA Sequence Encoding Function
#'
#' Converts an input vector of DNA strings into a factor-encoded data frame (each nucleotide position as a column).
#'
#' @param dna_strings Character vector where each element is a DNA sequence of consistent length (containing only A/T/C/G).
#' @return Data frame with columns equal to the length of the DNA sequences (column names: `nt_pos1` to `nt_posN`).
#'         All columns are factor-type with levels: A, T, C, G.
dna_encoding <- function(dna_strings){
  nn <- nchar( dna_strings[1] )
  seq_m <- matrix( unlist( strsplit(dna_strings, "") ), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

#' Batch Prediction of m6A Sites
#'
#' Performs batch prediction of m6A site probability and status for multiple samples using a pre-trained machine learning model.
#'
#' @import randomForest
#' @importFrom stats predict
#' @param ml_fit Pre-trained machine learning model (e.g., random forest model).
#' @param feature_df Data frame that must contain the following columns:
#'                   gc_content (numeric), RNA_type (character; allowed values: mRNA/lincRNA/lncRNA/pseudogene),
#'                   RNA_region (character; allowed values: CDS/intron/3'UTR/5'UTR), exon_length (numeric),
#'                   distance_to_junction (numeric), evolutionary_conservation (numeric),
#'                   DNA_5mer (character; 5-nucleotide DNA sequence).
#' @param positive_threshold Numeric value (0-1) for positive status cutoff, default = 0.5.
#' @return Data frame with two additional columns added to the input feature data:
#'         predicted_m6A_prob (predicted probability), predicted_m6A_status (predicted status: Positive/Negative).
#' @export
#' @examples
#' # 1. Load pre-trained model and example feature data from the package
#' trained_model <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#' example_features <- read.csv(
#'   system.file("extdata", "m6A_input_example.csv", package = "m6APrediction")
#' )
#'
#' # 2. Run batch prediction with default positive threshold (0.5)
#' batch_pred_results <- prediction_multiple(ml_fit = trained_model, feature_df = example_features)
#'
#' # 3. View key results (first 5 rows, focusing on prediction columns)
#' head(
#'   batch_pred_results[, c(
#'     "gc_content", "RNA_type", "RNA_region",
#'     "predicted_m6A_prob", "predicted_m6A_status"
#'   )]
#' )
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction", "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df))) #Check errors if incorrect column names of input data.frame
  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  onehot_encode <- dna_encoding(feature_df$DNA_5mer)
  feature_df <- cbind(feature_df, onehot_encode)

  # 显式指定stats::predict，双重保障（避免NAMESPACE同步问题）
  pred_result <- stats::predict(ml_fit, newdata = feature_df, type = "prob")
  pred_prob <- pred_result[,"Positive"]
  pred_status <- ifelse(pred_prob > positive_threshold,"Positive", "Negative")

  feature_df$predicted_m6A_prob <- pred_prob
  feature_df$predicted_m6A_status <- factor(pred_status, levels = c("Negative", "Positive"))
  return(feature_df)
}

#' Single-Sample Prediction of m6A Sites
#'
#' Predicts m6A site probability and status for a single sample using a pre-trained machine learning model.
#'
#' @param ml_fit Pre-trained machine learning model (e.g., random forest model).
#' @param gc_content Numeric value, representing GC content.
#' @param RNA_type Character string, specifying RNA type (allowed values: mRNA/lincRNA/lncRNA/pseudogene).
#' @param RNA_region Character string, specifying RNA region (allowed values: CDS/intron/3'UTR/5'UTR).
#' @param exon_length Numeric value, representing exon length.
#' @param distance_to_junction Numeric value, representing distance to junction.
#' @param evolutionary_conservation Numeric value, representing evolutionary conservation.
#' @param DNA_5mer Character string, 5-nucleotide DNA sequence (containing only A/T/C/G).
#' @param positive_threshold Numeric value (0-1) for positive status cutoff, default = 0.5.
#' @return Named vector with two elements: predicted_m6A_prob (predicted probability), predicted_m6A_status (predicted status).
#' @export
#' @examples
#' # 1. Load pre-trained random forest model from the package
#' trained_model <- readRDS(
#'   system.file("extdata", "rf_fit.rds", package = "m6APrediction")
#' )
#'
#' # 2. Run single-sample prediction (input valid feature values manually)
#' single_sample_pred <- prediction_single(
#'   ml_fit = trained_model,
#'   gc_content = 0.62,
#'   RNA_type = "mRNA",          # Must be one of: mRNA/lincRNA/lncRNA/pseudogene
#'   RNA_region = "CDS",         # Must be one of: CDS/intron/3'UTR/5'UTR
#'   exon_length = 10.8,
#'   distance_to_junction = 8.9,
#'   evolutionary_conservation = 0.75,
#'   DNA_5mer = "ATCGG",         # Must be a 5-nucleotide string (A/T/C/G only)
#'   positive_threshold = 0.5
#' )
#'
#' # 3. View prediction results (named vector: probability + status)
#' print(single_sample_pred)
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length, distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  single_df <- data.frame(
    gc_content = gc_content,
    RNA_type = RNA_type,
    RNA_region = RNA_region,
    exon_length = exon_length,
    distance_to_junction = distance_to_junction,
    evolutionary_conservation = evolutionary_conservation,
    DNA_5mer = DNA_5mer,
    stringsAsFactors = FALSE
  )

  result_df <- prediction_multiple(ml_fit, single_df, positive_threshold)

  returned_vector <- c(
    predicted_m6A_prob = result_df$predicted_m6A_prob,
    predicted_m6A_status = as.character(result_df$predicted_m6A_status)
  )
  return(returned_vector)
}
