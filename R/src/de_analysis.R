#' Run differential expression analysis using limma
#' @param counts A matrix of counts with genes as rows and samples as columns.
#' @param sample_ids A vector of sample IDs.
#' @param gene_ids A vector of gene IDs.
#' @param design_matrix A matrix or data frame with sample IDs as rows and covariates as columns.
#' @param contrast_groups A list of contrast groups, where each group is a vector of sample IDs.
#' @param reference_groups A list of reference groups, where each group is a vector of sample IDs.
#' @return A list containing the results of the differential expression analysis.
#'
#' @import limma
#' @export
run_de_analysis <- function(counts, sample_ids, gene_ids, design_matrix, contrast_groups, reference_groups) {
  # Convert counts to ExpressionSet object
  eset <- Biobase::ExpressionSet(assayData = counts,
                                 phenoData = Biobase::AnnotatedDataFrame(design_matrix),
                                 featureData = Biobase::AnnotatedDataFrame(data.frame(gene_id = gene_ids)))

  # Normalize data (you may want to adjust this based on your specific needs)
  normalized_counts <- limma::voom(eset)

  # Create design matrix
  design <- model.matrix(~ 0 + design_matrix)
  colnames(design) <- colnames(design_matrix)

  # Fit the model
  fit <- limma::lmFit(normalized_counts, design)

  # Create contrasts
  contrasts <- limma::makeContrasts(contrasts = sapply(seq_along(contrast_groups), function(i) {
    paste0(contrast_groups[[i]], " - ", reference_groups[[i]])
  }), levels = design)

  # Fit contrasts
  fit2 <- limma::contrasts.fit(fit, contrasts)

  # Empirical Bayes statistics
  fit2 <- limma::eBayes(fit2)

  # Get results
  results <- lapply(seq_along(contrast_groups), function(i) {
    limma::topTable(fit2, coef = i, n = Inf, sort.by = "P")
  })

  names(results) <- paste0("Contrast_", seq_along(contrast_groups))

  return(results)
}
