#' Test for enrichment in a dataset with the hypergeometric test
#'
#' @param go_sample_counts, data.frame with the following columns:
#' "goterm_id", character - Either the identifier or description of the GO term of interest to test for enrichment in sample.
#' "sample_id", character - The sample identifier.
#' "go_of_interest_in_sample", integer - Number of hits of the GO term of interest in the sample (q parameter in phyper), 'hit in sample'.
#' "go_of_interest_in_pop" - How many times the GO term of interest appears in the entire dataset (m parameter in phyper), 'hit in population'.
#' "other_gos_in_pop" - How many other GO terms (those not of interest) appear in the entire dataset (n parameter in phyper), 'fail in population'.
#' "total_gos_in_sample" - How many GOs are there in total in the sample (k parameter in phyper), sample size.
#'  
#' @param multitest_correction, character, either empty or one of the options of stats::p.adjust
#' "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param cutoff, numeric, default: 0.05 
#' @param verbose_output, bool, whether to store or not non-significative results
#'
#' @return data.frame with the following columns:
#' "goterm_id", character - Either the identifier or description of the GO term of interest to test for enrichment in sample.
#' "pvalue" - numeric
#' "adj_pvalue" - numeric
#' "sample_id", character - The sample identifier.#' "pvalue" -
#' "note", character - The numbers that respond to the wording 
#' 'There are k hits of the GO of interest in the sample (out of m hits of that GO term in the entire dataset) / out of a grand total of n GO terms annotated in the entire dataset'
#' "p_less_than_cutoff" - YES / NO (whether the adjusted pvalue < cuttoff)
#' @export
#'
#' @examples
#' \dontrun{
#' go_sample_counts <- read.csv("testdata.txt", sep = "\t", header = FALSE)
#' colnames(go_sample_counts) <- c("sample_id", "go_of_interest_in_sample","total_gos_in_sample", "go_of_interest_in_pop", "other_gos_in_pop", "goterm_id")
#' Bonferroni: The pvalues are multiplied by the number of comparisons
#' result <- hypergeom_test_enrichment(go_sample_counts, "bonferroni")
#' }
hypergeom_test_enrichment <- function(go_sample_counts, multitest_correction="", cutoff=0.05, verbose_output=TRUE){
  # q: hitInSample - The GO term of interest #go_of_interest_in_sample
  # m: hitInPop - How many times the GO of interest appears in the entire dataset #go_of_interest_in_pop
  # n: failInPop - How many GOs (not of interest) appear in the entire dataset #other_gos_in_pop
  # k: sampleSize - How many GOs are there in total in the category #total_gos_in_sample
  
  mandatory_columns <- c("goterm_id", "sample_id",
                         "go_of_interest_in_sample", "go_of_interest_in_pop",
                         "other_gos_in_pop", "total_gos_in_sample")
  if (!all(mandatory_columns %in% colnames(go_sample_counts))) {
    stop(paste("The input dataframe must have all of the following columns:",
         mandatory_columns)
         )
  }
  
  max_size <- nrow(go_sample_counts)
  res <- data.frame("goterm_id" = character(max_size),
                    "pvalue" = numeric(max_size),
                    "adj_pvalue" = numeric(max_size),
                    "sample_id" = character(max_size),
                    "note" = character(max_size),
                    "p_less_than_cutoff" = character(max_size))

  for (i in seq_len(nrow(go_sample_counts))) {
    res[i, "goterm_id"] <- as.character(go_sample_counts[i, "goterm_id"])
    res[i, "sample_id"] <- go_sample_counts[i, "sample_id"]
    q <- go_sample_counts[i, "go_of_interest_in_sample"]
    k <- go_sample_counts[i, "total_gos_in_sample"]
    m <- go_sample_counts[i, "go_of_interest_in_pop"]
    n <- go_sample_counts[i, "other_gos_in_pop"]
    res[i, "pvalue"] <- phyper(q - 1, m, n, k, lower.tail = FALSE) # p(X >= q)
    res[i, "note"] <- paste0(q,'(',m,')','/',k)
    res[i, "adj_pvalue"] <- res[i, "pvalue"]
  }
  
  if (multitest_correction != "") { 
    res[, "adj_pvalue"] <- p.adjust(res[, "pvalue"], "fdr")
  }
  
  for (i in seq_len(nrow(go_sample_counts))) {
    if(res[i, "adj_pvalue"] < cutoff){
      res[i, "note"] <- "YES"
    }
  }
  return(res)
}

