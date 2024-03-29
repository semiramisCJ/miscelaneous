% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hypergeom_multitest_enrichment.R
\name{hypergeom_test_enrichment}
\alias{hypergeom_test_enrichment}
\title{Test for enrichment in a dataset with the hypergeometric test}
\usage{
hypergeom_test_enrichment(
  go_sample_counts,
  multitest_correction = "",
  cutoff = 0.05,
  verbose_output = TRUE
)
}
\arguments{
\item{go_sample_counts, }{data.frame with the following columns:
"goterm_id", character - Either the identifier or description of the GO term of interest to test for enrichment in sample.
"sample_id", character - The sample identifier.
"go_of_interest_in_sample", integer - Number of hits of the GO term of interest in the sample (q parameter in phyper), 'hit in sample'.
"go_of_interest_in_pop" - How many times the GO term of interest appears in the entire dataset (m parameter in phyper), 'hit in population'.
"other_gos_in_pop" - How many other GO terms (those not of interest) appear in the entire dataset (n parameter in phyper), 'fail in population'.
"total_gos_in_sample" - How many GOs are there in total in the sample (k parameter in phyper), sample size.}

\item{multitest_correction, }{character, either empty or one of the options of stats::p.adjust
"holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"}

\item{cutoff, }{numeric, default: 0.05}

\item{verbose_output, }{bool, whether to store or not non-significative results}
}
\value{
data.frame with the following columns:
"goterm_id", character - Either the identifier or description of the GO term of interest to test for enrichment in sample.
"pvalue" - numeric
"adj_pvalue" - numeric
"sample_id", character - The sample identifier.#' "pvalue" -
"note", character - The numbers that respond to the wording
'There are k hits of the GO of interest in the sample (out of m hits of that GO term in the entire dataset) / out of a grand total of n GO terms annotated in the entire dataset'
"p_less_than_cutoff" - YES / NO (whether the adjusted pvalue < cuttoff)
}
\description{
Test for enrichment in a dataset with the hypergeometric test
}
\examples{
\dontrun{
go_sample_counts <- read.csv("testdata.txt", sep = "\t", header = FALSE)
colnames(go_sample_counts) <- c("sample_id", "go_of_interest_in_sample","total_gos_in_sample", "go_of_interest_in_pop", "other_gos_in_pop", "goterm_id")
Bonferroni: The pvalues are multiplied by the number of comparisons
result <- hypergeom_test_enrichment(go_sample_counts, "bonferroni")
}
}
