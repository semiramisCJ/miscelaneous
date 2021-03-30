testEnrichment<-function(go_sample_counts, multitest_correction="", cutoff=0.05, verbose_output=TRUE){
  
  # Multi-test correction methods available:
  # "holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
  # "fdr"
  # q: hitInSample - The GO term of interest #go_of_interest_in_sample
  # m: hitInPop - How many times the GO of interest appears in the entire dataset #go_of_interest_in_pop
  # n: failInPop - How many GOs (not of interest) appear in the entire dataset #other_gos_in_sample
  # k: sampleSize - How many GOs are there in total in the category #total_gos_in_sample
  
  # go_sample_counts - The dataset
  # *Index - The column number with the parameter
  # id* - The category name (HGT/nonHGT)
  # desc* - The description, for example, the GO name
  
  # multitest_correction - Multi-test correction
  # cutoff - Minimum accepted p-value

  # The note will be of the 
  # "annotSpot(annotGenome)/spotSize"
  
  mandatory_columns <- c("goterm_id", "sample_id",
                         "go_of_interest_in_sample", "go_of_interest_in_pop",
                         "other_gos_in_sample", "total_gos_in_sample")
  if(!all(mandatory_columns %in% colnames(go_sample_counts))){
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
    n <- go_sample_counts[i, "other_gos_in_sample"]
    res[i, "pvalue"] <- phyper(q - 1, m, n, k, lower.tail = FALSE) # p(X >= q)
    res[i, "note"] <- paste0(q,'(',m,')','/',k)
    res[i, "adj_pvalue"] <- pval
  }
  
  if(multitest_correction != ""){ 
    res[, "adj_pvalue"] <- p.adjust(res[, "pvalue"], "fdr")
  }
  
  for (i in seq_len(nrow(go_sample_counts))) {
    if(adjPval < cutoff){
      res[i, "note"] <- "YES"
    }
  }
  return(res)
}

go_sample_counts <- read.csv("testdata.txt", sep = "\t", header = FALSE)
colnames(go_sample_counts) <- c("sample_id", "go_of_interest_in_sample","total_gos_in_sample",
                                "go_of_interest_in_pop", "other_gos_in_sample", "goterm_id")
head(go_sample_counts)

# Bonferroni: The pvalues are multiplied by the number of comparisons
result <- testEnrichment(go_sample_counts, "bonferroni")
head(result)

