testEnrichment<-function(geneCounts, idIndex, qIndex, kIndex, mIndex, nIndex, 
                       descIndex, multitestCorrection="", cutoff=0.05){
  # Multi-test correction methods available:
  # "holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
  # "fdr"
  # q: hitInSample - The GO term of interest
  # m: hitInPop - How many times the GO of interest appears in the entire dataset
  # n: failInPop - How many GOs (not of interest) appear in the entire dataset
  # k: sampleSize - How many GOs are there in total in the category
  
  # geneCounts - The dataset
  # *Index - The column number with the parameter
  # id* - The category name (HGT/nonHGT)
  # desc* - The description, for example, the GO name
  
  # multitestCorrection - Multi-test correction
  # cutoff - Minimum accepted p-value

  # The note will be of the 
  # "annotSpot(annotGenome)/spotSize"
  
  res<-matrix(ncol = 6, nrow = nrow(geneCounts))
  colnames(res) <- c("ID","pval","adjPval","description", "note", "signif")
  for (i in seq_len(nrow(geneCounts))) {
    ID<-as.character(geneCounts[i,idIndex])
    q<-geneCounts[i,qIndex]
    k<-geneCounts[i,kIndex]
    m<-geneCounts[i,mIndex]
    n<-geneCounts[i,nIndex]
    description<-as.character(geneCounts[i,descIndex])
    pval<-phyper(q-1, m, n, k, lower.tail = FALSE) #p(X>=q)
    note<-paste0(q,'(',m,')','/',k)
    adjPval <- pval
    res[i,] <- c(ID,pval,adjPval,description, note, "")
  }
  
  if(multitestCorrection != ""){ 
    res[,3]<-p.adjust(res[,2], "fdr")
  }
  
  for (i in seq_len(nrow(geneCounts))) {
    if(adjPval < cutoff){
      res[i, 6] <- "YES"
    }
  }
  return(res)
}

geneCounts <- read.csv("testdata.txt", sep = "\t")
head(geneCounts)

# Bonferroni: The pvalues are multiplied by the number of comparisons
result <- testEnrichment(geneCounts, 1, 2, 3, 4, 5, 6, "bonferroni")
head(result)

