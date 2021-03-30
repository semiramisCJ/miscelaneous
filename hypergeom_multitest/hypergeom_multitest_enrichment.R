##########Function definition
testEnrichment<-function(geneCounts, idIndex, qIndex, kIndex, mIndex, nIndex, 
                         descIndex, multitestCorrection=TRUE, cutoff=0.05){
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
    
    res<-matrix(ncol = 5, nrow = 0)
    for (i in seq_len(nrow(geneCounts))) {
      ID<-as.character(geneCounts[i,idIndex])
      q<-geneCounts[i,qIndex]
      k<-geneCounts[i,kIndex]
      m<-geneCounts[i,mIndex]
      n<-geneCounts[i,nIndex]
      description<-as.character(geneCounts[i,descIndex])
      pval<-phyper(q-1, m, n, k, lower.tail = FALSE) #p(X>=q)
      
      if(multitestCorrection){ adjPval<-pval*nrow(geneCounts)/i }
      else{adjPval<-pval}
      
      if(adjPval < cutoff){
        note<-paste0(q,'(',m,')','/',k)
        res<-rbind(res, c(ID,pval,adjPval,description, note))
      }
    }
    
    return(res)
}


