# Calculate pval given the test statistic and df (chi2)

# Read working directory from command line options
args = commandArgs(trailingOnly=TRUE)
currCase<-args[1]
wDirBase<-args[2]
outDir<-args[3]

# Function definition; this includes contrasts for other models not compared in this script
M0_M3<-function(all_lnL){
  lnL_H0<-all_lnL["M0",1]
  lnL_H1<-all_lnL["M3",1]
  df_test<-4
  testStatistic<-2*(lnL_H1-lnL_H0)
  pval<-pchisq(testStatistic, df=df_test, lower.tail = FALSE)
  return(pval)
}


M1a_M2a<-function(all_lnL){
  lnL_H0<-all_lnL["M1",1]
  lnL_H1<-all_lnL["M2",1]
  df_test<-2
  testStatistic<-2*(lnL_H1-lnL_H0)
  pval<-pchisq(testStatistic, df=df_test, lower.tail = FALSE)
  return(pval)
}


M7_M8<-function(all_lnL){
  lnL_H0<-all_lnL["M7",1]
  lnL_H1<-all_lnL["M8",1]
  df_test<-2
  testStatistic<-2*(lnL_H1-lnL_H0)
  pval<-pchisq(testStatistic, df=df_test, lower.tail = FALSE)
  return(pval)
}



calcAll_05<-function(geneFam, inDir, outDir){
  cutoff<-0.05
  res<-matrix(nrow = 1, ncol = 6)
  
  colnames(res)<-c("contrast","pval", "H0", "H1","rejectH0", "geneFam")
  currFileName<-paste0(inDir, geneFam, ".lnL.tab")
  all_lnL<-read.table(currFileName, sep = "\t", row.names = 1)
  
  
  ##M1a_M2a
  answer<-"no"
  p<-M1a_M2a(all_lnL)
  if(!is.na(p) && p < cutoff){
    answer<-"yes"
  }
  res[1,]<-c("M1a_M2a", p, "neutral", "selection", answer, geneFam)
  
  ##Write results
  write.table(res, paste0(outDir, "LRTs_05.tab"), append = TRUE, sep = "\t", 
              row.names = FALSE, col.names = FALSE, quote = FALSE)

  return()
}

calcAll_001<-function(geneFam, inDir, outDir){
  cutoff<-0.001
  res<-matrix(nrow = 1, ncol = 6)
  
  colnames(res)<-c("contrast","pval", "H0", "H1","rejectH0", "geneFam")
  currFileName<-paste0(inDir, geneFam, ".lnL.tab")
  all_lnL<-read.table(currFileName, sep = "\t", row.names = 1)
  

  ##M1a_M2a
  answer<-"no"
  p<-M1a_M2a(all_lnL)
  if(!is.na(p) && p < cutoff){
    answer<-"yes"
  }
  res[1,]<-c("M1a_M2a", p, "neutral", "selection", answer, geneFam)
  
  ##Write results
  write.table(res, paste0(outDir, "LRTs_001.tab"), append = TRUE, sep = "\t", 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  return()
}


# Iterate through each gene family for each case
analyzeResults<-function(currCase, wDirBase, outDir){
  # Declare directories and location of the gene family list 
  wDir<-paste0(wDirBase, currCase, "/")
  gfFile<-paste0(wDirBase,currCase,"/noREC.",currCase,".gf.list")
  geneFamilies<-read.table(gfFile)
  apply(geneFamilies, 1, calcAll_05, wDir, outDir)
  apply(geneFamilies, 1, calcAll_001, wDir, outDir)
}

# Run function
analyzeResults(currCase, wDirBase, outDir)

