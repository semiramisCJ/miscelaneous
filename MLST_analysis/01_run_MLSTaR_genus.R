library(MLSTar)

#Find Cronobacter 
listPubmlst_orgs()[1:50] #cronobacter is number 32

lst <- listPubmlst_orgs() 
mlst_id<-lst[32]

listPubmlst_schemes(org = mlst_id)
#Scheme 1 has 7 genes

# Set genomes to query the database
working_dir<-"/my/path/to/working/directory/"
genomes <- list.files(path = working_dir, 
                      pattern = "\\.fasta$", 
                      full.names = TRUE)

# Send query
x <- doMLST(
  infiles = genomes, # The fasta files
  org = mlst_id, # The organism, in this case is "cronobacter"
  scheme = 1, # Scheme id number
  write = "all", # Only write fasta files for all alleles found
  ddir = paste0(working_dir,'pubmlst_cron_sch1', collapse = '/') # Write downloaded sequences and profiles here.
  )

# Explore object
attributes(x)

# Check result
x$result

# Format table and write result to outfile
tab_res<-cbind(rownames(x$result))
tab_res<-cbind(tab_res, x$result)
rownames(tab_res)<-NULL
colnames(tab_res)<-purrr::prepend(colnames(x$result), "strain")
write.table(tab_res, paste0(working_dir, "genus_MLST_result.tab"), 
			sep="\t", quote=FALSE, row.names=FALSE)


#Save object with results if further exploration is desired
saveRDS(x, paste0(working_dir, "genus_mlst_result.rds"))
