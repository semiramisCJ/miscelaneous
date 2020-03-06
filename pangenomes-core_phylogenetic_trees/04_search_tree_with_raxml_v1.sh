# Perform tree search with RAxML with a previously chosen model

# "allCoreGenes.fas" is the concantenated alignment file with which we will be working

# Go to working directory 
raxmlDir="/my/path/to/working/directory"
mkdir "$raxmlDir"
cd "$raxmlDir"

# Run multiple searches to get the best scoring tree
raxmlHPC-PTHREADS-SSE3 -s allCoreGenes.fas -n AllCoreGenes -m GTRGAMMAI -p 123456789 -N 50 -T 5

# Perform bootstrap search
raxmlHPC-PTHREADS-SSE3 -s allCoreGenes.fas -n AllCoreGenesBS -m GTRGAMMAI -p 123456789 -b 123456789 -N 100 -T 5

# Add bootstrap values on best tree
raxmlHPC-PTHREADS-SSE3 -m GTRGAMMAI -p 123456789 -f b -t RAxML_bestTree.AllCoreGenes -z RAxML_bootstrap.AllCoreGenesBS -n AllCoreGenes_wBootstrap -T 5

# Rename resulting tree
mv RAxML_bipartitionsBranchLabels.AllCoreGenes_wBootstrap AllCoreGenes_RAxML.tree


