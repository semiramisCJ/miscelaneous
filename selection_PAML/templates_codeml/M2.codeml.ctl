       seqfile = seqfile.txt          * sequence data filename
      treefile = treefile.txt         * SET THIS for tree file with ML branch lengths under M0
      outfile = results.txt           * main result file name

        noisy = 9                     * lots of rubbish on the screen
      verbose = 1                     * detailed output
      runmode = 0                     * user defined tree
      seqtype = 1                     * codons
    CodonFreq = 2                     * F3X4 for codon ferquencies
        model = 0                     * one omega ratio for all branches

      NSsites = 2                     * SET THIS for M2

        icode = 0                     * universal code
    fix_kappa = 0                     * kappa fixed
        kappa = 2                     * RESET THIS for all models

    fix_omega = 0                     * omega to be estimated 
        omega = 0.1                   * initial omega

  fix_blength = 0                     * fixed branch lengths from tree file
