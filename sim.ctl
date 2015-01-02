          seed = -1

       seqfile = seq.txt   * comment out this line if you don't want seqs
*      treefile = tree.tre   * comment out this line if you don't want trees
*      Imapfile = MyImap.txt

  species&tree = 10  A  B  C  D  E  F  G  H  I  J
                     1  1  1  1  1  1  1  1  1  1
   ((((((A #0.01,B #0.01):0.008 #0.01,C #0.01):0.018 #0.01,((D #0.01,E #0.01):0.004 #0.01,F #0.01):0.026 #0.01):0.032 #0.01,(G #0.01,H #0.01):0.022 #0.01):0.036 #0.01,I #0.01):0.038 #0.01,J #0.01):0.04 #0.01;

	loci&length = 1000 51 * number of loci & number of sites at each locus
	locusrate = 1.0     * alpha for gamma for locus rate
