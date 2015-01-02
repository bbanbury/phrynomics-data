#!/bin/sh

for((j=0; j<100; j++)) ; 

do  
../scripts/MCcoalROOT sim.ctl;		

grep -v "^ *$" concatenated.txt > concatenated.phy	

Rscript ../scripts/RemNonRandom.R concatenated.phy seq.txt 51 12 FALSE concat.phy 
 	
Rscript ../scripts/MakeSubData.R concat.phy 4		
 	
raxmlHPC-PTHREADS-AVX -T4 -s s4.txt -f a -m GTRCAT -V -x $RANDOM -\# 2 -p $RANDOM -n s4AllSites;
raxmlHPC-PTHREADS-AVX -T4 -s s6.txt -f a -m GTRCAT -V -x $RANDOM -\# 2 -p $RANDOM -n s6AllSites;
raxmlHPC-PTHREADS-AVX -T4 -s s8.txt -f a -m GTRCAT -V -x $RANDOM -\# 2 -p $RANDOM -n s8AllSites;
cat RAxML_bestTree.s4AllSites >> AllSites_s4.tre
cat RAxML_bestTree.s6AllSites >> AllSites_s6.tre
cat RAxML_bestTree.s8AllSites >> AllSites_s8.tre	

Rscript ../scripts/RemInvSites.R s4.txt s4v.phy stamatakis; 
Rscript ../scripts/RemInvSites.R s6.txt s6v.phy stamatakis; 
Rscript ../scripts/RemInvSites.R s8.txt s8v.phy stamatakis; 
	
raxmlHPC-PTHREADS-AVX -T4 -s s4v.phy -f a -m GTRCAT -V -x $RANDOM -\# 2 -p $RANDOM -n s4Uncorrected;
raxmlHPC-PTHREADS-AVX -T4 -s s6v.phy -f a -m GTRCAT -V -x $RANDOM -\# 2 -p $RANDOM -n s6Uncorrected;
raxmlHPC-PTHREADS-AVX -T4 -s s8v.phy -f a -m GTRCAT -V -x $RANDOM -\# 2 -p $RANDOM -n s8Uncorrected;
cat RAxML_bestTree.s4Uncorrected >> Uncorrected_s4.tre
cat RAxML_bestTree.s6Uncorrected >> Uncorrected_s6.tre
cat RAxML_bestTree.s8Uncorrected >> Uncorrected_s8.tre

raxmlHPC-PTHREADS-AVX -T4 -s s4v.phy -f a -m ASC_GTRCAT -V --asc-corr=lewis -x $RANDOM -\# 2 -p $RANDOM -n s4Conditional;
raxmlHPC-PTHREADS-AVX -T4 -s s6v.phy -f a -m ASC_GTRCAT -V --asc-corr=lewis -x $RANDOM -\# 2 -p $RANDOM -n s6Conditional;
raxmlHPC-PTHREADS-AVX -T4 -s s8v.phy -f a -m ASC_GTRCAT -V --asc-corr=lewis -x $RANDOM -\# 2 -p $RANDOM -n s8Conditional;
cat RAxML_bestTree.s4Conditional >> Conditional_s4.tre
cat RAxML_bestTree.s6Conditional >> Conditional_s6.tre
cat RAxML_bestTree.s8Conditional >> Conditional_s8.tre

raxmlHPC-PTHREADS-AVX -T4 -s s4v.phy -f a -m ASC_GTRCAT -V --asc-corr=stamatakis -q s4v.phy.part -x $RANDOM -\# 2 -p $RANDOM -n s4Reconstituted;
raxmlHPC-PTHREADS-AVX -T4 -s s6v.phy -f a -m ASC_GTRCAT -V --asc-corr=stamatakis -q s6v.phy.part -x $RANDOM -\# 2 -p $RANDOM -n s6Reconstituted;
raxmlHPC-PTHREADS-AVX -T4 -s s8v.phy -f a -m ASC_GTRCAT -V --asc-corr=stamatakis -q s8v.phy.part -x $RANDOM -\# 2 -p $RANDOM -n s8Reconstituted;
cat RAxML_bestTree.s4Reconstituted >> Reconstituted_s4.tre
cat RAxML_bestTree.s6Reconstituted >> Reconstituted_s6.tre
cat RAxML_bestTree.s8Reconstituted >> Reconstituted_s8.tre

Rscript ../scripts/RemInvSites.R s4.txt s4v.phy felsenstein; 
Rscript ../scripts/RemInvSites.R s6.txt s6v.phy felsenstein; 
Rscript ../scripts/RemInvSites.R s8.txt s8v.phy felsenstein; 

raxmlHPC-PTHREADS-AVX -T4 -s s4v.phy -f a -m ASC_GTRCAT -V --asc-corr=felsenstein -q s4v.phy.part -x $RANDOM -\# 2 -p $RANDOM -n s4ReconstitutedF;
raxmlHPC-PTHREADS-AVX -T4 -s s6v.phy -f a -m ASC_GTRCAT -V --asc-corr=felsenstein -q s6v.phy.part -x $RANDOM -\# 2 -p $RANDOM -n s6ReconstitutedF;
raxmlHPC-PTHREADS-AVX -T4 -s s8v.phy -f a -m ASC_GTRCAT -V --asc-corr=felsenstein -q s8v.phy.part -x $RANDOM -\# 2 -p $RANDOM -n s8ReconstitutedF;
cat RAxML_bestTree.s4ReconstitutedF >> ReconstitutedF_s4.tre
cat RAxML_bestTree.s6ReconstitutedF >> ReconstitutedF_s6.tre
cat RAxML_bestTree.s8ReconstitutedF >> ReconstitutedF_s8.tre

grep 'Tree-Length: ' RAxML_info.* >> results.txt

rm ./RAxML_* Imap* SeedUsed concatenated.txt 

done;

Rscript ../scripts/MakeTable.R results.txt output.txt

Rscript ../scripts/MakePlot.R output.txt plotLocusRates1k.pdf

