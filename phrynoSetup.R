##  ----------------------------------------  ##
##                                            ##
##               Phrynosomas                  ##
##           edited: 20 Nov 14                ##
##                                            ##
##  ----------------------------------------  ##

source("~/phrynomics-data/trunk/phrynomicsFunctions.R")
library(phrynomics)


##  ----------------------------------------  ##
##          Writing RAxML Data Files          ##
##  ----------------------------------------  ##


fullfiles <- system("ls ~/Dropbox/UWstuff/phrynomics/Analyses/DataForPaper/full/s*full.phy", intern=TRUE)

setwd("~/Dropbox/UWstuff/phrynomics/Analyses/DataForPaper/vstam")
for(i in sequence(length(fullfiles))){
  name <- strsplit(fullfiles[i], "/")[[1]][length(strsplit(fullfiles[i], "/")[[1]])]
  name <- sub("full", "noAmbigs", name)
  invoc <- paste("Rscript /Users/Barb/phrynomics-data/trunk/Rscripts/RemInvSites.R", fullfiles[i], name, "stam")
  system(invoc)
}
system("cp *noAmbigs.phy ../vlewis")
system("cp *noAmbigs.phy ../nonasc")







##  ----------------------------------------  ##
##            RAxML on the Cluster            ##
##  ----------------------------------------  ##


# Four runs on cluster 1) full analysis with s*full.phy, 2) snps only (no correction) with s*noAmbigs.phy, 3) snps only with lewis correction using s*noAmbigs.phy, and 4) snps only with felstam correction using s*noAmbigs.phy AND parts files

# ## full.pl 
# use strict;
# my $seed1=shift;
# my $seed2=shift;
# my $whichFile=shift;
# open(OUT, ">job.sh") or die;
  # print OUT "#!/bin/sh\n\n";
  # print OUT "#PBS -N \"full".$whichFile."\"\n";
  # print OUT "#PBS -l nodes=1:ppn=16,feature=16core,mem=22gb,walltime=500:00:00\n";
  # print OUT "#PBS -o /gscratch/leache/Barb/jobLogs/\n";
  # print OUT "#PBS -d /gscratch/leache/Barb/full/\n";
  # print OUT "/gscratch/leache/Barb/standard-RAxML-master/raxmlHPC-PTHREADS -T 16 -s ".$whichFile." -f a -m GTRCAT -V -x ".$seed1." -# autoMRE -p ".$seed2." -o GAWI1 -n full_out_".$whichFile."\n";
# close(OUT);

# ## full.sh 
# FILES="s*full.phy"
# RANDOM=3948
# for whichFile in $FILES
# do
  # seed1=$RANDOM
  # seed2=$RANDOM
  # echo $whichFile
  # perl -w full.pl $seed1 $seed2 $whichFile
  # qsub job.sh
# done

# ## nonsasc.pl 
# use strict;
# my $seed1=shift;
# my $seed2=shift;
# my $whichFile=shift;
# open(OUT, ">job.sh") or die;
  # print OUT "#!/bin/sh\n\n";
  # print OUT "#PBS -N \"nonasc".$whichFile."\"\n";
  # print OUT "#PBS -l nodes=1:ppn=16,feature=16core,mem=22gb,walltime=500:00:00\n";
  # print OUT "#PBS -o /gscratch/leache/Barb/jobLogs/\n";
  # print OUT "#PBS -d /gscratch/leache/Barb/nonasc/\n";
  # print OUT "/gscratch/leache/Barb/standard-RAxML-master/raxmlHPC-PTHREADS -T 16 -s ".$whichFile." -f a -m GTRCAT -V -x ".$seed1." -# autoMRE -p ".$seed2." -o GAWI1 -n nonasc_out_".$whichFile."\n";
# close(OUT);

# ## nonasc.sh
# FILES="s*noAmbigs.phy"
# RANDOM=135
# for whichFile in $FILES
# do
  # seed1=$RANDOM
  # seed2=$RANDOM
  # echo $whichFile
  # perl -w nonasc.pl $seed1 $seed2 $whichFile
  # qsub job.sh
# done

# ## lewis.pl
# use strict;
# my $seed1=shift;
# my $seed2=shift;
# my $whichFile=shift;
# open(OUT, ">job.sh") or die;
  # print OUT "#!/bin/sh\n\n";
  # print OUT "#PBS -N \"lewis".$whichFile."\"\n";
  # print OUT "#PBS -l nodes=1:ppn=16,feature=16core,mem=22gb,walltime=500:00:00\n";
  # print OUT "#PBS -o /gscratch/leache/Barb/jobLogs/\n";
  # print OUT "#PBS -d /gscratch/leache/Barb/vlewis/\n";
  # print OUT "/gscratch/leache/Barb/standard-RAxML-master/raxmlHPC-PTHREADS -T 16 -s ".$whichFile." -f a -m ASC_GTRCAT -V --asc-corr=lewis -x ".$seed1." -# autoMRE -p ".$seed2." -o GAWI1 -n lewis_out_".$whichFile."\n";
# close(OUT);

# ## lewis.sh
# FILES="s*noAmbigs.phy"
# RANDOM=456
# for whichFile in $FILES
# do
  # seed1=$RANDOM
  # seed2=$RANDOM
  # echo $whichFile
  # perl -w lewis.pl $seed1 $seed2 $whichFile
  # qsub job.sh
# done

# ## stam.pl
# use strict;
# my $seed1=shift;
# my $seed2=shift;
# my $whichFile=shift;
# open(OUT, ">job.sh") or die;
  # print OUT "#!/bin/sh\n\n";
  # print OUT "#PBS -N \"stam".$whichFile."\"\n";
  # print OUT "#PBS -l nodes=1:ppn=16,feature=16core,mem=22gb,walltime=500:00:00\n";
  # print OUT "#PBS -o /gscratch/leache/Barb/jobLogs/\n";
  # print OUT "#PBS -d /gscratch/leache/Barb/vstam/\n";
  # print OUT "/gscratch/leache/Barb/standard-RAxML-master/raxmlHPC-PTHREADS -T 16 -s ".$whichFile." -f a -m ASC_GTRCAT -V --asc-corr=stamatakis -q ".$whichFile.".part -x ".$seed1." -# autoMRE -p ".$seed2." -o GAWI1 -n stam_out_".$whichFile."\n";
# close(OUT);

# ## stam.sh
# FILES="s*noAmbigs.phy"
# RANDOM=6541
# for whichFile in $FILES
# do
  # seed1=$RANDOM
  # seed2=$RANDOM
  # echo $whichFile
  # perl -w stam.pl $seed1 $seed2 $whichFile
  # qsub job.sh
# done




##  ----------------------------------------  ##
##                     END                    ##
##  ----------------------------------------  ##


##  ----------------------------------------  ##
##      vfel/vstam/vlewis on the Cluster      ##
##  ----------------------------------------  ##

fullfiles <- system("ls ~/Dropbox/UWstuff/phrynomics/Analyses/DataForPaper/full/s*full.phy", intern=TRUE)


setwd("~/Dropbox/UWstuff/phrynomics/Analyses/DataForPaper/snpdata")
fullfiles <- system("ls s*.snps", intern=TRUE)
for(i in sequence(length(fullfiles))){
  name <- strsplit(fullfiles[i], "/")[[1]][length(strsplit(fullfiles[i], "/")[[1]])]
  name <- sub(".snps", "dataonly.snps", name)
  invoc <- paste("Rscript /Users/Barb/phrynomics-data/trunk/Rscripts/RemInvSites.R", fullfiles[i], name, "stam")
  system(invoc)
}
system("cp *noAmbigs.phy ../vlewis")
system("cp *noAmbigs.phy ../nonasc")




################################################
###################### END #####################
################################################































