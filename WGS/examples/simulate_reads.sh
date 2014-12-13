# *************************************************************
#
# $Source: $
# $Revision: $                                                                 
# $State: $                                                                     
# $Date: $                                                      
# $Author: $  
#
# $Log: $
#
#
# *************************************************************

# simulates 1000 paired-end reads from ../data/truncated_GCA_000001405.15_GRCh38_no_alt_analysis_set.fna 
# and writes them to ../data/sim_reads_1.fq.gz and to ../data/sim_reads_2.fq.gz
wgsim -N1000 -S1 ../data/truncated_GCA_000001405.15_GRCh38_no_alt_analysis_set.fna ../data/sim_reads_1.fq.gz ../data/sim_reads_2.fq.gz



