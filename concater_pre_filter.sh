#!/bin/bash
###NOTES
### argument 1 = where user has collated all their bam files they want to check for Hyaena
### argument 2 = where the user wants to output the intermediate files to

/home/mmeyer/perlscripts/solexa/analysis/count_informative_positions.pl -pos /mnt/scratch/thomas_harris_snell/hyena_mtDNA/hyena_diagnostic_positions/non_filtered_diag_pos/tom/dp_hgA.txt -mask 3 $1/*Crocuta_crocuta*.bam > $2/hgA_prefilt.txt
/home/mmeyer/perlscripts/solexa/analysis/count_informative_positions.pl -pos /mnt/scratch/thomas_harris_snell/hyena_mtDNA/hyena_diagnostic_positions/non_filtered_diag_pos/tom/dp_hgB.txt -mask 3 $1/*Crocuta_crocuta*.bam > $2/hgB_prefilt.txt
/home/mmeyer/perlscripts/solexa/analysis/count_informative_positions.pl -pos /mnt/scratch/thomas_harris_snell/hyena_mtDNA/hyena_diagnostic_positions/non_filtered_diag_pos/tom/dp_hgC.txt -mask 3 $1/*Crocuta_crocuta*.bam > $2/hgC_prefilt.txt
/home/mmeyer/perlscripts/solexa/analysis/count_informative_positions.pl -pos /mnt/scratch/thomas_harris_snell/hyena_mtDNA/hyena_diagnostic_positions/non_filtered_diag_pos/tom/dp_hgD.txt -mask 3 $1/*Crocuta_crocuta*.bam > $2/hgD_prefilt.txt

chmod +x $2/hgA_prefilt.txt
chmod +x $2/hgB_prefilt.txt
chmod +x $2/hgC_prefilt.txt
chmod +x $2/hgD_prefilt.txt

echo $1 > ./del_me.txt
echo $2
