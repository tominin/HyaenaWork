#!/bin/sh
dp_path=$(head -n 1 './del_me2.txt' | xargs)
bam_path=$(head -n 1 './del_me.txt' | xargs)

chmod +x "$dp_path"/filt_dp_hg*.txt

/home/mmeyer/perlscripts/solexa/analysis/count_informative_positions.pl -pos "$dp_path"/filt_dp_hgA.txt -mask 3 "$bam_path"/*Crocuta_crocuta*.bam > "$dp_path"/hgA_count_postfilt.txt
/home/mmeyer/perlscripts/solexa/analysis/count_informative_positions.pl -pos "$dp_path"/filt_dp_hgB.txt -mask 3 "$bam_path"/*Crocuta_crocuta*.bam > "$dp_path"/hgB_count_postfilt.txt
/home/mmeyer/perlscripts/solexa/analysis/count_informative_positions.pl -pos "$dp_path"/filt_dp_hgC.txt -mask 3 "$bam_path"/*Crocuta_crocuta*.bam > "$dp_path"/hgC_count_postfilt.txt
/home/mmeyer/perlscripts/solexa/analysis/count_informative_positions.pl -pos "$dp_path"/filt_dp_hgD.txt -mask 3 "$bam_path"/*Crocuta_crocuta*.bam > "$dp_path"/hgD_count_postfilt.txt

echo finished working, see final output files at: "$dp_path"

rm ./del_me.txt
rm ./del_me2.txt

