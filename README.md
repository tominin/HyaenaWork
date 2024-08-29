##HyaenaWork##
Full workings of creating my hyaena diagnostic positions, with some scripts that the user can use or alter to there needs.

##STEP 1##
Sequence Alignment:
Download all files in Hyena_datafile.docx (from NCBI):
Stack seq's in BioEdit, ensure NC_020671.1 is top (for coordinates later)
COMMAND LINE ARG: (install mafft prior)
mafft file > mafft_joint.fas

Open mafft_joint.fas
Cut and paste under/overlaps at the end and beginning of reads due to
mitochondria's circular nature.
Add insertions in genomes that dont align quite right (only need to do a few)

From this file get all the different groups (Haplogroup A, Haplogroup B, Haplogroup C, Haplogroup D, Proteles cristata, Striped Hyaena, Brown Hyaena) seperated into seperate files for
mmeyer script. All still need the NC_020671.1 at the top of the alignment.

##STEP 2##
UNIX command to define diagnostic pos for each haplogroup:
(Example)
/home/mmeyer/perlscripts/solexa/analysis/define_informative_positions_generic.pl -definitions A,Brown,1,/r1/people/thomas_harris_snell/hyena_mtDNA/brown_hyena/mafft_brown.fas:A,Striped,1,/r1/people/thomas_harris_snell/hyena_mtDNA/striped_hyena/mafft_striped.fas:A,Proteles,1,/r1/people/thomas_harris_snell/hyena_mtDNA/proteles/mafft_proteles.fas:A,Haplo_B,1,/r1/people/thomas_harris_snell/hyena_mtDNA/crocuta/mafft_hpgB:A,Haplo_C,1,/r1/people/thomas_harris_snell/hyena_mtDNA/crocuta/mafft_hpgC.fas:A,Haplo_D,1,/r1/people/thomas_harris_snell/hyena_mtDNA/crocuta/mafft_hpgD:B,Haplo_A,1,/r1/people/thomas_harris_snell/hyena_mtDNA/crocuta/mafft_hpgA.fas > /r1/people/thomas_harris_snell/hyena_mtDNA/hyena_dps/pre_edit/hpgA_prefilt.txt

Important to note that the derived group should be in B, not A!

##STEP 3##
Filtering out over-represented diagnostic positions
Create a bam directory with all selected bam files that you would like to run and analyse for Crocuta crocuta species. For this when cp the files from the sediment directory ensure to add a uniq letter for each run ID as there can be duplicate CapLibID/IndexLibID's therefore removing some of your desired data.
Then you are all ready to run the filtering scripts:

1ST COMMAND: bash concater_pre_filt.sh [path to bam files] [path to a directory that holds the intermediate files made] | python filt_dps_script.py [path to where the final output of your script goes]

2ND COMMAND: bash concater_pos_filt.sh

(I was having problems merging these 3 scripts together, a better programmer should be able to easily do it if you so wish)

This will output: The filtered diagnostic pos file of all 4 haplogroups The filtered count reads from the selected bam files

Please note that the filtered pos should be done every time you analyse a new set of data! (HYAENA ONLY)


##STEP 4##
Post running analysis, I have submitted my code for the user to read through, as long as they have the relevant datafiles with their output from mmeyer's scripts they should be able to quite easily alter my code to their need. Some changes with how layers are seperated may need to be done and different syntax readings may need to be taken into account.
Script = quicksand.automated.py

NOTES
Filters used in script:
Removes any samples that had lower than 3 hits for the group we want to investigate

Removes any samples that fall below the binomial distribution 95% confidence interval for the lowest bound at 10% coverage

Ensure to use time gaps from file #layerdates_ZJ.txt

Remove duplicate sample IDs

Count the number of samples (NOT THE NUMBER OF SEQUENCES IN EACH BAM FILE) that showed significant occurance rates

Obtain total percentage for the count of each haplogroups for each layer

To know which sample is in which layer use file #layers_ZJ_04012021.txt


