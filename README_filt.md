# HyaenaWork
Details for how to run my filtering code at the command line, as well as the logic used to create the graphic output if someone would ever like to recreate this procedure.

After running mmeyer's define_diagnostic positions script previously we need to remove any diagnostic positions that are over represented due to mitochondrial capture bias.

Important!
Ensure that you have created the derived state as group B and the ancestral as group A when you ran the code

Create a bam directory with all selected bam files that you would like to run and analyse for Crocuta crocuta species. For this when cp the files from the sediment directory ensure to add a uniq letter for each run ID as there can be duplicate CapLibID/IndexLibID's therefore removing some of your desired data.

Then you are all ready to run the filtering scripts:

1ST COMMAND:

bash concater_pre_filt.sh [path to bam files] [path to a directory that holds the intermediate files made] | python filt_dps_script.py [path to where the final output of your script goes]


2ND COMMAND:
bash concater_pos_filt.sh


(I was having problems merging these 3 scripts together, a better programmer should be able to easily do it if you so wish)


This will output:
The filtered diagnostic pos file of all 4 haplogroups
The filtered count reads from the selected bam files

Please note that the filtered pos should be done every time you analyse a new set of data!


