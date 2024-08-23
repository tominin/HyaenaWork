#HyaenaWork
Details showing the logic and filtering criteria used in my python script to hopefully help someone recreate this for future.

List of Files quicksand_automated.py reads in:
20240618_data_for_tom_DC_MC_SE-Profile_AA75_AllSamples_quicksandFilter.csv : File RunID, IndexID, CapLibID, SampleID and what family is found there from previous analysis (MEGAN/quicksand)

hom_curated_input19062020.txt : Helps remove any SampleID's that isn't from 1st screening (used for mammalia capture)

layers_ZJ_04012021.txt : Marker, Layer info for each sample ID

layerdates_ZJ.txt : what each is each markers time interval

hgA_count_postfilt.txt
hgB_count_postfilt.txt
hgC_count_postfilt.txt
hgD_count_postfilt.txt

The first four files will need to be changed for your specific cave of interest. As long as they have the appropiate information it should require little tinkering, depending on the information.
