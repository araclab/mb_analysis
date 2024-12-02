#!/bin/sh
#SBATCH --time 1-00:00:00
#SBATCH -p short
#SBATCH -o dada2_%j.out
#SBATCH -e dada2_%j.err

STARTTIMER="$(date +%s)"

# Set your conda enviroment location here as well
export scripts='YOUR_PATH/dada2_pipeline_scripts'
export rdp='YOUR_PATH/rdp_classifier_2.12/dist'
export v6species='YOUR_PATH/V6_species'
export database='YOUR_PATH/dada2_pipeline_databases'
export HDF5_USE_FILE_LOCKING='FALSE'

usage() {
	echo "
--------------------------------------------------------------------------------------------------------------------------------------

	USAGE: options available
		mandatory options:

				-i or --input : provide the name of your input folder (containing FASTQ files)
		                
		optional options:
				-p or --primer : path to fasta file with primer sequences (default:${database}/amplicon_primers_bactquant.fasta) 
				-r or --run : provide the name of your sequencing run e.g.: run17
				-t or --threads : give the number of threads to use (default is 16)
				-dt or --taxonomydb : give the name of your taxonomy database
				-ds or --speciesdb : give the name of your species database
	Help options:
				-h or --help : pop up the usage guide

		e.g. :

				sbatch dada2_pipeline_vP_v11_cleaned.sh -i raw_reads -r run17 -t 2 -dt rdp_train_set_16.fa.gz -ds rdp_species_assignment_16.fa.gz -p amplicon_primers_fadrosh.fasta
				raw_reads is the name of the input folder
				run17 is the run name
				2 is the total number of threads
				rdp_train_set_16.fa.gz is the name of your taxonomy database
				rdp_species_assignment_16.fa.gz is the name of your species database
				amplicon_primers_fadrosh.fasta is the name of your fadrosh amplicon primer sequences fasta file

--------------------------------------------------------------------------------------------------------------------------------------

	"
	}


# Give 2 variables at least
if [ $# -lt 2 ]
then
	echo "--------------AT LEAST give 2 arguments ( -i input_folder ) to represent input_name--------------"
	usage
	exit 1
fi


# Default settings
threads=16
run_name=run_unknown
db_taxonomy=rdp_train_set_16.fa.gz
db_species=rdp_species_assignment_16.fa.gz
primer_path="${database}/amplicon_primers_bactquant.fasta"

# Flag options
nopts=$#
for ((i=1 ; i <= nopts ; i++)); do
	case $1 in
	-r | --run)
	run_name=$2
			echo "------------------------------ This is run: $2 ------------------------------"
			shift 2
			;;

	-i | --input)
	input_folder=$2
			echo "------------------------------ Input_folder: $2 ------------------------------"
			shift 2
			;;
	-p | --primer)
	primer_path=$2
			echo "------------------------------ Primer fasta path: $2 ------------------------------"
			shift 2
			;;
	-t | --threads)
	threads=$2
         		echo "------------------------------ The number of our threads: $2 ------------------------------"
  	        	shift 2
			;;
	-dt | --taxonomydb)
	db_taxonomy=$2
			echo "------------------------------ Taxonomy Databases you are using: $2 ------------------------------"
			shift 2
			;;
	-ds | --speciesdb)
	db_species=$2
			echo "------------------------------ Species Databases you are using: $2 ------------------------------"
			shift 2
			;;
	-h | --help)
			usage
			exit 1
			;;
	esac
done


if [ -z "$input_folder" ]
then
	echo "------------------------------ Can not go on without your input ------------------------------"
	echo "------------------------------ AT LEAST give 2 arguments ( -i or --input ) to represent input_folder ------------------------------"
	usage
	exit 1
fi




export primer_path
my_output=my_output_$run_name
mkdir $my_output
mkdir $my_output/End_user_results
mkdir $my_output/RDS_files
mkdir trimmed_reads
mkdir primer_trimmed_reads



# Counting raw reads from unfiltered .gz files
echo "------------------------------------------------------------ STEP 1: ORIGINAL RAW DATA COUNTING ------------------------------------------------------------"
for i in $input_folder/*_R1_*.fastq.gz; do echo -n "$i" $'\t'; echo $(zcat $i | wc -l) / 4 | bc; done > $my_output/LP_${run_name}_readCount.txt

# Text editing: delete prefix and suffix
awk '{split($1,a,"/"); $1=a[length(a)]"\t";}1' $my_output/LP_${run_name}_readCount.txt > $my_output/LP_${run_name}_readCount.intermediate1.txt
awk '{split($1,a,"_"); $1=a[1]"\t";}1' $my_output/LP_${run_name}_readCount.intermediate1.txt > $my_output/LP_${run_name}_readCount.intermediate2.txt

# Text editing: add headers to each column
echo -e "Sample\tRaw_reads" | cat - $my_output/LP_${run_name}_readCount.intermediate2.txt > $my_output/LP_${run_name}_readCount.final.txt

rm $my_output/LP_${run_name}_readCount.intermediate1.txt
rm $my_output/LP_${run_name}_readCount.intermediate2.txt
rm $my_output/LP_${run_name}_readCount.txt



# Run Trimmomatic to remove low quality sequences from raw data (input_file_name)
echo "------------------------------------------------------------ STEP 2: TRIM LOW QUALITY SEQUENCES ------------------------------------------------------------"
perl $scripts/run_trimmomatic_amplicon_mode.pl $input_folder/*



# Move files that have "paired_trimmed" in their file name to trimmed_reads
echo "------------------------------------------------------------ STEP 3: MOVE FILES ------------------------------------------------------------"
mv $input_folder/*paired_trimmed* trimmed_reads
rm $input_folder/*single_trimmed*



# Primer sequences need to be removed so that DADA2 won't generate inaccurate data or chokes on ambigious bases
echo "------------------------------------------------------------ STEP 4: REMOVE PRIMER SEQUENCES ------------------------------------------------------------"
conda activate dada2

remove_adapter() {
f1="$(basename -- $1)"
f2="$(basename -- $2)"

echo "cutadapt --pair-filter any --no-indels -e 0.15 --discard-untrimmed -g file:${primer_path} -G file:${primer_path} -Z -o primer_trimmed_reads/${f1}.gz -p primer_trimmed_reads/${f2}.gz ${1} ${2} > primer_trimmed_reads/${f1%%.fastq*}_cutadapt_log.txt"

cutadapt --pair-filter any --no-indels -e 0.15 --discard-untrimmed -g file:${primer_path} -G file:${primer_path} -Z -o primer_trimmed_reads/${f1}.gz -p primer_trimmed_reads/${f2}.gz ${1} ${2} > primer_trimmed_reads/${f1%%.fastq*}_cutadapt_log.txt
} 
export -f remove_adapter;
parallel --link --jobs $threads remove_adapter ::: trimmed_reads/*_R1_*.fastq ::: trimmed_reads/*_R2_*.fastq


# Filter out reads with ambigious base calls (Ns)
echo "------------------------------------------------------------ STEP 5: FILTER OUT READS WITH AMBIGIOUS BASE CALLS (Ns) ------------------------------------------------------------"
$scripts/dada2_filter.R -f primer_trimmed_reads --maxN 0 --truncQ 0 --threads $threads --truncLen 0 --maxEE 2,2 --f_match _R1_.*fastq.gz --r_match _R2_.*fastq.gz -o $my_output/${run_name}_filtered_fastqs



# DADA2 ASV calling - creates seqtab.rds file - THIS TAKES A WHILE
echo "------------------------------------------------------------ STEP 6: CREATES SEQTAB.RDS FILE, HOLD ON ------------------------------------------------------------"
$scripts/dada2_inference.R -f $my_output/${run_name}_filtered_fastqs/ --seed 4124 -t $threads --verbose --plot_errors --minOverlap 12 -o $my_output/${run_name}_seqtab.rds



# Remove spurious features and assign taxonomy - creates seqtab_final.rds
echo "------------------------------------------------------------ STEP 7: NOW WE ARE FILTERING SEQTAB.RDS, THIS STEP TAKES THE LONGEST ------------------------------------------------------------"
$scripts/dada2_chimera_taxa.R -i $my_output/${run_name}_seqtab.rds -t 2 -r $database/$db_taxonomy -s $database/$db_species --allowMultiple --minBoot 80 --count_out $my_output/${run_name}_seqtab_final.rds --sp_out $my_output/${run_name}_species_final.rds --tax_out $my_output/${run_name}_tax_final.rds



# Convert RDS objects to TSV (in classic biom format) and FASTA - creates file.biom.tsv and file_species.fasta
echo "------------------------------------------------------------ STEP 8: CREATING MY_OUTPUT FILES . . . ALMOST DONE ------------------------------------------------------------"
$scripts/convert_dada2_out.R -i $my_output/${run_name}_seqtab_final.rds -b $my_output/${run_name}_biom.tsv -f $my_output/${run_name}_species.fasta --taxa_in $my_output/${run_name}_species_final.rds --taxa_out $my_output/${run_name}_taxa_metadata.txt



echo "------------------------------------------------------------ STEP 9: CONVERT TSV TO BIOM FILE ------------------------------------------------------------"
biom convert -i $my_output/${run_name}_biom.tsv -o $my_output/${run_name}_json.biom --to-json



# Create the summary of your biom file
echo "------------------------------------------------------------ STEP 10: CREATE THE SUMMARY OF YOUR BIOM FILE ------------------------------------------------------------"
biom summarize-table -i $my_output/${run_name}_json.biom -o $my_output/${run_name}_json_summary.txt



# Create log file, text editing: extract columns we need, calculate the summation
echo "------------------------------------------------------------ STEP 11: CREATE LOG FILE ------------------------------------------------------------"

awk 'FNR == 1 {next} NR==FNR {h[$1] = $2; next} {OFS = "\t" ; print $1,$2,h[$1]}' dada2_filter_read_counts.txt $my_output/LP_${run_name}_readCount.final.txt > $my_output/${run_name}_final1.txt
awk 'NR==FNR {h[$1] = $3; next} {OFS = "\t" ; print $0,h[$1]}' dada2_filter_read_counts.txt $my_output/${run_name}_final1.txt > $my_output/${run_name}_final2.txt
awk 'NR==FNR {h[$1] = $7; next} {OFS = "\t" ; print $0,h[$1]}' dada2_inferred_read_counts.txt $my_output/${run_name}_final2.txt > $my_output/${run_name}_final3.txt
awk 'NR==FNR {h[$1] = $3; next} {OFS = "\t" ; print $0,h[$1]}' dada2_nonchimera_counts.txt $my_output/${run_name}_final3.txt > $my_output/${run_name}_final4.txt
echo -e "Sample\tRaw_reads\tQuality_trimmed_reads\tReads_without_'N'_bases\tDenoised_dereplicated_reads\tNonchimera_reads" | cat - $my_output/${run_name}_final4.txt > $my_output/${run_name}_dada2_final5.txt

awk 'NR==1{$7="Proportion_of_saved_reads"}NR>1{$7 = (($6)/($2))*100"\%"} {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' $my_output/${run_name}_dada2_final5.txt > $my_output/${run_name}_dada2_final6.txt
awk '{print $0}''NR>1{sum2+=$2; sum3+=$3; sum4+=$4; sum5+=$5; sum6+=$6}END{print "Total_reads\t"sum2"\t"sum3"\t"sum4"\t"sum5"\t"sum6"\t"(sum6/sum2)*100"\%"}' $my_output/${run_name}_dada2_final6.txt > $my_output/${run_name}_dada2_QC_stats.txt

# Remove intermediate files
rm $my_output/${run_name}_dada2_final1.txt
rm $my_output/${run_name}_dada2_final2.txt
rm $my_output/${run_name}_dada2_final3.txt
rm $my_output/${run_name}_dada2_final4.txt
rm $my_output/${run_name}_dada2_final5.txt
rm $my_output/${run_name}_dada2_final6.txt


# Move QC files to main output folder
mv dada2_filter_read_counts.txt $my_output
mv dada2_inferred_read_counts.txt $my_output
mv dada2_nonchimera_counts.txt $my_output
mv estimated_forward_err.pdf $my_output/${run_name}_estimated_forward_err.pdf
mv estimated_reverse_err.pdf $my_output/${run_name}_estimated_reverse_err.pdf
mv $my_output/*.rds $my_output/RDS_files



# Summary
echo "------------------------------------------------------------ STEP 12: SUMMARY ------------------------------------------------------------"
# Genus level classification
java -Xmx12g -jar $rdp/classifier.jar classify -o $my_output/RDP212_GENUS_taxa_metadata.txt -q $my_output/${run_name}_species.fasta -f fixrank

# Species level classification using a custom classifier
java -Xmx12g -jar $rdp/classifier.jar classify -o $my_output/RDP212_V6species_taxa_metadata.txt -q $my_output/${run_name}_species.fasta -t $v6species/rRNAClassifier.properties

conda deactivate


perl $scripts/convert_default_settings_rdpclassifier_genus_taxa_to_biom_metadata_v043021.pl $my_output/RDP212_GENUS_taxa_metadata.txt
perl $scripts/convert_default_settings_rdpclassifier_speciesV6_taxa_to_biom_metadata_v043021.pl $my_output/RDP212_V6species_taxa_metadata.txt


conda activate dada2

biom add-metadata -i $my_output/${run_name}_json.biom -o $my_output/${run_name}_RDP212_GENUS_c80_tax.biom --observation-metadata-fp $my_output/RDP212_GENUS_taxa_metadata.txt.c80_BIOMversion.tsv --observation-header OTUID,taxonomy --sc-separated taxonomy --output-as-json

biom add-metadata -i $my_output/${run_name}_json.biom -o $my_output/${run_name}_RDP212_V6-SPECIES_c80_tax.biom --observation-metadata-fp $my_output/RDP212_V6species_taxa_metadata.txt.c80_BIOMversion.tsv --observation-header OTUID,taxonomy --sc-separated taxonomy --output-as-json

conda deactivate

conda activate qiime1

summarize_taxa.py -i $my_output/${run_name}_RDP212_GENUS_c80_tax.biom -L 2,3,4,5,6 -o $my_output/${run_name}_RDP212_GENUS_c80_tax_summarizeTaxa --absolute_abundance
summarize_taxa.py -i $my_output/${run_name}_RDP212_V6-SPECIES_c80_tax.biom -L 2,3,4,5,6,7 -o $my_output/${run_name}_RDP212_V6-SPECIES_c80_tax_summarizeTaxa --absolute_abundance

conda deactivate

conda activate dada2

biom convert -i $my_output/${run_name}_RDP212_GENUS_c80_tax.biom -o $my_output/${run_name}_RDP212_GENUS_c80_tax.tsv --to-tsv --header-key taxonomy
biom convert -i $my_output/${run_name}_RDP212_V6-SPECIES_c80_tax.biom -o $my_output/${run_name}_RDP212_V6-SPECIES_c80_tax.tsv --to-tsv --header-key taxonomy

conda deactivate



# Aureus merging
echo "------------------------------------------------------------ STEP 13: AUREUS MERGING ------------------------------------------------------------"
# Data summation of bacteria that has "Staphylococcus_aureus<0.8"
awk -F '\t' '$1 ~ /Staphylococcus_aureus/&&/<0.8/ {x=NF; for (i=1;i<x;i++) {sum[i] += $(i+1)} } END { for (i=1;i<x;i++) {print sum[i]} }' $my_output/${run_name}_RDP212_V6-SPECIES_c80_tax_summarizeTaxa/${run_name}_RDP212_V6-SPECIES_c80_tax_L7.txt > L7_intermediate1.txt

# Change my_output record separator to tab instead of defalut new line
awk 'BEGIN {ORS="\t"}; {print $0}' L7_intermediate1.txt > L7_intermediate2.txt

# Adding "Staphylococcus_aureus<0.8" before the data
awk '{print "Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus_aureus<0.8\t"$0}' L7_intermediate2.txt > L7_intermediate3.txt

# Data summation of bacteria that has "Staphylococcus_aureus>=0.8"
awk -F '\t' '$1 ~ /Staphylococcus_aureus/&&!/<0.8/ {x=NF; for (i=1;i<x;i++) {sum[i] += $(i+1)} } END { for (i=1;i<x;i++) {print sum[i]} }' $my_output/${run_name}_RDP212_V6-SPECIES_c80_tax_summarizeTaxa/${run_name}_RDP212_V6-SPECIES_c80_tax_L7.txt > L7_intermediate4.txt

# Change my_output record separator to tab instead of defalut new line
awk 'BEGIN {ORS="\t"}; {print $0}' L7_intermediate4.txt > L7_intermediate5.txt

# Adding "Staphylococcus_aureus<0.8" before the data
awk '{print "Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus_aureus\t"$0}' L7_intermediate5.txt > L7_intermediate6.txt

# Delete bacteria that has "Staphylococcus_aureus" from the original V6-SPECIES file
awk '!/Staphylococcus_aureus/{print $0}' $my_output/${run_name}_RDP212_V6-SPECIES_c80_tax_summarizeTaxa/${run_name}_RDP212_V6-SPECIES_c80_tax_L7.txt > L7_intermediate7.txt

# Add the generated Staphylococcus_aureus line(s) to the end of the original V6-SPECIES file
cat L7_intermediate7.txt L7_intermediate3.txt L7_intermediate6.txt > $my_output/${run_name}_RDP212_V6-SPECIES_c80_tax_summarizeTaxa/${run_name}_RDP212_V6-SPECIES_c80_tax_L7_aureus_final.txt

rm L7_intermediate1.txt
rm L7_intermediate2.txt
rm L7_intermediate3.txt
rm L7_intermediate4.txt
rm L7_intermediate5.txt
rm L7_intermediate6.txt
rm L7_intermediate7.txt


# Copy user files to separate folder
cp $my_output/${run_name}_RDP212_GENUS_c80_tax_summarizeTaxa/${run_name}_RDP212_GENUS_c80_tax_L6.txt $my_output/End_user_results
cp $my_output/${run_name}_RDP212_V6-SPECIES_c80_tax_summarizeTaxa/${run_name}_RDP212_V6-SPECIES_c80_tax_L7.txt $my_output/End_user_results
cp $my_output/${run_name}_RDP212_V6-SPECIES_c80_tax_summarizeTaxa/${run_name}_RDP212_V6-SPECIES_c80_tax_L7_aureus_final.txt $my_output/End_user_results
cp $my_output/${run_name}_json_summary.txt $my_output/End_user_results
cp $my_output/${run_name}_dada2_QC_stats.txt $my_output/End_user_results
cp $my_output/${run_name}_RDP212_GENUS_c80_tax.tsv $my_output/End_user_results
cp $my_output/${run_name}_species.fasta $my_output/End_user_results



echo "------------------------------------------------------------ STEP 14: CREATING DADA2 QC SUMMARY ------------------------------------------------------------"
# Create a DADA2 QC Summary (Additional QC Check)
qc_file_list=$(cat $my_output/${run_name}_dada2_QC_stats.txt | awk '{print $1}' | grep -Ev "Sample|Undetermined|Total_reads")
json_file_list=$(cat $my_output/${run_name}_json_summary.txt | awk '{print $1}' | tail -n +16 | tr -d : | grep -v Undetermined)

qc_file_num=$(cat $my_output/${run_name}_dada2_QC_stats.txt | awk '{print $1}' | grep -Ev "Sample|Undetermined|Total_reads" | wc -l)
json_file_num=$(cat $my_output/${run_name}_json_summary.txt | awk '{print $1}' | tail -n +16 | tr -d : | grep -v Undetermined | wc -l)

dada2_QC_summary=$my_output/End_user_results/${run_name}_dada2_QC_summary.txt

echo ------------DADA2_QC_Summary Table------------ >> $dada2_QC_summary
echo Total number of submitted samples: $qc_file_num >> $dada2_QC_summary
echo Total number of successfully processed samples: $json_file_num >> $dada2_QC_summary
echo ---------------------- >> $dada2_QC_summary
echo Total number of failed samples: `expr $qc_file_num - $json_file_num` >> $dada2_QC_summary
echo Failed samples listed below: >> $dada2_QC_summary
echo "$qc_file_list" | grep -v "$json_file_list" >> $dada2_QC_summary
echo ---------------------- >> $dada2_QC_summary
echo Samples with less than 1000 reads. >> $dada2_QC_summary
echo "Sample_ID Sample_Reads(Nonchimera_reads)" >> $dada2_QC_summary
awk '{if(($6 < 1000) && ($3 != 0)) {print $1, $6}}' $my_output/${run_name}_dada2_QC_stats.txt >> $dada2_QC_summary
echo ---------------------- >> $dada2_QC_summary
echo See $my_output/${run_name}_dada2_QC_stats.txt for more details. >> $dada2_QC_summary


# Script Timer
ENDTIMER="$(date +%s)"
DURATION=$[${ENDTIMER} - ${STARTTIMER}]
HOURS=$((${DURATION} / 3600)) 
MINUTES=$(((${DURATION} % 3600)/ 60))
SECONDS=$(((${DURATION} % 3600) % 60))

echo "RUNTIMER: $HOURS:$MINUTES:$SECONDS (hh:mm:ss)" >> $dada2_QC_summary



echo "------------------------------------------------------------ STEP 15: CLEANING UP DIRECTORY ------------------------------------------------------------"
# Directory clean up for Box transfer and End User Folder separation
mkdir Supporting-Files
mv $(ls | grep -v '\.err$\|\.out$') Supporting-Files

mkdir End-User-Results
cp Supporting-Files/$my_output/End_user_results/* End-User-Results


echo "------------------------------------------------------------ DADA2 PIPELINE COMPLETE, GOOD BYE ------------------------------------------------------------"
