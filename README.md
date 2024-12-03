# mb_analysis
Code repository used to process 16S sequencing data used in: Composition and Dynamics of the Adult Nasal Microbiome.


### Inhouse DADA2 pipeline requirements:
| Tool           | Version   		   |
| -------------- | ----------------------- |
| SLURM          | v23.02.4  		   | 
| Trimmomatic    | v0.39     		   |
| Cutadapt       | v2.10     		   |
| GNU parallel   | v20201122 		   |
| R	         | >= v3.5.0    	   |
| R-Bioconductor | >= v3.12   		   |
| DADA2 	 | v1.16  		   |
| BIOM 		 | v2.1.10 		   |
| QIIME 	 | v1.9.1 		   |
| RDP-Classifier | v2.12 		   |
| perl 		 | >= v5.26.2 		   |
| java 		 | build 10.0.2 2018-07-17 |

Note!!: Samplename (the first part of the file name, for example: B123-GWU-DNA) should avoid containing underscores (this applies to DADA2 and all Illumina systems). Use of dashes is recommended.

Please use the provided conda yaml files in `conda_enviroment_yml` folder to build the associated conda enviroments used in **dada2_pipeline_vP_v11_cleaned.sh**.  
Please update the path locations of the pipeline files and all dependency scripts to your own file paths:  
&nbsp;&nbsp;&nbsp;&nbsp; **scripts** = 'YOUR_PATH/dada2_pipeline_scripts'  
&nbsp;&nbsp;&nbsp;&nbsp; **rdp** = 'YOUR_PATH/rdp_classifier_2.12/dist'  
&nbsp;&nbsp;&nbsp;&nbsp; **v6species** = 'YOUR_PATH/V6_species'  
&nbsp;&nbsp;&nbsp;&nbsp; **database** = 'YOUR_PATH/dada2_pipeline_databases'  

General pipeline system call for SLURM environment (use qsub instead of sbatch for Sun Grid Engine)

---

## Script Usage

### dada2_pipeline_vP_v11_cleaned.sh
```
dada2_pipeline_vP_v11_cleaned.sh --h

USAGE: options available
                
                mandatory options:

                                -i or --input : provide the name of your input folder (containing FASTQ files)
                                
                optional options:
                                -p or --primer : path to fasta file with primer sequences (default:YOUR_PATH/dada2_pipeline_databases/amplicon_primers_bactquant.fasta) 
                                -r or --run : provide the name of your sequencing run e.g.:run17
                                -t or --threads : give the number of threads to use (default is 16)
                                -dt or --taxonomydb : give the name of your taxonomy database (rdp_train_set_16.fa.gz)
                                -ds or --speciesdb : give the name of your species database (rdp_species_assignment_16.fa.gz)
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
```
Note!!: When migrating this pipeline to environments with a Sun Grid Engine scheduler, make sure to alter the first five lines of the workflow script to make it compatible!!  
Two folders will be generated called **End-User-Results** and **Supporting-Files**. The End-User-Results folder will contain files that are fed into 16S_Sequencing_QC_v6.0.3.R. The Supporting-Files contains everything in the End-User-Results folder, as well as the intermediate product files.

### 16S_Sequencing_QC_v6.0.3.R
Please ensure the **16S_Sequencing_QC_16SAcceptedSpecies.txt** from the `database_files` folder is in the same working directory as **16S_Sequencing_QC_v6.0.3.R**. The 16S_Sequencing_QC_16SAcceptedSpecies.txt contains species that were specifically trained to be accruately classified by the Naïve Bayesian Classifier.  
After running the **dada2_pipeline_vP_v11_cleaned.sh** pipeline, several dada2 output files will be generated which are fed into the **16S_Sequencing_QC_v6.0.3.R** to generate read counts and ASVs.  
```
# -------- START OF USER INPUTS --------
# Set Working Directory
work_dir <- "Type_in_full_working_directory"                       # Example: /Users/edward.sung/Desktop/MyWorkingDirectory


# File Uploads
rawGenus_filename <- "Type_in_RDP212_GENUS_c80_tax_L6_txt"         # Example: "LiuPriceLab_Run53_GWU_RDP212_GENUS_c80_tax_L6.txt"
rawSpecies_filename <- "Type_in_RDP212_V6_SPECIES_c80_tax_L7_txt"  # Example: "LiuPriceLab_Run53_GWU_RDP212_V6-SPECIES_c80_tax_L7.txt"
rawASV_filename <- "Type_in_RDP212_GENUS_c80_tax_tsv"              # Example: "LiuPriceLab_Run53_GWU_RDP212_GENUS_c80_tax.tsv"
rawASV_Seq_filename <- "Type_in_species_fasta"                     # Example: "LiuPriceLab_Run53_GWU_species.fasta"


# Name output files (Do not include .xlsx extension in name)
download_data_filename <- "Type_in_output_file_name"               # Example: LiuPriceLab_Run53_20240513


# Resequencing Cutoff
reseq_cutoff_value <- 1000                                         # Default: 1000


# Sample filtering
# For sample filtering, please choose either filtering by samplelist or by regex, but not both.

# filter_samplelist_filename is a csv file with each row being a sample name that you want to keep.
# Make sure to convert any non-alphanumerical characters to .
# For example, any "-" dashes are converted to "." periods.
filter_samplelist_filename <- ""

# Please do not include in your regex filter: nec, ntc, pec, ptc, X.OTU.ID, taxonomy. It will remove key columns necessary, pending fix in future.
# Performs 1st sample filtering
filter_regex_negate_1 <- "yes"                                      # 'yes' - Apply filter as negate (Uses ! negation) or 'no' - Apply filters as is and keeps those that matches regex.
filter_regex_1 <- ""                                                # Uses the matches() helper function of dyplr for select function; Example: "Pre|MN|IN|blank"

# Performs 2nd sample filtering after 1st filter
filter_regex_negate_2 <- "yes"                                      # 'yes' - Apply filter as negate (Uses ! negation) or 'no' - Apply filters as is and keeps those that matches regex.
filter_regex_2 <- ""                                                # Uses the matches() helper function of dyplr for select function; Example: "Pre|MN|IN|blank"


# Species Data Filtering - Filters the MainData species tabs for only confidently classified species by RDP (rdp_classifier_version-speciesSummary10MAY2018)
species_data_filter <- "yes"                                        # Choices are `yes` or `no`


# Set and Search for ASV of interests
# Please captialize the first letter
# Options for ASV_#_ClassLevel = 'Phylum', 'Class', 'Order', 'Family', 'Genus'

ASV_1 <- "NA"                                                       # Example: Streptococcus
ASV_1_ClassLevel = "NA"                                             # Example: Genus

ASV_2 <- "NA"                                                       # Example: Staphylococcus
ASV_2_ClassLevel = "NA"                                             # Example: Genus

ASV_3 <- "NA"                                                       # Example: Corynebacterium
ASV_3_ClassLevel = "NA"                                             # Example: Genus

ASV_4 <- "NA"                                                       # Example: Moraxella
ASV_4_ClassLevel = "NA"                                             # Example: Genus

ASV_5 <- "NA"                                                       # Example: Propionibacterium
ASV_5_ClassLevel = "NA"                                             # Example: Genus

ASV_6 <- "NA"                                                       # Example: Dolosigranulum
ASV_6_ClassLevel = "NA"                                             # Example: Genus

ASV_7 <- "NA"
ASV_7_ClassLevel = "NA"

ASV_8 <- "NA"
ASV_8_ClassLevel = "NA"

ASV_9 <- "NA"
ASV_9_ClassLevel = "NA"

ASV_10 <- "NA"
ASV_10_ClassLevel = "NA"

# -------- END OF USER INPUTS --------
```

### (Name of CST Classifier Algo text file)
Small descriptions and code block of algo.


---

## Dada2 Pipeline Processing Description

### Data transfer to the cluster and check MD5 codes (perform manually):
Copy files to a common location (when you need to transfer compressed FASTQ files (.gz) from the MiSeq to your computer system).
```
cp **/*.gz raw_reads_directory
```
Example of the demultiplexed file naming format:  
B123-GWU-DNA_325_L001_R1_001.fastq.gz (contains suffix with R1 or R2 to designate forward and reverse for pair-end reads).


### Create input folders (performed by pipeline):
Our pipeline creates two directories for storing trimmed and primer-free reads and a folder for storing all the output data.


### Quality trim raw data (performed by pipeline):
We run Trimmomatic to trim and remove low quality sequences from our raw data.
```
perl $scripts/run_trimmomatic_amplicon_mode.pl $input_folder/*

	Selected settings (inside Perl script):
  	--leading 3		Remove leading low quality or 3 bases	
  	--trailing 3		Remove trailing low quality or 3 bases
  	--slidingwindow 4:15	Scan the read with a 4-base wide sliding window, trimming when the average quality per base drops below 15
  	--minlen 225		Drop the read if it is below 225 bp (for V3-V4 regions)
```


### Pre-processing raw data (performed by pipeline):
DADA2 infers error models from reads and will generate inaccurate data if primer sequences are included, so we first need to remove them. We use Cutadapt 2.1 on each FASTQ file pair (forward and reverse) with GNU parallel.
```
cutadapt --pair-filter any --no-indels -e 0.15 --discard-untrimmed -g file:${primer_path} -G file:${primer_path} -Z -o primer_trimmed_reads/${f1}.gz -p primer_trimmed_reads/${f2}.gz ${1} ${2} > primer_trimmed_reads/${f1%%.fastq*}_cutadapt_log.txt

export -f remove_adapter;

parallel --link --jobs $threads remove_adapter ::: trimmed_reads/*_R1_*.fastq ::: trimmed_reads/*_R2_*.fastq
```
file:amplicon_primers_f.fasta or file:amplicon_primers_r.fasta – make sure the file contains only the primers you have used for sequencing. 

Then filter out reads with ambigious base calls (Ns)
```
$scripts/dada2_filter.R -f primer_trimmed_reads --maxN 0 --truncQ 0 --threads $threads --truncLen 0 --maxEE 2,2 --f_match _R1_.*fastq.gz --r_match _R2_.*fastq.gz -o $my_output/${run_name}_filtered_fastqs
```


### DADA2 inference (performed by pipeline):
Merge clean paired-end reads and subsequently infer error models from the reads and identify Amplicon Sequence Variants (ASVs). 
```
$scripts/dada2_inference.R -f $my_output/${run_name}_filtered_fastqs/ --seed 4124 -t $threads --verbose --plot_errors --minOverlap 12 -o $my_output/${run_name}_seqtab.rds
```


### DADA2 chimera checking and taxonomy assignment (performed by pipeline): 
Remove chimeric variants and assign taxonomy to each non-chimeric variant. The dada2 package implements a method to make species level assignments based on exact matching between ASVs and sequenced reference strains. Recent analysis suggests that exact matching (or 100% identity) is the only appropriate way to assign species to 16S gene fragments. 
```
$scripts/dada2_chimera_taxa.R -i $my_output/${run_name}_seqtab.rds -t 2 -r $database/$db_taxonomy -s $database/$db_species --allowMultiple --minBoot 80 --count_out $my_output/${run_name}_seqtab_final.rds --sp_out $my_output/${run_name}_species_final.rds --tax_out $my_output/${run_name}_tax_final.rds
```


### Convert DADA2 R objects to BIOM and FASTA (performed by pipeline):
Convert RDS object to a variant count table (BIOM) and variant sequence file (FASTA).  
```
$scripts/convert_dada2_out.R -i $my_output/${run_name}_seqtab_final.rds -b $my_output/${run_name}_biom.tsv -f $my_output/${run_name}_species.fasta --taxa_in $my_output/${run_name}_species_final.rds --taxa_out $my_output/${run_name}_taxa_metadata.txt
```


### Convert BIOM TSV format to JSON and create summary (performed by pipeline):
JSON is an efficient and commonly used data format for storing bacterial count data. HDF5 is another format, however ColonialOne does not handle it well. 

```
biom convert -i $my_output/${run_name}_biom.tsv -o $my_output/${run_name}_json.biom --to-json

biom summarize-table -i $my_output/${run_name}_json.biom -o $my_output/${run_name}_json_summary.txt
```


### Classify variant sequences in FASTA file using RDP Classifier (performed by pipeline):
Assign taxonomy to all variants (second time) using RDP classification. The pipeline classifies bacterial lineages until the genus level, and performs a species level assignment using a custom classifier, called V6species.
```
# Genus level classification
java -Xmx12g -jar $rdp/classifier.jar classify -o $my_output/RDP212_GENUS_taxa_metadata.txt -q $my_output/${run_name}_species.fasta -f fixrank

# Species level classification using a custom classifier
java -Xmx12g -jar $rdp/classifier.jar classify -o $my_output/RDP212_V6species_taxa_metadata.txt -q $my_output/${run_name}_species.fasta -t $v6species/rRNAClassifier.properties
```


### Add RDP taxonomy labels to BIOM files (performed by pipeline):
Add taxonomy labels to the BIOM files. We add the V6species and default RDP taxonomy labels to separate BIOM tables. 
```
biom add-metadata -i $my_output/${run_name}_json.biom -o $my_output/${run_name}_RDP212_GENUS_c80_tax.biom --observation-metadata-fp $my_output/RDP212_GENUS_taxa_metadata.txt.c80_BIOMversion.tsv --observation-header OTUID,taxonomy --sc-separated taxonomy --output-as-json

biom add-metadata -i $my_output/${run_name}_json.biom -o $my_output/${run_name}_RDP212_V6-SPECIES_c80_tax.biom --observation-metadata-fp $my_output/RDP212_V6species_taxa_metadata.txt.c80_BIOMversion.tsv --observation-header OTUID,taxonomy --sc-separated taxonomy --output-as-json
```


### Summarize taxa and store results in new tables (performed by pipeline):
Provide summary information of the representation of taxonomic groups within each sample. The taxonomic level for which the summary information is provided is designated with -L option (L2 = Phylum, L7 = Species). Please be aware that this step removes the ‘strain’ level information associated with ASVs. For instance, if a sample contains five S. aureus ASVs, the count data for those variants will be merged/summed (in the L7 file). Those five variants will be merged with other Staph species in the L6 file to give one count value per sample etc. 
```
summarize_taxa.py -i $my_output/${run_name}_RDP212_GENUS_c80_tax.biom -L 2,3,4,5,6 -o $my_output/${run_name}_RDP212_GENUS_c80_tax_summarizeTaxa --absolute_abundance

summarize_taxa.py -i $my_output/${run_name}_RDP212_V6-SPECIES_c80_tax.biom -L 2,3,4,5,6,7 -o $my_output/${run_name}_RDP212_V6-SPECIES_c80_tax_summarizeTaxa --absolute_abundance
```


### 16S_Sequencing_QC_v6.0.3.R
After filling in the USER INPUTS section and running the script using the files generated by dada2 pipeline, the script will output several compiled files that can be further processed in R for analysis.
* _Groups.xlsx -- Aggregate the read counts at the taxonomic levels from Phylum to Family.
* _mainData.xlsx -- Primary excel that contains Genus and Species level read counts, along with presence and proportions.
* _QC.xlsx -- Filters out samples identified as NEC and PEC samples for QC checking, along with samples with less than the cutoff (default is 1000).
* _taxaData.csv -- Read counts with ASV seq_ID as rows. Utilized after blasting and identifying ASV seq_id species.
* _ASV.xlsx -- Excel tabs compiles ASV seq_IDs that are labeled as the taxas of interest (specified in User Inputs).
