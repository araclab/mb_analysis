# mb_analysis

**Inhouse DADA2 pipeline requirements:**
1.	SLURM
2.	Trimmomatic 0.38 or 0.39
3.	Cutadapt version 2.4
4.	GNU parallel 
5.	R > 3.5.0 and Bioconductor version > 3.8 
6.	DADA2 version 1.10
7.	BIOM version 2.1.6
8.	QIIME 1.9.1
9.	RDP classifier 2.12
10.	dada2
12.	Perl 


Note!!: Samplename (the first part of the file name, for example: B123-GWU-DNA) should avoid containing underscores (this applies to DADA2 and all Illumina systems). Use of dashes is recommended.

Please update the location of the pipeline and all dependency scripts in dada2_pipeline_vP_v10.sh :  
  **scripts** = 'YOUR_PATH/dada2_pipeline_scripts'  
  **rdp** = 'YOUR_PATH/rdp_classifier_2.12/dist'  
  **v6species** = 'YOUR_PATH/V6_species'  
  **database** = 'YOUR_PATH/dada2_pipeline_databases'  


General pipeline system call for SLURM environment (use qsub instead of sbatch for Sun Grid Engine)

USAGE: options available

sbatch dada2_pipeline_vP_v10.sh -i raw_reads -r runX
                
                mandatory options:

                                -i or --input : provide the name of your input folder (containing FASTQ files)
                                
                optional options:
                                -p or --primer : path to fasta file with primer sequences (default:/GWSPH/groups/liu_price_lab/pegasus_bin/Tools/dada2_pipeline_databases/amplicon_primers_bactquant.fasta) 
                                -r or --run : provide the name of your sequencing run e.g.:run17
                                -t or --threads : give the number of threads to use (default is 2)
                                -dt or --taxonomydb : give the name of your taxonomy database
                                -ds or --speciesdb : give the name of your species database
          Help options:
                                -h or --help : pop up the usage guide

                e.g. :

                                sbatch dada2_pipeline_vP_v10.sh -i raw_reads -r run17 -t 2 -dt rdp_train_set_16.fa.gz -ds rdp_species_assignment_16.fa.gz
                                raw_reads is the name of the input folder
                                run17 is the run name
                                2 is the total number of threads
                                rdp_train_set_16.fa.gz is the name of your taxonomy database
                                rdp_species_assignment_16.fa.gz is the name of your species database

Note!!: When migrating this pipeline to environments with a Sun Grid Engine scheduler, make sure to alter the first five lines of the workflow script to make it compatible!!

**Data transfer to the cluster and check MD5 codes (perform manually):**

Copy files to a common location (when you need to transfer compressed FASTQ files (.gz) from the MiSeq to your computer system)

cp **/*.gz raw_reads_directory

  Example of the demultiplexed file naming format:
  B123-GWU-DNA_325_L001_R1_001.fastq.gz :
  
**Create input folders (performed by pipeline):**

Our pipeline creates two directories for storing trimmed and primer-free reads and a folder for storing all the output data

**Quality trim raw data (performed by pipeline):** 

We run Trimmomatic to trim and remove low quality sequences from our raw data

perl run_trimmomatic_amplicon_mode.pl raw_reads/* 

  Selected settings (inside Perl script):
  --leading 3		Remove leading low quality or 3 bases	
  --trailing 3		Remove trailing low quality or 3 bases
  --slidingwindow 4:15	Scan the read with a 4-base wide sliding window, trimming when the average quality per base drops below 15
  --minlen 225		Drop the read if it is below 225 bp (for V3-V4 regions)

**Pre-processing raw data (performed by pipeline):**

DADA2 infers error models from reads and will generate inaccurate data if primer sequences are included, so we first need to remove them. We use Cutadapt 2.4 on each FASTQ file pair (forward and reverse) with GNU parallel. 
file:amplicon_primers_f.fasta or file:amplicon_primers_r.fasta – make sure the file contains only the primers you have used for sequencing. 

dada2_filter.R -f primer_trimmed_reads --maxN 0 --truncQ 0 --threads $threads --truncLen 0 --maxEE 2,2 --f_match _R1_.*fastq.gz --r_match _R2_.*fastq.gz -o $my_output/${run_name}_filtered_fastqs

**DADA2 inference (performed by pipeline):**

Merge clean paired-end reads and subsequently infer error models from the reads and identify Amplicon Sequence Variants (ASVs). 

dada2_inference.R -f $my_output/${run_name}_filtered_fastqs/ --seed 4124 -t $threads --verbose --plot_errors --minOverlap 12 -o $my_output/${run_name}_seqtab.rds

**DADA2 chimera checking and taxonomy assignment (performed by pipeline):** 

Remove chimeric variants and assign taxonomy to each non-chimeric variant. The dada2 package implements a method to make species level assignments based on exact matching between ASVs and sequenced reference strains. Recent analysis suggests that exact matching (or 100% identity) is the only appropriate way to assign species to 16S gene fragments. 

dada2_chimera_taxa.R -i $my_output/${run_name}_seqtab.rds -t 2 -r $database/$db_taxonomy -s $database/$db_species --allowMultiple --minBoot 80 --count_out $my_output/${run_name}_seqtab_final.rds --sp_out $my_output/${run_name}_species_final.rds --tax_out $my_output/${run_name}_tax_final.rds


**Convert DADA2 R objects to BIOM and FASTA (performed by pipeline):**

Convert RDS object to a variant count table (BIOM) and variant sequence file (FASTA).  

convert_dada2_out.R -i $my_output/${run_name}_seqtab_final.rds -b $my_output/${run_name}_biom.tsv -f $my_output/${run_name}_species.fasta --taxa_in $my_output/${run_name}_species_final.rds --taxa_out $my_output/${run_name}_taxa_metadata.txt

**Convert BIOM TSV format to JSON (performed by pipeline):**

JSON is an efficient and commonly used data format for storing bacterial count data. HDF5 is another format, however ColonialOne does not handle it well. 

biom convert -i $my_output/${run_name}_biom.tsv -o $my_output/${run_name}_json.biom --to-json


**Classify variant sequences in FASTA file using RDP Classifier (performed by pipeline):**

Assign taxonomy to all variants (second time) using RDP classification. The pipeline classifies bacterial lineages until the genus level, and performs a species level assignment using a custom classifier, called V6species.

java -Xmx12g -jar $rdp/classifier.jar classify -o $my_output/RDP212_GENUS_taxa_metadata.txt -q $my_output/${run_name}_species.fasta -f fixrank

java -Xmx12g -jar $rdp/classifier.jar classify -o $my_output/RDP212_V6species_taxa_metadata.txt -q $my_output/${run_name}_species.fasta -t $v6species/rRNAClassifier.properties

**Add RDP taxonomy labels to BIOM files (performed by pipeline):**

Add taxonomy labels to the BIOM files. We add the V6species and default RDP taxonomy labels to separate BIOM tables. 

biom add-metadata -i $my_output/${run_name}_json.biom -o $my_output/${run_name}_RDP212_GENUS_c80_tax.biom --observation-metadata-fp $my_output/RDP212_GENUS_taxa_metadata.txt.c80_BIOMversion.tsv --observation-header OTUID,taxonomy --sc-separated taxonomy --output-as-json

biom add-metadata -i $my_output/${run_name}_json.biom -o $my_output/${run_name}_RDP212_V6-SPECIES_c80_tax.biom --observation-metadata-fp $my_output/RDP212_V6species_taxa_metadata.txt.c80_BIOMversion.tsv --observation-header OTUID,taxonomy --sc-separated taxonomy --output-as-json

BIOM files with taxonomy labels can be converted to plain text (TSV), which is what the pipeline does, and can be used to extract ASV information from a given lineage. 

grep --line-buffered 'Staphylococcus_aureus' run18_test_RDP212_V6-SPECIES_c80_tax.tsv | awk '{print $1}' > saureus_ASVs.txt

perl get_seqs.pl saureus_ASVs.txt    $my_output/${run_name}_species.fasta > saureus.fasta

**Summarize taxa and store results in new tables (performed by pipeline):**

Provide summary information of the representation of taxonomic groups within each sample. The taxonomic level for which the summary information is provided is designated with -L option (L2 = Phylum, L7 = Species). Please be aware that this step removes the ‘strain’ level information associated with ASVs. For instance, if a sample contains five S. aureus ASVs, the count data for those variants will be merged/summed (in the L7 file). Those five variants will be merged with other Staph species in the L6 file to give one count value per sample etc. 

summarize_taxa.py -i $my_output/${run_name}_RDP212_GENUS_c80_tax.biom -L 2,3,4,5,6 -o $my_output/${run_name}_RDP212_GENUS_c80_tax_summarizeTaxa --absolute_abundance

summarize_taxa.py -i $my_output/${run_name}_RDP212_V6-SPECIES_c80_tax.biom -L 2,3,4,5,6,7 -o $my_output/${run_name}_RDP212_V6-SPECIES_c80_tax_summarizeTaxa --absolute_abundance
