# 16S_Sequencing_QC
# Original Script: 20220629_16S_Sequencing_QC_Template_v4.Rmd by Tony and Neelesh
# Modified and updated by: Edward Sung (edward.sung@gwu.edu)
#
# Current Version: v6.0.3
# Version Date: 8/08/24
#
# 16S Sequencing QC script is used to process 16S dada2 pipeline generated output files.



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



# Script Table of Contents
# 1. Libaries
# 2. Functions
# 3. Load in data
# 4. Sample Data Filtering / Selecting
# 5. Data QC, Cleaning and Subsetting
# 6. ASV Searching
# 7. Export Data



# 1. Libaries ----
# Load Packages
library("Biostrings")
library("stringr")
library("dplyr")
library("tidyr")
library("tidyverse")
library("ggplot2")
library("ggthemes")
library("ggpubr")
library("openxlsx")


# 2. Functions ----

# Global function variables
classLevel_ID <- list("phylum_index" = 2,
                      "class_index" = 3,
                      "order_index" = 4,
                      "family_index" = 5,
                      "genus_index" = 6,
                      "species_index" = 7)


# Converts ASV_#_ClassLevel user input into index version
asv_classlevel_to_index <- function(ASV_ClassLevel) {
  ASV_CLassLevel_Index <- list(
    Phylum = "phylum_index",
    Class = "class_index",
    Order = "order_index",
    Family = "family_index",
    Genus = "genus_index"
  )
  return(ASV_CLassLevel_Index[[ASV_ClassLevel]])
}


# Helper function that extracts the lowest target classification ID
class_ID_helper <- function(taxonomy_id, classlevel_index) {
  lineName <- unlist(strsplit(taxonomy_id, split = ";"))[1:classlevel_index]
  final_id <- list()

  for (i in 0:(length(lineName) - 1)) {
    pos_pointer <- length(lineName) - i
    if (lineName[pos_pointer] == "Unclassified") {
      final_id <- append(lineName[pos_pointer], final_id)
    } else {
      final_id <- append(lineName[pos_pointer], final_id)
      break
    }
  }
  return(paste(final_id, collapse = "_"))
}


# Function transposes the dataframe into standard viewing/processing format
tranpose_func <- function(df_input) {
  df_t <- t(df_input)
  colnames(df_t) <- df_t[1,]
  df_t <- df_t[-1,] %>%
    as.data.frame() %>%
    rownames_to_column("Sample_ID")

  # Ensures column id is characters and rest of the columns as numeric
  df_t[, 1] <- as.character(df_t[, 1])
  df_t[, -1] <- lapply(df_t[, -1], function(x) as.numeric(as.character(x)))

  return(df_t)
}


# Helper function that removes all zero taxas and obtains lowest target classification ID for data processing
clean_data_helper <- function(df_input, classlevel_index) {
  df <- df_input %>%
    dplyr::rename("taxonomy" = "X.OTU.ID") %>%
    select(-contains("Undetermined")) %>% # Remove Undetermined
    mutate(taxonomy = gsub("<", "_", taxonomy)) %>%
    mutate(taxonomy = gsub(">=", "_G", taxonomy)) %>%
    mutate(zerosum_check = rowSums(select(., -taxonomy))) %>%
    filter(zerosum_check > 0) %>%  # Select all non-zero taxa rows
    select(-zerosum_check) %>%
    rowwise() %>%
    mutate(class_ID = class_ID_helper(taxonomy, classLevel_ID[[classlevel_index]])) %>%
    relocate(class_ID, .before = taxonomy) %>%
    select(-taxonomy) %>%
    arrange(class_ID) %>%
    as.data.frame()

  return(df)
}


# Requires the output df from clean_data_helper function as the input_df
QC_NEC_PEC_func <- function(df_input) {

  # NEC Table
  nec_df <- df_input %>%
    select(starts_with(c("class_ID", "nec", "ntc"))) %>%
    mutate(zerosum_check = rowSums(select(., -class_ID))) %>%
    filter(zerosum_check > 0) %>% # Removes taxas that do not have reads
    select(-zerosum_check) %>%
    tranpose_func() %>%
    mutate(Total_Reads = rowSums(select(., -Sample_ID))) %>%
    relocate(Total_Reads, .after = Sample_ID)

  # PEC Table
  pec_df <- df_input %>%
    select(starts_with(c("class_ID", "pec", "ptc"))) %>%
    mutate(zerosum_check = rowSums(select(., -class_ID))) %>%
    filter(zerosum_check > 0) %>% # Removes taxas that do not have reads
    select(-zerosum_check) %>%
    tranpose_func() %>%
    mutate(Total_Reads = rowSums(select(., -Sample_ID)))

  pec_tmp <- (select(pec_df, -c(Sample_ID, Total_Reads)) / pec_df$Total_Reads) %>%
    round(digits = 2) %>%
    rename_with(~paste("prop", ., sep = "_"))

  pec_df <- pec_df %>%
    cbind(pec_tmp)

  return(list(nec_df, pec_df))
}


# Requires the output df from clean_data_helper function as the input_df
# Generates the reseq_data as well as finialized df table by removing the NEC and PEC and reseq samples
cleaned_data_func <- function(df_input, reseq_cutoff_value) {

  # reseq Table
  reseq_df <- df_input %>%
    select(-starts_with(c("nec", "ntc", "pec", "ptc"))) %>%
    tranpose_func() %>%
    mutate(Total_Reads = rowSums(select(., -Sample_ID))) %>%
    relocate(Total_Reads, .after = Sample_ID) %>%
    filter(Total_Reads < reseq_cutoff_value)

  # finalized Table
  final_df <- df_input %>%
    select(-starts_with(c("nec", "ntc", "pec", "ptc", reseq_df$Sample_ID))) %>%
    tranpose_func() %>%
    mutate(Total_Reads = rowSums(select(., -Sample_ID))) %>%
    relocate(Total_Reads, .after = Sample_ID)

  return(list(reseq_df, final_df))
}


# Takes in finalized Genus level data and marks samples that are contaminated with those in that list above threshold
contamination_check_func <- function(df_input) {
  contamination_list <- c("Cellulomonas",
                          "Chryseobacterium",
                          "Sphingomonas",
                          "Delftia",
                          "Pseudoxanthomonas")

  df <- df_input %>%
    select(Sample_ID, Total_Reads, any_of(contamination_list)) %>%
    # Sums up all reads in contamination list, divides by Total_Reads for that sample, then converts to 1 or 0, if greater than or equal to 0.25
    mutate(contaminated = ifelse((rowSums(select(., -c(Sample_ID, Total_Reads))) / Total_Reads) >= 0.25, 1, 0)) %>%
    select(Sample_ID, contaminated)

  return(df)
}


# Generates summary data based on the different class levels
grouped_class_summary_func <- function(df_input, classlevel_index, reseq_df) {
  output_df <- clean_data_helper(df_input, classLevel_ID[[classlevel_index]]) %>%
    select(-starts_with(c("nec", "ntc", "pec", "ptc", reseq_df$Sample_ID))) %>%
    group_by(class_ID) %>%
    summarize(across(everything(), sum)) %>%
    tranpose_func() %>%
    mutate(Total_Reads = rowSums(select(., -Sample_ID))) %>%
    relocate(Total_Reads, .after = Sample_ID)

  return(output_df)
}


# Function to search for ASV of interests from the ASV data
asv_of_interest_func <- function(df_asv_data_input, df_asv_seq_input, reseq_df, asv_search, classlevel_index) {
  df <- df_asv_data_input %>%
    dplyr::rename("ASV_ID" = "X.OTU.ID") %>%
    left_join(df_asv_seq_input, by = "ASV_ID") %>%
    select(-starts_with(c("nec", "ntc", "pec", "ptc", "undetermined", reseq_df$Sample_ID))) %>%
    mutate(taxonomy = gsub("<", "_", taxonomy)) %>%
    mutate(taxonomy = gsub(">=", "_G", taxonomy)) %>%
    mutate(Total_Seq_Reads = rowSums(select(., -c(taxonomy, Sequence, ASV_ID)))) %>%
    filter(Total_Seq_Reads > 0) %>%  # Select all non-zero taxa rows
    rowwise() %>%
    mutate(class_ID = class_ID_helper(taxonomy, classLevel_ID[[classlevel_index]])) %>%
    mutate(class_ID = gsub(" ", "", class_ID)) %>%
    filter(str_detect(class_ID, regex(paste0("^", asv_search, "$"), ignore_case = TRUE))) %>% # \\b needed to specify regex matchcase to be specific for end
    as.data.frame() %>%
    mutate(Prop_ASV_Total_Seq_Reads = Total_Seq_Reads / sum(Total_Seq_Reads)) %>%
    select(ASV_ID, class_ID, Total_Seq_Reads, Prop_ASV_Total_Seq_Reads, Sequence) %>%
    arrange(desc(Total_Seq_Reads))

  return(df)
}


# Write to xlsx documents / workbooks
xlsxWrite <- function(wb, sheet, data) {
  addWorksheet(wb, sheet)

  writeData(wb,
            sheet,
            data,
            withFilter = TRUE)

  setColWidths(wb, sheet, 1:ncol(data), widths = "Auto")
}


# 3. Load in data ----
# Set working directory
setwd(work_dir)

# Load in the accepted_16S_species text file (Currently textfile needs to be in the working directory)
accepted_16S_Species <- c("Sample_ID", "contaminated", "Total_Reads") # Includes these columns when filtering for selected species list.
accepted_16S_Species <- append(accepted_16S_Species, readLines("16S_Sequencing_QC_16SAcceptedSpecies.txt"))

# Raw data
rawGenus_data <- read.delim(rawGenus_filename, sep = "\t", skip = 1, header = TRUE)
rawSpecies_data <- read.delim(rawSpecies_filename, sep = "\t", skip = 1, header = TRUE)
rawASV_data <- read.delim(rawASV_filename, sep = "\t", skip = 1, header = TRUE)

# Raw ASV Data
rawASV_seq_data <- readDNAStringSet(rawASV_Seq_filename)
ASV_df <- data.frame(names(rawASV_seq_data), paste(rawASV_seq_data))
colnames(ASV_df) <- c("ASV_ID", "Sequence")

# Prints the Total_samples
# Use this number to compare to total expected number of samples
# This is used to catch any low/zero read samples that wouldn't show up in the QC - Resequence tab.
# Example would be dada2 pipeline removed all reads during the trimming process.

Total_samples = rawGenus_data %>%
  select(-contains(c("X.OTU.ID", "Undetermined"))) %>%
  ncol()
print(paste("Total Number of Samples:", Total_samples))


# 4. Sample Data Filtering / Selecting ----
# Columns to keep when using the regex sample name filtering
if (filter_samplelist_filename != "") {
  CoreColumns <- c("X.OTU.ID", "taxonomy", "nec", "ntc", "pec", "ptc")

  filter_samplelist_raw = read.csv(filter_samplelist_filename, header = FALSE)
  filter_samplelist = filter_samplelist_raw[, 1]

  rawGenus_data_filter = rawGenus_data %>%
    select(starts_with(c(CoreColumns, filter_samplelist)))
  rawSpecies_data_filter <- rawSpecies_data %>%
    select(starts_with(c(CoreColumns, filter_samplelist)))
  rawASV_data_filter <- rawASV_data %>%
    select(starts_with(c(CoreColumns, filter_samplelist)))
}


# First Regex Filter
if (filter_regex_1 != "") {
  CoreColumns <- paste("X.OTU.ID", "taxonomy", "nec", "ntc", "pec", "ptc", sep = "|")

  if (filter_regex_negate_1 == "yes") {
    rawGenus_data_filter <- rawGenus_data %>%
      select(!matches(c(paste(filter_regex_1, sep = "|"))))
    rawSpecies_data_filter <- rawSpecies_data %>%
      select(!matches(c(paste(filter_regex_1, sep = "|"))))
    rawASV_data_filter <- rawASV_data %>%
      select(!matches(c(paste(filter_regex_1, sep = "|"))))
  }
  else {
    rawGenus_data_filter <- rawGenus_data %>%
      select(matches(c(paste(CoreColumns, filter_regex_1, sep = "|"))))
    rawSpecies_data_filter <- rawSpecies_data %>%
      select(matches(c(paste(CoreColumns, filter_regex_1, sep = "|"))))
    rawASV_data_filter <- rawASV_data %>%
      select(matches(c(paste(CoreColumns, filter_regex_1, sep = "|"))))
  }
}

# Second Regex Filter 2 only works if Regex Filter 1 is used
if (filter_regex_1 != "" && filter_regex_2 != "") {
  if (filter_regex_negate_2 == "yes") {
    rawGenus_data_filter <- rawGenus_data_filter %>%
      select(!matches(c(paste(filter_regex_2, sep = "|"))))
    rawSpecies_data_filter <- rawSpecies_data_filter %>%
      select(!matches(c(paste(filter_regex_2, sep = "|"))))
    rawASV_data_filter <- rawASV_data_filter %>%
      select(!matches(c(paste(filter_regex_2, sep = "|"))))
  }
  else {
    rawGenus_data_filter <- rawGenus_data_filter %>%
      select(matches(c(paste(CoreColumns, filter_regex_2, sep = "|"))))
    rawSpecies_data_filter <- rawSpecies_data_filter %>%
      select(matches(c(paste(CoreColumns, filter_regex_2, sep = "|"))))
    rawASV_data_filter <- rawASV_data_filter %>%
      select(matches(c(paste(CoreColumns, filter_regex_2, sep = "|"))))
  }
}

# No sample filtering applied
if (filter_regex_1 == "" && filter_regex_2 == "" && filter_samplelist_filename == "") {
  rawGenus_data_filter <- rawGenus_data
  rawSpecies_data_filter <- rawSpecies_data
  rawASV_data_filter <- rawASV_data
}


# 5. Data QC, Cleaning and Subsetting ----
Genus_data <- clean_data_helper(rawGenus_data_filter, "genus_index")
Species_data <- clean_data_helper(rawSpecies_data_filter, "species_index")

# QC Data
Genus_data_nec_pec <- QC_NEC_PEC_func(Genus_data)
Genus_data_nec <- Genus_data_nec_pec[[1]]
Genus_data_pec <- Genus_data_nec_pec[[2]]

# Cleaned Data
Genus_data_reseq_mainData <- cleaned_data_func(Genus_data, reseq_cutoff_value)
Genus_data_reseq <- Genus_data_reseq_mainData[[1]]
Genus_data_mainData <- Genus_data_reseq_mainData[[2]]
Species_data_reseq_mainData <- cleaned_data_func(Species_data, reseq_cutoff_value)
Species_data_mainData <- Species_data_reseq_mainData[[2]]

# Contamination Check
Genus_data_contaminate <- contamination_check_func(Genus_data_mainData)

# Subsetted and Further Processing Data - Genus
Genus_data_mainData_reads <- Genus_data_mainData %>%
  left_join(Genus_data_contaminate, by = "Sample_ID") %>%
  relocate(contaminated, .after = Sample_ID)

Genus_data_mainData_presence <- Genus_data_mainData %>%
  select(-Total_Reads) %>%
  mutate(across(2:ncol(.), ~replace(., . > 0, 1))) %>%
  left_join(Genus_data_contaminate, by = "Sample_ID") %>%
  relocate(contaminated, .after = Sample_ID)

Genus_data_mainData_proportion <- Genus_data_mainData %>%
  rowwise() %>%
  mutate(across(3:ncol(.), ~replace(., . != 0, round(. / Total_Reads, 4)))) %>%
  left_join(Genus_data_contaminate, by = "Sample_ID") %>%
  relocate(contaminated, .after = Sample_ID)

# Subsetted and Further Processing Data - Species
Species_data_mainData_reads <- Species_data_mainData %>%
  left_join(Genus_data_contaminate, by = "Sample_ID") %>%
  relocate(contaminated, .after = Sample_ID) %>%
{ if (species_data_filter == "yes") select(., matches(paste0("^", paste(accepted_16S_Species, collapse = "|")))) else . }

Species_data_mainData_presence <- Species_data_mainData %>%
  select(-Total_Reads) %>%
  mutate(across(2:ncol(.), ~replace(., . > 0, 1))) %>%
  left_join(Genus_data_contaminate, by = "Sample_ID") %>%
  relocate(contaminated, .after = Sample_ID) %>%
{ if (species_data_filter == "yes") select(., matches(paste0("^", paste(accepted_16S_Species, collapse = "|")))) else . }

Species_data_mainData_proportion <- Species_data_mainData %>%
  rowwise() %>%
  mutate(across(3:ncol(.), ~replace(., . != 0, round(. / Total_Reads, 4)))) %>%
  left_join(Genus_data_contaminate, by = "Sample_ID") %>%
  relocate(contaminated, .after = Sample_ID) %>%
{ if (species_data_filter == "yes") select(., matches(paste0("^", paste(accepted_16S_Species, collapse = "|")))) else . }

# Subsetted and Further Processing Data - Phylum
grouped_phylum_data <<- grouped_class_summary_func(rawGenus_data_filter, "phylum_index", Genus_data_reseq) %>%
  left_join(Genus_data_contaminate, by = "Sample_ID") %>%
  relocate(contaminated, .before = Total_Reads)

# Subsetted and Further Processing Data - Class
grouped_class_data <<- grouped_class_summary_func(rawGenus_data_filter, "class_index", Genus_data_reseq) %>%
  left_join(Genus_data_contaminate, by = "Sample_ID") %>%
  relocate(contaminated, .before = Total_Reads)

# Subsetted and Further Processing Data - Order
grouped_order_data <<- grouped_class_summary_func(rawGenus_data_filter, "order_index", Genus_data_reseq) %>%
  left_join(Genus_data_contaminate, by = "Sample_ID") %>%
  relocate(contaminated, .before = Total_Reads)

# Subsetted and Further Processing Data - Family
grouped_family_data <<- grouped_class_summary_func(rawGenus_data_filter, "family_index", Genus_data_reseq) %>%
  left_join(Genus_data_contaminate, by = "Sample_ID") %>%
  relocate(contaminated, .before = Total_Reads)

# View QC Dataframes
view(Genus_data_nec)
view(Genus_data_pec)
view(Genus_data_reseq)


# 6. ASV Searching ----
# Initalize asv variables
asv_inputs <- c()
asv_table <- c()
asvList <- c()

# ASV_1
if (ASV_1 != "NA") {
  asvList <- c(asvList, "asv1")
  asv_inputs$asv1 <- ASV_1
  asv_table$asv1 <- asv_of_interest_func(rawASV_data_filter, ASV_df, Genus_data_reseq, ASV_1, asv_classlevel_to_index(ASV_1_ClassLevel))
}

# ASV_2
if (ASV_2 != "NA") {
  asvList <- c(asvList, "asv2")
  asv_inputs$asv2 <- ASV_2
  asv_table$asv2 <- asv_of_interest_func(rawASV_data_filter, ASV_df, Genus_data_reseq, ASV_2, asv_classlevel_to_index(ASV_2_ClassLevel))
}

# ASV_3
if (ASV_3 != "NA") {
  asvList <- c(asvList, "asv3")
  asv_inputs$asv3 <- ASV_3
  asv_table$asv3 <- asv_of_interest_func(rawASV_data_filter, ASV_df, Genus_data_reseq, ASV_3, asv_classlevel_to_index(ASV_3_ClassLevel))
}

# ASV_4
if (ASV_4 != "NA") {
  asvList <- c(asvList, "asv4")
  asv_inputs$asv4 <- ASV_4
  asv_table$asv4 <- asv_of_interest_func(rawASV_data_filter, ASV_df, Genus_data_reseq, ASV_4, asv_classlevel_to_index(ASV_4_ClassLevel))
}

# ASV_5
if (ASV_5 != "NA") {
  asvList <- c(asvList, "asv5")
  asv_inputs$asv5 <- ASV_5
  asv_table$asv5 <- asv_of_interest_func(rawASV_data_filter, ASV_df, Genus_data_reseq, ASV_5, asv_classlevel_to_index(ASV_5_ClassLevel))
}

# ASV_6
if (ASV_6 != "NA") {
  asvList <- c(asvList, "asv6")
  asv_inputs$asv6 <- ASV_6
  asv_table$asv6 <- asv_of_interest_func(rawASV_data_filter, ASV_df, Genus_data_reseq, ASV_6, asv_classlevel_to_index(ASV_6_ClassLevel))
}

# ASV_7
if (ASV_7 != "NA") {
  asvList <- c(asvList, "asv7")
  asv_inputs$asv7 <- ASV_7
  asv_table$asv7 <- asv_of_interest_func(rawASV_data_filter, ASV_df, Genus_data_reseq, ASV_7, asv_classlevel_to_index(ASV_7_ClassLevel))
}

# ASV_8
if (ASV_8 != "NA") {
  asvList <- c(asvList, "asv8")
  asv_inputs$asv8 <- ASV_8
  asv_table$asv8 <- asv_of_interest_func(rawASV_data_filter, ASV_df, Genus_data_reseq, ASV_8, asv_classlevel_to_index(ASV_8_ClassLevel))
}

# ASV_9
if (ASV_9 != "NA") {
  asvList <- c(asvList, "asv9")
  asv_inputs$asv9 <- ASV_9
  asv_table$asv9 <- asv_of_interest_func(rawASV_data_filter, ASV_df, Genus_data_reseq, ASV_9, asv_classlevel_to_index(ASV_9_ClassLevel))
}

# ASV_10
if (ASV_10 != "NA") {
  asvList <- c(asvList, "asv10")
  asv_inputs$asv10 <- ASV_10
  asv_table$asv10 <- asv_of_interest_func(rawASV_data_filter, ASV_df, Genus_data_reseq, ASV_10, asv_classlevel_to_index(ASV_10_ClassLevel))
}


# 7. Export Data ----
# mainData_wb contain sheets: 16S_Reads_Genus, 16S_Presence_Genus, 16S_Proportional_Abundance_Genus
mainData_wb <- createWorkbook()

# QC_wb contain sheets: NEC, PEC, Samples_To_Resequence
QC_wb <- createWorkbook()

# ASV_wb contain sheets: ASV_Searches
ASV_wb <- createWorkbook()

# Groups_wb contain sheets: Phylum, Class, Order, Family
Groups_wb <- createWorkbook()

# Genus Data mainData_wb
xlsxWrite(mainData_wb, "Genus_Reads", Genus_data_mainData_reads)
xlsxWrite(mainData_wb, "Genus_Presence", Genus_data_mainData_presence)
xlsxWrite(mainData_wb, "Genus_Proportional_Abundance", Genus_data_mainData_proportion)

# Species Data mainData_wb
xlsxWrite(mainData_wb, "Species_Reads", Species_data_mainData_reads)
xlsxWrite(mainData_wb, "Species_Presence", Species_data_mainData_presence)
xlsxWrite(mainData_wb, "Species_Proportional_Abundance", Species_data_mainData_proportion)

# QC_wb
xlsxWrite(QC_wb, "NEC", Genus_data_nec)
xlsxWrite(QC_wb, "PEC", Genus_data_pec)
xlsxWrite(QC_wb, "Samples_To_Resequence", Genus_data_reseq)

# ASV_wb
asv_save <- 0
if (length(asvList) != 0) {
  for (asv_num in asvList) {
    xlsxWrite(ASV_wb, paste("ASV", asv_inputs[[asv_num]], sep = "_"), asv_table[[asv_num]])
  }
  asv_save <- 1
}

# Groups_wb
xlsxWrite(Groups_wb, "Phylum,", grouped_phylum_data)
xlsxWrite(Groups_wb, "Class", grouped_class_data)
xlsxWrite(Groups_wb, "Order", grouped_order_data)
xlsxWrite(Groups_wb, "Family", grouped_family_data)

outputFiles <- c("mainData.xlsx", "QC.xlsx", "ASV.xlsx", "Groups.xlsx")
for (file_index in 1:length(outputFiles)) {
  outputFiles[file_index] <- paste(download_data_filename, outputFiles[file_index], sep = "_")
}

# Save xlsx to working directory
saveWorkbook(mainData_wb, outputFiles[1], overwrite = TRUE)
saveWorkbook(QC_wb, outputFiles[2], overwrite = TRUE)
if (asv_save == 1) {
  saveWorkbook(ASV_wb, outputFiles[3], overwrite = TRUE)
}
saveWorkbook(Groups_wb, outputFiles[4], overwrite = TRUE)

# Save a zip version of the files
zip(download_data_filename, outputFiles)
