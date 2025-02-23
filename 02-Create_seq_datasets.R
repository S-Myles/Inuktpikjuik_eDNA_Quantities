# Packages
library(phyloseq)     # Data structure and functions for seq data
library(tidyverse)    # Data handling and all
library(ape)          # phylogenetic trees
###############################################################################

# Note that depending on the analysis to follow (scripts with "01-" or "02-") prefixes
# You will want to generate phyloseq objects only up to a certain point in this processing file.
# i.e. no abundance filtering for alpha diversity metrics.

# Refer to the paper for description of data for each downstream analyses



###########
# Load 4 markers metabarcoding data and create phyloseq objects
###########

# Project metadata
metadata <-  read_csv("./data/metadata_qPCR_calculated.csv") %>%
  column_to_rownames(var = "sampleID")
metadata <- sample_data(metadata) 


####
# 12S MiFish-U2
####
# ASV table (Sample X Taxon, read counts)
S12_ASV_table <- read_tsv("./data/12S-ASV-table.tsv") %>%
  column_to_rownames('ASV_ID')     
# Taxonomy reference table
S12_taxonomy <- read_tsv("./data/12S-taxonomy.tsv") %>%
  column_to_rownames('ASV_ID') %>%
  as.matrix()
S12_tree <- phy_tree(read.tree("data/12S-filt_aligned_masked_tree_rooted.nwk"))

# Respective phyloseq objects
S12_ASV_table <- otu_table(S12_ASV_table, taxa_are_rows = TRUE)         
S12_taxonomy <- tax_table(S12_taxonomy)
# Merging into 1 global object
(S12_physeq_data <-  merge_phyloseq(S12_ASV_table, metadata, S12_taxonomy, S12_tree)) 

# Cleanup if you want
rm(S12_ASV_table, S12_taxonomy, S12_tree)
##################################


####
# 16S V4V5
####
# ASV table (Sample X Taxon, read counts)
S16_ASV_table <- read_tsv("./data/16S-ASV-table.tsv") %>%
  column_to_rownames('ASV_ID')      
# Taxonomy reference table
S16_taxonomy <- read_tsv("./data/16S-taxonomy.tsv") %>%
  column_to_rownames('ASV_ID') %>%
  as.matrix()
S16_tree <- phy_tree(read.tree("data/16S-filt_aligned_masked_tree_rooted.nwk"))

# Respective phyloseq objects
S16_ASV_table <- otu_table(S16_ASV_table, taxa_are_rows = TRUE)         
S16_taxonomy <- tax_table(S16_taxonomy)
# Merging into 1 global object
(S16_physeq_data <-  merge_phyloseq(S16_ASV_table, metadata, S16_taxonomy,S16_tree))

# Cleanup if you want
rm(S16_ASV_table, S16_taxonomy,S16_tree)
##################################


####
# 18S V4
####
# ASV table (Sample X Taxon, read counts)
S18_ASV_table <- read_tsv("./data/18S-ASV-table.tsv") %>%
  column_to_rownames('ASV_ID')      
# Taxonomy reference table
S18_taxonomy <- read_tsv("./data/18S-taxonomy.tsv") %>%
  column_to_rownames('ASV_ID') %>%
  as.matrix()
S18_tree <- phy_tree(read.tree("data/18S-filt_aligned_masked_tree_rooted.nwk"))


# Respective phyloseq objects
S18_ASV_table <- otu_table(S18_ASV_table, taxa_are_rows = TRUE)         
S18_taxonomy <- tax_table(S18_taxonomy)
# Merging into 1 global object
(S18_physeq_data <-  merge_phyloseq(S18_ASV_table, metadata, S18_taxonomy, S18_tree))

# Cleanup if you want
rm(S18_ASV_table, S18_taxonomy, S18_tree)
##################################




###############################################################################
# Setting up a pre-processing monitoring method
#########################################
# Initialize summary table
S12_summary <- data.frame(SampleID = sample_names(S12_physeq_data), stringsAsFactors = FALSE)
S16_summary <- data.frame(SampleID = sample_names(S16_physeq_data), stringsAsFactors = FALSE)
S18_summary <- data.frame(SampleID = sample_names(S18_physeq_data), stringsAsFactors = FALSE)


# Impact tracking function
track_per_sample_step <- function(physeq, step, summary_table) {
  # Capture the name of the summary table variable
  summary_table_name <- deparse(substitute(summary_table))

  reads_per_sample <- sample_sums(physeq)  # Total reads per sample
  asvs_per_sample <- colSums(otu_table(physeq) > 0)  # ASV counts
  
  # Ensure alignment of data with the summary table
  reads_aligned <- reads_per_sample[match(summary_table$SampleID, names(reads_per_sample))]
  asvs_aligned <- asvs_per_sample[match(summary_table$SampleID, names(asvs_per_sample))]
  
  # Update the summary table
  summary_table[[paste0(step, "_Reads")]] <- reads_aligned
  summary_table[[paste0(step, "_ASVs")]] <- asvs_aligned
  
  # Assign the updated table back to the original variable
  assign(summary_table_name, summary_table, envir = .GlobalEnv)
}


# Track initial dataset and then re-apply to each future step
track_per_sample_step(S12_physeq_data, "Loaded", S12_summary)
track_per_sample_step(S16_physeq_data, "Loaded", S16_summary)
track_per_sample_step(S18_physeq_data, "Loaded", S18_summary)
###############################################################################




#############################################################################
# Taxonomic Cleanup
#####################
S12_taxfilt_data <- subset_taxa(S12_physeq_data, !(Genus %in% c("Homo", "Bos", "Sus", "Anas", "Canis", "Blumeria", "Squalius"))) %>% 
  subset_taxa(!(Family == "NA")) %>% 
  subset_taxa(!is.na(Class)) %>% 
  subset_taxa(!(Phylum == "Pseudomonadota"))

S16_taxfilt_data <- subset_taxa(S16_physeq_data, !(Genus %in% c("Homo", "Bos", "Sus", "Anas", "Canis"))) %>% 
  subset_taxa(!(Phylum == "Eukaryota")) %>% 
  subset_taxa(!(Family %in% c("Mitochondria", "Chloroplast"))) %>% 
  subset_taxa(!is.na(Class))

S18_taxfilt_data <- subset_taxa(S18_physeq_data, !(Genus %in% c("Homo", "Bos", "Sus", "Anas", "Canis"))) %>% 
  subset_taxa(!(Phylum == "Unassigned")) %>% 
  subset_taxa(!is.na(Class))


# Track processing step impact
track_per_sample_step(S12_taxfilt_data, "Taxonomy_cleanup", S12_summary)
track_per_sample_step(S16_taxfilt_data, "Taxonomy_cleanup", S16_summary)
track_per_sample_step(S18_taxfilt_data, "Taxonomy_cleanup", S18_summary)

# Cleanup if you want
#rm(S12_physeq_data, S16_physeq_data, S18_physeq_data)
###################################################################################



###############################################################################
# Sample sequencing depth filtering
#########################################
# Assess
(sample_sum_df <- data.frame(sum = sample_sums(S12_taxfilt_data)) %>% 
  arrange(desc(sum)))

(sample_sum_df <- data.frame(sum = sample_sums(S16_taxfilt_data)) %>% 
    arrange(desc(sum)))

(sample_sum_df <- data.frame(sum = sample_sums(S18_taxfilt_data)) %>% 
    arrange(desc(sum)))


# Find out more
#out <- boxplot.stats(sample_sum_df$sum)$out
#out_ind <- which(sample_sum_df$sum %in% c(out))
#outliers <- row.names(metadata[out_ind])

# remove samples with read depths threshold
# 2000 read depth requirement would limit loss but clean up important bad quality
(S12_deep_taxfilt_data <- prune_samples(sample_sums(S12_taxfilt_data) >= 2000 , S12_taxfilt_data))
(S16_deep_taxfilt_data <- prune_samples(sample_sums(S16_taxfilt_data) >= 2000 , S16_taxfilt_data))
(S18_deep_taxfilt_data <- prune_samples(sample_sums(S18_taxfilt_data) >= 2000 , S18_taxfilt_data))


# Track processing step impact
track_per_sample_step(S12_deep_taxfilt_data, "Depth_cleanup", S12_summary)
track_per_sample_step(S16_deep_taxfilt_data, "Depth_cleanup", S16_summary)
track_per_sample_step(S18_deep_taxfilt_data, "Depth_cleanup", S18_summary)

# Cleanup if you want
rm(S12_taxfilt_data, S16_taxfilt_data, S18_taxfilt_data)
##################################################################################



###############################################################################
# Removing ASVs bellow bleed-through threshold
#########################################
(S12_deep_bt_taxfilt_data <- filter_taxa(S12_deep_taxfilt_data, function(x) sum(x) >= 0.001 * mean(sample_sums(S12_deep_taxfilt_data)), TRUE))
(S16_deep_bt_taxfilt_data <- filter_taxa(S16_deep_taxfilt_data, function(x) sum(x) >= 0.001 * mean(sample_sums(S16_deep_taxfilt_data)), TRUE))
(S18_deep_bt_taxfilt_data <- filter_taxa(S18_deep_taxfilt_data, function(x) sum(x) >= 0.001 * mean(sample_sums(S18_deep_taxfilt_data)), TRUE))


# Track processing step impact
track_per_sample_step(S12_deep_bt_taxfilt_data, "Bleed-through_cleanup", S12_summary)
track_per_sample_step(S16_deep_bt_taxfilt_data, "Bleed-through_cleanup", S16_summary)
track_per_sample_step(S18_deep_bt_taxfilt_data, "Bleed-through_cleanup", S18_summary)

# Cleanup if you want
rm(S12_deep_taxfilt_data, S16_deep_taxfilt_data, S18_deep_taxfilt_data)
##################################################################################



###############################################################################
# ASV prevalence filtering
#########################################
# Define prevalence of each taxon (in how many samples did each taxon appear)
prev12S <- apply(otu_table(S12_deep_bt_taxfilt_data), 1, function(x) sum(x > 0))
# Filter taxa: Keep taxa present in more than 1 sample
(S12_filt_data <- filter_taxa(S12_deep_bt_taxfilt_data, function(x) sum(x > 0) > 1, TRUE))

prev16S <- apply(otu_table(S16_deep_bt_taxfilt_data), 1, function(x) sum(x > 0))
(S16_filt_data <- filter_taxa(S16_deep_bt_taxfilt_data, function(x) sum(x > 0) > 1, TRUE))

prev18S <- apply(otu_table(S18_deep_bt_taxfilt_data), 1, function(x) sum(x > 0))
(S18_filt_data <- filter_taxa(S18_deep_bt_taxfilt_data, function(x) sum(x > 0) > 1, TRUE))



# Track processing step impact
track_per_sample_step(S12_filt_data, "Prevalence_cleanup", S12_summary)
track_per_sample_step(S16_filt_data, "Prevalence_cleanup", S16_summary)
track_per_sample_step(S18_filt_data, "Prevalence_cleanup", S18_summary)

# Cleanup if you want
rm(S12_deep_bt_taxfilt_data, S16_deep_bt_taxfilt_data, S18_deep_bt_taxfilt_data, COI_deep_bt_taxfilt_data, 
   sample_sum_df, prev12S, prev16S, prev18S, prevCOI)
#################################################################################

################################################################################
# Optional cleanup
#rm(list = ls(pattern = "^(COI_|S12_|S16_|S18_)(physeq_data|physeq_Ldata|physeq_Sdata|summary)"))
#rm(metadata, track_per_sample_step)


###############################################################################
# Finalize and print summary datasets
####################
# Calculate % retained for reads and ASVs
S12_summary$Percent_Reads_Retained <- (S12_summary$Prevalence_cleanup_Reads / S12_summary$Loaded_Reads) * 100
S12_summary$Percent_ASVs_Retained <- (S12_summary$Prevalence_cleanup_ASVs / S12_summary$Loaded_ASVs) * 100

S16_summary$Percent_Reads_Retained <- (S16_summary$Prevalence_cleanup_Reads / S16_summary$Loaded_Reads) * 100
S16_summary$Percent_ASVs_Retained <- (S16_summary$Prevalence_cleanup_ASVs / S16_summary$Loaded_ASVs) * 100

S18_summary$Percent_Reads_Retained <- (S18_summary$Prevalence_cleanup_Reads / S18_summary$Loaded_Reads) * 100
S18_summary$Percent_ASVs_Retained <- (S18_summary$Prevalence_cleanup_ASVs / S18_summary$Loaded_ASVs) * 100


write.csv(S12_summary, "outputs/S12_sample_summary.csv", row.names = FALSE)
write.csv(S16_summary, "outputs/S16_sample_summary.csv", row.names = FALSE)
write.csv(S18_summary, "outputs/S18_sample_summary.csv", row.names = FALSE)


