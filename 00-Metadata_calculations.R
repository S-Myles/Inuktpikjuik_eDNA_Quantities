library(tidyverse)
library(lubridate)
library(ggpubr)      # For data checks
library(GGally)      # For data checks

###################
### Data + metadata import, tidy, and format
###################
metadata <-  read_tsv("data/metadata_and_qPCR.txt")

# remove empty rows
if (all(is.na(tail(metadata, 1)))) {
  metadata <- metadata %>%
    slice(-n())
}

# Convert columns' data types accordingly and add rxn volumes
metadata <- metadata %>%
  mutate(across(c(latitude, longitude, `18S_DADA2_reads`, `18S_spike_main_ASV`, 
                  `16S_DADA2_reads`, `16S_spike_main_ASV`, `12S_DADA2_reads`), 
                ~ as.numeric(na_if(as.character(.), "N/A"))),
         Location = factor(Location, levels = c("Cooler blank 1",
                                               "Cooler blank 2",
                                               "Sackville River", 
                                               "DRDC barge", 
                                               "McKay bridge", 
                                               "Tufts Cove",
                                               "McNabs Island")),
         Date = ymd(Date),
         uL_per_qPCR_reaction = 5,
         uL_per_PCR_reaction = 2
  )
         


###################################################################################################
##      Gene quantities calculations
###################################################################################################
###################
### 16Sin qPCR calculated Extraction Efficiencies
###################
# Extr. efficiency coefficient = (16Sin_qPCRquantity for 5uL <- dilution <- elution vol) / spikes added
metadata <- metadata %>% 
  mutate(`16Sin extraction efficiency` = ((`16Sin_qPCR_quantity` / COI_16S_18S_qPCR_dilutions ) /   
           uL_per_qPCR_reaction * 
           elution_volume) / 
           `16S_spike_copies`
  )

# Check the data
hist(metadata$`16Sin extraction efficiency`, xlab = "Extraction Efficiency Coefficient")


###################
### Nanodrop DNA concentration method
###################
# [DNA] <- Elution volume <- Volume filtered = [DNA] in SW
metadata <- metadata %>% 
  mutate(`DNA (ng/L SW)` = (nanodrop_dna * elution_volume)   # DNA ng per filter
         / `Volume filtered (ml SW)` * 1000)                       # DNA ng / L Seawater

# Data check
hist(metadata$`DNA (ng/L SW)`)


#doubtful of bellow
#BBT_data <- BBT_data %>% 
#  mutate(`dna_per_L_SW_with_efficiencies` = `dna (ng/L SW)` * 1 / `16Sin extraction efficiency`)  #factoring in extraction efficiency
#hist(BBT_data$`dna_per_L_SW_with_efficiencies`)


###################
### qPCR gene counts method
###################

# qPCR quantities <- tmpl vol <- dilution <- Vol filt
### WITHOUT EXTRACTION EFFICIENCIES ###
metadata <- metadata %>% 
  mutate(`12S copies/mL SW` = (`12S_qPCR_quantity` / `12S_dilutions`) /   # gene copies per 5 uL considering tmpl dilution
           uL_per_qPCR_reaction *                                                 # copies per uL
           elution_volume /                                                       # copies on filter
           `Volume filtered (ml SW)`,                                                    # copies / mL of SW
         
         `16S copies/mL SW` = (`16S_qPCR_quantity` / COI_16S_18S_qPCR_dilutions) /   # gene copies per 5 uL considering tmpl dilution
           uL_per_qPCR_reaction *                                                 # copies per uL
           elution_volume /                                                       # copies on filter
           `Volume filtered (ml SW)`,                                                    # copies / mL of SW
         
         `18S copies/mL SW` = (`18S_qPCR_quantity` / COI_16S_18S_qPCR_dilutions) /   # gene copies per 5 uL considering tmpl dilution
           uL_per_qPCR_reaction *                                                 # copies per uL
           elution_volume /                                                       # copies on filter
           `Volume filtered (ml SW)`,                                                    # copies / mL of SW
         
         `COI copies/mL SW` = (`COI_qPCR_quantity` / COI_16S_18S_qPCR_dilutions) /   # gene copies per 5 uL considering tmpl dilution
           uL_per_qPCR_reaction *                                                 # copies per uL
           elution_volume /                                                       # copies on filter
           `Volume filtered (ml SW)`                                                    # copies / mL of SW
  )


# qPCR quantities <- tmpl vol <- dilution <- extr. eff. <- Vol filt
### WITH EXTRACTION EFFICIENCIES ###
metadata <- metadata %>% 
  mutate(`12S copies/mL SW * eff` = (`12S_qPCR_quantity` / `12S_dilutions`) /   # gene copies per 5 uL considering tmpl dilution
           uL_per_qPCR_reaction *                                                 # copies per uL
           elution_volume *                                                       # copies on filter
           1/`16Sin extraction efficiency` /                                      # still copies on filter (factoring: efficiency of extracting spike DNA)
           `Volume filtered (ml SW)`,                                                    # copies / mL of SW
         
         `16S copies/mL SW * eff` = (`16S_qPCR_quantity` / COI_16S_18S_qPCR_dilutions) /   # gene copies per 5 uL considering tmpl dilution
           uL_per_qPCR_reaction *                                                 # copies per uL
           elution_volume *                                                       # copies on filter
           1/`16Sin extraction efficiency` /                                      # still copies on filter (factoring: efficiency of extracting spike DNA)
           `Volume filtered (ml SW)`,                                                    # copies / mL of SW
         
         `18S copies/mL SW * eff` = (`18S_qPCR_quantity` / COI_16S_18S_qPCR_dilutions) /   # gene copies per 5 uL considering tmpl dilution
           uL_per_qPCR_reaction *                                                 # copies per uL
           elution_volume *                                                       # copies on filter
           1/`16Sin extraction efficiency` /                                      # still copies on filter (factoring: efficiency of extracting spike DNA)
           `Volume filtered (ml SW)`,                                                    # copies / mL of SW
         
         `COI copies/mL SW * eff` = (`COI_qPCR_quantity` / COI_16S_18S_qPCR_dilutions) /   # gene copies per 5 uL considering tmpl dilution
           uL_per_qPCR_reaction *                                                 # copies per uL
           elution_volume *                                                       # copies on filter
           1/`16Sin extraction efficiency` /                                      # still copies on filter (factoring: efficiency of extracting spike DNA)
           `Volume filtered (ml SW)`                                                    # copies / mL of SW
  )


# Check the data
hist(metadata$`12S copies/mL SW * eff`)
hist(metadata$`16S copies/mL SW * eff`)
hist(metadata$`18S copies/mL SW * eff`)
hist(metadata$`COI copies/mL SW * eff`)



###################
### ASV sequence counts method
###################

# Simplest
metadata <- metadata %>% 
  mutate(`16S_spike_ratios` = (`16S_spike_main_ASV` / `16S_DADA2_reads`), # Spike to total reads ratios
         
         `18S_spike_ratios` = (`18S_spike_main_ASV` / `18S_DADA2_reads`), # Spike to total reads ratios
         
         `16S_quant_reads` = (`16S_DADA2_reads` * `16S_spike_copies` /
                                          `16S_spike_main_ASV`),          # Spike corrected read counts
         
         `18S_quant_reads` = (`18S_DADA2_reads` * `18S_spike_copies` /
                                          `18S_spike_main_ASV`),          # Spike corrected read counts
         )


# Spike ASV counts <- total read counts <- ratio * Other ASV reads <- Extr. eff. <- volume filtered

metadata <- metadata %>% 
  mutate(`16S seq/mL SW` = (`16S_spike_copies` * `16S_DADA2_reads` /
                              `16S_spike_main_ASV`) *                # Spike corrected read counts
           1/`16Sin extraction efficiency` /                             # Extr. eff. corrected read counts             
           COI_16S_18S_PCR_dilutions /                                   # Considering PCR template dilutions
           `Volume filtered (ml SW)`,                                           # Seqs per mL seawater
         
          `18S seq/mL SW` = (`18S_spike_copies` * `18S_DADA2_reads` /
                                   `18S_spike_main_ASV`) *                # Spike corrected read counts
            1/`16Sin extraction efficiency` /                             # Extr. eff. corrected read counts             
            COI_16S_18S_PCR_dilutions /                                   # Considering PCR template dilutions
           `Volume filtered (ml SW)` ,                                            # Seqs per mL seawater

### Think here, should I consider PCR tmpl dilutions and equimolar pooling at Seq library prep??? ###
  )


write_csv(metadata, "data/metadata_qPCR_calculated.csv")


