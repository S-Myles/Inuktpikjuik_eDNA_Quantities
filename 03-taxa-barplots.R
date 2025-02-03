########
# Global Settings
########
# Libraries
library(phyloseq)
library(tidyverse)
# Global plot theme setting
theme_set(theme_bw())


########
# 12S
########
# -------------- Taxa barplots
########
# Aggregate taxa to the desired level and transform to relative abundances
S12_physeq_agg <- S12_taxfilt_data %>%
  #tax_glom(taxrank = "Species") %>%  # Aggregate at the Species level
  transform_sample_counts(function(x) x / sum(x))  # Convert to relative abundances

# Melt the data for ggplot
ps_melted <- psmelt(S12_physeq_agg) %>% 
  filter(Depth..m. < 2) %>% 
  mutate(Location = factor(Location, levels = c("Cooler blank 1",
                                                "Cooler blank 2",
                                                "Sackville River", 
                                                "DRDC barge", 
                                                "McKay bridge", 
                                                "Tufts Cove",
                                                "McNabs Island")),
         Date = format(Date, "%Y-%b-%d"),
         Date = factor(Date, levels = c("2023-Jun-22",
                                        "2023-Jul-25",
                                        "2023-Aug-31",
                                        "2023-Sep-25",
                                        "2023-Nov-03",
                                        "2023-Nov-24",
                                        "2023-Dec-14",
                                        "2024-Jan-19",
                                        "2024-Feb-16",
                                        "2024-Mar-28",
                                        "2024-Apr-18",
                                        "2024-May-16"
                                        )))

# Create the faceted bar plot with custom facet titles
(plot <- ggplot(ps_melted, aes(x = Location, y = Abundance, fill = Species)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(~ Date) +
    labs(x = "Sampling Location", y = "Relative Abundance") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 10),
      legend.position = "none"
    ))

# Save the plot
ggsave(
  filename = "outputs/metabarcoding_barplots/S12-sp_barplot.png",
  plot = plot,
  width = 14, height = 6,
  dpi = 300
)



## Make it quantitative w qPCR data
# Melt the data for ggplot
ps_melted <- psmelt(S12_physeq_agg) %>% 
  mutate(quantitative_abundance = Abundance * X12S.copies.mL.SW) %>% 
  filter(Depth..m. < 2) %>% 
  mutate(Location = factor(Location, levels = c("Cooler blank 1",
                                                "Cooler blank 2",
                                                "Sackville River", 
                                                "DRDC barge", 
                                                "McKay bridge", 
                                                "Tufts Cove",
                                                "McNabs Island")))

# Create the faceted bar plot with custom facet titles
(plot <- ggplot(ps_melted, aes(x = Location, y = quantitative_abundance, fill = Species)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(~ Date) +
    labs(x = "Sampling Location", y = "Quantitative Abundance") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 10),
      legend.position = "none"
    ))

# Save the plot
ggsave(
  filename = "outputs/metabarcoding_barplots/S12-Quantitative_sp_barplot.png",
  plot = plot,
  width = 14, height = 6,
  dpi = 300
)





########
# 16S
########
# -------------- Taxa barplots
########
# Aggregate taxa to the desired level and transform to relative abundances
S16_physeq_agg <- S16_taxfilt_data %>%
  #tax_glom(taxrank = "Phylum") %>%  # Aggregate at the Species level
  transform_sample_counts(function(x) x / sum(x))  # Convert to relative abundances

# Melt the data for ggplot
ps_melted <- psmelt(S16_physeq_agg) %>% 
  mutate(quantitative_abundance = Abundance * X16S.copies.mL.SW) %>% ## Make it quantitative w qPCR data
  filter(Depth..m. < 2) %>% 
  mutate(Location = factor(Location, levels = c("Cooler blank 1",
                                                "Cooler blank 2",
                                                "Sackville River", 
                                                "DRDC barge", 
                                                "McKay bridge", 
                                                "Tufts Cove",
                                                "McNabs Island")),
         Date = format(Date, "%Y-%b-%d"),
         Date = factor(Date, levels = c("2023-Jun-22",
                                        "2023-Jul-25",
                                        "2023-Aug-31",
                                        "2023-Sep-25",
                                        "2023-Nov-03",
                                        "2023-Nov-24",
                                        "2023-Dec-14",
                                        "2024-Jan-19",
                                        "2024-Feb-16",
                                        "2024-Mar-28",
                                        "2024-Apr-18",
                                        "2024-May-16"
         )))


# Create the faceted bar plot with custom facet titles
(plot <- ggplot(ps_melted, aes(x = Location, y = Abundance, fill = Phylum)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid( ~ Date) +  # Add the custom labeller here
    labs(x = "Sampling Location", y = "Relative Abundance") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 10),
      legend.position = "none"
    ))

# Save the plot
ggsave(
  filename = "outputs/metabarcoding_barplots/S16-phylum_barplot.png",
  plot = plot,
  width = 14, height = 6,
  dpi = 300
)



# Create the faceted bar plot with custom facet titles
(plot <- ggplot(ps_melted, aes(x = Location, y = quantitative_abundance, fill = Phylum)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(~ Date) +
    labs(x = "Sampling Location", y = "Quantitative Abundance") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 10),
      legend.position = "none"
    ))

# Save the plot
ggsave(
  filename = "outputs/metabarcoding_barplots/S16-Quantitative_sp_barplot.png",
  plot = plot,
  width = 14, height = 6,
  dpi = 300
)







########
# 18S
########
# -------------- Taxa barplots
########
# Aggregate taxa to the desired level and transform to relative abundances
S18_physeq_agg <- S18_taxfilt_data %>%
  #tax_glom(taxrank = "Phylum") %>%  # Aggregate at the Species level
  transform_sample_counts(function(x) x / sum(x))  # Convert to relative abundances

# Melt the data for ggplot
ps_melted <- psmelt(S18_physeq_agg) %>% 
  mutate(quantitative_abundance = Abundance * X18S.copies.mL.SW) %>% ## Make it quantitative w qPCR data
  filter(Depth..m. < 2) %>% 
  mutate(Location = factor(Location, levels = c("Cooler blank 1",
                                                "Cooler blank 2",
                                                "Sackville River", 
                                                "DRDC barge", 
                                                "McKay bridge", 
                                                "Tufts Cove",
                                                "McNabs Island")),
         Date = format(Date, "%Y-%b-%d"),
         Date = factor(Date, levels = c("2023-Jun-22",
                                        "2023-Jul-25",
                                        "2023-Aug-31",
                                        "2023-Sep-25",
                                        "2023-Nov-03",
                                        "2023-Nov-24",
                                        "2023-Dec-14",
                                        "2024-Jan-19",
                                        "2024-Feb-16",
                                        "2024-Mar-28",
                                        "2024-Apr-18",
                                        "2024-May-16"
         )))

# Create the faceted bar plot with custom facet titles
(plot <- ggplot(ps_melted, aes(x = Location, y = Abundance, fill = Phylum)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid( ~ Date) +  # Add the custom labeller here
    labs(x = "Sampling Location", y = "Relative Abundance") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 10),
      legend.position = "none"
    ))

# Save the plot
ggsave(
  filename = "outputs/metabarcoding_barplots/S18-phylum_barplot.png",
  plot = plot,
  width = 14, height = 6,
  dpi = 300
)


# Create the faceted bar plot with custom facet titles
(plot <- ggplot(ps_melted, aes(x = Location, y = quantitative_abundance, fill = Phylum)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(~ Date) +
    labs(x = "Sampling Location", y = "Quantitative Abundance") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 10),
      legend.position = "none"
    ))

# Save the plot
ggsave(
  filename = "outputs/metabarcoding_barplots/S16-Quantitative_sp_barplot.png",
  plot = plot,
  width = 14, height = 6,
  dpi = 300
)
