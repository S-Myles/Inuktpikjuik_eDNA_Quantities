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
  filter(Depth..m. < 2)

# Create the faceted bar plot with custom facet titles
(plot <- ggplot(ps_melted, aes(x = Date, y = Abundance, fill = Species)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(~ Location) +
    labs(x = "Samples", y = "Relative Abundance") +
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
  width = 8, height = 6,
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
  filter(Depth..m. < 2)

# Create the faceted bar plot with custom facet titles
(plot <- ggplot(ps_melted, aes(x = Date, y = Abundance, fill = Phylum)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid( ~ Location) +  # Add the custom labeller here
    labs(x = "Samples", y = "Relative Abundance") +
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
  width = 8, height = 6,
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
  filter(Depth..m. < 2)

# Create the faceted bar plot with custom facet titles
(plot <- ggplot(ps_melted, aes(x = Date, y = Abundance, fill = Phylum)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid( ~ Location) +  # Add the custom labeller here
    labs(x = "Samples", y = "Relative Abundance") +
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
  width = 8, height = 6,
  dpi = 300
)


