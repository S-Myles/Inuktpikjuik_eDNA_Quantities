########
# Global Settings
########
# Libraries
library(phyloseq)
library(tidyverse)
# Global plot theme setting
theme_set(theme_bw())

colorblind_palette <- c("#117733", "#44AA99", "#88CCEE", "#332288",
                        "#DDCC77", "#999933", "#CC6677", "#882255",
                        "#AA4499", "#DDDDDD") 

set3_palette <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
                  "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
                  "#D9D9D9", "#BC80BD")

nature_palette <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                    "#66A61E", "#E6AB02", "#A6761D", "#666666",
                    "#FF7F00", "#6A3D9A")

########
# Cooler blanks top 10
########
# -------------- 12S Taxa barplots
########
# Aggregate taxa to the desired level and transform to relative abundances
S12_cooler_blanks <- S12_taxfilt_data %>%
  subset_samples(Location %in% c("Cooler blank 1", "Cooler blank 2")) %>% 
  #tax_glom(taxrank = "Species") %>%  # Aggregate at the Species level
  transform_sample_counts(function(x) x / sum(x))  # Convert to relative abundances

# Melt the data for ggplot
ps_melted <- psmelt(S12_cooler_blanks) %>% 
  mutate(Location = factor(Location, levels = c("Cooler blank 1",
                                                "Cooler blank 2")),
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


# Identify the top 10 ASVs by cumulative relative abundance across the Cooler blank samples
top10_asvs <- ps_melted %>%
  group_by(Species) %>%                          # Group by ASV (Species)
  filter(Abundance > 0) %>%                         # take out zeros
  summarize(total_abundance = sum(Abundance)) %>% # Sum their relative abundances
  arrange(desc(total_abundance)) %>%            # Order descending
  slice_head(n = 10) %>%                          # Select the top 10
  pull(Species) %>%                                  # Extract their names
  factor(levels = c("Fundulus heteroclitus",
                    "Triglops murrayi",
                    "Pungitius pungitius",
                    "Pseudopleuronectes americanus", 
                    "Pollachius virens", 
                    "Pholis gunnellus",
                    "Hyperoplus lanceolatus",
                    "Clupea harengus",
                    "Scomber scombrus",
                    "NA"))


# Create the faceted bar plot
(plot <- ggplot(ps_melted, aes(x = Location, y = Abundance, fill = Species)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~ Date) +
  labs(x = "Sampling Location", y = "Relative Abundance") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10),
    legend.position = "right"
  ) +
  # Only include top 10 in the legend; assign grey for "Others"
  scale_fill_manual(values = c(setNames(colorblind_palette, top10_asvs))))




# -------------- 16S Taxa barplots
########
# Aggregate taxa to the desired level and transform to relative abundances
S16_cooler_blanks <- S16_taxfilt_data %>%
  subset_samples(Location %in% c("Cooler blank 1", "Cooler blank 2")) %>% 
  #tax_glom(taxrank = "Species") %>%  # Aggregate at the Species level
  transform_sample_counts(function(x) x / sum(x))  # Convert to relative abundances

# Melt the data for ggplot
ps_melted <- psmelt(S16_cooler_blanks) %>% 
  mutate(Location = factor(Location, levels = c("Cooler blank 1",
                                                "Cooler blank 2")),
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


# Identify the top 10 ASVs by cumulative relative abundance across the Cooler blank samples
top10_asvs <- ps_melted %>%
  group_by(Phylum) %>%                          # Group by ASV (Species)
  filter(Abundance > 0) %>%                         # take out zeros
  summarize(total_abundance = sum(Abundance)) %>% # Sum their relative abundances
  arrange(desc(total_abundance)) %>%            # Order descending
  slice_head(n = 10) %>%                          # Select the top 10
  pull(Phylum)  %>%                             # Extract their names
  factor(levels = c("Desulfobacterota",
                    "Cyanobacteria",
                    "Nitrospinota",
                    "Crenarchaeota",
                    "WPS-2",
                    "Firmicutes", 
                    "Actinobacteriota", 
                    "Planctomycetota",
                    "Bacteroidota",
                    "Proteobacteria"))




# Create the faceted bar plot
(plot <- ggplot(ps_melted, aes(x = Location, y = Abundance, fill = Phylum)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(~ Date) +
    labs(x = "Sampling Location", y = "Relative Abundance") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 10),
      legend.position = "right"
    ) +
    # Only include top 10 in the legend; assign grey for "Others"
    scale_fill_manual(values = c(setNames(set3_palette, top10_asvs))))




# -------------- 18S Taxa barplots
########
# Aggregate taxa to the desired level and transform to relative abundances
S18_cooler_blanks <- S18_taxfilt_data %>%
  subset_samples(Location %in% c("Cooler blank 1", "Cooler blank 2")) %>% 
  #tax_glom(taxrank = "Species") %>%  # Aggregate at the Species level
  transform_sample_counts(function(x) x / sum(x))  # Convert to relative abundances

# Melt the data for ggplot
ps_melted <- psmelt(S18_cooler_blanks) %>% 
  mutate(Location = factor(Location, levels = c("Cooler blank 1",
                                                "Cooler blank 2")),
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


# Identify the top 10 ASVs by cumulative relative abundance across the Cooler blank samples
top10_asvs <- ps_melted %>%
  group_by(Phylum) %>%                          # Group by ASV (Species)
  filter(Abundance > 0) %>%                         # take out zeros
  summarize(total_abundance = sum(Abundance)) %>% # Sum their relative abundances
  arrange(desc(total_abundance)) %>%            # Order descending
  slice_head(n = 10) %>%                          # Select the top 10
  pull(Phylum)  %>%                        # Extract their names
  factor(levels = c("Phragmoplastophyta",
                    "Rotifera",
                    "Ochrophyta",
                    "Protalveolata",
                    "Diatomea",
                    "Mollusca", 
                    "Ciliophora", 
                    "Arthropoda",
                    "Dinoflagellata",
                    "Ascomycota"))


# Create the faceted bar plot
(plot <- ggplot(ps_melted, aes(x = Location, y = Abundance, fill = Phylum)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(~ Date) +
    labs(x = "Sampling Location", y = "Relative Abundance") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 10),
      legend.position = "right"
    ) +
    # Only include top 10 in the legend; assign grey for "Others"
    scale_fill_manual(values = c(setNames(nature_palette, top10_asvs))))








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
  filter(Depth..m. == 1) %>% 
  mutate(Location = factor(Location, levels = c("Sackville River", 
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

# Identify the top 10 ASVs by cumulative relative abundance across the Cooler blank samples
top10_asvs <- ps_melted %>%
  group_by(Species) %>%                          # Group by ASV (Species)
  filter(Abundance > 0) %>%                         # take out zeros
  summarize(total_abundance = sum(Abundance)) %>% # Sum their relative abundances
  arrange(desc(total_abundance)) %>%            # Order descending
  slice_head(n = 10) %>%                          # Select the top 10
  pull(Species)

# Create the faceted bar plot with custom facet titles
(plot <- ggplot(ps_melted, aes(x = Location, y = Abundance, fill = Species)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(~ Date) +
    labs(x = "Sampling Location", y = "Relative Abundance") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 10)#,
      #legend.position = "none"
    ) +
  # Only include top 10 in the legend; assign grey for "Others"
  scale_fill_manual(values = c(setNames(colorblind_palette, top10_asvs))))

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
  mutate(quantitative_abundance = Abundance * X12S.copies.L.SW) %>% 
  filter(Depth..m. == 1) %>% 
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
    ) +
    # Only include top 10 in the legend; assign grey for "Others"
    scale_fill_manual(values = c(setNames(colorblind_palette, top10_asvs))))

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
  mutate(quantitative_abundance = Abundance * X16S.copies.L.SW) %>% ## Make it quantitative w qPCR data
  filter(Depth..m. == 1) %>% 
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


# Identify the top 10 ASVs by cumulative relative abundance across the Cooler blank samples
top10_asvs <- ps_melted %>%
  group_by(Phylum) %>%                          # Group by ASV (Species)
  filter(Abundance > 0) %>%                         # take out zeros
  summarize(total_abundance = sum(Abundance)) %>% # Sum their relative abundances
  arrange(desc(total_abundance)) %>%            # Order descending
  slice_head(n = 10) %>%                          # Select the top 10
  pull(Phylum)


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
    ) +
    # Only include top 10 in the legend; assign grey for "Others"
    scale_fill_manual(values = c(setNames(set3_palette, top10_asvs))))

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
    ) +
    # Only include top 10 in the legend; assign grey for "Others"
    scale_fill_manual(values = c(setNames(set3_palette, top10_asvs))))

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
  mutate(quantitative_abundance = Abundance * X18S.copies.L.SW) %>% ## Make it quantitative w qPCR data
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


# Identify the top 10 ASVs by cumulative relative abundance across the Cooler blank samples
top10_asvs <- ps_melted %>%
  group_by(Phylum) %>%                          # Group by ASV (Species)
  filter(Abundance > 0) %>%                         # take out zeros
  summarize(total_abundance = sum(Abundance)) %>% # Sum their relative abundances
  arrange(desc(total_abundance)) %>%            # Order descending
  slice_head(n = 10) %>%                          # Select the top 10
  pull(Phylum)


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
    ) +
    # Only include top 10 in the legend; assign grey for "Others"
    scale_fill_manual(values = c(setNames(nature_palette, top10_asvs))))

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
    ) +
    # Only include top 10 in the legend; assign grey for "Others"
    scale_fill_manual(values = c(setNames(nature_palette, top10_asvs))))

# Save the plot
ggsave(
  filename = "outputs/metabarcoding_barplots/S18-Quantitative_sp_barplot.png",
  plot = plot,
  width = 14, height = 6,
  dpi = 300
)
