########
# Global Settings
########
# Libraries
library(phyloseq)
library(tidyverse)
# Global plot theme setting
theme_set(theme_bw())

colorblind_palette <- c("#117733", "#44AA99", "#88CCEE", "red",
                        "#DDCC77", "#999933", "#CC6677", "#882255",
                        "#AA4499", "#DDDDDD", "grey15") 

set3_palette <- c("#8DD3C7", "#FFFFB3", "red", "#FB8072",
                  "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
                  "#D9D9D9", "#BC80BD", "grey15")

nature_palette <- c("#1B9E77", "red", "#7570B3", "#E7298A",
                    "#66A61E", "#E6AB02", "#A6761D", "#666666",
                    "#FF7F00", "#6A3D9A", "grey15")


########
# -------------- 12S Taxa barplots
########
# Aggregate taxa to the desired level and transform to relative abundances
S12_rel_abun <- S12_taxfilt_data %>%
  transform_sample_counts(function(x) x / sum(x))  # Convert to relative abundances

ps_melted <- psmelt(S12_rel_abun) %>% 
  filter(!is.na(Species))

# Top 10 taxa to plot
top10_asvs <- ps_melted %>%
  group_by(Species) %>% 
  summarize(total_abundance = sum(Abundance)) %>% # Sum relative abundances
  arrange(desc(total_abundance)) %>%         # Order in descending order
  slice_head(n = 10) %>%                      # Select the top 10 ASVs
  pull(Species)                              # Extract sp names



########
# Cooler blanks top 10
########
# Melt the data for ggplot
ps_melted <- psmelt(S12_rel_abun) %>% 
  filter(Location %in% c("Cooler blank 1", "Cooler blank 2")) %>% 
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
                                        )),
         Species = ifelse(Species %in% top10_asvs, Species, "Others"),
         Species = factor(Species, levels = c("Others", rev(top10_asvs))),
         quantitative_abundance = Abundance * X12S.copies.L.SW)


# Create the faceted bar plot
(plot <- ggplot(ps_melted, aes(x = Location, y = Abundance, fill = Species)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~ Date) +
  labs(x = "Sampling Location", y = "Relative Abundance") +
  theme_minimal() +
    theme(
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 14, angle = 45, hjust = 1, face = "bold"),  # Bold x-axis values
      axis.text.y = element_text(size = 14, face = "bold"),  # Bold y-axis values
      axis.title.x = element_text(size = 14, face = "bold"),    # Bold axis labels
      axis.title.y = element_text(size = 14, face = "bold"),    # Bold axis labels
      strip.text = element_text(size = 10, face = "bold")   
    ) +
  scale_fill_manual(values = c(setNames(colorblind_palette, c(top10_asvs, "Others")))))

# Save the plot
ggsave(
  filename = "outputs/metabarcoding_barplots/S12-BL_rel-abun_sp_barplot.png",
  plot = plot,
  width = 16, height = 6,
  dpi = 300
)



# Quantitative plot
(plot <- ggplot(ps_melted, aes(x = Location, y = quantitative_abundance, fill = Species)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(~ Date) +
    labs(x = "Sampling Location", y = "eDNA gene concentration in seawater (copies/L)") +
    theme_minimal() +
    theme(
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 14, angle = 45, hjust = 1, face = "bold"),  # Bold x-axis values
      axis.text.y = element_text(size = 14, face = "bold"),  # Bold y-axis values
      axis.title.x = element_text(size = 14, face = "bold"),    # Bold axis labels
      axis.title.y = element_text(size = 14, face = "bold"),    # Bold axis labels
      strip.text = element_text(size = 10, face = "bold")   
    ) +
    scale_fill_manual(values = c(setNames(colorblind_palette, c(top10_asvs, "Others")))))

# Save the plot
ggsave(
  filename = "outputs/metabarcoding_barplots/S12-BL_copies_sp_barplot.png",
  plot = plot,
  width = 16, height = 6,
  dpi = 300
)


########
# Environmental genes top 10
########

# Melt the data for ggplot
ps_melted <- psmelt(S12_rel_abun) %>% 
  filter(Depth..m. == 1) %>% 
  filter(Location %in% c("Sackville River", 
                         "DRDC barge", 
                         "McKay bridge", 
                         "Tufts Cove",
                         "McNabs Island")) %>% 
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
                                        "2024-May-16")),
         Species = ifelse(Species %in% top10_asvs, Species, "Others"),
         Species = factor(Species, levels = c("Others", rev(top10_asvs))),
         quantitative_abundance = Abundance * X12S.copies.L.SW)


# Create the faceted bar plot
(plot <- ggplot(ps_melted, aes(x = Location, y = Abundance, fill = Species)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(~ Date) +
    labs(x = "Sampling Location", y = "Relative Abundance") +
    theme_minimal() +
    theme(
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "bold"),  # Bold x-axis values
      axis.text.y = element_text(size = 14, face = "bold"),  # Bold y-axis values
      axis.title.x = element_text(size = 14, face = "bold"),    # Bold axis labels
      axis.title.y = element_text(size = 14, face = "bold"),    # Bold axis labels
      strip.text = element_text(size = 10, face = "bold")   
    ) +
    scale_fill_manual(values = c(setNames(colorblind_palette, c(top10_asvs, "Others")))))

# Save the plot
ggsave(
  filename = "outputs/metabarcoding_barplots/S12-sp_barplot.png",
  plot = plot,
  width = 16, height = 6,
  dpi = 300
)


# Quantitative plot
(plot <- ggplot(ps_melted, aes(x = Location, y = quantitative_abundance, fill = Species)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(~ Date) +
    labs(x = "Sampling Location", y = "eDNA gene concentration in seawater (copies/L)") +
    theme_minimal() +
    theme(
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "bold"),  # Bold x-axis values
      axis.text.y = element_text(size = 14, face = "bold"),  # Bold y-axis values
      axis.title.x = element_text(size = 14, face = "bold"),    # Bold axis labels
      axis.title.y = element_text(size = 14, face = "bold"),    # Bold axis labels
      strip.text = element_text(size = 10, face = "bold")   
    ) +
    scale_fill_manual(values = c(setNames(colorblind_palette, c(top10_asvs, "Others")))))

# Save the plot
ggsave(
  filename = "outputs/metabarcoding_barplots/S12-Quantitative_sp_barplot.png",
  plot = plot,
  width = 16, height = 6,
  dpi = 300
)









########
# -------------- 16S Taxa barplots
########
# Aggregate taxa to the desired level and transform to relative abundances
S16_rel_abun <- S16_taxfilt_data %>%
  transform_sample_counts(function(x) x / sum(x))  # Convert to relative abundances

ps_melted <- psmelt(S16_rel_abun) %>% 
  filter(!is.na(Phylum))

# Top 10 taxa to plot
top10_asvs <- ps_melted %>%
  group_by(Phylum) %>% 
  summarize(total_abundance = sum(Abundance)) %>% # Sum relative abundances
  arrange(desc(total_abundance)) %>%         # Order in descending order
  slice_head(n = 10) %>%                      # Select the top 10 ASVs
  pull(Phylum)                              # Extract sp names



########
# Cooler blanks top 10
########
# Melt the data for ggplot
ps_melted <- psmelt(S16_rel_abun) %>% 
  filter(Location %in% c("Cooler blank 1", "Cooler blank 2")) %>% 
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
         ))) %>% 
  mutate(Phylum = ifelse(Phylum %in% top10_asvs, Phylum, "Others"),
         Phylum = factor(Phylum, levels = c(rev(top10_asvs), "Others")),
         quantitative_abundance = Abundance * X16S.copies.L.SW)


# Create the faceted bar plot
(plot <- ggplot(ps_melted, aes(x = Location, y = Abundance, fill = Phylum)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(~ Date) +
    labs(x = "Sampling Location", y = "Relative Abundance") +
    theme_minimal() +
    theme(
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 14, angle = 45, hjust = 1, face = "bold"),  # Bold x-axis values
      axis.text.y = element_text(size = 14, face = "bold"),  # Bold y-axis values
      axis.title.x = element_text(size = 14, face = "bold"),    # Bold axis labels
      axis.title.y = element_text(size = 14, face = "bold"),    # Bold axis labels
      strip.text = element_text(size = 10, face = "bold")   
    ) +
    scale_fill_manual(values = c(setNames(set3_palette, c(top10_asvs, "Others")))))

# Save the plot
ggsave(
  filename = "outputs/metabarcoding_barplots/S16-BL_rel-abun_sp_barplot.png",
  plot = plot,
  width = 16, height = 6,
  dpi = 300
)



# Quantitative plot
(plot <- ggplot(ps_melted, aes(x = Location, y = quantitative_abundance, fill = Phylum)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(~ Date) +
    labs(x = "Sampling Location", y = "eDNA gene concentration in seawater (copies/L)") +
    theme_minimal() +
    theme(
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 14, angle = 45, hjust = 1, face = "bold"),  # Bold x-axis values
      axis.text.y = element_text(size = 14, face = "bold"),  # Bold y-axis values
      axis.title.x = element_text(size = 14, face = "bold"),    # Bold axis labels
      axis.title.y = element_text(size = 14, face = "bold"),    # Bold axis labels
      strip.text = element_text(size = 10, face = "bold")   
    ) +
    scale_fill_manual(values = c(setNames(set3_palette, c(top10_asvs, "Others")))))

# Save the plot
ggsave(
  filename = "outputs/metabarcoding_barplots/S16-BL_copies_sp_barplot.png",
  plot = plot,
  width = 16, height = 6,
  dpi = 300
)



########
# Environmental genes top 10
########

# Melt the data for ggplot
ps_melted <- psmelt(S16_rel_abun) %>% 
  filter(Depth..m. == 1) %>% 
  filter(Location %in% c("Sackville River", 
                         "DRDC barge", 
                         "McKay bridge", 
                         "Tufts Cove",
                         "McNabs Island")) %>% 
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
                                        "2024-May-16")),
         Phylum = ifelse(Phylum %in% top10_asvs, Phylum, "Others"),
         Phylum = factor(Phylum, levels = c("Others", rev(top10_asvs))),
         quantitative_abundance = Abundance * X16S.copies.L.SW)


# Create the faceted bar plot
(plot <- ggplot(ps_melted, aes(x = Location, y = Abundance, fill = Phylum)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(~ Date) +
    labs(x = "Sampling Location", y = "Relative Abundance") +
    theme_minimal() +
    theme(
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "bold"),  # Bold x-axis values
      axis.text.y = element_text(size = 14, face = "bold"),  # Bold y-axis values
      axis.title.x = element_text(size = 14, face = "bold"),    # Bold axis labels
      axis.title.y = element_text(size = 14, face = "bold"),    # Bold axis labels
      strip.text = element_text(size = 10, face = "bold")   
    ) +
    scale_fill_manual(values = c(setNames(set3_palette, c(top10_asvs, "Others")))))

# Save the plot
ggsave(
  filename = "outputs/metabarcoding_barplots/S16-sp_barplot.png",
  plot = plot,
  width = 16, height = 6,
  dpi = 300
)


# Quantitative plot
(plot <- ggplot(ps_melted, aes(x = Location, y = quantitative_abundance, fill = Phylum)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(~ Date) +
    labs(x = "Sampling Location", y = "eDNA gene concentration in seawater (copies/L)") +
    theme_minimal() +
    theme(
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "bold"),  # Bold x-axis values
      axis.text.y = element_text(size = 14, face = "bold"),  # Bold y-axis values
      axis.title.x = element_text(size = 14, face = "bold"),    # Bold axis labels
      axis.title.y = element_text(size = 14, face = "bold"),    # Bold axis labels
      strip.text = element_text(size = 10, face = "bold")   
    ) +
    scale_fill_manual(values = c(setNames(set3_palette, c(top10_asvs, "Others")))))

# Save the plot
ggsave(
  filename = "outputs/metabarcoding_barplots/S16-Quantitative_sp_barplot.png",
  plot = plot,
  width = 16, height = 6,
  dpi = 300
)






########
# -------------- 18S Taxa barplots
########
# Aggregate taxa to the desired level and transform to relative abundances
S18_rel_abun <- S18_taxfilt_data %>%
  transform_sample_counts(function(x) x / sum(x))  # Convert to relative abundances

ps_melted <- psmelt(S18_rel_abun) %>% 
  filter(!is.na(Phylum))

# Top 10 taxa to plot
top10_asvs <- ps_melted %>%
  group_by(Phylum) %>% 
  summarize(total_abundance = sum(Abundance)) %>% # Sum relative abundances
  arrange(desc(total_abundance)) %>%         # Order in descending order
  slice_head(n = 10) %>%                      # Select the top 10 ASVs
  pull(Phylum)                              # Extract sp names



########
# Cooler blanks top 10
########
# Melt the data for ggplot
ps_melted <- psmelt(S18_rel_abun) %>% 
  filter(Location %in% c("Cooler blank 1", "Cooler blank 2")) %>% 
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
         ))) %>% 
  mutate(Phylum = ifelse(Phylum %in% top10_asvs, Phylum, "Others"),
         Phylum = factor(Phylum, levels = c(rev(top10_asvs), "Others")),
         quantitative_abundance = Abundance * X18S.copies.L.SW)


# Create the faceted bar plot
(plot <- ggplot(ps_melted, aes(x = Location, y = Abundance, fill = Phylum)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(~ Date) +
    labs(x = "Sampling Location", y = "Relative Abundance") +
    theme_minimal() +
    theme(
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 14, angle = 45, hjust = 1, face = "bold"),  # Bold x-axis values
      axis.text.y = element_text(size = 14, face = "bold"),  # Bold y-axis values
      axis.title.x = element_text(size = 14, face = "bold"),    # Bold axis labels
      axis.title.y = element_text(size = 14, face = "bold"),    # Bold axis labels
      strip.text = element_text(size = 10, face = "bold")   
    ) +
    scale_fill_manual(values = c(setNames(nature_palette, c(top10_asvs, "Others")))))

# Save the plot
ggsave(
  filename = "outputs/metabarcoding_barplots/S18-BL_rel-abun_sp_barplot.png",
  plot = plot,
  width = 16, height = 6,
  dpi = 300
)



# Quantitative plot
(plot <- ggplot(ps_melted, aes(x = Location, y = quantitative_abundance, fill = Phylum)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(~ Date) +
    labs(x = "Sampling Location", y = "eDNA gene concentration in seawater (copies/L)") +
    theme_minimal() +
    theme(
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 14, angle = 45, hjust = 1, face = "bold"),  # Bold x-axis values
      axis.text.y = element_text(size = 14, face = "bold"),  # Bold y-axis values
      axis.title.x = element_text(size = 14, face = "bold"),    # Bold axis labels
      axis.title.y = element_text(size = 14, face = "bold"),    # Bold axis labels
      strip.text = element_text(size = 10, face = "bold")   
    ) +
    scale_fill_manual(values = c(setNames(nature_palette, c(top10_asvs, "Others")))))

# Save the plot
ggsave(
  filename = "outputs/metabarcoding_barplots/S18-BL_copies_sp_barplot.png",
  plot = plot,
  width = 16, height = 6,
  dpi = 300
)




########
# Environmental genes top 10
########

# Melt the data for ggplot
ps_melted <- psmelt(S18_rel_abun) %>% 
  filter(Depth..m. == 1) %>% 
  filter(Location %in% c("Sackville River", 
                         "DRDC barge", 
                         "McKay bridge", 
                         "Tufts Cove",
                         "McNabs Island")) %>% 
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
                                        "2024-May-16")),
         Phylum = ifelse(Phylum %in% top10_asvs, Phylum, "Others"),
         Phylum = factor(Phylum, levels = c("Others", rev(top10_asvs))),
         quantitative_abundance = Abundance * X18S.copies.L.SW)


# Create the faceted bar plot
(plot <- ggplot(ps_melted, aes(x = Location, y = Abundance, fill = Phylum)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(~ Date) +
    labs(x = "Sampling Location", y = "Relative Abundance") +
    theme_minimal() +
    theme(
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "bold"),  # Bold x-axis values
      axis.text.y = element_text(size = 14, face = "bold"),  # Bold y-axis values
      axis.title.x = element_text(size = 14, face = "bold"),    # Bold axis labels
      axis.title.y = element_text(size = 14, face = "bold"),    # Bold axis labels
      strip.text = element_text(size = 10, face = "bold")   
    ) +
    scale_fill_manual(values = c(setNames(nature_palette, c(top10_asvs, "Others")))))

# Save the plot
ggsave(
  filename = "outputs/metabarcoding_barplots/S18-sp_barplot.png",
  plot = plot,
  width = 16, height = 6,
  dpi = 300
)


# Quantitative plot
(plot <- ggplot(ps_melted, aes(x = Location, y = quantitative_abundance, fill = Phylum)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(~ Date) +
    labs(x = "Sampling Location", y = "eDNA gene concentration in seawater (copies/L)") +
    theme_minimal() +
    theme(
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "bold"),  # Bold x-axis values
      axis.text.y = element_text(size = 14, face = "bold"),  # Bold y-axis values
      axis.title.x = element_text(size = 14, face = "bold"),    # Bold axis labels
      axis.title.y = element_text(size = 14, face = "bold"),    # Bold axis labels
      strip.text = element_text(size = 10, face = "bold")   
    ) +
    scale_fill_manual(values = c(setNames(nature_palette, c(top10_asvs, "Others")))))

# Save the plot
ggsave(
  filename = "outputs/metabarcoding_barplots/S18-Quantitative_sp_barplot.png",
  plot = plot,
  width = 16, height = 6,
  dpi = 300
)




