library(tidyverse)

### Importing dataset
metadata_calc <-  read_csv("data/metadata_qPCR_calculated.csv") %>% 
  mutate(Location = factor(Location, levels = c("Cooler blank 1",
                                         "Cooler blank 2",
                                         "Sackville River", 
                                         "DRDC barge", 
                                         "McKay bridge", 
                                         "Tufts Cove",
                                         "McNabs Island"))
  )
metadata_calc$Date <- as.Date(metadata_calc$Date)




##########
## Volume filt, [DNA], 16S internal qPCR (raw), Extr. efficiency coefficient
##########
long_DNA_data <- metadata_calc %>%
  pivot_longer(cols = c("Volume filtered (ml SW)", 
                        "DNA (μg/L SW)"), 
               names_to = "value", 
               values_to = "quantities") %>% 
  mutate(value = factor(value, levels = c("Volume filtered (ml SW)", 
                                          "DNA (μg/L SW)")))

# Plot surface data
(plot <- long_DNA_data %>%   
  filter(`Depth (m)` < 2) %>% 
  ggplot(aes(x = Date, y = quantities, color = Location, group = Location)) +
  geom_line(linewidth = 1) +
  #geom_smooth(method = "loess", se = FALSE) +  # Add LOESS smoother without CI
  facet_wrap(~value, scales = "free_y") +
  labs(x = "Sampling date",
       y = "Value",
       color = "Sampling location") +
  theme_minimal() +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  scale_color_manual(values = c("#F0B8B2", "#C1CDCD", "#009E73", "#0072B2", "#9932CC", "#D55E00", "#F0E442")) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1, face = "bold"),  # Bold x-axis values
    axis.text.y = element_text(size = 16, face = "bold"),  # Bold y-axis values
    axis.title.x = element_text(size = 16, face = "bold"),    # Bold axis labels
    axis.title.y = element_text(size = 16, face = "bold"),    # Bold axis labels
    strip.text = element_text(size = 16, face = "bold")    # Bold facet titles
  ))

ggsave(
  filename = "outputs/metadata_TimeSeries/Surface_vol-dna-eff.png",
  plot = plot, width = 12, height = 6, dpi = 300
)

# Plot deep data
(plot <- long_DNA_data %>%   
  filter(`Depth (m)` > 2 | `Depth (m)` == 0) %>% 
  ggplot(aes(x = Date, y = quantities, color = Location, group = Location)) +
  geom_line(linewidth = 1) +
  #geom_smooth(method = "loess", se = FALSE) +  # Add LOESS smoother without CI
  facet_wrap(~value, scales = "free_y") +
  labs(x = "Date", 
       y = "Value",
       color = "Sampling location") +
  theme_minimal() +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  scale_color_manual(values = c("#F0B8B2", "#C1CDCD", "#009E73", "#0072B2", "#9932CC", "#D55E00", "#F0E442")) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 16, angle = 45, hjust = 1, face = "bold"),  # Bold x-axis values
      axis.text.y = element_text(size = 16, face = "bold"),  # Bold y-axis values
      axis.title.x = element_text(size = 16, face = "bold"),    # Bold axis labels
      axis.title.y = element_text(size = 16, face = "bold"),    # Bold axis labels
      strip.text = element_text(size = 16, face = "bold")    # Bold facet titles
    ))
ggsave(
  filename = "outputs/metadata_TimeSeries/Deep_vol-dna-eff.png",
  plot = plot, width = 12, height = 6, dpi = 300
)





##########
## qPCR quantities, no efficiency calculations
##########
long_qPCR_data <- metadata_calc %>%
  pivot_longer(cols = c("12S copies/L SW", 
                        "16S copies/L SW",
                        "18S copies/L SW",
                        "COI copies/L SW"), 
               names_to = "gene", 
               values_to = "gene_quantity")


# Plot surface data
(plot <- long_qPCR_data %>% 
  filter(`Depth (m)` < 2) %>% 
  ggplot(aes(x = Date, y = gene_quantity, color = Location, group = Location)) +
  geom_line(size = 1) +
  #geom_smooth(method = "loess", se = FALSE) +  # Add LOESS smoother without CI
  facet_wrap(~gene, scales = "free_y") +
  labs(x = "Sampling date", 
       y = "Gene Concentrations in Seawater",
       color = "Sampling location") +
  theme_minimal() +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  scale_y_continuous(labels = scales::scientific) +  # Force scientific notation on y-axis
  scale_color_manual(values = c("#F0B8B2", "#C1CDCD", "#009E73", "#0072B2", "#9932CC", "#D55E00", "#F0E442")) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  # Bold x-axis values
      axis.text.y = element_text(face = "bold"),  # Bold y-axis values
      axis.title.x = element_text(face = "bold"),    # Bold axis labels
      axis.title.y = element_text(face = "bold"),    # Bold axis labels
      strip.text = element_text(face = "bold")    # Bold facet titles
    ))

ggsave(
  filename = "outputs/metadata_TimeSeries/surf_qPCR_quants.png",
  plot = plot, width = 12, height = 6, dpi = 300
)

(plot <- long_qPCR_data %>% 
  filter(`Depth (m)` > 2 | `Depth (m)` == 0) %>% 
  ggplot(aes(x = Date, y = gene_quantity, color = Location, group = Location)) +
  geom_line(size = 1) +
  facet_wrap(~gene, scales = "free_y") +
  labs(x = "Sampling date", 
       y = "Gene Concentrations in Seawater",
       color = "Sampling location") +
  theme_minimal() +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  scale_y_continuous(labels = scales::scientific) +  # Force scientific notation on y-axis
  scale_color_manual(values = c("#F0B8B2", "#C1CDCD", "#009E73", "#0072B2", "#9932CC", "#D55E00", "#F0E442")) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  # Bold x-axis values
      axis.text.y = element_text(face = "bold"),  # Bold y-axis values
      axis.title.x = element_text(face = "bold"),    # Bold axis labels
      axis.title.y = element_text(face = "bold"),    # Bold axis labels
      strip.text = element_text(face = "bold")    # Bold facet titles
    ))

ggsave(
  filename = "outputs/metadata_TimeSeries/deep_qPCR_quants.png",
  plot = plot, width = 12, height = 6, dpi = 300
)


##########
## qPCR quantities, with efficiency calculations
##########
long_qPCR_eff_data <- metadata_calc %>%
  pivot_longer(cols = c("12S copies/L SW * eff", 
                        "16S copies/L SW * eff",
                        "18S copies/L SW * eff",
                        "COI copies/L SW * eff"), 
               names_to = "gene", 
               values_to = "gene_quantity")


# Plot surface data
(plot <- long_qPCR_eff_data %>% 
  filter(`Depth (m)` < 2) %>% 
  ggplot(aes(x = Date, y = gene_quantity, color = Location, group = Location)) +
  geom_line(size = 1) +
  #geom_smooth(method = "loess", se = FALSE) +  # Add LOESS smoother without CI
  facet_wrap(~gene, scales = "free_y") +
  labs(x = "Sampling date", 
       y = "Gene Concentrations in Seawater",
       color = "Sampling location") +
  theme_minimal() +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  scale_color_manual(values = c("#F0B8B2", "#C1CDCD", "#009E73", "#0072B2", "#9932CC", "#D55E00", "#F0E442")) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  # Bold x-axis values
    axis.text.y = element_text(face = "bold"),  # Bold y-axis values
    axis.title.x = element_text(face = "bold"),    # Bold axis labels
    axis.title.y = element_text(face = "bold"),    # Bold axis labels
    strip.text = element_text(face = "bold")    # Bold facet titles
  ))

ggsave(
  filename = "outputs/metadata_TimeSeries/surface_qPCR_quants_w_eff.png",
  plot = plot, width = 12, height = 6, dpi = 300
)



(plot <- long_qPCR_eff_data %>% 
  filter(`Depth (m)` > 2 | `Depth (m)` == 0) %>% 
  ggplot(aes(x = Date, y = gene_quantity, color = Location, group = Location)) +
  geom_line(size = 1) +
  facet_wrap(~gene, scales = "free_y") +
  labs(x = "Sampling date", 
       y = "Gene Concentrations in Seawater",
       color = "Sampling location") +
  theme_minimal() +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  scale_color_manual(values = c("#F0B8B2", "#C1CDCD", "#009E73", "#0072B2", "#9932CC", "#D55E00", "#F0E442")) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  # Bold x-axis values
    axis.text.y = element_text(face = "bold"),  # Bold y-axis values
    axis.title.x = element_text(face = "bold"),    # Bold axis labels
    axis.title.y = element_text(face = "bold"),    # Bold axis labels
    strip.text = element_text(face = "bold")    # Bold facet titles
  ))

ggsave(
  filename = "outputs/metadata_TimeSeries/deep_qPCR_quants_w_eff.png",
  plot = plot, width = 12, height = 6, dpi = 300
)




##########
## Spike to total reads ratios to transform total read counts
##########
long_seq_data <- metadata_calc %>%
  pivot_longer(cols = c("16S_DADA2_reads",
                        "18S_DADA2_reads",
                        "16S_spike_ratios",
                        "18S_spike_ratios",
                        "16S seq/L SW",
                        "18S seq/L SW"),
               names_to = "gene", 
               values_to = "reads") %>% 
  mutate(gene = factor(gene, levels = c("16S_DADA2_reads",
                                        "18S_DADA2_reads",
                                        "16S_spike_ratios",
                                        "18S_spike_ratios",
                                        "16S seq/L SW",
                                        "18S seq/L SW")))


(plot <- long_seq_data %>% 
  filter(`Depth (m)` == 1) %>% 
  ggplot(aes(x = Date, y = reads, color = Location, group = Location)) +
  geom_line(size = 1) +
  #geom_smooth(method = "loess", se = FALSE) +  # Add LOESS smoother without CI
  facet_wrap(~gene, scales = "free_y", ncol = 2) +
  labs(x = "Sampling date", 
       y = "Sequencing Values",
       color = "Sampling location") +
  theme_minimal() +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  scale_color_manual(values = c("#009E73", "#0072B2", "#9932CC", "#D55E00", "#F0E442")) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  # Bold x-axis values
    axis.text.y = element_text(face = "bold"),  # Bold y-axis values
    axis.title.x = element_text(face = "bold"),    # Bold axis labels
    axis.title.y = element_text(face = "bold"),    # Bold axis labels
    strip.text = element_text(face = "bold")    # Bold facet titles
  ))

ggsave(
  filename = "outputs/metadata_TimeSeries/surface_SEQ_quants.png",
  plot = plot, width = 12, height = 10, dpi = 300
)


(plot <- long_seq_data %>% 
  filter(`Depth (m)` > 2 ) %>% 
  ggplot(aes(x = Date, y = reads, color = Location, group = Location)) +
  geom_line(size = 1) +
  #geom_smooth(method = "loess", se = FALSE) +  # Add LOESS smoother without CI
  facet_wrap(~gene, scales = "free_y", ncol = 2) +
  labs(x = "Sampling date", 
       y = "Sequencing Values",
       color = "Sampling location") +
  theme_minimal() +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  scale_color_manual(values = c("#009E73", "#0072B2", "#9932CC", "#D55E00", "#F0E442")) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  # Bold x-axis values
    axis.text.y = element_text(face = "bold"),  # Bold y-axis values
    axis.title.x = element_text(face = "bold"),    # Bold axis labels
    axis.title.y = element_text(face = "bold"),    # Bold axis labels
    strip.text = element_text(face = "bold")    # Bold facet titles
  ))


ggsave(
  filename = "outputs/metadata_TimeSeries/deep_SEQ_quants.png",
  plot = plot, width = 12, height = 10, dpi = 300
)










########################################
#### PER GENE PLOTTING
########################################

###########
## 12S
###########
long_12S_data <- metadata_calc %>%
  pivot_longer(cols = c("DNA (μg/L SW)",
                        "16Sin extraction efficiency",
                        "12S copies/L SW",
                        "12S copies/L SW * eff"), 
               names_to = "gene", 
               values_to = "gene_quantity") %>% 
  mutate(gene = factor(gene, levels = c("DNA (μg/L SW)",
                                        "16Sin extraction efficiency",
                                        "12S copies/L SW",
                                        "12S copies/L SW * eff")))


# Plot surface data
(plot <- long_12S_data %>% 
  filter(`Depth (m)` < 2) %>% 
  ggplot(aes(x = Date, y = gene_quantity, color = Location, group = Location)) +
  geom_line(size = 1) +
  #geom_smooth(method = "loess", se = FALSE) +  # Add LOESS smoother without CI
  facet_wrap(~gene, scales = "free_y") +
  labs(x = "Sampling date", 
       y = "Gene Quantity",
       color = "Sampling location") +
  theme_minimal() +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  scale_color_manual(values = c("#F0B8B2", "#C1CDCD", "#009E73", "#0072B2", "#9932CC", "#D55E00", "#F0E442")) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  # Bold x-axis values
    axis.text.y = element_text(face = "bold"),  # Bold y-axis values
    axis.title.x = element_text(face = "bold"),    # Bold axis labels
    axis.title.y = element_text(face = "bold"),    # Bold axis labels
    strip.text = element_text(face = "bold")    # Bold facet titles
  ))


ggsave(
  filename = "outputs/metadata_TimeSeries/surface_12S_all_metrics.png",
  plot = plot, width = 12, height = 6, dpi = 300
)


(plot <- long_12S_data %>% 
  filter(`Depth (m)` > 2 | `Depth (m)` == 0) %>% 
  ggplot(aes(x = Date, y = gene_quantity, color = Location, group = Location)) +
  geom_line(size = 1) +
  facet_wrap(~gene, scales = "free_y") +
  labs(x = "Sampling date", 
       y = "Gene Quantity",
       color = "Sampling location") +
  theme_minimal() +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  scale_color_manual(values = c("#F0B8B2", "#C1CDCD", "#009E73", "#0072B2", "#9932CC", "#D55E00", "#F0E442")) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  # Bold x-axis values
    axis.text.y = element_text(face = "bold"),  # Bold y-axis values
    axis.title.x = element_text(face = "bold"),    # Bold axis labels
    axis.title.y = element_text(face = "bold"),    # Bold axis labels
    strip.text = element_text(face = "bold")    # Bold facet titles
  ))

ggsave(
  filename = "outputs/metadata_TimeSeries/deep_12S_all_metrics.png",
  plot = plot, width = 12, height = 6, dpi = 300
)


###########
## 16S
###########
long_16S_data <- metadata_calc %>%
  pivot_longer(cols = c("DNA (μg/L SW)",
                        "16Sin extraction efficiency",
                        "16S copies/L SW",
                        "16S copies/L SW * eff",                        
                        "16S seq/L SW"),
               names_to = "gene", 
               values_to = "gene_quantity") %>% 
  mutate(gene = factor(gene, levels = c("DNA (μg/L SW)",
                                        "16Sin extraction efficiency",
                                        "16S copies/L SW",
                                        "16S copies/L SW * eff",                        
                                        "16S seq/L SW")))

# Plot surface data
(plot <- long_16S_data %>% 
  filter(`Depth (m)` < 2) %>% 
  ggplot(aes(x = Date, y = gene_quantity, color = Location, group = Location)) +
  geom_line(size = 1) +
  #geom_smooth(method = "loess", se = FALSE) +  # Add LOESS smoother without CI
  facet_wrap(~gene, scales = "free_y", ncol = 2) +
  labs(x = "Sampling date", 
       y = "Gene Quantity",
       color = "Sampling location") +
  theme_minimal() +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  scale_color_manual(values = c("#F0B8B2", "#C1CDCD", "#009E73", "#0072B2", "#9932CC", "#D55E00", "#F0E442")) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  # Bold x-axis values
      axis.text.y = element_text(face = "bold"),  # Bold y-axis values
      axis.title.x = element_text(face = "bold"),    # Bold axis labels
      axis.title.y = element_text(face = "bold"),    # Bold axis labels
      strip.text = element_text(face = "bold")    # Bold facet titles
    ))

ggsave(
  filename = "outputs/metadata_TimeSeries/surface_16S_all_metrics.png",
  plot = plot, width = 12, height = 10, dpi = 300
)


(plot <- long_16S_data %>% 
  filter(`Depth (m)` > 2 | `Depth (m)` == 0) %>% 
  ggplot(aes(x = Date, y = gene_quantity, color = Location, group = Location)) +
  geom_line(size = 1) +
  facet_wrap(~gene, scales = "free_y", ncol = 2) +
  labs(x = "Sampling date", 
       y = "Gene Quantity",
       color = "Sampling location") +
  theme_minimal() +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  scale_color_manual(values = c("#F0B8B2", "#C1CDCD", "#009E73", "#0072B2", "#9932CC", "#D55E00", "#F0E442")) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  # Bold x-axis values
    axis.text.y = element_text(face = "bold"),  # Bold y-axis values
    axis.title.x = element_text(face = "bold"),    # Bold axis labels
    axis.title.y = element_text(face = "bold"),    # Bold axis labels
    strip.text = element_text(face = "bold")    # Bold facet titles
  ))

ggsave(
  filename = "outputs/metadata_TimeSeries/deep_16S_all_metrics.png",
  plot = plot, width = 12, height = 10, dpi = 300
)



###########
## 18S
###########
long_18S_data <- metadata_calc %>%
  pivot_longer(cols = c("DNA (μg/L SW)",
                        "16Sin extraction efficiency",
                        "18S copies/L SW",
                        "18S copies/L SW * eff",                        
                        "18S seq/L SW"), 
               names_to = "gene", 
               values_to = "gene_quantity") %>% 
  mutate(gene = factor(gene, levels = c("DNA (μg/L SW)",
                                        "16Sin extraction efficiency",
                                        "18S copies/L SW",
                                        "18S copies/L SW * eff",                        
                                        "18S seq/L SW")))


# Plot surface data
(plot <- long_18S_data %>% 
  filter(`Depth (m)` < 2) %>% 
  ggplot(aes(x = Date, y = gene_quantity, color = Location, group = Location)) +
  geom_line(size = 1) +
  #geom_smooth(method = "loess", se = FALSE) +  # Add LOESS smoother without CI
  facet_wrap(~gene, scales = "free_y", ncol = 2) +
  labs(x = "Sampling date", 
       y = "Gene Quantity",
       color = "Sampling location") +
  theme_minimal() +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  scale_color_manual(values = c("#F0B8B2", "#C1CDCD", "#009E73", "#0072B2", "#9932CC", "#D55E00", "#F0E442")) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  # Bold x-axis values
      axis.text.y = element_text(face = "bold"),  # Bold y-axis values
      axis.title.x = element_text(face = "bold"),    # Bold axis labels
      axis.title.y = element_text(face = "bold"),    # Bold axis labels
      strip.text = element_text(face = "bold")    # Bold facet titles
    ))

ggsave(
  filename = "outputs/metadata_TimeSeries/surface_18S_all_metrics.png",
  plot = plot, width = 12, height = 10, dpi = 300
)


(plot <- long_18S_data %>% 
  filter(`Depth (m)` > 2 | `Depth (m)` == 0) %>% 
  ggplot(aes(x = Date, y = gene_quantity, color = Location, group = Location)) +
  geom_line(size = 1) +
  facet_wrap(~gene, scales = "free_y", ncol = 2) +
  labs(x = "Sampling date", 
       y = "Gene Quantity",
       color = "Sampling location") +
  theme_minimal() +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  scale_color_manual(values = c("#F0B8B2", "#C1CDCD", "#009E73", "#0072B2", "#9932CC", "#D55E00", "#F0E442")) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  # Bold x-axis values
      axis.text.y = element_text(face = "bold"),  # Bold y-axis values
      axis.title.x = element_text(face = "bold"),    # Bold axis labels
      axis.title.y = element_text(face = "bold"),    # Bold axis labels
      strip.text = element_text(face = "bold")    # Bold facet titles
    ))


ggsave(
  filename = "outputs/metadata_TimeSeries/deep_18S_all_metrics.png",
  plot = plot, width = 12, height = 10, dpi = 300
)



###########
## COI
###########
long_COI_data <- metadata_calc %>%
  pivot_longer(cols = c("DNA (μg/L SW)",
                        "16Sin extraction efficiency",
                        "COI copies/L SW",
                        "COI copies/L SW * eff"), 
               names_to = "gene", 
               values_to = "gene_quantity") %>% 
  mutate(gene = factor(gene, levels = c("DNA (μg/L SW)",
                                        "16Sin extraction efficiency",
                                        "COI copies/L SW",
                                        "COI copies/L SW * eff")))

# Plot surface data
(plot <- long_COI_data %>% 
  filter(`Depth (m)` < 2) %>% 
  ggplot(aes(x = Date, y = gene_quantity, color = Location, group = Location)) +
  geom_line(size = 1) +
  #geom_smooth(method = "loess", se = FALSE) +  # Add LOESS smoother without CI
  facet_wrap(~gene, scales = "free_y") +
  labs(x = "Sampling date", 
       y = "Gene Quantity",
       color = "Sampling location") +
  theme_minimal() +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  scale_color_manual(values = c("#F0B8B2", "#C1CDCD", "#009E73", "#0072B2", "#9932CC", "#D55E00", "#F0E442")) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  # Bold x-axis values
      axis.text.y = element_text(face = "bold"),  # Bold y-axis values
      axis.title.x = element_text(face = "bold"),    # Bold axis labels
      axis.title.y = element_text(face = "bold"),    # Bold axis labels
      strip.text = element_text(face = "bold")    # Bold facet titles
    ))


ggsave(
  filename = "outputs/metadata_TimeSeries/surface_COI_all_metrics.png",
  plot = plot, width = 12, height = 6, dpi = 300
)

(plot <- long_COI_data %>% 
  filter(`Depth (m)` > 2 | `Depth (m)` == 0) %>% 
  ggplot(aes(x = Date, y = gene_quantity, color = Location, group = Location)) +
  geom_line(size = 1) +
  facet_wrap(~gene, scales = "free_y") +
  labs(x = "Sampling date", 
       y = "Gene Quantity",
       color = "Sampling location") +
  theme_minimal() +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  scale_color_manual(values = c("#F0B8B2", "#C1CDCD", "#009E73", "#0072B2", "#9932CC", "#D55E00", "#F0E442")) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  # Bold x-axis values
      axis.text.y = element_text(face = "bold"),  # Bold y-axis values
      axis.title.x = element_text(face = "bold"),    # Bold axis labels
      axis.title.y = element_text(face = "bold"),    # Bold axis labels
      strip.text = element_text(face = "bold")    # Bold facet titles
    ))

ggsave(
  filename = "outputs/metadata_TimeSeries/deep_COI_all_metrics.png",
  plot = plot, width = 12, height = 6, dpi = 300
)
