library(tidyverse)
library(GGally)

#######################
# Importing data
#######################
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


########################
# Correlation matrices
########################
# Transect data matrix 
(plot <- metadata_calc %>%
   ggpairs(columns = c('Date', 'Location', 
                       'Depth (m)', 'DNA (μg/L SW)',
                       '16Sin extraction efficiency'),
           aes(color = Location, alpha = 0.5),
           upper = list(continuous = wrap("cor", size = 5)),
           cardinality_threshold = NULL) +
   scale_color_manual(values = c("#D2B8A6", "#9BA8A8", "#009E73", "#0072B2", "#9932CC", "#D55E00", "#F0E442")) +
   scale_fill_manual(values = c("#D2B8A6", "#9BA8A8", "#009E73", "#0072B2", "#9932CC", "#D55E00", "#F0E442")) +
   theme(strip.text = element_text(size = 14),  # Column titles
         axis.text = element_text(size = 12)    # Axis text
   ))

# Save plot
ggsave(
  filename = "outputs/metadata_correlations/transect_matrix.png",
  plot = plot,
  width = 16, height = 16,
  dpi = 300
)


# 12S
(plot <- metadata_calc %>%
  ggpairs(columns = c('DNA (μg/L SW)', '16Sin extraction efficiency', 
                      '12S copies/L SW', '12S copies/L SW * eff', '12S_DADA2_reads'),
          aes(color = Location, alpha = 0.5),
          upper = list(continuous = wrap("cor", size = 5)),
          cardinality_threshold = NULL) +
  scale_color_manual(values = c("#D2B8A6", "#9BA8A8", "#009E73", "#0072B2", "#9932CC", "#D55E00", "#F0E442")) +
  scale_fill_manual(values = c("#D2B8A6", "#9BA8A8", "#009E73", "#0072B2", "#9932CC", "#D55E00", "#F0E442")) +
  theme(strip.text = element_text(size = 14),  # Column titles
        axis.text = element_text(size = 12)    # Axis text
  ))


# Save plot
ggsave(
  filename = "outputs/metadata_correlations/12S_matrix.png",
  plot = plot,
  width = 16, height = 16,
  dpi = 300
)


# 16S
(plot <- metadata_calc %>%
  ggpairs(columns = c('DNA (μg/L SW)', '16Sin extraction efficiency',
                      '16S copies/L SW', '16S copies/L SW * eff', 
                      '16S seq/L SW'),
          aes(color = Location, alpha = 0.5),
          upper = list(continuous = wrap("cor", size = 5)),
          cardinality_threshold = NULL) +
  scale_color_manual(values = c("#D2B8A6", "#9BA8A8", "#009E73", "#0072B2", "#9932CC", "#D55E00", "#F0E442")) +
  scale_fill_manual(values = c("#D2B8A6", "#9BA8A8", "#009E73", "#0072B2", "#9932CC", "#D55E00", "#F0E442")) +
  theme(strip.text = element_text(size = 14),  # Column titles
        axis.text = element_text(size = 12)    # Axis text
  ))

# Save plot
ggsave(
  filename = "outputs/metadata_correlations/16S_matrix.png",
  plot = plot,
  width = 16, height = 16,
  dpi = 300
)



# 18S
(plot <- metadata_calc %>%
    ggpairs(columns = c('DNA (μg/L SW)', '16Sin extraction efficiency',
                        '18S copies/L SW', '18S copies/L SW * eff', 
                        '18S seq/L SW'),
          aes(color = Location, alpha = 0.5),
          upper = list(continuous = wrap("cor", size = 5)),
          cardinality_threshold = NULL) +
  scale_color_manual(values = c("#D2B8A6", "#9BA8A8", "#009E73", "#0072B2", "#9932CC", "#D55E00", "#F0E442")) +
  scale_fill_manual(values = c("#D2B8A6", "#9BA8A8", "#009E73", "#0072B2", "#9932CC", "#D55E00", "#F0E442")) +
  theme(strip.text = element_text(size = 14),  # Column titles
        axis.text = element_text(size = 12)    # Axis text
  ))

# Save plot
ggsave(
  filename = "outputs/metadata_correlations/18S_matrix.png",
  plot = plot,
  width = 16, height = 16,
  dpi = 300
)



# COI
(plot <- metadata_calc %>%
    ggpairs(columns = c('DNA (μg/L SW)', '16Sin extraction efficiency',
                        'COI copies/L SW', 'COI copies/L SW * eff'),
          aes(color = Location, alpha = 0.5),
          upper = list(continuous = wrap("cor", size = 5)),
          cardinality_threshold = NULL) +
  scale_color_manual(values = c("#D2B8A6", "#9BA8A8", "#009E73", "#0072B2", "#9932CC", "#D55E00", "#F0E442")) +
  scale_fill_manual(values = c("#D2B8A6", "#9BA8A8", "#009E73", "#0072B2", "#9932CC", "#D55E00", "#F0E442")) +
  theme(strip.text = element_text(size = 14),  # Column titles
        axis.text = element_text(size = 12)    # Axis text
  ))

# Save plot
ggsave(
  filename = "outputs/metadata_correlations/COI_matrix.png",
  plot = plot,
  width = 16, height = 16,
  dpi = 300
)





# Cross-marker matrix
(plot <- metadata_calc %>%
    ggpairs(columns = c('12S copies/L SW','16S copies/L SW', 
                        '18S copies/L SW', 'COI copies/L SW'),
            aes(color = Location, alpha = 0.5),
            upper = list(continuous = wrap("cor", size = 5)),
            cardinality_threshold = NULL) +
    scale_color_manual(values = c("#D2B8A6", "#9BA8A8", "#009E73", "#0072B2", "#9932CC", "#D55E00", "#F0E442")) +
    scale_fill_manual(values = c("#D2B8A6", "#9BA8A8", "#009E73", "#0072B2", "#9932CC", "#D55E00", "#F0E442")) +
    theme(strip.text = element_text(size = 14),  # Column titles
          axis.text = element_text(size = 12)    # Axis text
    ))

# Save plot
ggsave(
  filename = "outputs/metadata_correlations/markers_matrix.png",
  plot = plot,
  width = 16, height = 16,
  dpi = 300
)

