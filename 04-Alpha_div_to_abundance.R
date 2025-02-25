library(phyloseq)
library(tidyverse)


# 12S Surface 
chao1_values <- estimate_richness(S12_taxfilt_data, measures = "Chao1")
rownames(chao1_values) <- gsub("\\.", "-", rownames(chao1_values))
chao1_values$SampleID <- rownames(chao1_values)

sample_metadata <- data.frame(sample_data(S12_taxfilt_data))
sample_metadata$SampleID <- rownames(sample_metadata)

sample_metadata <- merge(sample_metadata, chao1_values, by = "SampleID", all.x = TRUE) %>% 
  as.data.frame()

  
summary_data <- sample_metadata %>%
  filter(Depth..m. == 1) %>%
  group_by(Date) %>%
  summarise(
    mean_12S = mean(`X12S.copies.L.SW`, na.rm = TRUE),
    sd_12S = sd(`X12S.copies.L.SW`, na.rm = TRUE),
    n_12S = sum(!is.na(`X12S.copies.L.SW`)),  # Count non-missing values
    mean_Chao1 = mean(Chao1, na.rm = TRUE),
    sd_Chao1 = sd(Chao1, na.rm = TRUE),
    n_Chao1 = sum(!is.na(Chao1))  # Count non-missing values
  ) %>%
  mutate(
    se_12S = sd_12S / sqrt(n_12S),
    se_Chao1 = sd_Chao1 / sqrt(n_Chao1),
    ci_12S_lower = mean_12S - 1.96 * se_12S,
    ci_12S_upper = mean_12S + 1.96 * se_12S,
    ci_Chao1_lower = mean_Chao1 - 1.96 * se_Chao1,
    ci_Chao1_upper = mean_Chao1 + 1.96 * se_Chao1
  )

# Normalize for dual axis visualization
summary_data <- summary_data %>%
  mutate(
    norm_12S = scale(mean_12S, center = TRUE, scale = TRUE)[,1],
    norm_Chao1 = scale(mean_Chao1, center = TRUE, scale = TRUE)[,1],
    norm_ci_12S_lower = scale(ci_12S_lower, center = TRUE, scale = TRUE)[,1],
    norm_ci_12S_upper = scale(ci_12S_upper, center = TRUE, scale = TRUE)[,1],
    norm_ci_Chao1_lower = scale(ci_Chao1_lower, center = TRUE, scale = TRUE)[,1],
    norm_ci_Chao1_upper = scale(ci_Chao1_upper, center = TRUE, scale = TRUE)[,1]
  )


# Plot with confidence intervals
(plot <- ggplot(summary_data, aes(x = Date)) +
  
  # 12S copies line + confidence interval shading
  geom_ribbon(aes(ymin = norm_ci_12S_lower, ymax = norm_ci_12S_upper, fill = " "), alpha = 0.2) +
  geom_line(aes(y = norm_12S, color = "12S copies/Litre SW"), size = 1) +
  
  # Chao1 richness line + confidence interval shading
  geom_ribbon(aes(ymin = norm_ci_Chao1_lower, ymax = norm_ci_Chao1_upper, fill = "95% Confidence intervals"), alpha = 0.2) +
  geom_line(aes(y = norm_Chao1, color = "Chao1 richness"), size = 1, linetype = "dashed") +
  
  # Set axis labels and dual y-axis
  scale_y_continuous(
    name = "12S copies/Litre SW", 
    sec.axis = sec_axis(~ ., name = "Chao1 richness")  # Right y-axis
  ) +
  
  # Customize colors and theme
  scale_color_manual(values = c("12S copies/Litre SW" = "#4A4E9E", "Chao1 richness" = "#006D40")) +
  scale_fill_manual(values = c(" " = "#4A4E9E", "95% Confidence intervals" = "#006D40")) +
  
  labs(x = "Date", color = "Gene metric:") +
  theme_minimal() +
  theme(legend.position = "none",
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "bold"),  # Bold x-axis values
    axis.text.y = element_text(size = 14, face = "bold"),  # Bold y-axis values
    axis.title.x = element_text(size = 14, face = "bold"),    # Bold axis labels
    axis.title.y = element_text(size = 14, face = "bold")    # Bold axis labels
  ))


# Save the plot
ggsave(
  filename = "outputs/alpha/S12-surf_alpha_counts.png",
  plot = plot,
  width = 6, height = 4,
  dpi = 300
)


# 12S Deep

summary_data <- sample_metadata %>%
  filter(Depth..m. > 1) %>%
  group_by(Date) %>%
  summarise(
    mean_12S = mean(`X12S.copies.L.SW`, na.rm = TRUE),
    sd_12S = sd(`X12S.copies.L.SW`, na.rm = TRUE),
    n_12S = sum(!is.na(`X12S.copies.L.SW`)),  # Count non-missing values
    mean_Chao1 = mean(Chao1, na.rm = TRUE),
    sd_Chao1 = sd(Chao1, na.rm = TRUE),
    n_Chao1 = sum(!is.na(Chao1))  # Count non-missing values
  ) %>%
  mutate(
    se_12S = sd_12S / sqrt(n_12S),
    se_Chao1 = sd_Chao1 / sqrt(n_Chao1),
    ci_12S_lower = mean_12S - 1.96 * se_12S,
    ci_12S_upper = mean_12S + 1.96 * se_12S,
    ci_Chao1_lower = mean_Chao1 - 1.96 * se_Chao1,
    ci_Chao1_upper = mean_Chao1 + 1.96 * se_Chao1
  )

# Normalize for dual axis visualization
summary_data <- summary_data %>%
  mutate(
    norm_12S = scale(mean_12S, center = TRUE, scale = TRUE)[,1],
    norm_Chao1 = scale(mean_Chao1, center = TRUE, scale = TRUE)[,1],
    norm_ci_12S_lower = scale(ci_12S_lower, center = TRUE, scale = TRUE)[,1],
    norm_ci_12S_upper = scale(ci_12S_upper, center = TRUE, scale = TRUE)[,1],
    norm_ci_Chao1_lower = scale(ci_Chao1_lower, center = TRUE, scale = TRUE)[,1],
    norm_ci_Chao1_upper = scale(ci_Chao1_upper, center = TRUE, scale = TRUE)[,1]
  )


# Plot with confidence intervals
(plot <- ggplot(summary_data, aes(x = Date)) +
    
    # 12S copies line + confidence interval shading
    geom_ribbon(aes(ymin = norm_ci_12S_lower, ymax = norm_ci_12S_upper, fill = " "), alpha = 0.2) +
    geom_line(aes(y = norm_12S, color = "12S copies/Litre SW"), size = 1) +
    
    # Chao1 richness line + confidence interval shading
    geom_ribbon(aes(ymin = norm_ci_Chao1_lower, ymax = norm_ci_Chao1_upper, fill = "95% Confidence intervals"), alpha = 0.2) +
    geom_line(aes(y = norm_Chao1, color = "Chao1 richness"), size = 1, linetype = "dashed") +
    
    # Set axis labels and dual y-axis
    scale_y_continuous(
      name = "12S copies/Litre SW", 
      sec.axis = sec_axis(~ ., name = "Chao1 richness")  # Right y-axis
    ) +
    
    # Customize colors and theme
    scale_color_manual(values = c("12S copies/Litre SW" = "#4A4E9E", "Chao1 richness" = "#006D40")) +
    scale_fill_manual(values = c(" " = "#4A4E9E", "95% Confidence intervals" = "#006D40")) +
    
    labs(x = "Date", color = "Gene metric:") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "bold"),  # Bold x-axis values
          axis.text.y = element_text(size = 14, face = "bold"),  # Bold y-axis values
          axis.title.x = element_text(size = 14, face = "bold"),    # Bold axis labels
          axis.title.y = element_text(size = 14, face = "bold")    # Bold axis labels
    ))


# Save the plot
ggsave(
  filename = "outputs/alpha/S12-deep_alpha_counts.png",
  plot = plot,
  width = 6, height = 4,
  dpi = 300
)









# 16S Surface 
chao1_values <- estimate_richness(S16_taxfilt_data, measures = "Chao1")
rownames(chao1_values) <- gsub("\\.", "-", rownames(chao1_values))
chao1_values$SampleID <- rownames(chao1_values)

sample_metadata <- data.frame(sample_data(S16_taxfilt_data))
sample_metadata$SampleID <- rownames(sample_metadata)

sample_metadata <- merge(sample_metadata, chao1_values, by = "SampleID", all.x = TRUE) %>% 
  as.data.frame()


summary_data <- sample_metadata %>%
  filter(Depth..m. == 1) %>%
  group_by(Date) %>%
  summarise(
    mean_16S = mean(`X16S.copies.L.SW`, na.rm = TRUE),
    sd_16S = sd(`X16S.copies.L.SW`, na.rm = TRUE),
    n_16S = sum(!is.na(`X16S.copies.L.SW`)),  # Count non-missing values
    mean_Chao1 = mean(Chao1, na.rm = TRUE),
    sd_Chao1 = sd(Chao1, na.rm = TRUE),
    n_Chao1 = sum(!is.na(Chao1))  # Count non-missing values
  ) %>%
  mutate(
    se_16S = sd_16S / sqrt(n_16S),
    se_Chao1 = sd_Chao1 / sqrt(n_Chao1),
    ci_16S_lower = mean_16S - 1.96 * se_16S,
    ci_16S_upper = mean_16S + 1.96 * se_16S,
    ci_Chao1_lower = mean_Chao1 - 1.96 * se_Chao1,
    ci_Chao1_upper = mean_Chao1 + 1.96 * se_Chao1
  )

# Normalize for dual axis visualization
summary_data <- summary_data %>%
  mutate(
    norm_16S = scale(mean_16S, center = TRUE, scale = TRUE)[,1],
    norm_Chao1 = scale(mean_Chao1, center = TRUE, scale = TRUE)[,1],
    norm_ci_16S_lower = scale(ci_16S_lower, center = TRUE, scale = TRUE)[,1],
    norm_ci_16S_upper = scale(ci_16S_upper, center = TRUE, scale = TRUE)[,1],
    norm_ci_Chao1_lower = scale(ci_Chao1_lower, center = TRUE, scale = TRUE)[,1],
    norm_ci_Chao1_upper = scale(ci_Chao1_upper, center = TRUE, scale = TRUE)[,1]
  )


# Plot with confidence intervals
(plot <- ggplot(summary_data, aes(x = Date)) +
    
    # 12S copies line + confidence interval shading
    geom_ribbon(aes(ymin = norm_ci_16S_lower, ymax = norm_ci_16S_upper, fill = " "), alpha = 0.2) +
    geom_line(aes(y = norm_16S, color = "16S copies/Litre SW"), size = 1) +
    
    # Chao1 richness line + confidence interval shading
    geom_ribbon(aes(ymin = norm_ci_Chao1_lower, ymax = norm_ci_Chao1_upper, fill = "95% Confidence intervals"), alpha = 0.2) +
    geom_line(aes(y = norm_Chao1, color = "Chao1 richness"), size = 1, linetype = "dashed") +
    
    # Set axis labels and dual y-axis
    scale_y_continuous(
      name = "16S copies/Litre SW", 
      sec.axis = sec_axis(~ ., name = "Chao1 richness")  # Right y-axis
    ) +
    
    # Customize colors and theme
    scale_color_manual(values = c("16S copies/Litre SW" = "#4A4E9E", "Chao1 richness" = "#006D40")) +
    scale_fill_manual(values = c(" " = "#4A4E9E", "95% Confidence intervals" = "#006D40")) +
    
    labs(x = "Date", color = "Gene metric:") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "bold"),  # Bold x-axis values
          axis.text.y = element_text(size = 14, face = "bold"),  # Bold y-axis values
          axis.title.x = element_text(size = 14, face = "bold"),    # Bold axis labels
          axis.title.y = element_text(size = 14, face = "bold")    # Bold axis labels
    ))

# Save the plot
ggsave(
  filename = "outputs/alpha/S16-surf_alpha_counts.png",
  plot = plot,
  width = 6, height = 4,
  dpi = 300
)


# 16S Deep

summary_data <- sample_metadata %>%
  filter(Depth..m. > 1) %>%
  group_by(Date) %>%
  summarise(
    mean_16S = mean(`X16S.copies.L.SW`, na.rm = TRUE),
    sd_16S = sd(`X16S.copies.L.SW`, na.rm = TRUE),
    n_16S = sum(!is.na(`X16S.copies.L.SW`)),  # Count non-missing values
    mean_Chao1 = mean(Chao1, na.rm = TRUE),
    sd_Chao1 = sd(Chao1, na.rm = TRUE),
    n_Chao1 = sum(!is.na(Chao1))  # Count non-missing values
  ) %>%
  mutate(
    se_16S = sd_16S / sqrt(n_16S),
    se_Chao1 = sd_Chao1 / sqrt(n_Chao1),
    ci_16S_lower = mean_16S - 1.96 * se_16S,
    ci_16S_upper = mean_16S + 1.96 * se_16S,
    ci_Chao1_lower = mean_Chao1 - 1.96 * se_Chao1,
    ci_Chao1_upper = mean_Chao1 + 1.96 * se_Chao1
  )

# Normalize for dual axis visualization
summary_data <- summary_data %>%
  mutate(
    norm_16S = scale(mean_16S, center = TRUE, scale = TRUE)[,1],
    norm_Chao1 = scale(mean_Chao1, center = TRUE, scale = TRUE)[,1],
    norm_ci_16S_lower = scale(ci_16S_lower, center = TRUE, scale = TRUE)[,1],
    norm_ci_16S_upper = scale(ci_16S_upper, center = TRUE, scale = TRUE)[,1],
    norm_ci_Chao1_lower = scale(ci_Chao1_lower, center = TRUE, scale = TRUE)[,1],
    norm_ci_Chao1_upper = scale(ci_Chao1_upper, center = TRUE, scale = TRUE)[,1]
  )


# Plot with confidence intervals
(plot <- ggplot(summary_data, aes(x = Date)) +
    
    # 12S copies line + confidence interval shading
    geom_ribbon(aes(ymin = norm_ci_16S_lower, ymax = norm_ci_16S_upper, fill = " "), alpha = 0.2) +
    geom_line(aes(y = norm_16S, color = "16S copies/Litre SW"), size = 1) +
    
    # Chao1 richness line + confidence interval shading
    geom_ribbon(aes(ymin = norm_ci_Chao1_lower, ymax = norm_ci_Chao1_upper, fill = "95% Confidence intervals"), alpha = 0.2) +
    geom_line(aes(y = norm_Chao1, color = "Chao1 richness"), size = 1, linetype = "dashed") +
    
    # Set axis labels and dual y-axis
    scale_y_continuous(
      name = "16S copies/Litre SW", 
      sec.axis = sec_axis(~ ., name = "Chao1 richness")  # Right y-axis
    ) +
    
    # Customize colors and theme
    scale_color_manual(values = c("16S copies/Litre SW" = "#4A4E9E", "Chao1 richness" = "#006D40")) +
    scale_fill_manual(values = c(" " = "#4A4E9E", "95% Confidence intervals" = "#006D40")) +
    
    labs(x = "Date", color = "Gene metric:") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "bold"),  # Bold x-axis values
          axis.text.y = element_text(size = 14, face = "bold"),  # Bold y-axis values
          axis.title.x = element_text(size = 14, face = "bold"),    # Bold axis labels
          axis.title.y = element_text(size = 14, face = "bold")    # Bold axis labels
    ))


# Save the plot
ggsave(
  filename = "outputs/alpha/S16-deep_alpha_counts.png",
  plot = plot,
  width = 6, height = 4,
  dpi = 300
)









# 18S Surface 
chao1_values <- estimate_richness(S18_taxfilt_data, measures = "Chao1")
rownames(chao1_values) <- gsub("\\.", "-", rownames(chao1_values))
chao1_values$SampleID <- rownames(chao1_values)

sample_metadata <- data.frame(sample_data(S18_taxfilt_data))
sample_metadata$SampleID <- rownames(sample_metadata)

sample_metadata <- merge(sample_metadata, chao1_values, by = "SampleID", all.x = TRUE) %>% 
  as.data.frame()


summary_data <- sample_metadata %>%
  filter(Depth..m. == 1) %>%
  group_by(Date) %>%
  summarise(
    mean_18S = mean(`X18S.copies.L.SW`, na.rm = TRUE),
    sd_18S = sd(`X18S.copies.L.SW`, na.rm = TRUE),
    n_18S = sum(!is.na(`X18S.copies.L.SW`)),  # Count non-missing values
    mean_Chao1 = mean(Chao1, na.rm = TRUE),
    sd_Chao1 = sd(Chao1, na.rm = TRUE),
    n_Chao1 = sum(!is.na(Chao1))  # Count non-missing values
  ) %>%
  mutate(
    se_18S = sd_18S / sqrt(n_18S),
    se_Chao1 = sd_Chao1 / sqrt(n_Chao1),
    ci_18S_lower = mean_18S - 1.96 * se_18S,
    ci_18S_upper = mean_18S + 1.96 * se_18S,
    ci_Chao1_lower = mean_Chao1 - 1.96 * se_Chao1,
    ci_Chao1_upper = mean_Chao1 + 1.96 * se_Chao1
  )

# Normalize for dual axis visualization
summary_data <- summary_data %>%
  mutate(
    norm_18S = scale(mean_18S, center = TRUE, scale = TRUE)[,1],
    norm_Chao1 = scale(mean_Chao1, center = TRUE, scale = TRUE)[,1],
    norm_ci_18S_lower = scale(ci_18S_lower, center = TRUE, scale = TRUE)[,1],
    norm_ci_18S_upper = scale(ci_18S_upper, center = TRUE, scale = TRUE)[,1],
    norm_ci_Chao1_lower = scale(ci_Chao1_lower, center = TRUE, scale = TRUE)[,1],
    norm_ci_Chao1_upper = scale(ci_Chao1_upper, center = TRUE, scale = TRUE)[,1]
  )


# Plot with confidence intervals
(plot <- ggplot(summary_data, aes(x = Date)) +
    
    # 12S copies line + confidence interval shading
    geom_ribbon(aes(ymin = norm_ci_18S_lower, ymax = norm_ci_18S_upper, fill = " "), alpha = 0.2) +
    geom_line(aes(y = norm_18S, color = "18S copies/Litre SW"), size = 1) +
    
    # Chao1 richness line + confidence interval shading
    geom_ribbon(aes(ymin = norm_ci_Chao1_lower, ymax = norm_ci_Chao1_upper, fill = "95% Confidence intervals"), alpha = 0.2) +
    geom_line(aes(y = norm_Chao1, color = "Chao1 richness"), size = 1, linetype = "dashed") +
    
    # Set axis labels and dual y-axis
    scale_y_continuous(
      name = "18S copies/Litre SW", 
      sec.axis = sec_axis(~ ., name = "Chao1 richness")  # Right y-axis
    ) +
    
    # Customize colors and theme
    scale_color_manual(values = c("18S copies/Litre SW" = "#4A4E9E", "Chao1 richness" = "#006D40")) +
    scale_fill_manual(values = c(" " = "#4A4E9E", "95% Confidence intervals" = "#006D40")) +
    
    labs(x = "Date", color = "Gene metric:") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "bold"),  # Bold x-axis values
          axis.text.y = element_text(size = 14, face = "bold"),  # Bold y-axis values
          axis.title.x = element_text(size = 14, face = "bold"),    # Bold axis labels
          axis.title.y = element_text(size = 14, face = "bold")    # Bold axis labels
    ))


# Save the plot
ggsave(
  filename = "outputs/alpha/S18-surf_alpha_counts.png",
  plot = plot,
  width = 6, height = 4,
  dpi = 300
)


# 16S Deep

summary_data <- sample_metadata %>%
  filter(Depth..m. > 1) %>%
  group_by(Date) %>%
  summarise(
    mean_18S = mean(`X18S.copies.L.SW`, na.rm = TRUE),
    sd_18S = sd(`X18S.copies.L.SW`, na.rm = TRUE),
    n_18S = sum(!is.na(`X18S.copies.L.SW`)),  # Count non-missing values
    mean_Chao1 = mean(Chao1, na.rm = TRUE),
    sd_Chao1 = sd(Chao1, na.rm = TRUE),
    n_Chao1 = sum(!is.na(Chao1))  # Count non-missing values
  ) %>%
  mutate(
    se_18S = sd_18S / sqrt(n_18S),
    se_Chao1 = sd_Chao1 / sqrt(n_Chao1),
    ci_18S_lower = mean_18S - 1.96 * se_18S,
    ci_18S_upper = mean_18S + 1.96 * se_18S,
    ci_Chao1_lower = mean_Chao1 - 1.96 * se_Chao1,
    ci_Chao1_upper = mean_Chao1 + 1.96 * se_Chao1
  )

# Normalize for dual axis visualization
summary_data <- summary_data %>%
  mutate(
    norm_18S = scale(mean_18S, center = TRUE, scale = TRUE)[,1],
    norm_Chao1 = scale(mean_Chao1, center = TRUE, scale = TRUE)[,1],
    norm_ci_18S_lower = scale(ci_18S_lower, center = TRUE, scale = TRUE)[,1],
    norm_ci_18S_upper = scale(ci_18S_upper, center = TRUE, scale = TRUE)[,1],
    norm_ci_Chao1_lower = scale(ci_Chao1_lower, center = TRUE, scale = TRUE)[,1],
    norm_ci_Chao1_upper = scale(ci_Chao1_upper, center = TRUE, scale = TRUE)[,1]
  )


# Plot with confidence intervals
(plot <- ggplot(summary_data, aes(x = Date)) +
    
    # 12S copies line + confidence interval shading
    geom_ribbon(aes(ymin = norm_ci_18S_lower, ymax = norm_ci_18S_upper, fill = " "), alpha = 0.2) +
    geom_line(aes(y = norm_18S, color = "18S copies/Litre SW"), size = 1) +
    
    # Chao1 richness line + confidence interval shading
    geom_ribbon(aes(ymin = norm_ci_Chao1_lower, ymax = norm_ci_Chao1_upper, fill = "95% Confidence intervals"), alpha = 0.2) +
    geom_line(aes(y = norm_Chao1, color = "Chao1 richness"), size = 1, linetype = "dashed") +
    
    # Set axis labels and dual y-axis
    scale_y_continuous(
      name = "18S copies/Litre SW", 
      sec.axis = sec_axis(~ ., name = "Chao1 richness")  # Right y-axis
    ) +
    
    # Customize colors and theme
    scale_color_manual(values = c("18S copies/Litre SW" = "#4A4E9E", "Chao1 richness" = "#006D40")) +
    scale_fill_manual(values = c(" " = "#4A4E9E", "95% Confidence intervals" = "#006D40")) +
    
    labs(x = "Date", color = "Gene metric:") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "bold"),  # Bold x-axis values
          axis.text.y = element_text(size = 14, face = "bold"),  # Bold y-axis values
          axis.title.x = element_text(size = 14, face = "bold"),    # Bold axis labels
          axis.title.y = element_text(size = 14, face = "bold")    # Bold axis labels
    ))


# Save the plot
ggsave(
  filename = "outputs/alpha/S18-deep_alpha_counts.png",
  plot = plot,
  width = 6, height = 4,
  dpi = 300
)
