# Install once if needed
# install.packages(c("tidyverse","forcats","viridis","Cairo"))

library(tidyverse)
library(forcats)
library(viridis)
library(Cairo)

setwd("/home/ubuntu/NIV-projects/Mosquito-paleovirology/EEfinder/Seq-length-comaparison")

# Read data
df <- read.csv("TED-all-protein.length.csv", sep = "\t")

ref_len  <- 1084
mean_len <- mean(df$length, na.rm = TRUE)

# Order organisms by median length
df$Organism <- fct_reorder(df$Organism, df$length, .fun = median)

# Calculate mean and SD per organism (for error bars)
stats_df <- df %>%
  group_by(Organism) %>%
  summarise(
    mean_len_org = mean(length, na.rm = TRUE),
    sd_len_org   = sd(length, na.rm = TRUE),
    .groups = "drop"
  )

# Plot
p <- ggplot(df, aes(x = Organism, y = length, fill = Organism, color = Organism)) +
  
  # Colourful boxplots (width reflects sample size)
  geom_boxplot(outlier.shape = NA,
               varwidth = TRUE,
               width = 0.7,
               alpha = 0.75,
               linewidth = 0.5) +
  
  # Individual EVE sequences
  geom_jitter(width = 0.15, size = 1.6, alpha = 0.6) +
  
  # Error bars (mean ± SD)  ---- FIXED LAYER
  geom_errorbar(data = stats_df,
                aes(x = Organism,
                    ymin = mean_len_org - sd_len_org,
                    ymax = mean_len_org + sd_len_org),
                inherit.aes = FALSE,
                width = 0.25,
                linewidth = 1.2,
                color = "black") +
  
  # Mean point (black diamond) ---- FIXED LAYER
  geom_point(data = stats_df,
             aes(x = Organism, y = mean_len_org),
             inherit.aes = FALSE,
             color = "black",
             size = 3,
             shape = 18) +
  
  # Reference lines
  geom_hline(yintercept = ref_len,
             color = "blue",
             linetype = "dashed",
             linewidth = 1.3) +
  
  geom_hline(yintercept = mean_len,
             color = "red",
             linetype = "dotted",
             linewidth = 1.3) +
  
  # Color palette
  scale_fill_viridis_d(option = "turbo") +
  scale_color_viridis_d(option = "turbo") +
  
  # Theme
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 60, hjust = 1, size = 13),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.caption = element_text(size = 11)
  ) +
  
  labs(
    y = "Protein Length (aa)",
    x = "Host Organism",
    caption = paste0(
      "Black diamond = Mean ± SD per host; ",
      "Blue dashed = Reference (", ref_len, " aa); ",
      "Red dotted = Overall mean (", round(mean_len,1), " aa)"
    )
  )

# Preview
print(p)

# Save vector PDF
ggsave("TED-Seq-length-comparison.pdf",
       plot = p,
       width = 15,
       height = 10,
       device = cairo_pdf)

# Save high-resolution JPEG (600 dpi)
ggsave("TED-Seq-length-comparison.jpeg",
       plot = p,
       width = 15,
       height = 10,
       dpi = 600,
       units = "in",
       device = "jpeg")
       
#################################################################################################################################333
#for clado plot

df <- read.csv("Clado-all-protein.length.csv", sep = ",")

ref_len  <- 1045
mean_len <- mean(df$length, na.rm = TRUE)

# Order organisms by median length
df$Organism <- fct_reorder(df$Organism, df$length, .fun = median)

# Calculate mean and SD per organism (for error bars)
stats_df <- df %>%
    group_by(Organism) %>%
    summarise(
        mean_len_org = mean(length, na.rm = TRUE),
        sd_len_org   = sd(length, na.rm = TRUE),
        .groups = "drop"
    )

# Plot
p <- ggplot(df, aes(x = Organism, y = length, fill = Organism, color = Organism)) +
    
    # Colourful boxplots (width reflects sample size)
    geom_boxplot(outlier.shape = NA,
                 varwidth = TRUE,
                 width = 0.7,
                 alpha = 0.75,
                 linewidth = 0.5) +
    
    # Individual EVE sequences
    geom_jitter(width = 0.15, size = 1.6, alpha = 0.6) +
    
    # Error bars (mean ± SD)  ---- FIXED LAYER
    geom_errorbar(data = stats_df,
                  aes(x = Organism,
                      ymin = mean_len_org - sd_len_org,
                      ymax = mean_len_org + sd_len_org),
                  inherit.aes = FALSE,
                  width = 0.25,
                  linewidth = 1.2,
                  color = "black") +
    
    # Mean point (black diamond) ---- FIXED LAYER
    geom_point(data = stats_df,
               aes(x = Organism, y = mean_len_org),
               inherit.aes = FALSE,
               color = "black",
               size = 3,
               shape = 18) +
    
    # Reference lines
    geom_hline(yintercept = ref_len,
               color = "blue",
               linetype = "dashed",
               linewidth = 1.3) +
    
    geom_hline(yintercept = mean_len,
               color = "red",
               linetype = "dotted",
               linewidth = 1.3) +
    
    # Color palette
    scale_fill_viridis_d(option = "turbo") +
    scale_color_viridis_d(option = "turbo") +
    
    # Theme
    theme_bw(base_size = 14) +
    theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 60, hjust = 1, size = 13),
        axis.title = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.caption = element_text(size = 11)
    ) +
    
    labs(
        y = "Protein Length (aa)",
        x = "Host Organism",
        caption = paste0(
            "Black diamond = Mean ± SD per host; ",
            "Blue dashed = Reference (", ref_len, " aa); ",
            "Red dotted = Overall mean (", round(mean_len,1), " aa)"
        )
    )

# Preview
print(p)

# Save vector PDF
ggsave("Clado-Seq-length-comparison.pdf",
       plot = p,
       width = 15,
       height = 10,
       device = cairo_pdf)

# Save high-resolution JPEG (600 dpi)
ggsave("Clado-Seq-length-comparison.jpeg",
       plot = p,
       width = 15,
       height = 10,
       dpi = 600,
       units = "in",
       device = "jpeg")
