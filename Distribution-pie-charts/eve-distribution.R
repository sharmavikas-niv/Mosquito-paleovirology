setwd("/home/ubuntu/NIV-projects/Mosquito-paleovirology/github/Distribution-pie-charts")

# Load required libraries
library(ggpubr)
library(rstatix)
library(ggplot2)
library(plyr)
library(pheatmap)
library(tidyverse)
library(reshape2)
library(GenomicRanges)
library(dplyr)
library(stringr)
library(tidyr)
library(ggpubr)
library(scales)

# Read input files
df <- read.delim("mosquito-genomes-06-03-2025.csv", sep = "\t", header = TRUE)
merged_output <- read.delim("merged_all_columns_nonoverlapping.tsv", sep = "\t", header = TRUE)
# Ensure start and end are numeric
merged_output$start <- as.numeric(merged_output$start)
merged_output$end <- as.numeric(merged_output$end)

# Calculate region length
merged_output <- merged_output %>%
  mutate(region_bp = end - start + 1)

################################################################################################
#pie chart 
# Load libraries
library(dplyr)
library(ggplot2)
library(ggrepel)
library(rlang)
library(forcats)

# Function to make a pie chart with numbers + percentages in labels
make_pie <- function(df, column, top_n = 20, file_name = NULL, start_angle = -pi/2) {
    col_q <- enquo(column)
    
    # Prepare data
    plot_data <- df %>%
        count(!!col_q, sort = TRUE) %>%
        mutate(category = as.character(!!col_q)) %>%   
        mutate(category = ifelse(row_number() > top_n, "Others", category)) %>%
        group_by(category) %>%
        summarise(n = sum(n), .groups = "drop") %>%
        arrange(desc(n)) %>%
        mutate(
            pct   = n / sum(n) * 100,
            label = paste0(category, " (", n, ", ", round(pct, 1), "%)")
        )
    
    # Compute label positions safely
    plot_data <- plot_data %>%
        mutate(
            csum = rev(cumsum(rev(n))),
            pos  = n / 2 + lead(csum, 1),
            pos  = if_else(is.na(pos), n / 2, pos)
        ) %>%
        filter(!is.na(pos) & n > 0)
    
    # Base pie chart
    p <- ggplot(plot_data, aes(x = 1, y = n, fill = fct_inorder(category))) +
        geom_col(width = 1, color = "white") +
        coord_polar(theta = "y", start = start_angle) +
        scale_x_continuous(limits = c(0, 1.5)) +  
        theme_void() +
        theme(legend.position = "none") +
        labs(title = paste("Distribution of", as_name(col_q)))
    
    # Add external labels using computed positions
    p <- p +
        geom_label_repel(
            data = plot_data,
            aes(x = 1.3, y = pos, label = label, fill = category),
            inherit.aes = FALSE,
            size = 2,
            color = "black",
            box.padding = 0.35,
            point.padding = 0.3,
            segment.size = 0.35,
            segment.color = "grey30",
            max.overlaps = Inf,
            show.legend = FALSE
        )
    
    # Save plot if file_name is given
    if (!is.null(file_name)) {
        ggsave(file_name, p, width = 6, height = 6, device = "pdf")
    }
    
    return(p)
}

# Example usage
p1 <- make_pie(merged_output, Molecule_type, top_n = 20, file_name = "Molecule_type_pie_real-number.pdf")
p2 <- make_pie(merged_output, Family,        top_n = 20, file_name = "Family_pie_real-number.pdf")
p3 <- make_pie(merged_output, Genus,         top_n = 20, file_name = "Genus_pie_real-number.pdf")
p4 <- make_pie(merged_output, Species,       top_n = 20, file_name = "Species_pie_real-number.pdf")

