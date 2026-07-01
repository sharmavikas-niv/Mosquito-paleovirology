setwd("/home/ubuntu/NIV-projects/Mosquito-paleovirology/github/Neha/figure2")
library(readr)

#Read file
merged_output <- read.delim("merged_all_columns_nonoverlapping.tsv", stringsAsFactors = FALSE)

View(merged_output)

head(merged_output)
#check
length(merged_output$Organism.Name)
length(unique(merged_output$Organism.Name))

#count EVE
eve_counts <- as.data.frame(table(merged_output$Organism.Name))
colnames(eve_counts) <- c("Genome", "EVE_Count")

eve_counts <- eve_counts[order(eve_counts$EVE_Count, decreasing = TRUE), ]

View(eve_counts)

#Addition of genome data for normaliaziton 
df <- read.csv("mosquito-genomes-06-03-2025.csv", header = TRUE)


# STEP 1: Count EVEs per genome
library(dplyr)

eve_counts <- merged_output %>%
  count(Organism.Name, name = "EVE_Count")

# STEP 2: Select relevant genome size info
genome_info <- df %>%
  select(Organism.Name, GenomeSize = Assembly.Stats.Total.Sequence.Length) %>%
  distinct()

# STEP 3: Merge EVE count with genome size
merged_data <- left_join(eve_counts, genome_info, by = "Organism.Name")

# STEP 4: Calculate EVEs per megabase
merged_data <- merged_data %>%
  mutate(EVE_per_Mb = EVE_Count / (GenomeSize / 1e6))

View(merged_data)

head(merged_data)
#Step 4: Save to file
write.csv(merged_data, "genome_added_EVE_mb_count.csv", row.names = FALSE)

#########################################################################################
#plot
library(readr)

#EVE count
merged_data <- read.csv("genome_added_EVE_mb_count.csv", stringsAsFactors = FALSE)

# normalized plot
library(ggplot2)


P1 <- ggplot(merged_data,
             aes(x = reorder(Organism.Name, -EVE_per_Mb),
                 y = EVE_per_Mb)) +
  
  geom_bar(stat = "identity", width = 0.8, fill = "#f4a582", color = "black") +
  
  #  Inside bar (formula)
  geom_text(
    aes(
      y = EVE_per_Mb / 2,
      label = paste0(EVE_Count, "/", round(GenomeSize/1e6, 1)),
      size = ifelse(rank(EVE_per_Mb) <= 5, 2.4, 3)   #  last 5 bars smaller and it will be 2.4 in size
    ),
    angle = 90,
    color = "black"
  ) +
  scale_size_identity() +
  #  Top of bar (density value)
  geom_text(aes(label = sprintf("%.2f", EVE_per_Mb)),
            vjust = -0.6,
            angle = 45,
            hjust = 0,
            size = 3.5) +
  
  expand_limits(y = max(merged_data$EVE_per_Mb) * 1.1) +
  
  theme_minimal() +
  
  labs(
    x = "Mosquito Genomes",
    y = " Genome-normalized EVE density"
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 60, 
      hjust = 1, 
      color = "black",      # only black (no bold)
      face = "plain", 
      size = 10
    ),
    axis.text.y = element_text(
      color = "black",      # black color
      face = "bold",        # make bold
      size = 10
    ),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    plot.margin = margin(10, 10, 10, 10)
  )
print(P1)

######################################################################################
#file creation for for EVE percentage 

# Load required package
library(dplyr)



# 1. Read the TSV file
df <- read.table(
  "merged_all_columns_nonoverlapping.tsv",
  sep = "\t",
  header = TRUE,
  quote = "",
  fill = TRUE,
  stringsAsFactors = FALSE
)

# 2. Summarize total region_bp per organism
bp_per_organism <- df %>%
  group_by(Organism.Name) %>%
  summarise(total_region_bp = sum(region_bp, na.rm = TRUE))

print(head(bp_per_organism))

View(bp_per_organism)

# 3. Read the genome assembly CSV file
genomes <- read.csv(
  "mosquito-genomes-06-03-2025.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)

# Add Assembly-length to bp_per_organism
bp_per_organism <- bp_per_organism %>%
  left_join(genomes %>% select(Organism.Name, Assembly.length),
            by = "Organism.Name")

View(bp_per_organism)
print(bp_per_organism)

#percentage
bp_per_organism <- bp_per_organism %>%
  mutate(EVE_percentage = (total_region_bp / Assembly.length) * 100)

#descendig order
bp_per_organism <- bp_per_organism %>%
  arrange(desc(EVE_percentage))


View(bp_per_organism)
#save
write.table(bp_per_organism, "EVE_summary_by_organism.tsv", sep = "\t", quote = FALSE, row.names = FALSE)



########################################################################################

#EVE Percentage Plot
#Read file
bp_per_organism <- read.delim("EVE_summary_by_organism.tsv", stringsAsFactors = FALSE)
#plot with formula inside
library(ggplot2)

P2 <- ggplot(bp_per_organism,
             aes(x = reorder(Organism.Name, -EVE_percentage),
                 y = EVE_percentage)) +
  geom_bar(stat = "identity", width = 0.8, fill = "#b39ddb", color = "black") +
  geom_text(
    aes(
      label = paste0(total_region_bp, "/", Assembly.length),
      size = ifelse(rank(EVE_percentage) <= 2, 2.4, 3)   #  last 2 bars smaller
    ),
    position = position_stack(vjust = 0.5),
    angle = 90,
    color = "black"
  ) +
  scale_size_identity() +
  #  Top of bar (density value)
  geom_text(aes(label = sprintf("%.2f", EVE_percentage)),
            vjust = -0.6,
            angle = 45,
            hjust = 0,
            size = 3.5) +
  labs(x = "Mosquito Genomes", y = "Percentage of genome occupied by EVEs") +
  scale_y_continuous(expand = c(0, 0)) +   #remove padding
  coord_cartesian(clip = "off") + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 60, 
      hjust = 1, 
      color = "black",      # only black (no bold)
      face = "plain", 
      size = 10
    ),
    axis.text.y = element_text(
      color = "black",      # black color
      face = "bold",        # make bold
      size = 10
    ),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    plot.margin = margin(10, 10, 10, 10)
  )
print(P2)

#########################################################################
#collage
library(patchwork)


# Side-by-side collage
#labelling of plot
combined_plot <- (P1 | P2) +
  plot_annotation(tag_levels = "A")

combined_plot <- combined_plot &
  theme(plot.margin = margin(15, 10, 10, 10))

combined_plot
