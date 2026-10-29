
# ============================================================
# INTEGRATED FIGURE 1
# A. Overall classification of candidate regions
# B. Category composition across 19 mosquito genomes
# C. Top 30 viral-family annotations
#
# Input: TableS1.csv (tab-delimited)
#
# Outputs:
#   Figure1_candidate_classification_with_viral_families.pdf
#   Figure1_candidate_classification_with_viral_families.tiff
#   Figure1A_category_summary.csv
#   Figure1B_genome_category_summary.csv
#   Figure1C_top30_viral_families.csv
#   Figure1C_all_viral_families.csv
# ============================================================


# ============================================================
# 1. SET WORKING DIRECTORY AND LOAD PACKAGES
# ============================================================

setwd(
  "/home/ubuntu/NIV-projects/Mosquito-paleovirology/manuscript/scientific-report/comments/Figure-script/revision-codes/Figure.1"
)

input_file <- "TableS1.csv"

required_packages <- c(
  "tidyverse",
  "scales",
  "patchwork"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(tidyverse)
library(scales)
library(patchwork)

options(na.print = "NA")


# ============================================================
# 2. READ INPUT DATA
# ============================================================

dat <- read.delim(
  input_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = ""
)

cat("====================================================\n")
cat("DATA SUMMARY\n")
cat("====================================================\n")

cat("Total rows:", nrow(dat), "\n")
cat("Unique Element.ID values:", length(unique(dat$Element.ID)), "\n")
cat("Number of mosquito genomes:", length(unique(dat$Organism.Name)), "\n")

required_columns <- c(
  "Element.ID",
  "Organism.Name",
  "Family",
  "Genus",
  "Category"
)

missing_columns <- setdiff(required_columns, colnames(dat))

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

dat <- dat %>%
  mutate(
    Element.ID = trimws(as.character(Element.ID)),
    Organism.Name = trimws(as.character(Organism.Name)),
    Family = trimws(as.character(Family)),
    Category = trimws(as.character(Category))
  ) %>%
  filter(
    !is.na(Category),
    Category != "",
    !is.na(Element.ID),
    Element.ID != ""
  )


# ============================================================
# 3. DEFINE CATEGORY ORDER AND COLOURS
# ============================================================

category_levels <- c(
  "LTR_retrotransposon",
  "Viral_reference",
  "Other_TE",
  "DNA_or_other_transposon",
  "No_significant_hit"
)

category_labels <- c(
  "LTR retrotransposon-related",
  "Viral-reference-like",
  "Other TE-related",
  "DNA/other transposon-related",
  "No significant protein hit"
)

category_colors <- c(
  "LTR_retrotransposon"     = "#F8766D",
  "Viral_reference"         = "#00B0F6",
  "Other_TE"                = "#00BF7D",
  "DNA_or_other_transposon" = "#C49A00",
  "No_significant_hit"      = "#D966E8"
)

observed_categories <- unique(dat$Category)

unexpected_categories <- setdiff(
  observed_categories,
  category_levels
)

if (length(unexpected_categories) > 0) {
  stop(
    "Unexpected categories found in the input: ",
    paste(unexpected_categories, collapse = ", "),
    "\nPlease check the Category column before plotting."
  )
}

dat <- dat %>%
  mutate(
    Category = factor(
      Category,
      levels = category_levels
    )
  )


# ============================================================
# 4. PANEL A: OVERALL CLASSIFICATION PIE CHART
# ============================================================

category_summary <- dat %>%
  count(
    Category,
    name = "Count",
    .drop = FALSE
  ) %>%
  mutate(
    Percentage = 100 * Count / sum(Count),

    Category_label = category_labels[
      match(
        as.character(Category),
        category_levels
      )
    ],

    Legend_label = sprintf(
      "%s: %s (%.1f%%)",
      Category_label,
      comma(Count),
      Percentage
    ),

    Slice_label = if_else(
      Percentage >= 5,
      sprintf("%.1f%%", Percentage),
      ""
    )
  )

cat("\nPANEL A: OVERALL CATEGORY DISTRIBUTION\n")

base::print(
  as.data.frame(
    category_summary %>%
      select(Category, Count, Percentage)
  ),
  row.names = FALSE
)

p1a <- ggplot(
  category_summary,
  aes(
    x = "",
    y = Count,
    fill = Category
  )
) +
  geom_col(
    width = 1,
    colour = "white",
    linewidth = 0.7
  ) +
  coord_polar(theta = "y") +
  geom_text(
    aes(label = Slice_label),
    position = position_stack(vjust = 0.5),
    colour = "white",
    fontface = "bold",
    size = 4
  ) +
  scale_fill_manual(
    values = category_colors,
    breaks = category_levels,
    labels = category_summary$Legend_label[
      match(
        category_levels,
        as.character(category_summary$Category)
      )
    ],
    drop = FALSE
  ) +
  labs(
    title = "A",
    subtitle = "Overall classification of candidate regions",
    fill = "Category: count (percentage)"
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16
    ),
    plot.subtitle = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_text(
      face = "bold",
      size = 9
    ),
    legend.text = element_text(size = 8),
    legend.key.size = grid::unit(0.35, "cm"),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  guides(
    fill = guide_legend(
      ncol = 1,
      byrow = TRUE
    )
  )


# ============================================================
# 5. PANEL B: CATEGORY COMPOSITION ACROSS MOSQUITO GENOMES
# ============================================================

species_summary <- dat %>%
  count(
    Organism.Name,
    Category,
    name = "Count",
    .drop = FALSE
  ) %>%
  complete(
    Organism.Name,
    Category = factor(
      category_levels,
      levels = category_levels
    ),
    fill = list(Count = 0)
  ) %>%
  group_by(Organism.Name) %>%
  mutate(
    Total = sum(Count),
    Percentage = if_else(
      Total > 0,
      100 * Count / Total,
      0
    )
  ) %>%
  ungroup()

# Order genomes by total candidate-region count
species_order <- dat %>%
  count(
    Organism.Name,
    sort = TRUE
  ) %>%
  pull(Organism.Name)

species_summary <- species_summary %>%
  mutate(
    Organism.Name = factor(
      Organism.Name,
      levels = rev(species_order)
    )
  )

cat("\nPANEL B: CATEGORY COMPOSITION BY MOSQUITO GENOME\n")

base::print(
  as.data.frame(
    species_summary %>%
      select(
        Organism.Name,
        Category,
        Count,
        Total,
        Percentage
      )
  ),
  row.names = FALSE
)

p1b <- ggplot(
  species_summary,
  aes(
    x = Percentage,
    y = Organism.Name,
    fill = Category
  )
) +
  geom_col(
    width = 0.75,
    colour = "white",
    linewidth = 0.2
  ) +
  scale_fill_manual(
    values = category_colors,
    breaks = category_levels,
    labels = category_labels,
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = seq(0, 100, by = 20),
    limits = c(0, 100),
    expand = c(0, 0),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    title = "B",
    subtitle = "Category composition across mosquito genomes",
    x = "Candidate regions within each genome (%)",
    y = NULL,
    fill = "Candidate category"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16
    ),
    plot.subtitle = element_text(size = 10),
    axis.text.y = element_text(
      size = 8.5,
      colour = "black"
    ),
    axis.text.x = element_text(
      size = 9,
      colour = "black"
    ),
    axis.title.x = element_text(size = 10),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(
      face = "bold",
      size = 9
    ),
    legend.text = element_text(size = 8),
    legend.key.size = grid::unit(0.35, "cm"),
    plot.margin = margin(5, 10, 5, 5)
  ) +
  guides(
    fill = guide_legend(
      nrow = 2,
      byrow = TRUE
    )
  )


# ============================================================
# 6. PANEL C: TOP 30 VIRAL-FAMILY ANNOTATIONS
# ============================================================

# Restrict this analysis to Viral_reference candidates.
# Family labels are the original annotations from TableS1.csv.
# These annotations do not independently establish viral ancestry
# or confirm integration into the mosquito genome.

viral_dat <- dat %>%
  filter(
    as.character(Category) == "Viral_reference"
  ) %>%
  mutate(
    Family = na_if(trimws(Family), ""),
    Family = if_else(
      is.na(Family) |
        tolower(Family) %in% c(
          "undefined",
          "unknown",
          "na"
        ),
      "Unclassified",
      Family
    )
  ) %>%
  distinct(
    Element.ID,
    Family
  )

total_viral_loci <- n_distinct(viral_dat$Element.ID)

cat(
  "\nUnique viral-reference-like candidate loci:",
  total_viral_loci,
  "\n"
)

# Count unique Element.ID–Family combinations
family_summary_all <- viral_dat %>%
  count(
    Family,
    name = "Candidate_loci"
  ) %>%
  mutate(
    Percentage = 100 * Candidate_loci / total_viral_loci
  ) %>%
  arrange(
    desc(Candidate_loci),
    Family
  ) %>%
  mutate(
    Rank = row_number(),
    Percentage = round(Percentage, 2)
  ) %>%
  select(
    Rank,
    Family,
    Candidate_loci,
    Percentage
  )

# Export the complete family table
write.csv(
  family_summary_all,
  "Figure1C_all_viral_families.csv",
  row.names = FALSE,
  na = ""
)

# Exclude unclassified labels when selecting the 30 named families
top30_families <- family_summary_all %>%
  filter(
    !tolower(Family) %in% c(
      "unclassified",
      "undefined",
      "unknown",
      "na"
    )
  ) %>%
  slice_max(
    order_by = Candidate_loci,
    n = 30,
    with_ties = FALSE
  ) %>%
  arrange(
    Candidate_loci,
    Family
  ) %>%
  mutate(
    Family = factor(
      Family,
      levels = Family
    ),
    Label = paste0(
      comma(Candidate_loci),
      " (",
      sprintf("%.1f", Percentage),
      "%)"
    )
  )

# Export the top-30 table
write.csv(
  top30_families %>%
    mutate(Family = as.character(Family)) %>%
    select(
      Rank,
      Family,
      Candidate_loci,
      Percentage
    ),
  "Figure1C_top30_viral_families.csv",
  row.names = FALSE,
  na = ""
)

cat("\nPANEL C: TOP 30 VIRAL FAMILIES\n")

base::print(
  as.data.frame(
    top30_families %>%
      mutate(Family = as.character(Family)) %>%
      select(
        Rank,
        Family,
        Candidate_loci,
        Percentage
      )
  ),
  row.names = FALSE
)

p1c <- ggplot(
  top30_families,
  aes(
    x = Candidate_loci,
    y = Family
  )
) +
  geom_col(
    fill = "#08A9E6",
    colour = "black",
    linewidth = 0.2,
    width = 0.75
  ) +
  geom_text(
    aes(label = Label),
    hjust = -0.08,
    size = 2.8,
    colour = "black"
  ) +
  scale_x_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.23))
  ) +
  labs(
    title = "C",
    subtitle = paste0(
      "Top 30 family annotations among viral-reference-like candidates (N = ",
      comma(total_viral_loci),
      ")"
    ),
    x = "Number of unique candidate loci",
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16
    ),
    plot.subtitle = element_text(size = 10),
    axis.text.y = element_text(
      size = 8,
      colour = "black"
    ),
    axis.text.x = element_text(
      size = 9,
      colour = "black"
    ),
    axis.title.x = element_text(
      face = "bold",
      size = 10
    ),
    plot.margin = margin(5, 35, 5, 5)
  )


# ============================================================
# 7. INTEGRATE PANELS A, B AND C
# ============================================================

# Layout:
# -------------------------
# |     A     |     B     |
# -------------------------
# |           C           |
# -------------------------

top_row <- p1a | p1b

figure1 <- (
  top_row /
    p1c
) +
  plot_layout(
    heights = c(1, 1.35),
    guides = "keep"
  ) +
  plot_annotation(
    title = paste(
      "Comparative classification of candidate",
      "viral-homologous regions"
    ),
    subtitle = paste(
      "Initial EEFinder candidates across 19 Culicinae mosquito genomes;",
      "categories based on viral and transposable-element reference matches"
    ),
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 15
      ),
      plot.subtitle = element_text(size = 10)
    )
  )

print(figure1)


# ============================================================
# 8. EXPORT PUBLICATION-QUALITY FIGURES
# ============================================================

# Taller canvas accommodates all 30 family labels.

ggsave(
  filename = "Figure1_candidate_classification_with_viral_families.pdf",
  plot = figure1,
  width = 16,
  height = 17,
  units = "in",
  device = grDevices::cairo_pdf,
  limitsize = FALSE
)

ggsave(
  filename = "Figure1_candidate_classification_with_viral_families.tiff",
  plot = figure1,
  width = 16,
  height = 17,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white",
  limitsize = FALSE
)


# ============================================================
# 9. EXPORT PANELS A AND B SUMMARY TABLES
# ============================================================

write.csv(
  category_summary %>%
    mutate(Category = as.character(Category)) %>%
    select(
      Category,
      Count,
      Percentage
    ),
  file = "Figure1A_category_summary.csv",
  row.names = FALSE
)

write.csv(
  species_summary %>%
    mutate(
      Organism.Name = as.character(Organism.Name),
      Category = as.character(Category)
    ),
  file = "Figure1B_genome_category_summary.csv",
  row.names = FALSE
)


# ============================================================
# 10. FINAL CHECKS
# ============================================================

cat("\n====================================================\n")
cat("INTEGRATED FIGURE GENERATION COMPLETED\n")
cat("====================================================\n")

cat("Input rows used:", nrow(dat), "\n")
cat(
  "Mosquito genomes:",
  n_distinct(dat$Organism.Name),
  "\n"
)
cat(
  "Rows represented in Panel A:",
  sum(category_summary$Count),
  "\n"
)
cat(
  "Unique viral-reference-like candidate loci:",
  total_viral_loci,
  "\n"
)
cat(
  "Named families displayed in Panel C:",
  nrow(top30_families),
  "\n"
)
cat(
  "Total family annotations in exported table:",
  nrow(family_summary_all),
  "\n"
)

cat("\nOutput files:\n")
cat("1. Figure1_candidate_classification_with_viral_families.pdf\n")
cat("2. Figure1_candidate_classification_with_viral_families.tiff\n")
cat("3. Figure1A_category_summary.csv\n")
cat("4. Figure1B_genome_category_summary.csv\n")
cat("5. Figure1C_top30_viral_families.csv\n")
cat("6. Figure1C_all_viral_families.csv\n")

cat("\nAll done.\n")

