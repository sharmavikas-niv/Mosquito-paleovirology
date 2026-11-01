###############################################################
# MOSQUITO EVE ANALYSIS: Pfam07253 PUTATIVE ENVELOPE DOMAIN
# Genomic context: 5 kb upstream + 5 kb downstream
#
# Outputs:
# 1. Input QC summary
# 2. Classified domain-hit table
# 3. Domain prevalence across unique loci
# 4. Pfam07253 prevalence plot
# 5. Domain co-occurrence plot
# 6. Domain presence/absence heatmap
# 7. Pfam07253 genomic-position plot (if coordinates permit)
# 8. CSV tables for downstream analysis
#
# IMPORTANT:
# Pfam07253 detection alone does not establish a functional
# envelope protein, viral ancestry, or EVE integration.
###############################################################
setwd("/home/ubuntu/NIV-projects/Mosquito-paleovirology/manuscript/scientific-report/comments/Figure-script/revision-codes/5kb-up-down/TED")
rm(list = ls())

options(stringsAsFactors = FALSE)

# =============================================================
# 1. USER CONFIGURATION
# =============================================================

# Change this to your analysis directory.
workdir <- "."

# Change this to the actual mosquito EVE/flanking-region TSV.
input_file <- file.path(
  workdir,
  "Trichoplusia_ni_TED_virus_orfB-5KB.fasta-rpsblast-seg_mapped.tsv"
)

outdir <- file.path(
  workdir,
  "Mosquito_EVE_Pfam07253_results"
)

# Coordinate settings:
# Assumes a 10,000-bp query with the reference point at position 5000.
# Relative coordinates will range approximately from -5 kb to +5 kb.
WINDOW_BP <- 10000
ORIGIN_POSITION <- 5000

# Set TRUE only if qstart/qend are positions on the 10-kb
# mosquito flanking sequence and the reference point is correct.
PLOT_FLANK_COORDINATES <- TRUE

# Set to FALSE to skip the potentially large heatmap.
MAKE_HEATMAP <- TRUE

PNG_DPI <- 600

# Expected input columns:
# 1  locus
# 2  hit_id
# 3  pident
# 4  alignment_length
# 5  mismatch
# 6  gapopen
# 7  qstart
# 8  qend
# 9  sstart
# 10 send
# 11 evalue
# 12 bitscore
# 13 accession
# 14 domain
# 15 description
# 16 domain_length

# =============================================================
# 2. INSTALL AND LOAD PACKAGES
# =============================================================

packages <- c(
  "data.table",
  "dplyr",
  "tidyr",
  "stringr",
  "ggplot2",
  "pheatmap"
)

missing <- packages[
  !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing) > 0) {
  install.packages(missing)
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(pheatmap)
})

# =============================================================
# 3. CHECK INPUT AND CREATE OUTPUT DIRECTORY
# =============================================================

if (!file.exists(input_file)) {
  stop(
    "Input file not found: ", input_file,
    "\nPlease update workdir and input_file at the top of the script."
  )
}

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

message("Input: ", normalizePath(input_file))
message("Output directory: ", normalizePath(outdir))

# =============================================================
# 4. READ THE TAB-DELIMITED INPUT
# =============================================================

dat <- fread(
  input_file,
  sep = "\t",
  header = FALSE,
  quote = "",
  fill = TRUE,
  data.table = FALSE
)

if (ncol(dat) < 16) {
  stop(
    "Expected at least 16 columns, but found ",
    ncol(dat),
    ". Check the TSV format and column order."
  )
}

column_names <- c(
  "locus",
  "hit_id",
  "pident",
  "alignment_length",
  "mismatch",
  "gapopen",
  "qstart",
  "qend",
  "sstart",
  "send",
  "evalue",
  "bitscore",
  "accession",
  "domain",
  "description",
  "domain_length"
)

colnames(dat)[1:16] <- column_names

# Convert numeric columns.
numeric_cols <- c(
  "pident",
  "alignment_length",
  "mismatch",
  "gapopen",
  "qstart",
  "qend",
  "sstart",
  "send",
  "evalue",
  "bitscore",
  "domain_length"
)

for (nm in numeric_cols) {
  dat[[nm]] <- suppressWarnings(as.numeric(dat[[nm]]))
}

text_cols <- c(
  "locus",
  "hit_id",
  "accession",
  "domain",
  "description"
)

for (nm in text_cols) {
  dat[[nm]] <- trimws(as.character(dat[[nm]]))
}

dat <- dat %>%
  filter(
    !is.na(locus),
    locus != ""
  )

if (nrow(dat) == 0) {
  stop("No valid domain-hit records were found.")
}

# =============================================================
# 5. BASIC QUALITY CONTROL
# =============================================================

total_records <- nrow(dat)

total_loci <- n_distinct(dat$locus)

total_domain_labels <- n_distinct(
  dat$domain[!is.na(dat$domain) & dat$domain != ""]
)

qc_summary <- data.frame(
  Metric = c(
    "Input file",
    "Total domain-hit records",
    "Unique query/locus identifiers",
    "Unique domain labels"
  ),
  Value = c(
    basename(input_file),
    total_records,
    total_loci,
    total_domain_labels
  )
)

write.csv(
  qc_summary,
  file.path(outdir, "01_QC_summary.csv"),
  row.names = FALSE
)

# =============================================================
# 6. CLASSIFY DOMAIN HITS
# =============================================================

# Pfam07253 and Pfam12259 are handled separately.
# Check your original annotation table to confirm that accession
# strings and descriptions follow the expected naming conventions.

dat2 <- dat %>%
  mutate(
    search_text = str_to_lower(
      paste(
        coalesce(domain, ""),
        coalesce(accession, ""),
        coalesce(description, "")
      )
    ),

    # Exact Pfam07253 accession recognition.
    Envelope = str_detect(
      search_text,
      "pfam[:._ -]?07253|pfam07253"
    ),

    # Separate Baculo_F-related annotation.
    Baculo_F = str_detect(
      search_text,
      "pfam[:._ -]?12259|baculo[_ -]?f"
    ),

    RT = str_detect(
      search_text,
      "rvt_1|rvt_2|rt_rnaseh|reverse transcriptase|reverse-transcriptase"
    ),

    RNaseH = str_detect(
      search_text,
      "rnase[ _.-]?h|ribonuclease h"
    ),

    Integrase = str_detect(
      search_text,
      "integrase|\\brve\\b"
    ),

    Protease = str_detect(
      search_text,
      "protease|asp_protease"
    ),

    Gag = str_detect(
      search_text,
      "\\bgag\\b|retrotrans_gag"
    ),

    Gypsy = str_detect(
      search_text,
      "gypsy|ty3"
    )
  )

write.csv(
  dat2,
  file.path(outdir, "02_Classified_domain_hits.csv"),
  row.names = FALSE
)

# =============================================================
# 7. DOMAIN FREQUENCY ACROSS UNIQUE LOCI
# =============================================================

domain_frequency <- dat2 %>%
  filter(
    !is.na(domain),
    domain != ""
  ) %>%
  distinct(locus, domain) %>%
  count(domain, name = "Unique_loci") %>%
  mutate(
    Total_loci = total_loci,
    Percent_loci = 100 * Unique_loci / total_loci
  ) %>%
  arrange(desc(Unique_loci))

write.csv(
  domain_frequency,
  file.path(outdir, "03_Domain_frequency.csv"),
  row.names = FALSE
)

# Plot the 30 most frequent domain labels.
if (nrow(domain_frequency) > 0) {

  top_domains <- domain_frequency %>%
    slice_head(n = min(30, nrow(domain_frequency)))

  p_domains <- ggplot(
    top_domains,
    aes(
      x = Unique_loci,
      y = reorder(domain, Unique_loci)
    )
  ) +
    geom_col(
      width = 0.72,
      fill = "grey45"
    ) +
    geom_text(
      aes(label = Unique_loci),
      hjust = -0.15,
      size = 3.2
    ) +
    scale_x_continuous(
      expand = expansion(mult = c(0, 0.15))
    ) +
    labs(
      title = "Domain diversity across analysed loci",
      subtitle = "Counts represent unique loci",
      x = "Number of unique loci",
      y = "Domain annotation"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black")
    )

  plot_height <- max(
    5,
    0.25 * nrow(top_domains)
  )

  ggsave(
    file.path(outdir, "04_Top_domains.pdf"),
    p_domains,
    width = 9,
    height = plot_height,
    device = "pdf"
  )

  ggsave(
    file.path(outdir, "04_Top_domains.png"),
    p_domains,
    width = 9,
    height = plot_height,
    dpi = PNG_DPI,
    bg = "white"
  )
}

# =============================================================
# 8. CREATE LOCUS-LEVEL DOMAIN PRESENCE/ABSENCE TABLE
# =============================================================

domain_columns <- c(
  "RT",
  "RNaseH",
  "Protease",
  "Gag",
  "Integrase",
  "Gypsy",
  "Baculo_F",
  "Envelope"
)

locus_domains <- dat2 %>%
  group_by(locus) %>%
  summarise(
    across(
      all_of(domain_columns),
      ~ any(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

write.csv(
  locus_domains,
  file.path(outdir, "05_Locus_domain_presence_absence.csv"),
  row.names = FALSE
)

# =============================================================
# 9. SUMMARISE MAJOR DOMAIN CLASSES
# =============================================================

feature_labels <- c(
  RT = "Reverse transcriptase-associated",
  RNaseH = "RNase H",
  Protease = "Protease",
  Gag = "Gag",
  Integrase = "Integrase",
  Gypsy = "Gypsy/Ty3-associated",
  Baculo_F = "Putative envelope-associated domain / Baculo_F (Pfam12259)",
  Envelope = "Putative envelope-associated domain (Pfam07253)"
)

feature_summary <- data.frame(
  Feature = names(feature_labels),
  Description = unname(feature_labels),
  Unique_loci = vapply(
    names(feature_labels),
    function(x) {
      sum(locus_domains[[x]], na.rm = TRUE)
    },
    numeric(1)
  )
) %>%
  mutate(
    Total_loci = total_loci,
    Percent_loci = 100 * Unique_loci / total_loci
  )

write.csv(
  feature_summary,
  file.path(outdir, "06_Major_domain_summary.csv"),
  row.names = FALSE
)

# =============================================================
# 10. FOCUS ON PFAM07253
# =============================================================

envelope_hits <- dat2 %>%
  filter(Envelope)

envelope_loci <- unique(envelope_hits$locus)

n_envelope <- length(envelope_loci)

n_without_envelope <- total_loci - n_envelope

envelope_summary <- data.frame(
  Metric = c(
    "Total unique loci",
    "Pfam07253-positive loci",
    "Loci without detected Pfam07253",
    "Percentage positive",
    "Number of Pfam07253 hit records"
  ),
  Value = c(
    total_loci,
    n_envelope,
    n_without_envelope,
    100 * n_envelope / total_loci,
    nrow(envelope_hits)
  )
)

write.csv(
  envelope_summary,
  file.path(outdir, "07_Pfam07253_summary.csv"),
  row.names = FALSE
)

write.csv(
  envelope_hits,
  file.path(outdir, "08_Pfam07253_detailed_hits.csv"),
  row.names = FALSE
)

message("Total unique loci: ", total_loci)
message("Pfam07253-positive loci: ", n_envelope)

# =============================================================
# 11. PUBLICATION FIGURE A: PFAM07253 PREVALENCE
# =============================================================

prevalence_data <- data.frame(
  Category = c(
    "Pfam07253 detected",
    "Pfam07253 not detected"
  ),
  Count = c(
    n_envelope,
    n_without_envelope
  )
) %>%
  mutate(
    Percentage = 100 * Count / total_loci,
    Label = sprintf(
      "n = %d (%.1f%%)",
      Count,
      Percentage
    )
  )

p_prevalence <- ggplot(
  prevalence_data,
  aes(
    x = Count,
    y = reorder(Category, Count)
  )
) +
  geom_col(
    width = 0.62,
    fill = "grey45"
  ) +
  geom_text(
    aes(label = Label),
    hjust = -0.08,
    size = 3.8
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.30))
  ) +
  labs(
    title = "Prevalence of Pfam07253 in analysed loci",
    subtitle = paste0(
      "Unique loci analysed: N = ",
      total_loci
    ),
    x = "Number of unique loci",
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text = element_text(colour = "black"),
    axis.title.x = element_text(face = "bold"),
    plot.margin = margin(8, 25, 8, 8)
  )

ggsave(
  file.path(outdir, "09_Pfam07253_prevalence.pdf"),
  p_prevalence,
  width = 8.5,
  height = 4.5,
  device = "pdf"
)

ggsave(
  file.path(outdir, "09_Pfam07253_prevalence.png"),
  p_prevalence,
  width = 8.5,
  height = 4.5,
  dpi = PNG_DPI,
  bg = "white"
)

# =============================================================
# 12. PUBLICATION FIGURE B: DOMAIN CO-OCCURRENCE
# =============================================================

if (n_envelope > 0) {

  envelope_architecture <- locus_domains %>%
    filter(locus %in% envelope_loci) %>%
    summarise(
      across(
        all_of(domain_columns),
        ~ sum(.x, na.rm = TRUE)
      )
    ) %>%
    pivot_longer(
      cols = everything(),
      names_to = "Feature",
      values_to = "Unique_loci"
    ) %>%
    mutate(
      Percentage = 100 * Unique_loci / n_envelope,
      Description = unname(feature_labels[Feature])
    ) %>%
    arrange(Unique_loci)

  write.csv(
    envelope_architecture,
    file.path(outdir, "10_Pfam07253_domain_cooccurrence.csv"),
    row.names = FALSE
  )

  p_cooccurrence <- ggplot(
    envelope_architecture,
    aes(
      x = Unique_loci,
      y = reorder(Description, Unique_loci)
    )
  ) +
    geom_col(
      width = 0.68,
      fill = "grey45"
    ) +
    geom_text(
      aes(
        label = sprintf(
          "n = %d (%.1f%%)",
          Unique_loci,
          Percentage
        )
      ),
      hjust = -0.08,
      size = 3.1
    ) +
    scale_x_continuous(
      expand = expansion(mult = c(0, 0.35))
    ) +
    labs(
      title = "Domain co-occurrence in Pfam07253-positive loci",
      subtitle = paste0(
        "Genomic loci: ",
        n_envelope,
        " unique Pfam07253-positive loci"
      ),
      x = "Number of positive loci containing each domain",
      y = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(face = "bold")
    )

  ggsave(
    file.path(outdir, "10_Pfam07253_cooccurrence.pdf"),
    p_cooccurrence,
    width = 10,
    height = 5.5,
    device = "pdf"
  )

  ggsave(
    file.path(outdir, "10_Pfam07253_cooccurrence.png"),
    p_cooccurrence,
    width = 10,
    height = 5.5,
    dpi = PNG_DPI,
    bg = "white"
  )

} else {

  message(
    "No Pfam07253-positive loci found. ",
    "Co-occurrence plot skipped."
  )
}

# =============================================================
# 13. PUBLICATION FIGURE C: DOMAIN HEATMAP
# =============================================================

if (MAKE_HEATMAP && nrow(locus_domains) > 0) {

  heatmap_data <- locus_domains %>%
    select(
      locus,
      all_of(domain_columns)
    )

  heatmap_matrix <- as.matrix(
    heatmap_data[, domain_columns, drop = FALSE] * 1
  )

  rownames(heatmap_matrix) <- heatmap_data$locus

  # Hide individual locus labels for readability.
  # Very large datasets may produce a large PDF.
  pdf_height <- min(
    24,
    max(7, 0.12 * nrow(heatmap_matrix))
  )

  pdf(
    file.path(outdir, "11_Domain_presence_absence_heatmap.pdf"),
    width = 10,
    height = pdf_height,
    useDingbats = FALSE
  )

  pheatmap(
    heatmap_matrix,
    cluster_rows = nrow(heatmap_matrix) > 1,
    cluster_cols = FALSE,
    show_rownames = FALSE,
    border_color = NA,
    main = "Domain presence/absence across analysed loci",
    legend_breaks = c(0, 1),
    legend_labels = c("Absent", "Present")
  )

  dev.off()
}

# =============================================================
# 14. PUBLICATION FIGURE D: PFAM07253 POSITION ON ±5 KB FLANKS
# =============================================================

# This figure is valid only when:
# 1. qstart/qend are query coordinates on the 10-kb flank FASTA.
# 2. Position 5000 is the reference point.
# 3. The query represents mosquito genomic flanking sequence.
#
# If these conditions are not met, set
# PLOT_FLANK_COORDINATES <- FALSE.

if (PLOT_FLANK_COORDINATES && n_envelope > 0) {

  coordinate_hits <- envelope_hits %>%
    filter(
      !is.na(qstart),
      !is.na(qend),
      qstart > 0,
      qend > 0
    ) %>%
    mutate(
      query_start = pmin(qstart, qend),
      query_end = pmax(qstart, qend),

      relative_start = query_start - ORIGIN_POSITION,
      relative_end = query_end - ORIGIN_POSITION,

      start_kb = relative_start / 1000,
      end_kb = relative_end / 1000,

      midpoint_kb = (
        relative_start + relative_end
      ) / 2000
    ) %>%
    filter(
      relative_start >= -ORIGIN_POSITION,
      relative_end <= WINDOW_BP - ORIGIN_POSITION
    )

  if (nrow(coordinate_hits) > 0) {

    # Sort loci by the median position of their detected hits.
    locus_order <- coordinate_hits %>%
      group_by(locus) %>%
      summarise(
        median_position = median(midpoint_kb, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(median_position) %>%
      pull(locus)

    coordinate_hits$locus <- factor(
      coordinate_hits$locus,
      levels = locus_order
    )

    p_flank <- ggplot(
      coordinate_hits,
      aes(
        y = locus,
        xmin = start_kb,
        xmax = end_kb
      )
    ) +
      geom_rect(
        aes(
          ymin = as.numeric(locus) - 0.30,
          ymax = as.numeric(locus) + 0.30
        ),
        fill = "grey35",
        colour = "black",
        linewidth = 0.20
      ) +
      geom_vline(
        xintercept = 0,
        linetype = "dashed",
        linewidth = 0.5
      ) +
      scale_x_continuous(
        limits = c(-5, 5),
        breaks = seq(-5, 5, by = 1),
        expand = c(0, 0)
      ) +
      labs(
        title = "Pfam07253 hits in mosquito genomic flanks",
        subtitle = "Relative to the configured reference point",
        x = "Position relative to reference point (kb)",
        y = "Mosquito EVE locus"
      ) +
      theme_classic(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold"),
        axis.text = element_text(colour = "black"),
        axis.title = element_text(face = "bold")
      )

    plot_height <- min(
      24,
      max(4.5, 0.28 * n_distinct(coordinate_hits$locus))
    )

    ggsave(
      file.path(outdir, "12_Pfam07253_flank_positions.pdf"),
      p_flank,
      width = 10,
      height = plot_height,
      device = "pdf",
      limitsize = FALSE
    )

    ggsave(
      file.path(outdir, "12_Pfam07253_flank_positions.png"),
      p_flank,
      width = 10,
      height = plot_height,
      dpi = PNG_DPI,
      bg = "white",
      limitsize = FALSE
    )

    write.csv(
      coordinate_hits %>%
        select(
          locus,
          hit_id,
          accession,
          domain,
          description,
          pident,
          evalue,
          bitscore,
          query_start,
          query_end,
          relative_start,
          relative_end,
          start_kb,
          end_kb,
          midpoint_kb
        ),
      file.path(outdir, "12_Pfam07253_flank_coordinates.csv"),
      row.names = FALSE
    )

  } else {

    message(
      "No valid Pfam07253 coordinates were available ",
      "for the flank plot."
    )
  }
}

# =============================================================
# 15. EXPORT A SHORT ANALYSIS REPORT
# =============================================================

report <- c(
  "MOSQUITO EVE Pfam07253 ANALYSIS",
  "============================================",
  paste("Input:", basename(input_file)),
  paste("Total hit records:", total_records),
  paste("Unique query/locus identifiers:", total_loci),
  paste("Pfam07253-positive loci:", n_envelope),
  paste(
    "Percentage positive:",
    sprintf("%.2f%%", 100 * n_envelope / total_loci)
  ),
  "",
  "Interpretation notes:",
  "1. Pfam07253 is treated as a putative envelope-associated domain.",
  "2. Pfam07253 and Baculo_F/Pfam12259 are classified separately.",
  "3. Domain co-occurrence does not prove that domains belong to",
  "   a single functional protein or intact viral genome.",
  "4. Domain detection alone does not confirm viral ancestry or EVE status.",
  "5. Confirm the genomic integration context and compare candidate",
  "   sequences with viral and non-viral transposable-element references.",
  "6. The ±5 kb plot assumes qstart/qend are coordinates on the",
  "   10-kb mosquito flank sequence and that the reference point is correct."
)

writeLines(
  report,
  file.path(outdir, "13_Analysis_report.txt")
)

# =============================================================
# 16. COMPLETION
# =============================================================

message("")
message("============================================")
message("ANALYSIS COMPLETED")
message("============================================")
message("Unique loci: ", total_loci)
message("Pfam07253-positive loci: ", n_envelope)
message(
  "Percentage positive: ",
  sprintf("%.2f%%", 100 * n_envelope / total_loci)
)
message("Results saved in: ", normalizePath(outdir))
message("============================================")

######################################################################33
# =============================================================
# SECTION 12B. DOMAIN PREVALENCE ACROSS ALL LOCI
# =============================================================

message("Generating domain prevalence across all loci...")

# Count each domain class once per unique locus.
domain_all <- locus_domains %>%
  summarise(
    across(
      all_of(domain_columns),
      ~ sum(.x, na.rm = TRUE)
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Feature",
    values_to = "Unique_loci"
  ) %>%
  mutate(
    Total_loci = nrow(locus_domains),
    Percent_loci = 100 * Unique_loci / Total_loci,
    Description = unname(feature_labels[Feature])
  ) %>%
  arrange(Unique_loci)

write.csv(
  domain_all,
  file.path(outdir, "14_All_loci_domain_prevalence.csv"),
  row.names = FALSE
)

# Publication-quality horizontal bar plot.
p_all_domains <- ggplot(
  domain_all,
  aes(
    x = Unique_loci,
    y = reorder(Description, Unique_loci),
    fill = Feature == "Envelope"
  )
) +
  geom_col(width = 0.68) +
  geom_text(
    aes(
      label = sprintf(
        "n = %s (%.1f%%)",
        format(Unique_loci, big.mark = ","),
        Percent_loci
      )
    ),
    hjust = -0.08,
    size = 3.4
  ) +
  scale_fill_manual(
    values = c(
      "FALSE" = "grey55",
      "TRUE" = "#B2182B"
    ),
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(
      0,
      max(domain_all$Unique_loci) * 1.32
    ),
    expand = c(0, 0)
  ) +
  labs(
    subtitle = paste0(
      "Denominator: ",
      nrow(locus_domains),
      " unique loci"
    ),
    x = "Number of unique loci containing each domain",
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text = element_text(colour = "black"),
    axis.title.x = element_text(face = "bold"),
    plot.margin = margin(10, 40, 10, 10)
  )

ggsave(
  file.path(outdir, "14_All_loci_domain_prevalence.pdf"),
  p_all_domains,
  width = 11,
  height = 6.5,
  device = "pdf"
)

ggsave(
  file.path(outdir, "14_All_loci_domain_prevalence.png"),
  p_all_domains,
  width = 11,
  height = 6.5,
  dpi = PNG_DPI,
  bg = "white"
)


# =============================================================
# SECTION 12C. PFAM07253-POSITIVE VS NEGATIVE LOCI
# =============================================================

message("Comparing Pfam07253-positive and negative loci...")

# Label each unique locus by Pfam07253 status.
locus_groups <- locus_domains %>%
  mutate(
    Pfam07253_status = if_else(
      Envelope,
      "Pfam07253-positive",
      "Pfam07253-negative"
    )
  )

# Calculate the prevalence of every domain class in each group.
group_summary <- locus_groups %>%
  group_by(Pfam07253_status) %>%
  summarise(
    across(
      all_of(domain_columns),
      ~ sum(.x, na.rm = TRUE)
    ),
    Group_size = n(),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = all_of(domain_columns),
    names_to = "Feature",
    values_to = "Unique_loci"
  ) %>%
  mutate(
    Percent = 100 * Unique_loci / Group_size,
    Description = unname(feature_labels[Feature])
  )

write.csv(
  group_summary,
  file.path(outdir, "15_Pfam07253_positive_vs_negative.csv"),
  row.names = FALSE
)

# Grouped bar plot: compare percentages rather than raw counts.
p_group_comparison <- ggplot(
  group_summary,
  aes(
    x = Percent,
    y = reorder(Description, Percent),
    fill = Pfam07253_status
  )
) +
  geom_col(
    position = position_dodge(width = 0.78),
    width = 0.68
  ) +
  labs(
    title = "Domain prevalence by Pfam07253 status",
    subtitle = paste0(
      "Positive loci: ", n_envelope,
      "; negative loci: ", total_loci - n_envelope
    ),
    x = "Loci containing domain (%)",
    y = NULL,
    fill = "Locus group"
  ) +
  scale_fill_manual(
    values = c(
      "Pfam07253-positive" = "#B2182B",
      "Pfam07253-negative" = "grey65"
    )
  ) +
  scale_x_continuous(
    limits = c(0, 105),
    breaks = seq(0, 100, 20),
    expand = c(0, 0)
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text = element_text(colour = "black"),
    axis.title.x = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  file.path(outdir, "15_Pfam07253_group_comparison.pdf"),
  p_group_comparison,
  width = 11,
  height = 6.5,
  device = "pdf"
)

ggsave(
  file.path(outdir, "15_Pfam07253_group_comparison.png"),
  p_group_comparison,
  width = 11,
  height = 6.5,
  dpi = PNG_DPI,
  bg = "white"
)


# =============================================================
# SECTION 12D. PREVALENCE OF ALL INDIVIDUAL DOMAIN LABELS
# =============================================================

message("Generating prevalence plot for all individual domains...")

# Count each domain label once per locus.
# This differs from the nine manually defined domain classes.
all_pfam_frequency <- dat2 %>%
  filter(
    !is.na(domain),
    domain != ""
  ) %>%
  distinct(locus, domain) %>%
  count(domain, name = "Unique_loci") %>%
  mutate(
    Total_loci = total_loci,
    Percent_loci = 100 * Unique_loci / Total_loci
  ) %>%
  arrange(desc(Unique_loci), domain)

write.csv(
  all_pfam_frequency,
  file.path(outdir, "16_All_individual_domains.csv"),
  row.names = FALSE
)

# Plot the 30 most prevalent individual domain labels.
top_pfam <- all_pfam_frequency %>%
  slice_head(n = min(30, nrow(all_pfam_frequency)))

if (nrow(top_pfam) > 0) {

  p_all_pfam <- ggplot(
    top_pfam,
    aes(
      x = Unique_loci,
      y = reorder(domain, Unique_loci),
      fill = str_detect(
        str_to_lower(domain),
        "pfam.?07253|envelope"
      )
    )
  ) +
    geom_col(width = 0.70) +
    geom_text(
      aes(
        label = sprintf(
          "%s (%.2f%%)",
          format(Unique_loci, big.mark = ","),
          Percent_loci
        )
      ),
      hjust = -0.08,
      size = 3.1
    ) +
    scale_fill_manual(
      values = c(
        "FALSE" = "grey55",
        "TRUE" = "#B2182B"
      ),
      guide = "none"
    ) +
    scale_x_continuous(
      limits = c(
        0,
        max(top_pfam$Unique_loci) * 1.35
      ),
      expand = c(0, 0)
    ) +
    labs(
      title = "Most prevalent individual domain annotations",
      subtitle = paste0(
        "Top ", nrow(top_pfam),
        " domain labels across ", total_loci,
        " unique loci"
      ),
      x = "Number of unique loci containing each domain",
      y = "Domain annotation"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(face = "bold"),
      plot.margin = margin(10, 40, 10, 10)
    )

  ggsave(
    file.path(outdir, "16_All_individual_domains_top30.pdf"),
    p_all_pfam,
    width = 11,
    height = 8,
    device = "pdf"
  )

  ggsave(
    file.path(outdir, "16_All_individual_domains_top30.png"),
    p_all_pfam,
    width = 11,
    height = 8,
    dpi = PNG_DPI,
    bg = "white"
  )
}


# =============================================================
# SECTION 12E. DOMAIN ASSOCIATION TESTS
# =============================================================

message("Testing associations between domains and Pfam07253 status...")

# Fisher's exact test for each domain class, excluding Envelope itself.
# Envelope cannot be tested against its own positive/negative definition.
test_features <- setdiff(domain_columns, "Envelope")

association_results <- lapply(
  test_features,
  function(feature) {

    positive_has <- sum(
      locus_groups$Pfam07253_status == "Pfam07253-positive" &
        locus_groups[[feature]],
      na.rm = TRUE
    )

    positive_lacks <- sum(
      locus_groups$Pfam07253_status == "Pfam07253-positive" &
        !locus_groups[[feature]],
      na.rm = TRUE
    )

    negative_has <- sum(
      locus_groups$Pfam07253_status == "Pfam07253-negative" &
        locus_groups[[feature]],
      na.rm = TRUE
    )

    negative_lacks <- sum(
      locus_groups$Pfam07253_status == "Pfam07253-negative" &
        !locus_groups[[feature]],
      na.rm = TRUE
    )

    contingency <- matrix(
      c(
        positive_has,
        positive_lacks,
        negative_has,
        negative_lacks
      ),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(
        Group = c("Positive", "Negative"),
        Domain = c("Present", "Absent")
      )
    )

    ft <- fisher.test(contingency)

    data.frame(
      Feature = feature,
      Description = unname(feature_labels[feature]),
      Positive_with_domain = positive_has,
      Positive_without_domain = positive_lacks,
      Negative_with_domain = negative_has,
      Negative_without_domain = negative_lacks,
      Odds_ratio = unname(ft$estimate),
      CI_lower = ft$conf.int[1],
      CI_upper = ft$conf.int[2],
      P_value = ft$p.value
    )
  }
)

association_results <- bind_rows(association_results) %>%
  mutate(
    FDR = p.adjust(P_value, method = "BH")
  ) %>%
  arrange(FDR)

write.csv(
  association_results,
  file.path(outdir, "17_Pfam07253_domain_association_tests.csv"),
  row.names = FALSE
)

message("Additional analyses completed.")
