# =============================================================================
# PROJECT      : Genome-wide Characterization of Endogenous Viral Elements 
#                Reveals Extensive Viral Diversity in Culicinae Mosquitoes
# PUBLICATION  : https://doi.org/10.21203/rs.3.rs-10266412/v1
# SCRIPT       : LG EVE Reference Selection and Evolutionary Distance Analysis
#
# CODE AUTHOR  : Ankush Kumar Sharma
# GITHUB       : https://github.com/Ankush-Kumar-Sharma
#
# CONTRIBUTIONS
# Development and implementation of the R-based analytical workflow for
# evolutionary-distance assessment, reference-sequence evaluation, statistical
# summarisation, and publication-quality visualisation.
#
# PROJECT CONTEXT
# Developed as part of a collaborative research project. Project-level
# attribution, licensing, and reuse are subject to the repository's licence
# and applicable institutional and project agreements.
#
# CITATION
# Please cite the associated publication when using this workflow.
# =============================================================================

#!/usr/bin/env Rscript

# Publication scatter plot of reference sequences against the two EVE groups.
# Each point is one reference sequence:
#   x = median LG+G distance to TED-homologous EVEs
#   y = median LG+G distance to CladoV-homologous EVEs
#
# Point colour is taken from Comparision_Meta_Colour in the metadata.
# Point shape distinguishes LTR retroelement references from extant-virus references.
# Red cutoff lines are defined independently on each axis by the maximum median
# distance observed among references annotated as Retroviridae. All finite points
# remain visible; the lines mark the proposed selection boundary and do not delete
# or hide observations.

options(stringsAsFactors = FALSE)

MATRIX_FILE <- file.path("mds", "results_mds", "LG_pairwise_distance_matrix.tsv")
FASTA_FILE <- "350andabove_nr_tedclado_linsi_mafftadd_trimmed_python.fasta"
METADATA_FILE <- "350andabove_nr_tedclado_metadata_reftaxa_final.tsv"
OUTDIR <- file.path("mds", "results_mds")
OUTPUT_FILE <- file.path(OUTDIR, "LG_EVE_reference_scatter.svg")
SUMMARY_FILE <- file.path(OUTDIR, "LG_EVE_reference_category_distance_summary.tsv")

EVE_CLASSES <- c("TED_EVE", "CladoV_EVE")
REFERENCE_CLASSES <- c("LTR-Retroelements-pol", "viral_ref")
RETROVIRUS_CATEGORY <- "Retroviridae"
EVE_X_CLASS <- "TED_EVE"
EVE_Y_CLASS <- "CladoV_EVE"
UNASSIGNED_LABEL <- "Unassigned viruses"
FALLBACK_COLOUR <- "#999999"

if (!file.exists(MATRIX_FILE)) stop("Distance matrix not found: ", MATRIX_FILE)
if (!file.exists(FASTA_FILE)) stop("Alignment FASTA not found: ", FASTA_FILE)
if (!file.exists(METADATA_FILE)) stop("Metadata file not found: ", METADATA_FILE)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

clean <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x %in% c("", "NA")] <- NA_character_
  x
}

fasta_lines <- readLines(FASTA_FILE, warn = FALSE)
fasta_ids <- substring(fasta_lines[startsWith(fasta_lines, ">")], 2L)
if (!length(fasta_ids)) stop("No FASTA headers found.")
if (anyDuplicated(fasta_ids)) stop("FASTA headers are not unique.")

meta <- read.delim(
  METADATA_FILE, header = TRUE, sep = "\t", quote = "",
  check.names = FALSE, stringsAsFactors = FALSE
)
names(meta) <- trimws(names(meta))

required_meta <- c("sequence_class", "Comparision_Meta_Category",
                   "Comparision_Meta_Colour")
if (!all(required_meta %in% names(meta))) {
  stop("Metadata is missing required columns: ",
       paste(setdiff(required_meta, names(meta)), collapse = ", "))
}

id_candidates <- intersect(
  c("sequence_id", "original_header", "processed_header", "sequence_name"),
  names(meta)
)
if (!length(id_candidates)) {
  stop("Metadata must contain an identifier column matching the FASTA headers.")
}

meta_index <- NULL
meta_id_column <- NULL
for (nm in id_candidates) {
  idx <- match(fasta_ids, trimws(as.character(meta[[nm]])))
  if (all(!is.na(idx))) {
    meta_index <- idx
    meta_id_column <- nm
    break
  }
}
if (is.null(meta_index)) {
  stop("Could not match every FASTA header to metadata using: ",
       paste(id_candidates, collapse = ", "))
}
meta <- meta[meta_index, , drop = FALSE]

meta$sequence_class <- clean(meta$sequence_class)
meta$lineage <- clean(meta$Comparision_Meta_Category)
meta$lineage[is.na(meta$lineage)] <- "Unlabelled"
meta$lineage[meta$lineage == "Unassigned"] <- UNASSIGNED_LABEL
meta$metadata_colour <- clean(meta$Comparision_Meta_Colour)

valid_colour <- !is.na(meta$metadata_colour) &
  grepl("^#[[:xdigit:]]{6}$", meta$metadata_colour)
colour_map <- tapply(
  meta$metadata_colour[valid_colour],
  meta$lineage[valid_colour],
  function(x) x[[1L]]
)

D <- as.matrix(read.delim(
  MATRIX_FILE, header = TRUE, sep = "\t", quote = "",
  row.names = 1, check.names = FALSE, na.strings = c("NA", "")
))
storage.mode(D) <- "numeric"

if (is.null(rownames(D)) || is.null(colnames(D))) {
  stop("Distance matrix must contain sequence identifiers as row and column names.")
}
if (nrow(D) != ncol(D)) stop("Distance matrix is not square.")
if (!all(fasta_ids %in% rownames(D)) || !all(fasta_ids %in% colnames(D))) {
  stop("FASTA identifiers do not match the distance-matrix row/column names.")
}
D <- D[fasta_ids, fasta_ids, drop = FALSE]

eve_idx <- which(meta$sequence_class %in% EVE_CLASSES)
ref_idx <- which(meta$sequence_class %in% REFERENCE_CLASSES)
if (!length(eve_idx)) stop("No TED_EVE or CladoV_EVE sequences were found.")
if (!any(meta$sequence_class[eve_idx] == EVE_X_CLASS)) stop("No TED_EVE sequences found.")
if (!any(meta$sequence_class[eve_idx] == EVE_Y_CLASS)) stop("No CladoV_EVE sequences found.")
if (!length(ref_idx)) stop("No reference sequences were found.")

ted_idx <- eve_idx[meta$sequence_class[eve_idx] == EVE_X_CLASS]
cladov_idx <- eve_idx[meta$sequence_class[eve_idx] == EVE_Y_CLASS]

median_finite <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else stats::median(x)
}

points <- lapply(ref_idx, function(i) {
  dx <- D[i, ted_idx]
  dy <- D[i, cladov_idx]
  data.frame(
    reference_id = fasta_ids[i],
    reference_type = if (meta$sequence_class[i] == "viral_ref")
      "Extant virus" else "LTR retroelement",
    lineage = meta$lineage[i],
    colour = unname(colour_map[meta$lineage[i]]),
    x = median_finite(dx),
    y = median_finite(dy),
    n_TED_pairs = sum(is.finite(dx)),
    n_CladoV_pairs = sum(is.finite(dy)),
    stringsAsFactors = FALSE
  )
})
points <- do.call(rbind, points)
points$colour[is.na(points$colour) | !grepl("^#[[:xdigit:]]{6}$", points$colour)] <-
  FALLBACK_COLOUR

n_before <- nrow(points)
points <- points[is.finite(points$x) & is.finite(points$y), , drop = FALSE]
if (!nrow(points)) stop("No reference sequences have finite medians on both axes.")

retro <- points[tolower(points$lineage) == tolower(RETROVIRUS_CATEGORY), , drop = FALSE]
if (!nrow(retro)) {
  stop("No reference points are labelled '", RETROVIRUS_CATEGORY,
       "' in Comparision_Meta_Category; cannot define the requested cutoffs.")
}

# Transparent, data-derived reference envelope: the maximum of the observed
# Retroviridae median distances on each axis. All points, including points beyond
# the lines, are drawn. Category-wise Q1, median, and Q3 are also reported below.
x_cutoff <- max(retro$x, na.rm = TRUE)
y_cutoff <- max(retro$y, na.rm = TRUE)
points$within_cutoff <- points$x <= x_cutoff & points$y <= y_cutoff

# Summarise observed distance distributions by metadata category. These summaries
# support interpretation of the selected reference range; Q3 is not imposed as a
# preselected cutoff.
summarise_axis <- function(z) {
  z <- z[is.finite(z)]
  if (!length(z)) return(c(n = 0, Q1 = NA_real_, median = NA_real_,
                           Q3 = NA_real_, maximum = NA_real_))
  c(n = length(z),
    Q1 = unname(stats::quantile(z, 0.25, na.rm = TRUE, names = FALSE, type = 7)),
    median = stats::median(z),
    Q3 = unname(stats::quantile(z, 0.75, na.rm = TRUE, names = FALSE, type = 7)),
    maximum = max(z))
}
category_summaries <- lapply(split(points, points$lineage), function(d) {
  sx <- summarise_axis(d$x)
  sy <- summarise_axis(d$y)
  data.frame(
    Comparision_Meta_Category = d$lineage[1],
    n_reference_sequences = nrow(d),
    TED_n = unname(sx["n"]),
    TED_Q1 = unname(sx["Q1"]),
    TED_median = unname(sx["median"]),
    TED_Q3 = unname(sx["Q3"]),
    TED_maximum = unname(sx["maximum"]),
    CladoV_n = unname(sy["n"]),
    CladoV_Q1 = unname(sy["Q1"]),
    CladoV_median = unname(sy["median"]),
    CladoV_Q3 = unname(sy["Q3"]),
    CladoV_maximum = unname(sy["maximum"]),
    stringsAsFactors = FALSE
  )
})
category_summaries <- do.call(rbind, category_summaries)
category_summaries <- category_summaries[
  order(tolower(category_summaries$Comparision_Meta_Category)), , drop = FALSE
]
write.table(category_summaries, SUMMARY_FILE, sep = "\\t",
            quote = FALSE, row.names = FALSE, na = "NA")

# Assign shapes by reference type while retaining metadata-derived point colours.
points$pch <- ifelse(points$reference_type == "Extant virus", 24L, 21L)

# A small amount of padding ensures that all data remain visible.
x_range <- range(c(0, points$x, x_cutoff), finite = TRUE)
y_range <- range(c(0, points$y, y_cutoff), finite = TRUE)
x_pad <- max(diff(x_range) * 0.06, 0.05)
y_pad <- max(diff(y_range) * 0.06, 0.05)
x_lim <- c(max(0, x_range[1L] - x_pad), x_range[2L] + x_pad)
y_lim <- c(max(0, y_range[1L] - y_pad), y_range[2L] + y_pad)

# Publication layout: main scatter panel on the left; a metadata-colour key
# made of coloured tiles on the right. No individual reference names are plotted.
grDevices::svg(
  OUTPUT_FILE, width = 12.2, height = 8.0,
  pointsize = 10, onefile = TRUE, family = "sans"
)
old_par <- par(no.readonly = TRUE)
layout(matrix(c(1, 2), nrow = 1), widths = c(2.65, 2.35))

# Main scatter panel
par(mar = c(5.1, 5.2, 3.1, 0.8), mgp = c(3.0, 0.9, 0),
    tcl = -0.25, las = 1, bty = "l")
plot(
  NA, xlim = x_lim, ylim = y_lim, axes = FALSE,
  xlab = "Median LG+G distance to TED-homologous EVEs",
  ylab = "Median LG+G distance to CladoV-homologous EVEs",
  main = "Reference-sequence distance profiles"
)
x_ticks <- pretty(x_lim, n = 5)
y_ticks <- pretty(y_lim, n = 5)
axis(1, at = x_ticks, labels = format(x_ticks, trim = TRUE, nsmall = 1), cex.axis = 0.9)
axis(2, at = y_ticks, labels = format(y_ticks, trim = TRUE, nsmall = 1), las = 1, cex.axis = 0.9)
box()
abline(v = x_cutoff, h = y_cutoff, col = "#D62728", lwd = 1.7, lty = 2)
points(
  points$x, points$y, pch = points$pch,
  bg = points$colour, col = "grey15", lwd = 0.55, cex = 1.15
)
legend(
  "topleft",
  legend = c("LTR retroelement", "Extant virus"),
  pch = c(21, 24), pt.bg = "grey85", col = "grey15",
  pt.cex = 1.1, bty = "n", cex = 0.8, title = "Reference type"
)
# Numeric cutoff labels are explicit so the boundary is interpretable in print.
text(x = x_cutoff, y = y_lim[2] - 0.035 * diff(y_lim),
     labels = sprintf("TED cutoff = %.2f", x_cutoff),
     col = "#D62728", cex = 0.76, pos = 4, xpd = NA)
text(x = x_lim[2] - 0.02 * diff(x_lim), y = y_cutoff,
     labels = sprintf("CladoV cutoff = %.2f", y_cutoff),
     col = "#D62728", cex = 0.76, pos = 3, xpd = NA, adj = 1)
mtext(
  "Selection rationale: closer references ~1.9–3.0; Retroviridae extend to ~4.5; more distant groups often ≥5.",
  side = 3, line = 0.45, adj = 0, cex = 0.62
)

# Separate legend panel: tiles are coloured with Comparision_Meta_Colour and
# labels are the exact metadata category values. Only categories represented
# among plotted reference points are shown.
par(mar = c(1.0, 0.2, 2.5, 0.3), xaxs = "i", yaxs = "i")
legend_categories <- unique(points$lineage)
legend_categories <- legend_categories[order(tolower(legend_categories))]
legend_colours <- vapply(legend_categories, function(cat) {
  cc <- points$colour[match(cat, points$lineage)]
  if (length(cc) && !is.na(cc[1])) cc[1] else FALLBACK_COLOUR
}, character(1))

plot(
  NA, xlim = c(0, 2.65), ylim = c(0, max(1, length(legend_categories) + 1)),
  axes = FALSE, xlab = "", ylab = "", bty = "n",
  main = "Reference categories"
)
text(0.02, max(1, length(legend_categories) + 0.45),
     labels = "Tile colour follows metadata: Comparision_Meta_Colour",
     adj = 0, cex = 0.65, xpd = NA)
for (i in seq_along(legend_categories)) {
  yy <- length(legend_categories) - i + 1
  rect(0.02, yy - 0.22, 0.18, yy + 0.22,
       col = legend_colours[i], border = "grey25", lwd = 0.5)
  text(0.24, yy, labels = legend_categories[i], adj = 0,
       cex = 0.72, xpd = NA)
}
par(old_par)
invisible(dev.off())

cat("Wrote: ", OUTPUT_FILE, "\n", sep = "")
cat("References with finite medians on both axes: ", nrow(points), "\n", sep = "")
cat("References with unavailable medians omitted from plot: ", n_before - nrow(points), "\n", sep = "")
cat("Retroviridae-derived x cutoff (TED): ", signif(x_cutoff, 6), "\n", sep = "")
cat("Retroviridae-derived y cutoff (CladoV): ", signif(y_cutoff, 6), "\n", sep = "")
cat("References within both cutoffs: ", sum(points$within_cutoff), "\n", sep = "")
cat("References outside at least one cutoff: ", sum(!points$within_cutoff), "\n", sep = "")
cat("All finite reference points remain visible; cutoff lines do not remove points.\n")
