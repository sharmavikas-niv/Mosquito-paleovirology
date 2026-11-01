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

FASTA <- "350andabove_nr_tedclado_linsi_mafftadd_trimmed_python.fasta"
METADATA_FILE <- "350andabove_nr_tedclado_metadata_reftaxa_final.tsv"
OUTDIR <- file.path("mds", "results_mds")
MATRIX_FILE <- file.path(OUTDIR, "LG_pairwise_distance_matrix.tsv")

EVE_TYPES <- c("TED_EVE", "CladoV_EVE")
REFERENCE_TYPES <- c("LTR-Retroelements-pol", "viral_ref")
EVE_LABELS <- c(
  TED_EVE = "TED-homologous EVE",
  CladoV_EVE = "CladoV-homologous EVE"
)
FALLBACK_COLOUR <- "#B8B8B8"

if (!file.exists(MATRIX_FILE)) stop("Distance matrix not found: ", MATRIX_FILE)
if (!file.exists(FASTA)) stop("Alignment file not found: ", FASTA)
if (!file.exists(METADATA_FILE)) stop("Metadata file not found: ", METADATA_FILE)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

clean <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x %in% c("", "NA")] <- NA_character_
  x
}

fasta_headers <- {
  x <- readLines(FASTA, warn = FALSE)
  substring(x[startsWith(x, ">")], 2L)
}

metadata <- read.delim(
  METADATA_FILE, sep = "\t", header = TRUE, check.names = FALSE,
  stringsAsFactors = FALSE
)
names(metadata) <- trimws(names(metadata))

required <- c("sequence_class", "Comparision_Meta_Category",
              "Comparision_Meta_Colour")
if (!all(required %in% names(metadata))) {
  stop("Metadata is missing required columns: ",
       paste(setdiff(required, names(metadata)), collapse = ", "))
}

id_candidates <- intersect(
  c("sequence_id", "original_header", "processed_header", "sequence_name"),
  names(metadata)
)
if (!length(id_candidates)) stop("Metadata requires a sequence ID column.")

metadata_index <- NULL
for (nm in id_candidates) {
  idx <- match(fasta_headers, trimws(as.character(metadata[[nm]])))
  if (all(!is.na(idx))) {
    metadata_index <- idx
    break
  }
}
if (is.null(metadata_index)) stop("Every FASTA header must match a metadata ID.")
metadata <- metadata[metadata_index, , drop = FALSE]

sequence_class <- clean(metadata$sequence_class)
category <- clean(metadata$Comparision_Meta_Category)
category[is.na(category)] <- "Unlabelled"
category[category == "Unassigned"] <- "Unassigned viruses"
metadata_colour <- clean(metadata$Comparision_Meta_Colour)

valid_colour <- !is.na(metadata_colour) &
  grepl("^#([[:xdigit:]]{6}|[[:xdigit:]]{8})$", metadata_colour)
metadata_colour[!valid_colour] <- FALLBACK_COLOUR

# Use the metadata's category colour consistently. If a category has multiple
# colours, use its most frequent valid colour and report the inconsistency.
category_colours <- split(metadata_colour, category)
category_colours <- vapply(category_colours, function(x) {
  tab <- sort(table(x), decreasing = TRUE)
  names(tab)[1L]
}, character(1))

D <- as.matrix(read.delim(
  MATRIX_FILE, sep = "\t", header = TRUE, row.names = 1,
  check.names = FALSE, na.strings = "NA"
))
storage.mode(D) <- "numeric"
if (is.null(rownames(D)) || is.null(colnames(D)))
  stop("Distance matrix must have sequence names in rows and columns.")
if (!all(fasta_headers %in% rownames(D)) ||
    !all(fasta_headers %in% colnames(D))) {
  stop("Distance matrix sequence names do not match the alignment.")
}
D <- D[fasta_headers, fasta_headers, drop = FALSE]
if (nrow(D) != ncol(D)) stop("Distance matrix must be square.")

is_eve <- sequence_class %in% EVE_TYPES
ref_idx <- which(sequence_class %in% REFERENCE_TYPES)
eve_idx <- which(is_eve)
if (!length(eve_idx) || !length(ref_idx))
  stop("The alignment must contain EVE and reference sequences.")

ref_type <- ifelse(
  sequence_class[ref_idx] == "LTR-Retroelements-pol",
  "LTR retroelement", "Extant virus"
)
ref_category <- category[ref_idx]

# Build EVE-to-reference pair data from the saved distance matrix.
plot_rows <- list()
z <- 0L
for (r in ref_idx) {
  for (e in eve_idx) {
    value <- D[e, r]
    if (!is.finite(value)) next
    z <- z + 1L
    plot_rows[[z]] <- data.frame(
      eve_id = fasta_headers[e],
      eve_class = sequence_class[e],
      eve_category = category[e],
      eve_colour = metadata_colour[e],
      reference_id = fasta_headers[r],
      reference_lineage = category[r],
      reference_colour = metadata_colour[r],
      reference_type = ref_type[match(r, ref_idx)],
      distance = value,
      stringsAsFactors = FALSE
    )
  }
}
if (!length(plot_rows))
  stop("No finite EVE-to-reference distances are available.")
pairs <- do.call(rbind, plot_rows)

pairs$group_label <- paste(
  pairs$reference_lineage, pairs$reference_type,
  unname(EVE_LABELS[pairs$eve_class]), sep = " | "
)
group_medians <- aggregate(distance ~ group_label, pairs, median)
group_medians <- group_medians[
  order(group_medians$distance, group_medians$group_label), , drop = FALSE
]
group_order <- group_medians$group_label
pairs$group_label <- factor(pairs$group_label, levels = group_order)

# Figure 1: all unique pairwise distances, stacked by broad pair type.
# EVE-to-reference distances are highlighted using the EVE category colours.
all_i <- which(upper.tri(D), arr.ind = TRUE)
all_values <- D[all_i]
finite <- is.finite(all_values)
all_i <- all_i[finite, , drop = FALSE]
all_values <- all_values[finite]

if (length(all_values)) {
  pair_group <- rep("Other pairs", length(all_values))
  ci <- sequence_class[all_i[, 1L]]
  cj <- sequence_class[all_i[, 2L]]
  eve_i <- ci %in% EVE_TYPES
  eve_j <- cj %in% EVE_TYPES
  ref_i <- ci %in% REFERENCE_TYPES
  ref_j <- cj %in% REFERENCE_TYPES

  pair_group[eve_i & eve_j] <- "EVE–EVE"
  pair_group[eve_i & ref_j] <- paste0("EVE–reference: ", ci[eve_i & ref_j])
  pair_group[ref_i & eve_j] <- paste0("EVE–reference: ", cj[ref_i & eve_j])
  pair_group[ref_i & ref_j] <- "Reference–reference"

  pair_levels <- c(
    "Other pairs", "Reference–reference", "EVE–EVE",
    paste0("EVE–reference: ", EVE_TYPES)
  )
  pair_levels <- pair_levels[pair_levels %in% pair_group]
  pair_colours <- c(
    "Other pairs" = "#D9D9D9",
    "Reference–reference" = "#969696",
    "EVE–EVE" = "#636363"
  )
  for (eve_type in EVE_TYPES) {
    eve_rows <- which(sequence_class == eve_type)
    eve_cat <- if (length(eve_rows)) category[eve_rows[1L]] else NA_character_
    eve_col <- if (!is.na(eve_cat) && eve_cat %in% names(category_colours)) {
      category_colours[[eve_cat]]
    } else {
      FALLBACK_COLOUR
    }
    pair_colours[paste0("EVE–reference: ", eve_type)] <- eve_col
  }

  breaks <- hist(all_values, breaks = "FD", plot = FALSE)$breaks
  bin_id <- cut(all_values, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  counts <- table(
    factor(bin_id, levels = seq_len(length(breaks) - 1L)),
    factor(pair_group, levels = pair_levels)
  )
  grDevices::svg(
    file.path(OUTDIR, "LG_distance_distribution.svg"),
    width = 9.5, height = 6.2, pointsize = 11, onefile = TRUE, family = "sans"
  )
  old_par <- par(no.readonly = TRUE)
  par(mar = c(4.5, 4.7, 3.2, 1), mgp = c(2.6, 0.8, 0),
      tcl = -0.25, las = 1, bty = "l")
  mids <- (breaks[-1L] + breaks[-length(breaks)]) / 2
  barplot(
    t(as.matrix(counts)), beside = FALSE, space = 0,
    names.arg = rep("", length(mids)),
    col = pair_colours[pair_levels], border = NA,
    xlab = "LG+G evolutionary distance (amino-acid substitutions per site)",
    ylab = "Number of unique sequence pairs",
    main = "Distribution of pairwise LG+G evolutionary distances",
    xaxt = "n"
  )
  tick_positions <- pretty(seq_along(mids))
  tick_positions <- tick_positions[
    tick_positions >= 1 & tick_positions <= length(mids)
  ]
  tick_labels <- format(mids[pmax(1, pmin(length(mids), round(tick_positions)))],
                        trim = TRUE, digits = 3)
  axis(1, at = tick_positions, labels = tick_labels)
  legend(
    "topright", legend = pair_levels, fill = pair_colours[pair_levels],
    border = NA, bty = "n", cex = 0.72
  )
  mtext(paste0("Finite unique pairs: ", format(length(all_values), big.mark = ",")),
        side = 3, line = 0.4, adj = 1, cex = 0.8)
  par(old_par)
  invisible(dev.off())
}

# Figure 2: box fills use EVE category colours; x-axis labels use reference
# category colours. Groups are sorted by increasing median distance and medians
# are printed numerically above their boxes.
group_levels <- levels(pairs$group_label)
box_data <- lapply(group_levels, function(g) pairs$distance[pairs$group_label == g])
medians <- vapply(box_data, median, numeric(1))

eve_colour_for_group <- vapply(group_levels, function(g) {
  row <- pairs[as.character(pairs$group_label) == g, , drop = FALSE]
  category_colours[[row$eve_category[1L]]]
}, character(1))
ref_colour_for_group <- vapply(group_levels, function(g) {
  row <- pairs[as.character(pairs$group_label) == g, , drop = FALSE]
  category_colours[[row$reference_lineage[1L]]]
}, character(1))

labels <- vapply(group_levels, function(g) {
  bits <- strsplit(g, " \\| ", perl = TRUE)[[1]]
  paste(bits, collapse = "\n")
}, character(1))
median_labels <- format(signif(medians, 3), trim = TRUE, scientific = FALSE)
plot_width <- max(12, min(32, 8 + 0.32 * length(group_levels)))

grDevices::svg(
  file.path(OUTDIR, "LG_EVE_to_reference_distances.svg"),
  width = plot_width, height = 9, pointsize = 10, onefile = TRUE, family = "sans"
)
old_par <- par(no.readonly = TRUE)
par(mar = c(13, 5, 3, 2), mgp = c(2.8, 0.8, 0),
    tcl = -0.25, las = 1, bty = "l", cex.axis = 0.75)
boxplot(
  box_data, at = seq_along(box_data), names = FALSE, outline = FALSE,
  col = eve_colour_for_group, border = "grey25",
  ylab = "LG+G evolutionary distance (amino-acid substitutions per site)",
  xlab = "", main = "EVE-to-reference evolutionary distances",
  ylim = c(0, max(c(unlist(box_data), medians), na.rm = TRUE) * 1.14)
)
axis(1, at = seq_along(labels), labels = FALSE)
usr <- par("usr")
label_y <- usr[3] - 0.035 * (usr[4] - usr[3])
text(
  x = seq_along(labels), y = label_y, labels = labels,
  srt = 90, adj = 1, xpd = NA, cex = 0.62,
  col = ref_colour_for_group
)
text(
  x = seq_along(medians), y = medians,
  labels = median_labels, pos = 3, cex = 0.72, xpd = NA
)
mtext("Box fill: EVE category colour. X-axis label: reference category colour.",
      side = 1, line = 11.4, adj = 0, cex = 0.8)
mtext("Groups ordered by increasing median distance; numeric labels indicate medians.",
      side = 3, line = 0.4, adj = 0, cex = 0.8)
legend(
  "topleft",
  legend = unique(pairs$eve_category[match(EVE_TYPES, pairs$eve_class)]),
  fill = vapply(EVE_TYPES, function(et) {
    rr <- which(pairs$eve_class == et)[1L]
    if (is.na(rr)) FALLBACK_COLOUR else pairs$eve_colour[rr]
  }, character(1)),
  border = "grey25", bty = "n", cex = 0.75
)
par(old_par)
invisible(dev.off())

message("Wrote figures to: ", normalizePath(OUTDIR))
