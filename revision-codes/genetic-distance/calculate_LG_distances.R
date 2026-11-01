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

suppressPackageStartupMessages(library(phangorn))

FASTA <- "350andabove_nr_tedclado_linsi_mafftadd_trimmed_python.fasta"
METADATA_FILE <- "350andabove_nr_tedclado_metadata_reftaxa_final.tsv"
OUTDIR <- file.path("mds", "results_mds")

MODEL <- "LG"
GAMMA_K <- 4L
GAMMA_SHAPE <- NA_real_
GAMMA_SEED <- 1L
GAMMA_PER_GROUP <- 40L

EVE_TYPES <- c("TED_EVE", "CladoV_EVE")
REFERENCE_TYPES <- c("LTR-Retroelements-pol", "viral_ref")
EVE_LABELS <- c(
  TED_EVE = "TED-Homologous",
  CladoV_EVE = "CladoV_Homologous"
)

if (!file.exists(FASTA)) stop("Alignment file not found: ", FASTA)
if (!file.exists(METADATA_FILE)) stop("Metadata file not found: ", METADATA_FILE)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

read_fasta <- function(path) {
  x <- readLines(path, warn = FALSE)
  h <- which(startsWith(x, ">"))
  if (!length(h)) stop("No FASTA records found.")
  ends <- c(h[-1L] - 1L, length(x))
  headers <- substring(x[h], 2L)
  seqs <- vapply(seq_along(h), function(i) {
    first <- h[i] + 1L
    last <- ends[i]
    if (first > last) return("")
    lines <- x[first:last]
    paste0(gsub("\\s+", "", lines[!startsWith(lines, ">")]), collapse = "")
  }, character(1))
  if (anyDuplicated(headers)) stop("FASTA headers must be unique.")
  if (any(!nzchar(seqs))) stop("Empty sequence found in FASTA.")
  if (length(unique(nchar(seqs))) != 1L)
    stop("FASTA sequences must have equal aligned lengths.")
  list(headers = headers, sequences = toupper(seqs), length = nchar(seqs[1L]))
}

clean <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x %in% c("", "NA")] <- NA_character_
  x
}

fasta <- read_fasta(FASTA)
metadata <- read.delim(
  METADATA_FILE, sep = "\t", header = TRUE, check.names = FALSE,
  stringsAsFactors = FALSE
)
names(metadata) <- trimws(names(metadata))

required <- c("sequence_class", "Comparision_Meta_Category")
if (!all(required %in% names(metadata)))
  stop("Metadata is missing required columns: ",
       paste(setdiff(required, names(metadata)), collapse = ", "))

id_candidates <- intersect(
  c("sequence_id", "original_header", "processed_header", "sequence_name"),
  names(metadata)
)
if (!length(id_candidates)) stop("Metadata requires a sequence ID column.")

metadata_id <- NULL
metadata_index <- NULL
for (nm in id_candidates) {
  idx <- match(fasta$headers, trimws(as.character(metadata[[nm]])))
  if (all(!is.na(idx))) {
    metadata_id <- nm
    metadata_index <- idx
    break
  }
}
if (is.null(metadata_id)) stop("Every FASTA header must match a metadata ID.")
if (anyDuplicated(trimws(as.character(metadata[[metadata_id]]))))
  stop("Metadata IDs must be unique.")
metadata <- metadata[metadata_index, , drop = FALSE]

sequence_class <- clean(metadata$sequence_class)
category <- clean(metadata$Comparision_Meta_Category)
category[is.na(category)] <- "Unlabelled"
category[category == "Unassigned"] <- "Unassigned viruses"

class_map <- c(
  TED_EVE = "EVE", CladoV_EVE = "EVE",
  `LTR-Retroelements-pol` = "LTR retroelement",
  viral_ref = "Extant virus"
)
unknown <- setdiff(unique(sequence_class), names(class_map))
if (length(unknown))
  stop("Unrecognised sequence_class value(s): ", paste(unknown, collapse = ", "))

for (g in names(EVE_LABELS)) {
  bad <- sequence_class == g & category != EVE_LABELS[[g]]
  if (any(bad))
    stop("EVE sequence_class and Comparision_Meta_Category disagree for ",
         sum(bad), " record(s) in ", g, ".")
}

alignment <- phangorn::read.phyDat(file = FASTA, format = "fasta", type = "AA")
if (length(alignment) != length(fasta$headers))
  stop("Alignment parsing changed the sequence count.")
if (!identical(names(alignment), fasta$headers))
  stop("Parsed alignment labels differ from FASTA headers.")

if (is.na(GAMMA_SHAPE)) {
  set.seed(GAMMA_SEED)
  strata <- interaction(sequence_class, category, drop = TRUE, lex.order = TRUE)
  sample_idx <- unlist(lapply(split(seq_along(strata), strata), function(idx) {
    if (length(idx) > GAMMA_PER_GROUP) sample(idx, GAMMA_PER_GROUP) else idx
  }), use.names = FALSE)
  sample_idx <- sort(unique(sample_idx))
  if (length(sample_idx) < 4L)
    stop("Too few sequences to estimate gamma shape.")

  sample_alignment <- subset(alignment, subset = sample_idx)
  sample_dist <- as.matrix(phangorn::dist.ml(
    sample_alignment, model = MODEL, exclude = "pairwise"
  ))

  keep <- rowSums(!is.finite(sample_dist)) == 0L
  while (sum(keep) >= 4L && any(!keep)) {
    sample_dist <- sample_dist[keep, keep, drop = FALSE]
    sample_alignment <- subset(
      sample_alignment,
      subset = match(rownames(sample_dist), names(sample_alignment))
    )
    keep <- rowSums(!is.finite(sample_dist)) == 0L
  }
  if (sum(keep) < 4L)
    stop("Gamma-shape estimation failed; set GAMMA_SHAPE explicitly.")
  if (any(!keep)) {
    sample_dist <- sample_dist[keep, keep, drop = FALSE]
    sample_alignment <- subset(
      sample_alignment,
      subset = match(rownames(sample_dist), names(sample_alignment))
    )
  }

  tree <- phangorn::NJ(stats::as.dist(sample_dist))
  tree$edge.length[!is.finite(tree$edge.length) | tree$edge.length <= 0] <- 1e-6
  fit <- phangorn::pml(tree, sample_alignment, model = MODEL, k = GAMMA_K, shape = 1)
  fit <- phangorn::optim.pml(
    fit, optEdge = TRUE, optGamma = TRUE,
    control = phangorn::pml.control(trace = 0)
  )
  GAMMA_SHAPE <- fit$shape
}
if (!is.finite(GAMMA_SHAPE) || GAMMA_SHAPE <= 0)
  stop("GAMMA_SHAPE must be a positive finite number.")

message("Calculating pairwise LG+G distances...")
D <- as.matrix(phangorn::dist.ml(
  alignment, model = MODEL, exclude = "pairwise",
  k = GAMMA_K, shape = GAMMA_SHAPE
))
if (is.null(rownames(D))) rownames(D) <- fasta$headers
if (is.null(colnames(D))) colnames(D) <- fasta$headers
D <- D[fasta$headers, fasta$headers, drop = FALSE]

if (any(D < 0, na.rm = TRUE)) stop("Negative distances were returned.")
if (max(abs(D - t(D)), na.rm = TRUE) > 1e-8)
  stop("Distance matrix is not symmetric.")
if (max(abs(diag(D)), na.rm = TRUE) > 1e-8)
  stop("Distance matrix diagonal is not zero.")

matrix_file <- file.path(OUTDIR, "LG_pairwise_distance_matrix.tsv")
write.table(
  D, matrix_file, sep = "\t", quote = FALSE,
  row.names = TRUE, col.names = NA, na = "NA"
)

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
groups <- unique(data.frame(
  lineage = ref_category,
  reference_type = ref_type,
  stringsAsFactors = FALSE
))

rows <- list()
z <- 0L
for (i in seq_len(nrow(groups))) {
  refs <- ref_idx[
    ref_category == groups$lineage[i] &
      ref_type == groups$reference_type[i]
  ]
  for (eve_type in EVE_TYPES) {
    eves <- eve_idx[sequence_class[eve_idx] == eve_type]
    values <- as.vector(D[eves, refs, drop = FALSE])
    values <- values[is.finite(values)]
    if (!length(values)) next
    z <- z + 1L
    rows[[z]] <- data.frame(
      eve_class = eve_type,
      eve_label = unname(EVE_LABELS[eve_type]),
      reference_lineage = groups$lineage[i],
      reference_type = groups$reference_type[i],
      n_eve_sequences = length(eves),
      n_reference_sequences = length(refs),
      n_finite_distances = length(values),
      median_distance = median(values),
      q1_distance = unname(quantile(values, 0.25)),
      q3_distance = unname(quantile(values, 0.75)),
      mean_distance = mean(values),
      stringsAsFactors = FALSE
    )
  }
}
if (!length(rows)) stop("No finite EVE-to-reference distances were calculated.")
eve_ref_summary <- do.call(rbind, rows)
eve_ref_summary <- eve_ref_summary[
  order(eve_ref_summary$median_distance, eve_ref_summary$reference_lineage,
        eve_ref_summary$eve_class),
]

all_pair_values <- D[upper.tri(D)]
finite_values <- all_pair_values[is.finite(all_pair_values)]

summary_lines <- c(
  "Pairwise LG+G evolutionary distances",
  "",
  paste("Alignment:", basename(FASTA)),
  paste("Metadata:", basename(METADATA_FILE)),
  paste("Sequences:", length(fasta$headers)),
  paste("Alignment length (amino acids):", fasta$length),
  paste("Substitution model:", MODEL),
  paste0("Rate heterogeneity: discrete gamma; categories = ", GAMMA_K,
         "; shape = ", signif(GAMMA_SHAPE, 6)),
  "Gap and ambiguous-residue treatment: pairwise deletion (phangorn)",
  paste("Unique sequence pairs:", length(all_pair_values)),
  paste("Finite pairwise distances:", length(finite_values)),
  paste("Non-finite pairwise distances:", sum(!is.finite(all_pair_values))),
  if (length(finite_values))
    paste("Median of all finite pairwise distances:", signif(median(finite_values), 6))
  else "Median of all finite pairwise distances: unavailable",
  "",
  "EVE-to-reference distances by EVE class and reference group",
  "Distances are amino-acid substitutions per site under LG+G.",
  "Median, first quartile (Q1), and third quartile (Q3) are reported.",
  ""
)

format_num <- function(x) format(signif(x, 6), scientific = FALSE, trim = TRUE)
for (i in seq_len(nrow(eve_ref_summary))) {
  r <- eve_ref_summary[i, ]
  summary_lines <- c(
    summary_lines,
    paste0(
      r$eve_label, " | ", r$reference_lineage, " | ", r$reference_type,
      " | n=", r$n_finite_distances,
      " | median=", format_num(r$median_distance),
      " | Q1=", format_num(r$q1_distance),
      " | Q3=", format_num(r$q3_distance),
      " | mean=", format_num(r$mean_distance)
    )
  )
}
summary_lines <- c(
  summary_lines, "",
  "Output: LG_pairwise_distance_matrix.tsv",
  "Output: distance_summary.txt"
)
writeLines(summary_lines, file.path(OUTDIR, "distance_summary.txt"))

message("Wrote distance matrix and summary to: ", normalizePath(OUTDIR))
