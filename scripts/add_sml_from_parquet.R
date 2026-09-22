#!/usr/bin/env Rscript

# Add the Small Molecule (SML) table to each mzTab-M file in mztab/ from the
# matching Parquet file in parquet/, using RmzTabM.
#
# Assay order in the existing MTD section is preserved. Parquet columns that
# repeat a sample name are stored with a trailing "_1" (and so on); those
# columns are matched to the later assays that share the same sample name.
# Columns that are not assays (Lipid IS, standard class, and MS-DIAL
# annotation fields) are not written as abundance columns.

args <- commandArgs(trailingOnly = TRUE)
mztab_dir <- if (length(args) >= 1) args[[1]] else "mztab"
parquet_dir <- if (length(args) >= 2) args[[2]] else "parquet"

if (!requireNamespace("RmzTabM", quietly = TRUE)) {
    stop("RmzTabM is required. Install it from https://github.com/kozo2/RmzTabM")
}
if (!requireNamespace("nanoparquet", quietly = TRUE)) {
    stop("nanoparquet is required to read the Parquet files.")
}

library(RmzTabM)

# Parquet columns that are not quantification assays.
metadata_columns <- c(
    "Alignment ID", "Average Rt(min)", "Average Mz", "Metabolite name",
    "Adduct type", "Post curation result", "Fill %", "MS/MS assigned",
    "Reference RT", "Reference m/z", "Formula", "Ontology", "INCHIKEY",
    "SMILES", "Annotation tag (VS1.0)", "RT matched", "m/z matched",
    "MS/MS matched", "Comment", "Manually modified",
    "Manually modified for annotation", "Isotope tracking parent ID",
    "Isotope tracking weight number", "Total score", "RT similarity",
    "Dot product", "Reverse dot product", "Fragment presence %",
    "S/N average", "Spectrum reference file name", "MS1 isotopic spectrum",
    "MS/MS spectrum", "Average mobility", "Average CCS", "Reference CCS",
    "CCS similarity", "CCS matched", "Method", "Lipid IS", "standard class"
)

# Annotation columns copied into optional SML fields.
optional_columns <- c(
    "Alignment ID" = "alignment_id",
    "Average Rt(min)" = "average_rt_min",
    "Average Mz" = "average_mz",
    "Ontology" = "ontology",
    "Annotation tag (VS1.0)" = "annotation_tag",
    "Total score" = "total_score",
    "Fill %" = "fill_percent",
    "MS/MS assigned" = "msn_assigned",
    "Reference RT" = "reference_rt",
    "Reference m/z" = "reference_mz",
    "Comment" = "comment",
    "Post curation result" = "post_curation_result",
    "S/N average" = "sn_average",
    "Spectrum reference file name" = "spectrum_reference_file_name",
    "RT matched" = "rt_matched",
    "m/z matched" = "mz_matched",
    "MS/MS matched" = "msms_matched",
    "Dot product" = "dot_product",
    "Reverse dot product" = "reverse_dot_product",
    "Fragment presence %" = "fragment_presence_percent",
    "RT similarity" = "rt_similarity",
    "Isotope tracking parent ID" = "isotope_tracking_parent_id",
    "Isotope tracking weight number" = "isotope_tracking_weight_number",
    "Manually modified" = "manually_modified",
    "Manually modified for annotation" = "manually_modified_for_annotation",
    "Average CCS" = "average_ccs",
    "Average mobility" = "average_mobility",
    "CCS similarity" = "ccs_similarity",
    "Reference CCS" = "reference_ccs",
    "CCS matched" = "ccs_matched",
    "Method" = "method"
)

as_null_chr <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x <- gsub("[\r\n\t]", " ", x)
    x <- gsub("|", "/", x, fixed = TRUE)
    x[x %in% c("", "NA", "NaN", "null", "#N/A", "N/A")] <- NA_character_
    x
}

assay_base_name <- function(assay) {
    sub(" \\([A-Z]+\\)$", "", assay)
}

match_abundance_columns <- function(column_names, assay_names) {
    bases <- assay_base_name(assay_names)
    used <- rep(FALSE, length(column_names))
    matched <- character(length(bases))
    for (i in seq_along(bases)) {
        base <- bases[[i]]
        suffixes <- c(base, paste0(base, "_", seq_len(5L)))
        candidates <- which(column_names %in% suffixes & !used)
        candidates <- candidates[!column_names[candidates] %in% metadata_columns]
        if (!length(candidates)) {
            stop("No Parquet column for assay ", assay_names[[i]])
        }
        pick <- candidates[[1L]]
        used[[pick]] <- TRUE
        matched[[i]] <- column_names[[pick]]
    }
    matched
}

read_comments <- function(path) {
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    comments <- lines[startsWith(lines, "COM\t")]
    if (!length(comments)) {
        return(character())
    }
    sub("^COM\t", "", comments)
}

add_sml <- function(mztab_path, parquet_path) {
    pq <- nanoparquet::read_parquet(parquet_path)
    mz <- readMzTabM(mztab_path)
    meta <- mtd(mz)
    assays <- unname(getMtdField(meta, "^assay\\[[0-9]+\\]$"))
    if (!length(assays)) {
        stop("No assays in ", mztab_path)
    }
    abundance_columns <- match_abundance_columns(colnames(pq), assays)
    abundances <- vapply(
        pq[abundance_columns],
        function(column) {
            values <- as.character(column)
            # Excel exports use "#N/A" for missing quantification cells.
            values[values %in% c("", "NA", "NaN", "null", "#N/A", "N/A")] <- NA_character_
            as.numeric(values)
        },
        numeric(nrow(pq))
    )
    if (is.null(dim(abundances))) {
        abundances <- matrix(abundances, nrow = nrow(pq))
    }
    colnames(abundances) <- abundance_columns

    inchikey <- as_null_chr(pq[["INCHIKEY"]])
    database_identifier <- ifelse(
        is.na(inchikey),
        NA_character_,
        paste0("INCHIKEY:", inchikey)
    )
    present <- intersect(names(optional_columns), colnames(pq))
    optional <- lapply(present, function(column) as_null_chr(pq[[column]]))
    names(optional) <- unname(optional_columns[present])

    sml <- do.call(smlCreate, c(
        list(
            x = abundances,
            database_identifier = database_identifier,
            chemical_formula = as_null_chr(pq[["Formula"]]),
            smiles = as_null_chr(pq[["SMILES"]]),
            chemical_name = as_null_chr(pq[["Metabolite name"]]),
            adduct_ions = as_null_chr(pq[["Adduct type"]])
        ),
        optional
    ))
    sml <- smlAddStudyVariableColumns(sml, meta)
    numeric_columns <- vapply(sml, is.numeric, logical(1))
    sml[numeric_columns] <- lapply(sml[numeric_columns], function(column) {
        column[!is.finite(column)] <- NA_real_
        column
    })

    meta <- setMtdField(meta, field = "mzTab-profile", value = "M+S")
    description <- unname(getMtdField(meta, "^description$"))
    description <- sub(
        "Quantification values are not included\\.",
        "Quantification values are reported in the SML section.",
        description
    )
    if (!any(grepl("SML section", description))) {
        description <- paste(
            description,
            "Quantification values are reported in the SML section."
        )
    }
    meta <- setMtdField(meta, field = "description", value = description)

    comments <- read_comments(mztab_path)
    sml_comment <- paste(
        "SML abundances and annotations were added from",
        basename(parquet_path),
        "with RmzTabM."
    )
    if (!any(comments == sml_comment)) {
        comments <- c(comments, sml_comment)
    }
    writeMzTabM(list(MTD = meta, SML = sml), path = mztab_path, comments = comments)
    invisible(nrow(sml))
}

mztab_files <- list.files(mztab_dir, pattern = "\\.mzTab$", full.names = TRUE)
if (!length(mztab_files)) {
    stop("No mzTab-M files in ", mztab_dir)
}

for (mztab_path in mztab_files) {
    parquet_path <- file.path(
        parquet_dir,
        sub("\\.mzTab$", ".parquet", basename(mztab_path))
    )
    if (!file.exists(parquet_path)) {
        stop("Missing Parquet file for ", basename(mztab_path))
    }
    molecules <- add_sml(mztab_path, parquet_path)
    message(basename(mztab_path), ": ", molecules, " small molecules")
}
