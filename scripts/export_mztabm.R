#!/usr/bin/env Rscript

# Export one mzTab-M 2.1 metadata (profile M) file per [[datasets]] entry in
# sample_metadata.toml, using https://github.com/kozo2/RmzTabM.
#
# The TOML file records the sample-annotation block of each RIKEN LIPIDOMICS
# workbook, not the quantification matrix, so the files contain the MTD
# section only.

args <- commandArgs(trailingOnly = TRUE)
toml_path <- if (length(args) >= 1) args[[1]] else "sample_metadata.toml"
out_dir <- if (length(args) >= 2) args[[2]] else "mztab"

if (!requireNamespace("RmzTabM", quietly = TRUE)) {
    stop("RmzTabM is required. Install it from https://github.com/kozo2/RmzTabM")
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite is required to read the TOML file via Python.")
}

library(RmzTabM)

read_sample_metadata <- function(path) {
    py <- tempfile(fileext = ".py")
    writeLines(c(
        "import json, sys, tomllib",
        "with open(sys.argv[1], 'rb') as handle:",
        "    json.dump(tomllib.load(handle), sys.stdout)"
    ), py)
    json <- system2("python3", c(py, path), stdout = TRUE)
    if (!is.null(attr(json, "status")) && attr(json, "status") != 0) {
        stop("Failed to parse ", path)
    }
    jsonlite::fromJSON(paste(json, collapse = "\n"), simplifyVector = FALSE)
}

field <- function(sample, key) {
    value <- sample[[key]]
    if (is.null(value) || length(value) == 0 || is.na(value)) {
        return(NA_character_)
    }
    as.character(value)
}

or_null <- function(x) {
    x <- as.character(x)
    # The workbooks use "NA" for not applicable. mzTab uses "null", and a bare
    # NA token is read back as a missing value.
    x[is.na(x) | !nzchar(x) | x %in% c("NA", "N/A", "null")] <- "null"
    x
}

user_param <- function(name) {
    name <- gsub("[\r\n\t]", " ", name)
    name <- gsub(",", ";", name, fixed = TRUE)
    name <- gsub("[", "(", name, fixed = TRUE)
    name <- gsub("]", ")", name, fixed = TRUE)
    sprintf("[,, %s, ]", name)
}

is_binomial <- function(x) {
    grepl("^[A-Z][a-z]+ [a-z]", x)
}

split_taxon_tissue <- function(tissue) {
    parts <- trimws(strsplit(tissue, "/", fixed = TRUE)[[1]])
    if (length(parts) == 1) {
        if (is_binomial(parts[[1]])) {
            return(list(species = parts[[1]], organ = NA_character_))
        }
        return(list(species = NA_character_, organ = parts[[1]]))
    }
    if (is_binomial(parts[[1]]) && !is_binomial(parts[[2]])) {
        return(list(
            species = parts[[1]],
            organ = paste(parts[-1], collapse = "/")
        ))
    }
    if (length(parts) >= 2 && is_binomial(parts[[2]])) {
        return(list(species = parts[[2]], organ = parts[[1]]))
    }
    list(species = NA_character_, organ = tissue)
}

MOUSE <- "[NCBITaxon, NCBITaxon:10090, Mus musculus, ]"
HUMAN <- "[NCBITaxon, NCBITaxon:9606, Homo sapiens, ]"

organism_fields <- function(category, tissue) {
    species <- NULL
    organ <- NULL
    cell_type <- NULL
    if (is.na(category) || category %in% c("Blank", "QC")) {
        return(list(species = species, tissue = organ, cell_type = cell_type))
    }
    if (category %in% c("Mouse", "Mouse cultured cell")) {
        species <- MOUSE
    } else if (category %in% c("Human", "Human cultured cell")) {
        species <- HUMAN
    }
    if (grepl("cultured cell", category, fixed = TRUE)) {
        if (!is.na(tissue) && nzchar(tissue)) {
            cell_type <- user_param(tissue)
        }
        return(list(species = species, tissue = organ, cell_type = cell_type))
    }
    if (is.na(tissue) || !nzchar(tissue) || tissue %in% c("Blank", "QC")) {
        return(list(species = species, tissue = organ, cell_type = cell_type))
    }
    if (category %in% c("Plant", "Algae")) {
        parsed <- split_taxon_tissue(tissue)
        if (!is.na(parsed$species)) {
            species <- user_param(parsed$species)
        }
        if (!is.na(parsed$organ)) {
            organ <- user_param(parsed$organ)
        }
        return(list(species = species, tissue = organ, cell_type = cell_type))
    }
    list(species = species, tissue = user_param(tissue), cell_type = cell_type)
}

infer_polarity <- function(name) {
    low <- tolower(name)
    positive <- grepl("(^|[_. -])pos($|[_. -])|positive|posi_|lipp", low) |
        grepl("pp$", low)
    negative <- grepl("(^|[_. -])neg($|[_. -])|nega|negative|lipn", low) |
        grepl("nn$", low)
    polarity <- rep(NA_character_, length(name))
    polarity[negative & !positive] <- "negative"
    polarity[positive & !negative] <- "positive"
    polarity
}

quant_unit_param <- function(units) {
    units <- unique(units[units != "null"])
    if (length(units) == 1) {
        return(user_param(units))
    }
    "[PRIDE, PRIDE:0000330, Arbitrary quantification unit, ]"
}

# mtdSample() / mtdAssay() sort index columns lexicographically once a file
# has 10 or more entries. Restore numeric index order within each MTD block
# while keeping the section order produced by mtdSort().
reorder_mtd_indices <- function(mtd) {
    top <- sub("\\[.*", "", mtd[, 1L])
    runs <- rle(top)
    ends <- cumsum(runs$lengths)
    starts <- ends - runs$lengths + 1L
    rows <- integer()
    for (i in seq_along(starts)) {
        idx <- starts[[i]]:ends[[i]]
        keys <- vapply(mtd[idx, 1L], function(field) {
            nums <- regmatches(field, gregexpr("[0-9]+", field))[[1]]
            if (!length(nums)) {
                return("")
            }
            paste(sprintf("%06d", as.integer(nums)), collapse = ".")
        }, character(1))
        rows <- c(rows, idx[order(keys, seq_along(keys))])
    }
    mtd[rows, , drop = FALSE]
}

opt_chr <- function(values) {
    lapply(values, function(value) {
        if (is.null(value) || length(value) == 0 || is.na(value) || !nzchar(value)) {
            NULL
        } else {
            value
        }
    })
}

dataset_to_mztab <- function(dataset, path) {
    samples <- dataset$samples
    if (is.null(samples) || !length(samples)) {
        stop("Dataset ", dataset$name, " has no samples")
    }
    rows <- lapply(samples, function(sample) {
        name <- field(sample, "name")
        category <- field(sample, "category")
        if (is.na(name) || !nzchar(name) || is.na(category)) {
            return(NULL)
        }
        organism <- organism_fields(category, field(sample, "tissue_species"))
        data.frame(
            name = name,
            column = field(sample, "column"),
            public_private = or_null(field(sample, "public_private")),
            category = or_null(category),
            tissue_species = or_null(field(sample, "tissue_species")),
            genotype_background = or_null(field(sample, "genotype_background")),
            perturbation = or_null(field(sample, "perturbation")),
            diet_culture = or_null(field(sample, "diet_culture")),
            biological_replicate = or_null(field(sample, "biological_replicate")),
            technical_replicate = or_null(field(sample, "technical_replicate")),
            unit = or_null(field(sample, "unit")),
            species = if (is.null(organism$species)) NA_character_ else organism$species,
            tissue = if (is.null(organism$tissue)) NA_character_ else organism$tissue,
            cell_type = if (is.null(organism$cell_type)) NA_character_ else organism$cell_type,
            stringsAsFactors = FALSE
        )
    })
    df <- do.call(rbind, rows)
    if (is.null(df) || !nrow(df)) {
        stop("Dataset ", dataset$name, " has no annotated samples")
    }
    rownames(df) <- NULL

    df$polarity_from_name <- infer_polarity(df$name)
    known <- unique(df$polarity_from_name[!is.na(df$polarity_from_name)])
    df$scan_polarity <- df$polarity_from_name
    if (length(known) == 1) {
        df$scan_polarity[is.na(df$scan_polarity)] <- known
    }
    unspecified <- is.na(df$scan_polarity)
    df$scan_polarity[unspecified] <- "positive"

    ann <- paste(df$species, df$tissue, df$cell_type, sep = "|")
    conflict <- tapply(ann, df$name, function(values) length(unique(values)) > 1)
    conflict_names <- names(conflict)[as.logical(conflict)]
    df$sample_id <- df$name
    if (length(conflict_names)) {
        hit <- df$name %in% conflict_names
        df$sample_id[hit] <- paste(df$name[hit], df$column[hit], sep = " ")
    }
    df$assay <- paste0(df$name, " (", df$column, ")")
    if (anyDuplicated(df$assay)) {
        stop("Duplicated assay ids in ", dataset$name)
    }

    sample_order <- unique(df$sample_id)
    sample_index <- match(df$sample_id, sample_order)
    first <- !duplicated(df$sample_id)

    mtd <- mtdSkeleton(
        id = gsub("[^A-Za-z0-9._+-]", "_", sub("\\.xlsx$", "", dataset$name)),
        software = sprintf("[,, RmzTabM, %s]", as.character(packageVersion("RmzTabM"))),
        small_molecule_quantification_unit = quant_unit_param(df$unit),
        small_molecule_feature_quantification_unit = quant_unit_param(df$unit),
        mztab_profile = "M"
    )
    mtd <- setMtdField(mtd, field = "title", value = dataset$sheet)
    mtd <- setMtdField(
        mtd,
        field = "description",
        value = paste(
            "Sample metadata from RIKEN LIPIDOMICS workbook",
            dataset$name,
            "exported from sample_metadata.toml.",
            "Quantification values are not included."
        )
    )
    mtd <- setMtdField(
        mtd,
        field = "uri",
        value = "https://metabography.riken.jp/menta.cgi/lipidomics/download_data_set"
    )
    if (any(grepl("NCBITaxon", df$species))) {
        mtd <- rbind(mtd, rbind(
            c("cv[4]-label", "NCBITaxon"),
            c("cv[4]-full_name", "NCBI organismal classification"),
            c("cv[4]-version", "2026-07-12"),
            c("cv[4]-uri", "http://purl.obolibrary.org/obo/ncbitaxon.owl")
        ))
    }

    mtd_sample <- mtdSample(
        sample = sample_order,
        species = opt_chr(df$species[first][match(sample_order, df$sample_id[first])]),
        tissue = opt_chr(df$tissue[first][match(sample_order, df$sample_id[first])]),
        cell_type = opt_chr(df$cell_type[first][match(sample_order, df$sample_id[first])])
    )
    run_parameter <- lapply(seq_len(nrow(df)), function(i) {
        if (!unspecified[[i]]) {
            return(NULL)
        }
        "[,, scan polarity present in sample name, no]"
    })
    mtd_run <- mtdMsRun(
        location = rep("null", nrow(df)),
        scan_polarity = df$scan_polarity,
        parameters = run_parameter
    )
    mtd_assay <- mtdAssay(
        assay = df$assay,
        sample_ref = paste0("sample[", sample_index, "]"),
        ms_run_ref = paste0("ms_run[", seq_len(nrow(df)), "]")
    )
    groups <- c(
        "public_private", "category", "genotype_background", "perturbation",
        "diet_culture", "biological_replicate", "technical_replicate", "unit"
    )
    group_description <- c(
        "Public/Private", "Category", "Genotype/Background", "Perturbation",
        "Diet/Culture", "Biological replicate", "Technical replicate", "Unit"
    )
    keep <- vapply(groups, function(group) !all(df[[group]] == "null"), logical(1))
    mtd_svar <- mtdStudyVariables(
        df,
        groups = groups[keep],
        group_description = group_description[keep]
    )
    mtd <- reorder_mtd_indices(
        mtdSort(rbind(mtd, mtd_sample, mtd_run, mtd_assay, mtd_svar))
    )

    comments <- c(
        "Generated with RmzTabM from sample_metadata.toml.",
        "Source: RIKEN LIPIDOMICS, CC BY-NC 4.0."
    )
    if (any(unspecified)) {
        comments <- c(
            comments,
            paste(
                "Scan polarity was not encoded in the sample name for",
                sum(unspecified),
                "assay(s); those runs are reported as positive scan and",
                "carry ms_run-parameter [,, scan polarity present in sample name, no]."
            )
        )
    }
    writeMzTabM(list(MTD = mtd), path = path, comments = comments)
    invisible(df)
}

metadata <- read_sample_metadata(toml_path)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

written <- character()
for (dataset in metadata$datasets) {
    outfile <- file.path(out_dir, sub("\\.xlsx$", ".mzTab", dataset$name))
    dataset_to_mztab(dataset, outfile)
    written <- c(written, outfile)
}

message("Wrote ", length(written), " mzTab-M files to ", out_dir)
