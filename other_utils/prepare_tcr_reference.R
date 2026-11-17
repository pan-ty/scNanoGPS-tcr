#!/usr/bin/env Rscript
# Convert local OCM reference/sample metadata to portable exporter inputs.
# Source files are read-only. Do not commit the generated study metadata.

args <- commandArgs(trailingOnly = TRUE)
usage <- paste(
  "Usage: Rscript other_utils/prepare_tcr_reference.R",
  "--sample-paths samples.xlsx --pool-id POOL --output NEW_DIRECTORY",
  "[--ob-rda OB.rda]"
)
if (length(args) == 0L || "--help" %in% args) {
  cat(usage, "\n")
  quit(status = 0)
}
if (length(args) %% 2L) stop(usage)
keys <- args[seq.int(1L, length(args), 2L)]
if (anyDuplicated(keys)) stop("Repeated argument")
if (any(!keys %in% c("--sample-paths", "--pool-id", "--output", "--ob-rda"))) stop(usage)
opt <- setNames(as.list(args[seq.int(2L, length(args), 2L)]), keys)
required <- c("--sample-paths", "--pool-id", "--output")
if (!all(required %in% names(opt))) stop(usage)
if (!requireNamespace("readxl", quietly = TRUE)) stop("The readxl R package is required")

sheet <- readxl::read_xlsx(opt[["--sample-paths"]])
if (!all(c("sample_id", "NP_id", "NB_id") %in% names(sheet))) {
  stop("Spreadsheet requires sample_id, NP_id, NB_id columns")
}
selected <- unique(as.data.frame(sheet[!is.na(sheet$NP_id) & sheet$NP_id == opt[["--pool-id"]],
                                      c("sample_id", "NP_id", "NB_id")]))
if (!nrow(selected)) stop("Pool ID is absent from the spreadsheet")
if (any(is.na(selected$sample_id)) || any(!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", selected$sample_id))) {
  stop("Unsafe or empty sample ID")
}
has_ob <- !is.na(selected$NB_id) & nzchar(selected$NB_id)
if (any(has_ob) && !all(has_ob)) stop("Pool mixes OCM and non-OCM rows")
ocm <- all(has_ob)
if (ocm && (any(!selected$NB_id %in% paste0("OB", 1:4)) || anyDuplicated(selected$NB_id))) {
  stop("OCM pool must map each OB group to exactly one sample")
}
if (!ocm && nrow(selected) != 1L) stop("A non-OCM library requires exactly one sample")

ob <- NULL
if (ocm) {
  if (!"--ob-rda" %in% names(opt)) stop("OCM mode requires --ob-rda")
  environment <- new.env(parent = emptyenv())
  loaded <- load(opt[["--ob-rda"]], envir = environment)
  if (!"OB" %in% loaded) stop("RDA must contain an OB object")
  ob <- environment$OB
  if (!all(paste0("OB", 1:4) %in% names(ob))) stop("OB must contain OB1, OB2, OB3, OB4")
  all_barcodes <- unlist(ob[paste0("OB", 1:4)], use.names = FALSE)
  if (!is.character(all_barcodes) || anyNA(all_barcodes) ||
      any(!grepl("^[ACGT]{16}$", all_barcodes)) || anyDuplicated(all_barcodes)) {
    stop("OB lists must contain non-overlapping, unique, full 16-base DNA barcodes")
  }
  rm(all_barcodes)
}

output <- path.expand(opt[["--output"]])
if (grepl("onedrive", output, ignore.case = TRUE)) stop("Choose a local output outside OneDrive")
if (file.exists(output)) stop("Output directory already exists")
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
staging <- tempfile(pattern = ".tcr_reference.", tmpdir = dirname(output))
dir.create(staging)
tryCatch({
  write.table(selected, file.path(staging, "source_samples.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE, na = "")
  if (ocm) {
    map <- data.frame(ocm_group = selected$NB_id, sample_id = selected$sample_id)
    write.table(map, file.path(staging, "sample_map.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    connection <- gzfile(file.path(staging, "ocm_barcodes.tsv.gz"), open = "wt")
    tryCatch({
      writeLines("barcode\tocm_group", connection)
      for (group in paste0("OB", 1:4)) {
        writeLines(paste(ob[[group]], group, sep = "\t"), connection)
      }
    }, finally = close(connection))
    writeLines("ocm", file.path(staging, "mode.txt"))
  } else {
    writeLines("non_ocm", file.path(staging, "mode.txt"))
    writeLines(selected$sample_id, file.path(staging, "sample_id.txt"))
  }
  if (file.exists(output) || !file.rename(staging, output)) stop("Could not publish new output directory")
}, finally = {
  if (dir.exists(staging)) unlink(staging, recursive = TRUE)
})
cat("Mode:", if (ocm) "ocm" else "non_ocm", "\n")
print(selected, row.names = FALSE)
cat("Reference outputs:", normalizePath(output), "\n")
if (!ocm) cat("Supply a chemistry-matched --whitelist and use --sample-id with tcr_export.py.\n")
