# =============================================================================
# ensannot.R  -  Ensembl gene annotation framework (no biomaRt needed)
#
# Package layout: save as  R/ensannot.R  in a package skeleton, then run
#   devtools::document()      # writes NAMESPACE and man/*.Rd from the roxygen blocks
#
# DESCRIPTION needs:
#   Depends:  R (>= 4.0.0)                 # tools::R_user_dir()
#   Imports:  data.table, utils, tools
#   Suggests: R.utils                      # much faster .gz reading
#   Roxygen:  list(markdown = TRUE)
#
# Standalone use (no package): run library(data.table), then source("ensannot.R").
# =============================================================================

#' ensannot: Ensembl gene annotation without biomaRt
#'
#' @description
#' Builds a local, reusable gene annotation table straight from the Ensembl FTP
#' site (the core MySQL dump files, read as plain tab-delimited text, so no SQL
#' server or SQL knowledge is needed) and uses it to annotate your own tables.
#'
#' @section Workflow:
#' **Part 1 - build once, reuse forever**
#' \itemize{
#'   \item [ens_build_annotation()] downloads, parses and caches the table for one
#'     species / release / assembly.
#'   \item [ens_load_annotation()] reloads a cached object (or a saved file).
#'   \item [ens_list_cached()] shows what is cached.
#'   \item [ens_resolve_db_url()] turns species, assembly and release into the FTP folder.
#' }
#' **Part 2 - annotate your data**
#' \itemize{
#'   \item [ens_annotate_ids()]: details for Ensembl gene IDs.
#'   \item [ens_annotate_regions()]: all genes overlapping each region (returns a
#'     per-region table and a gene table).
#'   \item [ens_annotate_loci()]: genes overlapping each locus, or the closest
#'     genes within a distance.
#' }
#' [ens_example_annotation()] returns a tiny synthetic annotation so every
#' example in this package runs offline.
#'
#' @section Column names:
#' Output columns follow biomaRt attribute names (`ensembl_gene_id`,
#' `chromosome_name`, `start_position`, `end_position`, `strand`,
#' `external_gene_name`, `description`, `external_synonym`, `gene_biotype`,
#' `entrezgene_id`, `ensembl_transcript_id`) plus one extra column,
#' `custom_gene_name`. See [ens_annot] for details.
#'
#' @section Options:
#' `options(ensannot.cache_dir = "path")` sets where built objects are cached
#' (default: [ens_cache_dir()]).
#'
#' @keywords internal
#' @import data.table
#' @importFrom utils download.file head tail
#' @importFrom tools R_user_dir
"_PACKAGE"

# Column names used with data.table's non-standard evaluation (silences R CMD check).
utils::globalVariables(c(
  ".", ".N", ".SD", ".chr_num", "attrib_type_id", "biotype", "chr", "chromosome_name",
  "code", "coord_system_id", "custom_gene_name", "db_name", "dbprimary_acc",
  "description", "display_xref_id", "dist", "end_position", "ensembl_gene_id",
  "ensembl_id", "ensembl_object_type", "ensembl_transcript_id", "entrez",
  "entrezgene_id", "external_db_id", "external_gene_name", "external_synonym",
  "ge", "gene_id", "gi", "gs", "internal_id", "is_current", "mind", "name",
  "qe", "qi", "qs", "rel", "seq_region_id", "seq_region_start", "seq_region_end",
  "seq_region_strand", "start_position", "stable_id", "synonym", "transcript_id",
  "translation_id", "xref_id"
))

.ENS_FTP        <- "https://ftp.ensembl.org/pub"
.ENS_ATTR_COLS  <- c("ensembl_gene_id", "external_gene_name", "custom_gene_name",
                     "description", "external_synonym", "gene_biotype",
                     "entrezgene_id", "ensembl_transcript_id")
.ENS_COORD_COLS <- c("chromosome_name", "start_position", "end_position", "strand")

# ------------------------------------------------------------------ utilities
# (internal helpers: documented with plain comments, not exported)

.msg <- function(verbose, ...) {
  if (isTRUE(verbose)) message(format(Sys.time(), "%H:%M:%S "), ...)
}

# unique, non-missing values joined with sep; NA if nothing is left
.collapse_unique <- function(x, sep = "; ") {
  x <- as.character(x)
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (!length(x)) NA_character_ else paste(x, collapse = sep)
}

# "chr1" -> "1", "chrM"/"m" -> "MT", "x" -> "X"
.norm_chr <- function(x) {
  x <- sub("^chr", "", trimws(as.character(x)), ignore.case = TRUE)
  special <- tolower(x) %in% c("x", "y", "mt", "m")
  x[special] <- toupper(x[special])
  x[x == "M"] <- "MT"
  x
}

# lower-case, spaces/dots -> "_"   ("Capra hircus" -> "capra_hircus")
.norm_species <- function(x) tolower(gsub("[ .]+", "_", trimws(x)))

# comparison key that ignores case and punctuation ("Black_Bengal" == "blackbengal")
.key <- function(x) gsub("[^a-z0-9]", "", tolower(x))

#' Cache directory for built annotation objects
#'
#' @description
#' Returns the directory in which [ens_build_annotation()] stores its results
#' (`<db_name>.rds` and `<db_name>.tsv.gz`) and where [ens_load_annotation()]
#' and [ens_list_cached()] look by default.
#'
#' @details
#' The default is a per-user location chosen by [tools::R_user_dir()]. Override
#' it for the whole session (or in `.Rprofile`) with
#' `options(ensannot.cache_dir = "/path/to/dir")`, or per call through the
#' `cache_dir` argument of the functions that take one.
#'
#' @return A single character string (a path). The directory is not created
#'   by this function.
#' @seealso [ens_build_annotation()], [ens_list_cached()]
#' @examples
#' ens_cache_dir()
#'
#' old <- options(ensannot.cache_dir = file.path(tempdir(), "ens_cache"))
#' ens_cache_dir()
#' options(old)
#' @export
ens_cache_dir <- function() {
  getOption("ensannot.cache_dir", tools::R_user_dir("ensannot", "cache"))
}

# Ensembl core DB names:  <genus_species>[_<assembly>]_core_<release>_<number>
#   homo_sapiens_core_115_38
#   capra_hircus_blackbengal_core_115_1
#   arabidopsis_thaliana_core_62_115_11      (Ensembl Genomes: two release numbers)
.parse_db_name <- function(db_name) {
  m <- regmatches(db_name, regexec("^(.+)_core_([0-9_]+)$", db_name))[[1]]
  if (length(m) != 3) stop("Cannot parse database name: ", db_name, call. = FALSE)
  parts <- strsplit(m[3], "_", fixed = TRUE)[[1]]
  list(species = m[2],                                   # genus_species[_assembly]
       release = paste(head(parts, -1), collapse = "_"),
       number  = tail(parts, 1))
}

# ----------------------------------------------------------- FTP navigation
# names of the entries of an Apache-style directory listing
.ens_list_dir <- function(url) {
  url  <- sub("/*$", "/", url)
  html <- tryCatch(readLines(url, warn = FALSE),
                   error = function(e)
                     stop("Could not list ", url, ": ", conditionMessage(e), call. = FALSE))
  hits  <- regmatches(html, gregexpr('href="[^"]+"', html))
  hrefs <- sub('"$', "", sub('^href="', "", unlist(hits)))
  hrefs <- hrefs[!grepl("^[?/#]|^https?:|^ftp:|^\\.\\.", hrefs)]
  unique(sub("/$", "", hrefs))
}

#' Find the Ensembl core-database folder for a species / assembly / release
#'
#' @description
#' Lists the `mysql` directory of an Ensembl release on the FTP site and returns
#' the URL of the core database that matches the requested species (and
#' assembly). [ens_build_annotation()] calls this for you; use it directly to
#' check which database will be used before downloading anything.
#'
#' @details
#' Ensembl names its core databases
#' `<genus_species>[_<assembly>]_core_<release>_<number>`:
#'
#' \itemize{
#'   \item `homo_sapiens_core_115_38` - no assembly part: `assembly = NULL`.
#'   \item `capra_hircus_blackbengal_core_115_1` - assembly part
#'     `blackbengal`: `species = "capra_hircus", assembly = "blackbengal"`.
#' }
#'
#' With `assembly = NULL` only the plain `<genus_species>_core_...` database
#' matches; with an assembly given, only `<genus_species>_<assembly>_core_...`
#' matches. This guarantees you get the assembly you asked for. Matching
#' ignores case and punctuation, so `"Black_Bengal"` and `"blackbengal"` are
#' equivalent.
#'
#' **Finding assembly names:** browse
#' `https://ftp.ensembl.org/pub/release-<N>/mysql/` (or `current_mysql/`) and
#' read the folder names. The assembly name is the text between the species and
#' `_core`. If a lookup fails, the error message lists the databases that exist
#' for that species.
#'
#' @param species Organism as `genus_species`, e.g. `"homo_sapiens"`
#'   (`"Homo sapiens"` is accepted too).
#' @param release Ensembl release number (e.g. `115`) or `"current"` for the
#'   latest release.
#' @param assembly Assembly name exactly as it appears in the FTP folder name
#'   (case and punctuation ignored), or `NULL` (default) for databases without
#'   an assembly part.
#' @param base_url Root of the FTP/HTTPS site. Change for mirrors or other
#'   Ensembl divisions (e.g. Ensembl Genomes).
#'
#' @return A single character string: the URL of the database directory, with
#'   a trailing slash.
#' @seealso [ens_build_annotation()]
#' @examples
#' \dontrun{
#' # Plain species, no assembly part: homo_sapiens_core_115_38
#' ens_resolve_db_url("homo_sapiens", release = 115)
#'
#' # Species with an assembly part: capra_hircus_blackbengal_core_115_1
#' ens_resolve_db_url("capra_hircus", release = 115, assembly = "blackbengal")
#'
#' # Latest release
#' ens_resolve_db_url("mus_musculus", release = "current")
#' }
#' @export
ens_resolve_db_url <- function(species, release = "current", assembly = NULL,
                               base_url = .ENS_FTP) {
  sp <- .norm_species(species)
  rel_dir <- if (identical(tolower(as.character(release)), "current"))
    "current_mysql" else sprintf("release-%s/mysql", release)
  dir_url <- sprintf("%s/%s/", sub("/+$", "", base_url), rel_dir)

  entries <- .ens_list_dir(dir_url)
  dbs     <- grep("^.+_core_[0-9_]+$", entries, value = TRUE)
  prefix  <- sub("_core_[0-9_]+$", "", dbs)
  want    <- if (is.null(assembly)) sp else paste0(sp, "_", assembly)

  hit <- dbs[.key(prefix) == .key(want)]
  if (!length(hit)) {
    related <- dbs[startsWith(.key(prefix), .key(sp))]
    stop("No core database named '", want, "_core_<release>_<n>' in ", dir_url, "\n",
         if (length(related))
           paste0("Databases for this species:\n  ", paste(related, collapse = "\n  "),
                  "\nSet 'assembly' to the text between the species and '_core'.\n")
         else "No database starts with that species name; check the spelling.\n",
         "Browse ", dir_url, " for the exact names, or pass ftp_url= directly.",
         call. = FALSE)
  }
  if (length(hit) > 1)
    stop("Several databases match:\n  ", paste(hit, collapse = "\n  "),
         "\nPass ftp_url= directly.", call. = FALSE)
  paste0(dir_url, hit, "/")
}

# download to a temp file first (no half-written files); skips existing non-empty files
.ens_download <- function(url, dest, verbose = TRUE) {
  if (file.exists(dest) && file.size(dest) > 0) return(invisible(dest))
  tmp <- paste0(dest, ".part")
  old <- options(timeout = max(getOption("timeout"), 7200)); on.exit(options(old), add = TRUE)
  .msg(verbose, "downloading ", basename(dest))
  res <- tryCatch(utils::download.file(url, tmp, mode = "wb", quiet = !verbose),
                  error = function(e) e)
  if (inherits(res, "error") || !file.exists(tmp))
    stop("Download failed: ", url, if (inherits(res, "error")) paste0(" (", conditionMessage(res), ")"),
         call. = FALSE)
  file.rename(tmp, dest)
  invisible(dest)
}

# ------------------------------------------------- reading MySQL dump tables
# Column names are parsed from the <db>.sql.gz CREATE TABLE statements, so the
# code survives column additions between Ensembl releases.
.sql_columns <- function(sql_gz, table) {
  con <- gzfile(sql_gz, "rt"); on.exit(close(con))
  lines <- readLines(con, warn = FALSE)
  start <- grep(sprintf("^CREATE TABLE `%s` \\(", table), lines)
  if (!length(start)) stop("Table '", table, "' not found in ", basename(sql_gz), call. = FALSE)
  cols <- character()
  for (ln in lines[-seq_len(start[1])]) {
    if (grepl("^\\)", ln)) break
    m <- regmatches(ln, regexec("^\\s*`([^`]+)`", ln))[[1]]
    if (length(m) == 2) cols <- c(cols, m[2])
  }
  cols
}

.fread_gz <- function(file, select) {
  args <- list(header = FALSE, sep = "\t", quote = "", na.strings = "\\N",
               select = select, data.table = TRUE, showProgress = FALSE)
  if (requireNamespace("R.utils", quietly = TRUE)) {
    do.call(fread, c(list(file = file), args))
  } else {
    do.call(fread, c(list(cmd = paste("gzip -dc", shQuote(file))), args))
  }
}

.read_table <- function(raw_dir, db_name, table, cols) {
  all_cols <- .sql_columns(file.path(raw_dir, paste0(db_name, ".sql.gz")), table)
  miss <- setdiff(cols, all_cols)
  if (length(miss))
    stop("Table '", table, "' lacks expected column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  idx <- sort(match(cols, all_cols))
  dt  <- .fread_gz(file.path(raw_dir, paste0(table, ".txt.gz")), select = idx)
  setnames(dt, all_cols[idx])
  dt[]
}

# =============================================================================
# The annotation object
# =============================================================================

#' The `ens_annot` annotation object
#'
#' @description
#' The object returned by [ens_build_annotation()], [ens_load_annotation()] and
#' [ens_example_annotation()], and accepted by all `ens_annotate_*()`
#' functions. It is a list of class `"ens_annot"` with two elements.
#'
#' @section `$genes`:
#' A [data.table::data.table] with one row per current gene, sorted by
#' chromosome (1, 2, ..., X, Y, MT, scaffolds) and start. Column names follow
#' biomaRt attribute names:
#' \describe{
#'   \item{`ensembl_gene_id`}{Stable Ensembl gene ID, without version.}
#'   \item{`chromosome_name`}{Chromosome or scaffold name as in Ensembl
#'     (`"1"`, `"X"`, `"MT"`, ...).}
#'   \item{`start_position`, `end_position`}{1-based, inclusive gene coordinates.}
#'   \item{`strand`}{`1` or `-1` (integer), as returned by biomaRt.}
#'   \item{`external_gene_name`}{Gene symbol (the gene's display xref); `NA` if
#'     the gene has no symbol.}
#'   \item{`custom_gene_name`}{Never empty: the symbol if present, otherwise the
#'     first (alphabetical) synonym, otherwise the Ensembl gene ID.}
#'   \item{`description`}{Gene description. Keeps the trailing
#'     `[Source:...]` text unless built with `clean_description = TRUE`.}
#'   \item{`external_synonym`}{Synonyms from all external references linked to
#'     the gene, separated by `"; "`.}
#'   \item{`gene_biotype`}{e.g. `protein_coding`, `lncRNA`, `miRNA`.}
#'   \item{`entrezgene_id`}{Entrez Gene ID(s), separated by `"; "`.}
#'   \item{`ensembl_transcript_id`}{Transcript IDs, separated by `"; "`.}
#' }
#' Values that are multiple within one gene are joined with `"; "`. Missing
#' values are `NA`.
#'
#' @section `$meta`:
#' A list: `species` (the `genus_species[_assembly]` part of the database
#' name), `release`, `number`, `db_name`, `source_url`, `primary_only` and
#' `built` (timestamp).
#'
#' @seealso [ens_build_annotation()], [ens_load_annotation()],
#'   [ens_example_annotation()]
#' @examples
#' ann <- ens_example_annotation()
#' class(ann)
#' names(ann$genes)
#' ann$meta$db_name
#' ann
#' @name ens_annot
NULL

#' Tiny synthetic annotation for examples and tests
#'
#' @description
#' Returns a small, fully synthetic [ens_annot][ens_annot] object (five
#' fictional genes on chromosomes 1, 2 and X). It has the same structure as a
#' real build, so all examples and unit tests can run offline and instantly.
#' The gene IDs (`TOYG001`, ...), names and coordinates are invented and
#' do not correspond to any real organism.
#'
#' @details
#' The toy genes cover the situations the package has to handle: a gene with a
#' symbol and several synonyms and transcripts (`TOYG001`), a gene without a
#' symbol but with synonyms (`TOYG002`; `custom_gene_name` falls back to the
#' first synonym), a gene with neither (`TOYG003`; `custom_gene_name` is the
#' gene ID), genes on both strands, and genes on different chromosomes.
#'
#' @return An object of class `"ens_annot"`; see [ens_annot].
#' @examples
#' ann <- ens_example_annotation()
#' ann$genes[, c("ensembl_gene_id", "external_gene_name", "custom_gene_name")]
#' @export
ens_example_annotation <- function() {
  genes <- data.table(
    ensembl_gene_id    = c("TOYG001", "TOYG002", "TOYG003", "TOYG004", "TOYG005"),
    chromosome_name    = c("1", "1", "1", "2", "X"),
    start_position     = c(1000, 6000, 20000, 500, 100),
    end_position       = c(5000, 9000, 22000, 3000, 900),
    strand             = c(1L, -1L, 1L, 1L, -1L),
    external_gene_name = c("ALPHA1", NA, NA, "GAMMA3", "DELTA4"),
    custom_gene_name   = c("ALPHA1", "B2", "TOYG003", "GAMMA3", "DELTA4"),
    description        = c("toy alpha gene [Source:Toy;Acc:1]",
                           "toy beta-like gene [Source:Toy;Acc:2]", NA,
                           "toy gamma gene [Source:Toy;Acc:4]",
                           "toy delta gene [Source:Toy;Acc:5]"),
    external_synonym   = c("A1; ALP", "B2; BETA-LIKE", NA, NA, "D4"),
    gene_biotype       = c("protein_coding", "protein_coding", "lncRNA",
                           "protein_coding", "miRNA"),
    entrezgene_id      = c("101", "102", NA, "104", NA),
    ensembl_transcript_id = c("TOYT001; TOYT002", "TOYT003", "TOYT004",
                              "TOYT005", "TOYT006")
  )
  structure(list(
    genes = genes,
    meta  = list(species = "toy_species", release = "0", number = "0",
                 db_name = "toy_species_core_0_0", source_url = NA_character_,
                 primary_only = TRUE, built = Sys.time())
  ), class = "ens_annot")
}

# =============================================================================
# PART 1 - build / cache / load
# =============================================================================

#' Build (or reuse) a gene annotation table from the Ensembl FTP site
#'
#' @description
#' Downloads the core database dump files of one Ensembl species / assembly /
#' release, parses them into one tidy table (one row per gene) and caches the
#' result so it only has to be built once. The dump files are plain
#' tab-delimited text, so no SQL server is involved.
#'
#' Identify the database either with `species` (+ `release`, + `assembly`) or
#' with the exact folder URL in `ftp_url`.
#'
#' @details
#' **Which database is used.** Ensembl names core databases
#' `<genus_species>[_<assembly>]_core_<release>_<number>`.
#' By default (`assembly = NULL`) the database `<genus_species>_core_...` is
#' used. If `assembly` is given, `<genus_species>_<assembly>_core_...` is used,
#' which makes the assembly choice explicit. Look up assembly names in the
#' folder listing of
#' `https://ftp.ensembl.org/pub/release-<N>/mysql/` (see
#' [ens_resolve_db_url()]).
#'
#' **What is downloaded and parsed.** The table definition file plus the tables
#' `gene`, `transcript`, `translation`, `seq_region`, `coord_system`,
#' `seq_region_attrib`, `attrib_type`, `xref`, `external_db`,
#' `external_synonym` and `object_xref`. Column positions are read from the
#' table definitions, so the parser adapts to schema changes between releases.
#'
#' **How the columns are derived.**
#' \itemize{
#'   \item Only current genes (`is_current = 1`) are kept. With
#'     `primary_only = TRUE` genes on non-reference regions (patches,
#'     haplotypes) and LRG are dropped.
#'   \item `external_gene_name`: the gene's display xref label.
#'   \item `external_synonym`: synonyms of all external references linked to the
#'     gene (including its display xref), so genes without an HGNC-style symbol
#'     can still get synonyms (e.g. from Entrez or other external databases).
#'   \item `entrezgene_id`: `EntrezGene` references linked to the gene, its
#'     transcripts or its translations.
#'   \item `ensembl_transcript_id`: all current transcripts of the gene.
#'   \item `custom_gene_name`: symbol, else first synonym (alphabetical), else
#'     the Ensembl gene ID.
#' }
#'
#' **Caching.** Results are written to `cache_dir` as `<db_name>.rds` (used for
#' reloading) and `<db_name>.tsv.gz` (for use outside R). If the cache file
#' exists it is returned immediately unless `force = TRUE`. When `release` is a
#' number and `species` is given, the cache is checked without any network
#' access.
#'
#' **Resources.** The first build downloads the dumps (for human, several
#' hundred MB to a few GB, dominated by `object_xref`) and needs a few GB of
#' RAM while parsing; installing the \pkg{R.utils} package makes reading the
#' `.gz` files much faster. The raw files are deleted after the build unless
#' `keep_raw = TRUE`. Subsequent loads take seconds.
#'
#' @param species Organism as `genus_species` (e.g. `"homo_sapiens"`). Not
#'   needed if `ftp_url` is given.
#' @param release Ensembl release number (e.g. `115`) or `"current"`.
#' @param assembly Assembly name as it appears in the FTP folder name (the text
#'   between the species and `_core`), e.g. `"blackbengal"` for
#'   `capra_hircus_blackbengal_core_115_1`. `NULL` (default) selects the
#'   database without an assembly part.
#' @param ftp_url URL of the core database folder, e.g.
#'   `"https://ftp.ensembl.org/pub/release-115/mysql/homo_sapiens_core_115_38/"`.
#'   Overrides `species`, `release` and `assembly`.
#' @param cache_dir Directory for cached objects; default [ens_cache_dir()].
#' @param primary_only If `TRUE` (default), drop genes on patch/haplotype
#'   (`non_ref`) regions and LRG records. Unplaced scaffolds are kept.
#' @param clean_description If `FALSE` (default, as in biomaRt) `description`
#'   keeps its trailing `[Source:...]` tag; if `TRUE` the tag is removed.
#' @param force Rebuild even if a cached object exists.
#' @param keep_raw Keep the downloaded raw dump files (in
#'   `<cache_dir>/raw/<db_name>`) instead of deleting them after the build.
#' @param base_url FTP/HTTPS root used when resolving `species`; change for
#'   mirrors or Ensembl Genomes.
#' @param verbose Print progress messages.
#'
#' @return An object of class `"ens_annot"` (see [ens_annot]) holding
#'   `$genes`, the annotation table, and `$meta`.
#' @seealso [ens_load_annotation()], [ens_list_cached()],
#'   [ens_resolve_db_url()], [ens_annotate_ids()], [ens_annotate_regions()],
#'   [ens_annotate_loci()]
#' @examples
#' \dontrun{
#' # Default: homo_sapiens_core_115_38 (no assembly part in the folder name)
#' ann <- ens_build_annotation("homo_sapiens", release = 115)
#'
#' # A species whose folder name carries an assembly part:
#' # capra_hircus_blackbengal_core_115_1
#' ann <- ens_build_annotation("capra_hircus", release = 115,
#'                             assembly = "blackbengal")
#'
#' # Give the exact folder instead
#' ann <- ens_build_annotation(ftp_url =
#'   "https://ftp.ensembl.org/pub/release-115/mysql/homo_sapiens_core_115_38/")
#'
#' # Latest release, cleaned descriptions, custom cache location
#' ann <- ens_build_annotation("mus_musculus", release = "current",
#'                             clean_description = TRUE,
#'                             cache_dir = "~/ensembl_cache")
#'
#' # Second call is instant: reads the cache
#' ann <- ens_build_annotation("homo_sapiens", release = 115)
#'
#' ann                 # summary
#' head(ann$genes)
#' }
#' @export
ens_build_annotation <- function(species = NULL, release = "current", assembly = NULL,
                                 ftp_url = NULL, cache_dir = ens_cache_dir(),
                                 primary_only = TRUE, clean_description = FALSE,
                                 force = FALSE, keep_raw = FALSE,
                                 base_url = .ENS_FTP, verbose = TRUE) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  # cheap offline cache hit for explicit releases
  if (!force && is.null(ftp_url) && !is.null(species) &&
      !identical(tolower(as.character(release)), "current")) {
    hit <- .find_cached(species, release, assembly, cache_dir)
    if (!is.null(hit)) { .msg(verbose, "using cached ", basename(hit)); return(.read_annot_file(hit)) }
  }

  if (is.null(ftp_url)) {
    if (is.null(species)) stop("Provide 'species' or 'ftp_url'.", call. = FALSE)
    ftp_url <- ens_resolve_db_url(species, release, assembly, base_url)
  }
  ftp_url <- sub("/*$", "/", ftp_url)
  db_name <- basename(sub("/+$", "", ftp_url))
  info    <- .parse_db_name(db_name)
  rds     <- file.path(cache_dir, paste0(db_name, ".rds"))

  if (!force && file.exists(rds)) {
    .msg(verbose, "using cached ", basename(rds)); return(readRDS(rds))
  }

  raw_dir <- file.path(cache_dir, "raw", db_name)
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  tables <- c("seq_region", "coord_system", "seq_region_attrib", "attrib_type",
              "gene", "transcript", "translation", "xref", "external_db",
              "external_synonym", "object_xref")
  for (f in c(paste0(db_name, ".sql.gz"), paste0(tables, ".txt.gz")))
    .ens_download(paste0(ftp_url, f), file.path(raw_dir, f), verbose)

  rd <- function(tb, cols) { .msg(verbose, "reading ", tb); .read_table(raw_dir, db_name, tb, cols) }

  # --- genes + chromosome names -------------------------------------------
  gene <- rd("gene", c("gene_id", "biotype", "seq_region_id", "seq_region_start",
                       "seq_region_end", "seq_region_strand", "display_xref_id",
                       "description", "is_current", "stable_id"))
  gene <- gene[is_current == 1L]
  setnames(gene, c("gene_id", "stable_id"), c("internal_id", "ensembl_gene_id"))

  sr <- rd("seq_region", c("seq_region_id", "name", "coord_system_id"))
  sr[, name := as.character(name)]
  if (primary_only) {
    cs    <- rd("coord_system", c("coord_system_id", "name"))
    atype <- rd("attrib_type", c("attrib_type_id", "code"))
    sra   <- rd("seq_region_attrib", c("seq_region_id", "attrib_type_id"))
    nonref_ids <- sra[attrib_type_id %in% atype[code == "non_ref", attrib_type_id], unique(seq_region_id)]
    lrg_cs     <- cs[tolower(name) == "lrg", coord_system_id]
    drop_ids   <- union(nonref_ids, sr[coord_system_id %in% lrg_cs, seq_region_id])
    gene <- gene[!seq_region_id %in% drop_ids]
  }
  gene[, chromosome_name := sr$name[match(seq_region_id, sr$seq_region_id)]]

  # --- description ----------------------------------------------------------
  gene[, description := as.character(description)]
  if (clean_description)
    gene[, description := trimws(sub("\\s*\\[Source:[^]]*\\]\\s*$", "", description))]
  gene[!is.na(description) & !nzchar(description), description := NA_character_]

  # --- gene symbol (display xref) -----------------------------------------
  xref  <- rd("xref", c("xref_id", "external_db_id", "dbprimary_acc", "display_label"))
  ext   <- rd("external_db", c("external_db_id", "db_name"))
  synon <- rd("external_synonym", c("xref_id", "synonym"))

  gene[, external_gene_name := as.character(xref$display_label[match(display_xref_id, xref$xref_id)])]
  gene[!is.na(external_gene_name) & !nzchar(trimws(external_gene_name)),
       external_gene_name := NA_character_]

  # --- transcripts --------------------------------------------------------
  tr <- rd("transcript", c("transcript_id", "gene_id", "stable_id", "is_current"))
  tr <- tr[is_current == 1L & gene_id %in% gene$internal_id]
  tr_agg <- tr[, .(ensembl_transcript_id = .collapse_unique(sort(stable_id))), by = gene_id]
  gene[, ensembl_transcript_id := tr_agg$ensembl_transcript_id[match(internal_id, tr_agg$gene_id)]]

  # --- synonyms + Entrez (via object_xref) --------------------------------
  # Synonyms are collected from ALL xrefs linked to the gene (plus its display
  # xref), so genes without an HGNC-style symbol can still get synonyms
  # (e.g. from Entrez / other external databases).
  entrez_db <- ext[db_name == "EntrezGene", external_db_id]
  ex   <- xref[external_db_id %in% entrez_db, .(xref_id, entrez = as.character(dbprimary_acc))]
  syn_x <- unique(synon$xref_id)

  ox <- rd("object_xref", c("ensembl_id", "ensembl_object_type", "xref_id"))
  ox <- ox[xref_id %in% union(syn_x, ex$xref_id)]

  gl <- unique(rbindlist(list(
    ox[ensembl_object_type == "Gene" & xref_id %in% syn_x,
       .(internal_id = ensembl_id, xref_id)],
    gene[!is.na(display_xref_id) & display_xref_id %in% syn_x,
         .(internal_id, xref_id = display_xref_id)]
  )))
  sj  <- merge(gl, synon[, .(xref_id, synonym = as.character(synonym))],
               by = "xref_id", allow.cartesian = TRUE)
  syn <- sj[, .(external_synonym = .collapse_unique(sort(unique(synonym)))), by = internal_id]
  gene[, external_synonym := syn$external_synonym[match(internal_id, syn$internal_id)]]
  rm(gl, sj, syn)

  gene[, entrezgene_id := NA_character_]
  if (nrow(ex)) {
    oe <- ox[xref_id %in% ex$xref_id]
    oe[, entrez := ex$entrez[match(xref_id, ex$xref_id)]]
    g_link <- oe[ensembl_object_type == "Gene", .(internal_id = ensembl_id, entrez)]
    t_link <- oe[ensembl_object_type == "Transcript",
                 .(internal_id = tr$gene_id[match(ensembl_id, tr$transcript_id)], entrez)]
    tl     <- rd("translation", c("translation_id", "transcript_id"))
    p_link <- oe[ensembl_object_type == "Translation",
                 .(internal_id = tr$gene_id[match(tl$transcript_id[match(ensembl_id, tl$translation_id)],
                                                  tr$transcript_id)], entrez)]
    all_link <- unique(rbindlist(list(g_link, t_link, p_link))[!is.na(internal_id) & !is.na(entrez)])
    ent <- all_link[, .(entrezgene_id = .collapse_unique(entrez[order(suppressWarnings(as.numeric(entrez)))])),
                    by = internal_id]
    gene[, entrezgene_id := ent$entrezgene_id[match(internal_id, ent$internal_id)]]
    rm(oe, g_link, t_link, p_link, tl, all_link, ent)
  }
  rm(ox)

  # --- custom gene name: symbol -> first synonym -> Ensembl gene ID ---------
  first_syn <- trimws(sub(";.*$", "", gene$external_synonym))
  gene[, custom_gene_name := fifelse(!is.na(external_gene_name), external_gene_name,
                                     fifelse(!is.na(first_syn) & nzchar(first_syn), first_syn,
                                             ensembl_gene_id))]

  # --- final tidy table (biomaRt attribute names) ---------------------------
  genes <- gene[, .(
    ensembl_gene_id,
    chromosome_name,
    start_position = as.numeric(seq_region_start),
    end_position   = as.numeric(seq_region_end),
    strand         = as.integer(seq_region_strand),
    external_gene_name,
    custom_gene_name,
    description,
    external_synonym,
    gene_biotype   = biotype,
    entrezgene_id,
    ensembl_transcript_id
  )]
  genes[, .chr_num := suppressWarnings(as.integer(chromosome_name))]
  setorder(genes, .chr_num, chromosome_name, start_position, na.last = TRUE)  # 1..22, X, Y, MT, scaffolds
  genes[, .chr_num := NULL]

  obj <- structure(list(
    genes = genes,
    meta  = list(species = info$species, release = info$release, number = info$number,
                 db_name = db_name, source_url = ftp_url,
                 primary_only = primary_only, built = Sys.time())
  ), class = "ens_annot")

  saveRDS(obj, rds)
  fwrite(genes, file.path(cache_dir, paste0(db_name, ".tsv.gz")), sep = "\t", na = "")
  if (!keep_raw) unlink(raw_dir, recursive = TRUE)
  .msg(verbose, "done: ", nrow(genes), " genes -> ", rds)
  obj
}

# newest cached .rds matching species / release / assembly, or NULL
.find_cached <- function(species, release, assembly, cache_dir) {
  files <- list.files(cache_dir, "_core_.*\\.rds$", full.names = TRUE)
  if (!length(files)) return(NULL)
  want <- .norm_species(species)
  if (!is.null(assembly)) want <- paste0(want, "_", assembly)
  inf <- lapply(basename(sub("\\.rds$", "", files)), .parse_db_name)
  ok  <- vapply(inf, function(i) .key(i$species) == .key(want), logical(1))
  if (!identical(tolower(as.character(release)), "current"))
    ok <- ok & vapply(inf, function(i) i$release == as.character(release), logical(1))
  files <- files[ok]
  if (!length(files)) return(NULL)
  if (length(files) > 1) warning("Several cached matches; using the most recent.", call. = FALSE)
  files[which.max(file.mtime(files))]
}

.read_annot_file <- function(file) {
  if (grepl("\\.rds$", file, ignore.case = TRUE)) return(readRDS(file))
  genes <- fread(file, sep = "\t", data.table = TRUE, colClasses = list(character = "chromosome_name"))
  info <- tryCatch(.parse_db_name(sub("\\.tsv(\\.gz)?$", "", basename(file))), error = function(e) NULL)
  structure(list(genes = genes, meta = c(info, list(source_url = NA, file = file))), class = "ens_annot")
}

#' Load a previously built annotation
#'
#' @description
#' Loads an annotation built earlier by [ens_build_annotation()], either by
#' looking it up in the cache (by species / release / assembly) or from an
#' explicit file. No network access is needed.
#'
#' @details
#' If several cached objects match (for example `release = "current"` after
#' building two different releases), the most recently written one is used and
#' a warning is issued. Give an explicit `release` to be unambiguous.
#'
#' A `.tsv` / `.tsv.gz` file written by [ens_build_annotation()] can also be
#' loaded with `file =`; this is useful to share the table with colleagues who
#' do not have the cache. Column types are restored from the file, so the
#' chromosome column stays character.
#'
#' @param species Organism as `genus_species`; not needed if `file` is given.
#' @param release Ensembl release number or `"current"` (any cached release).
#' @param assembly Assembly name as in the FTP folder name, or `NULL` for
#'   databases without an assembly part.
#' @param file Path to a saved `.rds` or `.tsv`/`.tsv.gz` annotation. If given,
#'   the other lookup arguments are ignored.
#' @param cache_dir Cache directory; default [ens_cache_dir()].
#'
#' @return An object of class `"ens_annot"`; see [ens_annot].
#' @seealso [ens_build_annotation()], [ens_list_cached()]
#' @examples
#' # Round trip with the synthetic example data
#' ann <- ens_example_annotation()
#' f <- tempfile(fileext = ".rds")
#' saveRDS(ann, f)
#' ann2 <- ens_load_annotation(file = f)
#' identical(ann$genes, ann2$genes)
#'
#' \dontrun{
#' # From the cache, after ens_build_annotation() has been run once
#' ann <- ens_load_annotation("homo_sapiens", release = 115)
#' ann <- ens_load_annotation("capra_hircus", release = 115, assembly = "blackbengal")
#'
#' # From a shared table
#' ann <- ens_load_annotation(file = "homo_sapiens_core_115_38.tsv.gz")
#' }
#' @export
ens_load_annotation <- function(species = NULL, release = "current", assembly = NULL,
                                file = NULL, cache_dir = ens_cache_dir()) {
  if (!is.null(file)) return(.read_annot_file(file))
  if (is.null(species)) stop("Provide 'species' or 'file'.", call. = FALSE)
  hit <- .find_cached(species, release, assembly, cache_dir)
  if (is.null(hit)) stop("No cached annotation found; run ens_build_annotation() first.", call. = FALSE)
  .read_annot_file(hit)
}

#' List cached annotation objects
#'
#' @description
#' Shows which annotation objects have been built and cached, so you can see
#' what is available offline.
#'
#' @param cache_dir Cache directory; default [ens_cache_dir()].
#'
#' @return A [data.table::data.table] with one row per cached object and
#'   columns `species_assembly` (the `genus_species[_assembly]` part of the
#'   database name), `release`, `number` (trailing number of the database
#'   name), `file`, `size_mb` and `modified`. Empty (zero columns) if nothing is
#'   cached.
#' @seealso [ens_build_annotation()], [ens_load_annotation()]
#' @examples
#' # An empty temporary cache has nothing in it
#' ens_list_cached(cache_dir = tempfile())
#'
#' \dontrun{
#' ens_list_cached()
#' }
#' @export
ens_list_cached <- function(cache_dir = ens_cache_dir()) {
  files <- list.files(cache_dir, "_core_.*\\.rds$", full.names = TRUE)
  if (!length(files)) return(data.table())
  rbindlist(lapply(files, function(f) {
    i <- .parse_db_name(sub("\\.rds$", "", basename(f)))
    data.table(species_assembly = i$species, release = i$release, number = i$number,
               file = f, size_mb = round(file.size(f) / 1e6, 1), modified = file.mtime(f))
  }))
}

#' Print an `ens_annot` object
#'
#' @description
#' Prints a short summary: species/assembly, release, number of genes and
#' chromosomes, the column names and the most frequent gene biotypes.
#'
#' @param x An [ens_annot] object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @examples
#' ens_example_annotation()
#' @export
print.ens_annot <- function(x, ...) {
  m <- x$meta
  cat("<ens_annot>", m$species, "| release", m$release, "| db number", m$number, "\n")
  cat("  genes:", nrow(x$genes), " chromosomes/regions:", uniqueN(x$genes$chromosome_name), "\n")
  cat("  columns:", paste(names(x$genes), collapse = ", "), "\n")
  top <- head(sort(table(x$genes$gene_biotype), decreasing = TRUE), 5)
  cat("  top biotypes:", paste(names(top), top, sep = "=", collapse = ", "), "\n")
  invisible(x)
}

# =============================================================================
# PART 2 - annotating user tables
# =============================================================================

# data.frame / data.table / path to a delimited file -> data.table
.as_input <- function(x) {
  if (is.character(x) && length(x) == 1 && file.exists(x)) return(fread(x, data.table = TRUE))
  as.data.table(x)
}

# ens_annot / path to saved annotation -> validated ens_annot
.as_annot <- function(a) {
  if (is.character(a) && length(a) == 1 && file.exists(a)) a <- .read_annot_file(a)
  if (!inherits(a, "ens_annot"))
    stop("'annot' must be an ens_annot object or a path to a saved .rds/.tsv(.gz).", call. = FALSE)
  if (!"chromosome_name" %in% names(a$genes))
    stop("This annotation was built with an older version of ensannot.R ",
         "(old column names). Rebuild with ens_build_annotation(..., force = TRUE).", call. = FALSE)
  a
}

.check_cols <- function(annot, cols) {
  bad <- setdiff(cols, names(annot$genes))
  if (length(bad)) stop("Unknown annotation column(s): ", paste(bad, collapse = ", "),
                        "\nAvailable: ", paste(names(annot$genes), collapse = ", "), call. = FALSE)
  cols
}

# bind input and result; rename result columns that clash with input columns
.finish <- function(x, res, out_file = NULL) {
  dup <- names(res) %in% names(x)
  if (any(dup)) setnames(res, names(res)[dup], paste0(names(res)[dup], ".ens"))
  out <- cbind(x, res)
  if (!is.null(out_file)) fwrite(out, out_file, sep = "\t", na = "")
  out[]
}

.query_table <- function(x, chr_col, start_col, end_col) {
  for (cc in c(chr_col, start_col, end_col))
    if (!cc %in% names(x)) stop("Column '", cc, "' not found in input.", call. = FALSE)
  s <- suppressWarnings(as.numeric(x[[start_col]]))
  e <- suppressWarnings(as.numeric(x[[end_col]]))
  lo <- pmin(s, e); hi <- pmax(s, e)
  chr <- .norm_chr(x[[chr_col]])
  list(chr = chr, s = lo, e = hi, valid = !is.na(chr) & !is.na(lo) & !is.na(hi))
}

# all (query row, gene row) pairs whose intervals come within `pad` bp
.overlap_pairs <- function(q, genes, pad = 0) {
  G <- data.table(chr = .norm_chr(genes$chromosome_name), gs = as.numeric(genes$start_position),
                  ge = as.numeric(genes$end_position), gi = seq_len(nrow(genes)))
  setkey(G, chr, gs, ge)
  v <- which(q$valid)
  Q <- data.table(qi = v, chr = q$chr[v], qs = pmax(q$s[v] - pad, 1), qe = q$e[v] + pad)
  setkey(Q, chr, qs, qe)
  ov <- foverlaps(Q, G, by.x = c("chr", "qs", "qe"), by.y = c("chr", "gs", "ge"),
                  type = "any", nomatch = NULL)
  ov[, .(qi, gi)]
}

# na_string = NULL: missing values stay NA and are simply left out when collapsing
.gene_vals <- function(genes, gi, cols, inner_sep, na_string = NULL) {
  out <- lapply(cols, function(cn) {
    v <- as.character(genes[[cn]][gi])
    v <- gsub("; ", inner_sep, v, fixed = TRUE)
    if (!is.null(na_string)) v[is.na(v)] <- na_string
    v
  })
  names(out) <- cols
  out
}

# collapse per-pair values to one string per input row; NAs are dropped, and a
# row with nothing left stays NA (written as empty)
.collapse_by_row <- function(n, qi, vals, sep) {
  if (!length(qi)) {
    return(as.data.table(lapply(vals, function(v) rep(NA_character_, n))))
  }
  dt  <- as.data.table(c(list(qi = qi), vals))
  agg <- dt[, lapply(.SD, function(v) {
    v <- v[!is.na(v)]
    if (length(v)) paste(v, collapse = sep) else NA_character_
  }), by = qi]
  out <- agg[match(seq_len(n), agg$qi)]
  out[, qi := NULL]
  out[]
}

#' Annotate Ensembl gene IDs
#'
#' @description
#' Mode 1 of the annotation functions: looks up each Ensembl gene ID of your
#' input in the annotation and appends the requested columns. The output has
#' exactly one row per input row, in the input order.
#'
#' @details
#' IDs are matched after trimming whitespace and, by default, removing a
#' version suffix (`ENSG00000139618.15` matches `ENSG00000139618`). IDs that are
#' not in the annotation (a different release or species, retired genes,
#' typos) get `NA` in all added columns; the number of unmatched IDs is
#' reported with a message. Duplicated IDs in the input are fine.
#'
#' If an added column has the same name as a column already present in your
#' input, the added column gets the suffix `.ens`.
#'
#' @param x The input table: a data.frame / data.table, or the path to a
#'   delimited file (read with [data.table::fread()]).
#' @param annot An [ens_annot] object from [ens_build_annotation()] /
#'   [ens_load_annotation()], or the path to a saved `.rds` / `.tsv(.gz)`.
#' @param id_col Name of the column of `x` holding Ensembl gene IDs.
#' @param add_cols Annotation columns to append. Default: everything except
#'   `ensembl_gene_id` (so `external_gene_name`, `custom_gene_name`,
#'   `description`, `external_synonym`, `gene_biotype`, `entrezgene_id`,
#'   `ensembl_transcript_id`, `chromosome_name`, `start_position`,
#'   `end_position`, `strand`). See [ens_annot] for the full list.
#' @param strip_version Remove a trailing `.<number>` from the input IDs before
#'   matching.
#' @param out_file Optional path; if given, the result is also written there as
#'   a tab-delimited file.
#'
#' @return A [data.table::data.table] with all columns of `x` followed by the
#'   requested annotation columns.
#' @seealso [ens_annotate_regions()], [ens_annotate_loci()], [ens_annot]
#' @examples
#' ann <- ens_example_annotation()
#' x <- data.frame(gene = c("TOYG001", "TOYG003.2", "NOT_A_GENE"),
#'                 logFC = c(1.5, -2.1, 0.3))
#'
#' # everything
#' ens_annotate_ids(x, ann, id_col = "gene")
#'
#' # only selected columns
#' ens_annotate_ids(x, ann, id_col = "gene",
#'                  add_cols = c("custom_gene_name", "gene_biotype"))
#'
#' # write to file as well
#' ens_annotate_ids(x, ann, "gene", out_file = tempfile(fileext = ".tsv"))
#' @export
ens_annotate_ids <- function(x, annot, id_col,
                             add_cols = c(.ENS_ATTR_COLS[-1], .ENS_COORD_COLS),
                             strip_version = TRUE, out_file = NULL) {
  x <- .as_input(x); annot <- .as_annot(annot); add_cols <- .check_cols(annot, add_cols)
  if (!id_col %in% names(x)) stop("Column '", id_col, "' not found in input.", call. = FALSE)
  ids <- trimws(as.character(x[[id_col]]))
  if (strip_version) ids <- sub("\\.[0-9]+$", "", ids)
  m   <- match(ids, annot$genes$ensembl_gene_id)
  message(sum(is.na(m) & !is.na(ids)), " of ", length(ids), " IDs not found in annotation.")
  res <- annot$genes[m, add_cols, with = FALSE]
  .finish(x, res, out_file)
}

#' Annotate genomic regions with all overlapping genes
#'
#' @description
#' Mode 2 of the annotation functions: for every region (row) of your input,
#' finds **all** genes that overlap it. Returns a list with two tables:
#' \describe{
#'   \item{`regions`}{your input rows (one per input row, same order) with the
#'     gene information of all overlapping genes collapsed into the row.}
#'   \item{`genes`}{one row per distinct gene found in any of the regions, with
#'     the full gene details, in the layout of [ens_annotate_ids()].}
#' }
#'
#' @details
#' **Two tables.** `regions` is the compact per-region view (several genes
#' joined into one cell). `genes` is the long view: each overlapping gene
#' appears once, even if it falls in several regions, and has a `region_row`
#' column listing the input row number(s) it was found in, so the two tables
#' can be linked. The `genes` table is produced by calling [ens_annotate_ids()]
#' on the gene IDs found, so everything that function accepts can be passed
#' through `gene_args` (e.g. `add_cols`, `strip_version`, `out_file`).
#'
#' **Overlap.** A gene is reported when it overlaps the region by at least
#' `min_overlap` base pairs (both ends inclusive). Coordinates are 1-based, as
#' in Ensembl. If `start > end` in a row, the two are swapped. Rows with missing
#' or non-numeric coordinates get `n_genes = NA` and empty annotation.
#' Chromosome names are normalised: `chr1`, `Chr1` and `1` are equivalent, and
#' `chrM`/`M`/`MT` all match `MT`. Genes are matched on the Ensembl chromosome
#' name, so scaffold names must agree with the Ensembl ones.
#'
#' **Output layout.** Genes in a row are listed in genomic order and joined with
#' `sep`. Values that are multiple inside one gene (synonyms, Entrez IDs,
#' transcripts) are joined with `inner_sep` instead, so the two levels can be
#' told apart. Missing values (for example a gene with no synonym) are simply
#' left out rather than written as `"NA"`; because `ensembl_gene_id` and
#' `custom_gene_name` are never missing, they always list every gene in the
#' region. Set `na_string` (e.g. `"NA"`) to write a placeholder for missing
#' values, which keeps all columns aligned gene by gene. A row without any
#' overlapping gene has `n_genes = 0` and empty (NA) annotation columns.
#'
#' If an added column has the same name as a column already present in your
#' input, the added column gets the suffix `.ens`.
#'
#' @param x The input table: a data.frame / data.table, or the path to a
#'   delimited file.
#' @param annot An [ens_annot] object, or the path to a saved annotation.
#' @param chr_col,start_col,end_col Names of the columns of `x` holding
#'   chromosome, region start and region end.
#' @param add_cols Annotation columns to report. Default: `ensembl_gene_id`,
#'   `external_gene_name`, `custom_gene_name`, `description`,
#'   `external_synonym`, `gene_biotype`, `entrezgene_id`,
#'   `ensembl_transcript_id`. Gene coordinates (`chromosome_name`,
#'   `start_position`, `end_position`, `strand`) can be added here too.
#' @param min_overlap Minimum overlap in bp for a gene to count (default `1`,
#'   i.e. any overlap).
#' @param sep Separator between genes in the output (default `"; "`).
#' @param inner_sep Separator between multiple values belonging to one gene
#'   (default `","`).
#' @param na_string `NULL` (default): leave missing values out. A string such
#'   as `"NA"`: write that string for missing values.
#' @param gene_args Named list of extra arguments passed to
#'   [ens_annotate_ids()] when building the `genes` table, e.g.
#'   `list(add_cols = c("custom_gene_name", "start_position"),
#'   strip_version = FALSE, out_file = "genes.tsv")`. `x`, `annot` and `id_col`
#'   are set internally and cannot be overridden. By default the `genes` table
#'   has the [ens_annotate_ids()] default columns: `external_gene_name`,
#'   `custom_gene_name`, `description`, `external_synonym`, `gene_biotype`,
#'   `entrezgene_id`, `ensembl_transcript_id`, `chromosome_name`,
#'   `start_position`, `end_position`, `strand`.
#' @param region_rows If `TRUE` (default), the `genes` table gets a first column
#'   `region_row` with the row number(s) of `x` (joined with `sep`) in which
#'   each gene was found. `FALSE` omits it.
#' @param out_file Optional path; if given, the `regions` table is also written
#'   there as a tab-delimited file. (To write the `genes` table, use
#'   `gene_args = list(out_file = ...)`.)
#'
#' @return A list with two [data.table::data.table]s:
#' \describe{
#'   \item{`regions`}{all columns of `x`, then `n_genes` (number of overlapping
#'     genes), then the columns requested in `add_cols`.}
#'   \item{`genes`}{`region_row` (unless `region_rows = FALSE`),
#'     `ensembl_gene_id`, then the [ens_annotate_ids()] columns, one row per
#'     distinct overlapping gene, ordered by first appearance (input row, then
#'     genomic position). Zero rows if no region overlaps a gene.}
#' }
#' @seealso [ens_annotate_loci()], [ens_annotate_ids()], [ens_annot]
#' @examples
#' ann <- ens_example_annotation()
#' regions <- data.frame(
#'   chr   = c("1",  "chr1",  "3"),
#'   start = c(4000, 15000,   1),
#'   end   = c(7000, 25000, 100)
#' )
#'
#' # default: a list with the per-region table and the gene table
#' res <- ens_annotate_regions(regions, ann, "chr", "start", "end")
#' res$regions   # input rows + collapsed gene information
#' res$genes     # one row per distinct overlapping gene (with region_row)
#'
#' # just names and biotypes in the regions table
#' ens_annotate_regions(regions, ann, "chr", "start", "end",
#'                      add_cols = c("ensembl_gene_id", "custom_gene_name",
#'                                   "gene_biotype"))$regions
#'
#' # choose the columns of the genes table (passed to ens_annotate_ids)
#' ens_annotate_regions(regions, ann, "chr", "start", "end",
#'                      gene_args = list(add_cols = c("custom_gene_name",
#'                                                    "gene_biotype",
#'                                                    "start_position",
#'                                                    "end_position")))$genes
#'
#' # genes table without the region_row link column, written to a file
#' ens_annotate_regions(regions, ann, "chr", "start", "end",
#'                      region_rows = FALSE,
#'                      gene_args = list(out_file = tempfile(fileext = ".tsv")))$genes
#'
#' # require at least 1 kb overlap, keep placeholders so columns stay aligned
#' ens_annotate_regions(regions, ann, "chr", "start", "end",
#'                      add_cols = c("custom_gene_name", "external_synonym"),
#'                      min_overlap = 1000, na_string = "NA")$regions
#'
#' # include gene coordinates in the regions table, change separators
#' ens_annotate_regions(regions, ann, "chr", "start", "end",
#'                      add_cols = c("custom_gene_name", "start_position",
#'                                   "end_position"),
#'                      sep = " | ")$regions
#' @export
ens_annotate_regions <- function(x, annot, chr_col, start_col, end_col,
                                 add_cols = .ENS_ATTR_COLS, min_overlap = 1,
                                 sep = "; ", inner_sep = ",", na_string = NULL,
                                 gene_args = list(), region_rows = TRUE,
                                 out_file = NULL) {
  x <- .as_input(x); annot <- .as_annot(annot); add_cols <- .check_cols(annot, add_cols)
  q  <- .query_table(x, chr_col, start_col, end_col)
  ov <- .overlap_pairs(q, annot$genes, pad = 0)
  gs <- as.numeric(annot$genes$start_position[ov$gi]); ge <- as.numeric(annot$genes$end_position[ov$gi])
  olen <- pmin(q$e[ov$qi], ge) - pmax(q$s[ov$qi], gs) + 1
  ov <- ov[olen >= min_overlap]
  gs <- as.numeric(annot$genes$start_position[ov$gi])
  ov <- ov[order(ov$qi, gs)]

  n <- nrow(x)
  vals <- .gene_vals(annot$genes, ov$gi, add_cols, inner_sep, na_string)
  res  <- .collapse_by_row(n, ov$qi, vals, sep)
  n_genes <- tabulate(ov$qi, nbins = n)
  n_genes[!q$valid] <- NA_integer_
  res <- cbind(data.table(n_genes = n_genes), res)
  regions_df <- .finish(x, res, out_file)

  # ---- gene table: one row per distinct gene, built with ens_annotate_ids() ----
  ugi <- unique(ov$gi)                                  # order of first appearance
  gx  <- data.table(ensembl_gene_id = annot$genes$ensembl_gene_id[ugi])
  if (isTRUE(region_rows)) {
    rr <- ov[, .(region_row = paste(unique(qi), collapse = sep)), by = gi]
    gx <- cbind(data.table(region_row = as.character(rr$region_row[match(ugi, rr$gi)])), gx)
  }
  gene_args <- gene_args[setdiff(names(gene_args), c("x", "annot", "id_col"))]
  genes_df  <- suppressMessages(
    do.call(ens_annotate_ids,
            c(list(x = gx, annot = annot, id_col = "ensembl_gene_id"), gene_args)))

  list(regions = regions_df, genes = genes_df)
}

#' Annotate loci with overlapping or closest genes
#'
#' @description
#' Mode 3 of the annotation functions: for every locus (row) of your input,
#' reports the genes that overlap it; if none overlap, the closest gene(s)
#' within `max_dist` base pairs. The output has exactly one row per input row,
#' in the input order.
#'
#' @details
#' **Loci.** A locus is a single position (`pos_col` only) or a small interval
#' (`pos_col` = start and `end_col` = end). Coordinates are 1-based, as in
#' Ensembl. Chromosome names are normalised (`chr1` = `1`, `chrM` = `MT`).
#' Rows with missing or non-numeric coordinates get `n_genes = NA`.
#'
#' **Selection rule.**
#' \enumerate{
#'   \item If at least one gene overlaps the locus, **all** overlapping genes
#'     are reported (distance 0) and nothing else.
#'   \item Otherwise, genes lying within `max_dist` bp are considered. With
#'     `nonoverlap = "closest"` only the nearest one is reported (all genes
#'     are reported if several are equally near); with `"all_within"` every
#'     gene within `max_dist` is reported, nearest first.
#'   \item If no gene is within `max_dist`, `n_genes = 0` and the annotation
#'     columns are empty.
#' }
#' With the default `max_dist = 0` only overlapping genes are ever reported.
#'
#' **Distance and relation.** `gene_distance` is the number of bases from the
#' locus to the nearest base of the gene: `0` means overlap, `1` means the
#' locus is directly adjacent to the gene. `gene_relation` says where the
#' locus lies relative to the gene, taking the gene's strand into account:
#' `"overlap"`, `"upstream"` (the locus is before the gene's start in the
#' direction of transcription) or `"downstream"` (after its end).
#'
#' **Output layout.** Same conventions as [ens_annotate_regions()]: several
#' genes per row are joined with `sep`; multi-valued fields within one gene use
#' `inner_sep`; missing values are left out unless `na_string` is set;
#' `gene_relation` and `gene_distance` are listed in the same order as
#' `ensembl_gene_id`.
#'
#' @param x The input table: a data.frame / data.table, or the path to a
#'   delimited file.
#' @param annot An [ens_annot] object, or the path to a saved annotation.
#' @param chr_col Name of the chromosome column of `x`.
#' @param pos_col Name of the position column (or the interval start if
#'   `end_col` is given).
#' @param end_col Optional name of an interval end column. `NULL` (default)
#'   treats every row as a single position.
#' @param max_dist Maximum distance in bp to search when no gene overlaps the
#'   locus. Default `0` (overlap only).
#' @param nonoverlap `"closest"` (default) reports the nearest gene(s) when
#'   nothing overlaps; `"all_within"` reports all genes within `max_dist`.
#' @param add_cols Annotation columns to report; same default and options as
#'   in [ens_annotate_regions()].
#' @param sep Separator between genes in the output (default `"; "`).
#' @param inner_sep Separator between multiple values belonging to one gene
#'   (default `","`).
#' @param na_string `NULL` (default): leave missing values out. A string such
#'   as `"NA"`: write that string for missing values.
#' @param out_file Optional path; if given, the result is also written there as
#'   a tab-delimited file.
#'
#' @return A [data.table::data.table] with all columns of `x`, then `n_genes`,
#'   then the requested annotation columns, then `gene_relation` and
#'   `gene_distance`.
#' @seealso [ens_annotate_regions()], [ens_annotate_ids()], [ens_annot]
#' @examples
#' ann <- ens_example_annotation()
#' loci <- data.frame(
#'   chr = c("1",  "chr1", "2",    "X"),
#'   pos = c(2500, 5500,   10000,  50)
#' )
#'
#' # overlapping genes only (max_dist = 0)
#' ens_annotate_loci(loci, ann, "chr", "pos")
#'
#' # also report the closest gene within 1 kb when nothing overlaps
#' ens_annotate_loci(loci, ann, "chr", "pos", max_dist = 1000,
#'                   add_cols = c("ensembl_gene_id", "custom_gene_name"))
#'
#' # every gene within 10 kb (nearest first), instead of only the closest
#' ens_annotate_loci(loci, ann, "chr", "pos", max_dist = 10000,
#'                   nonoverlap = "all_within",
#'                   add_cols = c("custom_gene_name", "gene_biotype"))
#'
#' # loci that are intervals
#' loci2 <- data.frame(chr = "1", from = 4900, to = 6100)
#' ens_annotate_loci(loci2, ann, "chr", "from", end_col = "to",
#'                   add_cols = "custom_gene_name")
#' @export
ens_annotate_loci <- function(x, annot, chr_col, pos_col, end_col = NULL, max_dist = 0,
                              nonoverlap = c("closest", "all_within"),
                              add_cols = .ENS_ATTR_COLS,
                              sep = "; ", inner_sep = ",", na_string = NULL,
                              out_file = NULL) {
  nonoverlap <- match.arg(nonoverlap)
  x <- .as_input(x); annot <- .as_annot(annot); add_cols <- .check_cols(annot, add_cols)
  q  <- .query_table(x, chr_col, pos_col, if (is.null(end_col)) pos_col else end_col)
  g  <- annot$genes

  ov <- .overlap_pairs(q, g, pad = max_dist)
  qs <- q$s[ov$qi]; qe <- q$e[ov$qi]
  gs <- as.numeric(g$start_position[ov$gi]); ge <- as.numeric(g$end_position[ov$gi])
  d  <- ifelse(ge < qs, qs - ge, ifelse(gs > qe, gs - qe, 0))
  plus <- is.na(g$strand[ov$gi]) | g$strand[ov$gi] != -1L
  rel  <- ifelse(d == 0, "overlap",
                 ifelse(ge < qs, ifelse(plus, "downstream", "upstream"),
                        ifelse(plus, "upstream", "downstream")))
  ov[, `:=`(dist = d, rel = rel, gs = gs)]
  ov <- ov[dist <= max_dist]
  ov[, mind := min(dist), by = qi]
  ov <- if (nonoverlap == "closest") ov[dist == mind] else ov[mind > 0 | dist == 0]
  setorder(ov, qi, dist, gs)

  n <- nrow(x)
  vals <- c(.gene_vals(g, ov$gi, add_cols, inner_sep, na_string),
            list(gene_relation = ov$rel, gene_distance = as.character(ov$dist)))
  res  <- .collapse_by_row(n, ov$qi, vals, sep)
  n_genes <- tabulate(ov$qi, nbins = n)
  n_genes[!q$valid] <- NA_integer_
  res <- cbind(data.table(n_genes = n_genes), res)
  .finish(x, res, out_file)
}

