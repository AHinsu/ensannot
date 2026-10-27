# ensannot

**Ensembl gene annotation in R, without `biomaRt`.**

`ensannot` builds a local, reusable gene annotation table straight from the
[Ensembl FTP site](https://ftp.ensembl.org/pub/) (the core database dump files,
read as plain tab-delimited text, so no SQL server is needed) and uses it to
annotate your own data: gene lists, genomic regions (peaks, CNV segments) and
single loci (GWAS SNPs, variants).

Build once per species / release / assembly, cache it, and reuse it offline.

- [Why ensannot?](#why-ensannot)
- [Installation](#installation)
- [Quick start (offline, 1 minute)](#quick-start-offline-1-minute)
- [Step 1: build an annotation](#step-1-build-an-annotation)
- [Step 2: annotate your data](#step-2-annotate-your-data)
- [Real-case examples](#real-case-examples)
- [Output columns](#output-columns)
- [Function reference](#function-reference)
- [Notes and good practice](#notes-and-good-practice)
- [Troubleshooting](#troubleshooting)
- [Development](#development)
- [Data and license](#data-and-license)

---

## Why ensannot?

Ensembl changed how annotation is served, which broke or complicated
`biomaRt`-based workflows. The FTP site still has everything, and the files are
stable and versioned. `ensannot`:

- builds **one tidy table per Ensembl release and assembly**, so your
  annotation is reproducible and pinned to a version;
- caches it locally, so later sessions load in seconds and work **offline**;
- uses **biomaRt attribute names** (`ensembl_gene_id`, `external_gene_name`,
  `chromosome_name`, ...), so existing code needs few changes;
- covers three annotation tasks with one consistent interface: **IDs**,
  **regions** and **loci**;
- gives every gene a never-empty `custom_gene_name` (symbol, else synonym, else
  Ensembl ID), handy for plots and tables.

---

## Installation

### Requirements

- R >= 4.0.0
- Packages: [`data.table`](https://CRAN.R-project.org/package=data.table)
  (installed automatically)
- Recommended: [`R.utils`](https://CRAN.R-project.org/package=R.utils), which
  makes reading the large `.gz` files much faster
- Internet access for the first build of each annotation (not needed afterwards)

### From GitHub

```r
# install.packages("remotes")
remotes::install_github("ahinsu/ensannot")

# optional but recommended
install.packages("R.utils")
```

### From a local clone or source tarball

```r
# in a terminal:  git clone https://github.com/ahinsu/ensannot.git
install.packages("path/to/ensannot", repos = NULL, type = "source")
# or, with devtools:
devtools::install("path/to/ensannot")
```

Then load it:

```r
library(ensannot)
```

> **Windows note:** without `R.utils`, the `.gz` fallback calls `gzip` on the
> command line, which Windows may not have. Install `R.utils` to avoid this.

---

## Quick start (offline, 1 minute)

The package ships a tiny **synthetic** annotation (genes `TOYG001`-`TOYG005`,
invented coordinates), so you can try every function without downloading
anything:

```r
library(ensannot)

ann <- ens_example_annotation()
ann                    # summary
ann$genes              # the annotation table (a data.table)
```

Real data work the same way: just replace `ann` with a build from
`ens_build_annotation()` (next section).

---

## Step 1: build an annotation

```r
# Human, Ensembl release 115
ann <- ens_build_annotation("homo_sapiens", release = 115)

# Later sessions: load from the cache (fast, offline)
ann <- ens_load_annotation("homo_sapiens", release = 115)

# What is cached?
ens_list_cached()
```

The first build downloads the needed Ensembl dump files and parses them. For
large genomes this can take a while and needs a few GB of RAM and disk space;
the raw downloads are deleted after the build unless `keep_raw = TRUE`. The
result is cached as `<db_name>.rds` (for R) and `<db_name>.tsv.gz` (for
anything else).

### Choosing species, assembly and release

Ensembl names its core databases:

```
<genus_species>[_<assembly>]_core_<release>_<number>
```

| Folder on the FTP site | How to request it |
|---|---|
| `homo_sapiens_core_115_38` | `ens_build_annotation("homo_sapiens", release = 115)` |
| `capra_hircus_blackbengal_core_115_1` | `ens_build_annotation("capra_hircus", release = 115, assembly = "blackbengal")` |

- With **no** `assembly`, `ensannot` looks for `<genus_species>_core_...`.
- With `assembly = "x"`, it looks for `<genus_species>_x_core_...`. This makes
  sure you get exactly the assembly you asked for.
- To find assembly names, browse
  `https://ftp.ensembl.org/pub/release-<N>/mysql/` (or `current_mysql/`) and
  take the text between the species and `_core`. If a lookup fails, the error
  message lists what exists for that species.
- `release = "current"` uses the latest release.
- Prefer to give the exact folder? Use `ftp_url`:

```r
ann <- ens_build_annotation(
  ftp_url = "https://ftp.ensembl.org/pub/release-115/mysql/homo_sapiens_core_115_38/"
)
```

Useful options:

| Argument | Default | Meaning |
|---|---|---|
| `primary_only` | `TRUE` | Drop genes on patches/haplotypes and LRG records |
| `clean_description` | `FALSE` | Strip the trailing `[Source:...]` from `description` |
| `force` | `FALSE` | Rebuild even if a cached build exists |
| `keep_raw` | `FALSE` | Keep the downloaded dump files |
| `cache_dir` | per-user dir | Where builds are stored (see `ens_cache_dir()`) |

Set a cache location for the whole session (or in `.Rprofile`):

```r
options(ensannot.cache_dir = "~/ensembl_cache")
```

---

## Step 2: annotate your data

All three functions take a data frame **or a file path** (read with
`data.table::fread`), and return one result row per input row.

| Task | Function |
|---|---|
| Details for Ensembl gene IDs | `ens_annotate_ids()` |
| All genes overlapping each region | `ens_annotate_regions()` |
| Overlapping or closest genes for each locus | `ens_annotate_loci()` |

Using the toy annotation:

```r
ann <- ens_example_annotation()
```

### Gene IDs

```r
x <- data.frame(gene = c("TOYG001", "TOYG003.2", "NOT_A_GENE"),
                logFC = c(1.5, -2.1, 0.3))

ens_annotate_ids(x, ann, id_col = "gene")
```

Version suffixes (`.2`) are ignored, and unmatched IDs give `NA`.

### Regions

```r
regions <- data.frame(chr   = c("chr1", "1", "X"),
                      start = c(4000, 15000, 1),
                      end   = c(7000, 25000, 100))

res <- ens_annotate_regions(regions, ann, "chr", "start", "end",
                            add_cols = c("ensembl_gene_id", "custom_gene_name"),
                            gene_args = list(add_cols = c("custom_gene_name",
                                                          "gene_biotype")))
res$regions
res$genes
```

`ens_annotate_regions()` returns **a list of two tables**:

`res$regions`: your rows with all overlapping genes collapsed into each row

| chr | start | end | n_genes | ensembl_gene_id | custom_gene_name |
|---|---|---|---|---|---|
| chr1 | 4000 | 7000 | 2 | TOYG001; TOYG002 | ALPHA1; B2 |
| 1 | 15000 | 25000 | 1 | TOYG003 | TOYG003 |
| X | 1 | 100 | 1 | TOYG005 | DELTA4 |

`res$genes`: one row per distinct gene found in any region (long format);
`region_row` says which input row(s) it came from

| region_row | ensembl_gene_id | custom_gene_name | gene_biotype |
|---|---|---|---|
| 1 | TOYG001 | ALPHA1 | protein_coding |
| 1 | TOYG002 | B2 | protein_coding |
| 2 | TOYG003 | TOYG003 | lncRNA |
| 3 | TOYG005 | DELTA4 | miRNA |

`gene_args` is passed on to `ens_annotate_ids()` when building the `genes`
table, so you can choose its columns (`add_cols`), keep version suffixes
(`strip_version = FALSE`) or write it to a file (`out_file = "genes.tsv"`).

### Loci

```r
loci <- data.frame(chr = c("1", "chr1", "2", "X"),
                   pos = c(2500, 5500, 10000, 50))

ens_annotate_loci(loci, ann, "chr", "pos", max_dist = 1000,
                  add_cols = c("ensembl_gene_id", "custom_gene_name"))
```

| chr | pos | n_genes | ensembl_gene_id | custom_gene_name | gene_relation | gene_distance |
|---|---|---|---|---|---|---|
| 1 | 2500 | 1 | TOYG001 | ALPHA1 | overlap | 0 |
| chr1 | 5500 | 2 | TOYG001; TOYG002 | ALPHA1; B2 | downstream; downstream | 500; 500 |
| 2 | 10000 | 0 | | | | |
| X | 50 | 1 | TOYG005 | DELTA4 | downstream | 50 |

Rule: if any gene overlaps the locus, **all overlapping genes** are reported.
Otherwise the **closest** gene(s) within `max_dist` bp are reported (all of
them if tied; use `nonoverlap = "all_within"` for every gene within the
distance). `gene_relation` is relative to the gene's strand (`overlap`,
`upstream`, `downstream`) and `gene_distance` is in bp (`0` = overlap).

---

## Real-case examples

The examples below use human (GRCh38) and Ensembl release 115. Positions given
here are for illustration; always use the **same assembly and a matching
Ensembl release** as the coordinates in your data.

```r
library(ensannot)
ann <- ens_build_annotation("homo_sapiens", release = 115)   # or ens_load_annotation()
```

### 1) Add gene symbols and biotypes to differential expression results

Results from DESeq2/edgeR/limma often have versioned Ensembl IDs
(`ENSG00000141510.18`). Annotate them directly:

```r
res <- read.csv("deseq2_results.csv")        # has a column "gene_id"

res_annot <- ens_annotate_ids(
  res, ann, id_col = "gene_id",
  add_cols = c("custom_gene_name", "external_gene_name", "gene_biotype",
               "entrezgene_id", "chromosome_name"),
  out_file = "deseq2_results_annotated.tsv"
)

# keep only protein-coding genes, significant
subset(res_annot, gene_biotype == "protein_coding" & padj < 0.05)
```

A tiny runnable version:

```r
res <- data.frame(
  gene_id = c("ENSG00000139618.17", "ENSG00000141510.18",
              "ENSG00000146648.20", "ENSG00000136997.20"),   # BRCA2, TP53, EGFR, MYC
  log2FoldChange = c(1.2, -2.3, 0.8, 3.1)
)
ens_annotate_ids(res, ann, "gene_id",
                 add_cols = c("custom_gene_name", "gene_biotype", "entrezgene_id"))
```

`custom_gene_name` is always filled, so it is safe to use as a plot label even
for genes without an official symbol.

### 2) Annotate ChIP-seq / ATAC-seq peaks with the genes they overlap

```r
library(data.table)
peaks <- fread("peaks.narrowPeak", select = 1:3,
               col.names = c("chr", "start", "end"))
peaks[, start := start + 1L]     # BED is 0-based; ensannot uses 1-based coordinates

ann_peaks <- ens_annotate_regions(
  peaks, ann, "chr", "start", "end",
  add_cols = c("ensembl_gene_id", "custom_gene_name", "gene_biotype"),
  out_file = "peaks_annotated.tsv"
)

ann_peaks$regions                 # one row per peak
table(ann_peaks$genes$gene_biotype)   # which kinds of genes the peaks hit
```

Peaks that overlap no gene have `n_genes = 0` and empty annotation. To count
only substantial overlaps, add `min_overlap = 100` (bp).

### 3) Annotate copy-number segments and export the affected genes

```r
segs <- data.frame(chr = c("17", "8"),
                   start = c(7600000, 127700000),
                   end   = c(7700000, 127800000))

cnv <- ens_annotate_regions(
  segs, ann, "chr", "start", "end",
  gene_args = list(add_cols = c("custom_gene_name", "gene_biotype",
                                "chromosome_name", "start_position", "end_position"),
                   out_file = "cnv_genes.tsv")
)

cnv$regions    # compact: genes per segment
cnv$genes      # long: one row per gene, with region_row linking back to segs
```

### 4) Map GWAS lead SNPs to genes (overlap, otherwise nearest within 100 kb)

```r
snps <- data.frame(
  rsid = c("rs429358", "rs7903146"),
  chr  = c("19", "10"),
  pos  = c(44908684, 112998590)          # GRCh38 positions
)

ens_annotate_loci(snps, ann, "chr", "pos", max_dist = 100000,
                  add_cols = c("ensembl_gene_id", "custom_gene_name", "gene_biotype"))
```

Each SNP gets its overlapping gene(s), or the nearest gene(s) within 100 kb,
with `gene_relation` and `gene_distance`. Use `nonoverlap = "all_within"` to
list every gene in the window instead of only the closest.

### 5) Another species, with an assembly name in the Ensembl folder

Some Ensembl species have several assemblies. The folder name then contains the
assembly, e.g. `capra_hircus_blackbengal_core_115_1`:

```r
goat <- ens_build_annotation("capra_hircus", release = 115, assembly = "blackbengal")
goat
```

Browse `https://ftp.ensembl.org/pub/release-115/mysql/` to see the names that
exist for your species.

### 6) Share the table with someone who does not use R

```r
ens_list_cached()                       # shows the file paths
# <db_name>.tsv.gz is a plain tab-delimited table you can open anywhere

ann <- ens_load_annotation(file = "homo_sapiens_core_115_38.tsv.gz")   # reload later
```

### 7) Annotate straight from files

All annotate functions accept a file path:

```r
ens_annotate_loci("snps.tsv", ann, "chr", "pos", max_dist = 50000,
                  out_file = "snps_annotated.tsv")
```

---

## Output columns

The annotation table (`ann$genes`) has one row per current gene, with biomaRt
attribute names:

| Column | Content |
|---|---|
| `ensembl_gene_id` | Stable Ensembl gene ID (no version) |
| `chromosome_name` | Chromosome or scaffold name as in Ensembl (`1`, `X`, `MT`, ...) |
| `start_position`, `end_position` | 1-based, inclusive gene coordinates |
| `strand` | `1` or `-1` |
| `external_gene_name` | Gene symbol; `NA` if none |
| `custom_gene_name` | Symbol, else first synonym (alphabetical), else Ensembl gene ID. Never empty |
| `description` | Gene description (keeps `[Source:...]` unless `clean_description = TRUE`) |
| `external_synonym` | Synonyms, separated by `"; "` |
| `gene_biotype` | e.g. `protein_coding`, `lncRNA`, `miRNA` |
| `entrezgene_id` | Entrez Gene ID(s), separated by `"; "` |
| `ensembl_transcript_id` | Transcript IDs, separated by `"; "` |

Where the values come from:

- **Symbol:** the gene's display reference.
- **Synonyms:** all external references linked to the gene, so genes without an
  HGNC-style symbol can still get synonyms from other databases.
- **Entrez ID:** `EntrezGene` references linked to the gene, its transcripts or
  its translations.

In the region and locus outputs, several genes in one row are joined with
`"; "`, and multiple values belonging to one gene (synonyms, Entrez IDs,
transcripts) with `","`, so the two levels stay distinguishable. Missing values
are left out (not written as `NA`); use `na_string = "NA"` to write a
placeholder and keep columns aligned gene by gene.

---

## Function reference

| Function | Purpose |
|---|---|
| `ens_build_annotation()` | Download, parse and cache an annotation (reuses the cache) |
| `ens_load_annotation()` | Load a cached build, or a saved `.rds` / `.tsv(.gz)` |
| `ens_list_cached()` | Table of cached builds |
| `ens_resolve_db_url()` | Species + assembly + release -> Ensembl database folder URL |
| `ens_cache_dir()` | Current cache directory |
| `ens_example_annotation()` | Tiny synthetic annotation for examples and tests |
| `ens_annotate_ids()` | Details for Ensembl gene IDs |
| `ens_annotate_regions()` | All genes per region; returns `list(regions, genes)` |
| `ens_annotate_loci()` | Overlapping or closest genes per locus |

Every function has full documentation with examples: `?ens_annotate_loci`,
`?ens_build_annotation`, `?ens_annot` (the annotation object), and so on.

Key arguments of the annotate functions:

| Argument | Applies to | Meaning |
|---|---|---|
| `add_cols` | all | Which annotation columns to add |
| `min_overlap` | regions | Minimum overlap (bp) for a gene to count |
| `max_dist` | loci | Search distance (bp) when nothing overlaps |
| `nonoverlap` | loci | `"closest"` (default) or `"all_within"` |
| `end_col` | loci | Treat loci as small intervals instead of single positions |
| `gene_args` | regions | Extra arguments passed to `ens_annotate_ids()` for the `genes` table |
| `sep`, `inner_sep`, `na_string` | regions, loci | Output separators and missing-value handling |
| `out_file` | all | Also write the result as a tab-delimited file |

---

## Notes and good practice

- **Match assembly and release to your data.** Gene coordinates differ between
  assemblies (e.g. GRCh37 vs GRCh38). For human GRCh37, Ensembl keeps its
  databases in a separate part of the FTP site (under `pub/grch37/`); find the
  right folder there and pass it with `ftp_url`.
- **Coordinates are 1-based and inclusive.** BED files are 0-based; add 1 to the
  start for exactness (see example 2).
- **Chromosome names are normalised.** `chr1`, `Chr1` and `1` are the same;
  `chrM`/`M`/`MT` match `MT`. Scaffold names must agree with Ensembl's.
- **Gene-level only.** The table has genes with their transcript IDs, not
  exon-level coordinates.
- **Synonyms and Entrez IDs depend on the Ensembl cross-references** of the
  release and species; coverage varies, especially outside model organisms.
- **Patches and haplotypes** are dropped by default (`primary_only = TRUE`);
  unplaced scaffolds are kept.
- **Pin your release** (e.g. `release = 115`) in scripts for reproducibility.
  `release = "current"` changes over time.

---

## Troubleshooting

| Problem | What to do |
|---|---|
| `No core database named '...'` | Check the species spelling. If the species has several assemblies, set `assembly` (the text between species and `_core` in the FTP folder). The error lists what exists |
| `Several databases match` | Pass the exact folder with `ftp_url` |
| `Download failed` / timeouts | Retry (finished files are kept); check your network or proxy. Large files can need a long timeout; `ensannot` raises R's timeout itself |
| Build is slow or memory-hungry | Install `R.utils`; close other memory-heavy sessions; run the first build for large species on a machine with several GB of RAM |
| `built with an older version of ensannot.R` | The cached object predates the current column names; rebuild with `ens_build_annotation(..., force = TRUE)` |
| All IDs "not found" | Wrong species, assembly or release for your IDs. Retired genes are not in newer releases |
| Peaks/SNPs get no genes | Check the chromosome naming and that your coordinates match the annotation's assembly |
| Cache in the wrong place | Set `options(ensannot.cache_dir = "...")` or pass `cache_dir` |

---

## Development

Bug reports and suggestions are welcome via
[GitHub issues](https://github.com/ahinsu/ensannot/issues).

---

## Data and license

- Package code: MIT license (see `LICENSE`).
- Ensembl data are provided by EMBL-EBI under the
  [Ensembl data disclaimer](https://www.ensembl.org/info/about/legal/disclaimer.html);
  please cite Ensembl when you use annotations built with this package.
- `ensannot` is an independent tool and is not affiliated with or endorsed by
  Ensembl or EMBL-EBI.

## Citation

If use ensannot in your work, please cite the tools as: 
Ankit Hinsu. (2026). AHinsu/ensannot: ensannot v0.1.0 (Version v0.1.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.23688272
