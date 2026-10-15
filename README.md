# ensannot

`ensannot` builds reusable gene annotation tables from Ensembl FTP core database dumps (without `biomaRt`) and uses them to annotate your own datasets.

## What the code does

The package has two parts:

1. **Build and cache annotation**
   - `ens_build_annotation()` downloads/parses Ensembl core dump tables and creates one annotation table per species/release/assembly.
   - `ens_load_annotation()` reloads a cached or saved annotation object.
   - `ens_list_cached()` lists locally cached builds.
   - `ens_resolve_db_url()` resolves species/release/assembly to the Ensembl database folder URL.

2. **Annotate user data**
   - `ens_annotate_ids()` adds annotation fields for Ensembl gene IDs.
   - `ens_annotate_regions()` finds all genes overlapping each genomic region.
   - `ens_annotate_loci()` finds overlapping genes or nearest genes to loci/positions.

The annotation output follows biomaRt-style column names (`ensembl_gene_id`, `external_gene_name`, `description`, `gene_biotype`, etc.) plus `custom_gene_name`.

## Quick use cases / examples

```r
library(ensannot)
```

### 1) Work offline with built-in example annotation

```r
ann <- ens_example_annotation()
ann
head(ann$genes)
```

### 2) Annotate a gene-level results table by Ensembl ID

```r
x <- data.frame(
  gene  = c("TOYG001", "TOYG003.2", "NOT_A_GENE"),
  logFC = c(1.5, -2.1, 0.3)
)

id_annot <- ens_annotate_ids(x, ann, id_col = "gene")
id_annot
```

### 3) Annotate genomic regions with all overlapping genes

```r
regions <- data.frame(
  chr   = c("chr1", "1", "X"),
  start = c(4000, 15000, 1),
  end   = c(7000, 25000, 100)
)

region_annot <- ens_annotate_regions(regions, ann, "chr", "start", "end")
region_annot$regions
region_annot$genes
```

### 4) Annotate loci (position-level) with overlap/nearest-gene logic

```r
loci <- data.frame(
  chr = c("1", "chr1", "2", "X"),
  pos = c(2500, 5500, 10000, 50)
)

locus_annot <- ens_annotate_loci(loci, ann, "chr", "pos", max_dist = 1000)
locus_annot
```

### 5) Build and reuse a real Ensembl annotation

```r
# First run downloads and builds (can take time/resources for large species)
ann_hs <- ens_build_annotation("homo_sapiens", release = 115)

# Later runs are fast from cache
ann_hs2 <- ens_load_annotation("homo_sapiens", release = 115)

# See what is cached
ens_list_cached()
```
