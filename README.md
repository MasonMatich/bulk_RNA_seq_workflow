# Bulk RNA-seq Analysis Workflow

A reproducible, automated workflow for differential gene expression analysis of bulk RNA-seq data using DESeq2, with comprehensive gene ontology enrichment analysis and interactive HTML reports.

## Overview

This workflow performs:
- **Differential expression analysis** using DESeq2
- **Principal component analysis** (PCA) with batch effect visualization
- **Optional batch correction** using ComBat_seq
- **Gene ontology enrichment** (GO classification, over-representation, GSEA)
- **Reactome pathway analysis**
- **Interactive HTML reports** with downloadable tables and visualizations

Originally designed for *C. elegans* data but adaptable to other organisms.

---

## Features

✅ **Multi-sample support** - Analyze multiple datasets in one run
✅ **Batch effect correction** - Optional ComBat_seq integration
✅ **Comprehensive validation** - Input validation with clear error messages
✅ **Reproducible** - Session info and package versions recorded
✅ **Interactive outputs** - HTML tables with export to CSV/Excel
✅ **Publication-ready figures** - High-resolution volcano plots, MA plots, PCA
✅ **Robust path management** - Works from any directory using `here` package

---

## Requirements

### R Version
- R ≥ 4.0.0 recommended

### Required Packages

**CRAN Packages:**
```r
c("bookdown", "dplyr", "DT", "knitr", "kableExtra", "ashr", "htmltools",
  "ggplot2", "PCAtools", "knitrBootstrap", "stringr", "tibble", "readr",
  "patchwork", "tidyverse", "RColorBrewer", "here")
```

**Bioconductor Packages:**
```r
c("ReactomePA", "clusterProfiler", "DESeq2", "GOSemSim", "sva", "vsn",
  "org.Ce.eg.db", "EnhancedVolcano", "rrvgo", "PCAtools")
```

**Note:** For organisms other than *C. elegans*, replace `org.Ce.eg.db` with the appropriate organism database (e.g., `org.Hs.eg.db` for human, `org.Mm.eg.db` for mouse).

### Installation

All packages are automatically checked and installed when the workflow runs. Alternatively, install manually:

```r
# CRAN packages
install.packages(c("bookdown", "dplyr", "DT", "knitr", "here", ...))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("DESeq2", "clusterProfiler", "org.Ce.eg.db", ...))
```

---

## Directory Structure

```
bulk_RNA_seq_workflow/
├── index.Rmd                    # Main analysis script
├── _bookdown.yml               # Bookdown configuration
├── _output.yml                 # Output format settings
├── style.css                   # Custom CSS for reports
├── scripts/
│   ├── analysis.Rmd           # Child document for individual samples
│   ├── sampleInfo.csv         # Parameter file (REQUIRED)
│   ├── chen_phillips_metadata.tsv  # Example metadata
│   ├── compileReport.R        # Report generation functions
│   ├── translateBioIDs.R      # Gene ID conversion
│   ├── callClusterProfilerFunc.R  # GO/pathway analysis
│   ├── runDESeq2Analysis.R    # DESeq2 wrapper
│   ├── visualizePCsByCondition.R  # PCA plotting
│   ├── knitDataTable.R        # Interactive table rendering
│   ├── emptyTableMessage.R    # Empty result handling
│   └── reduceGORedundancy.R   # GO term simplification
├── raw_data/                   # Input count matrices and metadata
├── figures/                    # Output plots (auto-created)
└── results/                    # Output tables and RDS files (auto-created)
```

---

## Quick Start

### 1. Prepare Your Data

#### Count Matrix
- Format: TSV file with columns: `gene_id`, `gene_name`, sample columns
- Example:
  ```
  gene_id         gene_name    sample1  sample2  sample3
  WBGene00000001  aap-1        150      200      175
  WBGene00000002  aat-1        1200     1150     1300
  ```

#### Metadata File
- Format: TSV with sample conditions
- Must match count matrix column names
- Example (`chen_phillips_metadata.tsv`):
  ```
  strain    replicate
  wt        rep1
  wt        rep2
  hrde1     rep1
  hrde1     rep2
  ```

### 2. Configure Parameters

Edit `scripts/sampleInfo.csv`:

```csv
sample_title,sample_name,salmon_merged_gene_counts_file_path,metadata_file_path,ref_level_cond,ref_level_value,design,batch_correct,batch_cond,l2fc_filter,padj_filter
Chen & Phillips,chen_phillips,../raw_data/chen_phillips/salmon.merged.gene_counts.tsv,scripts/chen_phillips_metadata.tsv,strain,wt,~ replicate + strain,TRUE,replicate,0.8,0.1
```

**Parameter Descriptions:**

| Parameter | Description | Example |
|-----------|-------------|---------|
| `sample_title` | Display name for report | `"Chen & Phillips"` |
| `sample_name` | Internal identifier (no spaces) | `chen_phillips` |
| `salmon_merged_gene_counts_file_path` | Path to count matrix | `../raw_data/counts.tsv` |
| `metadata_file_path` | Path to sample metadata | `scripts/metadata.tsv` |
| `ref_level_cond` | Column name for experimental condition | `strain` |
| `ref_level_value` | Control/reference level | `wt` |
| `design` | DESeq2 design formula | `~ replicate + strain` |
| `batch_correct` | Enable batch correction? | `TRUE` or `FALSE` |
| `batch_cond` | Batch variable to regress out | `replicate` |
| `l2fc_filter` | Log2 fold-change threshold | `0.8` |
| `padj_filter` | Adjusted p-value threshold | `0.1` |

### 3. Run Analysis

**Option A: RStudio**
1. Open `bulk_RNA_seq_workflow.Rproj`
2. Open `index.Rmd`
3. Click "Knit" → "Knit to bookdown::html_document2"

**Option B: Command Line**
```r
# In R console
rmarkdown::render("index.Rmd", params = list(params_file = "scripts/sampleInfo.csv"))
```

**Option C: With Custom Parameters**
```r
bookdown::render_book("index.Rmd",
                      params = list(
                        dynamic_title = "My RNA-seq Analysis",
                        params_file = "scripts/sampleInfo.csv"
                      ))
```

### 4. View Results

Output: `bulk_RNA_seq_workflow.html` (interactive HTML report)

---

## Output Files

### Figures (`figures/`)
- `{sample_name}_PCA.pdf` - PCA plots
- `{sample_name}_{contrast}_volcano.png` - Volcano plots

### Results (`results/`)
- `{sample_name}_dds.rds` - DESeq2 dataset object
- `{sample_name}_rld.rds` - Regularized log-transformed counts
- `{sample_name}_{contrast}_DEGs.csv` - Differentially expressed genes
- `{sample_name}_{contrast}_GO_*.csv` - Gene ontology results
- `{sample_name}_{contrast}_Reactome_*.csv` - Pathway analysis results

**File naming convention:** `{sample_name}` from sampleInfo.csv, `{contrast}` format: `{condition}_{level1}_vs_{level2}`

---

## Understanding the Workflow

### Analysis Steps

1. **Input Validation** - Checks file paths, parameters, data structure
2. **Count Matrix Loading** - Reads and validates count data
3. **Metadata Alignment** - Matches samples to experimental conditions
4. **DESeq2 Analysis** - Differential expression with optional batch correction
5. **PCA Visualization** - Before/after batch correction
6. **Contrast Extraction** - All experimental vs. control comparisons
7. **Gene ID Translation** - WormBase → Entrez ID conversion
8. **Enrichment Analysis** - GO classification, over-representation, GSEA
9. **Report Generation** - Interactive HTML with figures and tables

### Design Formula

The `design` parameter follows DESeq2 formula syntax:

```r
~ batch_variable + condition_of_interest
```

**Examples:**
- Simple comparison: `~ treatment`
- With batch effects: `~ replicate + treatment`
- Multiple factors: `~ sex + treatment`
- Interaction term: `~ sex + treatment + sex:treatment`

**Important:** The last term is your primary variable of interest!

### Batch Correction

When `batch_correct = TRUE`:
- **Method:** ComBat_seq from `sva` package
- **What it does:** Removes technical variation while preserving biological signal
- **When to use:** Technical replicates, sequencing batches, different labs
- **Preserves:** All conditions NOT specified in `batch_cond`

---

## Customization

### Changing Organism

1. Replace organism database:
   ```r
   # In index.Rmd, line ~31
   bioconductor <- c(..., "org.Hs.eg.db", ...)  # Human
   # OR
   bioconductor <- c(..., "org.Mm.eg.db", ...)  # Mouse
   ```

2. Update semantic data:
   ```r
   # In index.Rmd, line ~83
   semantic_data <- GOSemSim::godata(org.Hs.eg.db, ont="BP", keytype = "ENTREZID")
   ```

3. Update gene ID translation:
   ```r
   # In translateBioIDs.R, modify bioID check and bitr() calls
   ```

### Adjusting Thresholds

Edit `scripts/sampleInfo.csv`:
- `l2fc_filter`: Minimum absolute log2 fold-change (default: 0.8 = 1.74-fold)
- `padj_filter`: Maximum adjusted p-value (default: 0.1)

Or modify globally in `analysis.Rmd` line ~258-261.

### Figure Sizes

**Volcano plots:** `scripts/compileReport.R`, line 46
```r
ggsave(volcano_path, volcano, height = 12, width = 10, units = "in", dpi = 300)
```

**Image display size:** `scripts/compileReport.R`, line 47
```r
cat(sprintf('<img src="%s" width="80%%"/>\n\n', volcano_rel_path))
```

---

## Troubleshooting

### Common Errors

**Error: Count matrix file not found**
```
❌ Count matrix file not found: ../raw_data/counts.tsv
```
**Solution:** Check path in `sampleInfo.csv` is correct relative to project root

**Error: Metadata dimensions don't match**
```
❌ number of rows entered as metadata doesn't match number of columns in count matrix
```
**Solution:** Ensure metadata has exactly one row per sample (column) in count matrix

**Error: Reference level not found**
```
❌ Reference level value 'control' not found in metadata column 'treatment'
```
**Solution:** Check spelling of `ref_level_value` matches values in metadata file

**Error: Invalid design formula**
```
❌ Invalid design formula for sample 'my_sample': ~treatment
```
**Solution:** Ensure formula syntax is correct and variables exist in metadata

### Package Installation Issues

If BiocManager fails with "type='source'" error:
```r
# Try binary installation instead
BiocManager::install(c("DESeq2", "clusterProfiler"), type="binary")
```

### Memory Issues

For very large datasets (>50,000 genes, >100 samples):
```r
# Increase memory (before running)
options(java.parameters = "-Xmx8g")  # 8GB RAM
```

---

## Citation

If you use this workflow, please cite:

**DESeq2:**
Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* 15:550.

**clusterProfiler:**
Yu, G., Wang, L.G., Han, Y., He, Q.Y. (2012) clusterProfiler: an R package for comparing biological themes among gene clusters. *OMICS* 16(5):284-287.

**Method reference:**
Tutorial adapted from NYU Gencore: https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/

---

## FAQ

**Q: Can I analyze multiple experiments at once?**
A: Yes! Add additional rows to `sampleInfo.csv`, each representing a separate experiment.

**Q: What if I don't have batch effects?**
A: Set `batch_correct = FALSE` and leave `batch_cond` empty.

**Q: How do I export tables from the HTML report?**
A: Click the "CSV" or "Excel" buttons above each interactive table.

**Q: Can I re-run analysis without re-computing everything?**
A: Load saved RDS files:
```r
dds <- readRDS("results/chen_phillips_dds.rds")
rld <- readRDS("results/chen_phillips_rld.rds")
```

**Q: How do I change the report title?**
A: Modify `dynamic_title` parameter when rendering:
```r
rmarkdown::render("index.Rmd",
                  params = list(dynamic_title = "My Custom Title"))
```

**Q: What's the difference between GO classification, over-representation, and GSEA?**
- **Classification:** Groups genes into broad GO categories
- **Over-representation:** Tests if specific GO terms are enriched beyond chance
- **GSEA:** Ranks all genes and finds coordinated changes in pathways

---

## Contributing

Issues and pull requests welcome! For major changes, please open an issue first to discuss.

---

## License

This project is provided as-is for academic research use.

---

## Contact

**Author:** Mason Matich (masonmatich@gmail.com)

For questions or issues, please:
1. Check the troubleshooting section above
2. Review error messages carefully (they include suggested fixes!)
3. Open a GitHub issue with your error message and `sessionInfo()` output

---

## Changelog

### Version 2.0 (Current)
- ✅ Added comprehensive parameter validation
- ✅ Implemented robust path management with `here` package
- ✅ Added file existence checks with clear error messages
- ✅ Auto-create output directories
- ✅ Include session info in reports for reproducibility
- ✅ Refactored DESeq2 code to reduce duplication
- ✅ Improved error messages with actionable suggestions

### Version 1.0
- Initial release with basic DESeq2 workflow
- GO and Reactome enrichment analysis
- Interactive HTML reports
