# SETUP

How to configure the R, Python and database environments needed to run the code in
`analysis/`, `figures/` and `tables/`.

> **Scope.** This document covers **only** `analysis/`, `figures/` and `tables/` —
> the code that produced the figures, tables and downstream analyses of the
> manuscript. It does **not** cover `legacy_codes/` or `preprocessing_pipeline/`.
> Those two directories were created separately, through the standardized GLASS
> consortium pipeline protocols previously described in the publication linked in
> [`README.md`](README.md), and they carry their own conda environment
> specifications. See *Directories this document does not cover* at the end.

---

## 1. What you need before anything else

| Requirement | Notes |
|---|---|
| **R ≥ 4.2** | The floor imposed by the Bioconductor release these packages come from. R 4.3 is what the code was last run under. |
| **Python ≥ 3.9** | Only `figures/python/` and the two pipeline drivers in `analysis/python/` use it. |
| **PostgreSQL client libraries** | `libpq` headers, for building `RPostgres`. On Linux: `libpq-dev` / `postgresql-devel`. |
| **unixODBC** | Only if you use the ODBC route in section 3. `unixodbc-dev` / `unixODBC-devel` plus `psqlodbc`. |
| **Access to the GLASS database** | Or the equivalent flat files from Synapse — see section 3. |
| **Optional: conda/mamba** | The two pipeline drivers in `analysis/python/` expect named conda environments. Everything else runs from a plain R/Python install. |

A compiler toolchain is required: several packages build from source
(`RPostgres`, `svglite`, `data.table` on some platforms).

---

## 2. R environment

### 2.1 CRAN packages

Everything the three directories need from CRAN, in one call:

```r
install.packages(c(
  # data handling and the tidyverse core
  "tidyverse", "data.table", "vroom", "fs", "magrittr", "broom", "tibble",
  "dplyr", "tidyr", "purrr", "stringr", "forcats",

  # database
  "DBI", "RPostgres", "odbc", "pool",

  # plotting
  "ggplot2", "ggpubr", "ggbeeswarm", "ggridges", "ggalluvial", "ggh4x",
  "ggrepel", "patchwork", "cowplot", "egg", "gridExtra", "scales",
  "RColorBrewer", "colorspace", "circlize", "corrplot", "svglite",

  # statistics and survival
  "survival", "survminer", "rstatix", "EnvStats",

  # tables
  "gt", "gtsummary",

  # command-line parsing used by the snakemake-facing scripts
  "optparse"
))
```

`httpgd` appears in a few figure scripts as an interactive graphics device. It is
**not required** — the scripts run headless without it. Install it only if you
want live plot preview:

```r
install.packages("httpgd")
```

### 2.2 Bioconductor packages

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "GenomicRanges",
  "VariantAnnotation",
  "ComplexHeatmap",
  "BSgenome.Hsapiens.UCSC.hg19",   # used by the kataegis / Palimpsest stage
  "BSgenome.Hsapiens.UCSC.hg38"    # used by the CN and chromothripsis figures
))
```

Both BSgenome packages are large (roughly 1 GB installed each). Install only the build you need:
the **hg19/GRCh37** build is used by `analysis/python/*.sh`, the **hg38** build by
`figures/R/`. See section 5 on why both appear.

### 2.3 Packages not on CRAN or Bioconductor

These five are installed from source repositories and are each required by a
specific script. Install `remotes` first.

```r
install.packages("remotes")

# dN/dS selection analysis  -> figures/R/GLASS-I_dNdScv.R
remotes::install_github("im3sanger/dndscv")

# chromothripsis calling    -> analysis/R/GLASS-I_ShatterSeek.R
remotes::install_github("parklab/ShatterSeek")

# Sankey QC diagram         -> figures/R/GLASS-I_Sankey_QC.R
remotes::install_github("davidsjoberg/ggsankey")

# allele-specific copy number -> analysis/R/ascat_analysis.R, analysis/snakemake/ascat.smk
remotes::install_github("VanLoo-lab/ascat", subdir = "ASCAT")

# mutual exclusivity        -> figures/R/GLASS-I_sSNV_Mutual_Exclusivity.R
# DISCOVER is distributed from the authors' own CRAN-style repository, not GitHub.
# Check https://github.com/NKI-CCB/DISCOVER for the current repository URL if this
# one has moved.
install.packages("discover", repos = "http://ccb.nki.nl/software/discover/repos/r")
```

**Palimpsest** is needed only by the signature deconvolution stage of
`analysis/python/mutational_signature_pipeline.sh`, which runs R inside a conda
environment (section 4.2), not from your main R library:

```r
remotes::install_github("FunGeST/Palimpsest")
```

### 2.4 Verifying the R install

```r
pkgs <- c("tidyverse","data.table","vroom","fs","magrittr","broom",
          "DBI","RPostgres","odbc","pool",
          "ggplot2","ggpubr","ggbeeswarm","ggridges","ggalluvial","ggh4x",
          "ggrepel","patchwork","cowplot","egg","gridExtra","scales",
          "RColorBrewer","colorspace","circlize","corrplot","svglite",
          "survival","survminer","rstatix","EnvStats","gt","gtsummary","optparse",
          "GenomicRanges","VariantAnnotation","ComplexHeatmap",
          "dndscv","ShatterSeek","ggsankey","ASCAT","discover")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) cat("MISSING:", paste(missing, collapse = ", "), "\n") else
  cat("All packages present.\n")
```

### 2.5 Which packages each directory needs

If you only want to run part of the repository:

| Directory | R packages beyond the tidyverse/DBI core |
|---|---|
| `tables/R/` | `gt`, `gtsummary`, `survival`, `survminer`, `RPostgres` |
| `figures/R/` | the full plotting set, plus `ComplexHeatmap`, `GenomicRanges`, `BSgenome.Hsapiens.UCSC.hg38`, `circlize`, `dndscv`, `discover`, `ggsankey`, `survival`, `survminer`, `rstatix`, `EnvStats`, `corrplot` |
| `analysis/R/` | `ShatterSeek`, `VariantAnnotation`, `GenomicRanges`, `ASCAT`, `optparse`, `data.table`, `vroom`, `fs`, `pool`, `odbc` (plus base `parallel`, which ships with R) |

---

## 3. Database connection

Every R script in `figures/R/` and `tables/R/` opens the GLASS PostgreSQL database
at the top. The underlying data is public. As noted in `README.md`, the GLASS release is on
Synapse — `syn17038081` (tables) and `syn26465623` (files). Most of the scripts
that read the database carry an inline comment naming the Synapse ID for the table
they query, so each `dbGetQuery()` can be replaced with a `read_csv()` of the
corresponding download. There is no automatic fallback in the code; this is a
manual substitution.

---

## 4. Python environments

### 4.1 Figure scripts

`figures/python/GLASS-I_World_Map_Cohort_Distribution.py` is the only Python
figure script. It needs a geospatial stack, which is much easier to install from
conda-forge than from pip:

```bash
conda create -n glass-figures -c conda-forge python=3.11 \
    geopandas shapely matplotlib numpy pandas
conda activate glass-figures
python figures/python/GLASS-I_World_Map_Cohort_Distribution.py
```

pip equivalent, if you prefer:

```bash
python -m venv .venv && source .venv/bin/activate
pip install geopandas shapely matplotlib numpy pandas
```

`tables/python/` contains no code at present (a `.gitkeep` placeholder only).

### 4.2 The two pipeline drivers in `analysis/python/`

`mutational_signature_pipeline.sh` and `kataegis_pipeline.sh` are SLURM-ready bash
drivers that orchestrate R and Python stages. They expect **two named conda
environments** and activate them by name. The defaults are `palimpsest` (R side)
and `sigprofiler` / `--cluster-env` (Python side); override with `--r-env` and
`--sigprofiler-env` / `--cluster-env`.

```bash
# Python side: SigProfiler toolkit
conda create -n sigprofiler -c conda-forge python=3.10 pip
conda activate sigprofiler
pip install SigProfilerMatrixGenerator SigProfilerExtractor \
            SigProfilerSimulator SigProfilerClusters
# install the reference genome(s) once, ~3 GB each
python -c "from SigProfilerMatrixGenerator import install as g; g.install('GRCh37')"

# R side: Palimpsest + database access
conda create -n palimpsest -c conda-forge r-base=4.3 r-tidyverse r-dbi r-odbc \
    r-devtools bioconductor-bsgenome.hsapiens.ucsc.hg19
conda activate palimpsest
R -e 'remotes::install_github("FunGeST/Palimpsest")'
```

Both scripts load a `miniconda` module before activating, which suits an HPC with
Lmod. On a workstation, pass `--no-module-load`, or `--no-conda` if your
environments are already active.

Both are resume-aware: they write marker files under
`WORKDIR/.pipeline_markers` and helper scripts under `WORKDIR/_pipeline_scripts`,
so a re-run skips completed stages unless you pass `--force`. Start with
`--dry-run` to see the planned actions.

### 4.3 The snakemake workflows in `analysis/snakemake/`

`ampsuite.smk`, `ascat.smk` and `sigprofilertoolkit.smk` activate their own conda
environments by name inside each rule (`ampsuite`, `ascat`, `sigprofilertoolkit`)
and assume `~/bin/myconda.sh` exists to make `conda activate` available inside a
non-interactive shell. They are cluster workflows over BAM files, not part of the
figure/table path, and they depend on external tools installed separately:

| Workflow | External tool |
|---|---|
| `ampsuite.smk` | [AmpliconSuite-pipeline](https://github.com/AmpliconSuite/AmpliconSuite-pipeline) (AmpliconArchitect + AmpliconClassifier), plus its reference data and a Mosek licence |
| `ascat.smk` | [ASCAT](https://github.com/VanLoo-lab/ascat) with its GC-content and replication-timing reference files |
| `sigprofilertoolkit.smk` | the SigProfiler toolkit as in 4.2, plus the custom COSMIC-style reference signature sets documented in the header of that file |
| `GLASS-I_SvABA.sbatch` | [SvABA](https://github.com/walaj/svaba), run from a conda env named `renv` |

---

## 5. Reference genome builds

Both GRCh37/hg19 and GRCh38/hg38 appear in this repository, and that is
deliberate, not an inconsistency:

- **GRCh37 / hg19** — the GLASS variant calls themselves. The signature and
  kataegis pipelines in `analysis/python/` default to `--genome GRCh37` and load
  `BSgenome.Hsapiens.UCSC.hg19`.
- **hg38** — the copy-number and chromothripsis figure scripts in `figures/R/`
  load `BSgenome.Hsapiens.UCSC.hg38` for cytoband and segment coordinates.

Check the `--genome` flag and the `BSgenome` import at the top of a script before
assuming a build.

---

## 6. Running the code

There is no master driver; each script is standalone and is run directly.

```bash
# a table
Rscript tables/R/GLASS-I_Baseline.R

# a figure
Rscript figures/R/GLASS-I_Cohort_Description.R

# a pipeline
bash analysis/python/mutational_signature_pipeline.sh \
     --workdir /path/to/work --sql analysis/SQL/GLASS-I_dNdScv_NHM_input.sql \
     --dry-run
```

Three things to expect on a first run.

**1. Output paths are not portable.** 11 scripts write to `~/`, and a few carry
absolute paths from the original working environment
(`/vast/palmer/pi/verhaak/...`, including other users' directories). Check the
`write_csv()` / `ggsave()` calls at the bottom of a script and the file paths at
the top, and adjust both before running.

**2. Some scripts expect their SQL in your home directory.** The `.sql` files in
`analysis/SQL/` are query definitions read by the R scripts, not run on their own —
but three scripts read them from `~/` rather than from the repository:

| Script | Expects |
|---|---|
| `tables/R/GLASS-I_Drivers_Table.R` | `~/GLASS-I_CNA_Data.sql` |
| `figures/R/GLASS-I_dNdScv.R` | `~/GLASS-I_dNdScv_NHM_input.sql` and `~/GLASS-I_dNdScv_HM_input.sql` |

Copy or symlink them before running, or edit the `read_file()` paths:

```bash
ln -s "$PWD/analysis/SQL/GLASS-I_CNA_Data.sql"         ~/GLASS-I_CNA_Data.sql
ln -s "$PWD/analysis/SQL/GLASS-I_dNdScv_NHM_input.sql" ~/GLASS-I_dNdScv_NHM_input.sql
ln -s "$PWD/analysis/SQL/GLASS-I_dNdScv_HM_input.sql"  ~/GLASS-I_dNdScv_HM_input.sql
```

(`analysis/R/prepare_gistic_input.R` reads `gistic_prepare.sql` from the working
directory instead, so run it from `analysis/SQL/` or adjust the path.)

**3. Some inputs come from an earlier step, not from the database.** A few scripts
read CSVs produced by a pipeline rather than by a query — for example the
mutational-signature proportion table read by `GLASS-I_IDH_Mutations.R`,
`GLASS-I_TMZ_Cycles.R`, `GLASS-I_Univariable_KM_ID8.R` and
`GLASS-I_Summarizing_Table_Molecular_Cohort.R`. Run the corresponding pipeline in
`analysis/python/` first. Where a script reads such a file, the provenance is
noted in a trailing comment on that line.

---

## Directories this document does not cover

`legacy_codes/` and `preprocessing_pipeline/` were **not** built with the
environment described above. They were created separately, through the
standardized GLASS consortium pipeline protocols **described in the manuscript
linked in [`README.md`](README.md)** and established in the prior GLASS releases
listed there, and are included in this repository for provenance and
reproducibility rather than as code to be re-run in this environment.

- **`preprocessing_pipeline/`** — the GLASS Snakemake preprocessing workflow that
  turns raw sequencing data into the database contents: alignment, variant
  calling, copy number, and the SQL that populates the PostgreSQL schema. It
  ships its **own** conda environment specifications in
  `preprocessing_pipeline/envs/` (one YAML per tool: `align.yaml`, `gatk4.yaml`,
  `sequenza.yaml`, `pyclone.yaml`, and so on) together with its own
  `Snakefile`, cluster configuration under `conf/`, and Java dependencies under
  `jar/`. Use those specifications, not section 2 or 4 of this document.
- **`legacy_codes/`** — R, Python and Julia code from the earlier GLASS
  publications listed under *Prior releases* in `README.md`. It is retained as a
  historical record. It is not maintained against the current data release, and
  its dependencies are not part of this setup.

If you need to re-run either of those, follow the protocols in the manuscript and
the environment files shipped inside `preprocessing_pipeline/`.
