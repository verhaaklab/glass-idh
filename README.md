![GLASS_logo_072219_08](https://user-images.githubusercontent.com/6731211/64618915-2ed59e80-d3af-11e9-8983-d41414379ad3.png)

## The GLASS consortium

### Overview
The Glioma Longitudinal AnalySiS (GLASS) consortium consists of clinical, bioinformaticians, and basic science researchers from leading institutions across the world striving to better understand glioma tumor evolution and to expose its therapeutic vulnerabilities. The code in this respository was used to generate the figures and perform the analyses described in our [2025 publication](https://www.biorxiv.org/content/10.1101/2025.07.11.664189v1.full) in BioRxiv.

R code used to make each of the figures can be found in the figures/R subdirectory.

### Data Release version 2025.
The data analyzed using this code was part of the fourth release of the GLASS dataset which was made available in 2025. This dataset is managed internally using the PostgreSQL database management system. These data are under active curation so future versions will include additional data as well as correct potential errors.

### Data Download
The GLASS data can be downloaded from the `Tables` page [here](https://www.synapse.org/#!Synapse:syn17038081/tables/) and the `Files` page [here](https://www.synapse.org/#!Synapse:syn26465623). It is also possible to query the data directly using the the API by using queries. You can read more about that [here](https://docs.synapse.org/articles/tables.html).

### Prior releases
Previous projects and associated repositories are listed below:

Barthel FP, Johnson KC, Varn FS, et al. Longitudinal molecular trajectories of diffuse glioma in adults. Nature. 2019; 576(7785): 112-120. [Paper.](https://www.nature.com/articles/s41586-019-1775-1) [Code.](https://github.com/TheJacksonLaboratory/GLASS)

Varn FS, Johnson KC, Martinek J, et al. Glioma progression is shaped by genetic evolution and microenvironment interactions. Cell. 2022; 185(12): 2184-2199.e16. [Paper.](https://www.sciencedirect.com/science/article/pii/S0092867422005360?via%3Dihub) [Code.](https://github.com/fsvarn/GLASSx)

---

## Repository layout

The repository separates **the code behind this manuscript** from **the standardized
GLASS infrastructure it was built on top of**. The first two directories below are
the former; the last two are the latter.

| Directory | What it holds |
|---|---|
| [`analysis/`](#analysis) | Analyses that generate derived results — copy number, structural variants, mutational signatures |
| [`figures/`](#figures) | One script per manuscript figure panel |
| [`tables/`](#tables) | One script per manuscript table, mostly cohort description and survival models |
| [`legacy_codes/`](#legacy_codes) | Code from the prior GLASS publications listed above, retained for provenance |
| [`preprocessing_pipeline/`](#preprocessing_pipeline) | The GLASS workflow that produced the database itself |

Environment and database configuration for `analysis/`, `figures/` and `tables/`
is documented separately in [`SETUP.md`](SETUP.md).

---

### `analysis/`

Derived analyses: everything that computes a result which the figures and tables
then read. 23 files in four subdirectories.

| Subdirectory | Files | Contents |
|---|---|---|
| `analysis/R/` | 12 | Chromothripsis calling (ShatterSeek), ecDNA/amplicon analysis (AmpliconSuite), allele-specific copy number (ASCAT), whole-genome doubling, fraction-of-genome-altered calculation, GISTIC input preparation, copy-number signature exploration, telomere-length normalisation |
| `analysis/SQL/` | 4 | Query definitions consumed by the scripts and pipelines — CNA extraction, the hypermutant and non-hypermutant dN/dS inputs, GISTIC preparation |
| `analysis/python/` | 4 | Two self-contained, resume-aware pipeline drivers — `mutational_signature_pipeline.sh` (repeat-region filtering → SigProfilerMatrixGenerator/Extractor → Palimpsest deconvolution) and `kataegis_pipeline.sh` (clustered-mutation / kataegis calling via SigProfilerSimulator + SigProfilerClusters) — plus an AmpliconSuite manifest builder and the SvABA SLURM script |
| `analysis/snakemake/` | 3 | Cluster workflows over BAM files: `ampsuite.smk`, `ascat.smk`, `sigprofilertoolkit.smk` |

Note that the `analysis/python/` drivers are bash, not Python — they orchestrate
R and Python stages and are named for the toolchain they drive.

### `figures/`

One script per figure panel, 28 files. Each is standalone: it opens the GLASS
database, builds its own plotting table, and writes the panel.

| Subdirectory | Files | Contents |
|---|---|---|
| `figures/R/` | 27 | Cohort description and QC Sankey; swimmer's plot; copy-number heatmap and per-subtype CN segment plots (astrocytoma and oligodendroglioma, primary / recurrence / delta); fraction-of-genome-altered, with and without hypermutants and by grade; chromosomal-instability landscape and CIN markers; IDH mutation figures; mutation burden; mutational-signature boxplots; dN/dS selection; sSNV mutual exclusivity; TMZ cycles; and Kaplan–Meier curves for CNA, drivers (GENIE), hypermutation, ID8 and sSNV burden |
| `figures/python/` | 1 | The world map of cohort contributing centres |

### `tables/`

One script per manuscript table, 16 files, all R.

| Group | Contents |
|---|---|
| Cohort description | Baseline table, molecular cohort summary, prior-treatment table |
| Drivers | Driver gene table, drivers by prior treatment, IDH mutation table, McNemar test on paired CNA changes |
| Survival | Univariable and multivariable models for three endpoints — overall survival (OS), post-recurrence survival (PRS) and time to recurrence (TTR) — each split into clinical and molecular covariates |

Output is formatted with `gt` and `gtsummary`.

### `legacy_codes/`

**Created separately, outside the scope of this manuscript's analysis code**, under
the same standardized GLASS pipeline protocols. 333 files of R (308), Python (19)
and Julia (5) from the earlier GLASS publications listed under *Prior releases*
above, retained for provenance and reproducibility of those papers.

Organised by topic under `legacy_codes/R/`: single-cell deconvolution
(CIBERSORTx), neoantigen analysis, SNV and CNV processing, manifest handling,
mutation timing, figures and tables for the 2019 and 2022 papers, and a Shiny
app. The Julia code covers subclonal selection analysis.

This directory is **not** maintained against the current data release and its
dependencies are not part of [`SETUP.md`](SETUP.md).

### `preprocessing_pipeline/`

**Created separately, through the standardized GLASS consortium pipeline protocols
described in the [2025 publication](https://www.biorxiv.org/content/10.1101/2025.07.11.664189v1.full)
linked above and established in the prior GLASS releases listed under *Prior
releases*.** 251 files. This is the workflow that produced the PostgreSQL database
the rest of the repository queries — it is not re-run as part of reproducing a
figure, and it is not covered by [`SETUP.md`](SETUP.md).

| Subdirectory | Files | Contents |
|---|---|---|
| `sql/` | 125 | The schema and queries that populate and serve the database, grouped by domain: `clinical`, `variants`, `snv`, `cnv`, `drivers`, `mut_freq`, `mut_sig`, `pyclone`, `timing`, `neutrality`, `neoag`, `immune`, `lohhla`, `expression`, `heatmap`, `metadata`, `set`, `dndscv`, `figures` |
| `snakemake/` | 30 | Per-stage workflow modules, driven by the top-level `Snakefile` |
| `envs/` | 44 | One conda environment specification per tool — alignment, GATK4, Sequenza, PyClone, Manta, Delly, LUMPY, OptiType, pVACseq, MiXCR, Kallisto, and others |
| `conf/` | — | Cluster and pipeline configuration |
| `bin/`, `jar/` | — | Bundled executables and Java dependencies |
| `dag/`, `dbm/` | — | Workflow DAGs and the database schema diagram |

Because this directory carries its own environment specifications in `envs/`, use
those rather than the instructions in [`SETUP.md`](SETUP.md) if you need to re-run
any part of it.

---

## Where to start

- **Reproducing a figure or table** → read [`SETUP.md`](SETUP.md), configure the
  database connection, then run the single script for that panel.
- **Reproducing a derived analysis** → the corresponding script or pipeline in
  `analysis/`; the two drivers in `analysis/python/` accept `--dry-run`.
- **Understanding how the data was generated** → `preprocessing_pipeline/` and the
  methods of the manuscript.
- **Looking for code from the 2019 or 2022 papers** → `legacy_codes/`, or the
  dedicated repositories linked under *Prior releases*.
