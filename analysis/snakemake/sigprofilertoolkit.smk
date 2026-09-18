## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Author: Emre Kocakavuk
## SigProfilerToolkit: unified toolkit for mutational signature analysis, combining tools for matrix generation, plotting, assignment, and extraction.
## - SigProfilerMatrixGenerator for SBS, DBS, INDEL: see https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-6041-2 (Wiki: https://osf.io/s93d5/wiki/home/)
## - SigProfilerMatrixGenerator for SVs, CNVs: see https://doi.org/10.1186/s12864-023-09584-y 
## - SigProfilerAssignment: see https://doi.org/10.1093/bioinformatics/btad756
## - SigProfilerExtractor: see https://doi.org/10.1016/j.xgen.2022.100179 
## Github repo: https://github.com/AlexandrovLab/SigProfilerToolkit
## We are using reference signatures from different publications: 
## - SBS/DBS/ID COSMIC: https://www.nature.com/articles/s41586-020-1943-3 (embedded in SigProfilerToolkit, COSMIC v3.4)
## - SBS Jin: https://www.nature.com/articles/s41588-024-01659-0; retreived from https://github.com/parklab/MuSiCal/blob/main/musical/data/COSMIC-MuSiCal_v3p2_SBS_WGS.csv & transformed into COSMIC style via /vast/palmer/pi/verhaak/shared/glass6/scripts/R/sigprofiler/reference_sig_transform.R
## - SBS Degasperi: https://www.science.org/doi/10.1126/science.abl9283; retreived from https://github.com/Nik-Zainal-Group/signature.tools.lib/blob/master/data/RefSigSBS_v2.03/RefSig_SBS_v2.03.tsv & transformed into COSMIC style via /vast/palmer/pi/verhaak/shared/glass6/scripts/R/sigprofiler/reference_sig_transform.R
## - ID Jin: https://www.nature.com/articles/s41588-024-01659-0; retreived from https://github.com/parklab/MuSiCal/blob/main/musical/data/MuSiCal_v4_Indel_WGS.csv & transformed into COSMIC style via /vast/palmer/pi/verhaak/shared/glass6/scripts/R/sigprofiler/reference_sig_transform.R
## - CN COSMIC: https://www.nature.com/articles/s41586-022-04738-6 (embedded in SigProfilerToolkit, COSMIC v3.4)
## - SV COSMIC: https://www.nature.com/articles/s41588-025-02474-x (embedded in SigProfilerToolkit, COSMIC v3.4)
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

#### SBS/DBS/ID
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Prepare input for SigProfilerMatrixGenerator for SBS, DBS, INDEL (SDI)
## From case-level geno file, take only ssMutect2 Filter = PASS + Force-include hotspot variants
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
rule sigp_matgen_sdi_prepare:
    input:
        geno = "results/mutect2/geno2db/{case_barcode}.geno.tsv"
    output:
        tsv = "results/sigprofiler/matrixgenerator/sbs_dbs_indel/input/{case_barcode}/{aliquot_barcode}/{aliquot_barcode}.SBS-DBS-INDEL.txt"
    params:
        hotspot_vcf = config["mutect2"]["given_alleles"],
        fasta       = config["reference_fasta"], 
        genome      = config["sigprofiler"]["genome_build"],
        project     = config["sigprofiler"]["project"]
    threads:
        CLUSTER_META["sigp_matgen_prepare"]["cpus-per-task"]
    log:
        "logs/sigprofiler/matgen/sdi/prepare/{case_barcode}/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/sigprofiler/matgen/sdi/prepare/{case_barcode}/{aliquot_barcode}.txt"
    shell:
        """
        set -euo pipefail
        source ~/bin/myconda.sh && conda activate sigprofilertoolkit

        Rscript scripts/R/snakemake/sigp_matgen_sdi_prepare.R \
            --geno_file    {input.geno} \
            --hotspot_vcf  {params.hotspot_vcf} \
            --fasta        {params.fasta} \
            --genome       {params.genome} \
            --project      {params.project} \
            --aliquot      {wildcards.aliquot_barcode} \
            --out_tsv      {output.tsv} \
            > {log} 2>&1
        """


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## SigProfilerMatrixGenerator for SBS, DBS, INDEL (SDI)
## Running on aliquot - level
## Complex input/output structure related to SigProfiler behavior of copying input files into separate directory and not being able to specify output directory
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
rule sigp_matgen_sdi_run:
    input:
        tsv = "results/sigprofiler/matrixgenerator/sbs_dbs_indel/input/{case_barcode}/{aliquot_barcode}/{aliquot_barcode}.SBS-DBS-INDEL.txt"
    output:
        directory("results/sigprofiler/matrixgenerator/sbs_dbs_indel/input/{case_barcode}/{aliquot_barcode}/output")
    params:
        refgenome = config["sigprofiler"]["genome_build"],
        project   = config["sigprofiler"]["project"],
        in_dir    = "results/sigprofiler/matrixgenerator/sbs_dbs_indel/input/{case_barcode}/{aliquot_barcode}"
    threads:
        CLUSTER_META["sigp_matgen_sdi_run"]["cpus-per-task"]
    log:
        "logs/sigprofiler/matgen/sdi/run/{case_barcode}/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/sigprofiler/matgen/sdi/run/{case_barcode}/{aliquot_barcode}.txt"
    shell:
        """
        set -euo pipefail
        source ~/bin/myconda.sh && conda activate sigprofilertoolkit

        SigProfilerToolkit matrix_generator matrix_generator \
            {params.project} \
            {params.refgenome} \
            {params.in_dir} \
            --plot TRUE \
            > {log} 2>&1
        """

# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# SigProfilerAssignment for SBS
# Use COSMIC (Alexandrov et al.) reference signature set
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
rule sigp_fit_sbs_cosmic:
    input:
        matrix_file = "results/sigprofiler/matrixgenerator/sbs_dbs_indel/input/{case_barcode}/{aliquot_barcode}/output/SBS/GLASS.SBS96.all"
    output:
        directory("results/sigprofiler/assignment/sbs/cosmic/{case_barcode}/{aliquot_barcode}")
    params:
        genome_build   = config["sigprofiler"]["genome_build"],
        cosmic_version = config["sigprofiler"]["cosmic_version"],
        exome_flag     = lambda wc: "True" if manifest.isExome(wc.aliquot_barcode) else "False",
        export_flag    = "True",
        context_type = "96"
    threads:
        CLUSTER_META["sigp_fit_sdi_cosmic"]["cpus-per-task"]
    log:
        "logs/sigprofiler/assignment/sbs/cosmic/{case_barcode}/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/sigprofiler/assignment/sbs/cosmic/{case_barcode}/{aliquot_barcode}.txt"
    shell:
        """
        set -euo pipefail
        source ~/bin/myconda.sh && conda activate sigprofilertoolkit

        SigProfilerAssignment cosmic_fit \
            "{input.matrix_file}" \
            "{output}" \
            --input_type matrix \
            --genome_build {params.genome_build} \
            --cosmic_version {params.cosmic_version} \
            --export_probabilities_per_mutation {params.export_flag} \
            --exome {params.exome_flag} \
            --context_type {params.context_type} \
            > {log} 2>&1
        """

# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# SigProfilerAssignment for SBS
# Use MuSiCal (Jin et al.) reference signature set
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
rule sigp_fit_sbs_jin:
    input:
        matrix_file = "results/sigprofiler/matrixgenerator/sbs_dbs_indel/input/{case_barcode}/{aliquot_barcode}/output/SBS/GLASS.SBS96.all"
    output:
        directory("results/sigprofiler/assignment/sbs/jin/{case_barcode}/{aliquot_barcode}")
    params:
        genome_build    = config["sigprofiler"]["genome_build"],
        signatures    = config["sigprofiler"]["jin_sbs_ref"],
        signature_db    = config["sigprofiler"]["jin_sbs_ref"],
        exome_flag      = lambda wc: "True" if manifest.isExome(wc.aliquot_barcode) else "False",
        export_flag     = "True",
        context_type    = "96"
    threads:
        CLUSTER_META["sigp_fit_sdi_cosmic"]["cpus-per-task"]
    log:
        "logs/sigprofiler/assignment/sbs/jin/{case_barcode}/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/sigprofiler/assignment/sbs/jin/{case_barcode}/{aliquot_barcode}.txt"
    shell:
        """
        set -euo pipefail
        source ~/bin/myconda.sh && conda activate sigprofilertoolkit

        SigProfilerAssignment decompose_fit \
            "{input.matrix_file}" \
            "{output}" \
            --input_type matrix \
            --genome_build {params.genome_build} \
            --signatures {params.signatures} \
            --signature_database {params.signature_db} \
            --export_probabilities_per_mutation {params.export_flag} \
            --exome {params.exome_flag} \
            --context_type {params.context_type} \
            > {log} 2>&1
        """

# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# SigProfilerAssignment for SBS
# Use Signal (Degasperi et al.) reference signature set
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
rule sigp_fit_sbs_degasperi:
    input:
        matrix_file = "results/sigprofiler/matrixgenerator/sbs_dbs_indel/input/{case_barcode}/{aliquot_barcode}/output/SBS/GLASS.SBS96.all"
    output:
        directory("results/sigprofiler/assignment/sbs/degasperi/{case_barcode}/{aliquot_barcode}")
    params:
        genome_build    = config["sigprofiler"]["genome_build"],
        signatures    = config["sigprofiler"]["degasperi_sbs_ref"],
        signature_db    = config["sigprofiler"]["degasperi_sbs_ref"],
        exome_flag      = lambda wc: "True" if manifest.isExome(wc.aliquot_barcode) else "False",
        export_flag     = "True",
        context_type    = "96"
    threads:
        CLUSTER_META["sigp_fit_sdi_cosmic"]["cpus-per-task"]
    log:
        "logs/sigprofiler/assignment/sbs/degasperi/{case_barcode}/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/sigprofiler/assignment/sbs/degasperi/{case_barcode}/{aliquot_barcode}.txt"
    shell:
        """
        set -euo pipefail
        source ~/bin/myconda.sh && conda activate sigprofilertoolkit

        SigProfilerAssignment decompose_fit \
            "{input.matrix_file}" \
            "{output}" \
            --input_type matrix \
            --genome_build {params.genome_build} \
            --signatures {params.signatures} \
            --signature_database {params.signature_db} \
            --export_probabilities_per_mutation {params.export_flag} \
            --exome {params.exome_flag} \
            --context_type {params.context_type} \
            > {log} 2>&1
        """

# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# SigProfilerAssignment for DBS
# Use COSMIC (Alexandrov et al.) reference signature set
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
rule sigp_fit_dbs_cosmic:
    input:
        matrix_file = "results/sigprofiler/matrixgenerator/sbs_dbs_indel/input/{case_barcode}/{aliquot_barcode}/output/DBS/GLASS.DBS78.all"
    output:
        directory("results/sigprofiler/assignment/dbs/cosmic/{case_barcode}/{aliquot_barcode}")
    params:
        genome_build   = config["sigprofiler"]["genome_build"],
        cosmic_version = config["sigprofiler"]["cosmic_version"],
        exome_flag     = lambda wc: "True" if manifest.isExome(wc.aliquot_barcode) else "False",
        export_flag    = "True",
        context_type = "DINUC", 
        collapse = "False"
    threads:
        CLUSTER_META["sigp_fit_sdi_cosmic"]["cpus-per-task"]
    log:
        "logs/sigprofiler/assignment/dbs/cosmic/{case_barcode}/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/sigprofiler/assignment/dbs/cosmic/{case_barcode}/{aliquot_barcode}.txt"
    shell:
        """
        set -euo pipefail
        source ~/bin/myconda.sh && conda activate sigprofilertoolkit

        SigProfilerAssignment cosmic_fit \
            "{input.matrix_file}" \
            "{output}" \
            --input_type matrix \
            --genome_build {params.genome_build} \
            --cosmic_version {params.cosmic_version} \
            --export_probabilities_per_mutation {params.export_flag} \
            --exome {params.exome_flag} \
            --context_type {params.context_type} \
            --collapse_to_SBS96 {params.collapse} \
            > {log} 2>&1
        """

# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# SigProfilerAssignment for INDEL
# Use COSMIC (Alexandrov et al.) reference signature set
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
rule sigp_fit_indel_cosmic:
    input:
        matrix_file = "results/sigprofiler/matrixgenerator/sbs_dbs_indel/input/{case_barcode}/{aliquot_barcode}/output/ID/GLASS.ID83.all"
    output:
        directory("results/sigprofiler/assignment/indel/cosmic/{case_barcode}/{aliquot_barcode}")
    params:
        genome_build   = config["sigprofiler"]["genome_build"],
        cosmic_version = config["sigprofiler"]["cosmic_version"],
        exome_flag     = lambda wc: "True" if manifest.isExome(wc.aliquot_barcode) else "False",
        export_flag    = "True",
        context_type = "ID", 
        collapse = "False"
    threads:
        CLUSTER_META["sigp_fit_sdi_cosmic"]["cpus-per-task"]
    log:
        "logs/sigprofiler/assignment/indel/cosmic/{case_barcode}/{aliquot_barcode}.log"
    benchmark:
        ("benchmarks/sigprofiler/assignment/indel/cosmic/{case_barcode}/{aliquot_barcode}.txt")
    shell:
        """
        set -euo pipefail
        source ~/bin/myconda.sh && conda activate sigprofilertoolkit

        SigProfilerAssignment cosmic_fit \
            "{input.matrix_file}" \
            "{output}" \
            --input_type matrix \
            --genome_build {params.genome_build} \
            --cosmic_version {params.cosmic_version} \
            --export_probabilities_per_mutation {params.export_flag} \
            --exome {params.exome_flag} \
            --context_type {params.context_type} \
            --collapse_to_SBS96 {params.collapse} \
            > {log} 2>&1
        """

# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# SigProfilerAssignment for INDEL
# Use MuSiCal (Jin et al.) reference signature set
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
rule sigp_fit_indel_jin:
    input:
        matrix_file = "results/sigprofiler/matrixgenerator/sbs_dbs_indel/input/{case_barcode}/{aliquot_barcode}/output/ID/GLASS.ID83.all"
    output:
        directory("results/sigprofiler/assignment/indel/jin/{case_barcode}/{aliquot_barcode}")
    params:
        genome_build    = config["sigprofiler"]["genome_build"],
        signatures    = config["sigprofiler"]["jin_indel_ref"],
        signature_db    = config["sigprofiler"]["jin_indel_ref"],
        exome_flag      = lambda wc: "True" if manifest.isExome(wc.aliquot_barcode) else "False",
        export_flag     = "True",
        context_type    = "ID",
        collapse = "False"
    threads:
        CLUSTER_META["sigp_fit_sdi_cosmic"]["cpus-per-task"]
    log:
        "logs/sigprofiler/assignment/indel/jin/{case_barcode}/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/sigprofiler/assignment/indel/jin/{case_barcode}/{aliquot_barcode}.txt"
    shell:
        """
        set -euo pipefail
        source ~/bin/myconda.sh && conda activate sigprofilertoolkit

        SigProfilerAssignment decompose_fit \
            "{input.matrix_file}" \
            "{output}" \
            --input_type matrix \
            --genome_build {params.genome_build} \
            --signatures {params.signatures} \
            --signature_database {params.signature_db} \
            --export_probabilities_per_mutation {params.export_flag} \
            --exome {params.exome_flag} \
            --context_type {params.context_type} \
            > {log} 2>&1
        """

### CNV
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Prepare input for SigProfilerMatrixGenerator for CNV
## From case-level RefPhase phased CNV seg file
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
rule sigp_matgen_cnv_prepare:
    input:
        seg_file = "results/refphase/{case_barcode}/{case_barcode}-phased-segments.tsv"
    output:
        "results/sigprofiler/matrixgenerator/cnv/input/{case_barcode}/{aliquot_barcode}/{aliquot_barcode}.CNV.txt"
    params:
        pair_barcode = lambda wc: next(p for p in manifest.getPairsByCase(wc.case_barcode) if manifest.getTumor(p) == wc.aliquot_barcode)
    threads:
        CLUSTER_META["sigp_matgen_prepare"]["cpus-per-task"]
    log:
        "logs/sigprofiler/matgen/cnv/prepare/{case_barcode}/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/sigprofiler/matgen/cnv/prepare/{case_barcode}/{aliquot_barcode}.txt"
    shell:
        """
        set -euo pipefail
        source ~/bin/myconda.sh && conda activate sigprofilertoolkit

        Rscript scripts/R/snakemake/sigp_matgen_cnv_prepare.R \
            --seg_file  "{input.seg_file}" \
            --pair      "{params.pair_barcode}" \
            --aliquot   "{wildcards.aliquot_barcode}" \
            --out_txt   "{output}" \
            > {log} 2>&1
        """

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## SigProfilerMatrixGenerator for CNV
## Running on aliquot - level
## Complex input/output structure related to SigProfiler behavior of copying input files into separate directory and not being able to specify output directory
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
rule sigp_matgen_cnv_run:
    input:
        tsv = "results/sigprofiler/matrixgenerator/cnv/input/{case_barcode}/{aliquot_barcode}/{aliquot_barcode}.CNV.txt"
    output:
        # tsv = "results/sigprofiler/matrixgenerator/cnv/input/{case_barcode}/{aliquot_barcode}/output/GLASS.CNV48.matrix.tsv"
        out_dir = directory("results/sigprofiler/matrixgenerator/cnv/input/{case_barcode}/{aliquot_barcode}/output")
    params:
        refgenome = config["sigprofiler"]["genome_build"],
        project   = config["sigprofiler"]["project"],
        file_type      = "ASCAT", 
        out_dir = directory("results/sigprofiler/matrixgenerator/cnv/input/{case_barcode}/{aliquot_barcode}/output")
    threads:
        CLUSTER_META["sigp_matgen_sdi_run"]["cpus-per-task"]
    log:
        "logs/sigprofiler/matgen/cnv/run/{case_barcode}/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/sigprofiler/matgen/cnv/run/{case_barcode}/{aliquot_barcode}.txt"
    shell:
        """
        set -euo pipefail
        source ~/bin/myconda.sh && conda activate sigprofilertoolkit

        SigProfilerToolkit matrix_generator cnv_matrix_generator \
            {params.file_type} \
            {input} \
            {params.project} \
            {params.out_dir} \
            > {log} 2>&1
        """

# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# SigProfilerAssignment for CNV
# Use COSMIC (Alexandrov et al.) reference signature set
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
rule sigp_fit_cnv_cosmic:
    input:
        matrix_file = "results/sigprofiler/matrixgenerator/cnv/input/{case_barcode}/{aliquot_barcode}/output/GLASS.CNV48.matrix.tsv"
    output:
        directory("results/sigprofiler/assignment/cnv/cosmic/{case_barcode}/{aliquot_barcode}")
    params:
        genome_build    = config["sigprofiler"]["genome_build"],
        cosmic_version = config["sigprofiler"]["cosmic_version"],
        exome_flag      = lambda wc: "True" if manifest.isExome(wc.aliquot_barcode) else "False",
        export_flag     = "True",
        collapse = "False"
    threads:
        CLUSTER_META["sigp_fit_sdi_cosmic"]["cpus-per-task"]
    log:
        "logs/sigprofiler/assignment/cnv/cosmic/{case_barcode}/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/sigprofiler/assignment/cnv/cosmic/{case_barcode}/{aliquot_barcode}.txt"
    shell:
        """
        set -euo pipefail
        source ~/bin/myconda.sh && conda activate sigprofilertoolkit

        SigProfilerAssignment cosmic_fit \
            {input.matrix_file} \
            {output} \
            --input_type matrix \
            --genome_build {params.genome_build} \
            --cosmic_version {params.cosmic_version} \
            --export_probabilities_per_mutation {params.export_flag} \
            --exome {params.exome_flag} \
            --collapse_to_SBS96 {params.collapse} \
            > {log} 2>&1
        """

### SV
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Prepare input for SigProfilerMatrixGenerator for SV
## From sample-level SvABA SV file (aso used for ShatterSeek input)
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
rule sigp_matgen_sv_prepare:
    input:
        sv_tsv = "results/shatterseek/input_sv/{case_barcode}/{aliquot_barcode}.sv.tsv"
    output:
        bedpe = "results/sigprofiler/matrixgenerator/sv/input/{case_barcode}/{aliquot_barcode}/{aliquot_barcode}.bedpe"
    threads:
        CLUSTER_META["sigp_matgen_prepare"]["cpus-per-task"]
    log:
        "logs/sigprofiler/matgen/sv/prepare/{case_barcode}/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/sigprofiler/matgen/sv/prepare/{case_barcode}/{aliquot_barcode}.txt"
    shell:
        """
        set -euo pipefail
        source ~/bin/myconda.sh && conda activate sigprofilertoolkit

        Rscript scripts/R/snakemake/sigp_matgen_sv_prepare.R \
          --sv_tsv {input.sv_tsv} \
          --aliquot {wildcards.aliquot_barcode} \
          --out_tsv {output.bedpe} \
          > {log} 2>&1
        """

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## SigProfilerMatrixGenerator for SV
## Running on aliquot - level
## Complex input/output structure related to SigProfiler behavior of copying input files into separate directory and not being able to specify output directory
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
rule sigp_matgen_sv_run:
    input:
        bedpe = "results/sigprofiler/matrixgenerator/sv/input/{case_barcode}/{aliquot_barcode}/{aliquot_barcode}.bedpe"
    output:
        out_dir = directory("results/sigprofiler/matrixgenerator/sv/input/{case_barcode}/{aliquot_barcode}/output")
    params:
        in_dir  = "results/sigprofiler/matrixgenerator/sv/input/{case_barcode}/{aliquot_barcode}",
        out_dir = "results/sigprofiler/matrixgenerator/sv/input/{case_barcode}/{aliquot_barcode}/output",
        project = config["sigprofiler"]["project"]
    threads:
        CLUSTER_META["sigp_matgen_sdi_run"]["cpus-per-task"]
    log:
        "logs/sigprofiler/matgen/sv/run/{case_barcode}/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/sigprofiler/matgen/sv/run/{case_barcode}/{aliquot_barcode}.txt"
    shell:
        """
        set -euo pipefail
        source ~/bin/myconda.sh && conda activate sigprofilertoolkit

        SigProfilerToolkit matrix_generator sv_matrix_generator \
            {params.in_dir} \
            {params.project} \
            {params.out_dir}
          > {log} 2>&1
        """

# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# SigProfilerAssignment for SV
# Use COSMIC (Alexandrov et al.) reference signature set
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
rule sigp_fit_sv_cosmic:
    input:
        matrix_file = "results/sigprofiler/matrixgenerator/sv/input/{case_barcode}/{aliquot_barcode}/output/GLASS.SV32.matrix.tsv"
    output:
        directory("results/sigprofiler/assignment/sv/cosmic/{case_barcode}/{aliquot_barcode}")
    params:
        genome_build = config["sigprofiler"]["genome_build"],
        cosmic_version = config["sigprofiler"]["cosmic_version"],
        exome_flag = lambda wc: "True" if manifest.isExome(wc.aliquot_barcode) else "False",
        export_flag = "True",
        collapse = "False"
    threads:
        CLUSTER_META["sigp_fit_sdi_cosmic"]["cpus-per-task"]
    log:
        "logs/sigprofiler/assignment/sv/cosmic/{case_barcode}/{aliquot_barcode}.log"
    benchmark:
        "benchmarks/sigprofiler/assignment/sv/cosmic/{case_barcode}/{aliquot_barcode}.txt"
    shell:
        """
        set -euo pipefail
        source ~/bin/myconda.sh && conda activate sigprofilertoolkit

        SigProfilerAssignment cosmic_fit \
            {input.matrix_file} \
            {output} \
            --input_type matrix \
            --genome_build {params.genome_build} \
            --cosmic_version {params.cosmic_version} \
            --export_probabilities_per_mutation {params.export_flag} \
            --exome {params.exome_flag} \
            --collapse_to_SBS96 {params.collapse} \
            > {log} 2>&1
        """