## Authors: Emre Kocakavuk
## New implementation of AmpliconSuite for detection of chromothripsis and extrachromosomal DNA (ecDNA) elements
## See: https://github.com/AmpliconSuite/AmpliconSuite-pipeline
## See: https://www.biorxiv.org/content/10.1101/2024.05.06.592768v1
## AmpliconSuite combines the following two steps:
## - AmpliconArchitect (AA)
## - AmpliconClassifier (AC)
## Note: Due to the longitudinal nature of GLASS, we are running AmpliconSuite in multisample mode starting from BAM files

## Python bad interpreter issue (also noted on Github), quick fix in our case:
# conda activate ampsuite
# for stub in "$CONDA_PREFIX"/bin/*AmpSuite*.py; do
#     sed -i '1c #!/usr/bin/env python3' "$stub"
# done

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Prepare AmpliconSuite manifest on a case-level basis
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
rule ampsuite_manifest:
    input:
        tumor  = lambda wc: ancient(expand("results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam", aliquot_barcode = manifest.getTumorByCase(wc.case_barcode))), 
        normal = lambda wc: ancient(expand("results/align/bqsr/{aliquot_barcode}.realn.mdup.bqsr.bam", aliquot_barcode = manifest.getNormalByCase(wc.case_barcode)))
    output:
        tsv = "results/ampsuite/manifest/{case_barcode}.tsv"
    params:
        tumor_bams  = lambda wc, input, **_: " ".join(input.tumor),
        normal_bams = lambda wc, input, **_: " ".join(input.normal)
    threads: 
        CLUSTER_META["ampsuite_manifest"]["cpus-per-task"]
    log: 
        "logs/ampsuite/manifest/{case_barcode}.log"
    benchmark: 
        "benchmarks/ampsuite/manifest/{case_barcode}.txt"
    message: 
        "Prepare Manifest for AmpliconSuite\n"
        "Case: {wildcards.case_barcode}"
    shell:
        """
        echo "$-"
        source ~/bin/myconda.sh && conda activate ampsuite
        echo "$-"
        python scripts/python/snakemake/ampsuite_manifest.py \
               --tumor "{params.tumor_bams}" \
               --normal "{params.normal_bams}" \
               --output {output.tsv} \
        > {log} 2>&1
        """

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Run AmpliconSuite with AmpliconArchitect & AmpliconClassifier
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
rule ampsuite_run:
    input:
        tsv = "results/ampsuite/manifest/{case_barcode}.tsv"
    output:
        done = "results/ampsuite/run/{case_barcode}/AA.done"
    params:
        outdir     = lambda wc: f"results/ampsuite/run/{wc.case_barcode}",
        reference  = "GRCh37"
    threads:
        CLUSTER_META["ampsuite_run"]["cpus-per-task"]
    log:
        "logs/ampsuite/run/{case_barcode}.log"
    benchmark:
        "benchmarks/ampsuite/run/{case_barcode}.txt"
    message:
        "Grouped AmpliconSuite\n"
        "Case: {wildcards.case_barcode}"
    shell:
        """
        echo "$-"
        source ~/bin/myconda.sh && conda activate ampsuite
        echo "$-"
        GroupedAnalysisAmpSuite.py \
            -i {input.tsv} \
            -o {params.outdir} \
            -t {threads} \
            --ref {params.reference} \
        > {log} 2>&1
        touch {output.done}
        """