## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Author: Emre Kocakavuk
## ASCAT (Allele-Specific Copy Number Analysis of Tumors) 
## - Preparation involves transformation of GATK CNV output into ASCAT-readable input via GATK2ASCAT
## - Creation of allele-specific copy number calls via ASCAT
## - Additional informative metrics: Purity, Ploidy, WGD
## See: https://academic.oup.com/bioinformatics/article/37/13/1909/5843787 
## Github repo: https://github.com/VanLoo-lab/ascat/tree/master 
## Modules partly adapted from: https://github.com/nf-core/modules/blob/733023d250311ee76c46d6863a4e056f9855eb5d/modules/nf-core/ascat/main.nf 
## Update EK 01/25: ASCAT comes with a multi-sample caller that has several advantages for longitudinal datasets (improved identification of shared vs. private breakpoints and segments)
## See: https://academic.oup.com/bioinformatics/article/37/13/1909/5843787 
## Github repo: https://github.com/VanLoo-lab/ascat/blob/master/ASCAT/vignettes/asmultipcf-vignette.Rmd (particularly adapting the asmultipcf function)
## ASCAT single sample run: `single_ascat`
## ASCAT case based multi sample run: `multi_ascat`
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
rule single_ascat: 
    input: 
        baft    = "results/gatk2ascat/{pair_barcode}/{pair_barcode}.BAF.tumor.txt",
        bafn    = "results/gatk2ascat/{pair_barcode}/{pair_barcode}.BAF.normal.txt",
        logrt   = "results/gatk2ascat/{pair_barcode}/{pair_barcode}.LogR.tumor.txt",
        logrn   = "results/gatk2ascat/{pair_barcode}/{pair_barcode}.LogR.normal.txt"
    output: 
        qc   = "results/ascat/single/qc/{pair_barcode}.QC_metrics.txt",
        seg  = "results/ascat/single/seg/{pair_barcode}.segments.txt"
    params: 
        gender  = lambda wildcards: manifest.getSex(manifest.getTumor(wildcards.pair_barcode)) if manifest.getSex(manifest.getTumor(wildcards.pair_barcode)) is not None else "NA",
        gcfile = lambda wc: config[f"ascat_{'wes' if manifest.isExome(manifest.getTumor(wc.pair_barcode)) else 'wgs'}"]["GCfile"],
        rtfile = lambda wc: config[f"ascat_{'wes' if manifest.isExome(manifest.getTumor(wc.pair_barcode)) else 'wgs'}"]["RTfile"],
        pdfdir = "results/ascat/single/pdf",
        mode = "single"
    threads: 
        CLUSTER_META["single_ascat"]["cpus-per-task"]
    log: 
        "logs/ascat/single/{pair_barcode}.log"
    benchmark: 
        "benchmarks/ascat/single/{pair_barcode}.txt"
    message: 
        "Run ASCAT\n"
        "Pair: {wildcards.pair_barcode}"
    shell: 
        """
        if [ ! -d "{params.pdfdir}" ]; then
            mkdir -p "{params.pdfdir}"
            echo "Directory created at: {params.pdfdir}" >> {log}
        else
            echo "Directory already exists: {params.pdfdir}" >> {log}
        fi
        echo "$-"
        source ~/bin/myconda.sh && conda activate ascat
        echo "$-"
        Rscript scripts/R/snakemake/ascat_analysis.R \
            {params.mode} \
            {input.baft} \
            {input.logrt} \
            {input.bafn} \
            {input.logrn} \
            {params.gender} \
            {params.gcfile} \
            {params.rtfile} \
            {params.pdfdir} \
            {output.qc} \
            {output.seg}> {log} 2>&1
        """

rule multi_ascat_merge:
    input:
        tumor_baf = lambda wildcards: expand("results/gatk2ascat/{pair_barcode}/{pair_barcode}.BAF.tumor.txt", pair_barcode = manifest.getPairsByCase(wildcards.case_barcode)),
        tumor_logr = lambda wildcards: expand("results/gatk2ascat/{pair_barcode}/{pair_barcode}.LogR.tumor.txt", pair_barcode = manifest.getPairsByCase(wildcards.case_barcode)),
        normal_baf = lambda wildcards: expand("results/gatk2ascat/{pair_barcode}/{pair_barcode}.BAF.normal.txt", pair_barcode = manifest.getPairsByCase(wildcards.case_barcode)),
        normal_logr = lambda wildcards: expand("results/gatk2ascat/{pair_barcode}/{pair_barcode}.LogR.normal.txt", pair_barcode = manifest.getPairsByCase(wildcards.case_barcode))
    output:
        tumor_baf_merged  = "results/ascat/multi/merge/{case_barcode}/{case_barcode}.merged.BAF.tumor.txt",
        tumor_logr_merged = "results/ascat/multi/merge/{case_barcode}/{case_barcode}.merged.LogR.tumor.txt",
        normal_baf_merged = "results/ascat/multi/merge/{case_barcode}/{case_barcode}.merged.BAF.normal.txt",
        normal_logr_merged= "results/ascat/multi/merge/{case_barcode}/{case_barcode}.merged.LogR.normal.txt"
    params:
        tbaf_tsv  = lambda wildcards, input: ",".join(input.tumor_baf),
        tlogr_tsv = lambda wildcards, input: ",".join(input.tumor_logr),
        nbaf_tsv  = lambda wildcards, input: input.normal_baf[0],
        nlogr_tsv = lambda wildcards, input: input.normal_logr[0]
    threads:
        CLUSTER_META["multi_ascat_merge"]["cpus-per-task"]
    log:
        "logs/ascat/multi/{case_barcode}.merge.log"
    shell:
        """
        echo "$-"
        source ~/bin/myconda.sh && conda activate ascat
        echo "$-"
        Rscript scripts/R/snakemake/ascat_merge.R \
            "{params.tbaf_tsv}" \
            "{params.tlogr_tsv}" \
            "{params.nbaf_tsv}" \
            "{params.nlogr_tsv}" \
            "{output.tumor_baf_merged}" \
            "{output.tumor_logr_merged}" \
            "{output.normal_baf_merged}" \
            "{output.normal_logr_merged}" \
            >> {log} 2>&1
        """

rule multi_ascat:
    input:
        tumor_baf_merged  = "results/ascat/multi/merge/{case_barcode}/{case_barcode}.merged.BAF.tumor.txt",
        tumor_logr_merged = "results/ascat/multi/merge/{case_barcode}/{case_barcode}.merged.LogR.tumor.txt",
        normal_baf_merged = "results/ascat/multi/merge/{case_barcode}/{case_barcode}.merged.BAF.normal.txt",
        normal_logr_merged= "results/ascat/multi/merge/{case_barcode}/{case_barcode}.merged.LogR.normal.txt"
    output:
        qc  = "results/ascat/multi/qc/{case_barcode}.QC_metrics.txt",
        seg = "results/ascat/multi/seg/{case_barcode}.segments.txt"
    params:
        gcfile = lambda wc: config[f"ascat_{'wes' if manifest.isExome(manifest.getTumor(manifest.getPairsByCase(wc.case_barcode)[0])) else 'wgs'}"]["GCfile"],
        rtfile = lambda wc: config[f"ascat_{'wes' if manifest.isExome(manifest.getTumor(manifest.getPairsByCase(wc.case_barcode)[0])) else 'wgs'}"]["RTfile"],
        pdfdir = "results/ascat/multi/pdf",
        gender = lambda wildcards: ",".join(manifest.getSex(manifest.getTumor(pb)) if manifest.getSex(manifest.getTumor(pb)) is not None else "NA" for pb in manifest.getPairsByCase(wildcards.case_barcode)),
        mode = "multi"
    threads:
        CLUSTER_META["multi_ascat"]["cpus-per-task"]
    log:
        "logs/ascat/multi/{case_barcode}.run.log"
    shell:
        """
        if [ ! -d "{params.pdfdir}" ]; then
            mkdir -p "{params.pdfdir}"
            echo "Directory created at: {params.pdfdir}" >> {log}
        else
            echo "Directory already exists: {params.pdfdir}" >> {log}
        fi
        echo "$-"
        source ~/bin/myconda.sh && conda activate ascat
        echo "$-"
        Rscript scripts/R/snakemake/ascat_analysis.R \
            {params.mode} \
            {input.tumor_baf_merged} \
            {input.tumor_logr_merged} \
            {input.normal_baf_merged} \
            {input.normal_logr_merged} \
            {params.gender} \
            {params.gcfile} \
            {params.rtfile} \
            {params.pdfdir} \
            {output.qc} \
            {output.seg} \
            >> {log} 2>&1
        """




# # ## Previous:
# # ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# # ## Author: Emre Kocakavuk
# # ## ASCAT (Allele-Specific Copy Number Analysis of Tumors) 
# # ## - Preparation involves transformation of GATK CNV output into ASCAT-readable input via GATK2ASCAT
# # ## - Creation of allele-specific copy number calls via ASCAT
# # ##- Additional informative metrics: Purity, Ploidy, WGD
# # ## See: https://academic.oup.com/bioinformatics/article/37/13/1909/5843787 
# # ## Github repo: https://github.com/VanLoo-lab/ascat/tree/master 
# # ## Modules partly adapted from: https://github.com/nf-core/modules/blob/733023d250311ee76c46d6863a4e056f9855eb5d/modules/nf-core/ascat/main.nf 
# # ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# # rule ascat_run: 
# #     input: 
# #         baft    = "results/gatk2ascat/{pair_barcode}/{pair_barcode}.BAF.tumor.txt",
# #         bafn    = "results/gatk2ascat/{pair_barcode}/{pair_barcode}.BAF.normal.txt",
# #         logrt   = "results/gatk2ascat/{pair_barcode}/{pair_barcode}.LogR.tumor.txt",
# #         logrn   = "results/gatk2ascat/{pair_barcode}/{pair_barcode}.LogR.normal.txt"
# #     output: 
# #         qc   = "results/ascat/qc/{pair_barcode}.QC_metrics.txt",
# #         seg  = "results/ascat/seg/{pair_barcode}.segments.txt"
# #     params: 
# #         gender  = lambda wildcards: manifest.getSex(manifest.getTumor(wildcards.pair_barcode)) if manifest.getSex(manifest.getTumor(wildcards.pair_barcode)) is not None else "NA",
# #         outdir  = "results/ascat", 
# #         gversion = config['ascat']['genomeVersion'], 
# #         gcfile  = config['ascat']['GCfile'], 
# #         rtfile  = config['ascat']['RTfile'],
# #         pdfdir = "results/ascat/pdf"
# #     threads: 
# #         CLUSTER_META["ascat_run"]["cpus-per-task"]
# #     log: 
# #         "logs/ascat/{pair_barcode}.log"
# #     benchmark: 
# #         "benchmarks/ascat/{pair_barcode}.txt"
# #     message: 
# #         "Run ASCAT\n"
# #         "Pair: {wildcards.pair_barcode}"
# #     shell: 
# #         """
# #         if [ ! -d "{params.pdfdir}" ]; then
# #             mkdir -p "{params.pdfdir}"
# #             echo "Directory created at: {params.pdfdir}" >> {log}
# #         else
# #             echo "Directory already exists: {params.pdfdir}" >> {log}
# #         fi
# #         echo "$-"
# #         source ~/bin/myconda.sh && conda activate ascat
# #         echo "$-"
# #         Rscript scripts/R/snakemake/ascat_analysis.R {input.baft} {input.bafn} {input.logrt} {input.logrn} {params.gender} {params.outdir} {params.gversion} {params.gcfile} {params.rtfile} {params.pdfdir} {output.qc} {output.seg}> {log} 2>&1
# #         """

# #         # """
# #         # echo "$-"
# #         # source ~/bin/myconda.sh && conda activate ascat
# #         # echo "$-"
# #         # Rscript scripts/R/snakemake/ascat_analysis.R \
# #         #     --tumor_bam {input.tumor}
# #         #     --normal_bam {input.normal}
# #         #     --output_prefix
# #         #     --allelecounter_exe {config[ascat][allelecounter_dir]} \
# #         #     --alleles_prefix {config[ascat][alleles_prefix]} \
# #         #     --loci_prefix {config[ascat][loci_prefix]} \
# #         #     --gender {params.gender}
# #         #     --genomeVersion {config[ascat][genomeVersion]}
# #         #     --nthreads
# #         #     > {log} 2>&1
# #         # """