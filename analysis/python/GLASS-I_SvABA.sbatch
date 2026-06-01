# GLASS-I: SvABA SV Calling
# Author: C.M.S. Tesileanu
# Date: 2026-05-27

#!/usr/bin/env bash
#SBATCH --job-name=svaba                
#SBATCH --chdir=~/results/svaba                         		
#SBATCH --output=~/logs/svaba/svaba%A_%a.out        
#SBATCH --error=~/logs/svaba/svaba%A_%a.out                   
#SBATCH --partition=day                              
#SBATCH --time=23:59:59                     	 
#SBATCH --mem=62g                              
#SBATCH --nodes=1                 		   
#SBATCH --cpus-per-task=12                 			
#SBATCH --mail-type=FAIL,END                     
#SBATCH --array=1-102%5 --requeue --export=all                          

## --- DO NOT EDIT from below here---- ##
## following will always overwrite previous output file, if any.
set +o noclobber

if [[ -z "${SLURM_ARRAY_TASK_ID}" ]]; then
    echo -e "\nERROR: Empty SLURM_ARRAY_TASK_ID variable. This script must be run in array job format.\n" >&2
    exit 1
fi

## sleep for n seconds before running any command
FORCEWAIT=$(shuf -i 10-200 -n 1)
printf "\nWaiting for %s seconds before starting workflow\n" "${FORCEWAIT:-10}"
sleep "${FORCEWAIT}"

echo "BEGIN at $(date)"

## notify slack when job starts if env variable PINGSTARTSLACK is set to YES
STARTMSG="MYJOB ID: ${SLURM_ARRAY_JOB_ID}: ${SLURM_JOB_NAME} starting at ${PWD} on $(hostname) for ${USER}"
echo -e "\n${STARTMSG}\n"

#### plate-channel manifest in tsv format ####
## This manifest instructs sbatch script to run each line as a separate job
manifest=~/svaba_manifest.txt

############################### START USER CMDS ################################
## parse manifest ##
## Extract each of N lines and process it as one of 1:N arrays specified above under
## SBATCH --array argument.

module load miniconda
conda activate renv

## extract plate and channel specific cellprofiler workflow command from a manifest
cp_cmd="$(sed -n ${SLURM_ARRAY_TASK_ID}p "${manifest}")"
TUMORS=$(ls ~/"${cp_cmd}"-TP*-WGS-*.realn.mdup.bqsr.bam 2>/dev/null \
               ~/"${cp_cmd}"-R*-WGS-*.realn.mdup.bqsr.bam 2>/dev/null)
NORMAL=$(ls ~/"${cp_cmd}"-NB*-WGS-*.realn.mdup.bqsr.bam)
OUTDIR=~/results/svaba

for TUMOR in $TUMORS; do
    TID=$(basename "$TUMOR" .realn.mdup.bqsr.bam)

    echo "Running SvABA for tumor: $TID"

    svaba run \
        -t "$TUMOR" \
        -n "$NORMAL" \
        -p 12 \
        -D ~/b37/dbsnp_138.b37.vcf \
        -a "${OUTDIR}/${TID}" \
        -G ~/b37/human_g1k_v37_decoy.fasta \
        -B ~/encode/hg19-blacklist.v2.bed \
        -z
done
svaba_exit=$?
svaba_end_time="$(date +%Y%m%d_%H%M%S_%Z)"

if [ "${svaba_exit}" -eq 0 ]; then
    printf "%s\tsvaba completed OK to %s\t%s\t%s\n" "${cp_cmd}" "${SLURM_ARRAY_JOB_ID}" "${svaba_end_time}" "${svaba_exit}"
    echo "svaba completed successfully for ${cp_cmd}"
else
    printf "%s\tsvaba completed  ERROR to %s\t%s\t%s\n" "${cp_cmd}" "${SLURM_ARRAY_JOB_ID}" "${svaba_end_time}" "${svaba_exit}"
    echo "svaba with an id: ${cp_cmd} failed. Exiting now without deleting the source file." >&2
    exit 1
fi

################################ END USER CMDS #################################

exitstat="${svaba_exit}"
echo "END at $(date)"

end_time="$(date +%d%b%y_%H%M%S_%Z)"
printf "%s_%s\t%s\t%s\tEND\t%s\n" "${SLURM_ARRAY_JOB_ID}" "${SLURM_ARRAY_TASK_ID}" "${cp_cmd}" "${end_time}" "${exitstat}"

echo "exit status was ${exitstat}"
exit ${exitstat}
## END ##
