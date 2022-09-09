#!/bin/bash

#SBATCH --job-name=masurca
#SBATCH --output=logs/%x-%u-%A-%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=12000
#SBATCH --mem=360G
##SBATCH --nodelist=omics-cn003

####################################################
# HELP                                             #
####################################################

usage="

Run masurca assembly using organellar illumina pair end reads and minion reads.

sbatch $(basename "$0") [-h|-f|-r|-l|-o]
where:
    -h  show this help text
    -f  path to illumina sequence_1 file in fastq.gz format
    -r  path to illumina sequence_2 file in fastq.gz format
    -l  long read path
    -o  output directory

example:
sbatch scripts/masurca.sh -f fnwi/202201-compgen-caulerpa-jongepier/results/organel-filt/\${spp}_R1_organel.fastq.gz -r fnwi/202201-compgen-caulerpa-jongepier/results/organel-filt/\${spp}_R2_organel.fastq.gz -l fnwi/202201-compgen-caulerpa-jongepier/results/organel-filt/\${spp}-minion-wgs-filtered-organel.fastq.gz -o fnwi/202203-compgen-chloroplast-jongepier/results/masurca/\${spp}
"



while getopts ':h:f:r:l:o:' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    f) FWDPATH="${OPTARG}"
       ;;
    r) REVPATH="${OPTARG}"
       ;;
    o) OUTDIR="${OPTARG}"
       ;;
    l) LRPATH="${OPTARG}"
       ;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
  esac
done
shift $((OPTIND - 1))


####################################################
# MAIN                                             #
####################################################

echo `date`" $SLURM_JOB_NAME started on node $SLURM_NODEID using $SLURM_CPUS_ON_NODE cpus."
echo "Command: sbatch $(basename "$0") -f $FWDPATH -r $REVPATH -o $OUTDIR -d $DB -p $PREFIX"



## =================================================
## Environment
## =================================================

#conda activate /zfs/omics/projects/caulerpa/202201-compgen-caulerpa-jongepier/envs/masurca

## =================================================
## Paths on scratch
## =================================================

DATE=`date +"%Y%m%dT%H%M%S"`
TMP=/scratch/$USER/$DATE
export TMPDIR=$TMP
srun mkdir -p $TMP



## =================================================
## Collect input data
## =================================================

echo `date`" Collecting input data"

srun mc cp $FWDPATH $TMP/
srun mc cp $REVPATH $TMP/
srun mc cp $LRPATH $TMP/

echo `date`" Collecting input data finshed"
## -------------------------------------------------



## =================================================
## Run masurca
## =================================================


echo `date`"  Masurca  started..."

cmd="srun masurca -i $TMP/$(basename $FWDPATH),$TMP/$(basename $REVPATH) -r $TMP/$(basename $LRPATH) -t $SLURM_CPUS_ON_NODE -o $TMP"
echo "Command: $cmd"
eval $cmd

echo `date`" Masurca finished"
## -------------------------------------------------



## =================================================
## Cleanup
## =================================================


echo `date`" All done!"
echo "$SLURM_JOB_NAME finished on node $SLURM_NODEID using $SLURM_CPUS_ON_NODE cpus."


####################################################
# THE END                                          #
####################################################



