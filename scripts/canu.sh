#!/bin/bash

#SBATCH --job-name=canu
#SBATCH --output=logs/%x-%u-%A-%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=12000
#SBATCH --mem=200G
##SBATCH --nodelist=omics-cn003

####################################################
# HELP                                             #
####################################################

usage="

Run canu assembly using organellar minion reads.

sbatch $(basename "$0") [-h|-l|-g|-o|-p]
where:
    -h  show this help text
    -l  long read path
    -g  genome size
    -o  output directory
    -p  prefix

example:
sbatch scripts/canu.sh -l fnwi/202201-compgen-caulerpa-jongepier/results/organel-filt/\${spp}-minion-wgs-filtered-organel.fastq.gz -g 345k -o fnwi/202203-compgen-chloroplast-jongepier/results/canu/\${spp} -p \${spp}-canu
"

while getopts ':h:l:g:o:p:' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    l) LRPATH="${OPTARG}"
       ;;
    g) GENOMESIZE="${OPTARG}"
       ;;
    o) OUTDIR="${OPTARG}"
       ;;
    p) PREFIX="${OPTARG}"
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

#conda activate /zfs/omics/projects/caulerpa/202201-compgen-caulerpa-jongepier/envs/canu

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

srun mc cp $LRPATH $TMP/

echo `date`" Collecting input data finshed"
## -------------------------------------------------



## =================================================
## Run canu
## =================================================


echo `date`"  Canu  started..."

cmd="srun canu -d $TMP -p $PREFIX genomeSize=$GENOMESIZE -nanopore-raw useGrid=false $TMP/$(basename $LRPATH)"
echo "Command: $cmd"
eval $cmd

echo `date`" Canu finished"
## -------------------------------------------------



## =================================================
## Cleanup
## =================================================


echo `date`" All done!"
echo "$SLURM_JOB_NAME finished on node $SLURM_NODEID using $SLURM_CPUS_ON_NODE cpus."


####################################################
# THE END                                          #
####################################################



