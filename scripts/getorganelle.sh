#!/bin/bash

#SBATCH --job-name=gorg
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

Run getorganelle assembly using illumina pair end reads.

sbatch $(basename "$0") [-h|-f|-r|-o|-d|-p]
where:
    -h  show this help text
    -f  path to illumina sequence_1 file in fastq.gz format
    -r  path to illumina sequence_2 file in fastq.gz format
    -o  output directory
    -d  database e.g. embplant_pt
    -p  prefix

example:
sbatch scripts/getorganelle.sh -f fnwi/202201-compgen-caulerpa-jongepier/results/qc/qc/sickle/\${spp}_R1_sickle.fq.gz -r fnwi/202201-compgen-caulerpa-jongepier/results/qc/qc/sickle/\${spp}_R2_sickle.fq.gz -o fnwi/202203-compgen-chloroplast-jongepier/results/assembly -d embplant_pt -p \${spp}-chloroplast
sbatch scripts/getorganelle.sh -f fnwi/202201-compgen-caulerpa-jongepier/results/qc/qc/sickle/\${spp}_R1_sickle.fq.gz -r fnwi/202201-compgen-caulerpa-jongepier/results/qc/qc/sickle/\${spp}_R2_sickle.fq.gz -o fnwi/202201-compgen-caulerpa-jongepier/results/getorganelle -d embplant_mt -p \${spp}-mitochondrium
"



while getopts ':h:f:r:o:d:p:' option; do
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
    d) DB="${OPTARG}"
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

source ~/personal/miniconda3/etc/profile.d/conda.sh
conda activate /zfs/omics/projects/caulerpa/202201-compgen-caulerpa-jongepier/envs/getorganelle



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

echo `date`" Collecting input data finshed"
## -------------------------------------------------



## =================================================
## Run getorganelle
## =================================================


echo `date`"  Getorganelle  started..."

mkdir -p $TMP/$PREFIX
#cmd="srun get_organelle_from_reads.py -1 $TMP/$(basename $FWDPATH) -2 $TMP/$(basename $REVPATH) -o $TMP/$PREFIX -t $SLURM_CPUS_ON_NODE -F $DB -R 15 -k 21,45,65,85,105 --overwrite"
cmd="srun get_organelle_from_reads.py -1 $TMP/$(basename $FWDPATH) -2 $TMP/$(basename $REVPATH) -o $TMP/$PREFIX -t $SLURM_CPUS_ON_NODE -F $DB -R 30 -k 21,45,65,85,105 --overwrite --reduce-reads-for-coverage inf --max-reads inf"
echo "Command: $cmd"
eval $cmd

echo `date`" Getorganelle finished"
## -------------------------------------------------



## =================================================
## Cleanup
## =================================================


echo `date`"  Archiving and moving data started..."

srun mv $TMP/${PREFIX} $TMP/${PREFIX}-getorganelle-run
srun tar --exclude='${PREFIX}-getorganelle-run/extended_*paired.fq' -czvf ${PREFIX}-getorganelle-run.tar.gz -C $TMP ${PREFIX}-getorganelle-run 
srun mc cp $TMP/${PREFIX}-getorganelle-run.tar.gz $OUTDIR/

echo `date`"  Archiving and moving data finished."
## -------------------------------------------------


echo `date`" All done!"
echo "$SLURM_JOB_NAME finished on node $SLURM_NODEID using $SLURM_CPUS_ON_NODE cpus."


####################################################
# THE END                                          #
####################################################



