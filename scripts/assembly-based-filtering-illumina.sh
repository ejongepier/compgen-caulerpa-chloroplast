#!/bin/bash

#SBATCH --job-name=ass-filt
#SBATCH --output=./logs/%x-%u-%A-%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=1260
#SBATCH --mem=72G

#######################################################
# HELP                                                #
#######################################################

usage="

Filtering of quality filtered illumina reads based on an assembly.

sbatch scritps/$(basename "$0") [-h|-f|-r|-d|-p|-o]

where:
    -h  show this help text
    -f  path to illumina fwd reads
    -r 	path to	illumina rev reads
    -d  path to assembly database in fasta.gz format
    -p  prefix used to label output files (e.g. species abbreviation)
    -o  output directory

example:
sbatch scripts/assembly-based-filtering-illumina.sh -f fnwi/202201-compgen-caulerpa-jongepier/results/qc/qc/sickle/\${spp}_R1_sickle.fq.gz -r fnwi/202201-compgen-caulerpa-jongepier/results/qc/qc/sickle/\${spp}_R2_sickle.fq.gz -d fnwi/202203-compgen-chloroplast-jongepier/results/assembly/\$spp-chloroplast-getorganelle-doubled -p \${spp} -o fnwi/202203-compgen-chloroplast-jongepier/results/readfilt
"



while getopts ':h:f:r:d:p:o:' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    f) FWDPATH=${OPTARG}
       ;;
    r) REVPATH=${OPTARG}
       ;;
    d) DB=${OPTARG}
       ;;
    p) PREFIX=${OPTARG}
       ;;
    o) OUTDIR=${OPTARG}
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



#######################################################
# MAIN                                                #
#######################################################

echo `date`" $SLURM_JOB_NAME started on node $SLURM_NODEID using $SLURM_CPUS_ON_NODE cpus."
echo "Command: sbatch $(basename "$0") -i $INPUT -p $PREFIX -o $OUTDIR"
echo "-------------------------------------------------"


## ====================================================
## Environment
## ====================================================




## ====================================================
## Paths etc
## ====================================================

DATE=`date +"%Y%m%dT%H%M%S"`
TMP=/scratch/$USER/$DATE
export TMPDIR=$TMP
srun mkdir -p $TMP




## ====================================================
## Collect input data
## ====================================================

echo `date`"  Collecting input data..."
srun mc cp $FWDPATH $TMP/
srun mc cp $REVPATH $TMP/
srun mc find $(dirname $DB) --name "$(basename $DB).*.bt2" --exec "mc cp {} $TMP/"
echo `date`"  Collecting input data finished"
echo "-------------------------------------------------"




## ====================================================
## Run bowtie2
## ====================================================

echo `date`"  Running bowtie2..."

cmd="srun bowtie2 -x $TMP/$(basename $DB) -p $SLURM_CPUS_ON_NODE -X 1000 -1 $TMP/$(basename $FWDPATH) -2 $TMP/$(basename $REVPATH) -S $TMP/${PREFIX}-mappings.sam --quiet --al-conc-gz $TMP/${PREFIX}_R%_$(basename $DB)-mapped.fastq.gz"
echo "Command: $cmd"
eval $cmd

echo `date`"  Running bowtie2 finished."
echo "-------------------------------------------------"



## ====================================================
## Cleanup
## ====================================================

echo `date`"  moving output files to $OUTDIR..."

srun mc cp $TMP/${PREFIX}_R1_$(basename $DB)-mapped.fastq.gz $OUTDIR/
srun mc cp $TMP/${PREFIX}_R2_$(basename $DB)-mapped.fastq.gz $OUTDIR/
srun rm -fr $TMP

echo `date`"  moving output files to $OUTDIR finished"
echo "-------------------------------------------------"





