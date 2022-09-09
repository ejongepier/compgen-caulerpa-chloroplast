#!/bin/bash

#SBATCH --job-name=depth-pileup
#SBATCH --output=./logs/%x-%u-%A-%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=1260
#SBATCH --mem=72G

#######################################################
# HELP                                                #
#######################################################

usage="

Compute sequencing depth and pileup of illumina reads based on an assembly.

sbatch scritps/$(basename "$0") [-h|-f|-r|-d|-p|-o]

where:
    -h  show this help text
    -f  path to illumina fwd reads
    -r 	path to	illumina rev reads
    -d  path to assembly database in fasta.gz format
    -p  prefix used to label output files (e.g. species abbreviation)
    -o  output directory

example:
sbatch scripts/depth-pileup-illumina.sh -f fnwi/202201-compgen-caulerpa-jongepier/results/assembly-based-filt/clen_R1_clen-chloroplast-noncirc-genome-mapped.fastq.gz -r fnwi/202201-compgen-caulerpa-jongepier/results/assembly-based-filt/clen_R2_clen-chloroplast-noncirc-genome-mapped.fastq.gz -d fnwi/202201-compgen-caulerpa-jongepier/results/masurca/clen-masurca-organel/clen-chloroplast-genome -p clen -o fnwi/202203-compgen-chloroplast-jongepier/results/pileup
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

#conda activate /zfs/omics/projects/caulerpa/202201-compgen-caulerpa-jongepier/envs/qc-caulerpa



## ====================================================
## Paths etc
## ====================================================

DATE=`date +"%Y%m%dT%H%M%S"`
TMP=/scratch/$USER/tmp-${DATE}
export TMPDIR=$TMP
srun mkdir -p $TMP




## ====================================================
## Collect input data
## ====================================================

echo `date`"  Collecting input data..."
srun mc cp $FWDPATH $TMP
srun mc cp $REVPATH $TMP
srun mc find $(dirname $DB) --name "$(basename $DB).*.bt2" --exec "mc cp {} $TMP/"
echo `date`"  Collecting input data finished"
echo "-------------------------------------------------"




## ====================================================
## Run bowtie2
## ====================================================

echo `date`"  Running bowtie2..."

cmd="srun bowtie2 -x $TMP/$(basename $DB) -p $SLURM_CPUS_ON_NODE -X 1000 -p $SLURM_CPUS_ON_NODE -1 $TMP/$(basename $FWDPATH) -2 $TMP/$(basename $REVPATH) -S $TMP/${PREFIX}-$(basename $DB)-mappings.sam --quiet"
echo "Command: $cmd"
eval $cmd

echo `date`"  Running bowtie2 finished."
echo "-------------------------------------------------"



## ====================================================
## Run sam to sorted bam
## ====================================================

echo `date`"  Running samtools sort..."

cmd="srun samtools view -@ $SLURM_CPUS_ON_NODE -S -b $TMP/${PREFIX}-$(basename $DB)-mappings.sam | samtools sort - -@ $SLURM_CPUS_ON_NODE -o $TMP/${PREFIX}-$(basename $DB)-mappings.bam"
echo "Command: $cmd"
eval $cmd

echo `date`"  Running samtools sort finished"
echo "-------------------------------------------------"



## ====================================================
## Computing coverage
## ====================================================

echo `date`"  Running samtools depth..."

## samtools depth automatically skips secondary alignments
cmd="srun samtools depth -d 1000000 -a $TMP/${PREFIX}-$(basename $DB)-mappings.bam >  $TMP/${PREFIX}-$(basename $DB)-mappings.depth"                               
echo "Command: $cmd"
eval $cmd

echo `date`"  Running samtools depth finished."
echo "-------------------------------------------------"



## ====================================================
## Run samtools pileup
## ====================================================

echo `date`"  Running samtools pileup..."

cmd="srun samtools mpileup -d 1000000 -o $TMP/${PREFIX}-$(basename $DB)-mappings.pileup $TMP/${PREFIX}-$(basename $DB)-mappings.bam"
echo "Command: $cmd"
eval $cmd

echo `date`"  Running samtools pileup finished."
echo "-------------------------------------------------"




## ====================================================
## Cleanup
## ====================================================

echo `date`"  moving output files to $OUTDIR..."

srun mc cp $TMP/${PREFIX}-$(basename $DB)-mappings.bam $OUTDIR/
srun mc cp $TMP/${PREFIX}-$(basename $DB)-mappings.depth $OUTDIR/
srun mc cp $TMP/${PREFIX}-$(basename $DB)-mappings.pileup $OUTDIR/
srun rm -fr $TMP

echo `date`"  moving output files to $OUTDIR finished"
echo "-------------------------------------------------"





