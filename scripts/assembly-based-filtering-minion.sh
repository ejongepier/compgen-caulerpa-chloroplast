#!/bin/bash 

#SBATCH --job-name=ass-filt
#SBATCH --output=./logs/%x-%u-%A-%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=4000
#SBATCH --mem=100G

#######################################################
# HELP                                                #
#######################################################

usage="

Assembly-based filtering of quality filtered minion reads.

sbatch scritps/$(basename "$0") [-h|-i|-d|-p|-o]

where:
    -h  show this help text
    -i  path to minion reads in fastq.gz format
    -d  path to assembly database in fna format
    -p  prefix used to label output files (e.g. species abbreviation)
    -o  output directory

example:
sbatch scripts/assembly-based-filtering-minion.sh -i fnwi/202011-metagenomics-caulerpa-anastasiabarilo/data/wgs/\${spp}-minion-wgs-filtered.fastq.gz -d  fnwi/202203-compgen-chloroplast-jongepier/results/assembly/\$spp-chloroplast-getorganelle-doubled -p \${spp} -o fnwi/202203-compgen-chloroplast-jongepier/results/readfilt
"


while getopts ':h:i:d:p:o:' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    i) INPUT=${OPTARG}
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

source ~/personal/miniconda3/etc/profile.d/conda.sh
conda activate /zfs/omics/projects/caulerpa/202201-compgen-caulerpa-jongepier/envs/qc-caulerpa




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
srun mc cp $INPUT $TMP/
srun mc cp ${DB}.fa.gz $TMP/
echo `date`"  Collecting input data finished"
echo "-------------------------------------------------"




## ====================================================
## Run minimap
## ====================================================

echo `date`"  Running minimap..."

cmd="srun minimap2 -ax map-ont -t $SLURM_CPUS_ON_NODE $TMP/$(basename $DB).fa.gz $TMP/$(basename $INPUT) | samtools view -bS -@ $SLURM_CPUS_ON_NODE -o $TMP/${PREFIX}-mappings.bam"
echo "Command: $cmd"
eval $cmd

echo `date`"  Running minimap finished."
echo "-------------------------------------------------"


echo `date`"  Get reads that mapped to the assembly..."

#Use samtools -F 4 to extract all mapped reads for single end
cmd="srun samtools fastq -F 4 -@ $SLURM_CPUS_ON_NODE $TMP/${PREFIX}-mappings.bam | gzip -c - > $TMP/${PREFIX}-minion-wgs-filtered-$(basename $DB)-mapped.fastq.gz"
echo "Command: $cmd"
eval $cmd

echo `date`"  Get reads that mapped to the assembly finished"
echo "-------------------------------------------------"



## ====================================================
## Cleanup
## ====================================================

echo `date`"  moving output files to $OUTDIR..."

srun mc cp $TMP/${PREFIX}-minion-wgs-filtered-$(basename $DB)-mapped.fastq.gz $OUTDIR/
srun rm -fr $TMP

echo `date`"  moving output files to $OUTDIR finished"
echo "-------------------------------------------------"

