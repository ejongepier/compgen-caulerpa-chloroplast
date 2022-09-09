#!/bin/bash 

#SBATCH --job-name=strainberry
#SBATCH --output=./logs/%x-%u-%A-%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=4000
#SBATCH --mem=100G
#SBATCH --nodelist=omics-cn004

#######################################################
# HELP                                                #
#######################################################

usage="

Strain-leel haplotype phased assembly using stainberry.

sbatch scritps/$(basename "$0") [-h|-i|-a|-p|-o]

where:
    -h  show this help text
    -i  path to minion reads in fastq.gz format
    -a  path to assembly in fna format
    -p  prefix used to label output files (e.g. species abbreviation)
    -o  output directory

example:
sbatch scripts/strainberry.sh -i fnwi/202011-metagenomics-caulerpa-anastasiabarilo/data/wgs/\${spp}-minion-wgs-filtered.fastq.gz -a fnwi/202201-compgen-caulerpa-jongepier/results/getorganelle/cpro-chloroplast-genome.fa.gz -p \${spp} -o fnwi/202203-compgen-chloroplast-jongepier/results/assembly/
"


while getopts ':h:i:a:p:o:' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    i) INPUT=${OPTARG}
       ;;
    a) ASSEMBLY=${OPTARG}
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
echo "Command: sbatch $(basename "$0") -i $INPUT -a $ASSEMBLY -p $PREFIX -o $OUTDIR"
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
#srun mc cp $INPUT $TMP/
#srun mc cp $ASSEMBLY $TMP/
echo `date`"  Collecting input data finished"
echo "-------------------------------------------------"




## ====================================================
## Run minimap
## ====================================================

echo `date`"  Running minimap..."

#cmd="srun minimap2 -ax map-ont -t $SLURM_CPUS_ON_NODE $TMP/$(basename $ASSEMBLY) $TMP/$(basename $INPUT) | samtools sort -@ $SLURM_CPUS_ON_NODE > $TMP/${PREFIX}-mappings.bam"
#echo "Command: $cmd"
#eval $cmd

#gzip -d $TMP/$(basename $ASSEMBLY)
ASS=$(basename $ASSEMBLY)
ASSIN=${ASS/.gz/}
#cmd="srun samtools faidx $TMP/$ASSIN"
#echo "Command: $cmd"
#eval $cmd


#cmd="srun samtools index $TMP/${PREFIX}-mappings.bam"
#echo "Command: $cmd"
#eval $cmd

echo `date`"  Running minimap finished."
echo "-------------------------------------------------"



## ====================================================
## Run strainberry
## ====================================================

echo `date`"  Running strainberry..."

source ~/personal/miniconda3/etc/profile.d/conda.sh
conda activate /zfs/omics/projects/caulerpa/202203-compgen-chloroplast-jongepier/envs/strainberry

cmd="srun strainberry -r $TMP/$ASSIN -b $TMP/${PREFIX}-mappings.bam -o $TMP/${PREFIX}-strainberry-out -c $SLURM_CPUS_ON_NODE --nanopore --max-strains 100"
echo "Command: $cmd"
eval $cmd

echo `date`"  Running strainberry finished"
echo "-------------------------------------------------"



## ====================================================
## Cleanup
## ====================================================

#echo `date`"  moving output files to $OUTDIR..."

#srun mc cp $OUT/${PREFIX}-minion-wgs-filtered-$(basename $DB)-mapped.fastq.gz $OUTDIR/
#srun mc cp $OUT/${PREFIX}-minion-wgs-filtered-$(basename $DB)-mapping.stats $OUTDIR/
#srun rm -R $OUT
#srun rm -fr $TMP

echo `date`"  moving output files to $OUTDIR finished"
echo "-------------------------------------------------"

