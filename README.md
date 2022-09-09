# README to Caulerpa chloroplast comparative genomics project

This repo contains the code and documentation for the comparative genomics of the Caulerpa chloroplast genomes.
How these genomes were obtained is detailed elsewhere (https://github.com/ejongepier/compgen-caulerpa-fin).
Results are stored in minio bucket fnwi/202203-compgen-chloroplast-jongepier/.

## Usage

### Preliminary chloroplast assembly

The purpose is to get a quick assembly to use for read mapping and filtering. The filtered reads will be used in the final assembly

```bash
sbatch scripts/getorganelle.sh -f fnwi/202201-compgen-caulerpa-jongepier/results/qc/qc/sickle/${spp}_R1_sickle.fq.gz -r fnwi/202201-compgen-caulerpa-jongepier/results/qc/qc/sickle/${spp}_R2_sickle.fq.gz -o fnwi/202203-compgen-chloroplast-jongepier/results/assembly -d embplant_pt -p ${spp}-chloroplast
```

Post process by changing the fasta header, doubling the sequence and building the bowtei2 db.
The purpose of doubling the seq is to allow for concordant mapping accross the circularized junction.


### Filter chloroplast reads

The purpose is to get the chloroplast reads as input to metaspades

```bash
sbatch scripts/assembly-based-filtering-illumina.sh -f fnwi/202201-compgen-caulerpa-jongepier/results/qc/qc/sickle/${spp}_R1_sickle.fq.gz -r fnwi/202201-compgen-caulerpa-jongepier/results/qc/qc/sickle/${spp}_R2_sickle.fq.gz -d fnwi/202203-compgen-chloroplast-jongepier/results/assembly/$spp-chloroplast-getorganelle-doubled -p ${spp} -o fnwi/202203-compgen-chloroplast-jongepier/results/readfilt
sbatch scripts/assembly-based-filtering-minion.sh -i fnwi/202011-metagenomics-caulerpa-anastasiabarilo/data/wgs/${spp}-minion-wgs-filtered.fastq.gz -d  fnwi/202203-compgen-chloroplast-jongepier/results/assembly/$spp-chloroplast-getorganelle-doubled -p ${spp} -o fnwi/202203-compgen-chloroplast-jongepier/results/readfilt
```

### Cloroplast meta-assembly



### Identify chloroplast ORFs

ORFs were identified with EMBOSS getorfs utility, and annotated using blastp against NCBI nr.
This analyses was run on blobfish.

Get chloroplast genome assemblies:

```bash
mkdir -p data/assemblies
mc cp fnwi/202201-compgen-caulerpa-jongepier/results/masurca/clen-masurca-organel/clen-chloroplast-genome.fa.gz data/assemblies/
mc cp fnwi/202201-compgen-caulerpa-jongepier/results/getorganelle/cpro-chloroplast-genome.fa.gz data/assemblies/
```

Get ORFs:

```bash
SPP=clen
mkdir -p results/getorfs/
gzip -d data/assemblies/${SPP}-chloroplast-genome.fa.gz
getorf -circular Y -reverse T -sequence data/assemblies/${SPP}-chloroplast-genome.fa -outseq results/getorfs/${SPP}-chloroplast-orfs.faa
```

Get NCBI db:

```bash
cd ~/Databases/ncbi-nr
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
gzip -d nr.gz
mv nr nr.faa
makeblastdb -dbtype 'prot' -in nr.faa
echo `date` > VERSION.txt
cd -
```

Blast orfs against nr:

```bash
blastp -db ~/Databases/ncbi-nr/nr -query results/getorfs/${SPP}-chloroplast-orfs.faa -outfmt 6
```

Cleanup

```bash
gzip results/getorfs/${SPP}-chloroplast-orfs.faa
mc cp results/getorfs/${SPP}-chloroplast-orfs.faa fnwi/202203-compgen-chloroplast-jongepier/results/getorfs/
rm -r data results
```


### Annotate chloroplast genomes

Chloroplast genomes were annotated with the webutility of ``GeSeq``.


## Poolseq variant calling

```bash 
sbatch scripts/depth-pileup-illumina.sh -f fnwi/202201-compgen-caulerpa-jongepier/results/assembly-based-filt/clen_R1_clen-chloroplast-noncirc-genome-mapped.fastq.gz -r fnwi/202201-compgen-caulerpa-jongepier/results/assembly-based-filt/clen_R2_clen-chloroplast-noncirc-genome-mapped.fastq.gz -d fnwi/202201-compgen-caulerpa-jongepier/results/masurca/clen-masurca-organel/clen-chloroplast-genome -p clen -o fnwi/202203-compgen-chloroplast-jongepier/results/pileup
```

## Authors

* Evelien Jongepier (e.jongepier@uva.nl)

