# README to Caulerpa chloroplast comparative genomics project

This repo contains the code and documentation for the comparative genomics of the Caulerpa chloroplast genomes.
How these genomes were obtained is detailed elsewhere ().

## Usage

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
mc cp results/getorfs/${SPP}-chloroplast-orfs.faa fnwi/202203-compgen-chloroplast-jongepier/results/getorfs/
rm -r data results
```


### Annotate chloroplast genomes

Chloroplast genomes were annotated with the webutility of ``GeSeq``.



## Authors

* Evelien Jongepier (e.jongepier@uva.nl)

