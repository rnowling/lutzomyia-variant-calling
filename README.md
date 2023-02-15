# lutzomyia-variant-calling
Pipeline for calling variants for Lutzoymia samples

The pipeline is organized into four stages:

* Alignment: trims reads, aligns reads to genome
* Haplotyping: Runs the GATK HaplotypeCaller on each BAM file
* Genotyping: Combines the resulting GVCFs and performs genotype calling
* Variant Filtering: Performs quality filtering, separates SNPs from indels

The pipeline expects to find a genome FASTA file and pairs of FASTQ files in the `input` directory.  The name of the genome file can be set in the `config.yaml` file.
The other files are found by searching for pairs of `*_1.fastq.gz.` and `*_2.fastq.gz` in the `input/samples` directory.

Parameters for trimming and haplotyping can be set in the `config.yaml` file.  The variant filtering parameters are currently hard-coded in `variant_filtering.Snakefile`, but they will be moved to `config.yaml` in the future. Parameters were taken primarily from the [1000 Anopheles genomes supplemental materials](https://www.nature.com/articles/nature24995#Sec13).

You can run the pipeline stages as follows:

```bash
$ snakemake --cores X --snakefile alignments.Snakefile run_pipeline
$ snakemake --cores X --snakefile haplotyping.Snakefile run_pipeline
$ snakemake --cores X --snakefile genotyping.Snakefile run_pipeline
$ snakemake --cores X --snakefile variant_filtering.Snakefile run_pipeline
```
where `X` is the number of cores to use.  The pipeline is set up to run on a single, high-memory machine.

Note that GATK version 4 will be needed.  GATK 4 depends on JVM version 8, so that might also need to be installed.
