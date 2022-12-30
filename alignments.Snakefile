configfile: "config.yaml"

from itertools import chain
import glob
from os.path import basename
from os.path import join


def get_samples(fastq_paths):
    samples = set()
    
    for path in fastq_paths:
        flname = basename(path)
        sample_end = flname.find("_")
        sample_name = flname[:sample_end]
        samples.add(sample_name)

    return list(samples)

fastq_paths = glob.glob(join(config["samples_dir"], "*.fastq.gz"))
all_samples = get_samples(fastq_paths)

# copy appropriate files to data directory
rule link_reads:
    input:
        src=lambda wildcards: "%s/{filename}.fastq.gz" % config["samples_dir"]
    output:
        "data/raw_reads/{filename}.fastq.gz"
    threads:
        1
    shell:
        "ln -s ../../{input.src} {output}"

rule copy_genome:
    input:
        config["genome_path"]
    output:
        "data/genome/genome.fa"
    threads:
        1
    shell:
        "cp {input} {output}"

rule clean_reads:
    input:
        reads1="data/raw_reads/{sample}_1.fastq.gz",
        reads2="data/raw_reads/{sample}_2.fastq.gz",
        adapters=config["trimmomatic_adapters"]
    params:
        log=lambda w: "data/cleaned_reads/{}.log".format(w.sample)
    output:
        "data/cleaned_reads/{sample}_1p.fastq.gz",
        "data/cleaned_reads/{sample}_1u.fastq.gz",
        "data/cleaned_reads/{sample}_2p.fastq.gz",
        "data/cleaned_reads/{sample}_2u.fastq.gz"
    threads:
        4
    shell:
        "TrimmomaticPE -threads {threads} {input.reads1} {input.reads2} {output} ILLUMINACLIP:{input.adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 HEADCROP:10 MINLEN:%(trimmomatic_minlen)s &> {params.log}" % config

rule index_genome:
    input:
        "data/genome/genome.fa"
    output:
        "data/genome/genome.fa.bwt",
    threads:
        1
    shell:
        "bwa index {input}"

rule bwa_mem:
    input:
        genome="data/genome/genome.fa",
        index="data/genome/genome.fa.bwt",
        fastq1="data/cleaned_reads/{sample}_1p.fastq.gz",
        fastq2="data/cleaned_reads/{sample}_2p.fastq.gz"
    output:
        "data/alignments/{sample}.bam"
    threads:
        4
    shell:
        "bwa mem -t {threads} {input.genome} {input.fastq1} {input.fastq2} | samtools sort -@ {threads} -O bam -o {output}"

rule mapping_stats:
    input:
        "data/alignments/{sample}.bam"
    output:
        "data/alignments/{sample}_flagstat.txt"
    threads:
        2
    shell:
        "samtools flagstat {input} > {output}"
        
# filter value (-F 260) came directly from Krystal's scripts
rule filter_bam:
    input:
        bam="data/alignments/{sample}.bam"
    output:
        "data/filtered_alignments/{sample}.filtered.bam"
    threads:
        6
    shell:
        "samtools view -b -@{threads} -F 260 {input.bam} > {output}"

rule index_filtered_alignments:
    input:
        "data/filtered_alignments/{sample}.filtered.bam"
    output:
        "data/filtered_alignments/{sample}.filtered.bai"
    threads:
        6
    shell:
        "samtools index -b {input} {output}"

rule mark_duplicates:
    input:
        bam="data/filtered_alignments/{sample}.filtered.bam",
        index="data/filtered_alignments/{sample}.filtered.bai"
    output:
        bam="data/deduped_alignments/{sample}.deduped.bam",
        metrics="data/deduped_alignments/{sample}.deduped.metrics"
    threads:
        6
    shell:
        "PicardCommandLine MarkDuplicates I={input.bam} O={output.bam} M={output.metrics} REMOVE_DUPLICATES=true"
        
rule cov_stats:
    input:
        "data/deduped_alignments/{sample}.deduped.bam"
    output:
        "data/deduped_alignments/{sample}_cov.txt"
    threads:
        2
    shell:
        "samtools coverage {input} > {output}"

rule run_pipeline:
    input:
        mapping_stats=expand("data/alignments/{sample}_flagstat.txt",
                             sample=all_samples),
        cov_stats=expand("data/deduped_alignments/{sample}_cov.txt",
                         sample=all_samples)
