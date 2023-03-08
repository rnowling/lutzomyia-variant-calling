configfile: "config.yaml"

from itertools import chain
import glob
from os.path import basename
from os.path import join


def get_samples(bam_paths):
    samples = set()
    
    for path in bam_paths:
        flname = basename(path)
        sample_end = flname.find(".")
        sample_name = flname[:sample_end]
        samples.add(sample_name)

    return list(samples)

bam_paths = glob.glob("data/deduped_alignments/*.deduped.bam")
all_samples = get_samples(bam_paths)

rule index_genome:
    input:
        "data/genome/genome.fa"
    output:
        "data/genome/genome.fa.bwt",
    threads:
        1
    shell:
        "bwa index {input}"
        
rule index_deduped_alignments:
    input:
        "data/deduped_alignments/{sample}.deduped.bam"
    output:
        "data/deduped_alignments/{sample}.deduped.bai"
    threads:
        6
    shell:
        "samtools index -b {input} {output}"

rule add_read_groups:
    input:
        bam="data/deduped_alignments/{sample}.deduped.bam",
        bai="data/deduped_alignments/{sample}.deduped.bai"
    output:
        bam="data/alignments_with_read_groups/{sample}.read_groups.bam"
    threads:
        6
    shell:
        "PicardCommandLine AddOrReplaceReadGroups I={input.bam} O={output.bam} RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM={wildcards.sample}"

rule genome_dict:
    input:
        fa="data/genome/genome.fa"
    params:
        gatk=config["gatk_jar"],
        java=config["java"]
    output:
        dict_path="data/genome/genome.dict"
    threads:
        4
    shell:
        "{params.java} -jar {params.gatk} CreateSequenceDictionary -R {input.fa} -O {output.dict_path}"

rule genome_index:
    input:
        fa="data/genome/genome.fa"
    output:
        index="data/genome/genome.fa.fai"
    threads:
        4
    shell:
        "samtools faidx {input.fa}"

rule index_read_groups:
    input:
        "data/alignments_with_read_groups/{sample}.read_groups.bam"
    output:
        "data/alignments_with_read_groups/{sample}.read_groups.bai"
    threads:
        6
    shell:
        "samtools index -b {input} {output}"

rule call_haplotypes:
    input:
        bam="data/alignments_with_read_groups/{sample}.read_groups.bam",
        bai="data/alignments_with_read_groups/{sample}.read_groups.bai",
        genome="data/genome/genome.fa",
        genome_idx="data/genome/genome.fa.fai",
        genome_dict="data/genome/genome.dict"
    params:
        gatk=config["gatk_jar"],
        java=config["java"],
        hetero=config["heterozygosity"],
        indel_hetero=config["indel_heterozygosity"],
        min_base_quality_score=config["min_base_quality_score"],
        contamination_fraction_to_filter=config["contamination_fraction_to_filter"]
    output:
        g_vcf="data/haplotypes/{sample}.g.vcf.gz"
    threads:
        4
    shell:
        "{params.java} -jar {params.gatk} HaplotypeCaller -R {input.genome} -I {input.bam} -O {output.g_vcf} -ERC GVCF --min-base-quality-score {params.min_base_quality_score} --contamination-fraction-to-filter {params.contamination_fraction_to_filter} --native-pair-hmm-threads {threads} --heterozygosity {params.hetero} --indel-heterozygosity {params.indel_hetero}"

rule run_pipeline:
    input:
        gvcfs=expand("data/haplotypes/{sample}.g.vcf.gz",
                     sample=all_samples)
