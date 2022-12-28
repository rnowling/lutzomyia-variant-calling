configfile: "config.yaml"

from itertools import chain
import glob
from os.path import basename
from os.path import join


def get_samples():
    samples = set()
    
    for path in glob.glob("data/haplotypes/*.g.vcf.gz"):
        flname = basename(path)
        sample_end = flname.find(".")
        sample_name = flname[:sample_end]
        samples.add(sample_name)

    return list(samples)

all_samples = get_samples()

rule combine_gvcfs:
    input:
        genome="data/genome/genome.fa",
        gvcfs=lambda w: expand("data/haplotypes/{sample}.g.vcf.gz", sample=all_samples)
    params:
        gatk=config["gatk_jar"],
        java=config["java"]
    output:
        "data/genotypes/combined.g.vcf.gz"
    threads:
        3
    run:
        cmd = "{params.java} -jar {params.gatk} CombineGVCFs -R {input.genome} "
        cmd += "--variant " + " --variant ".join(input.gvcfs)
        cmd += " -O {output}"
        shell(cmd)

rule genotype_vcfs:
    input:
        genome="data/genome/genome.fa",
        gvcf="data/genotypes/combined.g.vcf.gz"
    params:
        gatk=config["gatk_jar"],
        java=config["java"]
    output:
        "data/genotypes/genotypes.vcf.gz"
    threads:
        3
    shell:
        "{params.java} -jar {params.gatk} GenotypeGVCFs -R {input.genome} -V {input.gvcf} -O {output}"

rule run_pipeline:
    input:
        vcfs="data/genotypes/genotypes.vcf.gz"
