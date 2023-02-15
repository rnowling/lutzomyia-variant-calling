configfile: "config.yaml"

from itertools import chain
import glob
from os.path import basename
from os.path import join

rule annotate:
    input:
        "data/genotypes/genotypes.vcf.gz"
    params:
        gatk=config["gatk_jar"],
        java=config["java"]
    output:
        "data/genotypes/genotypes.annotated.vcf.gz"
    threads:
        4
    shell:
        """
        {params.java} -jar {params.gatk} VariantFiltration -V {input} \
        -filter "QD < 5.0" --filter-name "QD5" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "FS > 60.0" --filter-name "FS60" \
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
        -O {output}
        """
rule filter:
    input:
        "data/genotypes/genotypes.annotated.vcf.gz"
    output:
        gz="data/genotypes/genotypes.filtered.vcf.gz"
    threads:
        4
    shell:
        "vcftools --gzvcf {input} --min-meanDP 14 --max-meanDP 120 --max-missing-count 7 --maf 0.05 --remove-filtered-all --stdout --recode --recode-INFO-all | bgzip -c -@ {threads} > {output.gz}"

rule index:
    input:
        "data/genotypes/genotypes.filtered.vcf.gz"
    output:
        "data/genotypes/genotypes.filtered.vcf.gz.tbi"
    threads:
        4
    shell:
        "tabix -p vcf {input}"

rule select_snps:
    input:
        gz="data/genotypes/genotypes.filtered.vcf.gz",
        tbi="data/genotypes/genotypes.filtered.vcf.gz.tbi"
    params:
        gatk=config["gatk_jar"],
        java=config["java"]
    output:
        "data/genotypes/snps.filtered.vcf.gz"
    threads:
        3
    shell:
        """
        {params.java} -jar {params.gatk} SelectVariants -select-type SNP \
        -V {input.gz} -O {output}
        """

rule biallelic_snps:
    input:
        "data/genotypes/snps.filtered.vcf.gz"
    output:
        "data/genotypes/snps.filtered.biallelic.vcf.gz"
    threads:
        4
    shell:
        "vcftools --gzvcf {input} --min-alleles 2 --max-alleles 2 --stdout --recode --recode-INFO-all | bgzip -c -@ {threads} > {output}"
        
rule run_pipeline:
    input:
        filtered_snps="data/genotypes/snps.filtered.biallelic.vcf.gz"
