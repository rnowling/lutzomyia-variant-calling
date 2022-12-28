configfile: "config.yaml"

from itertools import chain
import glob
from os.path import basename
from os.path import join

rule select_snps:
    input:
        "data/genotypes/genotypes.vcf.gz"
    params:
        gatk=config["gatk_jar"],
        java=config["java"]
    output:
        "data/genotypes/snps.vcf.gz"
    threads:
        3
    shell:
        """
        {params.java} -jar {params.gatk} SelectVariants -select-type SNP \
        -V {input} -O {output}
        """

rule annotate_snps:
    input:
        "data/genotypes/snps.vcf.gz"
    params:
        gatk=config["gatk_jar"],
        java=config["java"]
    output:
        "data/genotypes/snps.annotated.vcf.gz"
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

rule filter_snps:
    input:
        "data/genotypes/snps.annotated.vcf.gz"
    output:
        "data/genotypes/snps.filtered.vcf.gz"
    threads:
        1
    shell:
        "vcftools --gzvcf {input} --min-meanDP 14 --max-meanDP 60 --max-missing-count 7 --remove-filtered-all --stdout --recode --recode-INFO-all --mac 2 | gzip -c > {output}" 

rule biallelic_snps:
    input:
        "data/genotypes/snps.filtered.vcf.gz"
    output:
        "data/genotypes/snps.filtered.biallelic.vcf.gz"
    threads:
        1
    shell:
        "vcftools --gzvcf {input} --min-alleles 2 --max-alleles 2 --stdout --recode --recode-INFO-all | gzip -c > {output}"
        
rule run_pipeline:
    input:
        filtered_snps="data/genotypes/snps.filtered.biallelic.vcf.gz"
