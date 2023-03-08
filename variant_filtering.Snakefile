configfile: "config.yaml"

from itertools import chain
import glob
from os.path import basename
from os.path import join

rule separate_by_chrom:
    input:
        "data/genotypes/genotypes.vcf.gz"
    output:
        "data/genotypes/{chrom}.genotypes.vcf.gz"
    threads:
        4
    shell:
        "vcftools --gzvcf {input} --chr {wildcards.chrom} --stdout --recode --recode-INFO-all | bgzip -c -@ {threads} > {output}"

rule missingness:
    input:
        "data/genotypes/{chrom}.genotypes.vcf.gz"
    output:
        "data/genotypes/{chrom}.genotypes.imiss"
    threads:
        4
    shell:
        "vcftools --gzvcf {input} --missing-indv --out data/genotypes/{wildcards.chrom}.genotypes"

rule index_by_chrom:
    input:
        "data/genotypes/{chrom}.genotypes.vcf.gz"
    output:
        "data/genotypes/{chrom}.genotypes.vcf.gz.tbi"
    threads:
        4
    shell:
        "tabix -p vcf {input}"

rule annotate:
    input:
        gz="data/genotypes/{chrom}.genotypes.vcf.gz",
        tbi="data/genotypes/{chrom}.genotypes.vcf.gz.tbi"
    params:
        gatk=config["gatk_jar"],
        java=config["java"],
        filters=" ".join(["-filter '{}' --filter-name '{}'".format(filter["expr"], filter["name"]) for filter in config["gatk_filters"]])
    output:
        "data/genotypes/{chrom}.genotypes.annotated.vcf.gz"
    threads:
        4
    shell:
        "{params.java} -jar {params.gatk} VariantFiltration -V {input.gz} {params.filters} -O {output}"

rule filter:
    input:
        "data/genotypes/{chrom}.genotypes.annotated.vcf.gz"
    params:
    	min_meanDP=config["min_meanDP"],
    	max_meanDP=config["max_meanDP"],
    	max_missing_count=config["max_missing_count"],
    	maf=config["maf"]
    output:
        "data/genotypes/{chrom}.genotypes.filtered.vcf.gz"
    threads:
        4
    shell:
        "vcftools --gzvcf {input} --min-meanDP {params.min_meanDP} --max-meanDP {params.max_meanDP} --max-missing-count {params.max_missing_count} --maf {params.maf} --remove-filtered-all --stdout --recode --recode-INFO-all | bgzip -c -@ {threads} > {output}"

rule index_filtered:
    input:
        "data/genotypes/{chrom}.genotypes.filtered.vcf.gz"
    output:
        "data/genotypes/{chrom}.genotypes.filtered.vcf.gz.tbi"
    threads:
        4
    shell:
        "tabix -p vcf {input}"

rule select_snps:
    input:
        gz="data/genotypes/{chrom}.genotypes.filtered.vcf.gz",
        tbi="data/genotypes/{chrom}.genotypes.filtered.vcf.gz.tbi"
    params:
        gatk=config["gatk_jar"],
        java=config["java"]
    output:
        "data/genotypes/{chrom}.snps.filtered.vcf.gz"
    threads:
        3
    shell:
        """
        {params.java} -jar {params.gatk} SelectVariants -select-type SNP \
        -V {input.gz} -O {output}
        """

rule biallelic_snps:
    input:
        "data/genotypes/{chrom}.snps.filtered.vcf.gz"
    output:
        "data/genotypes/{chrom}.snps.biallelic.vcf.gz"
    threads:
        4
    shell:
        "vcftools --gzvcf {input} --min-alleles 2 --max-alleles 2 --stdout --recode --recode-INFO-all | bgzip -c -@ {threads} > {output}"
        
rule run_pipeline:
    input:
        filtered_snps=expand("data/genotypes/{chrom}.snps.biallelic.vcf.gz",
        	  	     chrom=config["chromosomes"]),
        missingness=expand("data/genotypes/{chrom}.genotypes.imiss",
                           chrom=config["chromosomes"])
