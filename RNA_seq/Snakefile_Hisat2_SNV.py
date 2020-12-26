# Rule 0: Build rule of all rules_
sample_bee = ['6107','6154','6246_A','7454','7538','9540_A','332014','338391']

rule all:
    input:
        expand("10_lofreq_hisat2_results/{sample}/{sample}_lofreq.vcf", sample=sample_bee),
        expand("10_lofreq_hisat2_results/{sample}/{sample}_lofreq.vcf.gz", sample=sample_bee),
        expand("10_lofreq_hisat2_results/{sample}/{sample}_lofreq.vcf.gz.tbi", sample=sample_bee),
        expand("10_lofreq_hisat2_results/{sample}/{sample}_lofreq_stat.vchk", sample=sample_bee),
        expand("10_quasitool_hisat2_results/{sample}/{sample}_quasitool.vcf", sample=sample_bee),
        expand("10_quasitool_hisat2_results/{sample}/{sample}_quasitool.vcf.gz", sample=sample_bee),
        expand("10_quasitool_hisat2_results/{sample}/{sample}_quasitool.vcf.gz.tbi", sample=sample_bee),
        expand("10_quasitool_hisat2_results/{sample}/{sample}_quasitool_stat.vchk", sample=sample_bee)

        #expand("CliqueSNV_output/{sample}/", sample=sample_bee)


# Rule 1: Calling variants by LoFreq
rule lofreq_call:
    input:
        "00_reference/BQCV.fasta",
        "06_hisat2_results/BQCV/{sample}/{sample}_sorted_mapped.bam"
    output:
        "10_lofreq_hisat2_results/{sample}/{sample}_lofreq.vcf"
    log:
        "log/hisat2_lofreq/{sample}_lofreq"
    shell:
        "lofreq call-parallel --pp-threads 8 -f {input[0]} {input[1]} -o {output} 2>{log}"

# Rule 2: Compress the VCF file
rule bgzip_lofreq:
    input:
        "10_lofreq_hisat2_results/{sample}/{sample}_lofreq.vcf"
    output:
        "10_lofreq_hisat2_results/{sample}/{sample}_lofreq.vcf.gz"
    shell:
        "bgzip -c {input} > {output}"

# Rule 3: Create an index file for the VCF file
rule bcftools_index_lofreq:
    input:
        "10_lofreq_hisat2_results/{sample}/{sample}_lofreq.vcf.gz"
    output:
        "10_lofreq_hisat2_results/{sample}/{sample}_lofreq.vcf.gz.tbi"
    shell:
        "bcftools index -t {input} -o {output}"

# Rule 4: Generate a text file with different stats from the VCF file
rule bcftools_stats_lofreq:
    input:
        "10_lofreq_hisat2_results/{sample}/{sample}_lofreq.vcf.gz"
    output:
        "10_lofreq_hisat2_results/{sample}/{sample}_lofreq_stat.vchk"
    shell:
        "bcftools stats {input} > {output}"


# Rule 5: Calling variants by quasitools
rule quasitools_ntvar:
    input:
        "00_reference/BQCV.fasta",
        "06_hisat2_results/BQCV/{sample}/{sample}_sorted_mapped.bam"
    output:
        "10_quasitool_hisat2_results/{sample}/{sample}_quasitool.vcf"
    shell:
        "quasitools call ntvar {input[1]} {input[0]} -o {output}"


# Rule 6: Compress the VCF file
rule bgzip_quasi:
    input:
        "10_quasitool_hisat2_results/{sample}/{sample}_quasitool.vcf"
    output:
        "10_quasitool_hisat2_results/{sample}/{sample}_quasitool.vcf.gz"
    shell:
        "bgzip -c {input} > {output}"

# Rule 7: Create an index file for the VCF file
rule bcftools_index_quasi:
    input:
        "10_quasitool_hisat2_results/{sample}/{sample}_quasitool.vcf.gz"
    output:
        "10_quasitool_hisat2_results/{sample}/{sample}_quasitool.vcf.gz.tbi"
    shell:
        "bcftools index -t {input} -o {output}"

# Rule 8: Generate a text file with different stats from the VCF file
rule bcftools_stats_quasi:
    input:
        "10_quasitool_hisat2_results/{sample}/{sample}_quasitool.vcf.gz"
    output:
        "10_quasitool_hisat2_results/{sample}/{sample}_quasitool_stat.vchk"
    shell:
        "bcftools stats {input} > {output}"


# Rule2: Quasispecies reconstruction
#rule clique_snv:
#    input:
#        "06_hisat2_results/{sample}/{sample}_sorted_mapped.bam"
#    output:
#        "CliqueSNV_output/{sample}/"
#    log:
#        "log/hisat2_clique/{sample}_clique"
#    shell:
#        "java -jar -Xmx5120m clique-snv.jar -m snv-illumina -in {input} -outDir {output} -threads 8 2>{log}"
