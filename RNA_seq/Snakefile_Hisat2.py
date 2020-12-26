# Rule 0: Build rule of all rules_
sample_bee = ['330710','331014_B','333031','334300','334303','334304','334310','334705','334706','334712']

rule all:
    input:
        expand("02_fastqc_results/{sample}/{sample}_1_fastqc.zip", sample=sample_bee),
        expand("02_fastqc_results/{sample}/{sample}_2_fastqc.zip", sample=sample_bee),
        expand("02_fastqc_results/{sample}/{sample}_1_fastqc.html", sample=sample_bee),
        expand("02_fastqc_results/{sample}/{sample}_2_fastqc.html", sample=sample_bee),
        expand("04_fastq_good/{sample}_1_good.fastq.gz", sample=sample_bee),
        expand("04_fastq_good/{sample}_2_good.fastq.gz", sample=sample_bee),
        expand("03_fastp_results/{sample}/{sample}.html", sample=sample_bee),
        expand("03_fastp_results/{sample}/{sample}.json", sample=sample_bee),
        expand("06_hisat2_results/{sample}/{sample}.sam", sample=sample_bee),
        expand("06_hisat2_results/{sample}/{sample}.bam", sample=sample_bee),
        expand("06_hisat2_results/{sample}/{sample}_sorted_stats/genome_results.txt", sample=sample_bee),
        expand("06_hisat2_results/{sample}/{sample}_sorted_stats/report.pdf", sample=sample_bee),
        expand("06_hisat2_results/{sample}/{sample}_sorted.bam", sample=sample_bee),
        expand("06_hisat2_results/{sample}/{sample}_sorted_mapped.bam", sample=sample_bee),
        expand("06_hisat2_results/{sample}/{sample}_sorted_mapped.bai", sample=sample_bee),
        expand("07_hisat2_coverage/{sample}/{sample}.txt", sample=sample_bee),
        expand("07_hisat2_coverage/{sample}/{sample}_plot.txt", sample=sample_bee),
        expand("08_hisat2_depth/{sample}/{sample}.txt", sample=sample_bee),
        expand("08_hisat2_depth/{sample}/{sample}_meandepth.txt", sample=sample_bee)
        #expand("09_hisat2_stringtie_results/{sample}.gtf", sample=sample_bee),
        #expand("09_stringtie_for_ballgown/{sample}/{sample}.gtf", sample=sample_bee)


# Rule 1: fastqc
rule fastqc:
    input:
        "01_fastq/{sample}_1.fastq.gz",
        "01_fastq/{sample}_2.fastq.gz"

    output:
        "02_fastqc_results/{sample}/{sample}_1_fastqc.zip",
        "02_fastqc_results/{sample}/{sample}_2_fastqc.zip",
        "02_fastqc_results/{sample}/{sample}_1_fastqc.html",
        "02_fastqc_results/{sample}/{sample}_2_fastqc.html"

    log:
        "log/fastqc/{sample}_fastqc"

    shell:
        "fastqc {input} -o 02_fastqc_results/{wildcards.sample} --noextract {input[0]} {input[1]} 2>{log}"

# RUle 2: fastp
rule fastp:
    input:
        "01_fastq/{sample}_1.fastq.gz",
        "01_fastq/{sample}_2.fastq.gz"

    output:
        "04_fastq_good/{sample}_1_good.fastq.gz",
        "04_fastq_good/{sample}_2_good.fastq.gz",
        "03_fastp_results/{sample}/{sample}.html",
        "03_fastp_results/{sample}/{sample}.json"

    log:
        "log/fastp/{sample}_fastp"

    shell:
        "fastp -h 03_fastp_results/{wildcards.sample}/{wildcards.sample}.html -j 03_fastp_results/{wildcards.sample}/{wildcards.sample}.json -i {input[0]} -I {input[1]} -o {output[0]} -O {output[1]} 2>{log}"


# Rule 3: Map fastq to fasta index
rule hisat2_mapping:
    input:
        "04_fastq_good/{sample}_1_good.fastq.gz",
        "04_fastq_good/{sample}_2_good.fastq.gz"

    output:
        "06_hisat2_results/{sample}/{sample}.sam"

    log:
        "log/hisat2_map/{sample}_HISAT2"

    shell:
        "hisat2 -x 05_hisat2_index/BQCV/BQCV -1 {input[0]} -2 {input[1]} -S {output} 2>{log}"


# Rule 4: Use samtools convert sam to bam
rule samtools_view:
    input:
        "06_hisat2_results/{sample}/{sample}.sam"
    output:
        "06_hisat2_results/{sample}/{sample}.bam"
    shell:
        "samtools view -Sb {input} > {output}"


# Rule 5: Sort the bam
rule samtools_sort:
    input:
        "06_hisat2_results/{sample}/{sample}.bam"
    output:
        "06_hisat2_results/{sample}/{sample}_sorted.bam"
    shell:
        "samtools sort {input} -o {output}"


# Rule 6: Examines sequencing alignment data
rule qulimap_bampc:
    input:
        "06_hisat2_results/{sample}/{sample}_sorted.bam"
    output:
        "06_hisat2_results/{sample}/{sample}_sorted_stats/genome_results.txt",
        "06_hisat2_results/{sample}/{sample}_sorted_stats/report.pdf"
    shell:
        "qualimap bamqc -bam {input} -outformat pdf"


# Rule 7: Extrect mapped reads
rule samtools_view_mapped:
    input:
        "06_hisat2_results/{sample}/{sample}_sorted.bam"
    output:
        "06_hisat2_results/{sample}/{sample}_sorted_mapped.bam"
    shell:
        "samtools view -b -F4 {input} > {output}"

# Rule 8: Build the sorted bam hisat2_index
rule samtools_index:
    input:
        "06_hisat2_results/{sample}/{sample}_sorted_mapped.bam"
    output:
        "06_hisat2_results/{sample}/{sample}_sorted_mapped.bai"
    shell:
        "samtools index {input} {output}"

# Rule 9: Calculate coverage
rule samtools_coverage:
    input:
        "06_hisat2_results/{sample}/{sample}_sorted_mapped.bam"

    output:
        "07_hisat2_coverage/{sample}/{sample}.txt"

    shell:
        "samtools coverage {input} -o {output}"

# Rule 10: coverage plot
rule samtools_coverage_A:
    input:
        "06_hisat2_results/{sample}/{sample}_sorted_mapped.bam"
    output:
        "07_hisat2_coverage/{sample}/{sample}_plot.txt"
    shell:
        "samtools coverage -A {input} -o {output}"


# Rule 11: depth
rule samtools_depth:
    input:
        "06_hisat2_results/{sample}/{sample}_sorted_mapped.bam"
    output:
        "08_hisat2_depth/{sample}/{sample}.txt",
    shell:
        "samtools depth {input} > {output}"

# Rule 12: mean samtools_depth
rule samtools_meandepth:
    input:
        "08_hisat2_depth/{sample}/{sample}.txt"
    output:
        "08_hisat2_depth/{sample}/{sample}_meandepth.txt"
    shell:
        """cat {input} | awk "{{sum+=\$3}}END{{print sum/NR}}" > {output}"""


# Rule 13: assembly
#rule stringtie_assembly:
#    input:
#        "06_hisat2_results/{sample}/{sample}_sorted_mapped.bam",
#        "00_reference/DWV.gff"
#    output:
#        "09_hisat2_stringtie_results/{sample}.gtf"
#    shell:
#        "stringtie -p 8 -G {input[1]} {input[0]} -o {output}"

# Rule 14 preparing for Ballgown
#rule stringtie_ballgown:
#    input:
#        "06_hisat2_results/{sample}/{sample}_sorted_mapped.bam",
#        "09_hisat2_stringtie_results/stringtie_merge.gtf"
#    output:
#        "09_stringtie_for_ballgown/{sample}/{sample}.gtf"
#    shell:
#        "stringtie {input[0]} -e -B -p 8 -G {input[1]} -o {output}"
