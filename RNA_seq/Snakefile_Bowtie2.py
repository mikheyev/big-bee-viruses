# Rule 0: Build rule of all rules_
sample_bee = ['1275','1500_A','6107','6154','6246_A','7454','7538','9465','9540_A','332014','338391','338479']

rule all:
    input:
        expand("12_bowtie2_results/{sample}/{sample}.sam", sample=sample_bee),
        expand("12_bowtie2_results/{sample}/{sample}.bam", sample=sample_bee),
        expand("12_bowtie2_results/{sample}/{sample}_sorted.bam", sample=sample_bee),
        expand("12_bowtie2_results/{sample}/{sample}_sorted_stats/genome_results.txt", sample=sample_bee),
        expand("12_bowtie2_results/{sample}/{sample}_sorted_stats/report.pdf", sample=sample_bee),
        expand("12_bowtie2_results/{sample}/{sample}_sorted_mapped.bam", sample=sample_bee),
        expand("12_bowtie2_results/{sample}/{sample}_sorted_mapped.bai", sample=sample_bee),
        expand("13_bowtie2_coverage/{sample}/{sample}.txt", sample=sample_bee),
        expand("13_bowtie2_coverage/{sample}/{sample}_plot.txt", sample=sample_bee),
        expand("14_bowtie2_depth/{sample}/{sample}.txt", sample=sample_bee),
        expand("14_bowtie2_depth/{sample}/{sample}_meandepth.txt", sample=sample_bee),
        expand("15_bowtie2_stringtie_results/{sample}.gtf", sample=sample_bee),
        expand("15_stringtie_for_ballgown/{sample}/{sample}.gtf", sample=sample_bee)


# Rule 1: Map fastq to fasta index
rule bowtie2_mapping:
    input:
        "04_fastq_good/{sample}_1_good.fastq.gz",
        "04_fastq_good/{sample}_2_good.fastq.gz"
    output:
        "12_bowtie2_results/{sample}/{sample}.sam"
    log:
        "log/bowtie2_map/{sample}_bowtie2"
    shell:
        "bowtie2 -x 11_bowtie2_index/DWV/DWV_index -1 {input[0]} -2 {input[1]} -S {output} 2>{log}"

# Rule 2: Use samtools convert sam to bam
rule samtools_view:
    input:
        "12_bowtie2_results/{sample}/{sample}.sam"
    output:
        "12_bowtie2_results/{sample}/{sample}.bam"
    shell:
        "samtools view -Sb {input} > {output}"


# Rule 3: Sort the bam
rule samtools_sort:
    input:
        "12_bowtie2_results/{sample}/{sample}.bam"
    output:
        "12_bowtie2_results/{sample}/{sample}_sorted.bam"
    shell:
        "samtools sort {input} -o {output}"


# Rule 4: Examines sequencing alignment data
rule qulimap_bampc:
    input:
        "12_bowtie2_results/{sample}/{sample}_sorted.bam"
    output:
        "12_bowtie2_results/{sample}/{sample}_sorted_stats/genome_results.txt",
        "12_bowtie2_results/{sample}/{sample}_sorted_stats/report.pdf"
    shell:
        "qualimap bamqc -bam {input} -outformat pdf"


# Rule 5: Extrect mapped reads
rule samtools_view_mapped:
    input:
        "12_bowtie2_results/{sample}/{sample}_sorted.bam"
    output:
        "12_bowtie2_results/{sample}/{sample}_sorted_mapped.bam"
    shell:
        "samtools view -b -F4 {input} > {output}"


# Rule 6: Build the sorted bam hisat2_index
rule samtools_index:
    input:
        "12_bowtie2_results/{sample}/{sample}_sorted_mapped.bam"
    output:
        "12_bowtie2_results/{sample}/{sample}_sorted_mapped.bai"
    shell:
        "samtools index {input} {output}"


# Rule 7: Calculate coverage
rule samtools_coverage:
    input:
        "12_bowtie2_results/{sample}/{sample}_sorted_mapped.bam"
    output:
        "13_bowtie2_coverage/{sample}/{sample}.txt"
    shell:
        "samtools coverage {input} -o {output}"

# Rule 8: coverage plot
rule samtools_coverage_A:
    input:
        "12_bowtie2_results/{sample}/{sample}_sorted_mapped.bam"
    output:
        "13_bowtie2_coverage/{sample}/{sample}_plot.txt"
    shell:
        "samtools coverage -A {input} -o {output}"


# Rule 9: depth
rule samtools_depth:
    input:
        "12_bowtie2_results/{sample}/{sample}_sorted_mapped.bam"
    output:
        "14_bowtie2_depth/{sample}/{sample}.txt"
    shell:
        "samtools depth {input} > {output}"

# Rule 10: mean samtools_depth
rule samtools_meandepth:
    input:
        "14_bowtie2_depth/{sample}/{sample}.txt"
    output:
        "14_bowtie2_depth/{sample}/{sample}_meandepth.txt"
    shell:
        """cat {input} | awk "{{sum+=\$3}}END{{print sum/NR}}" > {output}"""

# Rule 11: assembly
rule stringtie_assembly:
    input:
        "12_bowtie2_results/{sample}/{sample}_sorted_mapped.bam",
        "00_reference/DWV.gff"
    output:
        "15_bowtie2_stringtie_results/{sample}.gtf"
    shell:
        "stringtie -p 8 -G {input[1]} {input[0]} -o {output}"


# Rule 12 preparing for Ballgown
rule stringtie_ballgown:
    input:
        "12_bowtie2_results/{sample}/{sample}_sorted_mapped.bam",
        "15_bowtie2_stringtie_results/stringtie_merge.gtf"
    output:
        "15_stringtie_for_ballgown/{sample}/{sample}.gtf"
    shell:
        "stringtie {input[0]} -e -B -p 8 -G {input[1]} -o {output}"
