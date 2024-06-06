
rule merge_reads:
    input:
        read1 = "{sample}.R1.fastq.gz",
        read2 = "{sample}.R2.fastq.gz"
    output:
        merge_gz = "{sample}.merged.fastq.gz",
        not_combined_1 = "{sample}.notCombined_1.fastq",
        not_combined_2 = "{sample}.notCombined_2.fastq"
    shell:
        """
        flash {input.read1} {input.read2} -o {wildcards.sample} -z
        """