
rule merge_reads:
    input:
        read1 = "{sample}.R1.fastq.gz",
        read2 = "{sample}.R2.fastq.gz"
    params:
        script = "src/merge_reads.py"
    output:
        merge_gz = "{sample}.merged.fastq.gz",
        not_combined_1 = "{sample}.notCombined_1.fastq",
        not_combined_2 = "{sample}.notCombined_2.fastq"
    shell:
        """
        set -eux
        python {params.script} --read1 {input.read1} --read2 {input.read2} --out {output.merge_gz} 
        """