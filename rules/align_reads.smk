
rule build_genome_index:
    input:
        ref_genome = ""
    output:
        pac = "ref_genome.pac",
        amb = "ref_genome.amb",
        ann = "ref_genome.ann",
        bwt = "ref_genome.bwt",
        sa = "ref_genome.sa"
    params:
        bwa = "",
    threads: 2
    shell:
        """
        for extension in .pac .amb .ann .bwt .sa; do
            if [ ! -f $HG19_path$extension ]; then
                genome_indexed=False
                break
            fi
        done

        if [ "$genome_indexed" = False ]; then
            echo 'Genome index files not detected. Running BWA to generate indices.'
            {params.bwa} index {input.ref_genome}
            echo 'BWA genome index generated'
        else
            echo 'BWA genome index found.'
        fi
        """


rule align_reads:
    input:
        fq1 = "read1.fastq",
        fq2 = "read2.fastq",
        pac = "ref_genome.pac",
        amb = "ref_genome.amb",
        ann = "ref_genome.ann",
        bwt = "ref_genome.bwt",
        sa = "ref_genome.sa"
    output:
        bam = "aligned.bam"
    params:
        bwa = "",
        index = "ref_genome"
    threads: 5
    shell:
        """
        set -eux
        {params.bwa} mem -t {threads} {input.index} {input.fq1} {input.fq2} | samtools view -bS - > {output.bam}
        echo 'Paired end mapping completed.'
        """


rule sort_bam:
    input:
        bam = "aligned.bam"
    output:
        sort_bam = "sorted.bam"
    threads: 5
    shell:
        """
        set -eux
        samtools sort -@ {threads} -o {output.sort_bam} {input.bam}
        samtools index {output.sort_bam}
        echo 'Sorting by coordinate position complete.'
        """
