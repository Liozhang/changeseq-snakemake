rule all:
    input:
        "sorted.bam"

rule check_genome_index:
    output:
        touch("genome_indexed.txt")
    shell:
        """
        BWA_path={BWA_path}
        HG19_path={HG19_path}
        index_files_extensions=".pac .amb .ann .bwt .sa"
        genome_indexed=True
        for extension in $index_files_extensions
        do
            if [ ! -f $HG19_path$extension ]; then
                genome_indexed=False
                break
            fi
        done
        if [ "$genome_indexed" = False ]; then
            echo 'Genome index files not detected. Running BWA to generate indices.'
            bwa_index_command="$BWA_path index $HG19_path"
            echo 'Running bwa command: $bwa_index_command'
            $bwa_index_command
            echo 'BWA genome index generated'
        else
            echo 'BWA genome index found.'
        fi
        """

rule align_reads:
    input:
        "genome_indexed.txt",
        "read1.fastq",
        "read2.fastq"
    output:
        "aligned.sam"
    shell:
        """
        BWA_path={BWA_path}
        HG19_path={HG19_path}
        bwa_alignment_command="$BWA_path mem $HG19_path {input[1]} {input[2]} > {output}"
        echo $bwa_alignment_command
        $bwa_alignment_command
        echo 'Paired end mapping completed.'
        """

rule sam_to_bam:
    input:
        "aligned.sam"
    output:
        "aligned.bam"
    shell:
        """
        samtools_sam_to_bam_command="samtools sort -o {output} {input}"
        echo $samtools_sam_to_bam_command
        $samtools_sam_to_bam_command
        echo 'Sorting by coordinate position complete.'
        """

rule index_bam:
    input:
        "aligned.bam"
    output:
        "aligned.bam.bai"
    shell:
        """
        samtools_index_command="samtools index {input}"
        echo $samtools_index_command
        $samtools_index_command
        echo 'Indexing complete.'
        """

rule sort_bam:
    input:
        "aligned.bam"
    output:
        "sorted.bam"
    shell:
        """
        samtools_sort_by_name_command="samtools sort -o {output} -n {input}"
        echo $samtools_sort_by_name_command
        $samtools_sort_by_name_command
        echo 'Sorting by name complete.'
        """