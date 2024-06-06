import os
import argparse
import gzip
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import reverse_complement, Seq


def read_fq(file):
    if re.search('.gz$', file):
        fastq = SeqIO.parse(gzip.open(file, 'rt'), 'fastq')
    else:
        fastq = SeqIO.parse(open(file, 'r'), 'fastq')
    for record in fastq:
        yield record

            
def mergeReads(fq1_file, fq2_file, out):
    fq1_records = read_fq(fq1_file)
    fq2_records = read_fq(fq2_file)
    records = []
    for r1, r2 in zip(fq1_records, fq2_records):
        merged_seq = str(reverse_complement(r1.seq)) + str(r2.seq)
        merged_qual = r1.letter_annotations['phred_quality'][::-1] + r2.letter_annotations['phred_quality']
        records.append(
            SeqRecord(
                seq = Seq(merged_seq),
                id = r1.id,
                description = r1.description,
                letter_annotations = {'phred_quality': merged_qual}
            )
        )
        # Write to file every 1000 records
        if len(records) >= 10000:
            with gzip.open(out, 'wt') as f:
                SeqIO.write(records, f, 'fastq')
            records = []
    # Write remaining records
    with gzip.open(out, 'wt') as f:
        SeqIO.write(records, f, 'fastq')


def main():
    parser = argparse.ArgumentParser(description='Merge CIRCLE-seq reads for alignment.')
    parser.add_argument('--read1', help='Read 1 filename', required=True)
    parser.add_argument('--read2', help='Read 2 filename', required=True)
    parser.add_argument('--out', help='Output filename', required=True)

    args = parser.parse_args()
    
    # check input files
    if not os.path.exists(args.read1):
        print('Read 1 file does not exist')
        exit(1)
    
    if not os.path.exists(args.read2):
        print('Read 2 file does not exist')
        exit(1)
    
    # check output file
    if not args.out.endswith('gz') and args.out.endswith('fastq'):
        print('Output file must be a .fastq.gz or .fastq file')
        exit(1)

    try:
        mergeReads(args.read1, args.read2, args.out)
        print('Merging reads complete')
    except Exception as e:
        raise

if __name__ == "__main__":
    main()



