# count_kmers_in_reads

To run kmers on one fastq sample:
- `./count_kmers_in_reads.sh <fastq> <kmers_tsv>`
- kmers_tsv must have two columns (kmer_name, kmer_sequence) without header

To collect results of multiple run:
- `python3 count_kmers_in_reads.sh --kmers <kmers_tsv>`

