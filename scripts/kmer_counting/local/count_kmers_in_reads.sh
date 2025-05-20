#!/bin/bash
# Usage: ./count_kmers_in_reads.sh <fastq> <kmers_tsv>
# kmers_tsv must have two columns without header:
#   1. kmer_name
#   2. kmer_sequence

#SBATCH --job-name=kmer_counter
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00

reads_path=$1
kmers_tsv=$2

fastq=`echo $reads_path | awk '{n=split($0, name, "/"); print name[n]}'`
script_dir="$(dirname "$(readlink -f "$0")")"

echo $fastq

for fv in `cut -f 2 $kmers_tsv`; do

  rc=`echo $fv | tr ACGTacgt TGCAtgca | rev`
  echo -n "$fv $rc "; date +"%D %T"

  # fastq.gz
  if (echo $reads_path | grep -q ".fastq.gz" ); then
    zgrep -oP "$fv|$rc" $reads_path | wc -l >> "$fastq.cnt"

  # fastq
  elif (echo $reads_path | grep -q ".fastq" ); then
    grep -oP "$fv|$rc" $reads_path | wc -l >> "$fastq.cnt"

  # bam
  elif (echo $reads_path | grep -q ".bam" ); then
    samtools view $reads_path | grep -oP "$fv|$rc" | wc -l >> "$fastq.cnt"
  
  # wrng file
  else
    echo "Wrong input file format (.fastq/.fastq.gz/.bam required)"
    echo "Usage: ./count_kmers_in_reads.sh kmer.tsv reads" 

  fi

done
