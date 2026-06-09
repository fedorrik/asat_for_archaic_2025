#!/bin/bash

stv_bed=$1         # Path to stv.bed file with all the assemblies
input_fasta=$2     # Path to fasta file with all the assemblies
hors_to_sample=$3  # Number of HORs to be sampled in echo stv of each cenhap
                   #     More HORs align longer

mkdir -p bed fa

# create bed files for each frequent stv of each cenhap
python3 select_stv.py $stv_bed $hors_to_sample
echo "Frequent stv bed files created"
echo

# process each bed file
for i in `ls bed`; do
  base=$(basename $i .bed);
  echo "Processing $base..."
  # transform to fasta and align
  bedtools getfasta -bed bed/$base.bed -fi $input_fasta -s | muscle > fa/$base.fa 2> current_muscle_progress.txt
  # make hmm
  hmmbuild $base.hmm fa/$base.fa > /dev/null
  # append to common hmm
  cat $base.hmm >> dumb_horhap.hmm
  # remove hmm
  rm $base.hmm
done

rm current_muscle_progress.txt

# RUN HMMER
#hmmer-run.sh cen12_fastas dumb_horhap.hmm
