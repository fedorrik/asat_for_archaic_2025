import argparse
import pandas as pd
from Bio.SeqIO import parse
from Bio.Seq import Seq

## Call with
# python3 .py \
#      --stv "S3C11H1L.1-2_3/4_5" \
#      -- HOR-SF_table.tsv \

## Notes
# How HOR-SF_table.tsv was created: 
#   HOR monomer consensus sequences were obtained from AS-HORs-hmmer3.4-071024.hmm
#   and SF-HMMER classifier were run on these consensuses to determine the SF class 
#   of each HOR monomer. The sequence of each SF class were obteined from 
#   databaseS13_SF_cons_monomers.aln.txt
#   (supplimentary files of [Altemos et all, 2022])

###############################################################################
##                          Function Definitions                             ##
###############################################################################

def read_hor_sf_table(hor_sf_table_file):
    hor_sf_table = pd.read_csv(hor_sf_table_file, sep='\t').set_index('HOR')
    return hor_sf_table


def decompress_stv_name(stv_name):
    hor_name, input_mon_numbers = stv_name.split('.')
    input_mon_numbers = input_mon_numbers.split('_')
    output_mon_numbers = []
    for mon_number in input_mon_numbers:
        # expand collapsed line of natural numbers
        if '-' in mon_number:
            start, end = mon_number.split('-')
            # positive strand (e.g. 1-8)
            if end > start:
                expanded_numbers = [str(i) for i in range(int(start), int(end)+1)]
                decompress_stv_name.to_reverse = False
            # negative strand (e.g. 8-1)
            else:
                expanded_numbers = [str(i) for i in range(int(start), int(end)-1, -1)]
                decompress_stv_name.to_reverse = True
            output_mon_numbers += expanded_numbers
        else:
            output_mon_numbers.append(mon_number)
    output_monomers = ['{}.{}'.format(hor_name, i) for i in output_mon_numbers]
    return output_monomers


def reconstruct_ancestral_hor(hor_monomers, hor_sf_table):
    ancestral_hor_sequence = ''
    for hor_monomer in hor_monomers:
        sf_mon_seq = hor_sf_table.loc[hor_monomer, 'SF_seq']
        # reverse compliment if needed
        if decompress_stv_name.to_reverse:
            sf_mon_seq = str(Seq(sf_mon_seq).reverse_complement())
        ancestral_hor_sequence += sf_mon_seq
    return ancestral_hor_sequence


def main():
    parser = argparse.ArgumentParser(description='Reconstruct ancestral HOR for a given StV')

    parser.add_argument('--stv', '-s', type=str, action='store', 
                         help='StV name string (e.g. S3C11H1L.1-2_3/4_5)')
    parser.add_argument('--hor_sf_table', '-t', type=str, action='store',
                         help='Table with SF sequence for every HOR monomer')
    
    args = parser.parse_args()

    
    # main stuff
    hor_sf_table = read_hor_sf_table(args.hor_sf_table)
    decompressed_stv = decompress_stv_name(args.stv)
    seq = reconstruct_ancestral_hor(decompressed_stv, hor_sf_table)
    print('>anc\n{}'.format(seq))

if __name__ == '__main__':
    main()
