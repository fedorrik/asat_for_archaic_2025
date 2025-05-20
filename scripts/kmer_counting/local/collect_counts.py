import argparse
import pandas as pd
from os import listdir


def collect_counts(kmers_file):

    datadir = [i for i in listdir('./') if '.cnt' in i]
    
    kmers = pd.read_csv(kmers_file, sep='\t', names=['kmer_name', 'kmer_seq'])
    kmers = list(kmers['kmer_name'])

    df = pd.DataFrame(columns=kmers)

    for i in datadir:
        print(i)
        df_sample = pd.read_csv(i, names=[i]).T.set_axis(kmers, axis=1)
        df = pd.concat([df, df_sample], axis=0)

    df.T.to_csv('cnttable', sep='\t')


def main():
    parser = argparse.ArgumentParser(description='Collect all .cnt files in current directory to one cnttable')
    parser.add_argument('--kmers', '-k', type=str, action='store', help='tsv file with kmer names and sequences')
    args = parser.parse_args()

    collect_counts(args.kmers)
    

if __name__ == '__main__':
    main()
