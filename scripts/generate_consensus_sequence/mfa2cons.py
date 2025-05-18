from collections import Counter
from Bio.SeqIO import parse
from sys import argv


mfa = argv[1]
# create list with seqs
seqs = [record.seq for record in parse(mfa, 'fasta')]
# create consensus sequence
cons = ''
for i in range(len(seqs[0])):
    position = []
    for seq in seqs:
        base = seq[i]
        position.append(base)
    most_freq = Counter(position).most_common(2)
    # put "N" if there are two most freq bases (not "-")
    if len(most_freq) > 1 and most_freq[0][1] == most_freq[1][1] and most_freq[0][0] != '-' and most_freq[1][0] != '-':
        cons += 'N'
    # append nothing to cons if gap is most freq
    elif most_freq[0][0] == '-':
        continue
    # append most freq
    else:
        cons += most_freq[0][0]
# print cons
name = '.'.join(mfa.split('.')[:-1]).split('/')[-1]
print('>CONS-{}\n{}'.format(name, cons))

