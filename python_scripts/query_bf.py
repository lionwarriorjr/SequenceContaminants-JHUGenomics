import sys
import os
import pybloom_live
import argparse
import numpy as np
from tqdm import tqdm

# Code from genome-sketching
def kmers(seq, k):
    kms = set()
    num_k = len(seq) - k + 1
    for i in range(num_k):
        km = seq[i:i + k]
        kms.add(km)
    return kms

# Code from lecture
def parse_fastq(fh):
    """ Parse reads from a FASTQ filehandle.  For each read, we
        return a name, nucleotide-string, quality-string triple. """
    reads = []
    q_score = []
    while True:
        first_line = fh.readline()
        if len(first_line) == 0:
            break  # end of file
        name = first_line[1:].rstrip()
        seq = fh.readline().rstrip()
        fh.readline()  # ignore line starting with +
        qual = fh.readline().rstrip()
        reads.append(seq)
        q_score.append(qual)
    return reads


def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)


def main():
    parse = argparse.ArgumentParser(description='Query Naive Bloom Filter')
    parse.add_argument('-bf', type=dir_path, help='Input reference bloom filters directory')
    parse.add_argument('-f', help='Input query fastq')
    args = parse.parse_args()

    reads_kms = []
    bf_files = ['e_coli_bf.txt', 's_maltophilia_bf.txt', 'a_xylosoxidans_bf.txt',
                'd_acidovorans_bf.txt', 'p_acnes_bf.txt', 's_epidermidis_bf.txt',
                'a_baumannii_bf.txt', 's_cerevisiae_bf.txt', 'w_anomalus_bf.txt',
                'm_restricta_bf.txt', 'p_notatum_bf.txt']
    bf_path = [os.path.join(args.bf, p) for p in bf_files]
    with open(args.f) as f:
        fq = parse_fastq(f)
    for read in fq:
        reads_kms.append(kmers(read, 20))

    output = np.zeros((len(fq), len(bf_path)), dtype=int)
    for j, bf in enumerate(bf_path):
        with open(bf, 'rb') as b:
            sbf_file = pybloom_live.ScalableBloomFilter.fromfile(b)
            for i, read_k in enumerate(reads_kms):
                read_kmer_len = len(read_k)
                mis = 0
                match = 0
                for k in read_k:
                    # if k in sbf_file:
                    #     output[i][j] = 1
                    #     break
                    if k in sbf_file:
                        match += 1
                        if match/read_kmer_len >= 0.2:
                            output[i][j] = 1
                            break
                        # if mis/read_kmer_len > 0.5:
                        #     break
                # if mis/read_kmer_len < 0.8:
                #     output[i][j] = 1

    score = np.zeros(len(bf_files), dtype=int)
    x, y = np.nonzero(output)
    for r in y:
        score[r] += 1

    # for out in output:
    #     sys.stdout.write(' '.join([str(x) for x in out]) + '\n')

    for b, s in zip(bf_files, score):
        print(str(b) + ': ' + str(s))


if __name__ == '__main__':
    main()
