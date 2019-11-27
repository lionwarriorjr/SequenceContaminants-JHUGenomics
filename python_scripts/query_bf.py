import sys
import pybloom_live
import argparse
import numpy as np
from tqdm import tqdm


def kmers(seq, k):
    kms = set()
    num_k = len(seq) - k + 1
    for i in range(num_k):
        km = seq[i:i + k]
        kms.add(km)
    return kms


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


def main():
    parse = argparse.ArgumentParser(description='Query Naive Bloom Filter')
    parse.add_argument('-ibf', nargs='+', help='Input reference bloom filters')
    parse.add_argument('-f', help='Input query fastq')
    args = parse.parse_args()

    reads_kms = []
    with open(args.f) as f:
        fq = parse_fastq(f)
    for read in fq:
        reads_kms.append(kmers(read, 20))

    output = np.zeros((len(fq), len(args.ibf)))
    for j, bf in enumerate(args.ibf):
        with open(bf, 'rb') as b:
            sbf_file = pybloom_live.ScalableBloomFilter.fromfile(b)
            for i, read_k in enumerate(reads_kms):
                output[i][j] = 1
                for k in read_k:
                    if k not in sbf_file:
                        output[i][j] = 0
                        break

    for out in output:
        sys.stdout.write(' '.join([str(x) for x in out]) + '\n')


if __name__ == '__main__':
    main()
