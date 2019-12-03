import pybloom_live
import pickle
import argparse
import os
from tqdm import tqdm


def kmers(seq, k):
    kms = []
    num_k = len(seq) - k + 1
    for i in range(num_k):
        km = seq[i:i + k]
        kms.append(km)
    return kms


def store_data(dict, filename):
    file = open(filename, 'ab')
    pickle.dump(dict)
    file.close()


def load_data(filename):
    file = open(filename, 'rb')
    d = pickle.load(file)
    file.close()
    return d


def main():
    parse = argparse.ArgumentParser(description='Build Naive Bloom Filter')
    parse.add_argument('-i', help='Input reference fasta directory')
    parse.add_argument('-o', help='Output bloom filter directory')
    args = parse.parse_args()
    SPECIES = {'GCF_000005845.2_ASM584v2_genomic.fna': 'e_coli_bf.txt',
               'GCF_000072485.1_ASM7248v1_genomic.fna': 's_maltophilia_bf.txt',
               'GCF_001457475.1_NCTC10807_genomic.fna': 'a_xylosoxidans_bf.txt',
               'GCF_000018665.1_ASM1866v1_genomic.fna': 'd_acidovorans_bf.txt',
               'GCF_000008345.1_ASM834v1_genomic.fna': 'p_acnes_bf.txt',
               'GCF_000007645.1_ASM764v1_genomic.fna': 's_epidermidis_bf.txt',
               'GCF_000746645.1_ASM74664v1_genomic.fna': 'a_baumannii_bf.txt',
               'GCF_000146045.2_R64_genomic.fna': 's_cerevisiae_bf.txt',
               'GCF_001661255.1_Wican1_genomic.fna': 'w_anomalus_bf.txt',
               'GCF_003290485.1_ASM329048v1_genomic.fna': 'm_restricta_bf.txt',
               'GCA_000710275.1_ASM71027v1_genomic.fna': 'p_notatum_bf.txt'}

    for ip, op in SPECIES.items():
        seq_file = os.path.join(args.i, ip)
        with open(seq_file, 'r') as f:
            seq = ''
            # l = 0
            for line in f:
                # if l > 1000:
                #     break
                if line.startswith('>'):
                    continue
                seq = seq + line.rstrip()
                # l += 1
        kms = kmers(seq, 20)
        print(op)
        print('number of kmers', len(kms))
        sbf = pybloom_live.ScalableBloomFilter(mode=pybloom_live.ScalableBloomFilter.LARGE_SET_GROWTH,
                                               error_rate=0.0005)

        for i, x in zip(kms, tqdm(range(len(kms)))):
            _ = sbf.add(i)

        if next(iter(kms)) in sbf:
            print('kms in sbf')
        out_file = os.path.join(args.o, op)
        with open(out_file, 'wb') as fh:
            sbf.tofile(fh)
        with open(out_file, 'rb') as b:
            sbf_file = pybloom_live.ScalableBloomFilter.fromfile(b)
        if next(iter(kms)) in sbf_file:
            print('kms in file')


if __name__ == '__main__':
    main()
