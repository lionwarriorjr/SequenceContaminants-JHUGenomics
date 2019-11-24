import pybloom_live
import pickle
import argparse
from tqdm import tqdm


def kmers(seq, k):
    kms = set()
    num_k = len(seq) - k + 1
    for i in range(num_k):
        km = seq[i:i + k]
        kms.add(km)
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
    parse.add_argument('-i', nargs='+', help='Input reference fasta')
    parse.add_argument('-o', nargs='+', help='Output bloom filter path')
    args = parse.parse_args()
    if len(args.i) != len(args.o):
        parse.error('Numbers of Input and Output do not match')
    for ip, op in zip(args.i, args.o):
        seq_file = ip
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
        print(ip)
        print('number of kmers', len(kms))
        sbf = pybloom_live.ScalableBloomFilter(mode=pybloom_live.ScalableBloomFilter.LARGE_SET_GROWTH)

        for i, x in zip(kms, tqdm(range(len(kms)))):
            _ = sbf.add(i)

        if next(iter(kms)) in sbf:
            print('kms in sbf')
        with open(op, 'wb') as fh:
            sbf.tofile(fh)
        with open(op, 'rb') as b:
            sbf_file = pybloom_live.ScalableBloomFilter.fromfile(b)
        if next(iter(kms)) in sbf_file:
            print('kms in file')


if __name__ == '__main__':
    main()
