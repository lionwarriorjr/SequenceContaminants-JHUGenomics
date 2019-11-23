import sys
import pybloom_live
import pickle


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
    seq_file = sys.argv[1]
    with open(seq_file, 'r') as f:
        seq = ''
        l = 0
        for line in f:
            if l > 1000:
                break
            if line.startswith('>'):
                continue
            seq = seq + line.rstrip()
            l += 1
    kms = kmers(seq, 20)
    print(kms[1])
    sbf = pybloom_live.ScalableBloomFilter(mode=pybloom_live.ScalableBloomFilter.LARGE_SET_GROWTH)
    count = 1
    for i in kms:
        _ = sbf.add(i)
        count += 1


    if kms[1] in sbf:
        print('kms in sbf')
    with open('test_bloom.txt', 'wb') as fh:
        sbf.tofile(fh)
    with open('test_bloom.txt', 'rb') as b:
        sbf_file = pybloom_live.ScalableBloomFilter.fromfile(b)
    if kms[1] in sbf_file:
        print('kms in file')


if __name__ == '__main__':
    main()
