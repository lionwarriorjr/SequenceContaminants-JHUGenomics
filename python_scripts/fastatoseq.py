import sys
def fastatoseqperline(filename):
    name = []
    sequence = []
    with open(filename) as fh:
        while True:
            first_line = fh.readline()
            if len(first_line) == 0:
                break  # end of file
            nm = first_line[1:].rstrip()
            seq = fh.readline().rstrip()
            name.append(nm)
            sequence.append(seq)
    return sequence

fastafile = sys.argv[1]
fna_seq = fastatoseqperline(fastafile)
for k in range(0,len(fna_seq)):
    print(fna_seq[k])