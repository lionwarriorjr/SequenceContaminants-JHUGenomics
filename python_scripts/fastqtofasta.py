import sys
def fastqtofasta(filename):
    """ Parse reads from a FASTQ filehandle.  For each read, we
        return a name, nucleotide-string, quality-string triple. """
    name = []
    sequence = []
    with open(filename) as fh:
        while True:
            first_line = fh.readline()
            if len(first_line) == 0:
                break  # end of file
            nm = first_line[1:].rstrip()
            seq = fh.readline().rstrip()
            fh.readline()  # ignore line starting with +
            fh.readline()
            name.append(nm)
            sequence.append(seq)
            # quality.append(qual)
    return name,sequence

fastqfile = sys.argv[1]
fna_name,fna_seq = fastqtofasta(fastqfile)
for k in range(0,len(fna_name)):
    print('@' + fna_name[k])
    print(fna_seq[k])