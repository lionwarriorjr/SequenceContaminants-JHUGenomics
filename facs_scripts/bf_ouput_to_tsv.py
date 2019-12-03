import os
import csv


def parse_fastq(fh):
    """ Parse reads from a FASTQ filehandle.  For each read, we
        return a name, nucleotide-string, quality-string triple. """
    names = []
    reads = []
    q_score = []
    while True:
        first_line = fh.readline()
        if len(first_line) == 0:
            break  # end of file
        name = first_line[1:].rstrip()
        names.append(name)
        seq = fh.readline().rstrip()
        fh.readline()  # ignore line starting with +
        qual = fh.readline().rstrip()
        reads.append(seq)
        q_score.append(qual)
    return reads, names


con_path = './contaminated_fastq'
output = {}
SPECIES = {'e_coli.bloom.fastq': 0,
           's_maltophilia.bloom.fastq': 1,
           'a_xylosoxidans.bloom.fastq': 2,
           'd_acidovorans.bloom.fastq': 3,
           'p_acne.bloom.fastq': 4,
           's_epidermidis.bloom.fastq': 5,
           'a_baumannii.bloom.fastq': 6,
           's_cerevisiae.bloom.fastq': 7,
           'w_anomalus.bloom.fastq': 8,
           'm_restricta.bloom.fastq': 9,
           'p_notatum.bloom.fastq': 10}
for file in os.listdir(con_path):
    con_file_path = os.path.join(con_path, file)
    if os.path.getsize(con_file_path) and file.endswith('fastq'):
        species_name = SPECIES[file]
        with open(con_file_path, 'r') as f:
            reads, names = parse_fastq(f)
            for r, n in zip(reads, names):
                if n in output:
                    output[n][1].append(species_name)
                else:
                    output[n] = [r, [species_name]]

with open(os.path.join(con_path, 'facs_output.tsv'), 'wt') as t:
    tsv_writer = csv.writer(t, delimiter='\t')
    for na, re in output.items():
        if len(re[1]) == 1:
            tsv_writer.writerow([na, re[0], re[1][0]])
        else:
            sn = ''
            for i in re[1]:
                if sn:
                    sn = str(sn) + ',' + str(i)
                else:
                    sn = str(i)
            tsv_writer.writerow([na, re[0], sn])
