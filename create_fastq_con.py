import sys
from itertools import islice
from shutil import copyfile


def append_to_human_fastq(human, bacteria):
    copyfile(human, 'hu_con.fastq')
    with open('hu_con.fastq', 'a') as f:
        num_lines = int(sum(1 for line in open(human))/4) * 0.01
        with open(bacteria, "r") as b:
            head = list(islice(b, int(num_lines*4+1)))
            f.write('\n')
            for item in head:
                f.write(item.rstrip())
                f.write('\n')

    return head


def main():
    human_fastq = sys.argv[1]
    bac_fastq = sys.argv[2]
    con_reads = append_to_human_fastq(human_fastq, bac_fastq)

    for o in con_reads:
        sys.stdout.write(''.join([str(x) for x in o]))

if __name__ == '__main__':
    main()
