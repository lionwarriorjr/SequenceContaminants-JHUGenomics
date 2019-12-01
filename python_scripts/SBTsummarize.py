import sys

def parse_fasta(fastafile):
    name = []
    sequence = []
    with open(fastafile) as fh:
        while True:
            first_line = fh.readline()
            if len(first_line) == 0:
                break  # end of file
            nm = first_line[1:].rstrip()
            seq = fh.readline().rstrip()
            name.append(nm)
            sequence.append(seq)
    return name,sequence

def parse_sbtresult(sbtfile):
    ## read in sbt result file
    fh = open(sbtfile).readlines()
    id = []
    species = []
    ## matched id and its species
    k = -1
    for i in range(0, len(fh)):
        line = fh[i]
        if line.startswith("*"):
            k = k+1
            if int(line.split()[1]) == 1:
                spec = fh[i + 1].split("/")[-1].replace(".bf.bv", "").rstrip()
                id.append(k)
                species.append(spec)
    return id,species



fastafile = sys.argv[1]#'unmapped_mix3species.fasta'
sbtfile = sys.argv[2] #'output-20.txt'

name,seq = parse_fasta(fastafile)
id,species = parse_sbtresult(sbtfile)


for j in range(0,len(id)):
    print(name[id[j]]+'\t'+species[j])










