import sys
import csv
import re


def parse_fastq(filename):
    """ Parse reads from a FASTQ filehandle.  For each read, we
        return a name, nucleotide-string, quality-string triple. """
    name = []
    sequence = []
    quality = []
    prefixname = re.sub('.fastq','',filename)
    with open(filename) as fh:
        while True:
            first_line = fh.readline()
            if len(first_line) == 0:
                break  # end of file
            nm = prefixname + '.' + first_line[1:].rstrip()
            seq = fh.readline().rstrip()
            fh.readline()  # ignore line starting with +
            qual = fh.readline().rstrip()
            name.append(nm)
            sequence.append(seq)
            quality.append(qual)
    return name,sequence,quality


def parse_source(sourcefile):
    with open(sourcefile) as f:
        rd = csv.reader(f, delimiter='\t')
        otherspecies = []
        mixprop = []  # * 100
        for row in rd:
            otherspecies.append(row[0])
            mixprop.append(row[1])
    return otherspecies,mixprop

def human_mix_multispecies(humanfile,sourcefile): # assume \t separated; 1st col: otherspecies's file name, 2nd col: mixing proportion
    humanName,humanSeq,humanQual = parse_fastq(humanfile)
    sourceSpecies,sourceProp = parse_source(sourcefile)
    humanLen = len(humanName)
    for i in range(0,len(sourceSpecies)):
        tmpName,tmpSeq,tmpQual = parse_fastq(sourceSpecies[i])
        tmpInsertions = round(float(sourceProp[i]) * len(tmpName)) # round to nearest integer
        for j in range(0,tmpInsertions):
            humanName.append(tmpName[j])
            humanSeq.append(tmpSeq[j])
            humanQual.append(tmpQual[j])
    return humanName,humanSeq,humanQual

## Example
HumanFile = sys.argv[1]
MixSourceFile = sys.argv[2]
# OutputFile = sys.argv[3]

# x1,x2,x3 = parse_fastq('chr21_sim.fastq')
# command line: python3 mixture.py chr21_sim.fastq test.txt > chr21_mix3species.fastq
y1,y2,y3 = human_mix_multispecies(HumanFile,MixSourceFile)
for k in range(0,len(y1)):
    print('@' + y1[k])
    print(y2[k])
    print('+')
    print(y3[k])






