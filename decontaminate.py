import os
import csv
import sys
import numpy as np
import random
import subprocess
from subprocess import DEVNULL, STDOUT
import pandas as pd
import string
import torch
import torch.nn as nn
import time
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from collections import Counter, defaultdict
from dna2vec.dna2vec.multi_k_model import MultiKModel


# Constants for Simulation

# length of reads for RNN embeddings
READ_LENGTH = 90
# how many bases should be written per line of a fasta file
LINE_LENGTH = 90

# default header for written fasta file
DEFAULT_HEADER = ">reads.fna"
# header in written fasta files from unmapped reads in overlap analysis
OVERLAP_HEADER = ">Overlapped from Contaminant Genomes"
# temp file that overlapped reads are written to
OVERLAP_OUTPUT_FILE = "overlap.fna"

# Constants for Bowtie2 alignment
UNMAPPED_OUTPUT_FILE = "unmapped.fastq"
MAPPED_OUTPUT_FILE = "mapped.fastq"
TABLE_LEN = 13
UNMAPPED_STATUS_COL = 1
UNMAPPED_FLAG = 4
MAPPED_FLAG = 0
NAME_COL = 0
READ_COL = 9
# initialize names of pre-built Bowtie indexes in current directory
indexes = [
    'coli',
    'maltophilia',
    'xylosoxidans',
    'acidovorans',
    'acnes',
    'epidermidis',
    'baumannii',
    'cerevisiae',
    'anomalus',
    'restricta',
    'notatum'
]

# Constants for RNN evaluation

# path to trained DL model
MODEL_PATH = "rnn.pt"
# threshold probability of max target class in RNN prediction
RNN_THRESHOLD = 0.5
# identities of contaminant sources with pre-built Bowtie indexes
all_categories = [
    'E. coli',
    'S. maltophilia',
    'A. xylosoxidans',
    'D. acidovorans',
    'P. acnes',
    'S. epidermidis',
    'A. baumannii',
    'S. cerevisiae',
    'W. anomalus',
    'M. restricta',
    'P. notatum'
]
n_categories = len(all_categories)
# length of reads during RNN training
read_length = 90
# k-mer length during RNN training
k = 3
n_kmers = read_length - k + 1
# load pretrained DNA word embeddings from disk
filepath = 'dna2vec/pretrained/dna2vec-20161219-0153-k3to8-100d-10c-29320Mbp-sliding-Xat.w2v'
# load pretrained word embeddings as global variable
mk_model = MultiKModel(filepath)

# Pytorch RNN module definition
class RNN(nn.Module):
    def __init__(self, input_size, hidden_size, output_size):
        super(RNN, self).__init__()
        self.hidden_size = hidden_size
        # 2 hidden layers for simple internal structure
        # hidden layer is recurrent - feeds in as input to next layer
        self.i2h = nn.Linear(input_size + hidden_size, hidden_size)
        self.i2o = nn.Linear(input_size + hidden_size, output_size)
        # softmax logits for probability estimates in output layer
        self.softmax = nn.LogSoftmax(dim=1)

    def forward(self, input, hidden):
        combined = torch.cat((input, hidden), 1)
        hidden = self.i2h(combined)
        output = self.i2o(combined)
        output = self.softmax(output)
        return output, hidden

    def initHidden(self):
        return torch.zeros(1, self.hidden_size)

# load trained RNN model from disk
rnn = torch.load(MODEL_PATH)

# Constants for Bloom Filter

# path to written Bloom Filter fastq files
CON_PATH = './build/facs_build/contaminated_fastq'
# species map from Bloom filter name to species id
BLOOM = {
    'e_coli.bloom.fastq': 0,
    's_maltophilia.bloom.fastq': 1,
    'a_xylosoxidans.bloom.fastq': 2,
    'd_acidovorans.bloom.fastq': 3,
    'p_acne.bloom.fastq': 4,
    's_epidermidis.bloom.fastq': 5,
    'a_baumannii.bloom.fastq': 6,
    's_cerevisiae.bloom.fastq': 7,
    'w_anomalus.bloom.fastq': 8,
    'm_restricta.bloom.fastq': 9,
    'p_notatum.bloom.fastq': 10
}
# output file for Bloom Filter predictions
FACS_OUTPUT = "facs.tsv"


# read fasta file into list of reads in-memory
def parse_fasta(filename, read_length, limit=-1):
    reads = []
    read, n_char = "", "N"
    with open(filename, "r") as fa:
        read = ""
        for line in fa:
            if len(reads) == limit:
                return reads
            # ignore first line of each fasta input file
            if line.startswith(">"):
                continue
            # for consistency retain only upper-cases bases
            line = line.upper()
            end = min(read_length - len(read), len(line))
            offset = 0
            while offset+end < len(line):
                # extend read when part of same line
                read += line[offset:offset+end]
                if n_char not in read:
                    reads.append(read)
                    if limit != -1 and len(reads) >= limit:
                        return reads
                offset += end
                end = read_length
                read = ""
            read += line[offset:-1]
            if len(read) == READ_LENGTH:
                if n_char not in read:
                    reads.append(read)
                    if limit != -1 and len(reads) >= limit:
                        return reads
                # start the next read
                read = ""
        if read:
            if n_char not in read:
                reads.append(read)
    return reads


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


# write reads to fasta file on disk 
def write_fasta(filename, reads, names):
    with open(filename, "w") as fw:
        for i in range(len(reads)):
            fw.write('>' + names[i] + '\n')
            fw.write(reads[i] + '\n')


def contaminate(src_reads, contaminant):
    # introduce contaminants into source reads during simuulation
    n_contaminants = int(CONTAMINATION_RATE * len(src_reads))
    print("introducing %d contaminants from contaminant %s"%(n_contaminants, contaminant))
    contaminant_reads = parse_fasta(contaminant, READ_LENGTH, n_contaminants)
    src_reads_copy = src_reads.copy()
    src_reads_copy.extend(contaminant_reads)
    return src_reads_copy


def write_facs_output():
    output = {}
    for f in os.listdir(CON_PATH):
        con_file_path = os.path.join(CON_PATH, f)
        if os.path.getsize(con_file_path) and f.endswith('fastq'):
            species_name = BLOOM[f]
            with open(con_file_path, 'r') as f:
                reads, names = parse_fastq(f)
                for r, n in zip(reads, names):
                    if n in output:
                        output[n][1].append(species_name)
                    else:
                        output[n] = [r, [species_name]]

    with open("temp/" + FACS_OUTPUT, 'wt') as t:
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


# apply word embeddings to DNA reads
def lineToTensor(line):
    tensor = torch.zeros(len(line)-k+1, 1, 100)
    for i in range(len(line)-k+1):
        # embed reads through pre-trained DNA word embeddings
        tensor[i][0] = torch.FloatTensor(mk_model.vector(line[i:i+k].upper()))
    return tensor


# forward propagation through network
def evaluate(line_tensor):
    hidden = rnn.initHidden()
    for i in range(line_tensor.size()[0]):
        # forward propagate through trained RNN
        output, hidden = rnn(line_tensor[i], hidden)
    return output


# generate estimates for target class and associated target probability
def predict(input_line, n_predictions=3):
    with torch.no_grad():
        # generate output logits from forward propagation
        output = evaluate(lineToTensor(input_line))
        # softmax the logits to evaluate target class probabilities
        prob = torch.softmax(output, dim=1)
        # Get top N categories
        topv, topi = output.topk(n_predictions, 1, True)
        predictions = []
        for i in range(n_predictions):
            value = topv[0][i].item()
            # extract class for estimate
            category_index = topi[0][i].item()
            # extract target propbability for the max estimate
            max_prob = prob[0,category_index].numpy()
            predictions.append((category_index, max_prob))
    return predictions


def main():

    # name of pre-built index for source genome in current directory
    index = sys.argv[1]
    # name of file with intermixed contaminant reads
    contaminated_file = sys.argv[2]
    # how many of the unmapped reads should be verified by Bowtie2
    n_args = len(sys.argv)-1
    # optional third command-line argument
    final_output = sys.argv[3]
    
    # run bowtie alignment to find mapped/unmapped reads
    process = subprocess.Popen( args=["bowtie2",
                                "-x", index,
                                "-U", contaminated_file,
                                "--un", "temp/" + UNMAPPED_OUTPUT_FILE ],
                                stdout=DEVNULL, stderr=STDOUT )
    # block until bowtie subprocess completes for initial alignment
    process.wait()

    # extract unmapped reads from Bowtie2 output
    print("Logging unmapped reads from Bowtie2 alignment of contaminated genome")
    unmapped_reads = []
    extract_read = False
    with open("temp/" + UNMAPPED_OUTPUT_FILE, 'r') as fr:
        for line in fr:
            if extract_read:
                unmapped_reads[len(unmapped_reads)-1][1] = line.strip()
                current = unmapped_reads[len(unmapped_reads)-1]
                print('unmapped read (name, text): %s, %s'%(current[0], current[1]))
            extract_read = line.startswith("@")
            if extract_read:
                unmapped_reads.append([line[1:].strip(), ""])
                
    print("loaded pretrained DNA word embeddings and trained RNN model from disk\n")
    print(rnn)

    # forward propagate through RNN to get predictions for each read
    rnn_pred = {}
    for (name, unmapped_read) in unmapped_reads:
        # run the unmapped read through the RNN
        if 'N' not in unmapped_read:
            results = predict(unmapped_read)
            categories = [ret[0] for ret in results]
            probabilities = [ret[1] for ret in results]
            rnn_pred[name] = [unmapped_read, categories, probabilities]

    # get Bloom Filter predictions
    process = subprocess.Popen( args = ["./facs_remove.sh"], stdout=DEVNULL, stderr=STDOUT )
    process.wait()
    # write sketching predictions to current directory
    write_facs_output()
    overlapped_reads = rnn_pred
    with open("temp/" + FACS_OUTPUT, 'r') as fr:
        for line in fr:
            values = line.split('\t')
            name, text, categories = values[0], values[1], values[2]
            categories = categories.strip().split(',')
            categories = [int(c) for c in categories]
            if name in overlapped_reads:
                both = set(rnn_pred[name][1]).intersection(categories)
                overlap_indices = [rnn_pred[name][1].index(x) for x in both]
                overlap_categories, overlap_probs = [], []
                for i in overlap_indices:
                    overlap_categories.append(rnn_pred[name][1][i])
                    overlap_probs.append(rnn_pred[name][2][i])
                for category in categories:
                    if category not in both:
                        overlap_categories.append(category)
                overlapped_reads[name][1] = overlap_categories
                overlapped_reads[name][2] = overlap_probs

    # overlap classifications from DL and Bloom Filter are our predictions
    print('predicted set of unmapped reads mapped to contaminant genomes:')
    class_preds = [list() for i in range(n_categories)]
    class_names = [list() for i in range(n_categories)]
    with open(final_output, 'w') as fw:
        for name in overlapped_reads:
            results = overlapped_reads[name]
            text, categories, probabilities = results[0], results[1], results[2]
            for i in range(len(categories)-len(probabilities)):
                probabilities.append('none')
            predicted = []
            for category in categories:
                predicted.append(all_categories[category])
                class_names[category].append(name)
                class_preds[category].append(text)
            write_pred, write_prob = str(predicted), str(probabilities)
            print('read %s introduced by contaminant genome(s) %s with probabilities %s\n'%(name, write_pred, write_prob))
            fw.write("%s\t%s\t%s\t%s\n"%(name, text, write_pred, write_prob))

    for i in range(n_categories):
        write_fasta(filename="temp/" + indexes[i] + ".fna", 
                    reads=class_preds[i],
                    names=class_names[i])
    
    # verify claimed reads from overlap analysis with Bowtie2
    print("Verifying unmapped reads against contaminant genome with Bowtie2 alignment")
    reads = set()
    for i in range(n_categories):
        if len(class_preds[i]) > 0:
            # run the verification alignment
            process = subprocess.Popen( args=["bowtie2",
                                        "-x", "indexes/" + indexes[i],
                                        "-f", "temp/" + indexes[i] + ".fna",
                                        "--al", "temp/" + MAPPED_OUTPUT_FILE],
                                        stdout=DEVNULL, stderr=STDOUT )
            process.wait()
            extract_read = False
            with open("temp/" + MAPPED_OUTPUT_FILE, 'r') as fr:
                for line in fr:
                    if extract_read:
                        print('%s mapped to contaminant genome %s'%(line.strip(), all_categories[i]))
                    extract_read = line.startswith(">")
                    if extract_read:
                        reads.add(line[1:].strip())
    print("Summary report: %d/%d of discovered unmapped reads were successfully mapped to a contaminant"%(len(reads), len(unmapped_reads))) 

if __name__ == "__main__":
    main()
