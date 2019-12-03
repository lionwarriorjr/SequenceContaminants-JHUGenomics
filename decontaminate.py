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
from dna2vec.dna2vec.multi_k_model import MultiKModel
from collections import Counter, defaultdict


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
BOWTIE_ALIGNMENTS_INPUT = "cont.sam"
TABLE_LEN = 13
UNMAPPED_STATUS_COL = 1
UNMAPPED_FLAG = 4
MAPPED_FLAG = 0
READ_COL = 9
# how many of unmapped reads should be sampled for verification with Bowtie
ALIGNMENT_VALIDATION_LIMIT = 1000

# Constants for RNN evaluation

# threshold probability of max target class in RNN prediction
RNN_THRESHOLD = 0.5
# identities of contaminant sources with pre-built Bowtie indexes
all_categories = ['D. acidovorans', 'S. epidermidis', 'W. anomalus']
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
rnn = torch.load('./rnn.pt')


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


# write reads to fasta file on disk 
def write_fasta(filename, reads, header=None, label=DEFAULT_HEADER, limit=150, append=False):
    index = 0
    mode = "a" if append else "w"
    with open(filename, mode) as fw:
        if header:
            fw.write(header + "\n")
        start = 0
        for read in reads:
            offset = 0
            while offset + limit - start < len(read):
                # if next substring of read terminates the line, write it
                fw.write(read[offset:offset+limit-start] + "\n")
                offset += limit-start
                index += 1
                if header:
                    fw.write(label + " read %d:\n"%(index))
                start = 0
            if start > 0 and limit - start < len(read):
                # start the next line written to file
                fw.write(read[:limit-start] + "\n")
                index += 1
                if header:
                    fw.write(label + " read %d:\n"%(index))
                offset, start = limit-start, 0
            # continue writes for current line being written
            fw.write(read[offset:])
            start += len(read)-offset


def contaminate(src_reads, contaminant):
    # introduce contaminants into source reads during simuulation
    n_contaminants = int(CONTAMINATION_RATE * len(src_reads))
    print("introducing %d contaminants from contaminant %s"%(n_contaminants, contaminant))
    contaminant_reads = parse_fasta(contaminant, READ_LENGTH, n_contaminants)
    src_reads_copy = src_reads.copy()
    src_reads_copy.extend(contaminant_reads)
    return src_reads_copy


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
def predict(input_line, n_predictions=1):
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
            # extract max class for estimate
            category_index = topi[0][i].item()
            # extract target propbability for the max estimate
            max_prob = prob[0,category_index].numpy()
            predictions.append((all_categories[category_index], max_prob))
    return predictions


def main():

    # name of pre-built index for source genome in current directory
    index = sys.argv[1]
    # name of file with intermixed contaminant reads
    contaminated_file = sys.argv[2]
    # how many of the unmapped reads should be verified by Bowtie2
    n_args = len(sys.argv)-1
    # optional third command-line argument
    alignment_limit = sys.argv[3] if n_args == 3 else ALIGNMENT_VALIDATION_LIMIT
    
    # run bowtie alignment to find mapped/unmapped reads
    process = subprocess.Popen( args=["bowtie2", 
                                  "-x", index, 
                                  "-f", contaminated_file, 
                                  "-S", BOWTIE_ALIGNMENTS_INPUT ] )
    # block until bowtie subprocess completes for initial alignment
    process.wait()
    
    # extract unmapped reads from Bowtie2 output
    print("Logging unmapped reads from Bowtie2 alignment of contaminated genome")
    unmapped_reads = []
    with open(BOWTIE_ALIGNMENTS_INPUT, 'r') as fr:
        for line in fr:
            values = line.split("\t")
            if len(values) > READ_COL:
                if int(values[UNMAPPED_STATUS_COL]) == UNMAPPED_FLAG:
                    unmapped_read = values[READ_COL]
                    print('unmapped read: %s'%(unmapped_read))
                    unmapped_reads.append(unmapped_read.strip())

    print("loaded pretrained DNA word embeddings and trained RNN model from disk\n")
    print(rnn)

    # forward propagate through RNN to get predictions for each read
    rnn_pred = defaultdict(list)
    for unmapped_read in unmapped_reads:
        # run the unmapped read through the RNN
        category, prob = predict(unmapped_read)[0]
        # only include predictions for contaminants with
        # a probability estimate from the RNN > RNN_THRESHOLD
        if prob > RNN_THRESHOLD:
            rnn_pred[unmapped_read].append(category)
    for read in rnn_pred:
        # get the best estimate for the contaminant source
        max_category = Counter(rnn_pred[read]).most_common(1)[0]
        rnn_pred[read] = max_category[0]

    # initialize names of pre-built Bowtie indexes in current directory
    indexes = {
        'D. acidovorans': 'acidovorans',
        'S. epidermidis': 'epidermidis',
        'W. anomalus': 'anomalus'
    }

    overlapped_reads = rnn_pred
    # overlapped classifications from DL and Bloom Filter are our predictions
    print('predicted set of unmapped reads mapped to contaminant genomes:')
    for read in overlapped_reads:
        print('%s introduced by contaminant genome %s'%(read, overlapped_reads[read]))

    # sample verification of first ALIGNMENT_VALIDATION_LIMIT unmapped reads
    print("Verifying unmapped reads against contaminant genome with Bowtie2 alignment")
    count = 0
    reads, mapped = [], set()
    for read in overlapped_reads:
        reads.append(read)
        # extract the pre-built Bowtie index for the predicted contaminant
        index = indexes[overlapped_reads[read]]
        write_fasta(filename=OVERLAP_OUTPUT_FILE, 
                    reads=[read], 
                    header=OVERLAP_HEADER, 
                    label=">%s"%(OVERLAP_OUTPUT_FILE),
                    limit=LINE_LENGTH)
        # run the verification alignment
        process = subprocess.Popen( args=["bowtie2", 
                                    "-x", index, 
                                    "-f", OVERLAP_OUTPUT_FILE, 
                                    "-S", BOWTIE_ALIGNMENTS_INPUT ],
                                    stdout=DEVNULL, stderr=STDOUT )
        process.wait()
        with open(BOWTIE_ALIGNMENTS_INPUT, 'r') as fr:
            for line in fr:
                values = line.split("\t")
                if len(values) > READ_COL:
                    # check if alignment is a match against the contaminant genome
                    if int(values[UNMAPPED_STATUS_COL]) == MAPPED_FLAG:
                        mapped.add(read)
                        print('%s mapped to %s'%(values[READ_COL], overlapped_reads[read]))
                        break
        count += 1
        if count == alignment_limit:
            break

    # find which of the first ALIGNMENT_VALIDATION_LIMIT unmapped reads actually 
    aligned = 0
    for read in reads:
        if read in mapped:
            aligned += 1
    print("%d/%d of the first %d unmapped reads generated from a contaminant genome"%(aligned, len(reads), len(reads)))

if __name__ == "__main__":
    main()
