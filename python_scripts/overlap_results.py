#!/usr/bin/env python3
import sys

SPECIES = {0: 'E.coli',
           1: 'S.maltophilia',
           2: 'A.xylosoxidans',
           3: 'D.acidovorans',
           4: 'P.acnes',
           5: 'S.epidermidis',
           6: 'A.baumannii',
           7: 'S.cerevisiae',
           8: 'W.anomalus',
           9: 'M.restricta',
           10: 'P.notatum'}

def main():
    '''
    1. READ IN THE FILE WITH UNMAPPED READS AND STORE THEM IN A LIST
    '''
    file_name = sys.argv[1] # Path to the file with unmapped reads
    with open(file_name) as in_file:
        lines = in_file.readlines()
    names = [] # Fasta name line
    reads = [] # Fasta sequence line
    for line in lines:
        if line.startswith('>'):
            names.append(line.rstrip())
        else:
            reads.append(line.rstrip())

    '''
    2-1. GET RNN PREDICTION FOR EACH READ
         ex. [[1, 3, 6], [4, 5], ..., [1, 3, 6]]
    '''
    rnn_pred = []
    for read in reads:
        # Get results from RNN
        # Append results to rnn_pred

    '''
    2-2. GET BLOOM FILTER PREDICITON FOR EACH READ
         ex. [[1, 3, 6], [4, 5], ..., [1, 3, 6]]
    '''
    bloom_pred = []
    for read in reads:
        # Get results from RNN
        # Append results to bloom_pred

    '''
    3. TAKE INTERSECTIONS OF THE PREDICTIONS
    '''
    intersections = {}
    for i in range(len(reads)):
        intersection = set(rnn_pred[i]).intersection(set(bloom_pred))
        if intersection not in intersections:
            intersections[intersection] = [i]
        else:
            intersections[intersection].append(i)

    '''
    4. WRITE READS WITH THE SAME INTERSECTIONS TO RESPECTIVE FILES
    '''
    for intersection, members in intersections.items():
        with open('_'.join(SPECIES[i] for i in intersection) + '.fasta', 'w') as out_file:
            for member in members:
                out_file.write(names[member] + '\n')
                out_file.write(reads[member] + '\n')

if __name__ == '__main__':
    main()
