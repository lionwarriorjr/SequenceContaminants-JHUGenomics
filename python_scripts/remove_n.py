import sys


with open(sys.argv[1],'r') as ref:
    with open('chr21_no_n.fa', 'a') as cp:
        n = 'N'
        n_string = ''.join([char*30 for char in n])
        for line in ref:
            if n_string in line.rstrip():
                continue
            else:
                cp.write(line)