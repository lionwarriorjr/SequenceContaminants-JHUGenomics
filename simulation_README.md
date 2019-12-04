# Simulate fastq with reference fasta

1. We want to first remove lines that consists all N bases
from human chr21 fasta. Run the following python script.
```shell script
python3 remove_n.py /input/file/path/ /output/file/path/
```  
2. Now we can simulate fastq from reference fasta for
each species using mason. 
Run the following script for each species fasta.  
```shell script
./mason2-2.0.9-Mac-x86_64/bin/mason_simulator -ir /path/to/fasta/ -n 100000 -o /output/fastq --illumina-read-lengh 150
```
This will simulate fastq of 100,000 reads from reference fasta file with the
following default variations:
- Insertion rate: 0.00005
- Deletion rate: 0.00005
- Average Probability for mismatch: 0.004  

3. Now we can create simulated contaminated human fastq
by appending a set proportion of bacteria fastq reads to 
the human fastq file. We set 1% of contamination for each species.
This means that for 100,000 human reads, we append 1000 reads 
for each bacteria. The mixture.py script takes in the human fastq
file as the first input and a txt file of the contaminated proportion 
for each bacteria then output the contaminated fastq.
Here we have a sameprop.txt that has 3 bacteria species of the 
same 1% proportion. The diffprop.txt has 3 bacteria with 1%, 
2%, 3% proportion.
```shell script
python3 mixture.py /path/to/human/fastq ../../data/sameprop.txt > /path/to/sim_contaminated_output_fastq
python3 mixture.py /path/to/human/fastq ../../data/diffprop.txt > /path/to/sim_contaminated_output_fastq
```

4. The next step is to build the bowtie2 index for the
 reference human fasta. This will create index files with 
 prefix of ref_chr21
```shell script
bowtie2-build path/to/human_ref_chr21.fasta ref_chr21
```

5. Align the contaminated fastq from step 3 to human reference 
using bowtie2 and extract the unmapped reads.
```shell script
bowtie2 -x path/to/ref_chr21 -U path/to/chr21_mix3species.fastq --un path/to/unmapped_mix3species.fastq -S path/to/out_mix3species.sam
bowtie2 -x path/to/ref_chr21 -U path/to/chr21_mix3species_diffprop.fastq --un path/to/unmapped_mix3species_diffprop.fastq -S path/to/out_mix3species_diffprop.sam
```
