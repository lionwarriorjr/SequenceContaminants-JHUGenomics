# Sequence Bloom Tree (SBT)

### 1. Dependencies:
* 1) Install jellyfish-2.3.0: https://github.com/gmarcais/Jellyfish
* 2) Install sdsl-lite: https://github.com/simongog/sdsl-lite
* 3) Install gcc (Version 4.9.1 or later): gcc (GCC) 8.3.1 20190223 (Red Hat 8.3.1-2) on CS Ungrad cluster
* 4) Install the latest version of SBT: https://github.com/Kingsford-Group/bloomtree


### 2. Specify package configuration path

Note: Suppose all the dependencies are installed under folder: ./resource/, please write/modify the following into ~/.bashrc file.

vim ~/.bashrc
```
#jellyfish
export path=$PATH:/home/bzhang34/genomics/resource/jellyfish-2.3.0/bin:
## package config
export PKG_CONFIG_PATH=/path/to/final/resource/jellyfish-2.3.0:/path/to/final/resource/sdsl-lite/build  
```

### 3. Build Sequence Bloom Tree based on 11 species
```
./code/sbt/sbt_build.sh
```
We can build the SBT for potential contamination sources, which is 11 species in our case. The following are the key steps in building a bloomtree:

* 1) cd into bloomtree directory to run ./bt and specify relative paths
```
SBT_DIR='./resource/bloomtree/src/'
DATA_DIR='../../../data/raw'
SIM_DIR='../../../data/sim'
BUILD_DIR='../../../build/sbt'
CODE_DIR='../../../code/sbt'
OUTPUT_DIR='../../../result/sbt'

cd $SBT_DIR
```

* 2) Initialize hash functions 
```
./bt hashes --k 20 $BUILD_DIR/myhash.hh 1
```
In our analysis, we compared different `--k` parameter. It specifies k-mer length. By comparing different kmer lengths, ranging from 20 to 40, we found that the best performance is using default 20-mer. User can adjust the parameter value by replacing `--k 20` by `--k Number`  (i.e. set Number = 40).

* 3) Build Bloom Filter for each species using a for loop
```
for f in  $DATA_DIR/*.fna ; do
         ./bt count --k 20 --cutoff 1 $BUILD_DIR/myhash.hh 2000000000  $f  $BUILD_DIR/`basename $f .fna`.bf.bv
done
```
In our analysis, we compared different `--cutoff` parameter. It specifies the minimum count required for a kmer to be added to the bloom filter. In our comparison, we found that the best choice is to include any element with one count or greater. 

* 4) Create a list for all species
```
find $BUILD_DIR/*.bf.bv -maxdepth 1 -type f > $BUILD_DIR/list.txt
```

* 5) Build a Sequence Bloom Tree
```
./bt build --k 20 $BUILD_DIR/myhash.hh $BUILD_DIR/list.txt $BUILD_DIR/allref.bloomtree
```

* 6) Draw the graph of bloom tree structure
```
./bt draw $BUILD_DIR/allref.bloomtree $BUILD_DIR/allrefstructure.dot
```
User can use the online tools or install graphic package to draw the .dot relationship of all nodes and edges.

### 4. Query the unmapped reads obtained from Bowtie2 using SBT and summarize the results
```
./code/sbt/sbt_run.sh
```
After building the bloomtree, we could run the shell script `./sbt_run.sh` to query the unmapped reads to find the potential contaminating source(s). The following steps are included in our pipeline as well.

* 1) cd into bloomtree directory to run ./bt and specify relative paths
```
SBT_DIR='./resource/bloomtree/src/'
DATA_DIR='../../../data/raw'
SIM_DIR='../../../data/sim'
BUILD_DIR='../../../build/sbt'
CODE_DIR='../../../code/sbt'
OUTPUT_DIR='../../../result/sbt'
UNMAPPED_FILE='unmapped_mix3species'

cd $SBT_DIR
```
The `UNMAPPED_FILE` specifies the file name for query set of sequences.

* 2) Preprocess
```
# convert faastq to fasta
python3 $CODE_DIR/fastqtofasta.py $SIM_DIR/$UNMAPPED_FILE.fastq > $SIM_DIR/$UNMAPPED_FILE.fasta
# convert fasta to line-separated reads
python3 $CODE_DIR/fastatoseq.py $SIM_DIR/$UNMAPPED_FILE.fasta > $SIM_DIR/$UNMAPPED_FILE.txt
```

* 3) Query the unmapped reads
```
./bt query -t 0.8 --k 20 $BUILD_DIR/allref.bloomtree $SIM_DIR/$UNMAPPED_FILE.txt $OUTPUT_DIR/output.txt
```
In our analysis, we compared different choice of threshold (`--t`). It defines a hit as the proportion of query kmers must be present in any bloom filter, which sets a tolerance for "approximate matching". We used the default `--t 0.8` as it produced reasonable result in our simulation setting.

* 4) Summarize the result 
```
python3 $CODE_DIR/SBTsummarize.py $SIM_DIR/$UNMAPPED_FILE.fasta $OUTPUT_DIR/output.txt > $OUTPUT_DIR/summary.txt
```
The `SBTsummarize.py` script takes two inputs: first as query file in .fasta format and second as the query output files obtainted by applying SBT. It returns a summary file that contains 3 columns: sequence name, sequence read and 'matched' contaminated species.






