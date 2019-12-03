#!/bin/bash
read -p "Enter fastq path: " fastq_path
read -p "Output path: " output_path
cd facs/facs
for filename in ../../facs_bf/*.bloom; do
  echo "$filename"
  (./facs remove -r "$filename" -q "$fastq_path") 2> "../../contaminated_fastq/$(basename "$filename" .fna).fastq"
done