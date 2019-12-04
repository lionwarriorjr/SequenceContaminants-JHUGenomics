#!/bin/bash
for filename in build/facs_build/facs_bf/*.bloom; do
  echo "$filename"
  (./facs remove -r "$filename" -q "unmapped_mix3species.fastq") 2> "build/facs_build/contaminated_fastq/$(basename "$filename" .fna).fastq"
done