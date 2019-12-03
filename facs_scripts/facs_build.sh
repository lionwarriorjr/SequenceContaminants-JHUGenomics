#!/bin/bash
cd facs/facs
make
for filename in ../../../ref/*.fna; do
  echo "$filename"
  ./facs build -r "$filename" -o "../../facs_bf/$(basename "$filename" .fna).bloom"
doner