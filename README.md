# cpup

forked from mpileup2readcounts

## Install

```bash
mkdir build
cd build
cmake ..
make
```

## Usage

```bash
samtools mpileup -f ref.fa -l regions.bed -d 0 -Q 0 --reverse-del alignments.bam | cpup sample1
```
