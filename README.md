# cpup

## Install

```bash
make
```

## Usage

```bash
samtools mpileup -d 0 -Q 0 --reverse-del -l <.bed> -f <.fa> <.bam> | cpup
```
