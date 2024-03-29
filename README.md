# cpup

**Convert samtools mpileup result into base count table**

> Feature:

- multiple bam files are supported
- filter sites before output
- output by strandness
- group and count indels by pattern

## Install

```bash
git clone https://github.com/y9c/cpup.git
cd ./cpup
make
```

## Usage

Pipe the output of samtools mpileup into this tools!

```
Usage:
  samtools mpileup -d 0 -Q 10 --reverse-del -l <.bed> -f <.fa> <.bam> | cpup

  -h, --help          show help
  -H, --headerless    hide header
  -S, --strandless    ignore strand information
  -s, --by-strand     output by strand
  -m, --major-strand  output major strand only
  -i, --indel         append indel count
  -e, --ends          append read ends (5' | 3') count
  -c, --count []      select count columns
  -f, --filter []     filter sites
  -F, --drop []       drop sites

```

`samtools mpileup` can mpileup the mapping result site by site in the format
below. The 5th (8th, 11st, 14th, ...) column report the observed bases in each site.

```
XII     455422  C       16      <<<<<<<<<<<<<<,,        FFFFFFFFFFFFFFFF        4       <<,,    FFFF
XII     455423  T       16      <<<<<<<<<<<<<<,,        FFFFFFFFFFFFFFFF        4       <<,,    FFFF
XII     455424  C       17      <<<<<<<<<<<<<<,,^$,     FFFFFFFFFFFFFFFFE       4       <<,,    FFFF
XII     455425  A       18      <<<<<<<<<<<<<<,,,^$,    FFFFFFFFFFFFFFFFFF      4       <<,,    FFFF
XII     455426  A       18      <<<<<<<<<<<<<<,,,,      FFFFFFFFFFFFFFFFFF      4       <<,,    FFFF
XII     455427  A       20      <<<<<<<<<<<<<<,,,,^$,^$,        FFFFFFFFFFFFFFFFFFEE    4       <<,,    FFFF
XII     455428  C       20      <<<<<<<<<<<<<<,,,,,,    FFFFFFFFFFFFFFFFFFFF    4       <<,,    FFFF
XII     455429  A       20      <<<<<<<<<<<<<<,,,,,,    FFFFFFFFFFFFFFFFFFFF    4       <<,,    FFFF
XII     455430  G       22      <<<<<<<<<<<<<<,,,,,,^$,^$,      FFFFFFFFFFFFFFFFFFFFBB  4       <<,,    FFFF
XII     455431  G       23      <<<<<<<<<<<<<<,,,,,,,,^$,       FFFFFFFFFFFFFFFFFFFFEED 4       <<,,    FFFF
...
```

`cpup` can convert the mpileup output into table below by parameters `-i -s -f mut:3`.

```
chr     pos     ref_base        strand  depth,a,c,g,t,n,gap,insert,delete,istat,dstat           depth,a,c,g,t,n,gap,insert,delete,istat,dstat
XII     455496  C       -       1409,17,0,0,0,0,2,0,2,,c:1|caa:1        139,3,0,0,0,0,1,0,2,,c:1|caa:1
XII     455498  T       -       1454,0,6,2,0,0,407,0,244,,at:220|att:24 144,0,0,1,0,0,37,0,27,,at:23|att:4
XII     455499  T       -       1496,0,1,12,0,0,250,0,12,,t:7|tt:5      151,0,0,2,0,0,29,0,1,,t:1
XII     455500  A       -       1455,0,2,2,0,0,256,0,0,,        145,0,1,0,0,0,28,0,0,,
XII     729179  T       +       223,0,3,0,0,0,33,0,0,,  56,0,0,0,0,0,14,0,0,,
XII     729182  C       +       235,1,0,0,4,0,2,0,0,,   56,0,0,0,0,0,1,0,0,,
```

> NOTE:

- `-Q` in samtools mpileup should not be set to zero, which might bring bug in counting overlapping reads.
- `-f` to check **any** (max) value greater than cutoff.
  eg: `mut:3` for filter sites with **more than or equal to(>=)** 3 mutations in any sample.
  You can set multile filters, like `mut:3,delete:2` to filter sites with >= 3 mutations in any samples, meanwhile, are should be >=2 delete events.
- `-F` to check **all** (min) value greater than cutoff.

## Q&A?

- filter input base by its quality?

`samtools` can do it

- filter output site by cutoff?

`awk` can do it

- samtools mpileup is too slow?

`bcftools mpileup` might a better choice, if you are not confusing with tags (`-t AD`, `-t ADF`, `-t ADR`).
