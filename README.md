# cpup

- convert samtools mpileup result into base count tsv
- multiple bam files are supported
- `bcftools mpileup` might a better choice, if you are not confusing with tags
  (`-t AD`, `-t ADF`, `-t ADR`).

`samtools mpileup` can mpileup the mapping result site by site in the format
below. The 5th (8th, 11st, 14th, ...) column report the observed bases in each site.

```
XII	455422	C	16	<<<<<<<<<<<<<<,,	FFFFFFFFFFFFFFFF	4	<<,,	FFFF
XII	455423	T	16	<<<<<<<<<<<<<<,,	FFFFFFFFFFFFFFFF	4	<<,,	FFFF
XII	455424	C	17	<<<<<<<<<<<<<<,,^$,	FFFFFFFFFFFFFFFFE	4	<<,,	FFFF
XII	455425	A	18	<<<<<<<<<<<<<<,,,^$,	FFFFFFFFFFFFFFFFFF	4	<<,,	FFFF
XII	455426	A	18	<<<<<<<<<<<<<<,,,,	FFFFFFFFFFFFFFFFFF	4	<<,,	FFFF
XII	455427	A	20	<<<<<<<<<<<<<<,,,,^$,^$,	FFFFFFFFFFFFFFFFFFEE	4	<<,,	FFFF
```

`cpup` can convert the mpileup output into table below.

```
chr     pos     ref_base        depth,ref,A,C,G,T,N,Gap,Insert,Delete,a,c,g,t,n,gap,insert,delete       depth,ref,A,C,G,T,N,Gap,Insert,Delete,a,c,g,t,n,gap,insert,delete
XII     455422  C       16,2,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0    4,2,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0
XII     455423  T       16,2,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0    4,2,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0
XII     455424  C       17,3,0,0,0,0,0,0,0,0,0,3,0,0,0,0,0,0    4,2,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0
XII     455425  A       18,4,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0    4,2,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0
XII     455426  A       18,4,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0    4,2,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0
XII     455427  A       20,6,0,0,0,0,0,0,0,0,6,0,0,0,0,0,0,0    4,2,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0
```

- `cpup -H` to hide header.
- `cpup -i` to append indel count in base count sequence.
- `cpup -f` to filter sites.

## Install

```bash
make
```

## Usage

linux pipeline:

```bash
samtools mpileup -d 0 -Q 0 --reverse-del -l <.bed> -f <.fa> <.bam> | cpup
```

## Q&A?

- filter input base by its quality?

`samtools` can do it

- filter output site by cutoff?

`awk` can do it
