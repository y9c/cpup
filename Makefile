#
# Makefile
# Ye Chang, 2020-12-22 01:14
#

all: cpup
	@echo "Done!"

cpup: cpup.cpp
	@g++ -O3 -o $@ $<

.PHONY : test
test: cpup
	@samtools mpileup -d 0 -Q 0 --reverse-del -l ./test/yeast.bed -f ./test/yeast.fa ./test/sample1.bam ./test/sample2.bam | ./$<
