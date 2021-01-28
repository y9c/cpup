#
# Makefile
# Ye Chang, 2020-12-22 01:14
#

all: cpup
	@echo "Done!"

cpup: cpup.cpp
	@g++ -o $@ $<

test: cpup
	@samtools mpileup -d 0 -Q 0 --reverse-del -f ./test/yeast.fa ./test/yeast.bam ./test/yeast.bam | ./$<
