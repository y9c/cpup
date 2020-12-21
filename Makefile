#
# Makefile
# Ye Chang, 2020-12-22 01:14
#

all: cpup
	@echo "Makefile needs your attention"

cpup: mpileup2readcounts.cc
	@gcc -o $@ $<
