.PHONY: check all

check:
	Rscript run_all.R --mode=check

all:
	Rscript run_all.R --mode=full

