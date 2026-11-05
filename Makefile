.PHONY: install test reference smoke full check

install:
	R CMD INSTALL .

test: install
	Rscript tests/test-core.R
	Rscript examples/quickstart.R
	Rscript examples/client-server.R
	Rscript examples/tuning.R

reference:
	Rscript reproduction/reproduce.R --mode=reference

smoke:
	Rscript reproduction/reproduce.R --mode=smoke

full:
	Rscript reproduction/reproduce.R --mode=full --workers=4

check:
	R CMD build .
	R CMD check --no-manual distkre_0.1.0.tar.gz
