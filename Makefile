# Get the version info for later
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)

all: docs check clean

docs:
	R -q -e 'library("roxygen2"); roxygenise(".")'

build: docs
	cd ..;\
	R CMD build temporalEF

check: build
	cd ..;\
	R CMD check temporalEF_$(PKGVERS).tar.gz

check-cran: build
	cd ..;\
	R CMD check --as-cran temporalEF_$(PKGVERS).tar.gz

install: build
	cd ..;\
	R CMD INSTALL temporalEF_$(PKGVERS).tar.gz

move: check
	cp ../temporalEF.Rcheck/temporalEF-Ex.Rout ./tests/Examples/temporalEF-Ex.Rout.save

clean:
	cd ..;\
	rm -r temporalEF.Rcheck/

purl:
	cd ..;\
	R -q -e 'library("knitr"); purl("./temporalEF/vignettes/temporalEF.Rnw")'
