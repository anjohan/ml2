sources = $(shell find -name "*.f90")
SHELL := /usr/bin/bash
deps = sources2.bib data/J_Ridge1.png data/J_ols1.png data/J_LASSO1.png

all:
	mkdir -p data
	make build
	make -j $(deps)
	make report.pdf

build: $(sources)
	mkdir -p build && cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && make

%.pdf: %.tex $(deps)
	latexmk -pdflua -time -g -shell-escape $*

%.pdf: %.asy
	asy -maxtile "(400,400)" -o $@ $<

sources2.bib: sources.bib
	betterbib -l $< $@

data/J_%1.png: ./programs/matrix_to_png.py data/J_%.dat
	python $< $*

data/J_ols.dat: build/linreg
	./$<

data/J_Ridge.dat: data/J_ols.dat

data/J_LASSO.dat: programs/lasso.py
	python $<

build/%: build programs/%.f90

clean:
	latexmk -c
	rm -rf *.run.xml *.bbl build
	find -name "*.mod" -delete
